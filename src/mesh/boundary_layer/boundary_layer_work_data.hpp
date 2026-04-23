// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#pragma once

#include "sqmesh/mesh/api.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace sqmesh::mesh {

// -- Vec3 utilities ------------------------------------------------------
using Vec3 = std::array<double, 3>;

inline Vec3 vec3_sub(const Vec3 &a, const Vec3 &b) noexcept
{
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

inline Vec3 vec3_add(const Vec3 &a, const Vec3 &b) noexcept
{
  return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

inline Vec3 vec3_scale(const Vec3 &a, double s) noexcept
{
  return {a[0] * s, a[1] * s, a[2] * s};
}

inline double vec3_dot(const Vec3 &a, const Vec3 &b) noexcept
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline Vec3 vec3_cross(const Vec3 &a, const Vec3 &b) noexcept
{
  return {
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0]
  };
}

inline double vec3_length(const Vec3 &a) noexcept
{
  return std::sqrt(vec3_dot(a, a));
}

inline Vec3 vec3_normalized(const Vec3 &a) noexcept
{
  const double len = vec3_length(a);
  if(len < 1e-30) return {0.0, 0.0, 0.0};
  return vec3_scale(a, 1.0 / len);
}

// -- AABB for spatial queries --------------------------------------------
struct AABB {
  Vec3 lo {1e30, 1e30, 1e30};
  Vec3 hi {-1e30, -1e30, -1e30};

  void expand(const Vec3 &p) noexcept
  {
    for(int i = 0; i < 3; ++i) {
      lo[i] = std::min(lo[i], p[i]);
      hi[i] = std::max(hi[i], p[i]);
    }
  }

  void expand(const AABB &other) noexcept
  {
    for(int i = 0; i < 3; ++i) {
      lo[i] = std::min(lo[i], other.lo[i]);
      hi[i] = std::max(hi[i], other.hi[i]);
    }
  }

  void pad(double margin) noexcept
  {
    for(int i = 0; i < 3; ++i) {
      lo[i] -= margin;
      hi[i] += margin;
    }
  }

  [[nodiscard]] bool overlaps(const AABB &other) const noexcept
  {
    return lo[0] <= other.hi[0] && hi[0] >= other.lo[0] &&
           lo[1] <= other.hi[1] && hi[1] >= other.lo[1] &&
           lo[2] <= other.hi[2] && hi[2] >= other.lo[2];
  }

  [[nodiscard]] Vec3 center() const noexcept
  {
    return {(lo[0]+hi[0])*0.5, (lo[1]+hi[1])*0.5, (lo[2]+hi[2])*0.5};
  }
};

// -- Simple BVH for face intersection queries ----------------------------
struct BVHNode {
  AABB box;
  std::uint32_t left = 0;   // child index or face index if leaf
  std::uint32_t right = 0;
  bool is_leaf = false;
};

class FaceBVH {
public:
  void build(const std::vector<Vec3> &nodes,
             const std::vector<std::array<std::uint32_t, 4>> &faces,
             const std::vector<int> &face_nv)
  {
    nodes_ = &nodes;
    faces_ = &faces;
    face_nv_ = &face_nv;
    bvh_nodes_.clear();

    const auto nf = static_cast<std::uint32_t>(faces.size());
    if(nf == 0) return;

    std::vector<std::uint32_t> indices(nf);
    for(std::uint32_t i = 0; i < nf; ++i) indices[i] = i;

    bvh_nodes_.reserve(nf * 2);
    build_recursive(indices.data(), nf);
  }

  // Query: find all faces whose AABB overlaps query_box
  void query_overlap(const AABB &query_box,
                     std::vector<std::uint32_t> &result) const
  {
    result.clear();
    if(bvh_nodes_.empty()) return;
    query_recursive(0, query_box, result);
  }

  // Ray-triangle intersection: returns distance or negative if no hit.
  // Skips faces in skip_set (topological neighbors).
  [[nodiscard]] double ray_nearest_intersection(
    const Vec3 &origin,
    const Vec3 &direction,
    double max_dist,
    const std::vector<bool> &skip_face
  ) const
  {
    if(bvh_nodes_.empty()) return -1.0;
    double best = max_dist;
    ray_query_recursive(0, origin, direction, best, skip_face);
    return (best < max_dist) ? best : -1.0;
  }

private:
  std::uint32_t build_recursive(std::uint32_t *indices, std::uint32_t count)
  {
    const auto node_idx = static_cast<std::uint32_t>(bvh_nodes_.size());
    bvh_nodes_.push_back({});

    // Compute AABB
    for(std::uint32_t i = 0; i < count; ++i) {
      const auto fi = indices[i];
      const auto &f = (*faces_)[fi];
      const int nv = (*face_nv_)[fi];
      for(int v = 0; v < nv; ++v) {
        bvh_nodes_[node_idx].box.expand((*nodes_)[f[v]]);
      }
    }

    if(count <= 2) {
      bvh_nodes_[node_idx].is_leaf = true;
      bvh_nodes_[node_idx].left = indices[0];
      bvh_nodes_[node_idx].right = (count > 1) ? indices[1] : indices[0];
      return node_idx;
    }

    // Split along longest axis. Use a median split instead of a simple
    // midpoint partition so the recursion depth stays bounded even when many
    // faces share the same centroid coordinate on the chosen axis.
    int axis = 0;
    double span = bvh_nodes_[node_idx].box.hi[0] - bvh_nodes_[node_idx].box.lo[0];
    for(int a = 1; a < 3; ++a) {
      const double s = bvh_nodes_[node_idx].box.hi[a] - bvh_nodes_[node_idx].box.lo[a];
      if(s > span) { span = s; axis = a; }
    }

    auto centroid_on_axis = [&](std::uint32_t fi) noexcept {
      const auto &f = (*faces_)[fi];
      const int nv = (*face_nv_)[fi];
      double sum = 0.0;
      for(int v = 0; v < nv; ++v) {
        sum += (*nodes_)[f[v]][axis];
      }
      return sum / static_cast<double>(nv);
    };

    const auto left_count = count / 2;
    std::nth_element(
      indices,
      indices + left_count,
      indices + count,
      [&](std::uint32_t lhs, std::uint32_t rhs) noexcept {
        return centroid_on_axis(lhs) < centroid_on_axis(rhs);
      });

    const auto left_child = build_recursive(indices, left_count);
    const auto right_child = build_recursive(indices + left_count, count - left_count);
    // Child recursion may reallocate bvh_nodes_, so only write back by index.
    bvh_nodes_[node_idx].left = left_child;
    bvh_nodes_[node_idx].right = right_child;
    return node_idx;
  }

  void query_recursive(std::uint32_t idx, const AABB &query,
                       std::vector<std::uint32_t> &result) const
  {
    const auto &node = bvh_nodes_[idx];
    if(!node.box.overlaps(query)) return;

    if(node.is_leaf) {
      result.push_back(node.left);
      if(node.right != node.left) result.push_back(node.right);
      return;
    }

    query_recursive(node.left, query, result);
    query_recursive(node.right, query, result);
  }

  // Möller–Trumbore ray-triangle intersection
  static double ray_tri_intersect(
    const Vec3 &orig, const Vec3 &dir,
    const Vec3 &v0, const Vec3 &v1, const Vec3 &v2
  ) noexcept
  {
    constexpr double eps = 1e-12;
    auto e1 = vec3_sub(v1, v0);
    auto e2 = vec3_sub(v2, v0);
    auto h = vec3_cross(dir, e2);
    double a = vec3_dot(e1, h);
    if(std::abs(a) < eps) return -1.0;
    double f = 1.0 / a;
    auto s = vec3_sub(orig, v0);
    double u = f * vec3_dot(s, h);
    if(u < 0.0 || u > 1.0) return -1.0;
    auto q = vec3_cross(s, e1);
    double v = f * vec3_dot(dir, q);
    if(v < 0.0 || u + v > 1.0) return -1.0;
    double t = f * vec3_dot(e2, q);
    return (t > eps) ? t : -1.0;
  }

  void ray_query_recursive(
    std::uint32_t idx,
    const Vec3 &origin,
    const Vec3 &dir,
    double &best,
    const std::vector<bool> &skip_face
  ) const
  {
    const auto &node = bvh_nodes_[idx];
    // Simple AABB-ray test (slab method)
    double tmin = 0.0, tmax = best;
    for(int a = 0; a < 3; ++a) {
      if(std::abs(dir[a]) < 1e-30) {
        if(origin[a] < node.box.lo[a] || origin[a] > node.box.hi[a])
          return;
      } else {
        double inv = 1.0 / dir[a];
        double t1 = (node.box.lo[a] - origin[a]) * inv;
        double t2 = (node.box.hi[a] - origin[a]) * inv;
        if(t1 > t2) std::swap(t1, t2);
        tmin = std::max(tmin, t1);
        tmax = std::min(tmax, t2);
        if(tmin > tmax) return;
      }
    }

    if(node.is_leaf) {
      auto test_face = [&](std::uint32_t fi) {
        if(fi < skip_face.size() && skip_face[fi]) return;
        const auto &f = (*faces_)[fi];
        const int nv = (*face_nv_)[fi];
        // Tri 0
        double t = ray_tri_intersect(origin, dir,
          (*nodes_)[f[0]], (*nodes_)[f[1]], (*nodes_)[f[2]]);
        if(t > 0 && t < best) best = t;
        // Tri 1 (quad second half)
        if(nv == 4) {
          t = ray_tri_intersect(origin, dir,
            (*nodes_)[f[0]], (*nodes_)[f[2]], (*nodes_)[f[3]]);
          if(t > 0 && t < best) best = t;
        }
      };
      test_face(node.left);
      if(node.right != node.left) test_face(node.right);
      return;
    }

    ray_query_recursive(node.left, origin, dir, best, skip_face);
    ray_query_recursive(node.right, origin, dir, best, skip_face);
  }

  const std::vector<Vec3> *nodes_ = nullptr;
  const std::vector<std::array<std::uint32_t, 4>> *faces_ = nullptr;
  const std::vector<int> *face_nv_ = nullptr;
  std::vector<BVHNode> bvh_nodes_;
};

// -- Node type classification --------------------------------------------
enum class NodeTypeBL { FLAT, CONCAVE, CONVEX, CONCAVE_CONVEX, OTHER };

// -- Surface face (tri or quad) ------------------------------------------
struct SurfaceFace {
  std::array<std::uint32_t, 4> vertices {0, 0, 0, 0};
  int num_vertices = 3; // 3=tri, 4=quad
  std::uint32_t topology_owner = 0; // Source face entity_group ID
};

struct LayerSurfaceFace {
  std::array<std::uint32_t, 4> vertices {0, 0, 0, 0};
  int num_vertices = 3;
  int layer = 0;
  std::uint32_t source_face = std::numeric_limits<std::uint32_t>::max();
};

// -- Boundary layer working data ----------------------------------------─
struct BLWorkData {
  // -- Input surface --
  std::vector<Vec3> surface_nodes;
  std::vector<SurfaceFace> surface_faces;
  std::vector<EntityRef> surface_node_refs; // Map back to Domain node refs

  // -- Node state --
  std::vector<bool> node_active;
  std::vector<bool> is_boundary_node;
  std::vector<NodeTypeBL> node_type;
  std::vector<Vec3> march_direction;
  std::vector<int> max_node_level;
  std::vector<double> node_thickness;      // Accumulated thickness per node
  std::vector<std::vector<double>> node_dihedral_angles;
  std::vector<double> node_manifold_angle; // Manifold opening angle

  // -- Face state --
  std::vector<bool> face_active;
  std::vector<Vec3> face_normal;
  std::vector<int> max_face_level;

  // -- Connectivity --
  std::vector<std::vector<std::uint32_t>> node_to_faces;  // VF
  std::vector<std::vector<std::uint32_t>> node_to_nodes;  // VV
  std::vector<std::vector<std::uint32_t>> face_to_faces;  // FF (via shared edge)

  // -- Output: layer positions --
  // layer_positions[layer][node_index] = 3D position
  // layer 0 = original surface
  std::vector<std::vector<Vec3>> layer_positions;
  std::vector<LayerSurfaceFace> top_cap_faces;
  std::vector<LayerSurfaceFace> transition_side_faces;

  // -- BVH for collision detection --
  FaceBVH face_bvh;

  // -- Non-boundary face indices (for collision queries) --
  std::vector<std::uint32_t> non_bl_face_indices;
  std::vector<Vec3> non_bl_nodes;
  std::vector<EntityRef> non_bl_node_refs;
  std::vector<std::array<std::uint32_t, 4>> non_bl_face_verts;
  std::vector<int> non_bl_face_nv;
  FaceBVH non_bl_bvh;

  void clear()
  {
    surface_nodes.clear();
    surface_faces.clear();
    surface_node_refs.clear();
    node_active.clear();
    is_boundary_node.clear();
    node_type.clear();
    march_direction.clear();
    max_node_level.clear();
    node_thickness.clear();
    node_dihedral_angles.clear();
    node_manifold_angle.clear();
    face_active.clear();
    face_normal.clear();
    max_face_level.clear();
    node_to_faces.clear();
    node_to_nodes.clear();
    face_to_faces.clear();
    layer_positions.clear();
    top_cap_faces.clear();
    transition_side_faces.clear();
    non_bl_face_indices.clear();
    non_bl_nodes.clear();
    non_bl_node_refs.clear();
    non_bl_face_verts.clear();
    non_bl_face_nv.clear();
  }

  void build_connectivity()
  {
    const auto num_nodes = surface_nodes.size();
    const auto num_faces = surface_faces.size();

    node_to_faces.assign(num_nodes, {});
    node_to_nodes.assign(num_nodes, {});

    for(std::uint32_t fi = 0; fi < num_faces; ++fi) {
      const auto &face = surface_faces[fi];
      for(int vi = 0; vi < face.num_vertices; ++vi) {
        const auto v = face.vertices[vi];
        node_to_faces[v].push_back(fi);
        const auto v_next = face.vertices[(vi + 1) % face.num_vertices];
        auto &nn = node_to_nodes[v];
        if(std::find(nn.begin(), nn.end(), v_next) == nn.end()) {
          nn.push_back(v_next);
        }
        auto &nn2 = node_to_nodes[v_next];
        if(std::find(nn2.begin(), nn2.end(), v) == nn2.end()) {
          nn2.push_back(v);
        }
      }
    }

    // Build FF connectivity (faces sharing an edge)
    face_to_faces.assign(num_faces, {});
    for(std::uint32_t fi = 0; fi < num_faces; ++fi) {
      const auto &face = surface_faces[fi];
      for(int ei = 0; ei < face.num_vertices; ++ei) {
        const auto v0 = face.vertices[ei];
        const auto v1 = face.vertices[(ei + 1) % face.num_vertices];
        for(const auto neighbor_fi : node_to_faces[v0]) {
          if(neighbor_fi == fi) continue;
          const auto &nf = surface_faces[neighbor_fi];
          for(int nej = 0; nej < nf.num_vertices; ++nej) {
            const auto nv0 = nf.vertices[nej];
            const auto nv1 = nf.vertices[(nej + 1) % nf.num_vertices];
            if((nv0 == v0 && nv1 == v1) || (nv0 == v1 && nv1 == v0)) {
              auto &ff = face_to_faces[fi];
              if(std::find(ff.begin(), ff.end(), neighbor_fi) == ff.end()) {
                ff.push_back(neighbor_fi);
              }
            }
          }
        }
      }
    }
  }
};

} // namespace sqmesh::mesh
