// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "region_detector.hpp"
#include "../boundary_layer/boundary_layer_work_data.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <queue>
#include <unordered_map>
#include <unordered_set>

namespace sqmesh::mesh {

namespace {

using Vec3 = std::array<double, 3>;

Vec3 v3sub(const Vec3 &a, const Vec3 &b) { return {a[0]-b[0], a[1]-b[1], a[2]-b[2]}; }
Vec3 v3add(const Vec3 &a, const Vec3 &b) { return {a[0]+b[0], a[1]+b[1], a[2]+b[2]}; }
Vec3 v3scale(const Vec3 &a, double s) { return {a[0]*s, a[1]*s, a[2]*s}; }
Vec3 v3cross(const Vec3 &a, const Vec3 &b) {
  return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
}
double v3dot(const Vec3 &a, const Vec3 &b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
double v3len(const Vec3 &a) { return std::sqrt(v3dot(a, a)); }
Vec3 v3norm(const Vec3 &a) { double l = v3len(a); return l > 1e-30 ? v3scale(a, 1.0/l) : Vec3{0,0,0}; }

// Möller-Trumbore ray-triangle, returns t or -1
double ray_tri(const Vec3 &orig, const Vec3 &dir,
               const Vec3 &v0, const Vec3 &v1, const Vec3 &v2)
{
  constexpr double eps = 1e-12;
  auto e1 = v3sub(v1, v0);
  auto e2 = v3sub(v2, v0);
  auto h = v3cross(dir, e2);
  double a = v3dot(e1, h);
  if(std::abs(a) < eps) return -1.0;
  double f = 1.0 / a;
  auto s = v3sub(orig, v0);
  double u = f * v3dot(s, h);
  if(u < 0.0 || u > 1.0) return -1.0;
  auto q = v3cross(s, e1);
  double v = f * v3dot(dir, q);
  if(v < 0.0 || u + v > 1.0) return -1.0;
  double t = f * v3dot(e2, q);
  return (t > eps) ? t : -1.0;
}

// Count all ray-triangle intersections (for even-odd test)
int ray_count_all(const Vec3 &origin, const Vec3 &dir,
                  const std::vector<Vec3> &nodes,
                  const std::vector<std::array<std::uint32_t, 3>> &tris,
                  const FaceBVH &bvh)
{
  // Simple brute-force for correctness (BVH acceleration can be added later)
  int count = 0;
  for(const auto &tri : tris) {
    double t = ray_tri(origin, dir, nodes[tri[0]], nodes[tri[1]], nodes[tri[2]]);
    if(t > 0) count++;
  }
  return count;
}

// Determine if a point is inside any enclosed region (even-odd rule)
// Cast ray in +X direction, count crossings
bool point_inside(const Vec3 &point,
                  const std::vector<Vec3> &nodes,
                  const std::vector<std::array<std::uint32_t, 3>> &tris)
{
  // Use 3 ray directions to be robust against edge/vertex hits
  Vec3 dirs[3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
  int votes = 0;
  for(int d = 0; d < 3; ++d) {
    int count = 0;
    for(const auto &tri : tris) {
      double t = ray_tri(point, dirs[d], nodes[tri[0]], nodes[tri[1]], nodes[tri[2]]);
      if(t > 0) count++;
    }
    if(count % 2 == 1) votes++;
  }
  return votes >= 2; // Majority vote
}

struct ExtractedMesh {
  std::vector<Vec3> nodes;
  std::vector<std::array<std::uint32_t, 3>> tris; // Triangle vertex indices
  std::vector<EntityRef> face_refs;                // Back to Domain
  std::vector<Vec3> face_normals;
  std::vector<bool> face_flipped;                  // True if winding was flipped for consistency

  // Adjacency
  std::vector<std::vector<std::uint32_t>> node_to_faces;
  std::vector<std::vector<std::uint32_t>> face_to_faces;
};

void extract_mesh(const Domain &domain, ExtractedMesh &mesh)
{
  std::unordered_map<std::uint64_t, std::uint32_t> node_map;
  auto pack_ref = [](EntityRef ref) -> std::uint64_t {
    return (static_cast<std::uint64_t>(ref.entity_group) << 32U) | ref.index;
  };

  for(const auto &entity_group : domain.entity_groups()) {
    if(entity_group.order() != EntityOrder::face ||
       entity_group.role() != EntityGroupRole::computational) continue;

    for(std::uint32_t fi = 0; fi < entity_group.faces().size(); ++fi) {
      EntityRef face_ref = {entity_group.id(), fi};
      auto fnodes = domain.face_nodes(face_ref);
      if(fnodes.size != 3) continue; // Only triangles

      std::array<std::uint32_t, 3> tri_idx;
      for(std::size_t ni = 0; ni < 3; ++ni) {
        auto key = pack_ref(fnodes[ni]);
        auto it = node_map.find(key);
        if(it == node_map.end()) {
          const auto &n = domain.node(fnodes[ni]);
          auto idx = static_cast<std::uint32_t>(mesh.nodes.size());
          mesh.nodes.push_back({n.coordinates[0], n.coordinates[1], n.coordinates[2]});
          node_map[key] = idx;
          tri_idx[ni] = idx;
        } else {
          tri_idx[ni] = it->second;
        }
      }
      mesh.tris.push_back(tri_idx);
      mesh.face_refs.push_back(face_ref);
    }
  }
}

void compute_normals(ExtractedMesh &mesh)
{
  mesh.face_normals.resize(mesh.tris.size());
  for(std::size_t i = 0; i < mesh.tris.size(); ++i) {
    const auto &t = mesh.tris[i];
    auto e1 = v3sub(mesh.nodes[t[1]], mesh.nodes[t[0]]);
    auto e2 = v3sub(mesh.nodes[t[2]], mesh.nodes[t[0]]);
    mesh.face_normals[i] = v3norm(v3cross(e1, e2));
  }
}

void build_adjacency(ExtractedMesh &mesh)
{
  const auto nn = mesh.nodes.size();
  const auto nf = mesh.tris.size();

  mesh.node_to_faces.assign(nn, {});
  for(std::uint32_t fi = 0; fi < nf; ++fi) {
    for(int vi = 0; vi < 3; ++vi) {
      mesh.node_to_faces[mesh.tris[fi][vi]].push_back(fi);
    }
  }

  mesh.face_to_faces.assign(nf, {});
  for(std::uint32_t fi = 0; fi < nf; ++fi) {
    const auto &tri = mesh.tris[fi];
    for(int ei = 0; ei < 3; ++ei) {
      auto v0 = tri[ei];
      auto v1 = tri[(ei + 1) % 3];
      for(auto nfi : mesh.node_to_faces[v0]) {
        if(nfi == fi) continue;
        const auto &ntri = mesh.tris[nfi];
        for(int nej = 0; nej < 3; ++nej) {
          auto nv0 = ntri[nej];
          auto nv1 = ntri[(nej + 1) % 3];
          if((nv0 == v0 && nv1 == v1) || (nv0 == v1 && nv1 == v0)) {
            auto &ff = mesh.face_to_faces[fi];
            if(std::find(ff.begin(), ff.end(), nfi) == ff.end()) {
              ff.push_back(nfi);
            }
          }
        }
      }
    }
  }
}

// Step 3: Propagate normal consistency via FF flood fill
// For adjacent faces sharing edge (v0,v1): in face A the edge is v0->v1,
// in face B it should be v1->v0 (opposite). If both have same direction,
// flip face B.
void propagate_normal_consistency(ExtractedMesh &mesh)
{
  const auto nf = mesh.tris.size();
  mesh.face_flipped.assign(nf, false);
  std::vector<bool> visited(nf, false);

  for(std::uint32_t seed = 0; seed < nf; ++seed) {
    if(visited[seed]) continue;

    std::queue<std::uint32_t> bfs;
    bfs.push(seed);
    visited[seed] = true;

    while(!bfs.empty()) {
      auto fi = bfs.front(); bfs.pop();
      const auto &tri_a = mesh.tris[fi];

      for(auto nfi : mesh.face_to_faces[fi]) {
        if(visited[nfi]) continue;

        // Find shared edge
        const auto &tri_b = mesh.tris[nfi];
        bool need_flip = false;

        for(int ea = 0; ea < 3; ++ea) {
          auto a0 = tri_a[ea];
          auto a1 = tri_a[(ea + 1) % 3];
          for(int eb = 0; eb < 3; ++eb) {
            auto b0 = tri_b[eb];
            auto b1 = tri_b[(eb + 1) % 3];
            if(a0 == b0 && a1 == b1) {
              // Same edge direction in both faces -> normals inconsistent -> flip B
              need_flip = true;
              goto found_edge;
            }
            if(a0 == b1 && a1 == b0) {
              // Opposite edge direction -> normals consistent -> no flip
              need_flip = false;
              goto found_edge;
            }
          }
        }
        found_edge:

        if(need_flip) {
          // Flip face B winding
          std::swap(mesh.tris[nfi][1], mesh.tris[nfi][2]);
          mesh.face_normals[nfi] = v3scale(mesh.face_normals[nfi], -1.0);
          mesh.face_flipped[nfi] = !mesh.face_flipped[nfi];
        }

        visited[nfi] = true;
        bfs.push(nfi);
      }
    }
  }
}

} // anonymous namespace

RegionDetectionResult detect_regions(const Domain &domain)
{
  RegionDetectionResult result;

  // Step 1: Extract mesh
  ExtractedMesh mesh;
  extract_mesh(domain, mesh);

  if(mesh.tris.empty()) {
    std::fprintf(stderr, "[RegionDetector] No triangles found in Domain\n");
    return result;
  }
  std::fprintf(stderr, "[RegionDetector] Extracted %zu nodes, %zu triangles\n",
    mesh.nodes.size(), mesh.tris.size());

  // Step 2: Compute face normals
  compute_normals(mesh);

  // Step 3: Build adjacency
  build_adjacency(mesh);

  // Step 4: Propagate normal consistency
  propagate_normal_consistency(mesh);
  std::fprintf(stderr, "[RegionDetector] Normal consistency propagated\n");

  // Step 5+6: Half-face flood fill first, then one ray test per component
  // This avoids O(N*M) ray tests — only O(components * M) instead.
  const auto nf = mesh.tris.size();
  const auto nhf = nf * 2;

  double bbox_diag = 0;
  {
    Vec3 lo = mesh.nodes[0], hi = mesh.nodes[0];
    for(const auto &n : mesh.nodes) {
      for(int k = 0; k < 3; ++k) {
        lo[k] = std::min(lo[k], n[k]);
        hi[k] = std::max(hi[k], n[k]);
      }
    }
    bbox_diag = v3len(v3sub(hi, lo));
  }
  const double eps = bbox_diag * 1e-6;

  // First: flood fill half-face connected components (ignoring inside/outside)
  // Each component is a contiguous volume on one side of the surface.
  std::vector<int> hf_component(nhf, -1);
  int next_comp = 0;

  struct CompInfo {
    std::uint32_t seed_fi;
    int seed_side;
    std::vector<std::uint32_t> face_ids; // faces in this component
  };
  std::vector<CompInfo> components;

  for(std::uint32_t hf = 0; hf < nhf; ++hf) {
    if(hf_component[hf] != -1) continue;

    int comp_id = next_comp++;
    CompInfo comp;
    comp.seed_fi = hf / 2;
    comp.seed_side = hf % 2;

    std::queue<std::uint32_t> bfs;
    bfs.push(hf);
    hf_component[hf] = comp_id;

    while(!bfs.empty()) {
      auto cur_hf = bfs.front(); bfs.pop();
      auto cur_fi = cur_hf / 2;
      auto cur_side = cur_hf % 2;
      comp.face_ids.push_back(cur_fi);

      for(auto adj_fi : mesh.face_to_faces[cur_fi]) {
        std::uint32_t adj_hf = adj_fi * 2 + cur_side;
        if(hf_component[adj_hf] != -1) continue;
        hf_component[adj_hf] = comp_id;
        bfs.push(adj_hf);
      }
    }

    // Deduplicate
    std::sort(comp.face_ids.begin(), comp.face_ids.end());
    comp.face_ids.erase(std::unique(comp.face_ids.begin(), comp.face_ids.end()), comp.face_ids.end());

    components.push_back(std::move(comp));
  }

  std::fprintf(stderr, "[RegionDetector] Found %zu half-face components\n", components.size());

  // Now: for each component, test ONE representative point to determine inside/outside
  std::vector<bool> comp_inside(components.size(), false);
  for(std::size_t ci = 0; ci < components.size(); ++ci) {
    const auto &comp = components[ci];
    const auto &t = mesh.tris[comp.seed_fi];
    Vec3 center = v3scale(v3add(v3add(mesh.nodes[t[0]], mesh.nodes[t[1]]), mesh.nodes[t[2]]), 1.0/3.0);
    double sign = (comp.seed_side == 0) ? 1.0 : -1.0;
    Vec3 test_pt = v3add(center, v3scale(mesh.face_normals[comp.seed_fi], sign * eps));

    comp_inside[ci] = point_inside(test_pt, mesh.nodes, mesh.tris);
  }

  int inside_comps = 0;
  for(auto v : comp_inside) if(v) inside_comps++;
  std::fprintf(stderr, "[RegionDetector] Inside components: %d / %zu\n",
    inside_comps, components.size());

  // Build regions from inside components
  std::vector<int> hf_region(nhf, -1);
  int next_region_id = 0;

  for(std::size_t ci = 0; ci < components.size(); ++ci) {
    if(!comp_inside[ci]) continue;

    int region_id = next_region_id++;
    MeshRegion region;
    region.id = region_id;
    region.boundary_face_ids = components[ci].face_ids;

    // Mark half-faces
    for(auto fi : region.boundary_face_ids) {
      hf_region[fi * 2 + components[ci].seed_side] = region_id;
    }

    // Interior point
    const auto &t = mesh.tris[components[ci].seed_fi];
    Vec3 center = v3scale(v3add(v3add(mesh.nodes[t[0]], mesh.nodes[t[1]]), mesh.nodes[t[2]]), 1.0/3.0);
    double sign = (components[ci].seed_side == 0) ? 1.0 : -1.0;
    region.interior_point = v3add(center, v3scale(mesh.face_normals[components[ci].seed_fi], sign * eps * 10));

    result.regions.push_back(std::move(region));
  }

  // Step 7: Build face_to_regions
  result.face_to_regions.resize(nf, {-1, -1});
  for(std::uint32_t fi = 0; fi < nf; ++fi) {
    result.face_to_regions[fi][0] = hf_region[fi * 2];     // Positive side
    result.face_to_regions[fi][1] = hf_region[fi * 2 + 1]; // Negative side
  }

  result.face_refs = mesh.face_refs;

  std::fprintf(stderr, "[RegionDetector] Detected %zu regions:\n", result.regions.size());
  for(const auto &r : result.regions) {
    std::fprintf(stderr, "  Region %d: %zu boundary faces, interior=(%.2f, %.2f, %.2f)\n",
      r.id, r.boundary_face_ids.size(),
      r.interior_point[0], r.interior_point[1], r.interior_point[2]);
  }

  return result;
}

} // namespace sqmesh::mesh
