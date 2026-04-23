// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "proxy_mesh_geometry.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace sqmesh::mesh::detail {
namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kEpsilon = 1.0e-30;

// Up to 2 triangle indices sharing an edge (manifold mesh).
struct EdgeNeighbors final {
  std::uint32_t tri[2] = {UINT32_MAX, UINT32_MAX};
  std::uint8_t count = 0;
};

[[nodiscard]] geo::Vector3 vec_sub(
  const geo::Point3 &a, const geo::Point3 &b
) noexcept
{
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

[[nodiscard]] geo::Vector3 vec_cross(
  const geo::Vector3 &a, const geo::Vector3 &b
) noexcept
{
  return {
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0],
  };
}

[[nodiscard]] double vec_dot(
  const geo::Vector3 &a, const geo::Vector3 &b
) noexcept
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

[[nodiscard]] double vec_norm(const geo::Vector3 &v) noexcept
{
  return std::sqrt(vec_dot(v, v));
}

[[nodiscard]] geo::Vector3 vec_normalized(const geo::Vector3 &v) noexcept
{
  const double len = vec_norm(v);
  if(len < kEpsilon) {
    return {0.0, 0.0, 0.0};
  }
  const double inv = 1.0 / len;
  return {v[0] * inv, v[1] * inv, v[2] * inv};
}

/// Pack two vertex indices into a single uint64 key (order-independent).
[[nodiscard]] std::uint64_t edge_key(std::uint32_t v0, std::uint32_t v1) noexcept
{
  if(v0 > v1) {
    std::swap(v0, v1);
  }
  return (static_cast<std::uint64_t>(v0) << 32U) | v1;
}

} // anonymous namespace

// =============================================================================
// compute_normals_and_areas
// =============================================================================

void ProxyMeshGeometry::compute_normals_and_areas()
{
  assert(nodes != nullptr && triangles != nullptr);
  const auto &N = *nodes;
  const auto &T = *triangles;
  const auto nv = N.size();
  const auto nt = T.size();

  triangle_normals.resize(nt);
  triangle_areas.resize(nt);
  node_normals.assign(nv, {0.0, 0.0, 0.0});

  for(std::size_t fi = 0; fi < nt; ++fi) {
    const auto &tri = T[fi];
    const auto e01 = vec_sub(N[tri[1]], N[tri[0]]);
    const auto e02 = vec_sub(N[tri[2]], N[tri[0]]);
    const auto cross_vec = vec_cross(e01, e02);
    const double cross_len = vec_norm(cross_vec);

    triangle_areas[fi] = 0.5 * cross_len;

    if(cross_len > kEpsilon) {
      const double inv = 1.0 / cross_len;
      triangle_normals[fi] = {cross_vec[0] * inv, cross_vec[1] * inv, cross_vec[2] * inv};
    } else {
      triangle_normals[fi] = {0.0, 0.0, 0.0};
    }

    // Accumulate area-weighted normal to vertices.
    for(int k = 0; k < 3; ++k) {
      auto &nn = node_normals[tri[k]];
      nn[0] += cross_vec[0];  // magnitude = 2*area, so area-weighted
      nn[1] += cross_vec[1];
      nn[2] += cross_vec[2];
    }
  }

  // Normalize node normals.
  for(std::size_t i = 0; i < nv; ++i) {
    node_normals[i] = vec_normalized(node_normals[i]);
  }
}

// =============================================================================
// compute_chord_curvature
// =============================================================================
//
// Chord-based two-centroid curvature: two triangles sharing an edge approximate
// two points on a circle — arc angle = dihedral θ, chord length = sum of
// centroid-to-edge perpendicular distances. κ = 2·sin(θ/2)/chord.

void ProxyMeshGeometry::compute_chord_curvature(double feature_angle_degrees)
{
  assert(nodes != nullptr && triangles != nullptr);
  const auto &N = *nodes;
  const auto &T = *triangles;
  const auto nt = T.size();

  const double threshold_angle = feature_angle_degrees * kPi / 180.0;

  triangle_curvatures.assign(nt, 0.0);
  std::vector<double> &tri_kmax = triangle_curvatures;

  // Per-face adjacency ranges (across_bound=0 semantics).
  std::vector<model::detail::ProxyTriangleRange> default_ranges;
  const std::vector<model::detail::ProxyTriangleRange> *ranges_ptr =
    face_triangle_ranges;
  if(ranges_ptr == nullptr || ranges_ptr->empty()) {
    default_ranges.push_back({0U, nt});
    ranges_ptr = &default_ranges;
  }

  // Precompute centroids.
  std::vector<geo::Point3> centroids(nt, {0.0, 0.0, 0.0});
  for(std::size_t fi = 0; fi < nt; ++fi) {
    const auto &tri = T[fi];
    const auto &p0 = N[tri[0]];
    const auto &p1 = N[tri[1]];
    const auto &p2 = N[tri[2]];
    centroids[fi] = {
      (p0[0] + p1[0] + p2[0]) / 3.0,
      (p0[1] + p1[1] + p2[1]) / 3.0,
      (p0[2] + p1[2] + p2[2]) / 3.0,
    };
  }

  for(const auto &range : *ranges_ptr) {
    if(range.count == 0U) {
      continue;
    }
    const std::size_t range_end = range.offset + range.count;

    std::unordered_map<std::uint64_t, EdgeNeighbors> local_adj;
    local_adj.reserve(range.count * 3U);
    for(std::size_t fi = range.offset; fi < range_end; ++fi) {
      const auto &tri = T[fi];
      for(int e = 0; e < 3; ++e) {
        const auto key = edge_key(tri[e], tri[(e + 1) % 3]);
        auto &en = local_adj[key];
        if(en.count < 2) {
          en.tri[en.count] = static_cast<std::uint32_t>(fi);
          ++en.count;
        }
      }
    }

    for(std::size_t fi = range.offset; fi < range_end; ++fi) {
      const auto &tri = T[fi];
      const auto &n_f = triangle_normals[fi];
      if(triangle_areas[fi] < kEpsilon) {
        continue;
      }
      const auto &Cf = centroids[fi];

      double max_kappa = 0.0;

      for(int e = 0; e < 3; ++e) {
        const auto v0 = tri[e];
        const auto v1 = tri[(e + 1) % 3];
        const auto key = edge_key(v0, v1);

        const auto adj_it = local_adj.find(key);
        if(adj_it == local_adj.end() || adj_it->second.count < 2) {
          continue;  // mesh or CAD-face boundary
        }

        const auto &en = adj_it->second;
        const std::uint32_t fj = (en.tri[0] == fi) ? en.tri[1] : en.tri[0];
        const auto &n_fj = triangle_normals[fj];

        const double cos_angle =
          std::clamp(vec_dot(n_f, n_fj), -1.0, 1.0);
        const double angle = std::acos(cos_angle);
        if(angle <= kEpsilon || angle > threshold_angle) {
          continue;  // near-flat or feature edge
        }

        // Edge tangent, normalized.
        const auto edge_vec = vec_sub(N[v1], N[v0]);
        const double elen = vec_norm(edge_vec);
        if(elen < kEpsilon) {
          continue;
        }
        const geo::Vector3 etan = {
          edge_vec[0] / elen, edge_vec[1] / elen, edge_vec[2] / elen
        };

        // Perpendicular distance from v0 to Cf (projected out of etan).
        const auto df_vec = vec_sub(Cf, N[v0]);
        const double df_dot = vec_dot(df_vec, etan);
        const double df_len_sq = vec_dot(df_vec, df_vec);
        const double df_perp_sq = df_len_sq - df_dot * df_dot;
        const double dist_f = (df_perp_sq > 0.0) ? std::sqrt(df_perp_sq) : 0.0;

        // Perpendicular distance from v0 to Cfj.
        const auto dj_vec = vec_sub(centroids[fj], N[v0]);
        const double dj_dot = vec_dot(dj_vec, etan);
        const double dj_len_sq = vec_dot(dj_vec, dj_vec);
        const double dj_perp_sq = dj_len_sq - dj_dot * dj_dot;
        const double dist_fj = (dj_perp_sq > 0.0) ? std::sqrt(dj_perp_sq) : 0.0;

        const double dist = dist_f + dist_fj;
        if(dist <= kEpsilon) {
          continue;
        }

        const double kappa = 2.0 * std::sin(angle * 0.5) / dist;
        if(kappa > max_kappa) {
          max_kappa = kappa;
        }
      }

      tri_kmax[fi] = max_kappa;
    }
  }
}

} // namespace sqmesh::mesh::detail
