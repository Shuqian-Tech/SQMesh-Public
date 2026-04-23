// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#pragma once

#include "model/geometry_proxy_mesh.hpp"
#include "sqmesh/geo/api.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace sqmesh::mesh::detail {

/// Discrete geometry computed from a proxy mesh (triangle soup).
/// Provides triangle normals, areas, and chord curvature — all from raw
/// vertex/triangle arrays without any OCC/parametric surface dependency.
struct ProxyMeshGeometry final {
  // ---- input references (NOT owned) ----
  const std::vector<geo::Point3> *nodes = nullptr;
  const std::vector<std::array<std::uint32_t, 3>> *triangles = nullptr;
  /// Per-CAD-face triangle ranges. When set, compute_chord_curvature
  /// restricts edge-adjacency queries to within the range so kmax never
  /// leaks across CAD face boundaries.
  const std::vector<model::detail::ProxyTriangleRange> *face_triangle_ranges
    = nullptr;

  // ---- computed results ----
  std::vector<geo::Vector3> triangle_normals {};
  std::vector<double> triangle_areas {};
  std::vector<geo::Vector3> node_normals {};
  std::vector<double> triangle_curvatures {};  // κ_max per triangle

  // ---- methods ----

  /// Compute triangle normals, areas, and area-weighted node normals.
  void compute_normals_and_areas();

  /// Per-triangle chord curvature. For each of the triangle's 3 edges:
  ///   dist = |(centroid_f  - v1) ⊥ edge_tangent|
  ///        + |(centroid_fj - v1) ⊥ edge_tangent|
  ///   κ_edge = 2·sin(θ/2) / dist
  /// tri_kmax = MAX over valid edges (dihedral θ ∈ (0, feature_angle]).
  ///
  /// When face_triangle_ranges is set, adjacency queries are restricted to
  /// within the same CAD face so kmax never leaks across CAD face boundaries.
  ///
  /// @param feature_angle_degrees  Edges with dihedral > this are treated as
  ///        feature creases and skipped.
  void compute_chord_curvature(double feature_angle_degrees);
};

} // namespace sqmesh::mesh::detail
