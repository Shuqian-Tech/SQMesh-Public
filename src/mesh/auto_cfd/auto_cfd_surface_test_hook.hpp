// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#pragma once

#include "sqmesh/base/api.hpp"
#include "sqmesh/geo/api.hpp"
#include "sqmesh/mesh/api.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace sqmesh::mesh::testing {

struct AutoCfdSurfaceQualityGate final {
  double minimum_min_angle_degrees = 0.0;
  double maximum_max_angle_degrees = 180.0;
  double maximum_aspect_ratio = std::numeric_limits<double>::infinity();
  double minimum_radius_ratio = 0.0;
  double maximum_skewness = 1.0;
};

struct AutoCfdSurfaceSourceFaceFailureCounts final {
  std::size_t node_reprojection_failure_count = 0U;
  std::size_t source_face_containment_failure_count = 0U;
  std::size_t orientation_flip_count = 0U;

  [[nodiscard]] std::size_t total_failure_count() const noexcept
  {
    return node_reprojection_failure_count +
           source_face_containment_failure_count +
           orientation_flip_count;
  }

  [[nodiscard]] bool has_failures() const noexcept
  {
    return total_failure_count() > 0U;
  }
};

struct AutoCfdSurfaceMeshingInstabilityDiagnostics final {
  AutoCfdSurfaceSourceFaceFailureCounts seed_source_face_failures {};
  AutoCfdSurfaceSourceFaceFailureCounts front_candidate_source_face_failures {};
  AutoCfdSurfaceSourceFaceFailureCounts repair_candidate_source_face_failures {};
  AutoCfdSurfaceSourceFaceFailureCounts final_filter_source_face_failures {};
  std::size_t repair_flip_count = 0U;

  [[nodiscard]] std::size_t total_failure_count() const noexcept
  {
    return seed_source_face_failures.total_failure_count() +
           front_candidate_source_face_failures.total_failure_count() +
           repair_candidate_source_face_failures.total_failure_count() +
           final_filter_source_face_failures.total_failure_count();
  }

  [[nodiscard]] bool has_failures() const noexcept
  {
    return total_failure_count() > 0U;
  }
};

struct AutoCfdSurfaceFaceMeshingStats final {
  geo::TopologyEntityId face {};
  std::size_t seed_vertex_count = 0U;
  std::size_t seed_triangle_count = 0U;
  std::size_t final_vertex_count = 0U;
  std::size_t accepted_triangle_count = 0U;
  std::size_t inserted_vertex_count = 0U;
  std::size_t front_iteration_count = 0U;
  bool inserted_vertex_limit_reached = false;
  bool front_iteration_limit_reached = false;
  AutoCfdSurfaceMeshingInstabilityDiagnostics instability {};
};

struct AutoCfdSurfaceCandidateStats final {
  std::size_t face_count = 0U;
  std::size_t total_accepted_triangle_count = 0U;
  std::size_t total_inserted_vertex_count = 0U;
  std::size_t total_front_iteration_count = 0U;
  std::size_t max_inserted_vertex_count = 0U;
  std::size_t max_front_iteration_count = 0U;
  std::size_t inserted_vertex_cap_per_face = 0U;
  std::size_t front_iteration_cap_per_face = 0U;
  bool any_inserted_vertex_limit_reached = false;
  bool any_front_iteration_limit_reached = false;
  std::size_t faces_with_seed_source_face_instability = 0U;
  std::size_t faces_with_front_source_face_instability = 0U;
  std::size_t faces_with_repair_source_face_instability = 0U;
  std::size_t faces_with_final_filter_source_face_instability = 0U;
  std::size_t faces_with_any_meshing_instability = 0U;
  AutoCfdSurfaceMeshingInstabilityDiagnostics instability {};
  std::vector<AutoCfdSurfaceFaceMeshingStats> faces {};
};

enum class AutoCfdSurfaceFinalScreenFailureKind : std::uint8_t {
  none = 0,
  non_triangle_face = 1,
  duplicate_triangle = 2,
  non_manifold_edge = 3,
  degenerate_triangle = 4,
  quality_gate_rejection = 5,
  missing_face_owner = 6,
  missing_face_view = 7,
  node_reprojection_failure = 8,
  centroid_left_source_face = 9,
  orientation_flip = 10,
  no_surface_triangles = 11,
  internal_quality_guardrail_rejection = 12,
  boundary_topology_mismatch = 13,
};

struct AutoCfdSurfaceFinalScreenDiagnostics final {
  std::size_t total_triangle_count = 0;
  std::size_t non_triangle_face_count = 0;
  std::size_t duplicate_triangle_count = 0;
  std::size_t non_manifold_edge_count = 0;
  std::size_t degenerate_triangle_count = 0;
  std::size_t quality_gate_rejection_count = 0;
  std::size_t missing_face_owner_count = 0;
  std::size_t missing_face_view_count = 0;
  std::size_t node_reprojection_failure_count = 0;
  std::size_t centroid_left_source_face_count = 0;
  std::size_t orientation_flip_count = 0;
  std::size_t internal_quality_guardrail_rejection_count = 0;
  std::size_t boundary_topology_mismatch_count = 0;
  std::size_t minimum_angle_rejection_count = 0;
  std::size_t maximum_angle_rejection_count = 0;
  std::size_t aspect_ratio_rejection_count = 0;
  std::size_t radius_ratio_rejection_count = 0;
  std::size_t skewness_rejection_count = 0;
  double rejected_min_angle_minimum = std::numeric_limits<double>::infinity();
  double rejected_max_angle_maximum = 0.0;
  double rejected_aspect_ratio_maximum = 0.0;
  double rejected_radius_ratio_minimum = std::numeric_limits<double>::infinity();
  double rejected_skewness_maximum = 0.0;
  AutoCfdSurfaceFinalScreenFailureKind first_failure_kind =
    AutoCfdSurfaceFinalScreenFailureKind::none;

  [[nodiscard]] bool has_failures() const noexcept
  {
    return first_failure_kind != AutoCfdSurfaceFinalScreenFailureKind::none;
  }
};

// Internal review/test seam for the Auto CFD surface mesh quality screens.
// This exposes both the broader development guardrail and the tighter public
// delivered gate for reviewer-visible diagnostics only. This header
// intentionally lives under `src/sdk` and is not part of the installed SDK.
[[nodiscard]] AutoCfdSurfaceQualityGate
auto_cfd_surface_internal_quality_guardrail() noexcept;
[[nodiscard]] AutoCfdSurfaceQualityGate auto_cfd_surface_quality_gate() noexcept;

[[nodiscard]] base::StatusCode build_auto_cfd_surface_candidate_domain(
  geo::ModelHandle model_handle,
  const ParameterDictionary &parameters,
  Domain &domain,
  base::ContextHandle context
);

[[nodiscard]] base::StatusCode build_auto_cfd_surface_boundary_domain(
  geo::ModelHandle model_handle,
  const ParameterDictionary &parameters,
  Domain &domain,
  base::ContextHandle context
);

[[nodiscard]] base::StatusCode inspect_auto_cfd_surface_candidate_stats(
  geo::ModelHandle model_handle,
  const ParameterDictionary &parameters,
  AutoCfdSurfaceCandidateStats &stats,
  base::ContextHandle context
);

[[nodiscard]] base::StatusCode inspect_auto_cfd_surface_quality_gate(
  const Domain &domain
) noexcept;

[[nodiscard]] base::StatusCode inspect_auto_cfd_surface_final_screen(
  geo::ModelHandle model_handle,
  const ParameterDictionary &parameters,
  const Domain &domain,
  AutoCfdSurfaceFinalScreenDiagnostics &diagnostics,
  base::ContextHandle context
);

} // namespace sqmesh::mesh::testing
