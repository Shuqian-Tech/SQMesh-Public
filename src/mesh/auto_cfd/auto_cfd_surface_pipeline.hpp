// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#pragma once

#include "../sizing/mesh_size_controls.hpp"
#include "../sizing/size_function.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <vector>

namespace sqmesh::mesh::detail {

inline constexpr std::uint32_t invalid_auto_cfd_surface_part_index =
  std::numeric_limits<std::uint32_t>::max();
inline constexpr std::uint32_t invalid_auto_cfd_surface_octree_node_index =
  std::numeric_limits<std::uint32_t>::max();

struct AutoCfdSurfaceParameters final {
  double minimum_length = 0.0;
  double maximum_length = 0.0;
  double distortion_angle = 0.0;
  double growth_rate = 0.0;
  double proximity_maximum_normals_angle = 90.0;
  double proximity_length_to_gap_ratio = 1.0 / 3.0;
  double proximity_minimum_length = 0.0;
  double spacing_growth_rate = 0.0;
  double spacing_feature_angle = 0.0;
  double spacing_minimum_length = 0.0;
  double spacing_maximum_length = 0.0;
  double spacing_sharp_angle_limit = 0.0;
  double spacing_sharp_angle_length = 0.0;
  double spacing_maximum_normals_angle = 60.0;
  double spacing_length_to_gap_ratio = 1.0 / 3.0;
  double spacing_proximity_minimum_length = 0.0;
  bool spacing_feature_angle_specified = false;
  bool spacing_sharp_angle_specified = false;
  bool proximity = false;
  bool self_proximity = false;
  bool pid_proximity = false;
  bool spacing_proximity = false;
  bool spacing_self_proximity = false;
  bool spacing_pid_proximity = false;
  bool allow_quality_gate_failure = false;
  bool allow_final_screen_failure = false;
};

enum class AutoCfdSurfaceAnchorSource : std::uint8_t {
  none = 0,
  hard_point = 1,
  non_tangent_feature_intersection = 2,
};

struct AutoCfdSurfaceAnchorCandidate final {
  geo::TopologyEntityId vertex {};
  AutoCfdSurfaceAnchorSource source = AutoCfdSurfaceAnchorSource::none;
  bool frozen = false;
  bool tangent_continuous = false;
  std::vector<geo::TopologyEntityId> incident_edges {};
  std::vector<geo::TopologyEntityId> incident_feature_edges {};
};

struct AutoCfdSurfaceCurveWorkItem final {
  geo::TopologyEntityId edge {};
  std::vector<geo::TopologyEntityId> face_ids {};
  std::vector<geo::TopologyEntityId> vertex_ids {};
  geo::EdgeCurveInfo curve {};
  double target_size = 0.0;
};

struct AutoCfdSurfaceFaceWorkItem final {
  geo::TopologyEntityId face {};
  geo::FaceBoundaryLoops boundary {};
  double target_size = 0.0;
  std::uint32_t proximity_part_index = invalid_auto_cfd_surface_part_index;
};

enum class AutoCfdSurfaceSizingSourceKind : std::uint8_t {
  edge_curvature = 0,
  face_curvature = 1,
  self_proximity = 2,
  pid_proximity = 3,
};

struct AutoCfdSurfaceSizingSource final {
  geo::TopologyEntityId owner {};
  AutoCfdSurfaceSizingSourceKind kind = AutoCfdSurfaceSizingSourceKind::edge_curvature;
  geo::TopologyEntityId hit_owner {};
  geo::Point3 position {0.0, 0.0, 0.0};
  double primary_parameter = 0.0;
  double secondary_parameter = 0.0;
  double hit_distance = 0.0;
  double hit_normal_angle_degrees = 0.0;
  double raw_target_size = 0.0;
  double target_size = 0.0;
  double relaxed_target_size = 0.0;
};

/// A triangle face source for the adaptive octree size function.
/// Stores three vertex positions and their corresponding curvature-derived
/// target sizes.  Used by rebuild_auto_cfd_surface_background_field to feed
/// add_face_source() for continuous coverage.
struct AutoCfdSurfaceCurvatureFaceSource final {
  std::array<geo::Point3, 3> positions {};
  std::array<double, 3> sizes {0.0, 0.0, 0.0};
};

/// A line-segment edge source for the adaptive octree size function.
/// Emitted per consecutive sample pair along a GEdge so curve curvature
/// drives continuous coverage along the edge instead of only point seeds.
struct AutoCfdSurfaceCurvatureEdgeSource final {
  std::array<geo::Point3, 2> positions {};
  std::array<double, 2> sizes {0.0, 0.0};
};

struct AutoCfdSurfaceBackgroundBounds final {
  geo::Point3 minimum {0.0, 0.0, 0.0};
  geo::Point3 maximum {0.0, 0.0, 0.0};
  bool defined = false;
};

struct AutoCfdSurfaceBackgroundGrid final {
  geo::Point3 origin {0.0, 0.0, 0.0};
  geo::Point3 minimum {0.0, 0.0, 0.0};
  geo::Point3 maximum {0.0, 0.0, 0.0};
  std::array<double, 3> spacing {0.0, 0.0, 0.0};
  std::array<std::size_t, 3> node_counts {0U, 0U, 0U};
  std::vector<double> node_values {};
  std::vector<AutoCfdSurfaceSizingSource> seed_sources {};
  std::vector<std::vector<std::uint32_t>> cell_source_indices {};
  std::size_t seed_count = 0U;
  bool built = false;
};

struct AutoCfdSurfaceBackgroundOctreeNode final {
  geo::Point3 bounds_min {0.0, 0.0, 0.0};
  geo::Point3 bounds_max {0.0, 0.0, 0.0};
  double minimum_target_size = std::numeric_limits<double>::infinity();
  std::uint32_t first_source = 0U;
  std::uint32_t source_count = 0U;
  std::array<std::uint32_t, 8> children {
    invalid_auto_cfd_surface_octree_node_index,
    invalid_auto_cfd_surface_octree_node_index,
    invalid_auto_cfd_surface_octree_node_index,
    invalid_auto_cfd_surface_octree_node_index,
    invalid_auto_cfd_surface_octree_node_index,
    invalid_auto_cfd_surface_octree_node_index,
    invalid_auto_cfd_surface_octree_node_index,
    invalid_auto_cfd_surface_octree_node_index,
  };
  bool leaf = true;
};

struct AutoCfdSurfaceBackgroundOctree final {
  geo::Point3 minimum {0.0, 0.0, 0.0};
  geo::Point3 maximum {0.0, 0.0, 0.0};
  std::vector<AutoCfdSurfaceSizingSource> seed_sources {};
  std::vector<std::uint32_t> source_indices {};
  std::vector<AutoCfdSurfaceBackgroundOctreeNode> nodes {};
  std::size_t leaf_count = 0U;
  std::size_t maximum_depth = 0U;
  bool built = false;
};

struct AutoCfdSurfaceProximityStats final {
  std::size_t shell_count = 0U;
  std::size_t part_count = 0U;
  std::size_t proxy_triangle_count = 0U;
  std::size_t bvh_node_count = 0U;
  std::size_t bvh_leaf_count = 0U;
  std::size_t bvh_max_leaf_triangle_count = 0U;
  std::size_t face_sample_count = 0U;
  std::size_t edge_sample_count = 0U;
  std::size_t ray_hit_count = 0U;
};

struct AutoCfdSurfaceBoundaryNode final {
  geo::Point3 position {0.0, 0.0, 0.0};
  geo::TopologyEntityId topology_vertex {};
  AutoCfdSurfaceAnchorSource anchor_source = AutoCfdSurfaceAnchorSource::none;
  bool topology_endpoint = false;
  bool frozen_anchor = false;
};

struct AutoCfdSurfaceEdgePoint final {
  std::uint32_t node_index = invalid_index;
  double parameter = 0.0;
  double target_size = 0.0;
};

struct AutoCfdSurfaceEdgeDiscretization final {
  geo::TopologyEntityId edge {};
  geo::TopologyEntityId start_vertex {};
  geo::TopologyEntityId end_vertex {};
  bool short_edge_fallback = false;
  double approximate_length = 0.0;
  std::vector<AutoCfdSurfaceEdgePoint> points {};
};

struct AutoCfdSurfaceBoundaryRetryRequest final {
  geo::TopologyEntityId edge {};
  std::uint32_t first_node_index = invalid_index;
  std::uint32_t second_node_index = invalid_index;
};

struct AutoCfdSurfaceBoundaryState final {
  std::vector<AutoCfdSurfaceBoundaryNode> nodes {};
  std::vector<std::uint32_t> vertex_node_indices {};
  std::vector<AutoCfdSurfaceEdgeDiscretization> edge_discretizations {};
};

enum class AutoCfdSurfaceSeamSupportKind : std::uint8_t {
  none = 0,
 
  supported_subset = 1,
  
  periodic = 2,
  
};

enum class AutoCfdSurfaceMetricFallbackKind : std::uint8_t {
  none = 0,
  isotropic_from_jacobian_average = 1,
  isotropic_from_uv_identity = 2,
};

enum class AutoCfdSurfaceFacePreprocessDisposition : std::uint8_t {
  uv_ready = 0,
  uv_ready_with_metric_fallback = 1,
  fallback_only = 2,
};

struct AutoCfdSurfaceFaceMetricTensor final {
  geo::TopologyEntityId face {};
  double u = 0.0;
  double v = 0.0;
  double target_size = 0.0;
  std::array<double, 3> first_fundamental_form {0.0, 0.0, 0.0};
  std::array<double, 3> tensor {0.0, 0.0, 0.0};
  std::array<double, 2> raw_eigenvalues {0.0, 0.0};
  std::array<double, 2> clamped_eigenvalues {0.0, 0.0};
  double determinant = 0.0;
  bool first_fundamental_form_defined = false;
  bool usable = false;
  bool clamped = false;
  AutoCfdSurfaceMetricFallbackKind fallback_kind =
    AutoCfdSurfaceMetricFallbackKind::none;
};

struct AutoCfdSurfaceFaceLoopPoint final {
  std::uint32_t node_index = invalid_index;
  double parameter = 0.0;
  double boundary_target_size = 0.0;
  geo::Point3 position {0.0, 0.0, 0.0};
  std::array<double, 2> uv {0.0, 0.0};
  bool uv_defined = false;
  AutoCfdSurfaceFaceMetricTensor metric {};
};

struct AutoCfdSurfaceFaceLoopSegment final {
  geo::FaceBoundaryEdgeUse edge_use {};
  std::vector<AutoCfdSurfaceFaceLoopPoint> points {};
};

struct AutoCfdSurfaceFaceLoopState final {
  geo::FaceBoundaryLoopKind kind = geo::FaceBoundaryLoopKind::unknown;
  bool closed = false;
  bool continuous = false;
  std::size_t seam_edge_use_count = 0U;
  std::size_t degenerate_edge_use_count = 0U;
  std::vector<AutoCfdSurfaceFaceLoopPoint> points {};
  std::vector<AutoCfdSurfaceFaceLoopSegment> segments {};
};


struct AutoCfdSurfaceFacePreprocessState final {
  geo::TopologyEntityId face {};
  geo::FaceBoundaryLoops boundary {};
  geo::FaceUvBounds uv_bounds {};
  std::array<double, 2> recovered_uv_min {0.0, 0.0};
  std::array<double, 2> recovered_uv_max {0.0, 0.0};
  std::vector<AutoCfdSurfaceFaceLoopState> loops {};
  AutoCfdSurfaceFaceMetricTensor reference_metric {};
  AutoCfdSurfaceSeamSupportKind seam_support =
    AutoCfdSurfaceSeamSupportKind::none;
  AutoCfdSurfaceFacePreprocessDisposition disposition =
    AutoCfdSurfaceFacePreprocessDisposition::fallback_only;
  bool uv_bounds_defined = false;
  bool uv_reconstruction_available = false;
  bool seam_unwrap_applied = false;
  bool metric_fallback_used = false;
};

struct AutoCfdSurfaceSizingFieldState final {
  double minimum_length = 0.0;
  double maximum_length = 0.0;
  double distortion_angle = 0.0;
  double growth_rate = 0.0;
  bool proximity_enabled = false;
  bool self_proximity_enabled = false;
  bool pid_proximity_enabled = false;
  std::size_t proximity_layer_count = 0U;
  double proximity_normal_angle_threshold_degrees = 0.0;
  double proximity_length_to_gap_ratio = 0.0;
  double proximity_minimum_length = 0.0;
  std::vector<double> face_target_sizes {};
  std::vector<double> edge_target_sizes {};
  std::vector<std::uint32_t> face_proximity_shell_ids {};
  std::vector<std::uint32_t> face_proximity_part_ids {};
  AutoCfdSurfaceProximityStats proximity_stats {};
  std::vector<AutoCfdSurfaceSizingSource> curvature_sources {};
  std::vector<AutoCfdSurfaceCurvatureFaceSource> curvature_face_sources {};
  std::vector<AutoCfdSurfaceCurvatureEdgeSource> curvature_edge_sources {};
  std::vector<AutoCfdSurfaceSizingSource> proximity_sources {};
  AutoCfdSurfaceBackgroundBounds background_bounds {};
  AutoCfdSurfaceBackgroundGrid background_grid {};
  AutoCfdSurfaceBackgroundOctree background_octree {};
  SizeFunction size_function {};
};

struct AutoCfdSurfacePipelineState final {
  std::uint64_t topology_revision = 0;
  AutoCfdSurfaceParameters parameters {};
  AutoCfdSurfaceSizingFieldState sizing_field {};
  AutoCfdSurfaceBoundaryState boundary {};
  std::vector<AutoCfdSurfaceAnchorCandidate> anchor_candidates {};
  std::vector<AutoCfdSurfaceCurveWorkItem> curve_work_items {};
  std::vector<AutoCfdSurfaceFaceWorkItem> face_work_items {};
  std::vector<AutoCfdSurfaceFacePreprocessState> face_preprocess_states {};
};

[[nodiscard]] base::StatusCode resolve_auto_cfd_surface_parameters(
  const ParameterDictionary &parameters,
  AutoCfdSurfaceParameters &resolved
);

[[nodiscard]] base::StatusCode build_auto_cfd_surface_pipeline_state(
  const geo::ModelView &model_view,
  const ResolvedMeshSizeControls &size_controls,
  const AutoCfdSurfaceParameters &parameters,
  AutoCfdSurfacePipelineState &state
);

[[nodiscard]] base::StatusCode rebuild_auto_cfd_surface_face_preprocess_states(
  const geo::ModelView &model_view,
  AutoCfdSurfacePipelineState &state
);

[[nodiscard]] base::StatusCode split_auto_cfd_surface_boundary_edges_for_retry(
  const geo::ModelView &model_view,
  const std::vector<AutoCfdSurfaceBoundaryRetryRequest> &requests,
  AutoCfdSurfacePipelineState &state
);

[[nodiscard]] const AutoCfdSurfaceFacePreprocessState *
find_auto_cfd_surface_face_preprocess_state(
  const AutoCfdSurfacePipelineState &state,
  geo::TopologyEntityId face
) noexcept;

[[nodiscard]] double clamp_auto_cfd_surface_size(
  const AutoCfdSurfaceSizingFieldState &field,
  double size
) noexcept;

[[nodiscard]] double effective_auto_cfd_surface_owner_target_size(
  const AutoCfdSurfaceSizingFieldState &field,
  geo::TopologyEntityId owner
) noexcept;

[[nodiscard]] base::StatusCode rebuild_auto_cfd_surface_background_field(
  AutoCfdSurfaceSizingFieldState &field
);

[[nodiscard]] double query_auto_cfd_surface_sizing_field(
  const AutoCfdSurfaceSizingFieldState &field,
  geo::TopologyEntityId owner,
  const geo::Point3 &position
) noexcept;

[[nodiscard]] double query_auto_cfd_surface_sizing_field(
  const AutoCfdSurfaceSizingFieldState &field,
  const geo::FaceView &face_view,
  double u,
  double v,
  const geo::Point3 &position
) noexcept;

[[nodiscard]] base::StatusCode sample_auto_cfd_surface_face_metric(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceSizingFieldState &field,
  double u,
  double v,
  AutoCfdSurfaceFaceMetricTensor &metric
);

} // namespace sqmesh::mesh::detail
