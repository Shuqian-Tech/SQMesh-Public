// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#pragma once

#include "sqmesh/base/api.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string_view>
#include <vector>

namespace sqmesh::geo {

using ModelHandle = sqmesh::Handle;

inline constexpr std::uint32_t invalid_topology_index =
  std::numeric_limits<std::uint32_t>::max();
inline constexpr std::size_t invalid_boundary_loop_index =
  std::numeric_limits<std::size_t>::max();

enum class IgesWriteMode : std::uint8_t {
  faces = 0,
  brep = 1,
};

enum class TopologyDimension : std::uint8_t {
  vertex = 0,
  edge = 1,
  face = 2,
  region = 3,
};

struct TopologyEntityId final {
  TopologyDimension dimension = TopologyDimension::vertex;
  std::uint32_t index = invalid_topology_index;
};

// `TopologyEntityId` values are stable only within one stored model topology
// revision. Operations that rebuild the underlying topology, such as an OCC
// `topo()` repair that modifies the shape, can advance that revision and
// require callers to reacquire ids/views from a fresh snapshot or model view.

[[nodiscard]] constexpr bool operator==(
  TopologyEntityId lhs,
  TopologyEntityId rhs
) noexcept
{
  return lhs.dimension == rhs.dimension && lhs.index == rhs.index;
}

[[nodiscard]] constexpr bool operator!=(
  TopologyEntityId lhs,
  TopologyEntityId rhs
) noexcept
{
  return !(lhs == rhs);
}

[[nodiscard]] constexpr bool is_valid(TopologyEntityId entity) noexcept
{
  return entity.index != invalid_topology_index;
}

struct TopologyEntityInfo final {
  TopologyEntityId entity {};
  std::size_t parent_count = 0;
  std::size_t child_count = 0;
};

struct TopologySnapshot final {
  std::uint64_t topology_revision = 0;
  std::vector<TopologyEntityInfo> regions {};
  std::vector<TopologyEntityInfo> faces {};
  std::vector<TopologyEntityInfo> edges {};
  std::vector<TopologyEntityInfo> vertices {};

  [[nodiscard]] const std::vector<TopologyEntityInfo> &entities(
    TopologyDimension dimension
  ) const noexcept
  {
    switch(dimension) {
    case TopologyDimension::region:
      return regions;
    case TopologyDimension::face:
      return faces;
    case TopologyDimension::edge:
      return edges;
    case TopologyDimension::vertex:
    default:
      return vertices;
    }
  }

  [[nodiscard]] std::size_t entity_count(TopologyDimension dimension) const noexcept
  {
    return entities(dimension).size();
  }
};

using Point3 = std::array<double, 3>;
using Vector3 = std::array<double, 3>;

struct FaceUvBounds final {
  TopologyEntityId face {};
  double u_min = 0.0;
  double u_max = 0.0;
  double v_min = 0.0;
  double v_max = 0.0;
};

struct FaceSample final {
  TopologyEntityId face {};
  double u = 0.0;
  double v = 0.0;
  Point3 position {0.0, 0.0, 0.0};
  Vector3 normal {0.0, 0.0, 0.0};
  bool normal_defined = false;
};

struct FaceCurvatureSample final {
  TopologyEntityId face {};
  double u = 0.0;
  double v = 0.0;
  double min_curvature = 0.0;
  double max_curvature = 0.0;
  double mean_curvature = 0.0;
  double gaussian_curvature = 0.0;
  bool curvature_defined = false;
};

struct FaceDerivatives final {
  TopologyEntityId face {};
  double u = 0.0;
  double v = 0.0;
  Point3 position {0.0, 0.0, 0.0};
  Vector3 du {0.0, 0.0, 0.0};
  Vector3 dv {0.0, 0.0, 0.0};
  Vector3 duu {0.0, 0.0, 0.0};
  Vector3 duv {0.0, 0.0, 0.0};
  Vector3 dvv {0.0, 0.0, 0.0};
  Vector3 normal {0.0, 0.0, 0.0};
  bool first_derivatives_defined = false;
  bool second_derivatives_defined = false;
  bool normal_defined = false;
};

struct FaceProjection final {
  TopologyEntityId face {};
  Point3 input_point {0.0, 0.0, 0.0};
  Point3 projected_point {0.0, 0.0, 0.0};
  double u = 0.0;
  double v = 0.0;
  double distance = 0.0;
  Vector3 normal {0.0, 0.0, 0.0};
  bool normal_defined = false;
};

struct FaceUvMapping final {
  TopologyEntityId face {};
  Point3 input_point {0.0, 0.0, 0.0};
  Point3 mapped_point {0.0, 0.0, 0.0};
  double u = 0.0;
  double v = 0.0;
  double distance = 0.0;
};

struct EdgeCurveInfo final {
  TopologyEntityId edge {};
  double parameter_min = 0.0;
  double parameter_max = 0.0;
  Point3 start_point {0.0, 0.0, 0.0};
  Point3 end_point {0.0, 0.0, 0.0};
  double approximate_length = 0.0;
};

struct EdgeTangentSample final {
  TopologyEntityId edge {};
  double parameter = 0.0;
  Point3 position {0.0, 0.0, 0.0};
  Vector3 derivative {0.0, 0.0, 0.0};
  Vector3 tangent {0.0, 0.0, 0.0};
  double speed = 0.0;
  bool tangent_defined = false;
};

struct EdgeCurveSamplingOptions final {
  // When positive, the edge approximate length is used to choose at least
  // `ceil(approximate_length / target_segment_length)` parameter-space
  // segments.
  double target_segment_length = 0.0;
  // The helper always returns `min_segment_count + 1` or more endpoint-inclusive
  // samples.
  std::size_t min_segment_count = 1;
};

struct EdgeCurveSamples final {
  EdgeCurveInfo curve {};
  std::vector<EdgeTangentSample> samples {};
};

enum class FaceBoundaryLoopKind : std::uint8_t {
  unknown = 0,
  outer = 1,
  inner = 2,
};

struct FaceBoundaryEdgeUse final {
  TopologyEntityId edge {};
  TopologyEntityId start_vertex {};
  TopologyEntityId end_vertex {};
  bool same_orientation_as_edge = false;
  // Trim-aware backends can mark face-local seam or degenerate edge uses
  // directly. Discrete backends keep these flags false.
  bool is_seam = false;
  bool is_degenerate = false;
};

struct FaceBoundaryLoop final {
  FaceBoundaryLoopKind kind = FaceBoundaryLoopKind::unknown;
  bool closed = false;
  bool continuous = false;
  std::size_t seam_edge_use_count = 0;
  std::size_t degenerate_edge_use_count = 0;
  std::size_t repeated_edge_use_count = 0;
  std::vector<FaceBoundaryEdgeUse> edge_uses {};
  std::vector<TopologyEntityId> vertex_ids {};
};

struct FaceBoundaryLoops final {
  TopologyEntityId face {};
  std::vector<FaceBoundaryLoop> loops {};
  std::size_t primary_outer_loop_index = invalid_boundary_loop_index;
  std::size_t outer_loop_count = 0;
  std::size_t inner_loop_count = 0;
  std::size_t unknown_loop_count = 0;
  std::size_t closed_loop_count = 0;
  std::size_t open_loop_count = 0;
  std::size_t non_continuous_loop_count = 0;
  std::size_t seam_loop_count = 0;
  std::size_t seam_edge_use_count = 0;
  std::size_t degenerate_loop_count = 0;
  std::size_t degenerate_edge_use_count = 0;
  std::size_t repeated_edge_loop_count = 0;
  std::size_t repeated_edge_use_count = 0;
  bool has_holes = false;
  bool has_seams = false;
};

// Stable high-level geometry access is a lightweight OCC-free value/view graph
// over `TopologyEntityId`, not a backend-owned `GEntity` clone.
struct VertexView final {
  ModelHandle model_handle = sqmesh::invalid_handle;
  base::ContextHandle context_handle = sqmesh::invalid_handle;
  TopologyEntityId entity {};
  std::vector<TopologyEntityId> edge_ids {};
};

struct EdgeView final {
  ModelHandle model_handle = sqmesh::invalid_handle;
  base::ContextHandle context_handle = sqmesh::invalid_handle;
  TopologyEntityId entity {};
  std::vector<TopologyEntityId> face_ids {};
  std::vector<TopologyEntityId> vertex_ids {};
};

struct FaceView final {
  ModelHandle model_handle = sqmesh::invalid_handle;
  base::ContextHandle context_handle = sqmesh::invalid_handle;
  TopologyEntityId entity {};
  std::vector<TopologyEntityId> region_ids {};
  std::vector<TopologyEntityId> edge_ids {};
  FaceBoundaryLoops ordered_boundary {};
};

struct RegionView final {
  ModelHandle model_handle = sqmesh::invalid_handle;
  base::ContextHandle context_handle = sqmesh::invalid_handle;
  TopologyEntityId entity {};
  std::vector<TopologyEntityId> face_ids {};
};

struct ModelView final {
  ModelHandle model_handle = sqmesh::invalid_handle;
  base::ContextHandle context_handle = sqmesh::invalid_handle;
  TopologySnapshot snapshot {};
  std::vector<RegionView> regions {};
  std::vector<FaceView> faces {};
  std::vector<EdgeView> edges {};
  std::vector<VertexView> vertices {};
  std::vector<TopologyEntityId> root_face_ids {};
  std::vector<TopologyEntityId> root_edge_ids {};
  std::vector<TopologyEntityId> root_vertex_ids {};

  [[nodiscard]] std::size_t entity_count(TopologyDimension dimension) const noexcept
  {
    return snapshot.entity_count(dimension);
  }

  [[nodiscard]] const RegionView *find_region(TopologyEntityId entity) const noexcept
  {
    if(entity.dimension != TopologyDimension::region || entity.index >= regions.size()) {
      return nullptr;
    }

    return &regions[entity.index];
  }

  [[nodiscard]] const FaceView *find_face(TopologyEntityId entity) const noexcept
  {
    if(entity.dimension != TopologyDimension::face || entity.index >= faces.size()) {
      return nullptr;
    }

    return &faces[entity.index];
  }

  [[nodiscard]] const EdgeView *find_edge(TopologyEntityId entity) const noexcept
  {
    if(entity.dimension != TopologyDimension::edge || entity.index >= edges.size()) {
      return nullptr;
    }

    return &edges[entity.index];
  }

  [[nodiscard]] const VertexView *find_vertex(TopologyEntityId entity) const noexcept
  {
    if(entity.dimension != TopologyDimension::vertex || entity.index >= vertices.size()) {
      return nullptr;
    }

    return &vertices[entity.index];
  }
};

[[nodiscard]] constexpr const FaceBoundaryLoops &ordered_boundary_loops(
  const FaceView &face_view
) noexcept
{
  return face_view.ordered_boundary;
}

[[nodiscard]] inline const FaceBoundaryLoop *primary_outer_boundary_loop(
  const FaceBoundaryLoops &boundary
) noexcept
{
  if(boundary.primary_outer_loop_index >= boundary.loops.size()) {
    return nullptr;
  }

  return &boundary.loops[boundary.primary_outer_loop_index];
}

struct FeatureEdgeOptions final {
  double feature_angle_degrees = 30.0;
  bool include_boundary_edges = true;
  bool include_non_manifold_edges = true;
};

struct FeatureEdgeReport final {
  std::vector<TopologyEntityId> edges {};
  std::size_t sharp_edge_count = 0;
  std::size_t boundary_edge_count = 0;
  std::size_t non_manifold_edge_count = 0;
};

struct ModelSummary final {
  std::size_t entity_count = 0;
  std::size_t compound_count = 0;
  std::size_t compsolid_count = 0;
  std::size_t solid_count = 0;
  std::size_t shell_count = 0;
  std::size_t face_count = 0;
  std::size_t wire_count = 0;
  std::size_t edge_count = 0;
  std::size_t vertex_count = 0;
};

struct StepImportOptions final {
  double system_length_unit = 0.0;
  // Healing defaults: all healing ON during import so that the geometry
  // is clean before proxy mesh construction and surface meshing.
  double tolerance = 1.0e-8;
  double min_tolerance = 0.0;
  double max_tolerance = 0.0;
  bool fix_degenerated = true;
  bool fix_small_edges = true;
  bool fix_small_faces = true;
  bool sew_faces = true;
  bool make_solids = true;
};

struct IgesImportOptions final {
  bool read_visible_only = false;
};

struct StlImportOptions final {
  // Coincident STL vertices are welded with a tolerance scaled by the model
  // bounding-box diagonal.
  double relative_merge_tolerance = 1.0e-9;
};

struct StepExportOptions final {
  double linear_tolerance = 0.0;
  std::string_view unit_name {};
  std::string_view schema {};
};

struct IgesExportOptions final {
  std::string_view unit_name {};
  IgesWriteMode write_mode = IgesWriteMode::brep;
};

struct TopologyCheckReport final {
  bool is_valid = false;
  std::size_t free_edge_count = 0;
  std::size_t contiguous_edge_count = 0;
  std::size_t multiple_edge_count = 0;
  ModelSummary model_summary {};
};

struct TopoOptions final {
  double tolerance = 1.0e-6;
  double min_tolerance = 0.0;
  double max_tolerance = 0.0;
  bool fix_degenerated = true;
  bool fix_small_edges = true;
  bool fix_small_faces = true;
  bool sew_faces = true;
  bool make_solids = true;
};

struct TopoReport final {
  TopologyCheckReport before {};
  TopologyCheckReport after {};
  std::uint64_t topology_revision_before = 0;
  std::uint64_t topology_revision_after = 0;
  bool modified = false;
  bool topology_identity_changed = false;
  bool free_edges_reduced = false;
  bool sewing_performed = false;
  bool sewing_modified = false;
  bool shape_fix_performed = false;
  bool shape_fix_modified = false;
};

[[nodiscard]] std::string_view module_name() noexcept;
[[nodiscard]] bool cad_io_available() noexcept;
[[nodiscard]] base::StatusCode create_placeholder_model(
  ModelHandle &model_handle,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode import_step(
  std::string_view path,
  ModelHandle &model_handle,
  const StepImportOptions &options = {},
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode import_iges(
  std::string_view path,
  ModelHandle &model_handle,
  const IgesImportOptions &options = {},
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
// STL import uses a discrete surface backend:
// - supported: topology_snapshot/topology_children/topology_parents,
//   edge_curve_info, sample_edge_tangent, face_boundary_loops, feature_edges
// - unsupported: parametric face queries such as face_uv_bounds,
//   sample_face, sample_face_curvature, sample_face_derivatives,
//   project_point_to_face, recover_face_uv, and OCC-only STEP/IGES
//   export/topology-repair operations
//
// Discrete face entities represent connected triangle patches. Their boundary
// loops expose ordered open-boundary chains/loops when present, and loop kind
// remains `unknown` because STL does not carry trim-wire classification. The
// richer outer/inner loop counts therefore remain zero unless the backend can
// classify trim wires explicitly. Seam and degenerate trim metadata likewise
// remain false/zero on discrete-backed boundaries.
[[nodiscard]] base::StatusCode import_stl(
  std::string_view path,
  ModelHandle &model_handle,
  const StlImportOptions &options = {},
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode export_step(
  ModelHandle model_handle,
  std::string_view path,
  const StepExportOptions &options = {},
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode export_iges(
  ModelHandle model_handle,
  std::string_view path,
  const IgesExportOptions &options = {},
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode model_summary(
  ModelHandle model_handle,
  ModelSummary &summary,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
// Returns the shared proxy-backed mesh handle when the model storage exposes a
// coarse proxy mesh. `mesh_summary()` / `Domain::summary()` on that handle
// count only `EntityGroupRole::computational` entities, so a proxy-only handle
// reports zero node/edge/face/cell counts until meshing algorithms append
// computational entity_groups. Use `domain_snapshot()` and entity_group-role inspection to
// inspect retained `geometric_proxy` content.
[[nodiscard]] base::StatusCode model_proxy_mesh(
  ModelHandle model_handle,
  sqmesh::Handle &mesh_handle,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode topology_snapshot(
  ModelHandle model_handle,
  TopologySnapshot &snapshot,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode model_view(
  ModelHandle model_handle,
  ModelView &view,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode topology_children(
  ModelHandle model_handle,
  TopologyEntityId entity,
  std::vector<TopologyEntityId> &children,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode topology_parents(
  ModelHandle model_handle,
  TopologyEntityId entity,
  std::vector<TopologyEntityId> &parents,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode face_uv_bounds(
  ModelHandle model_handle,
  TopologyEntityId face_entity,
  FaceUvBounds &bounds,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode face_uv_bounds(
  const FaceView &face_view,
  FaceUvBounds &bounds
) noexcept;
[[nodiscard]] base::StatusCode sample_face(
  ModelHandle model_handle,
  TopologyEntityId face_entity,
  double u,
  double v,
  FaceSample &sample,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode sample_face(
  const FaceView &face_view,
  double u,
  double v,
  FaceSample &sample
) noexcept;
[[nodiscard]] base::StatusCode sample_face_curvature(
  ModelHandle model_handle,
  TopologyEntityId face_entity,
  double u,
  double v,
  FaceCurvatureSample &sample,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode sample_face_curvature(
  const FaceView &face_view,
  double u,
  double v,
  FaceCurvatureSample &sample
) noexcept;
[[nodiscard]] base::StatusCode sample_face_derivatives(
  ModelHandle model_handle,
  TopologyEntityId face_entity,
  double u,
  double v,
  FaceDerivatives &sample,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode sample_face_derivatives(
  const FaceView &face_view,
  double u,
  double v,
  FaceDerivatives &sample
) noexcept;
[[nodiscard]] base::StatusCode project_point_to_face(
  ModelHandle model_handle,
  TopologyEntityId face_entity,
  const Point3 &point,
  FaceProjection &projection,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode project_point_to_face(
  const FaceView &face_view,
  const Point3 &point,
  FaceProjection &projection
) noexcept;
[[nodiscard]] base::StatusCode recover_face_uv(
  ModelHandle model_handle,
  TopologyEntityId face_entity,
  const Point3 &point,
  FaceUvMapping &mapping,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode recover_face_uv(
  const FaceView &face_view,
  const Point3 &point,
  FaceUvMapping &mapping
) noexcept;
[[nodiscard]] base::StatusCode recover_face_uv_from_edge(
  ModelHandle model_handle,
  TopologyEntityId face_entity,
  TopologyEntityId edge_entity,
  double edge_parameter,
  FaceUvMapping &mapping,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode recover_face_uv_from_edge(
  const FaceView &face_view,
  TopologyEntityId edge_entity,
  double edge_parameter,
  FaceUvMapping &mapping
) noexcept;
[[nodiscard]] base::StatusCode recover_face_uv_from_edge_use(
  ModelHandle model_handle,
  TopologyEntityId face_entity,
  const FaceBoundaryEdgeUse &edge_use,
  double edge_parameter,
  FaceUvMapping &mapping,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode recover_face_uv_from_edge_use(
  const FaceView &face_view,
  const FaceBoundaryEdgeUse &edge_use,
  double edge_parameter,
  FaceUvMapping &mapping
) noexcept;
[[nodiscard]] base::StatusCode edge_curve_info(
  ModelHandle model_handle,
  TopologyEntityId edge_entity,
  EdgeCurveInfo &info,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode edge_curve_info(
  const EdgeView &edge_view,
  EdgeCurveInfo &info
) noexcept;
[[nodiscard]] base::StatusCode sample_edge_tangent(
  ModelHandle model_handle,
  TopologyEntityId edge_entity,
  double parameter,
  EdgeTangentSample &sample,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode sample_edge_tangent(
  const EdgeView &edge_view,
  double parameter,
  EdgeTangentSample &sample
) noexcept;
[[nodiscard]] base::StatusCode sample_edge_curve(
  ModelHandle model_handle,
  TopologyEntityId edge_entity,
  const EdgeCurveSamplingOptions &options,
  EdgeCurveSamples &samples,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode sample_edge_curve(
  const EdgeView &edge_view,
  const EdgeCurveSamplingOptions &options,
  EdgeCurveSamples &samples
) noexcept;
[[nodiscard]] base::StatusCode face_boundary_loops(
  ModelHandle model_handle,
  TopologyEntityId face_entity,
  FaceBoundaryLoops &boundary,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode face_boundary_loops(
  const FaceView &face_view,
  FaceBoundaryLoops &boundary
) noexcept;
[[nodiscard]] base::StatusCode feature_edges(
  ModelHandle model_handle,
  FeatureEdgeReport &report,
  const FeatureEdgeOptions &options = {},
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode check_topology(
  ModelHandle model_handle,
  TopologyCheckReport &report,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode free_edge_count(
  ModelHandle model_handle,
  std::size_t &count,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode topo(
  ModelHandle model_handle,
  TopoReport &report,
  const TopoOptions &options = {},
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] ModelSummary placeholder_model_summary() noexcept;

} // namespace sqmesh::geo
