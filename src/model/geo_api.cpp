// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "../core/runtime_registry.hpp"
#include "../cad/occ/adapter.hpp"

#include "stl_import.hpp"

#include "sqmesh/geo/api.hpp"

#include <algorithm>
#include <cmath>
#include <string_view>
#include <unordered_set>

namespace sqmesh::geo {
namespace {

void sort_topology_entities(std::vector<TopologyEntityId> &entities)
{
  std::sort(
    entities.begin(),
    entities.end(),
    [](TopologyEntityId lhs, TopologyEntityId rhs) {
      if(lhs.dimension != rhs.dimension) {
        return static_cast<std::uint8_t>(lhs.dimension) <
               static_cast<std::uint8_t>(rhs.dimension);
      }

      return lhs.index < rhs.index;
    }
  );
}

template <typename ViewType>
void set_view_context(
  ModelHandle model_handle,
  base::ContextHandle context_handle,
  ViewType &view
) noexcept
{
  view.model_handle = model_handle;
  view.context_handle = context_handle;
}

[[nodiscard]] base::StatusCode invalid_face_view(std::string_view message) noexcept
{
  return core::detail::publish_error(base::StatusCode::invalid_argument, message);
}

[[nodiscard]] base::StatusCode invalid_edge_view(std::string_view message) noexcept
{
  return core::detail::publish_error(base::StatusCode::invalid_argument, message);
}

[[nodiscard]] base::StatusCode invalid_edge_sampling_options(
  std::string_view message
) noexcept
{
  return core::detail::publish_error(base::StatusCode::invalid_argument, message);
}

void finalize_face_boundary_loop(FaceBoundaryLoop &loop) noexcept
{
  loop.vertex_ids.clear();
  loop.continuous = false;
  loop.seam_edge_use_count = 0U;
  loop.degenerate_edge_use_count = 0U;
  loop.repeated_edge_use_count = 0U;

  if(loop.edge_uses.empty()) {
    loop.closed = false;
    return;
  }

  std::unordered_set<std::uint64_t> seen_edges;
  seen_edges.reserve(loop.edge_uses.size());

  loop.closed = loop.edge_uses.front().start_vertex == loop.edge_uses.back().end_vertex;
  loop.continuous = true;
  loop.vertex_ids.reserve(loop.edge_uses.size() + (loop.closed ? 0U : 1U));
  loop.vertex_ids.push_back(loop.edge_uses.front().start_vertex);

  for(std::size_t edge_index = 0U; edge_index < loop.edge_uses.size(); ++edge_index) {
    const auto &edge_use = loop.edge_uses[edge_index];
    if(edge_use.is_seam) {
      ++loop.seam_edge_use_count;
    }
    if(edge_use.is_degenerate) {
      ++loop.degenerate_edge_use_count;
    }

    const auto edge_key =
      (static_cast<std::uint64_t>(static_cast<std::uint8_t>(edge_use.edge.dimension)) << 32U) |
      static_cast<std::uint64_t>(edge_use.edge.index);
    if(!seen_edges.insert(edge_key).second) {
      ++loop.repeated_edge_use_count;
    }

    if(edge_index > 0U) {
      const auto &previous_edge_use = loop.edge_uses[edge_index - 1U];
      if(previous_edge_use.end_vertex != edge_use.start_vertex) {
        loop.continuous = false;
      }
    }

    const bool closes_cycle =
      loop.closed &&
      edge_index + 1U == loop.edge_uses.size() &&
      edge_use.end_vertex == loop.edge_uses.front().start_vertex;
    if(!closes_cycle) {
      loop.vertex_ids.push_back(edge_use.end_vertex);
    }
  }
}

void finalize_face_boundary_loops(FaceBoundaryLoops &boundary) noexcept
{
  boundary.primary_outer_loop_index = invalid_boundary_loop_index;
  boundary.outer_loop_count = 0U;
  boundary.inner_loop_count = 0U;
  boundary.unknown_loop_count = 0U;
  boundary.closed_loop_count = 0U;
  boundary.open_loop_count = 0U;
  boundary.non_continuous_loop_count = 0U;
  boundary.seam_loop_count = 0U;
  boundary.seam_edge_use_count = 0U;
  boundary.degenerate_loop_count = 0U;
  boundary.degenerate_edge_use_count = 0U;
  boundary.repeated_edge_loop_count = 0U;
  boundary.repeated_edge_use_count = 0U;
  boundary.has_holes = false;
  boundary.has_seams = false;

  for(std::size_t loop_index = 0U; loop_index < boundary.loops.size(); ++loop_index) {
    auto &loop = boundary.loops[loop_index];
    finalize_face_boundary_loop(loop);

    if(loop.closed) {
      ++boundary.closed_loop_count;
    } else {
      ++boundary.open_loop_count;
    }

    if(!loop.continuous) {
      ++boundary.non_continuous_loop_count;
    }

    if(loop.seam_edge_use_count > 0U) {
      ++boundary.seam_loop_count;
      boundary.seam_edge_use_count += loop.seam_edge_use_count;
    }

    if(loop.degenerate_edge_use_count > 0U) {
      ++boundary.degenerate_loop_count;
      boundary.degenerate_edge_use_count += loop.degenerate_edge_use_count;
    }

    if(loop.repeated_edge_use_count > 0U) {
      ++boundary.repeated_edge_loop_count;
      boundary.repeated_edge_use_count += loop.repeated_edge_use_count;
    }

    switch(loop.kind) {
    case FaceBoundaryLoopKind::outer:
      if(boundary.primary_outer_loop_index == invalid_boundary_loop_index) {
        boundary.primary_outer_loop_index = loop_index;
      }
      ++boundary.outer_loop_count;
      break;
    case FaceBoundaryLoopKind::inner:
      ++boundary.inner_loop_count;
      break;
    case FaceBoundaryLoopKind::unknown:
    default:
      ++boundary.unknown_loop_count;
      break;
    }
  }

  boundary.has_holes = boundary.inner_loop_count > 0U;
  boundary.has_seams = boundary.seam_edge_use_count > 0U;
}

[[nodiscard]] std::size_t edge_curve_segment_count(
  const EdgeCurveInfo &info,
  const EdgeCurveSamplingOptions &options
) noexcept
{
  std::size_t segment_count = std::max<std::size_t>(1U, options.min_segment_count);
  if(options.target_segment_length > 0.0 &&
     std::isfinite(info.approximate_length) &&
     info.approximate_length > 0.0) {
    segment_count =
      std::max<std::size_t>(
        segment_count,
        static_cast<std::size_t>(std::ceil(
          info.approximate_length / options.target_segment_length
        ))
      );
  }
  return segment_count;
}

} // namespace

std::string_view module_name() noexcept
{
  return "sqmesh::geo";
}

bool cad_io_available() noexcept
{
  return cad::occ::cad_io_available();
}

base::StatusCode create_placeholder_model(
  ModelHandle &model_handle,
  base::ContextHandle context_handle
) noexcept
{
  return core::detail::create_placeholder_model(model_handle, context_handle);
}

base::StatusCode import_step(
  std::string_view path,
  ModelHandle &model_handle,
  const StepImportOptions &options,
  base::ContextHandle context_handle
) noexcept
{
  return cad::occ::import_step(path, model_handle, options, context_handle);
}

base::StatusCode import_iges(
  std::string_view path,
  ModelHandle &model_handle,
  const IgesImportOptions &options,
  base::ContextHandle context_handle
) noexcept
{
  return cad::occ::import_iges(path, model_handle, options, context_handle);
}

base::StatusCode import_stl(
  std::string_view path,
  ModelHandle &model_handle,
  const StlImportOptions &options,
  base::ContextHandle context_handle
) noexcept
{
  return model::detail::import_stl(path, model_handle, options, context_handle);
}

base::StatusCode export_step(
  ModelHandle model_handle,
  std::string_view path,
  const StepExportOptions &options,
  base::ContextHandle context_handle
) noexcept
{
  return cad::occ::export_step(model_handle, path, options, context_handle);
}

base::StatusCode export_iges(
  ModelHandle model_handle,
  std::string_view path,
  const IgesExportOptions &options,
  base::ContextHandle context_handle
) noexcept
{
  return cad::occ::export_iges(model_handle, path, options, context_handle);
}

base::StatusCode model_summary(
  ModelHandle model_handle,
  ModelSummary &summary,
  base::ContextHandle context_handle
) noexcept
{
  return core::detail::model_summary(model_handle, summary, context_handle);
}

base::StatusCode model_proxy_mesh(
  ModelHandle model_handle,
  sqmesh::Handle &mesh_handle,
  base::ContextHandle context_handle
) noexcept
{
  return core::detail::model_proxy_mesh(model_handle, mesh_handle, context_handle);
}

base::StatusCode topology_snapshot(
  ModelHandle model_handle,
  TopologySnapshot &snapshot,
  base::ContextHandle context_handle
) noexcept
{
  return core::detail::topology_snapshot(model_handle, snapshot, context_handle);
}

base::StatusCode model_view(
  ModelHandle model_handle,
  ModelView &view,
  base::ContextHandle context_handle
) noexcept
{
  ModelView built;

  auto status =
    topology_snapshot(model_handle, built.snapshot, context_handle);
  if(status != base::StatusCode::ok) {
    view = {};
    return status;
  }

  built.model_handle = model_handle;
  built.context_handle = context_handle;
  built.regions.resize(built.snapshot.regions.size());
  built.faces.resize(built.snapshot.faces.size());
  built.edges.resize(built.snapshot.edges.size());
  built.vertices.resize(built.snapshot.vertices.size());

  for(const auto &region_info : built.snapshot.regions) {
    auto &region = built.regions[region_info.entity.index];
    set_view_context(model_handle, context_handle, region);
    region.entity = region_info.entity;

    status = topology_children(
      model_handle,
      region.entity,
      region.face_ids,
      context_handle
    );
    if(status != base::StatusCode::ok) {
      view = {};
      return status;
    }

    sort_topology_entities(region.face_ids);
  }

  for(const auto &face_info : built.snapshot.faces) {
    auto &face = built.faces[face_info.entity.index];
    set_view_context(model_handle, context_handle, face);
    face.entity = face_info.entity;

    status = topology_parents(
      model_handle,
      face.entity,
      face.region_ids,
      context_handle
    );
    if(status != base::StatusCode::ok) {
      view = {};
      return status;
    }
    sort_topology_entities(face.region_ids);

    status = topology_children(
      model_handle,
      face.entity,
      face.edge_ids,
      context_handle
    );
    if(status != base::StatusCode::ok) {
      view = {};
      return status;
    }
    sort_topology_entities(face.edge_ids);

    status = face_boundary_loops(
      model_handle,
      face.entity,
      face.ordered_boundary,
      context_handle
    );
    if(status != base::StatusCode::ok) {
      view = {};
      return status;
    }

    if(face.region_ids.empty()) {
      built.root_face_ids.push_back(face.entity);
    }
  }

  for(const auto &edge_info : built.snapshot.edges) {
    auto &edge = built.edges[edge_info.entity.index];
    set_view_context(model_handle, context_handle, edge);
    edge.entity = edge_info.entity;

    status = topology_parents(
      model_handle,
      edge.entity,
      edge.face_ids,
      context_handle
    );
    if(status != base::StatusCode::ok) {
      view = {};
      return status;
    }
    sort_topology_entities(edge.face_ids);

    status = topology_children(
      model_handle,
      edge.entity,
      edge.vertex_ids,
      context_handle
    );
    if(status != base::StatusCode::ok) {
      view = {};
      return status;
    }
    sort_topology_entities(edge.vertex_ids);

    if(edge.face_ids.empty()) {
      built.root_edge_ids.push_back(edge.entity);
    }
  }

  for(const auto &vertex_info : built.snapshot.vertices) {
    auto &vertex = built.vertices[vertex_info.entity.index];
    set_view_context(model_handle, context_handle, vertex);
    vertex.entity = vertex_info.entity;

    status = topology_parents(
      model_handle,
      vertex.entity,
      vertex.edge_ids,
      context_handle
    );
    if(status != base::StatusCode::ok) {
      view = {};
      return status;
    }
    sort_topology_entities(vertex.edge_ids);

    if(vertex.edge_ids.empty()) {
      built.root_vertex_ids.push_back(vertex.entity);
    }
  }

  sort_topology_entities(built.root_face_ids);
  sort_topology_entities(built.root_edge_ids);
  sort_topology_entities(built.root_vertex_ids);
  view = std::move(built);
  return core::detail::clear_error_state();
}

base::StatusCode topology_children(
  ModelHandle model_handle,
  TopologyEntityId entity,
  std::vector<TopologyEntityId> &children,
  base::ContextHandle context_handle
) noexcept
{
  return core::detail::topology_children(
    model_handle,
    entity,
    children,
    context_handle
  );
}

base::StatusCode topology_parents(
  ModelHandle model_handle,
  TopologyEntityId entity,
  std::vector<TopologyEntityId> &parents,
  base::ContextHandle context_handle
) noexcept
{
  return core::detail::topology_parents(
    model_handle,
    entity,
    parents,
    context_handle
  );
}

base::StatusCode face_uv_bounds(
  ModelHandle model_handle,
  TopologyEntityId face_entity,
  FaceUvBounds &bounds,
  base::ContextHandle context_handle
) noexcept
{
  return core::detail::face_uv_bounds(
    model_handle,
    face_entity,
    bounds,
    context_handle
  );
}

base::StatusCode face_uv_bounds(
  const FaceView &face_view,
  FaceUvBounds &bounds
) noexcept
{
  if(face_view.entity.dimension != TopologyDimension::face || !is_valid(face_view.entity)) {
    bounds = {};
    return invalid_face_view("Face UV bounds lookup requires a valid face view.");
  }

  return face_uv_bounds(
    face_view.model_handle,
    face_view.entity,
    bounds,
    face_view.context_handle
  );
}

base::StatusCode sample_face(
  ModelHandle model_handle,
  TopologyEntityId face_entity,
  double u,
  double v,
  FaceSample &sample,
  base::ContextHandle context_handle
) noexcept
{
  return core::detail::sample_face(
    model_handle,
    face_entity,
    u,
    v,
    sample,
    context_handle
  );
}

base::StatusCode sample_face(
  const FaceView &face_view,
  double u,
  double v,
  FaceSample &sample
) noexcept
{
  if(face_view.entity.dimension != TopologyDimension::face || !is_valid(face_view.entity)) {
    sample = {};
    return invalid_face_view("Face sampling requires a valid face view.");
  }

  return sample_face(
    face_view.model_handle,
    face_view.entity,
    u,
    v,
    sample,
    face_view.context_handle
  );
}

base::StatusCode sample_face_curvature(
  ModelHandle model_handle,
  TopologyEntityId face_entity,
  double u,
  double v,
  FaceCurvatureSample &sample,
  base::ContextHandle context_handle
) noexcept
{
  return core::detail::sample_face_curvature(
    model_handle,
    face_entity,
    u,
    v,
    sample,
    context_handle
  );
}

base::StatusCode sample_face_curvature(
  const FaceView &face_view,
  double u,
  double v,
  FaceCurvatureSample &sample
) noexcept
{
  if(face_view.entity.dimension != TopologyDimension::face || !is_valid(face_view.entity)) {
    sample = {};
    return invalid_face_view("Face curvature sampling requires a valid face view.");
  }

  return sample_face_curvature(
    face_view.model_handle,
    face_view.entity,
    u,
    v,
    sample,
    face_view.context_handle
  );
}

base::StatusCode sample_face_derivatives(
  ModelHandle model_handle,
  TopologyEntityId face_entity,
  double u,
  double v,
  FaceDerivatives &sample,
  base::ContextHandle context_handle
) noexcept
{
  return core::detail::sample_face_derivatives(
    model_handle,
    face_entity,
    u,
    v,
    sample,
    context_handle
  );
}

base::StatusCode sample_face_derivatives(
  const FaceView &face_view,
  double u,
  double v,
  FaceDerivatives &sample
) noexcept
{
  if(face_view.entity.dimension != TopologyDimension::face || !is_valid(face_view.entity)) {
    sample = {};
    return invalid_face_view("Face derivative sampling requires a valid face view.");
  }

  return sample_face_derivatives(
    face_view.model_handle,
    face_view.entity,
    u,
    v,
    sample,
    face_view.context_handle
  );
}

base::StatusCode project_point_to_face(
  ModelHandle model_handle,
  TopologyEntityId face_entity,
  const Point3 &point,
  FaceProjection &projection,
  base::ContextHandle context_handle
) noexcept
{
  return core::detail::project_point_to_face(
    model_handle,
    face_entity,
    point,
    projection,
    context_handle
  );
}

base::StatusCode project_point_to_face(
  const FaceView &face_view,
  const Point3 &point,
  FaceProjection &projection
) noexcept
{
  if(face_view.entity.dimension != TopologyDimension::face || !is_valid(face_view.entity)) {
    projection = {};
    return invalid_face_view("Face projection requires a valid face view.");
  }

  return project_point_to_face(
    face_view.model_handle,
    face_view.entity,
    point,
    projection,
    face_view.context_handle
  );
}

base::StatusCode recover_face_uv(
  ModelHandle model_handle,
  TopologyEntityId face_entity,
  const Point3 &point,
  FaceUvMapping &mapping,
  base::ContextHandle context_handle
) noexcept
{
  return core::detail::recover_face_uv(
    model_handle,
    face_entity,
    point,
    mapping,
    context_handle
  );
}

base::StatusCode recover_face_uv(
  const FaceView &face_view,
  const Point3 &point,
  FaceUvMapping &mapping
) noexcept
{
  if(face_view.entity.dimension != TopologyDimension::face || !is_valid(face_view.entity)) {
    mapping = {};
    return invalid_face_view("Face UV recovery requires a valid face view.");
  }

  return recover_face_uv(
    face_view.model_handle,
    face_view.entity,
    point,
    mapping,
    face_view.context_handle
  );
}

base::StatusCode recover_face_uv_from_edge(
  ModelHandle model_handle,
  TopologyEntityId face_entity,
  TopologyEntityId edge_entity,
  double edge_parameter,
  FaceUvMapping &mapping,
  base::ContextHandle context_handle
) noexcept
{
  return core::detail::recover_face_uv_from_edge(
    model_handle,
    face_entity,
    edge_entity,
    edge_parameter,
    mapping,
    context_handle
  );
}

base::StatusCode recover_face_uv_from_edge(
  const FaceView &face_view,
  TopologyEntityId edge_entity,
  double edge_parameter,
  FaceUvMapping &mapping
) noexcept
{
  if(face_view.entity.dimension != TopologyDimension::face || !is_valid(face_view.entity)) {
    mapping = {};
    return invalid_face_view("Face UV recovery from edge requires a valid face view.");
  }

  return recover_face_uv_from_edge(
    face_view.model_handle,
    face_view.entity,
    edge_entity,
    edge_parameter,
    mapping,
    face_view.context_handle
  );
}

base::StatusCode recover_face_uv_from_edge_use(
  ModelHandle model_handle,
  TopologyEntityId face_entity,
  const FaceBoundaryEdgeUse &edge_use,
  double edge_parameter,
  FaceUvMapping &mapping,
  base::ContextHandle context_handle
) noexcept
{
  return core::detail::recover_face_uv_from_edge_use(
    model_handle,
    face_entity,
    edge_use,
    edge_parameter,
    mapping,
    context_handle
  );
}

base::StatusCode recover_face_uv_from_edge_use(
  const FaceView &face_view,
  const FaceBoundaryEdgeUse &edge_use,
  double edge_parameter,
  FaceUvMapping &mapping
) noexcept
{
  if(face_view.entity.dimension != TopologyDimension::face || !is_valid(face_view.entity)) {
    mapping = {};
    return invalid_face_view("Face UV recovery from edge use requires a valid face view.");
  }

  return recover_face_uv_from_edge_use(
    face_view.model_handle,
    face_view.entity,
    edge_use,
    edge_parameter,
    mapping,
    face_view.context_handle
  );
}

base::StatusCode edge_curve_info(
  ModelHandle model_handle,
  TopologyEntityId edge_entity,
  EdgeCurveInfo &info,
  base::ContextHandle context_handle
) noexcept
{
  return core::detail::edge_curve_info(
    model_handle,
    edge_entity,
    info,
    context_handle
  );
}

base::StatusCode edge_curve_info(
  const EdgeView &edge_view,
  EdgeCurveInfo &info
) noexcept
{
  if(edge_view.entity.dimension != TopologyDimension::edge || !is_valid(edge_view.entity)) {
    info = {};
    return invalid_edge_view("Edge curve lookup requires a valid edge view.");
  }

  return edge_curve_info(
    edge_view.model_handle,
    edge_view.entity,
    info,
    edge_view.context_handle
  );
}

base::StatusCode sample_edge_tangent(
  ModelHandle model_handle,
  TopologyEntityId edge_entity,
  double parameter,
  EdgeTangentSample &sample,
  base::ContextHandle context_handle
) noexcept
{
  return core::detail::sample_edge_tangent(
    model_handle,
    edge_entity,
    parameter,
    sample,
    context_handle
  );
}

base::StatusCode sample_edge_tangent(
  const EdgeView &edge_view,
  double parameter,
  EdgeTangentSample &sample
) noexcept
{
  if(edge_view.entity.dimension != TopologyDimension::edge || !is_valid(edge_view.entity)) {
    sample = {};
    return invalid_edge_view("Edge tangent sampling requires a valid edge view.");
  }

  return sample_edge_tangent(
    edge_view.model_handle,
    edge_view.entity,
    parameter,
    sample,
    edge_view.context_handle
  );
}

base::StatusCode sample_edge_curve(
  ModelHandle model_handle,
  TopologyEntityId edge_entity,
  const EdgeCurveSamplingOptions &options,
  EdgeCurveSamples &samples,
  base::ContextHandle context_handle
) noexcept
{
  samples = {};

  if(edge_entity.dimension != TopologyDimension::edge || !is_valid(edge_entity)) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Edge curve sampling requires a valid edge topology entity id."
    );
  }
  if(options.min_segment_count == 0U) {
    return invalid_edge_sampling_options(
      "Edge curve sampling requires min_segment_count to be at least one."
    );
  }
  if(options.target_segment_length < 0.0) {
    return invalid_edge_sampling_options(
      "Edge curve sampling requires a non-negative target_segment_length."
    );
  }

  auto status = edge_curve_info(
    model_handle,
    edge_entity,
    samples.curve,
    context_handle
  );
  if(status != base::StatusCode::ok) {
    samples = {};
    return status;
  }
  if(samples.curve.parameter_max < samples.curve.parameter_min) {
    samples = {};
    return core::detail::publish_error(
      base::StatusCode::unsupported,
      "Edge curve sampling requires a non-decreasing edge parameter range."
    );
  }

  const auto segment_count = edge_curve_segment_count(samples.curve, options);
  samples.samples.reserve(segment_count + 1U);
  for(std::size_t segment_index = 0U; segment_index <= segment_count; ++segment_index) {
    const double alpha =
      static_cast<double>(segment_index) / static_cast<double>(segment_count);
    const double parameter = samples.curve.parameter_min +
                             alpha * (samples.curve.parameter_max - samples.curve.parameter_min);

    EdgeTangentSample sample;
    status = sample_edge_tangent(
      model_handle,
      edge_entity,
      parameter,
      sample,
      context_handle
    );
    if(status != base::StatusCode::ok) {
      samples = {};
      return status;
    }

    samples.samples.push_back(sample);
  }

  return core::detail::clear_error_state();
}

base::StatusCode sample_edge_curve(
  const EdgeView &edge_view,
  const EdgeCurveSamplingOptions &options,
  EdgeCurveSamples &samples
) noexcept
{
  if(edge_view.entity.dimension != TopologyDimension::edge || !is_valid(edge_view.entity)) {
    samples = {};
    return invalid_edge_view("Edge curve sampling requires a valid edge view.");
  }

  return sample_edge_curve(
    edge_view.model_handle,
    edge_view.entity,
    options,
    samples,
    edge_view.context_handle
  );
}

base::StatusCode face_boundary_loops(
  ModelHandle model_handle,
  TopologyEntityId face_entity,
  FaceBoundaryLoops &boundary,
  base::ContextHandle context_handle
) noexcept
{
  const auto status = core::detail::face_boundary_loops(
    model_handle,
    face_entity,
    boundary,
    context_handle
  );
  if(status != base::StatusCode::ok) {
    boundary = {};
    return status;
  }

  finalize_face_boundary_loops(boundary);
  return core::detail::clear_error_state();
}

base::StatusCode face_boundary_loops(
  const FaceView &face_view,
  FaceBoundaryLoops &boundary
) noexcept
{
  if(face_view.entity.dimension != TopologyDimension::face || !is_valid(face_view.entity)) {
    boundary = {};
    return invalid_face_view("Face boundary lookup requires a valid face view.");
  }

  boundary = face_view.ordered_boundary;
  finalize_face_boundary_loops(boundary);
  return core::detail::clear_error_state();
}

base::StatusCode feature_edges(
  ModelHandle model_handle,
  FeatureEdgeReport &report,
  const FeatureEdgeOptions &options,
  base::ContextHandle context_handle
) noexcept
{
  return core::detail::feature_edges(
    model_handle,
    report,
    options,
    context_handle
  );
}

base::StatusCode check_topology(
  ModelHandle model_handle,
  TopologyCheckReport &report,
  base::ContextHandle context_handle
) noexcept
{
  return cad::occ::check_topology(model_handle, report, context_handle);
}

base::StatusCode free_edge_count(
  ModelHandle model_handle,
  std::size_t &count,
  base::ContextHandle context_handle
) noexcept
{
  return cad::occ::free_edge_count(model_handle, count, context_handle);
}

base::StatusCode topo(
  ModelHandle model_handle,
  TopoReport &report,
  const TopoOptions &options,
  base::ContextHandle context_handle
) noexcept
{
  return cad::occ::topo(model_handle, report, options, context_handle);
}

ModelSummary placeholder_model_summary() noexcept
{
  return {};
}

} // namespace sqmesh::geo
