// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "mesh_size_controls.hpp"

#include "core/runtime_registry.hpp"

#include <algorithm>
#include <cmath>

namespace sqmesh::mesh::detail {
namespace {

constexpr double kMinimumTargetSize = 1.0e-9;

void assign_minimum_target(
  std::vector<double> &targets,
  std::uint32_t index,
  double target_size
)
{
  auto &stored_target_size = targets[index];
  if(stored_target_size > 0.0) {
    stored_target_size = std::min(stored_target_size, target_size);
    return;
  }

  stored_target_size = target_size;
}

[[nodiscard]] double apply_target_size(double current, double candidate) noexcept
{
  if(candidate <= 0.0) {
    return current;
  }
  return std::min(current, candidate);
}

} // namespace

base::StatusCode resolve_mesh_size_controls(
  const geo::ModelView &geometry_view,
  const MeshSizeControls &controls,
  double global_target_size,
  ResolvedMeshSizeControls &resolved
)
{
  resolved = {};
  resolved.global_target_size = std::max(global_target_size, kMinimumTargetSize);
  resolved.topology_revision = geometry_view.snapshot.topology_revision;
  resolved.face_target_sizes.assign(geometry_view.faces.size(), 0.0);
  resolved.edge_target_sizes.assign(geometry_view.edges.size(), 0.0);

  if(controls.empty()) {
    return core::detail::clear_error_state();
  }

  if(controls.topology_revision != 0U &&
     controls.topology_revision != geometry_view.snapshot.topology_revision) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Local mesh size controls reference a different geometry topology revision."
    );
  }

  for(const auto &local_size : controls.local_sizes) {
    if(!geo::is_valid(local_size.entity)) {
      return core::detail::publish_error(
        base::StatusCode::invalid_argument,
        "Local mesh size controls require valid geometry topology entity ids."
      );
    }
    if(local_size.target_size <= 0.0 || !std::isfinite(local_size.target_size)) {
      return core::detail::publish_error(
        base::StatusCode::invalid_argument,
        "Local mesh size controls require a positive finite target size."
      );
    }

    switch(local_size.entity.dimension) {
    case geo::TopologyDimension::face:
      if(local_size.entity.index >= geometry_view.faces.size()) {
        return core::detail::publish_error(
          base::StatusCode::invalid_argument,
          "A local mesh size control references a non-existent geometry face."
        );
      }
      assign_minimum_target(
        resolved.face_target_sizes,
        local_size.entity.index,
        local_size.target_size
      );
      break;
    case geo::TopologyDimension::edge:
      if(local_size.entity.index >= geometry_view.edges.size()) {
        return core::detail::publish_error(
          base::StatusCode::invalid_argument,
          "A local mesh size control references a non-existent geometry edge."
        );
      }
      assign_minimum_target(
        resolved.edge_target_sizes,
        local_size.entity.index,
        local_size.target_size
      );
      break;
    case geo::TopologyDimension::vertex:
    case geo::TopologyDimension::region:
    default:
      return core::detail::publish_error(
        base::StatusCode::invalid_argument,
        "Local mesh size controls currently support only geometry face and edge ids."
      );
    }
  }

  return core::detail::clear_error_state();
}

double effective_face_target_size(
  const ResolvedMeshSizeControls &controls,
  geo::TopologyEntityId face
) noexcept
{
  double target_size = controls.global_target_size;
  if(face.dimension == geo::TopologyDimension::face &&
     face.index < controls.face_target_sizes.size()) {
    target_size =
      apply_target_size(target_size, controls.face_target_sizes[face.index]);
  }
  return std::max(target_size, kMinimumTargetSize);
}

double effective_edge_target_size(
  const ResolvedMeshSizeControls &controls,
  const geo::EdgeView &edge_view
) noexcept
{
  double target_size = controls.global_target_size;

  if(edge_view.entity.dimension == geo::TopologyDimension::edge &&
     edge_view.entity.index < controls.edge_target_sizes.size()) {
    target_size =
      apply_target_size(target_size, controls.edge_target_sizes[edge_view.entity.index]);
  }

  for(const auto face_id : edge_view.face_ids) {
    target_size = apply_target_size(
      target_size,
      effective_face_target_size(controls, face_id)
    );
  }

  return std::max(target_size, kMinimumTargetSize);
}

} // namespace sqmesh::mesh::detail
