// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#pragma once

#include "sqmesh/base/api.hpp"
#include "sqmesh/geo/api.hpp"
#include "sqmesh/mesh/api.hpp"

#include <functional>
#include <memory>
#include <string_view>

namespace sqmesh::model::detail {

class GeometryModelStorage;
using GeometryModelStoragePtr = std::shared_ptr<GeometryModelStorage>;

} // namespace sqmesh::model::detail

namespace sqmesh::core::detail {

[[nodiscard]] base::StatusCode initialize(base::ContextHandle &context_handle) noexcept;
[[nodiscard]] base::StatusCode shutdown(base::ContextHandle context_handle) noexcept;
[[nodiscard]] base::StatusCode shutdown_all() noexcept;
[[nodiscard]] bool is_initialized() noexcept;
[[nodiscard]] base::ContextHandle current_context() noexcept;
[[nodiscard]] base::SessionHandle current_session(base::ContextHandle context_handle) noexcept;
[[nodiscard]] base::StatusCode last_error_code() noexcept;
[[nodiscard]] std::string_view last_error_message() noexcept;
[[nodiscard]] const char *status_code_name(base::StatusCode code) noexcept;
[[nodiscard]] base::StatusCode clear_error_state() noexcept;
[[nodiscard]] base::StatusCode publish_error(
  base::StatusCode code,
  std::string_view message
) noexcept;
[[nodiscard]] base::StatusCode create_placeholder_model(
  geo::ModelHandle &model_handle,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode store_model(
  model::detail::GeometryModelStoragePtr storage,
  geo::ModelSummary summary,
  geo::ModelHandle &model_handle,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode update_model(
  geo::ModelHandle model_handle,
  model::detail::GeometryModelStoragePtr storage,
  geo::ModelSummary summary,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
using ModelStorageBorrowCallback =
  std::function<base::StatusCode(const model::detail::GeometryModelStorage &storage)>;
using MeshDomainBorrowCallback =
  std::function<base::StatusCode(const mesh::Domain &domain)>;
using MeshDomainMutableBorrowCallback =
  std::function<base::StatusCode(mesh::Domain &domain)>;

// Internal-only borrow helpers keep the public handle/snapshot contract intact.
// Borrowed references remain valid only for the callback invocation and must
// not escape that scope, including async use or storing their address.
[[nodiscard]] base::StatusCode with_model_storage(
  geo::ModelHandle model_handle,
  ModelStorageBorrowCallback callback,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] model::detail::GeometryModelStoragePtr lookup_model_storage(
  geo::ModelHandle model_handle,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode model_summary(
  geo::ModelHandle model_handle,
  geo::ModelSummary &summary,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode topology_snapshot(
  geo::ModelHandle model_handle,
  geo::TopologySnapshot &snapshot,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode topology_children(
  geo::ModelHandle model_handle,
  geo::TopologyEntityId entity,
  std::vector<geo::TopologyEntityId> &children,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode topology_parents(
  geo::ModelHandle model_handle,
  geo::TopologyEntityId entity,
  std::vector<geo::TopologyEntityId> &parents,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode face_uv_bounds(
  geo::ModelHandle model_handle,
  geo::TopologyEntityId face_entity,
  geo::FaceUvBounds &bounds,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode sample_face(
  geo::ModelHandle model_handle,
  geo::TopologyEntityId face_entity,
  double u,
  double v,
  geo::FaceSample &sample,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode sample_face_curvature(
  geo::ModelHandle model_handle,
  geo::TopologyEntityId face_entity,
  double u,
  double v,
  geo::FaceCurvatureSample &sample,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode sample_face_derivatives(
  geo::ModelHandle model_handle,
  geo::TopologyEntityId face_entity,
  double u,
  double v,
  geo::FaceDerivatives &sample,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode project_point_to_face(
  geo::ModelHandle model_handle,
  geo::TopologyEntityId face_entity,
  const geo::Point3 &point,
  geo::FaceProjection &projection,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode recover_face_uv(
  geo::ModelHandle model_handle,
  geo::TopologyEntityId face_entity,
  const geo::Point3 &point,
  geo::FaceUvMapping &mapping,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode recover_face_uv_from_edge(
  geo::ModelHandle model_handle,
  geo::TopologyEntityId face_entity,
  geo::TopologyEntityId edge_entity,
  double edge_parameter,
  geo::FaceUvMapping &mapping,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode recover_face_uv_from_edge_use(
  geo::ModelHandle model_handle,
  geo::TopologyEntityId face_entity,
  const geo::FaceBoundaryEdgeUse &edge_use,
  double edge_parameter,
  geo::FaceUvMapping &mapping,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode edge_curve_info(
  geo::ModelHandle model_handle,
  geo::TopologyEntityId edge_entity,
  geo::EdgeCurveInfo &info,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode sample_edge_tangent(
  geo::ModelHandle model_handle,
  geo::TopologyEntityId edge_entity,
  double parameter,
  geo::EdgeTangentSample &sample,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode face_boundary_loops(
  geo::ModelHandle model_handle,
  geo::TopologyEntityId face_entity,
  geo::FaceBoundaryLoops &boundary,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode feature_edges(
  geo::ModelHandle model_handle,
  geo::FeatureEdgeReport &report,
  const geo::FeatureEdgeOptions &options = {},
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode validate_model_handle(
  geo::ModelHandle model_handle,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode store_mesh(
  geo::ModelHandle model_handle,
  std::string_view algorithm_name,
  mesh::MeshSummary summary,
  std::shared_ptr<mesh::Domain> domain,
  mesh::MeshHandle &mesh_handle,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode mesh_summary(
  mesh::MeshHandle mesh_handle,
  mesh::MeshSummary &summary,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode mesh_domain_snapshot(
  mesh::MeshHandle mesh_handle,
  mesh::Domain &domain,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode with_mesh_domain(
  mesh::MeshHandle mesh_handle,
  MeshDomainBorrowCallback callback,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode with_mesh_domain_mutable(
  mesh::MeshHandle mesh_handle,
  MeshDomainMutableBorrowCallback callback,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode refresh_mesh_summary(
  mesh::MeshHandle mesh_handle,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode model_proxy_mesh(
  geo::ModelHandle model_handle,
  mesh::MeshHandle &mesh_handle,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode mesh_node_count(
  mesh::MeshHandle mesh_handle,
  std::size_t &count,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode mesh_cell_count(
  mesh::MeshHandle mesh_handle,
  std::size_t &count,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;

} // namespace sqmesh::core::detail
