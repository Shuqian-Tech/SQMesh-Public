// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#pragma once

#include "geometry_proxy_mesh.hpp"

#include "sqmesh/base/api.hpp"
#include "sqmesh/geo/api.hpp"

#include <cstdint>
#include <memory>
#include <vector>

namespace sqmesh::model::detail {

enum class GeometryKernel : std::uint8_t {
  none = 0,
  occ = 1,
  discrete = 2,
};

class GeometryModelStorage {
public:
  virtual ~GeometryModelStorage() = default;

  [[nodiscard]] virtual GeometryKernel kernel() const noexcept = 0;
  // Internal-only proxy cache access. Implementations copy the cached proxy
  // mesh into `proxy_mesh`, keeping the public handle/borrow contract free of
  // raw-pointer lifetime coupling.
  [[nodiscard]] virtual base::StatusCode coarse_proxy_mesh(
    GeometryCoarseProxyMesh &proxy_mesh
  ) const noexcept = 0;
  [[nodiscard]] virtual base::StatusCode topology_snapshot(
    geo::TopologySnapshot &snapshot
  ) const noexcept = 0;
  [[nodiscard]] virtual base::StatusCode topology_children(
    geo::TopologyEntityId entity,
    std::vector<geo::TopologyEntityId> &children
  ) const noexcept = 0;
  [[nodiscard]] virtual base::StatusCode topology_parents(
    geo::TopologyEntityId entity,
    std::vector<geo::TopologyEntityId> &parents
  ) const noexcept = 0;
  [[nodiscard]] virtual base::StatusCode face_uv_bounds(
    geo::TopologyEntityId face_entity,
    geo::FaceUvBounds &bounds
  ) const noexcept = 0;
  [[nodiscard]] virtual base::StatusCode sample_face(
    geo::TopologyEntityId face_entity,
    double u,
    double v,
    geo::FaceSample &sample
  ) const noexcept = 0;
  [[nodiscard]] virtual base::StatusCode sample_face_curvature(
    geo::TopologyEntityId face_entity,
    double u,
    double v,
    geo::FaceCurvatureSample &sample
  ) const noexcept = 0;
  [[nodiscard]] virtual base::StatusCode sample_face_derivatives(
    geo::TopologyEntityId face_entity,
    double u,
    double v,
    geo::FaceDerivatives &sample
  ) const noexcept = 0;
  [[nodiscard]] virtual base::StatusCode project_point_to_face(
    geo::TopologyEntityId face_entity,
    const geo::Point3 &point,
    geo::FaceProjection &projection
  ) const noexcept = 0;
  [[nodiscard]] virtual base::StatusCode recover_face_uv(
    geo::TopologyEntityId face_entity,
    const geo::Point3 &point,
    geo::FaceUvMapping &mapping
  ) const noexcept = 0;
  [[nodiscard]] virtual base::StatusCode recover_face_uv_from_edge(
    geo::TopologyEntityId face_entity,
    geo::TopologyEntityId edge_entity,
    double edge_parameter,
    geo::FaceUvMapping &mapping
  ) const noexcept = 0;
  [[nodiscard]] virtual base::StatusCode recover_face_uv_from_edge_use(
    geo::TopologyEntityId face_entity,
    const geo::FaceBoundaryEdgeUse &edge_use,
    double edge_parameter,
    geo::FaceUvMapping &mapping
  ) const noexcept = 0;
  [[nodiscard]] virtual base::StatusCode edge_curve_info(
    geo::TopologyEntityId edge_entity,
    geo::EdgeCurveInfo &info
  ) const noexcept = 0;
  [[nodiscard]] virtual base::StatusCode sample_edge_tangent(
    geo::TopologyEntityId edge_entity,
    double parameter,
    geo::EdgeTangentSample &sample
  ) const noexcept = 0;
  [[nodiscard]] virtual base::StatusCode face_boundary_loops(
    geo::TopologyEntityId face_entity,
    geo::FaceBoundaryLoops &boundary
  ) const noexcept = 0;
  [[nodiscard]] virtual base::StatusCode feature_edges(
    geo::FeatureEdgeReport &report,
    const geo::FeatureEdgeOptions &options
  ) const noexcept = 0;
};

using GeometryModelStoragePtr = std::shared_ptr<GeometryModelStorage>;

} // namespace sqmesh::model::detail
