// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#pragma once

#include "sqmesh/geo/api.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace sqmesh::model::detail {

struct ProxyTriangleRange final {
  std::size_t offset = 0U;
  std::size_t count = 0U;
};

struct GeometryCoarseProxyMesh final {
  std::uint64_t topology_revision = 0U;
  // Counts the sum of per-face triangulation nodes before topology-aware
  // cross-face merging compacts shared edge/vertex nodes.
  std::size_t source_local_node_count = 0U;
  std::vector<geo::Point3> nodes {};
  std::vector<std::array<double, 2>> node_uv {};       // UV parameters
  std::vector<geo::TopologyEntityId> node_uv_face {};   // face the UV belongs to
  std::vector<geo::TopologyEntityId> node_topology_owner {};
  std::vector<std::array<std::uint32_t, 3>> triangles {};
  std::vector<geo::TopologyEntityId> triangle_face_owner {};
  std::vector<ProxyTriangleRange> face_triangle_ranges {};

  // Perfect boundaries: exact sequence of global node indices for each OCC feature edge.
  std::vector<std::vector<std::uint32_t>> edge_nodes {};
  std::vector<geo::TopologyEntityId> edge_topology_owner {};

  [[nodiscard]] bool empty() const noexcept
  {
    return nodes.empty() || triangles.empty();
  }
};

} // namespace sqmesh::model::detail
