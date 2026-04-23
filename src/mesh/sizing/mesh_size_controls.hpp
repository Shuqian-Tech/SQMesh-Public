// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#pragma once

#include "sqmesh/mesh/api.hpp"

#include <vector>

namespace sqmesh::mesh::detail {

struct ResolvedMeshSizeControls final {
  double global_target_size = 0.0;
  std::uint64_t topology_revision = 0;
  std::vector<double> face_target_sizes {};
  std::vector<double> edge_target_sizes {};
};

[[nodiscard]] base::StatusCode resolve_mesh_size_controls(
  const geo::ModelView &geometry_view,
  const MeshSizeControls &controls,
  double global_target_size,
  ResolvedMeshSizeControls &resolved
);
[[nodiscard]] double effective_face_target_size(
  const ResolvedMeshSizeControls &controls,
  geo::TopologyEntityId face
) noexcept;
[[nodiscard]] double effective_edge_target_size(
  const ResolvedMeshSizeControls &controls,
  const geo::EdgeView &edge_view
) noexcept;

} // namespace sqmesh::mesh::detail
