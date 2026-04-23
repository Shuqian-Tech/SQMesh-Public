// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#pragma once

#include "../framework/meshing_framework.hpp"

#include <string_view>

namespace sqmesh::mesh::detail {

constexpr const char *kBoundaryLayerMesherAlias = "boundary_layer_mesher";

[[nodiscard]] MeshingAlgorithmPtr create_boundary_layer_mesher();

} // namespace sqmesh::mesh::detail
