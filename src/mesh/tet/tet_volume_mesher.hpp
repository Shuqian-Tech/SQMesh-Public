// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#pragma once

#include "../framework/meshing_framework.hpp"

namespace sqmesh::mesh::detail {

[[nodiscard]] MeshingAlgorithmPtr create_native_tet_volume_mesher();

} // namespace sqmesh::mesh::detail
