// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#pragma once

#include "sqmesh/mesh/api.hpp"
#include <array>
#include <cstdint>
#include <vector>

namespace sqmesh::mesh {

struct MeshRegion {
  int id = -1;
  std::vector<std::uint32_t> boundary_face_ids; // Global face indices in the extracted list
  std::array<double, 3> interior_point {0, 0, 0};
};

struct RegionDetectionResult {
  std::vector<MeshRegion> regions;
  // face_to_regions[face_id] = {positive_side_region, negative_side_region}
  // -1 = external (unbounded)
  std::vector<std::array<int, 2>> face_to_regions;

  // Mapping: extracted face index -> (entity_group_id, face_index_in_entity_group)
  std::vector<EntityRef> face_refs;
};


// Detect enclosed regions from the surface mesh in the Domain.
// Uses face connectivity flood-fill + ray casting (even-odd rule).
RegionDetectionResult detect_regions(const Domain &domain);

} // namespace sqmesh::mesh
