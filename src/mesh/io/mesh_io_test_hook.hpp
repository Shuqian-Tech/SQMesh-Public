// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#pragma once

#include "sqmesh/base/api.hpp"
#include "sqmesh/mesh/api.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace sqmesh::mesh::testing {

struct ObjReviewOptionalIndex final {
  bool present = false;
  std::size_t value = 0U;
};

struct ObjReviewCoordinate final {
  std::array<double, 3> values {0.0, 0.0, 0.0};
  std::uint8_t component_count = 0U;
};

struct ObjReviewMaterialLibrary final {
  std::vector<std::string> entries {};
};

struct ObjReviewFaceCorner final {
  std::size_t vertex_index = 0U;
  ObjReviewOptionalIndex texture_index {};
  ObjReviewOptionalIndex normal_index {};
};

struct ObjReviewEdge final {
  std::array<std::size_t, 2> nodes {};
  std::string material_name {};
};

struct ObjReviewFace final {
  std::array<ObjReviewFaceCorner, 3> corners {};
  std::string material_name {};
};

struct ObjImportReview final {
  std::string name {};
  std::vector<ObjReviewCoordinate> texture_coordinates {};
  std::vector<std::array<double, 3>> normals {};
  std::vector<ObjReviewCoordinate> parameter_vertices {};
  std::vector<ObjReviewMaterialLibrary> material_libraries {};
  std::vector<ObjReviewEdge> edges {};
  std::vector<ObjReviewFace> faces {};
};

// Internal review/test seam for bounded OBJ parser fidelity checks. This header
// intentionally lives under `src/sdk` and is not part of the installed SDK.
[[nodiscard]] base::StatusCode inspect_obj_import_review(
  std::string_view path,
  ObjImportReview &review,
  const ObjImportOptions &options = {}
) noexcept;

}  // namespace sqmesh::mesh::testing
