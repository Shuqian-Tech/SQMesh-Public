// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#pragma once

#include "sqmesh/base/api.hpp"
#include "sqmesh/geo/api.hpp"

#include <string_view>

namespace sqmesh::model::detail {

[[nodiscard]] base::StatusCode import_stl(
  std::string_view path,
  geo::ModelHandle &model_handle,
  const geo::StlImportOptions &options,
  base::ContextHandle context_handle
) noexcept;

} // namespace sqmesh::model::detail
