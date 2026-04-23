// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#pragma once

#include "sqmesh/base/api.hpp"
#include "sqmesh/geo/api.hpp"

#include <string_view>

class TopoDS_Shape;

namespace sqmesh::cad::occ {

[[nodiscard]] bool cad_io_available() noexcept;
[[nodiscard]] base::StatusCode store_shape(
  const TopoDS_Shape &shape,
  geo::ModelHandle &model_handle,
  base::ContextHandle context_handle
) noexcept;
[[nodiscard]] base::StatusCode import_step(
  std::string_view path,
  geo::ModelHandle &model_handle,
  const geo::StepImportOptions &options,
  base::ContextHandle context_handle
) noexcept;
[[nodiscard]] base::StatusCode import_iges(
  std::string_view path,
  geo::ModelHandle &model_handle,
  const geo::IgesImportOptions &options,
  base::ContextHandle context_handle
) noexcept;
[[nodiscard]] base::StatusCode export_step(
  geo::ModelHandle model_handle,
  std::string_view path,
  const geo::StepExportOptions &options,
  base::ContextHandle context_handle
) noexcept;
[[nodiscard]] base::StatusCode export_iges(
  geo::ModelHandle model_handle,
  std::string_view path,
  const geo::IgesExportOptions &options,
  base::ContextHandle context_handle
) noexcept;
[[nodiscard]] base::StatusCode check_topology(
  geo::ModelHandle model_handle,
  geo::TopologyCheckReport &report,
  base::ContextHandle context_handle
) noexcept;
[[nodiscard]] base::StatusCode free_edge_count(
  geo::ModelHandle model_handle,
  std::size_t &count,
  base::ContextHandle context_handle
) noexcept;
[[nodiscard]] base::StatusCode topo(
  geo::ModelHandle model_handle,
  geo::TopoReport &report,
  const geo::TopoOptions &options,
  base::ContextHandle context_handle
) noexcept;

} // namespace sqmesh::cad::occ
