// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "adapter.hpp"

#include "../../core/runtime_registry.hpp"

namespace sqmesh::cad::occ {
namespace {

base::StatusCode unavailable(std::string_view operation) noexcept
{
  if(operation == "STEP import") {
    return core::detail::publish_error(
      base::StatusCode::unsupported,
      "STEP import requires OpenCASCADE support."
    );
  }
  if(operation == "IGES import") {
    return core::detail::publish_error(
      base::StatusCode::unsupported,
      "IGES import requires OpenCASCADE support."
    );
  }
  if(operation == "STEP export") {
    return core::detail::publish_error(
      base::StatusCode::unsupported,
      "STEP export requires OpenCASCADE support."
    );
  }
  if(operation == "IGES export") {
    return core::detail::publish_error(
      base::StatusCode::unsupported,
      "IGES export requires OpenCASCADE support."
    );
  }

  return core::detail::publish_error(
    base::StatusCode::unsupported,
    "CAD IO requires OpenCASCADE support."
  );
}

} // namespace

bool cad_io_available() noexcept
{
  return false;
}

base::StatusCode store_shape(
  const TopoDS_Shape &,
  geo::ModelHandle &model_handle,
  base::ContextHandle
) noexcept
{
  model_handle = sqmesh::invalid_handle;
  return unavailable("CAD IO");
}

base::StatusCode import_step(
  std::string_view,
  geo::ModelHandle &model_handle,
  const geo::StepImportOptions &,
  base::ContextHandle
) noexcept
{
  model_handle = sqmesh::invalid_handle;
  return unavailable("STEP import");
}

base::StatusCode import_iges(
  std::string_view,
  geo::ModelHandle &model_handle,
  const geo::IgesImportOptions &,
  base::ContextHandle
) noexcept
{
  model_handle = sqmesh::invalid_handle;
  return unavailable("IGES import");
}

base::StatusCode export_step(
  geo::ModelHandle,
  std::string_view,
  const geo::StepExportOptions &,
  base::ContextHandle
) noexcept
{
  return unavailable("STEP export");
}

base::StatusCode export_iges(
  geo::ModelHandle,
  std::string_view,
  const geo::IgesExportOptions &,
  base::ContextHandle
) noexcept
{
  return unavailable("IGES export");
}

base::StatusCode check_topology(
  geo::ModelHandle,
  geo::TopologyCheckReport &report,
  base::ContextHandle
) noexcept
{
  report = {};
  return unavailable("CAD topology");
}

base::StatusCode free_edge_count(
  geo::ModelHandle,
  std::size_t &count,
  base::ContextHandle
) noexcept
{
  count = 0U;
  return unavailable("CAD topology");
}

base::StatusCode topo(
  geo::ModelHandle,
  geo::TopoReport &report,
  const geo::TopoOptions &,
  base::ContextHandle
) noexcept
{
  report = {};
  return unavailable("CAD topology");
}

} // namespace sqmesh::cad::occ
