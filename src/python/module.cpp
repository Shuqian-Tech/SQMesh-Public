// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "bindings.hpp"

#include "sqmesh/base/api.hpp"

PYBIND11_MODULE(_sqmesh, module)
{
  module.doc() = "Bounded Python bindings for the SQMesh C++ SDK.";
  module.attr("invalid_handle") = py::int_(sqmesh::invalid_handle);

  register_sqmesh_exceptions(module);
  bind_base(module);
  bind_geo(module);
  bind_mesh(module);
}
