// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "bindings.hpp"
#include "exceptions.hpp"

#include "sqmesh/base/api.hpp"

void bind_base(py::module_ &module)
{
  auto base_module = module.def_submodule("base", "SQMesh runtime bindings.");
  base_module.attr("invalid_handle") = py::int_(sqmesh::invalid_handle);

  py::enum_<sqmesh::base::HandleKind>(base_module, "HandleKind")
    .value("invalid", sqmesh::base::HandleKind::invalid)
    .value("context", sqmesh::base::HandleKind::context)
    .value("session", sqmesh::base::HandleKind::session)
    .value("model", sqmesh::base::HandleKind::model)
    .value("mesh", sqmesh::base::HandleKind::mesh);

  py::enum_<sqmesh::base::StatusCode>(base_module, "StatusCode")
    .value("ok", sqmesh::base::StatusCode::ok)
    .value("invalid_argument", sqmesh::base::StatusCode::invalid_argument)
    .value("not_initialized", sqmesh::base::StatusCode::not_initialized)
    .value("invalid_handle", sqmesh::base::StatusCode::invalid_handle)
    .value("owner_mismatch", sqmesh::base::StatusCode::owner_mismatch)
    .value("io_error", sqmesh::base::StatusCode::io_error)
    .value("unsupported", sqmesh::base::StatusCode::unsupported)
    .value("internal_error", sqmesh::base::StatusCode::internal_error);

  base_module.def("module_name", []() {
    return std::string(sqmesh::base::module_name());
  });
  base_module.def(
    "status_code_name",
    [](sqmesh::base::StatusCode status) {
      return std::string(sqmesh::base::status_code_name(status));
    },
    py::arg("status")
  );
  base_module.def("initialize", []() {
    sqmesh::base::ContextHandle context_handle = sqmesh::invalid_handle;
    throw_on_status(sqmesh::base::initialize(context_handle));
    return context_handle;
  });
  base_module.def(
    "shutdown",
    [](sqmesh::base::ContextHandle context_handle) {
      throw_on_status(sqmesh::base::shutdown(context_handle));
    },
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  base_module.def("shutdown_all", []() {
    throw_on_status(sqmesh::base::shutdown_all());
  });
  base_module.def("is_initialized", &sqmesh::base::is_initialized);
  base_module.def("current_context", &sqmesh::base::current_context);
  base_module.def(
    "current_session",
    &sqmesh::base::current_session,
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  base_module.def("last_error_code", &sqmesh::base::last_error_code);
  base_module.def("last_error_message", []() {
    return std::string(sqmesh::base::last_error_message());
  });
}
