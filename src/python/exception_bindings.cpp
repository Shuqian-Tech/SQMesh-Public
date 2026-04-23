// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "bindings.hpp"
#include "exceptions.hpp"

#include "sqmesh/base/api.hpp"

namespace {

std::string last_error_message_or_name(sqmesh::base::StatusCode status)
{
  const auto message = sqmesh::base::last_error_message();
  if(!message.empty()) {
    return std::string(message);
  }

  const auto *name = sqmesh::base::status_code_name(status);
  if(name != nullptr && name[0] != '\0') {
    return std::string(name);
  }

  return "sqmesh error";
}

} // namespace

void register_sqmesh_exceptions(py::module_ &module)
{
  auto sqmesh_error =
    py::register_exception<SQMeshError>(module, "SQMeshError");
  py::register_exception<InvalidArgumentError>(
    module,
    "InvalidArgumentError",
    sqmesh_error.ptr()
  );
  py::register_exception<NotInitializedError>(
    module,
    "NotInitializedError",
    sqmesh_error.ptr()
  );
  py::register_exception<InvalidHandleError>(
    module,
    "InvalidHandleError",
    sqmesh_error.ptr()
  );
  py::register_exception<OwnerMismatchError>(
    module,
    "OwnerMismatchError",
    sqmesh_error.ptr()
  );
  py::register_exception<IoError>(module, "IoError", sqmesh_error.ptr());
  py::register_exception<UnsupportedError>(
    module,
    "UnsupportedError",
    sqmesh_error.ptr()
  );
  py::register_exception<InternalError>(
    module,
    "InternalError",
    sqmesh_error.ptr()
  );
}

[[noreturn]] void throw_sqmesh_error(sqmesh::base::StatusCode status)
{
  const auto message = last_error_message_or_name(status);
  switch(status) {
  case sqmesh::base::StatusCode::invalid_argument:
    throw InvalidArgumentError(message);
  case sqmesh::base::StatusCode::not_initialized:
    throw NotInitializedError(message);
  case sqmesh::base::StatusCode::invalid_handle:
    throw InvalidHandleError(message);
  case sqmesh::base::StatusCode::owner_mismatch:
    throw OwnerMismatchError(message);
  case sqmesh::base::StatusCode::io_error:
    throw IoError(message);
  case sqmesh::base::StatusCode::unsupported:
    throw UnsupportedError(message);
  case sqmesh::base::StatusCode::internal_error:
    throw InternalError(message);
  case sqmesh::base::StatusCode::ok:
  default:
    throw SQMeshError(message);
  }
}

void throw_on_status(sqmesh::base::StatusCode status)
{
  if(status == sqmesh::base::StatusCode::ok) {
    return;
  }

  throw_sqmesh_error(status);
}
