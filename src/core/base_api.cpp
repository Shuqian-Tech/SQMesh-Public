// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "runtime_registry.hpp"

#include "sqmesh/base/api.hpp"
#include "sqmesh/version.hpp"

#include <string_view>

namespace sqmesh {

const char *version_string() noexcept
{
  return SQMESH_VERSION_STRING;
}

} // namespace sqmesh

namespace sqmesh::base {

std::string_view module_name() noexcept
{
  return "sqmesh::base";
}

const char *status_code_name(StatusCode code) noexcept
{
  return core::detail::status_code_name(code);
}

StatusCode initialize(ContextHandle &context_handle) noexcept
{
  return core::detail::initialize(context_handle);
}

StatusCode shutdown(ContextHandle context_handle) noexcept
{
  return core::detail::shutdown(context_handle);
}

StatusCode shutdown_all() noexcept
{
  return core::detail::shutdown_all();
}

bool is_initialized() noexcept
{
  return core::detail::is_initialized();
}

ContextHandle current_context() noexcept
{
  return core::detail::current_context();
}

SessionHandle current_session(ContextHandle context_handle) noexcept
{
  return core::detail::current_session(context_handle);
}

StatusCode last_error_code() noexcept
{
  return core::detail::last_error_code();
}

std::string_view last_error_message() noexcept
{
  return core::detail::last_error_message();
}

} // namespace sqmesh::base
