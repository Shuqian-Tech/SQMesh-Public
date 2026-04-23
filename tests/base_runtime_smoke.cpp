// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "sqmesh/sqmesh.hpp"

#include <cstdio>
#include <cstdlib>
#include <string_view>
#include <type_traits>

namespace {

bool expect(bool condition, const char *message)
{
  if(condition) {
    return true;
  }

  std::fprintf(stderr, "base_runtime_smoke: %s\n", message);
  return false;
}

} // namespace

int main()
{
  static_assert(sizeof(sqmesh::Handle) == 8U, "sqmesh handles must be 64-bit.");
  static_assert(
    std::is_same_v<sqmesh::base::ContextHandle, sqmesh::Handle>,
    "ContextHandle must stay on the shared handle type."
  );
  static_assert(
    std::is_same_v<sqmesh::base::SessionHandle, sqmesh::Handle>,
    "SessionHandle must stay on the shared handle type."
  );

  if(!expect(
       sqmesh::base::status_code_name(sqmesh::base::StatusCode::ok) ==
         std::string_view("ok"),
       "status_code_name(ok) should be stable"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::base::module_name() == std::string_view("sqmesh::base"),
       "module_name() should report the sqmesh::base namespace"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::base::last_error_code() == sqmesh::base::StatusCode::ok,
       "initial last_error_code should be ok"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::base::last_error_message().empty(),
       "initial last_error_message should be empty"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::base::ContextHandle context_one = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::base::initialize(context_one) == sqmesh::base::StatusCode::ok,
       "initialize should create the first context"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       context_one != sqmesh::invalid_handle,
       "initialize should return a non-zero first context handle"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::base::is_initialized(),
       "runtime should report initialized after the first initialize"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::base::current_context() == context_one,
       "first initialize should bind the first context to the calling thread"
     )) {
    return EXIT_FAILURE;
  }

  const auto session_one = sqmesh::base::current_session(context_one);
  if(!expect(
       session_one != sqmesh::invalid_handle,
       "current_session should resolve the first context's default session"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::base::ContextHandle context_two = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::base::initialize(context_two) == sqmesh::base::StatusCode::ok,
       "initialize should create a second context"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       context_two != sqmesh::invalid_handle && context_two != context_one,
       "second initialize should return a distinct context handle"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::base::current_context() == context_two,
       "latest initialize on the thread should become current"
     )) {
    return EXIT_FAILURE;
  }

  const auto session_two = sqmesh::base::current_session(context_two);
  if(!expect(
       session_two != sqmesh::invalid_handle && session_two != session_one,
       "each context should own a distinct default session"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::base::current_session(context_one) == session_one,
       "explicit current_session lookup should keep older contexts reachable"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::base::last_error_code() == sqmesh::base::StatusCode::ok,
       "successful current_session lookups should clear last_error_code"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::base::last_error_message().empty(),
       "successful current_session lookups should clear last_error_message"
     )) {
    return EXIT_FAILURE;
  }

  if(!expect(
       sqmesh::base::shutdown(context_one) == sqmesh::base::StatusCode::ok,
       "shutdown should destroy an explicit context"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::base::current_context() == context_two,
       "destroying a non-current context should not disturb the thread current context"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::base::current_session(context_one) == sqmesh::invalid_handle,
       "destroyed contexts should no longer resolve to a session"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::base::last_error_code() == sqmesh::base::StatusCode::invalid_handle,
       "destroyed contexts should report invalid_handle"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       !sqmesh::base::last_error_message().empty(),
       "destroyed contexts should publish a diagnostic message"
     )) {
    return EXIT_FAILURE;
  }

  if(!expect(
       sqmesh::base::shutdown() == sqmesh::base::StatusCode::ok,
       "shutdown() without arguments should destroy the current context"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::base::current_context() == sqmesh::invalid_handle,
       "implicit shutdown should clear the thread current context"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       !sqmesh::base::is_initialized(),
       "implicit shutdown should leave the runtime uninitialized when it destroys the last context"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::base::current_session(context_two) == sqmesh::invalid_handle,
       "implicit shutdown should invalidate the destroyed current context"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::base::last_error_code() == sqmesh::base::StatusCode::invalid_handle,
       "implicit shutdown should leave stale-handle lookups reporting invalid_handle"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       !sqmesh::base::last_error_message().empty(),
       "implicit shutdown should preserve a diagnostic message for stale-handle lookups"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::base::ContextHandle context_three = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::base::initialize(context_three) == sqmesh::base::StatusCode::ok,
       "initialize should recreate the runtime after implicit shutdown"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       context_three != sqmesh::invalid_handle,
       "reinitialize should return a non-zero context handle"
     )) {
    return EXIT_FAILURE;
  }

  if(!expect(
       sqmesh::base::shutdown_all() == sqmesh::base::StatusCode::ok,
       "shutdown_all should tear down the remaining runtime state"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::base::current_context() == sqmesh::invalid_handle,
       "shutdown_all should clear the thread current context"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       !sqmesh::base::is_initialized(),
       "shutdown_all should leave the runtime uninitialized"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::base::current_session(context_three) == sqmesh::invalid_handle,
       "stale session lookup after shutdown_all should fail"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::base::last_error_code() == sqmesh::base::StatusCode::invalid_handle,
       "stale session lookup after shutdown_all should report invalid_handle"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       !sqmesh::base::last_error_message().empty(),
       "stale session lookup after shutdown_all should preserve an error message"
     )) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
