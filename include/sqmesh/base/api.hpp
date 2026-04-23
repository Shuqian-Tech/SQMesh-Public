// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>

namespace sqmesh {

using Handle = std::uint64_t;
inline constexpr Handle invalid_handle = 0;

} // namespace sqmesh

namespace sqmesh::base {

using ContextHandle = sqmesh::Handle;
using SessionHandle = sqmesh::Handle;

enum class HandleKind : std::uint8_t {
  invalid = 0,
  context = 1,
  session = 2,
  model = 3,
  mesh = 4,
};

enum class StatusCode : std::uint32_t {
  ok = 0,
  invalid_argument = 1,
  not_initialized = 2,
  invalid_handle = 3,
  owner_mismatch = 4,
  io_error = 5,
  unsupported = 6,
  internal_error = 7,
};

[[nodiscard]] std::string_view module_name() noexcept;
[[nodiscard]] const char *status_code_name(StatusCode code) noexcept;

// Create a context and bind it as current for the calling thread.
[[nodiscard]] StatusCode initialize(ContextHandle &context_handle) noexcept;
// Destroy a context; `invalid_handle` targets the calling thread's current context.
[[nodiscard]] StatusCode shutdown(
  ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
// Tear down every registered context and session in the runtime; other threads
// lazily observe their bound current contexts becoming invalid on the next API call.
[[nodiscard]] StatusCode shutdown_all() noexcept;
// Return true if the process-wide runtime still owns any active contexts.
[[nodiscard]] bool is_initialized() noexcept;
// Return the context currently bound to the calling thread, if any.
[[nodiscard]] ContextHandle current_context() noexcept;
[[nodiscard]] SessionHandle current_session(
  ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] StatusCode last_error_code() noexcept;
// The returned view is valid until the next sqmesh::base call on the same thread.
[[nodiscard]] std::string_view last_error_message() noexcept;

} // namespace sqmesh::base
