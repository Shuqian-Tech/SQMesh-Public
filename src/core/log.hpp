// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#pragma once

// SQMesh logging facade.
//
// Usage:
//   #include "core/log.hpp"
//
//   // Log messages (fmt syntax):
//   SQMESH_LOG_INFO("imported {} faces", count);
//   SQMESH_LOG_DEBUG("vertex ({}, {}, {})", x, y, z);
//   SQMESH_LOG_WARN("degenerate triangle at face {}", id);
//   SQMESH_LOG_ERROR("failed to open {}", path);
//
//   // Configure at startup:
//   sqmesh::log::init();                              // console only, level=info
//   sqmesh::log::init("sqmesh.log");                  // console + file
//   sqmesh::log::init("sqmesh.log", "debug");         // console + file, level=debug
//   sqmesh::log::init("", "trace");                   // console only, level=trace
//
//   // Change level at runtime:
//   sqmesh::log::set_level("debug");
//
// Log levels (low to high): trace, debug, info, warn, error, critical, off.

#ifdef SQMESH_HAS_SPDLOG

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>

#include <memory>
#include <string>
#include <string_view>
#include <vector>

#define SQMESH_LOG_TRACE(...)    SPDLOG_TRACE(__VA_ARGS__)
#define SQMESH_LOG_DEBUG(...)    SPDLOG_DEBUG(__VA_ARGS__)
#define SQMESH_LOG_INFO(...)     SPDLOG_INFO(__VA_ARGS__)
#define SQMESH_LOG_WARN(...)     SPDLOG_WARN(__VA_ARGS__)
#define SQMESH_LOG_ERROR(...)    SPDLOG_ERROR(__VA_ARGS__)
#define SQMESH_LOG_CRITICAL(...) SPDLOG_CRITICAL(__VA_ARGS__)

namespace sqmesh::log {

inline spdlog::level::level_enum parse_level(std::string_view name) {
  if(name == "trace")    return spdlog::level::trace;
  if(name == "debug")    return spdlog::level::debug;
  if(name == "info")     return spdlog::level::info;
  if(name == "warn")     return spdlog::level::warn;
  if(name == "error")    return spdlog::level::err;
  if(name == "critical") return spdlog::level::critical;
  if(name == "off")      return spdlog::level::off;
  return spdlog::level::info;
}

// Initialize the logging system.
//   file_path: if non-empty, also write to this file.
//   level:     "trace", "debug", "info", "warn", "error", "critical", "off".
inline void init(
  std::string_view file_path = "",
  std::string_view level = "info"
) {
  std::vector<spdlog::sink_ptr> sinks;

  // Console sink (colored).
  auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
  sinks.push_back(console_sink);

  // File sink (optional).
  if(!file_path.empty()) {
    auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(
      std::string(file_path), /*truncate=*/true
    );
    sinks.push_back(file_sink);
  }

  auto logger = std::make_shared<spdlog::logger>(
    "sqmesh", sinks.begin(), sinks.end()
  );
  logger->set_level(parse_level(level));
  logger->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] %v");
  spdlog::set_default_logger(logger);
}

// Change log level at runtime.
inline void set_level(std::string_view level) {
  spdlog::set_level(parse_level(level));
}

// Flush all pending log messages.
inline void flush() {
  spdlog::default_logger()->flush();
}

} // namespace sqmesh::log

#else // fallback: no spdlog

#include <cstdio>
#include <string_view>

#define SQMESH_LOG_TRACE(...)    ((void)0)
#define SQMESH_LOG_DEBUG(...)    ((void)0)
#define SQMESH_LOG_INFO(fmt, ...)  std::fprintf(stderr, "[info] " fmt "\n", ##__VA_ARGS__)
#define SQMESH_LOG_WARN(fmt, ...)  std::fprintf(stderr, "[warn] " fmt "\n", ##__VA_ARGS__)
#define SQMESH_LOG_ERROR(fmt, ...) std::fprintf(stderr, "[error] " fmt "\n", ##__VA_ARGS__)
#define SQMESH_LOG_CRITICAL(fmt, ...) std::fprintf(stderr, "[critical] " fmt "\n", ##__VA_ARGS__)

namespace sqmesh::log {
inline void init(std::string_view = "", std::string_view = "") {}
inline void set_level(std::string_view) {}
inline void flush() {}
} // namespace sqmesh::log

#endif
