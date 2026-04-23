cmake_minimum_required(VERSION 3.18)

project(
  SQMesh
  VERSION 0.1.0
  DESCRIPTION "Milestone 1 baseline for the sqmesh SDK"
  LANGUAGES CXX
)

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(SQMESH_GENERATED_INCLUDE_DIR "${PROJECT_BINARY_DIR}/generated/include")
file(MAKE_DIRECTORY "${SQMESH_GENERATED_INCLUDE_DIR}/sqmesh")
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/include/sqmesh/version.hpp.in"
  "${SQMESH_GENERATED_INCLUDE_DIR}/sqmesh/version.hpp"
  @ONLY
)

include(SQMeshOptions)
include(SQMeshDependencies)
include(SQMeshTargets)

message(STATUS "Configuring ${PROJECT_NAME} ${PROJECT_VERSION}")
message(STATUS "  SQMESH_BUILD_TESTS      = ${SQMESH_BUILD_TESTS}")
message(STATUS "  SQMESH_BUILD_BENCHMARKS = ${SQMESH_BUILD_BENCHMARKS}")
message(STATUS "  SQMESH_BUILD_PYTHON_BINDINGS = ${SQMESH_BUILD_PYTHON_BINDINGS}")
message(STATUS "  SQMESH_BUILD_LONG_BENCHMARKS = ${SQMESH_BUILD_LONG_BENCHMARKS}")
message(STATUS "  SQMESH_ENABLE_OCC       = ${SQMESH_ENABLE_OCC}")
message(STATUS "  SQMESH_OCC_FOUND        = ${SQMESH_OCC_FOUND}")
message(STATUS "  SQMESH_OCC_PACKAGE      = ${SQMESH_OCC_PACKAGE}")
message(STATUS "  SQMESH_OCC_VERSION      = ${SQMESH_OCC_VERSION}")
message(STATUS "  SQMESH_ENABLE_CGNS      = ${SQMESH_ENABLE_CGNS}")
message(STATUS "  SQMESH_CGNS_FOUND       = ${SQMESH_CGNS_FOUND}")
message(STATUS "  SQMESH_CGNS_PACKAGE     = ${SQMESH_CGNS_PACKAGE}")
message(STATUS "  SQMESH_CGNS_VERSION     = ${SQMESH_CGNS_VERSION}")

add_subdirectory(src)
add_subdirectory(examples)
if(SQMESH_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

if(SQMESH_BUILD_BENCHMARKS)
  add_subdirectory(benchmarks)
endif()
