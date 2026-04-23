include_guard(GLOBAL)

option(SQMESH_BUILD_TESTS "Build smoke tests for the baseline skeleton" ON)
option(SQMESH_BUILD_BENCHMARKS "Build benchmark placeholders" OFF)
option(
  SQMESH_BUILD_PYTHON_BINDINGS
  "Build the bounded Phase 2 Python bindings and package facade"
  OFF
)
option(
  SQMESH_BUILD_LONG_BENCHMARKS
  "Register long-running benchmark gates such as the 1M-cell Milestone 1 check"
  OFF
)
option(SQMESH_ENABLE_OCC "Build the optional OCC adapter boundary" OFF)
option(SQMESH_ENABLE_CGNS "Build the optional CGNS mesh IO boundary" OFF)
