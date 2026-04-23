set(SQMESH_PYTHON_PACKAGE_ROOT "${PROJECT_BINARY_DIR}/python")
set(SQMESH_PYTHON_PACKAGE_DIR "${SQMESH_PYTHON_PACKAGE_ROOT}/sqmesh")

file(MAKE_DIRECTORY "${SQMESH_PYTHON_PACKAGE_DIR}")

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG QUIET)

if(NOT pybind11_FOUND)
  include(FetchContent)
  FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG v2.13.6
    GIT_SHALLOW TRUE
  )
  FetchContent_MakeAvailable(pybind11)
endif()

set(sqmesh_python_sources
  bindings.hpp
  exceptions.hpp
  module.cpp
  exception_bindings.cpp
  base_bindings.cpp
  geo_bindings.cpp
  mesh_bindings.cpp
)

pybind11_add_module(_sqmesh MODULE ${sqmesh_python_sources})
sqmesh_configure_target(_sqmesh)
get_filename_component(
  SQMESH_PYTHON_INCLUDE_PARENT
  "${Python3_INCLUDE_DIR}"
  DIRECTORY
)
target_include_directories(
  _sqmesh
  PRIVATE
    ${Python3_INCLUDE_DIRS}
    "${SQMESH_PYTHON_INCLUDE_PARENT}"
)
target_link_libraries(_sqmesh PRIVATE sqmesh_sdk Python3::Module)

if(SQMESH_ENABLE_OCC AND SQMESH_OCC_FOUND)
  target_include_directories(_sqmesh SYSTEM PRIVATE ${SQMESH_OCC_INCLUDE_DIRS})
  target_link_libraries(_sqmesh PRIVATE ${SQMESH_OCC_LIBRARIES})
endif()

if(SQMESH_ENABLE_CGNS AND SQMESH_CGNS_FOUND)
  target_include_directories(_sqmesh SYSTEM PRIVATE ${SQMESH_CGNS_INCLUDE_DIRS})
  target_link_libraries(_sqmesh PRIVATE ${SQMESH_CGNS_LIBRARIES})
endif()

set_target_properties(
  _sqmesh
  PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${SQMESH_PYTHON_PACKAGE_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${SQMESH_PYTHON_PACKAGE_DIR}"
    ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/python"
)

foreach(_python_file IN ITEMS __init__.py base.py geo.py mesh.py)
  configure_file(
    "${PROJECT_SOURCE_DIR}/python/sqmesh/${_python_file}"
    "${SQMESH_PYTHON_PACKAGE_DIR}/${_python_file}"
    COPYONLY
  )
endforeach()

set(SQMESH_PYTHON_PACKAGE_ROOT "${SQMESH_PYTHON_PACKAGE_ROOT}" CACHE INTERNAL "")
set(SQMESH_PYTHON_EXECUTABLE "${Python3_EXECUTABLE}" CACHE INTERNAL "")
