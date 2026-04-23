function(sqmesh_add_smoke_test target_name source_file)
  add_executable(${target_name} ${source_file})
  target_link_libraries(
    ${target_name}
    PRIVATE
      sqmesh_sdk
      sqmesh_model
      sqmesh_core
      sqmesh_sdk
  )
  add_test(NAME ${target_name} COMMAND $<TARGET_FILE:${target_name}>)
endfunction()

function(sqmesh_enable_warnings_as_errors target_name)
  if(MSVC)
    target_compile_options(${target_name} PRIVATE /WX /utf-8)
  else()
    target_compile_options(${target_name} PRIVATE -Werror)
  endif()
endfunction()

function(sqmesh_add_python_test test_name script_name)
  add_test(
    NAME ${test_name}
    COMMAND "${SQMESH_PYTHON_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/python/${script_name}"
  )
  set_tests_properties(
    ${test_name}
    PROPERTIES
      WORKING_DIRECTORY "${PROJECT_BINARY_DIR}"
      ENVIRONMENT "PYTHONPATH=${SQMESH_PYTHON_PACKAGE_ROOT}"
  )
endfunction()

sqmesh_add_smoke_test(sqmesh_header_base_smoke header_base_smoke.cpp)
sqmesh_enable_warnings_as_errors(sqmesh_header_base_smoke)

sqmesh_add_smoke_test(sqmesh_header_geo_smoke header_geo_smoke.cpp)
sqmesh_enable_warnings_as_errors(sqmesh_header_geo_smoke)

sqmesh_add_smoke_test(sqmesh_header_mesh_smoke header_mesh_smoke.cpp)
sqmesh_enable_warnings_as_errors(sqmesh_header_mesh_smoke)

sqmesh_add_smoke_test(sqmesh_header_smoke header_smoke.cpp)
sqmesh_enable_warnings_as_errors(sqmesh_header_smoke)

sqmesh_add_smoke_test(sqmesh_base_runtime_smoke base_runtime_smoke.cpp)
sqmesh_enable_warnings_as_errors(sqmesh_base_runtime_smoke)
sqmesh_add_smoke_test(sqmesh_mesh_core_smoke mesh_core_smoke.cpp)
sqmesh_enable_warnings_as_errors(sqmesh_mesh_core_smoke)
sqmesh_add_smoke_test(sqmesh_mesh_quality_smoke mesh_quality_smoke.cpp)
sqmesh_enable_warnings_as_errors(sqmesh_mesh_quality_smoke)
sqmesh_add_smoke_test(sqmesh_mesh_io_smoke mesh_io_smoke.cpp)
sqmesh_enable_warnings_as_errors(sqmesh_mesh_io_smoke)
target_include_directories(
  sqmesh_mesh_io_smoke
  PRIVATE ${PROJECT_SOURCE_DIR}/src/sdk
)
sqmesh_add_smoke_test(sqmesh_dummy_mesher_smoke dummy_mesher_smoke.cpp)
sqmesh_enable_warnings_as_errors(sqmesh_dummy_mesher_smoke)
sqmesh_add_smoke_test(sqmesh_bottom_up_pipeline_smoke bottom_up_pipeline_smoke.cpp)
sqmesh_enable_warnings_as_errors(sqmesh_bottom_up_pipeline_smoke)
sqmesh_add_smoke_test(sqmesh_auto_cfd_surface_smoke auto_cfd_surface_smoke.cpp)
sqmesh_enable_warnings_as_errors(sqmesh_auto_cfd_surface_smoke)
target_include_directories(
  sqmesh_auto_cfd_surface_smoke
  PRIVATE ${PROJECT_SOURCE_DIR}/src/sdk
)
sqmesh_add_smoke_test(sqmesh_surface_mesher_smoke surface_mesher_smoke.cpp)
sqmesh_enable_warnings_as_errors(sqmesh_surface_mesher_smoke)
sqmesh_add_smoke_test(sqmesh_volume_mesher_smoke volume_mesher_smoke.cpp)
sqmesh_enable_warnings_as_errors(sqmesh_volume_mesher_smoke)
sqmesh_add_smoke_test(sqmesh_error_paths_smoke error_paths_smoke.cpp)
sqmesh_enable_warnings_as_errors(sqmesh_error_paths_smoke)
sqmesh_add_smoke_test(sqmesh_geo_io_smoke geo_io_smoke.cpp)
sqmesh_enable_warnings_as_errors(sqmesh_geo_io_smoke)

set(_sqmesh_occ_step_sample_candidates
  "${CMAKE_SOURCE_DIR}/examples/cube/cube.step"
  "/usr/share/opencascade/data/step/screw.step"
)

foreach(_candidate IN LISTS _sqmesh_occ_step_sample_candidates)
  if(EXISTS "${_candidate}")
    target_compile_definitions(
      sqmesh_geo_io_smoke
      PRIVATE sqmesh_TEST_OCC_STEP_SAMPLE_PATH="${_candidate}"
    )
    break()
  endif()
endforeach()

sqmesh_add_smoke_test(sqmesh_geo_topo_smoke geo_topo_smoke.cpp)
sqmesh_enable_warnings_as_errors(sqmesh_geo_topo_smoke)
sqmesh_add_smoke_test(sqmesh_geo_proxy_cache_smoke geo_proxy_cache_smoke.cpp)
sqmesh_enable_warnings_as_errors(sqmesh_geo_proxy_cache_smoke)
target_link_libraries(
  sqmesh_geo_proxy_cache_smoke
  PRIVATE sqmesh_core sqmesh_model sqmesh_sdk
)
sqmesh_add_smoke_test(sqmesh_geo_topology_query_smoke geo_topology_query_smoke.cpp)
sqmesh_enable_warnings_as_errors(sqmesh_geo_topology_query_smoke)
sqmesh_add_smoke_test(sqmesh_geo_model_view_smoke geo_model_view_smoke.cpp)
sqmesh_enable_warnings_as_errors(sqmesh_geo_model_view_smoke)
sqmesh_add_smoke_test(sqmesh_geo_stl_smoke geo_stl_smoke.cpp)
sqmesh_enable_warnings_as_errors(sqmesh_geo_stl_smoke)

if(SQMESH_ENABLE_OCC AND SQMESH_OCC_FOUND)
  target_include_directories(
    sqmesh_geo_topo_smoke
    SYSTEM PRIVATE ${SQMESH_OCC_INCLUDE_DIRS}
  )
  target_link_libraries(sqmesh_geo_topo_smoke PRIVATE ${SQMESH_OCC_LIBRARIES})
  target_compile_definitions(sqmesh_geo_topo_smoke PRIVATE SQMESH_TEST_OCC_ENABLED=1)

  target_include_directories(
    sqmesh_geo_proxy_cache_smoke
    SYSTEM PRIVATE ${SQMESH_OCC_INCLUDE_DIRS}
  )
  target_link_libraries(
    sqmesh_geo_proxy_cache_smoke
    PRIVATE ${SQMESH_OCC_LIBRARIES}
  )
  target_compile_definitions(
    sqmesh_geo_proxy_cache_smoke
    PRIVATE SQMESH_TEST_OCC_ENABLED=1
  )

  target_include_directories(
    sqmesh_geo_topology_query_smoke
    SYSTEM PRIVATE ${SQMESH_OCC_INCLUDE_DIRS}
  )
  target_link_libraries(
    sqmesh_geo_topology_query_smoke
    PRIVATE ${SQMESH_OCC_LIBRARIES}
  )
  target_compile_definitions(
    sqmesh_geo_topology_query_smoke
    PRIVATE SQMESH_TEST_OCC_ENABLED=1
  )

  target_include_directories(
    sqmesh_geo_model_view_smoke
    SYSTEM PRIVATE ${SQMESH_OCC_INCLUDE_DIRS}
  )
  target_link_libraries(
    sqmesh_geo_model_view_smoke
    PRIVATE ${SQMESH_OCC_LIBRARIES}
  )
  target_compile_definitions(
    sqmesh_geo_model_view_smoke
    PRIVATE SQMESH_TEST_OCC_ENABLED=1
  )

  target_include_directories(
    sqmesh_surface_mesher_smoke
    SYSTEM PRIVATE ${SQMESH_OCC_INCLUDE_DIRS}
  )
  target_link_libraries(
    sqmesh_surface_mesher_smoke
    PRIVATE ${SQMESH_OCC_LIBRARIES}
  )
  target_compile_definitions(
    sqmesh_surface_mesher_smoke
    PRIVATE SQMESH_TEST_OCC_ENABLED=1
  )

  target_include_directories(
    sqmesh_volume_mesher_smoke
    SYSTEM PRIVATE ${SQMESH_OCC_INCLUDE_DIRS}
  )
  target_link_libraries(
    sqmesh_volume_mesher_smoke
    PRIVATE ${SQMESH_OCC_LIBRARIES}
  )
  target_compile_definitions(
    sqmesh_volume_mesher_smoke
    PRIVATE SQMESH_TEST_OCC_ENABLED=1
  )

  target_include_directories(
    sqmesh_auto_cfd_surface_smoke
    SYSTEM PRIVATE ${SQMESH_OCC_INCLUDE_DIRS}
  )
  target_link_libraries(
    sqmesh_auto_cfd_surface_smoke
    PRIVATE ${SQMESH_OCC_LIBRARIES}
  )
  target_compile_definitions(
    sqmesh_auto_cfd_surface_smoke
    PRIVATE SQMESH_TEST_OCC_ENABLED=1
  )

endif()

if(SQMESH_ENABLE_CGNS AND SQMESH_CGNS_FOUND)
  target_include_directories(
    sqmesh_mesh_io_smoke
    SYSTEM PRIVATE ${SQMESH_CGNS_INCLUDE_DIRS}
  )
  target_link_libraries(
    sqmesh_mesh_io_smoke
    PRIVATE ${SQMESH_CGNS_LIBRARIES}
  )
  target_compile_definitions(
    sqmesh_mesh_io_smoke
    PRIVATE SQMESH_TEST_CGNS_ENABLED=1
  )
endif()

if(SQMESH_BUILD_PYTHON_BINDINGS)
  sqmesh_add_python_test(sqmesh_python_base_runtime test_base_runtime.py)
  sqmesh_add_python_test(sqmesh_python_geo_queries test_geo_queries.py)
  sqmesh_add_python_test(sqmesh_python_mesh_workflow test_mesh_workflow.py)
endif()
