add_executable(volume_mesh_example volume_mesh_example.cpp)
target_link_libraries(volume_mesh_example PRIVATE sqmesh_sdk)
if(SQMESH_SPDLOG_FOUND)
  target_compile_definitions(volume_mesh_example PRIVATE
    SQMESH_HAS_SPDLOG=1
    SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_DEBUG
  )
  target_link_libraries(volume_mesh_example PRIVATE spdlog::spdlog)
endif()

add_executable(surface_mesh_example surface_mesh_example.cpp)
target_link_libraries(surface_mesh_example PRIVATE sqmesh_sdk)
if(SQMESH_SPDLOG_FOUND)
  target_compile_definitions(surface_mesh_example PRIVATE
    SQMESH_HAS_SPDLOG=1
    SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_DEBUG
  )
  target_link_libraries(surface_mesh_example PRIVATE spdlog::spdlog)
endif()
