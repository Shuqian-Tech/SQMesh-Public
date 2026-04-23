set(sqmesh_cad_occ_sources adapter_stub.cpp)

if(SQMESH_ENABLE_OCC AND SQMESH_OCC_FOUND)
  set(sqmesh_cad_occ_sources adapter_occ.cpp)
endif()

add_library(sqmesh_cad_occ STATIC ${sqmesh_cad_occ_sources})
add_library(sqmesh::cad_occ ALIAS sqmesh_cad_occ)

sqmesh_configure_target(sqmesh_cad_occ)
target_link_libraries(sqmesh_cad_occ PUBLIC sqmesh_core)

if(SQMESH_ENABLE_OCC AND SQMESH_OCC_FOUND)
  target_include_directories(sqmesh_cad_occ SYSTEM PRIVATE ${SQMESH_OCC_INCLUDE_DIRS})
  target_link_libraries(sqmesh_cad_occ PRIVATE ${SQMESH_OCC_LIBRARIES})
endif()
