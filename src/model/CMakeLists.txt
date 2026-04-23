add_library(
  sqmesh_model
  STATIC
    geo_api.cpp
    stl_import.cpp
)
add_library(sqmesh::model ALIAS sqmesh_model)

sqmesh_configure_target(sqmesh_model)
target_link_libraries(sqmesh_model PUBLIC sqmesh_core sqmesh_cad_occ)
