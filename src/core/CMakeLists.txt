add_library(sqmesh_core STATIC base_api.cpp runtime_registry.cpp)
add_library(sqmesh::core ALIAS sqmesh_core)

sqmesh_configure_target(sqmesh_core)
