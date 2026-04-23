# Internal target names follow the implementation layering, while the public SDK
# surface they back is `sqmesh::base`, `sqmesh::geo`, and `sqmesh::mesh`.
add_subdirectory(core)
add_subdirectory(cad/occ)
add_subdirectory(model)
add_subdirectory(sdk)

if(SQMESH_BUILD_PYTHON_BINDINGS)
  add_subdirectory(python)
endif()
