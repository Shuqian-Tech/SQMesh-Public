add_library(
  sqmesh_sdk
  STATIC
    # mesh/api
    ../mesh/api/mesh_api.cpp
    ../mesh/api/mesh_sdk.cpp
    # mesh/io
    ../mesh/io/mesh_io.cpp
    # mesh/framework
    ../mesh/framework/meshing_framework.cpp
    # mesh/sizing
    ../mesh/sizing/adaptive_octree.cpp
    ../mesh/sizing/mesh_size_controls.cpp
    ../mesh/sizing/proxy_mesh_geometry.cpp
    ../mesh/sizing/size_function.cpp
    # core/predicates — Shewchuk's public-domain robust geometric predicates,
    # shared by the tet mesher and auto_cfd's CDT Delaunay path.
    ../core/predicates.cxx
    # mesh/auto_cfd
    ../mesh/auto_cfd/auto_cfd_surface_mesher.cpp
    ../mesh/auto_cfd/auto_cfd_surface_pipeline.cpp
    ../mesh/auto_cfd/auto_cfd_surface_proximity_index.cpp
    # mesh/tet/core (tet mesher engine, split from a single monolith along
    # its `//== <topic>_cxx ==//` banners; behavior unchanged)
    ../mesh/tet/core/mempool.cxx
    ../mesh/tet/core/geom.cxx
    ../mesh/tet/core/flip.cxx
    ../mesh/tet/core/delaunay.cxx
    ../mesh/tet/core/surface.cxx
    ../mesh/tet/core/constrained.cxx
    ../mesh/tet/core/steiner.cxx
    ../mesh/tet/core/reconstruct.cxx
    ../mesh/tet/core/refine.cxx
    ../mesh/tet/core/optimize.cxx
    ../mesh/tet/core/meshstat.cxx
    ../mesh/tet/core/output.cxx
    ../mesh/tet/core/run_tet_mesh_core.cxx
    # mesh/tet (SQMesh-owned tet API + adapter)
    ../mesh/tet/tet_native_algorithm.cpp
    # mesh/boundary_layer
    ../mesh/boundary_layer/boundary_layer_solver.cpp
    ../mesh/boundary_layer/boundary_layer_mesher.cpp
    # mesh/region
    ../mesh/region/region_detector.cpp
)
add_library(sqmesh::sdk ALIAS sqmesh_sdk)

sqmesh_configure_target(sqmesh_sdk)
# src/ as include root so mesh/* files can use "core/..." and "model/..." paths.
target_include_directories(sqmesh_sdk PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(sqmesh_sdk PUBLIC sqmesh_model)
target_link_libraries(sqmesh_core PUBLIC sqmesh_sdk)

if(SQMESH_SPDLOG_FOUND)
  target_compile_definitions(sqmesh_sdk PRIVATE
    SQMESH_HAS_SPDLOG=1
    SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_DEBUG
  )
  target_link_libraries(sqmesh_sdk PRIVATE spdlog::spdlog)
endif()

if(SQMESH_ENABLE_CGNS AND SQMESH_CGNS_FOUND)
  target_compile_definitions(sqmesh_sdk PRIVATE SQMESH_HAVE_CGNS=1)
  target_include_directories(sqmesh_sdk SYSTEM PRIVATE ${SQMESH_CGNS_INCLUDE_DIRS})
  target_link_libraries(sqmesh_sdk PRIVATE ${SQMESH_CGNS_LIBRARIES})
endif()
