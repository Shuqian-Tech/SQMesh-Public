// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "sqmesh/sqmesh.hpp"

#include <cstdio>
#include <cstdlib>
#include <unordered_set>

#if defined(SQMESH_TEST_OCC_ENABLED)
#include "../src/cad/occ/adapter.hpp"
#include "../src/core/runtime_registry.hpp"
#include "../src/model/geometry_model_storage.hpp"

#include <BRep_Builder.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Shape.hxx>
#endif

namespace {

bool expect(bool condition, const char *message)
{
  if(condition) {
    return true;
  }

  std::fprintf(stderr, "geo_proxy_cache_smoke: %s\n", message);
  return false;
}

bool expect_with_last_error(bool condition, const char *message)
{
  if(condition) {
    return true;
  }

  std::fprintf(
    stderr,
    "geo_proxy_cache_smoke: %s (last_error=%.*s)\n",
    message,
    static_cast<int>(sqmesh::base::last_error_message().size()),
    sqmesh::base::last_error_message().data()
  );
  return false;
}

#if defined(SQMESH_TEST_OCC_ENABLED)
TopoDS_Shape make_unstitched_box_faces()
{
  const TopoDS_Shape box = BRepPrimAPI_MakeBox(10.0, 20.0, 30.0).Shape();

  BRep_Builder builder;
  TopoDS_Compound compound;
  builder.MakeCompound(compound);

  for(TopExp_Explorer face_it(box, TopAbs_FACE); face_it.More(); face_it.Next()) {
    BRepBuilderAPI_Copy face_copy(face_it.Current());
    builder.Add(compound, face_copy.Shape());
  }

  return compound;
}

bool valid_proxy_ranges(const sqmesh::model::detail::GeometryCoarseProxyMesh &proxy)
{
  std::size_t expected_offset = 0U;
  for(std::size_t face_index = 0U; face_index < proxy.face_triangle_ranges.size(); ++face_index) {
    const auto &range = proxy.face_triangle_ranges[face_index];
    if(range.offset != expected_offset) {
      return false;
    }
    if(range.offset + range.count > proxy.triangle_face_owner.size()) {
      return false;
    }

    for(std::size_t offset = 0U; offset < range.count; ++offset) {
      const auto owner = proxy.triangle_face_owner[range.offset + offset];
      if(owner.dimension != sqmesh::geo::TopologyDimension::face ||
         owner.index != face_index) {
        return false;
      }
    }
    expected_offset += range.count;
  }

  return expected_offset == proxy.triangle_face_owner.size();
}

bool valid_node_owners(const sqmesh::model::detail::GeometryCoarseProxyMesh &proxy)
{
  if(proxy.node_topology_owner.size() != proxy.nodes.size()) {
    return false;
  }

  for(const auto owner : proxy.node_topology_owner) {
    if(!sqmesh::geo::is_valid(owner)) {
      return false;
    }
    switch(owner.dimension) {
    case sqmesh::geo::TopologyDimension::vertex:
    case sqmesh::geo::TopologyDimension::edge:
    case sqmesh::geo::TopologyDimension::face:
      break;
    case sqmesh::geo::TopologyDimension::region:
      return false;
    }
  }

  return true;
}

bool valid_edge_owners(const sqmesh::model::detail::GeometryCoarseProxyMesh &proxy)
{
  if(proxy.edge_topology_owner.size() != proxy.edge_nodes.size()) {
    return false;
  }

  for(std::size_t edge_index = 0U; edge_index < proxy.edge_topology_owner.size(); ++edge_index) {
    const auto owner = proxy.edge_topology_owner[edge_index];
    if(!sqmesh::geo::is_valid(owner) ||
       owner.dimension != sqmesh::geo::TopologyDimension::edge ||
       owner.index != edge_index) {
      return false;
    }
  }

  return true;
}

bool proxy_domain_edge_owners_are_valid(
  const sqmesh::mesh::Domain &domain,
  std::size_t expected_distinct_owner_count
)
{
  std::size_t proxy_edge_count = 0U;
  std::unordered_set<std::uint32_t> distinct_owner_ids;

  for(const auto &entity_group : domain.entity_groups()) {
    if(entity_group.role() != sqmesh::mesh::EntityGroupRole::geometric_proxy ||
       entity_group.order() != sqmesh::mesh::EntityOrder::edge) {
      continue;
    }

    for(std::uint32_t edge_index = 0U; edge_index < entity_group.edges().size(); ++edge_index) {
      ++proxy_edge_count;
      const auto owner = domain.edge_topology_owner({entity_group.id(), edge_index});
      if(!sqmesh::geo::is_valid(owner) ||
         owner.dimension != sqmesh::geo::TopologyDimension::edge) {
        return false;
      }
      distinct_owner_ids.insert(owner.index);
    }
  }

  return proxy_edge_count > 0U && distinct_owner_ids.size() == expected_distinct_owner_count;
}

bool capture_proxy_mesh(
  sqmesh::geo::ModelHandle model_handle,
  sqmesh::base::ContextHandle context,
  sqmesh::model::detail::GeometryCoarseProxyMesh &proxy
)
{
  proxy = {};
  return sqmesh::core::detail::with_model_storage(
           model_handle,
           [&](const sqmesh::model::detail::GeometryModelStorage &storage) {
             const auto status = storage.coarse_proxy_mesh(proxy);
             if(status != sqmesh::base::StatusCode::ok) {
               return status;
             }
             return sqmesh::base::StatusCode::ok;
           },
           context
         ) == sqmesh::base::StatusCode::ok;
}
#endif

} // namespace

int main()
{
  sqmesh::base::ContextHandle context = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::base::initialize(context) == sqmesh::base::StatusCode::ok,
       "initialize should create a runtime context"
     )) {
    return EXIT_FAILURE;
  }

  if(!sqmesh::geo::cad_io_available()) {
    return EXIT_SUCCESS;
  }

#if !defined(SQMESH_TEST_OCC_ENABLED)
  std::fprintf(
    stderr,
    "geo_proxy_cache_smoke: OCC runtime is enabled but the smoke target was not built with OCC.\n"
  );
  return EXIT_FAILURE;
#else
  sqmesh::geo::ModelHandle model_handle = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::cad::occ::store_shape(make_unstitched_box_faces(), model_handle, context) ==
         sqmesh::base::StatusCode::ok,
       "the test should be able to seed an OCC-backed geometry model"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::TopologySnapshot topology_before;
  if(!expect(
       sqmesh::geo::topology_snapshot(model_handle, topology_before, context) ==
         sqmesh::base::StatusCode::ok,
       "topology_snapshot should expose the pre-repair topology revision"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::model::detail::GeometryCoarseProxyMesh proxy_before;
  if(!expect(
       capture_proxy_mesh(model_handle, context, proxy_before),
       "the internal OCC proxy cache should be available through the storage borrow path"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       !proxy_before.empty(),
       "the OCC proxy cache should contain coarse proxy nodes and triangles"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       proxy_before.topology_revision == topology_before.topology_revision,
       "the OCC proxy cache should be stamped with the current topology revision"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       proxy_before.source_local_node_count == proxy_before.nodes.size(),
       "the pre-repair proxy for disconnected box faces should not merge face-local nodes across faces"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       valid_node_owners(proxy_before),
       "the OCC proxy cache should preserve a valid typed topology owner per global proxy node"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       valid_edge_owners(proxy_before),
       "the OCC proxy cache should preserve a valid typed topology owner per feature edge"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       proxy_before.triangle_face_owner.size() == proxy_before.triangles.size() &&
         proxy_before.face_triangle_ranges.size() == topology_before.faces.size() &&
         valid_proxy_ranges(proxy_before),
       "the OCC proxy cache should preserve per-triangle source face ownership and consistent face ranges"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle proxy_mesh_handle = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::geo::model_proxy_mesh(model_handle, proxy_mesh_handle, context) ==
         sqmesh::base::StatusCode::ok,
       "model_proxy_mesh should expose the shared proxy-backed mesh handle"
     )) {
    return EXIT_FAILURE;
  }
  sqmesh::mesh::MeshSummary proxy_mesh_summary;
  sqmesh::mesh::Domain proxy_domain;
  if(!expect(
       sqmesh::mesh::mesh_summary(proxy_mesh_handle, proxy_mesh_summary, context) ==
         sqmesh::base::StatusCode::ok &&
       sqmesh::mesh::domain_snapshot(proxy_mesh_handle, proxy_domain, context) ==
         sqmesh::base::StatusCode::ok,
       "proxy mesh handles should resolve through the public mesh summary and snapshot APIs"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       proxy_mesh_summary.node_count == 0U &&
         proxy_mesh_summary.edge_count == 0U &&
         proxy_mesh_summary.face_count == 0U &&
         proxy_mesh_summary.cell_count == 0U &&
         proxy_mesh_summary.source_topology_revision == topology_before.topology_revision,
       "proxy-only mesh summaries should report zero computational entities while preserving the source topology revision"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       proxy_domain.entity_group_count(sqmesh::mesh::EntityGroupRole::computational) == 0U &&
         proxy_domain.entity_group_count(sqmesh::mesh::EntityGroupRole::geometric_proxy) == 3U &&
         proxy_domain.source_topology_revision() == topology_before.topology_revision,
       "domain snapshots should still expose the retained geometric proxy entity_groups on a proxy-only mesh handle"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       proxy_domain_edge_owners_are_valid(proxy_domain, proxy_before.edge_nodes.size()),
       "proxy mesh snapshots should preserve per-edge topology ownership for geometric proxy segments"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::ParameterDictionary surface_parameters;
  surface_parameters.set_number("minimum_length", 4.0);
  surface_parameters.set_number("maximum_length", 8.0);
  surface_parameters.set_number("distortion_angle", 18.0);
  surface_parameters.set_number("growth_rate", 1.2);
  surface_parameters.set_text("element_type", "tri");
  sqmesh::mesh::MeshHandle meshed_proxy_handle = sqmesh::invalid_handle;
  if(!expect_with_last_error(
       sqmesh::mesh::create_surface_mesh(
         model_handle,
         "Auto CFD Surface Mesher",
         surface_parameters,
         meshed_proxy_handle,
         context
       ) == sqmesh::base::StatusCode::ok,
       "surface meshing should be able to append computational entity_groups onto the proxy-backed mesh handle"
     )) {
    return EXIT_FAILURE;
  }
  sqmesh::mesh::MeshSummary meshed_proxy_summary;
  sqmesh::mesh::Domain meshed_proxy_domain;
  if(!expect(
       meshed_proxy_handle == proxy_mesh_handle &&
         sqmesh::mesh::mesh_summary(meshed_proxy_handle, meshed_proxy_summary, context) ==
           sqmesh::base::StatusCode::ok &&
         sqmesh::mesh::domain_snapshot(meshed_proxy_handle, meshed_proxy_domain, context) ==
           sqmesh::base::StatusCode::ok,
       "surface meshing should keep using the stable proxy-backed mesh handle"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       meshed_proxy_summary.face_count > 0U &&
         meshed_proxy_summary.source_topology_revision == topology_before.topology_revision &&
         meshed_proxy_domain.entity_group_count(sqmesh::mesh::EntityGroupRole::computational) == 3U &&
         meshed_proxy_domain.entity_group_count(sqmesh::mesh::EntityGroupRole::geometric_proxy) == 3U,
       "surface meshing should append computational entity_groups while preserving the proxy entity_groups before topology repair"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::TopoOptions topo_options;
  topo_options.tolerance = 1.0e-6;
  sqmesh::geo::TopoReport topo_report;
  if(!expect(
       sqmesh::geo::topo(model_handle, topo_report, topo_options, context) ==
         sqmesh::base::StatusCode::ok,
       "topo should stitch the synthetic bad model before reusing the OCC proxy cache"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       topo_report.topology_revision_after > topo_report.topology_revision_before,
       "topo should advance the topology revision after stitching the synthetic bad model"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::TopologySnapshot topology_after;
  if(!expect(
       sqmesh::geo::topology_snapshot(model_handle, topology_after, context) ==
         sqmesh::base::StatusCode::ok,
       "topology_snapshot should expose the post-repair topology revision"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::model::detail::GeometryCoarseProxyMesh proxy_after;
  if(!expect(
       capture_proxy_mesh(model_handle, context, proxy_after),
       "the internal OCC proxy cache should rebuild after topology repair"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       proxy_after.topology_revision == topology_after.topology_revision &&
         proxy_after.topology_revision > proxy_before.topology_revision,
       "the OCC proxy cache should invalidate and rebuild against the advanced topology revision"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       proxy_after.source_local_node_count > proxy_after.nodes.size(),
       "the stitched OCC proxy should merge shared per-face nodes instead of staying disconnected polygon soup"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       proxy_after.nodes.size() < proxy_before.nodes.size(),
       "the stitched OCC proxy should compact shared global proxy nodes across adjacent faces"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       valid_node_owners(proxy_after),
       "the rebuilt OCC proxy cache should preserve a valid typed topology owner per global proxy node"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       valid_edge_owners(proxy_after),
       "the rebuilt OCC proxy cache should preserve a valid typed topology owner per feature edge"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       proxy_after.triangle_face_owner.size() == proxy_after.triangles.size() &&
         proxy_after.face_triangle_ranges.size() == topology_after.faces.size() &&
         valid_proxy_ranges(proxy_after),
       "the rebuilt OCC proxy cache should preserve per-triangle source face ownership and consistent face ranges"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshSummary proxy_mesh_summary_after;
  sqmesh::mesh::Domain proxy_domain_after;
  if(!expect(
       sqmesh::mesh::mesh_summary(proxy_mesh_handle, proxy_mesh_summary_after, context) ==
         sqmesh::base::StatusCode::ok &&
       sqmesh::mesh::domain_snapshot(proxy_mesh_handle, proxy_domain_after, context) ==
         sqmesh::base::StatusCode::ok,
       "the original proxy mesh handle should auto-refresh through mesh summary and snapshot after topology repair"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       proxy_mesh_summary_after.source_topology_revision == topology_after.topology_revision &&
         proxy_domain_after.source_topology_revision() == topology_after.topology_revision,
       "proxy mesh handle refresh should advance the runtime proxy mesh revision to the repaired topology revision"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       proxy_mesh_summary_after.node_count == 0U &&
         proxy_mesh_summary_after.edge_count == 0U &&
         proxy_mesh_summary_after.face_count == 0U &&
         proxy_mesh_summary_after.cell_count == 0U &&
         proxy_domain_after.entity_group_count(sqmesh::mesh::EntityGroupRole::computational) == 0U &&
         proxy_domain_after.entity_group_count(sqmesh::mesh::EntityGroupRole::geometric_proxy) == 3U,
       "proxy handle refresh should preserve the proxy-only summary semantics while rebuilding the retained proxy entity_groups"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       proxy_domain_edge_owners_are_valid(proxy_domain_after, proxy_after.edge_nodes.size()),
       "runtime proxy mesh refresh should rebuild proxy edge ownership onto the retained handle"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle proxy_mesh_handle_again = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::geo::model_proxy_mesh(model_handle, proxy_mesh_handle_again, context) ==
         sqmesh::base::StatusCode::ok &&
         proxy_mesh_handle_again == proxy_mesh_handle,
       "model_proxy_mesh should keep returning the same stable proxy handle after refresh"
     )) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
#endif
}
