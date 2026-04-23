// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "sqmesh/sqmesh.hpp"

#include <cstdio>
#include <cstdlib>

#if defined(SQMESH_TEST_OCC_ENABLED)
#include "../src/cad/occ/adapter.hpp"

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

  std::fprintf(stderr, "geo_topo_smoke: %s\n", message);
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
    sqmesh::geo::ModelHandle model_handle = sqmesh::invalid_handle;
    if(!expect(
         sqmesh::geo::create_placeholder_model(model_handle, context) ==
           sqmesh::base::StatusCode::ok,
         "placeholder geometry should still allocate a model handle"
       )) {
      return EXIT_FAILURE;
    }

    std::size_t free_edges = 1U;
    if(!expect(
         sqmesh::geo::free_edge_count(model_handle, free_edges, context) ==
           sqmesh::base::StatusCode::unsupported,
         "free_edge_count should report unsupported when OCC is disabled"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::geo::TopoReport report;
    if(!expect(
         sqmesh::geo::topo(model_handle, report, {}, context) ==
           sqmesh::base::StatusCode::unsupported,
         "topo should report unsupported when OCC is disabled"
       )) {
      return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
  }

#if !defined(SQMESH_TEST_OCC_ENABLED)
  std::fprintf(
    stderr,
    "geo_topo_smoke: OCC runtime is enabled but the smoke target was not built with OCC.\n"
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

  std::size_t free_edges_before = 0U;
  if(!expect(
       sqmesh::geo::free_edge_count(model_handle, free_edges_before, context) ==
         sqmesh::base::StatusCode::ok,
       "free_edge_count should inspect the bad OCC test shape"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       free_edges_before > 0U,
       "the synthetic bad model should expose free edges before repair"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::TopologyCheckReport check_before;
  if(!expect(
       sqmesh::geo::check_topology(model_handle, check_before, context) ==
         sqmesh::base::StatusCode::ok,
       "check_topology should report the pre-repair topology"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       check_before.free_edge_count == free_edges_before,
       "check_topology and free_edge_count should agree before repair"
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
  if(!expect(
       topology_before.topology_revision > 0U,
       "OCC-backed models should expose a non-zero topology revision"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::TopoOptions topo_options;
  topo_options.tolerance = 1.0e-6;

  sqmesh::geo::TopoReport topo_report;
  if(!expect(
       sqmesh::geo::topo(model_handle, topo_report, topo_options, context) ==
         sqmesh::base::StatusCode::ok,
       "topo should repair the synthetic bad OCC shape"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       topo_report.before.free_edge_count == free_edges_before,
       "topo should report the same free-edge count it observed before repair"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       topo_report.after.free_edge_count < topo_report.before.free_edge_count,
       "topo should reduce the number of free edges"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       topo_report.after.free_edge_count == 0U,
       "topo should stitch the synthetic bad model closed"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       topo_report.modified && topo_report.free_edges_reduced,
       "topo should report that it modified the model and reduced free edges"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       topo_report.topology_revision_after > 0U &&
         topo_report.topology_identity_changed,
       "topo should report that repair advanced the topology revision and invalidated entity ids"
     )) {
    return EXIT_FAILURE;
  }

  std::size_t free_edges_after = free_edges_before;
  if(!expect(
       sqmesh::geo::free_edge_count(model_handle, free_edges_after, context) ==
         sqmesh::base::StatusCode::ok,
       "free_edge_count should inspect the repaired OCC test shape"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       free_edges_after == topo_report.after.free_edge_count,
       "free_edge_count and topo should agree after repair"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       free_edges_after == 0U,
       "the repaired OCC test shape should not retain free edges"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::TopologySnapshot topology_after;
  if(!expect(
       sqmesh::geo::topology_snapshot(model_handle, topology_after, context) ==
         sqmesh::base::StatusCode::ok &&
         topology_after.topology_revision == topo_report.topology_revision_after,
       "topology_snapshot should expose the advanced topology revision after repair"
     )) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
#endif
}
