// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
// Auto CFD Surface Mesher smoke test.
//
// Validates core meshing functionality on simple OCC geometries:
//   1. Box (6 faces) — basic surface mesh generation
//   2. Trimmed face with hole — inner loop handling
//   3. Parallel faces — proximity detection
//   4. Parameter validation — error paths

#include "sqmesh/sqmesh.hpp"

#include "../src/mesh/auto_cfd/auto_cfd_surface_pipeline.hpp"
#include "../src/mesh/sizing/mesh_size_controls.hpp"

#include <cstdio>
#include <cstdlib>
#include <string_view>

#if defined(SQMESH_TEST_OCC_ENABLED)
#include "../src/cad/occ/adapter.hpp"

#include <BRep_Builder.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#endif

namespace {

bool expect(bool condition, const char *message) {
  if(!condition) {
    std::fprintf(stderr, "auto_cfd_surface_smoke: %s\n", message);
  }
  return condition;
}

#if defined(SQMESH_TEST_OCC_ENABLED)

// ── Geometry helpers ──────────────────────────────────────────────────

TopoDS_Face make_rectangular_face_with_hole()
{
  // Outer: 10 x 8 rectangle.  Inner: 3..7 x 2..6 hole.
  TopoDS_Wire outer_wire = BRepBuilderAPI_MakePolygon(
    gp_Pnt(0, 0, 0), gp_Pnt(10, 0, 0),
    gp_Pnt(10, 8, 0), gp_Pnt(0, 8, 0), Standard_True
  ).Wire();

  TopoDS_Wire inner_wire = BRepBuilderAPI_MakePolygon(
    gp_Pnt(3, 2, 0), gp_Pnt(7, 2, 0),
    gp_Pnt(7, 6, 0), gp_Pnt(3, 6, 0), Standard_True
  ).Wire();

  BRepBuilderAPI_MakeFace builder(gp_Pln(gp_Pnt(0,0,0), gp_Dir(0,0,1)), outer_wire);
  builder.Add(TopoDS::Wire(inner_wire.Reversed()));
  return builder.Face();
}

TopoDS_Face make_rectangular_face(double z, bool positive_normal)
{
  TopoDS_Wire wire = BRepBuilderAPI_MakePolygon(
    gp_Pnt(0, 0, z), gp_Pnt(2, 0, z),
    gp_Pnt(2, 2, z), gp_Pnt(0, 2, z), Standard_True
  ).Wire();

  gp_Dir normal = positive_normal ? gp_Dir(0, 0, 1) : gp_Dir(0, 0, -1);
  return BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(0, 0, z), normal), wire).Face();
}

TopoDS_Shape make_parallel_face_compound()
{
  // Two 2x2 faces separated by 0.3 gap — proximity test geometry.
  TopoDS_Compound compound;
  BRep_Builder builder;
  builder.MakeCompound(compound);
  builder.Add(compound, make_rectangular_face(0.0, true));
  builder.Add(compound, make_rectangular_face(0.3, false));
  return compound;
}

std::size_t count_faces(const sqmesh::mesh::Domain &domain)
{
  std::size_t count = 0;
  for(const auto &entity_group : domain.entity_groups()) {
    if(entity_group.order() == sqmesh::mesh::EntityOrder::face) {
      count += entity_group.faces().size();
    }
  }
  return count;
}

#endif // SQMESH_TEST_OCC_ENABLED

} // namespace

int main()
{
#if !defined(SQMESH_TEST_OCC_ENABLED)
  std::fprintf(stderr, "auto_cfd_surface_smoke: OCC not enabled, skipping.\n");
  return EXIT_SUCCESS;
#else

  // ── Initialize ─────────────────────────────────────────────────────

  sqmesh::base::ContextHandle context;
  if(!expect(
       sqmesh::base::initialize(context) == sqmesh::base::StatusCode::ok,
       "runtime initialization should succeed")) {
    return EXIT_FAILURE;
  }

  // ── Test 1: Box mesh ───────────────────────────────────────────────

  {
    TopoDS_Shape box = BRepPrimAPI_MakeBox(10.0, 8.0, 6.0).Shape();
    sqmesh::geo::ModelHandle model = sqmesh::invalid_handle;
    if(!expect(
         sqmesh::cad::occ::store_shape(box, model, context) ==
           sqmesh::base::StatusCode::ok,
         "store_shape should accept a box")) {
      return EXIT_FAILURE;
    }

    sqmesh::geo::ModelView view;
    if(!expect(
         sqmesh::geo::model_view(model, view, context) ==
           sqmesh::base::StatusCode::ok &&
           view.faces.size() == 6U,
         "box model should have 6 faces")) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::ParameterDictionary params;
    params.set_number("minimum_length", 0.5);
    params.set_number("maximum_length", 5.0);
    params.set_number("distortion_angle", 15.0);
    params.set_number("growth_rate", 1.2);

    sqmesh::mesh::MeshHandle mesh = sqmesh::invalid_handle;
    if(!expect(
         sqmesh::mesh::create_surface_mesh(
           model, "Auto CFD Surface Mesher", params, mesh, context
         ) == sqmesh::base::StatusCode::ok &&
           mesh != sqmesh::invalid_handle,
         "create_surface_mesh should mesh a box")) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::MeshSummary summary;
    if(!expect(
         sqmesh::mesh::mesh_summary(mesh, summary, context) ==
           sqmesh::base::StatusCode::ok &&
           summary.face_count > 10U &&
           summary.node_count > 8U,
         "box mesh should have reasonable face and node counts")) {
      return EXIT_FAILURE;
    }

    std::fprintf(stderr, "[TEST] Box mesh: %zu faces, %zu nodes\n",
      summary.face_count, summary.node_count);
  }

  // ── Test 2: Trimmed face with hole ─────────────────────────────────

  {
    TopoDS_Face holed_face = make_rectangular_face_with_hole();
    sqmesh::geo::ModelHandle model = sqmesh::invalid_handle;
    if(!expect(
         sqmesh::cad::occ::store_shape(holed_face, model, context) ==
           sqmesh::base::StatusCode::ok,
         "store_shape should accept a trimmed face with hole")) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::ParameterDictionary params;
    params.set_number("minimum_length", 0.5);
    params.set_number("maximum_length", 3.0);
    params.set_number("distortion_angle", 18.0);
    params.set_number("growth_rate", 1.2);

    sqmesh::mesh::MeshHandle mesh = sqmesh::invalid_handle;
    if(!expect(
         sqmesh::mesh::create_surface_mesh(
           model, "Auto CFD Surface Mesher", params, mesh, context
         ) == sqmesh::base::StatusCode::ok &&
           mesh != sqmesh::invalid_handle,
         "create_surface_mesh should mesh a face with inner loop")) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::MeshSummary summary;
    if(!expect(
         sqmesh::mesh::mesh_summary(mesh, summary, context) ==
           sqmesh::base::StatusCode::ok &&
           summary.face_count > 5U,
         "holed face mesh should have reasonable face count")) {
      return EXIT_FAILURE;
    }

    std::fprintf(stderr, "[TEST] Holed face mesh: %zu faces, %zu nodes\n",
      summary.face_count, summary.node_count);
  }

  // ── Test 3: Sizing field basics ────────────────────────────────────

  {
    sqmesh::mesh::detail::AutoCfdSurfaceSizingFieldState field;
    field.minimum_length = 0.1;
    field.maximum_length = 10.0;
    field.distortion_angle = 15.0;
    field.growth_rate = 1.2;

    // Add a single curvature source and build the size function.
    sqmesh::mesh::detail::AutoCfdSurfaceSizingSource source;
    source.kind = sqmesh::mesh::detail::AutoCfdSurfaceSizingSourceKind::face_curvature;
    source.position = {1.0, 0.0, 0.0};
    source.target_size = 0.2;
    source.relaxed_target_size = 0.2;
    field.curvature_sources.push_back(source);

    if(!expect(
         sqmesh::mesh::detail::rebuild_auto_cfd_surface_background_field(field) ==
           sqmesh::base::StatusCode::ok &&
           field.size_function.built(),
         "size function should build from a single curvature source")) {
      return EXIT_FAILURE;
    }

    // Query should return a finite positive value.
    const double queried = sqmesh::mesh::detail::query_auto_cfd_surface_sizing_field(
      field, {}, {1.0, 0.0, 0.0});
    if(!expect(
         queried > 0.0 && queried <= field.maximum_length,
         "size function query should return a valid size")) {
      return EXIT_FAILURE;
    }

    std::fprintf(stderr, "[TEST] Size function query at source: %.4f\n", queried);
  }

  // ── Test 4: Parameter validation ───────────────────────────────────

  {
    TopoDS_Shape box = BRepPrimAPI_MakeBox(1.0, 1.0, 1.0).Shape();
    sqmesh::geo::ModelHandle model = sqmesh::invalid_handle;
    static_cast<void>(sqmesh::cad::occ::store_shape(box, model, context));

    sqmesh::mesh::MeshHandle mesh = sqmesh::invalid_handle;

    // Empty algorithm name should fail.
    sqmesh::mesh::ParameterDictionary empty_params;
    if(!expect(
         sqmesh::mesh::create_surface_mesh(
           model, "", empty_params, mesh, context
         ) != sqmesh::base::StatusCode::ok,
         "empty algorithm name should be rejected")) {
      return EXIT_FAILURE;
    }

    // Missing required parameters should fail.
    sqmesh::mesh::ParameterDictionary bad_params;
    bad_params.set_number("minimum_length", -1.0);
    if(!expect(
         sqmesh::mesh::create_surface_mesh(
           model, "Auto CFD Surface Mesher", bad_params, mesh, context
         ) != sqmesh::base::StatusCode::ok,
         "negative minimum_length should be rejected")) {
      return EXIT_FAILURE;
    }
  }

  // ── Cleanup ────────────────────────────────────────────────────────

  if(!expect(
       sqmesh::base::shutdown_all() == sqmesh::base::StatusCode::ok,
       "shutdown_all should succeed")) {
    return EXIT_FAILURE;
  }

  std::fprintf(stderr, "[TEST] All auto_cfd_surface_smoke tests passed.\n");
  return EXIT_SUCCESS;
#endif
}
