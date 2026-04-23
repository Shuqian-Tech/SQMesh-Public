// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "sqmesh/sqmesh.hpp"

#include "../src/mesh/tet/tet_volume_mesher.hpp"

#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string_view>
#include <unordered_set>

#if defined(SQMESH_TEST_OCC_ENABLED)
#include "../src/cad/occ/adapter.hpp"

#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <TopoDS_Shape.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#endif

namespace {

constexpr std::string_view kVolumeMesherName = "Tetrahedral Volume Mesher";
constexpr std::string_view kVolumeMesherAlias = "tetrahedral_volume_mesher";
constexpr std::string_view kExpVolumeMesherName = "Experimental Tetrahedral Volume Mesher";
constexpr std::string_view kExpVolumeMesherAlias = "experimental_volume_mesher";
constexpr std::string_view kSurfaceMesherName = "Auto CFD Surface Mesher";
constexpr double kProjectionTolerance = 1.0e-6;

bool expect(bool condition, const char *message)
{
  if(condition) {
    return true;
  }

  std::fprintf(stderr, "volume_mesher_smoke: %s\n", message);
  return false;
}

bool expect_with_last_error(bool condition, const char *message)
{
  if(condition) {
    return true;
  }

  const auto error = sqmesh::base::last_error_message();
  std::fprintf(
    stderr,
    "volume_mesher_smoke: %s (last_error=%.*s)\n",
    message,
    static_cast<int>(error.size()),
    error.data()
  );
  return false;
}

[[nodiscard]] std::uint64_t pack_ref(sqmesh::mesh::EntityRef ref) noexcept
{
  return (static_cast<std::uint64_t>(ref.entity_group) << 32U) | ref.index;
}

[[nodiscard]] const sqmesh::mesh::KindQualitySummary *find_kind_summary(
  const sqmesh::mesh::MeshQualityReport &report,
  sqmesh::mesh::EntityKind kind
) noexcept
{
  const auto it = std::find_if(
    report.kinds.begin(),
    report.kinds.end(),
    [kind](const auto &entry) { return entry.kind == kind; }
  );
  return it == report.kinds.end() ? nullptr : &(*it);
}

sqmesh::mesh::MeshingOptions make_surface_options(
  double minimum_length,
  double maximum_length
)
{
  sqmesh::mesh::MeshingOptions options;
  options.parameters.set_number("minimum_length", minimum_length);
  options.parameters.set_number("maximum_length", maximum_length);
  options.parameters.set_number("distortion_angle", 15.0);
  options.parameters.set_number("growth_rate", 1.2);
  options.parameters.set_boolean("proximity", false);
  options.parameters.set_boolean("allow_quality_gate_failure", true);
  options.parameters.set_boolean("allow_final_screen_failure", true);
  return options;
}

sqmesh::mesh::ParameterDictionary make_tet_parameters(
  double minimum_length,
  double maximum_length
)
{
  sqmesh::mesh::ParameterDictionary parameters;
  parameters.set_number("minimum_length", minimum_length);
  parameters.set_number("maximum_length", maximum_length);
  parameters.set_number("growth_rate", 1.2);
  parameters.set_number("quality_ratio", 1.414);
  return parameters;
}

sqmesh::base::StatusCode create_surface_then_volume_mesh(
  sqmesh::geo::ModelHandle model_handle,
  const sqmesh::mesh::MeshingOptions &surface_options,
  const sqmesh::mesh::ParameterDictionary &tet_parameters,
  sqmesh::mesh::MeshHandle &surface_mesh,
  sqmesh::mesh::MeshHandle &volume_mesh,
  sqmesh::base::ContextHandle context
)
{
  surface_mesh = sqmesh::invalid_handle;
  volume_mesh = sqmesh::invalid_handle;

  auto status = sqmesh::mesh::create_surface_mesh(
    model_handle,
    kSurfaceMesherName,
    surface_options,
    surface_mesh,
    context
  );
  if(status != sqmesh::base::StatusCode::ok) {
    return status;
  }

  return sqmesh::mesh::create_volume_mesh(
    model_handle,
    kVolumeMesherName,
    tet_parameters,
    volume_mesh,
    context
  );
}

[[nodiscard]] std::size_t count_boundary_entity_groups(
  const sqmesh::mesh::Domain &domain
) noexcept
{
  std::size_t count = 0U;
  for(const auto &entity_group : domain.entity_groups()) {
    if(entity_group.role() == sqmesh::mesh::EntityGroupRole::computational &&
       entity_group.order() == sqmesh::mesh::EntityOrder::face &&
       entity_group.is_boundary()) {
      ++count;
    }
  }
  return count;
}

bool all_boundary_nodes_project_to_model_faces(
  sqmesh::geo::ModelHandle model_handle,
  const sqmesh::mesh::Domain &domain,
  sqmesh::base::ContextHandle context
)
{
  sqmesh::geo::ModelView geometry_view;
  if(sqmesh::geo::model_view(model_handle, geometry_view, context) !=
     sqmesh::base::StatusCode::ok) {
    return false;
  }

  std::unordered_set<std::uint64_t> boundary_nodes;
  for(const auto &entity_group : domain.entity_groups()) {
    if(entity_group.order() != sqmesh::mesh::EntityOrder::face || !entity_group.is_boundary()) {
      continue;
    }

    for(std::uint32_t face_index = 0; face_index < entity_group.faces().size(); ++face_index) {
      const auto face_ref = sqmesh::mesh::EntityRef {entity_group.id(), face_index};
      for(const auto node_ref : domain.face_nodes(face_ref)) {
        boundary_nodes.insert(pack_ref(node_ref));
      }
    }
  }

  if(boundary_nodes.empty()) {
    return false;
  }

  for(const auto packed_ref : boundary_nodes) {
    const auto node_ref = sqmesh::mesh::EntityRef {
      static_cast<sqmesh::mesh::EntityGroupIndex>(packed_ref >> 32U),
      static_cast<std::uint32_t>(packed_ref & 0xFFFFFFFFU),
    };
    const auto &node = domain.node(node_ref);

    bool projected = false;
    double best_distance = std::numeric_limits<double>::infinity();
    for(const auto &face_view : geometry_view.faces) {
      sqmesh::geo::FaceProjection projection;
      const auto status =
        sqmesh::geo::project_point_to_face(face_view, node.coordinates, projection);
      if(status != sqmesh::base::StatusCode::ok) {
        continue;
      }

      projected = true;
      best_distance = std::min(best_distance, projection.distance);
    }

    if(!projected || best_distance > kProjectionTolerance) {
      return false;
    }
  }

  return true;
}

#if defined(SQMESH_TEST_OCC_ENABLED)
TopoDS_Shape make_triangular_prism()
{
  BRepBuilderAPI_MakePolygon base;
  base.Add(gp_Pnt(0.0, 0.0, 0.0));
  base.Add(gp_Pnt(14.0, 0.0, 0.0));
  base.Add(gp_Pnt(4.0, 9.0, 0.0));
  base.Close();

  const TopoDS_Face base_face = BRepBuilderAPI_MakeFace(base.Wire()).Face();
  return BRepPrimAPI_MakePrism(base_face, gp_Vec(0.0, 0.0, 12.0)).Shape();
}

bool verify_tetra_volume_mesh(
  sqmesh::geo::ModelHandle model_handle,
  const char *fixture_name,
  sqmesh::base::ContextHandle context
)
{
  const auto surface_options = make_surface_options(10.0, 200.0);
  const auto tet_parameters = make_tet_parameters(10.0, 200.0);

  sqmesh::mesh::MeshHandle surface_mesh = sqmesh::invalid_handle;
  sqmesh::mesh::MeshHandle volume_mesh = sqmesh::invalid_handle;
  if(!expect_with_last_error(
       create_surface_then_volume_mesh(
         model_handle,
         surface_options,
         tet_parameters,
         surface_mesh,
         volume_mesh,
         context
       ) == sqmesh::base::StatusCode::ok,
       fixture_name
     )) {
    return false;
  }

  sqmesh::mesh::MeshSummary summary;
  sqmesh::mesh::Domain domain;
  sqmesh::mesh::MeshQualityReport quality_report;
  if(!expect(
       sqmesh::mesh::mesh_summary(volume_mesh, summary, context) ==
         sqmesh::base::StatusCode::ok &&
       sqmesh::mesh::domain_snapshot(volume_mesh, domain, context) ==
         sqmesh::base::StatusCode::ok &&
       sqmesh::mesh::mesh_quality_report(volume_mesh, quality_report, context) ==
         sqmesh::base::StatusCode::ok,
       "generated tetra mesh should expose summary, snapshot, and quality data"
     )) {
    return false;
  }

  if(!expect(
       summary.node_count > 0U &&
         summary.face_count > 0U &&
         summary.cell_count > 0U,
       "generated tetra mesh should contain nodes, faces, and cells"
     )) {
    return false;
  }

  const auto *tetra_summary =
    find_kind_summary(quality_report, sqmesh::mesh::EntityKind::cell_tetra);
  if(!expect(
       tetra_summary != nullptr &&
         tetra_summary->supported_element_count == summary.cell_count &&
         tetra_summary->valid_element_count == summary.cell_count &&
         tetra_summary->degenerate_element_count == 0U &&
         tetra_summary->inverted_element_count == 0U &&
         tetra_summary->jacobian.minimum > 0.0 &&
         tetra_summary->radius_ratio.minimum > 0.0,
       "generated tetra mesh should be valid and positively oriented"
     )) {
    return false;
  }

  if(!expect(
       count_boundary_entity_groups(domain) > 0U,
       "generated tetra mesh should expose boundary face entity_groups"
     )) {
    return false;
  }

  if(!expect(
       all_boundary_nodes_project_to_model_faces(model_handle, domain, context),
       "generated tetra mesh boundary nodes should project back to source geometry"
     )) {
    return false;
  }

  return true;
}

bool verify_coplanar_tet_failure_publishes_diagnostics(
  sqmesh::base::ContextHandle context
)
{
  sqmesh::mesh::Domain domain;

  sqmesh::mesh::EntityGroupDefinition node_definition;
  node_definition.order = sqmesh::mesh::EntityOrder::node;
  node_definition.name = "surface_nodes";
  node_definition.role = sqmesh::mesh::EntityGroupRole::computational;
  const auto node_entity_group = domain.create_entity_group(std::move(node_definition));

  sqmesh::mesh::EntityGroupDefinition face_definition;
  face_definition.order = sqmesh::mesh::EntityOrder::face;
  face_definition.name = "surface_faces";
  face_definition.role = sqmesh::mesh::EntityGroupRole::computational;
  face_definition.boundary = true;
  face_definition.default_kind = sqmesh::mesh::EntityKind::face_triangle;
  face_definition.zone_id = 7;
  const auto face_entity_group = domain.create_entity_group(std::move(face_definition));

  const auto n0 = domain.add_node(node_entity_group, {0.0, 0.0, 0.0});
  const auto n1 = domain.add_node(node_entity_group, {1.0, 0.0, 0.0});
  const auto n2 = domain.add_node(node_entity_group, {1.0, 1.0, 0.0});
  const auto n3 = domain.add_node(node_entity_group, {0.0, 1.0, 0.0});

  static_cast<void>(domain.add_triangle_face(face_entity_group, {n0, n1, n2}));
  static_cast<void>(domain.add_triangle_face(face_entity_group, {n0, n2, n3}));

  auto algorithm = sqmesh::mesh::detail::create_native_tet_volume_mesher();
  sqmesh::mesh::detail::MeshingRequest request;
  request.context_handle = context;
  request.target_dimension = sqmesh::mesh::detail::MeshingDimension::volume;

  const auto status = algorithm->execute(
    request,
    make_tet_parameters(10.0, 200.0),
    domain
  );
  if(!expect(
       status == sqmesh::base::StatusCode::invalid_argument,
       "coplanar tet-engine failures should preserve invalid_argument status"
     )) {
    return false;
  }

  const auto error = sqmesh::base::last_error_message();
  if(!expect(
       !error.empty(),
       "tet-engine failures should publish a diagnostic message"
     )) {
    return false;
  }

  return expect(
    error.find("build_initial_delaunay") != std::string_view::npos,
    "tet-engine diagnostics should include the last failing stage"
  );
}

bool verify_initial_delaunay_reject_summary_publishes_diagnostics(
  sqmesh::base::ContextHandle context
)
{
  sqmesh::mesh::Domain domain;

  sqmesh::mesh::EntityGroupDefinition node_definition;
  node_definition.order = sqmesh::mesh::EntityOrder::node;
  node_definition.name = "reject_surface_nodes";
  node_definition.role = sqmesh::mesh::EntityGroupRole::computational;
  const auto node_entity_group = domain.create_entity_group(std::move(node_definition));

  sqmesh::mesh::EntityGroupDefinition face_definition;
  face_definition.order = sqmesh::mesh::EntityOrder::face;
  face_definition.name = "reject_surface_faces";
  face_definition.role = sqmesh::mesh::EntityGroupRole::computational;
  face_definition.boundary = true;
  face_definition.default_kind = sqmesh::mesh::EntityKind::face_triangle;
  face_definition.zone_id = 9;
  const auto face_entity_group = domain.create_entity_group(std::move(face_definition));

  const auto n0 = domain.add_node(node_entity_group, {0.0, 0.0, 0.0});
  const auto n1 = domain.add_node(node_entity_group, {1.0, 0.0, 0.0});
  const auto n2 = domain.add_node(node_entity_group, {0.0, 1.0, 0.0});
  const auto n3 = domain.add_node(node_entity_group, {0.0, 0.0, 1.0});
  const auto n4 = domain.add_node(node_entity_group, {1.0e308, 1.0e308, 1.0e308});
  const auto n5 = domain.add_node(
    node_entity_group,
    {-0.27446388611905248, -0.37235110771364455, 0.0}
  );
  const auto n6 = domain.add_node(
    node_entity_group,
    {-626012511736.9292, 10695329907.956909, 55260990680.331665}
  );

  static_cast<void>(domain.add_triangle_face(face_entity_group, {n0, n1, n2}));
  static_cast<void>(domain.add_triangle_face(face_entity_group, {n0, n2, n3}));
  static_cast<void>(domain.add_triangle_face(face_entity_group, {n0, n3, n4}));
  static_cast<void>(domain.add_triangle_face(face_entity_group, {n0, n4, n5}));
  static_cast<void>(domain.add_triangle_face(face_entity_group, {n0, n5, n6}));

  auto algorithm = sqmesh::mesh::detail::create_native_tet_volume_mesher();
  sqmesh::mesh::detail::MeshingRequest request;
  request.context_handle = context;
  request.target_dimension = sqmesh::mesh::detail::MeshingDimension::volume;

  const auto status = algorithm->execute(
    request,
    make_tet_parameters(10.0, 200.0),
    domain
  );
  if(!expect(
       status == sqmesh::base::StatusCode::internal_error,
       "initial Delaunay reject failures should surface as internal_error"
     )) {
    return false;
  }

  const auto error = sqmesh::base::last_error_message();
  if(!expect(
       error.find("build_initial_delaunay") != std::string_view::npos,
       "initial Delaunay reject diagnostics should include the failing stage"
     )) {
    return false;
  }
  if(!expect(
       error.find("Initial Delaunay:") != std::string_view::npos,
       "initial Delaunay reject diagnostics should include the reject summary"
     )) {
    return false;
  }

  const bool has_nonzero_locate_failed =
    error.find("locate_failed=") != std::string_view::npos &&
    error.find("locate_failed=0") == std::string_view::npos;
  return expect(
    has_nonzero_locate_failed,
    "initial Delaunay reject diagnostics should include a nonzero locate_failed count"
  );
}
#endif

} // namespace

int main()
{
  if(!expect(
       sqmesh::mesh::is_algorithm_registered(kVolumeMesherName),
       "Tetrahedral Volume Mesher should be registered in the algorithm registry"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::mesh::is_algorithm_registered(kVolumeMesherAlias),
       "tetrahedral_volume_mesher alias should be registered in the algorithm registry"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::mesh::is_algorithm_registered(kExpVolumeMesherName),
       "Experimental Tetrahedral Volume Mesher should be registered"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::mesh::is_algorithm_registered(kExpVolumeMesherAlias),
       "experimental_volume_mesher alias should be registered"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::base::ContextHandle context = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::base::initialize(context) == sqmesh::base::StatusCode::ok,
       "initialize should create a runtime context"
     )) {
    return EXIT_FAILURE;
  }

  if(!verify_coplanar_tet_failure_publishes_diagnostics(context)) {
    return EXIT_FAILURE;
  }
  if(!verify_initial_delaunay_reject_summary_publishes_diagnostics(context)) {
    return EXIT_FAILURE;
  }

  const auto tet_parameters = make_tet_parameters(10.0, 200.0);

  if(!sqmesh::geo::cad_io_available()) {
    sqmesh::geo::ModelHandle model_handle = sqmesh::invalid_handle;
    if(!expect(
         sqmesh::geo::create_placeholder_model(model_handle, context) ==
           sqmesh::base::StatusCode::ok,
         "placeholder geometry should still allocate a model handle"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::MeshHandle mesh_handle = sqmesh::invalid_handle;
    if(!expect(
         sqmesh::mesh::create_volume_mesh(
           model_handle,
           kVolumeMesherName,
           tet_parameters,
           mesh_handle,
           context
         ) == sqmesh::base::StatusCode::invalid_argument,
         "volume meshing should reject placeholder models without a seeded surface domain"
       )) {
      return EXIT_FAILURE;
    }

    return sqmesh::base::shutdown_all() == sqmesh::base::StatusCode::ok
             ? EXIT_SUCCESS
             : EXIT_FAILURE;
  }

#if !defined(SQMESH_TEST_OCC_ENABLED)
  std::fprintf(
    stderr,
    "volume_mesher_smoke: OCC runtime is enabled but the smoke target was not built with OCC.\n"
  );
  return EXIT_FAILURE;
#else
  sqmesh::mesh::MeshHandle invalid_mesh = sqmesh::invalid_handle;
  if(!expect_with_last_error(
       sqmesh::mesh::create_volume_mesh(
         static_cast<sqmesh::geo::ModelHandle>(0xDEADBEEFULL),
         kVolumeMesherName,
         tet_parameters,
         invalid_mesh,
         context
       ) == sqmesh::base::StatusCode::invalid_handle,
       "volume meshing should reject invalid source model handles"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::ModelHandle box_model = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::cad::occ::store_shape(
         BRepPrimAPI_MakeBox(10.0, 20.0, 30.0).Shape(),
         box_model,
         context
       ) == sqmesh::base::StatusCode::ok,
       "store_shape should publish the OCC box fixture"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_tetra_volume_mesh(
       box_model,
       "box fixture should generate a tetrahedral volume mesh after surface meshing",
       context
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::ModelHandle prism_model = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::cad::occ::store_shape(make_triangular_prism(), prism_model, context) ==
         sqmesh::base::StatusCode::ok,
       "store_shape should publish the triangular prism fixture"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_tetra_volume_mesh(
       prism_model,
       "triangular prism fixture should generate a tetrahedral volume mesh after surface meshing",
       context
     )) {
    return EXIT_FAILURE;
  }

  if(!expect(
       sqmesh::base::shutdown_all() == sqmesh::base::StatusCode::ok,
       "shutdown_all should release the runtime after the volume mesher smoke test"
     )) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
#endif
}
