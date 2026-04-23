// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "sqmesh/sqmesh.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#if defined(SQMESH_TEST_OCC_ENABLED)
#include "../src/cad/occ/adapter.hpp"

#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <gp_Dir.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#endif

namespace {

constexpr double kProjectionTolerance = 1.0e-6;

bool expect(bool condition, const char *message)
{
  if(condition) {
    return true;
  }

  std::fprintf(stderr, "surface_mesher_smoke: %s\n", message);
  return false;
}

bool expect_with_last_error(bool condition, const char *message)
{
  if(condition) {
    return true;
  }

  std::fprintf(
    stderr,
    "surface_mesher_smoke: %s (last_error=%.*s)\n",
    message,
    static_cast<int>(sqmesh::base::last_error_message().size()),
    sqmesh::base::last_error_message().data()
  );
  return false;
}

double axis_span(double minimum, double maximum)
{
  return maximum - minimum;
}

#if defined(SQMESH_TEST_OCC_ENABLED)
TopoDS_Face make_rectangular_face_with_hole()
{
  BRepBuilderAPI_MakePolygon outer;
  outer.Add(gp_Pnt(0.0, 0.0, 0.0));
  outer.Add(gp_Pnt(10.0, 0.0, 0.0));
  outer.Add(gp_Pnt(10.0, 8.0, 0.0));
  outer.Add(gp_Pnt(0.0, 8.0, 0.0));
  outer.Close();

  BRepBuilderAPI_MakePolygon inner;
  inner.Add(gp_Pnt(3.0, 2.0, 0.0));
  inner.Add(gp_Pnt(7.0, 2.0, 0.0));
  inner.Add(gp_Pnt(7.0, 6.0, 0.0));
  inner.Add(gp_Pnt(3.0, 6.0, 0.0));
  inner.Close();

  BRepBuilderAPI_MakeFace builder(
    gp_Pln(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, 0.0, 1.0)),
    outer.Wire(),
    true
  );
  builder.Add(inner.Wire());
  return builder.Face();
}

TopoDS_Face make_rectangular_face_with_two_holes()
{
  BRepBuilderAPI_MakePolygon outer;
  outer.Add(gp_Pnt(0.0, 0.0, 0.0));
  outer.Add(gp_Pnt(12.0, 0.0, 0.0));
  outer.Add(gp_Pnt(12.0, 8.0, 0.0));
  outer.Add(gp_Pnt(0.0, 8.0, 0.0));
  outer.Close();

  BRepBuilderAPI_MakePolygon first_inner;
  first_inner.Add(gp_Pnt(2.0, 2.0, 0.0));
  first_inner.Add(gp_Pnt(4.0, 2.0, 0.0));
  first_inner.Add(gp_Pnt(4.0, 4.0, 0.0));
  first_inner.Add(gp_Pnt(2.0, 4.0, 0.0));
  first_inner.Close();

  BRepBuilderAPI_MakePolygon second_inner;
  second_inner.Add(gp_Pnt(8.0, 3.0, 0.0));
  second_inner.Add(gp_Pnt(10.0, 3.0, 0.0));
  second_inner.Add(gp_Pnt(10.0, 6.0, 0.0));
  second_inner.Add(gp_Pnt(8.0, 6.0, 0.0));
  second_inner.Close();

  BRepBuilderAPI_MakeFace builder(
    gp_Pln(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, 0.0, 1.0)),
    outer.Wire(),
    true
  );
  builder.Add(first_inner.Wire());
  builder.Add(second_inner.Wire());
  return builder.Face();
}

TopoDS_Face find_seam_face(
  const TopoDS_Shape &shape,
  sqmesh::base::ContextHandle context
)
{
  for(TopExp_Explorer face_it(shape, TopAbs_FACE); face_it.More(); face_it.Next()) {
    const TopoDS_Face face = TopoDS::Face(face_it.Current());

    sqmesh::geo::ModelHandle face_handle = sqmesh::invalid_handle;
    if(sqmesh::cad::occ::store_shape(face, face_handle, context) !=
       sqmesh::base::StatusCode::ok) {
      continue;
    }

    sqmesh::geo::ModelView face_view;
    if(sqmesh::geo::model_view(face_handle, face_view, context) !=
         sqmesh::base::StatusCode::ok ||
       face_view.faces.size() != 1U) {
      continue;
    }

    if(sqmesh::geo::ordered_boundary_loops(face_view.faces.front()).has_seams) {
      return face;
    }
  }

  return {};
}
#endif

std::uint64_t pack_entity_ref(sqmesh::mesh::EntityRef ref)
{
  return (static_cast<std::uint64_t>(ref.entity_group) << 32U) |
         static_cast<std::uint64_t>(ref.index);
}

bool all_nodes_project_to_model_faces(
  sqmesh::geo::ModelHandle model_handle,
  const sqmesh::mesh::Domain &domain,
  sqmesh::base::ContextHandle context
)
{
  sqmesh::geo::TopologySnapshot topology_snapshot;
  if(sqmesh::geo::topology_snapshot(model_handle, topology_snapshot, context) !=
     sqmesh::base::StatusCode::ok) {
    return false;
  }

  for(const auto &entity_group : domain.entity_groups()) {
    if(entity_group.order() != sqmesh::mesh::EntityOrder::node) {
      continue;
    }

    for(const auto &node : entity_group.nodes()) {
      bool projected = false;
      double best_distance = std::numeric_limits<double>::infinity();

      for(const auto &face_info : topology_snapshot.faces) {
        sqmesh::geo::FaceProjection projection;
        const auto status = sqmesh::geo::project_point_to_face(
          model_handle,
          face_info.entity,
          node.coordinates,
          projection,
          context
        );
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
  }

  return true;
}

bool face_views_share_edge(
  const sqmesh::geo::FaceView &lhs,
  const sqmesh::geo::FaceView &rhs
)
{
  for(const auto lhs_edge : lhs.edge_ids) {
    if(std::find(rhs.edge_ids.begin(), rhs.edge_ids.end(), lhs_edge) != rhs.edge_ids.end()) {
      return true;
    }
  }
  return false;
}

sqmesh::geo::TopologyEntityId find_disjoint_face(
  const sqmesh::geo::ModelView &geometry_view,
  sqmesh::geo::TopologyEntityId face
)
{
  const auto *target_face = geometry_view.find_face(face);
  if(target_face == nullptr) {
    return {};
  }

  for(const auto &candidate_face : geometry_view.faces) {
    if(candidate_face.entity == face) {
      continue;
    }
    if(!face_views_share_edge(*target_face, candidate_face)) {
      return candidate_face.entity;
    }
  }

  return {};
}

std::vector<std::size_t> count_surface_faces_per_geometry_face(
  sqmesh::geo::ModelHandle model_handle,
  const sqmesh::mesh::Domain &domain,
  sqmesh::base::ContextHandle context
)
{
  sqmesh::geo::ModelView geometry_view;
  if(sqmesh::geo::model_view(model_handle, geometry_view, context) !=
     sqmesh::base::StatusCode::ok) {
    return {};
  }

  std::vector<std::size_t> counts(geometry_view.faces.size(), 0U);
  for(const auto &entity_group : domain.entity_groups()) {
    if(entity_group.order() != sqmesh::mesh::EntityOrder::face) {
      continue;
    }

    for(std::uint32_t face_index = 0; face_index < entity_group.faces().size(); ++face_index) {
      const auto face_ref = sqmesh::mesh::EntityRef {entity_group.id(), face_index};
      const auto face_nodes = domain.face_nodes(face_ref);
      if(face_nodes.size != 3U) {
        continue;
      }

      sqmesh::geo::Point3 centroid {0.0, 0.0, 0.0};
      for(const auto node_ref : face_nodes) {
        const auto &node = domain.node(node_ref);
        centroid[0] += node.coordinates[0];
        centroid[1] += node.coordinates[1];
        centroid[2] += node.coordinates[2];
      }
      centroid[0] /= 3.0;
      centroid[1] /= 3.0;
      centroid[2] /= 3.0;

      double best_distance = std::numeric_limits<double>::infinity();
      std::uint32_t best_face_index = sqmesh::geo::invalid_topology_index;
      for(const auto &face_view : geometry_view.faces) {
        sqmesh::geo::FaceProjection projection;
        const auto status =
          sqmesh::geo::project_point_to_face(face_view, centroid, projection);
        if(status != sqmesh::base::StatusCode::ok || projection.distance > best_distance) {
          continue;
        }

        best_distance = projection.distance;
        best_face_index = face_view.entity.index;
      }

      if(best_face_index != sqmesh::geo::invalid_topology_index &&
         best_distance <= kProjectionTolerance) {
        ++counts[best_face_index];
      }
    }
  }

  return counts;
}

bool all_surface_faces_have_topology_owner(
  const sqmesh::mesh::Domain &domain
)
{
  for(const auto &entity_group : domain.entity_groups()) {
    if(entity_group.order() != sqmesh::mesh::EntityOrder::face) {
      continue;
    }

    for(std::uint32_t face_index = 0U; face_index < entity_group.faces().size(); ++face_index) {
      const auto face_ref = sqmesh::mesh::EntityRef {entity_group.id(), face_index};
      const auto owner = domain.face_topology_owner(face_ref);
      if(!sqmesh::geo::is_valid(owner) ||
         owner.dimension != sqmesh::geo::TopologyDimension::face) {
        return false;
      }
    }
  }

  return true;
}

std::size_t boundary_edge_component_count(const sqmesh::mesh::Domain &domain)
{
  std::unordered_map<std::uint64_t, std::vector<std::uint64_t>> adjacency;
  std::unordered_set<std::uint64_t> unvisited;

  for(const auto &entity_group : domain.entity_groups()) {
    if(entity_group.order() != sqmesh::mesh::EntityOrder::edge) {
      continue;
    }
    if(entity_group.role() == sqmesh::mesh::EntityGroupRole::geometric_proxy) {
      continue;
    }
    for(std::uint32_t edge_index = 0; edge_index < entity_group.edges().size(); ++edge_index) {
      const sqmesh::mesh::EntityRef edge_ref {entity_group.id(), edge_index};
      if(sqmesh::mesh::is_valid(
           domain.adjacent_face(edge_ref, sqmesh::mesh::FaceSide::right)
         )) {
        continue;
      }

      const auto edge_nodes = domain.edge_nodes(edge_ref);
      if(edge_nodes.size != 2U) {
        continue;
      }

      const auto first = pack_entity_ref(edge_nodes[0]);
      const auto second = pack_entity_ref(edge_nodes[1]);
      adjacency[first].push_back(second);
      adjacency[second].push_back(first);
      unvisited.insert(first);
      unvisited.insert(second);
    }
  }

  std::size_t component_count = 0U;
  std::vector<std::uint64_t> stack;
  while(!unvisited.empty()) {
    ++component_count;
    stack.push_back(*unvisited.begin());
    unvisited.erase(unvisited.begin());

    while(!stack.empty()) {
      const auto current = stack.back();
      stack.pop_back();
      const auto adjacency_it = adjacency.find(current);
      if(adjacency_it == adjacency.end()) {
        continue;
      }

      for(const auto neighbor : adjacency_it->second) {
        const auto unvisited_it = unvisited.find(neighbor);
        if(unvisited_it == unvisited.end()) {
          continue;
        }
        stack.push_back(neighbor);
        unvisited.erase(unvisited_it);
      }
    }
  }

  return component_count;
}

} // namespace

int main()
{
  if(!expect(
       sqmesh::mesh::is_algorithm_registered("Auto CFD Surface Mesher"),
       "Auto CFD Surface Mesher should be registered in the algorithm registry"
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

  if(!sqmesh::geo::cad_io_available()) {
    sqmesh::geo::ModelHandle model_handle = sqmesh::invalid_handle;
    if(!expect(
         sqmesh::geo::create_placeholder_model(model_handle, context) ==
           sqmesh::base::StatusCode::ok,
         "placeholder geometry should still allocate a model handle"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::ParameterDictionary parameters;
    parameters.set_number("minimum_length", 1.0);
    parameters.set_number("maximum_length", 2.0);
    parameters.set_number("distortion_angle", 18.0);
    parameters.set_number("growth_rate", 1.2);
    parameters.set_text("element_type", "tri");

    sqmesh::mesh::MeshHandle mesh_handle = sqmesh::invalid_handle;
    if(!expect(
         sqmesh::mesh::create_surface_mesh(
           model_handle,
           "Auto CFD Surface Mesher",
           parameters,
           mesh_handle,
           context
         ) == sqmesh::base::StatusCode::invalid_argument,
         "surface meshing should reject placeholder models without real geometry faces"
       )) {
      return EXIT_FAILURE;
    }
    if(!expect(
         mesh_handle == sqmesh::invalid_handle,
         "failed surface meshing should not allocate a mesh handle"
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
    "surface_mesher_smoke: OCC runtime is enabled but the smoke target was not built with OCC.\n"
  );
  return EXIT_FAILURE;
#else
  sqmesh::geo::ModelHandle model_handle = sqmesh::invalid_handle;
  const TopoDS_Shape box = BRepPrimAPI_MakeBox(10.0, 20.0, 30.0).Shape();
  if(!expect(
       sqmesh::cad::occ::store_shape(box, model_handle, context) ==
         sqmesh::base::StatusCode::ok,
       "store_shape should publish the OCC box through the public model handle pipeline"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::ModelView box_view;
  if(!expect(
       sqmesh::geo::model_view(model_handle, box_view, context) ==
         sqmesh::base::StatusCode::ok,
       "model_view should expose the OCC box geometry used by the local sizing smoke checks"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       !box_view.faces.empty(),
       "the OCC box should expose geometry faces for local sizing controls"
     )) {
    return EXIT_FAILURE;
  }
  const auto target_face = box_view.faces.front().entity;
  const auto opposite_face = find_disjoint_face(box_view, target_face);
  if(!expect(
       sqmesh::geo::is_valid(opposite_face),
       "the OCC box should expose at least one face disjoint from the targeted local-size face"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::ParameterDictionary coarse_parameters;
  coarse_parameters.set_number("minimum_length", 15.0);
  coarse_parameters.set_number("maximum_length", 30.0);
  coarse_parameters.set_number("distortion_angle", 18.0);
  coarse_parameters.set_number("growth_rate", 1.2);
  coarse_parameters.set_text("element_type", "tri");

  sqmesh::mesh::ParameterDictionary fine_parameters;
  fine_parameters.set_number("minimum_length", 5.0);
  fine_parameters.set_number("maximum_length", 15.0);
  fine_parameters.set_number("distortion_angle", 18.0);
  fine_parameters.set_number("growth_rate", 1.2);
  fine_parameters.set_text("element_type", "triangle");

  sqmesh::mesh::MeshingOptions local_face_options;
  local_face_options.parameters = coarse_parameters;
  local_face_options.size_controls.topology_revision =
    box_view.snapshot.topology_revision;
  local_face_options.size_controls.add_local_size(target_face, 5.0);

  // -- Unified Domain: generate → snapshot → assert → generate next --
  sqmesh::mesh::MeshHandle mesh_handle = sqmesh::invalid_handle;

  // 1. Coarse: generate then snapshot
  if(!expect(
       sqmesh::mesh::create_surface_mesh(
         model_handle, "Auto CFD Surface Mesher", coarse_parameters, mesh_handle, context
       ) == sqmesh::base::StatusCode::ok,
       "create_surface_mesh should generate a coarse surface mesh through the SDK"
     )) {
    return EXIT_FAILURE;
  }
  sqmesh::mesh::MeshSummary coarse_summary;
  sqmesh::mesh::Domain coarse_domain;
  if(!expect(
       sqmesh::mesh::mesh_summary(mesh_handle, coarse_summary, context) ==
         sqmesh::base::StatusCode::ok,
       "mesh_summary should resolve the coarse surface mesh handle"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::mesh::domain_snapshot(mesh_handle, coarse_domain, context) ==
         sqmesh::base::StatusCode::ok,
       "domain_snapshot should resolve the coarse surface mesh domain"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       coarse_summary.node_count > 0U && coarse_summary.edge_count > 0U &&
         coarse_summary.face_count > 0U && coarse_summary.cell_count == 0U,
       "coarse surface meshing should produce a non-empty triangle mesh without volume cells"
     )) {
    return EXIT_FAILURE;
  }

  // 2. Fine: generate, snapshot, compare with coarse snapshot
  sqmesh::mesh::MeshHandle fine_mesh = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::create_surface_mesh(
         model_handle, "Auto CFD Surface Mesher", fine_parameters, fine_mesh, context
       ) == sqmesh::base::StatusCode::ok,
       "create_surface_mesh should generate a finer surface mesh through the SDK"
     )) {
    return EXIT_FAILURE;
  }
  sqmesh::mesh::MeshSummary fine_summary;
  sqmesh::mesh::Domain fine_domain;
  if(!expect(
       sqmesh::mesh::mesh_summary(fine_mesh, fine_summary, context) ==
         sqmesh::base::StatusCode::ok,
       "mesh_summary should resolve the fine surface mesh handle"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::mesh::domain_snapshot(fine_mesh, fine_domain, context) ==
         sqmesh::base::StatusCode::ok,
       "domain_snapshot should resolve the fine surface mesh domain"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       fine_summary.node_count > coarse_summary.node_count &&
         fine_summary.edge_count > coarse_summary.edge_count &&
         fine_summary.face_count > coarse_summary.face_count &&
         fine_summary.cell_count == 0U,
       "smaller minimum_length values should produce a finer surface mesh"
     )) {
    return EXIT_FAILURE;
  }

  // 3. Local face sizing: generate, snapshot
  sqmesh::mesh::MeshHandle local_face_mesh = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::create_surface_mesh(
         model_handle, "Auto CFD Surface Mesher", local_face_options, local_face_mesh, context
       ) == sqmesh::base::StatusCode::ok,
       "create_surface_mesh should accept explicit local face size controls through the public mesh API"
     )) {
    return EXIT_FAILURE;
  }
  sqmesh::mesh::MeshSummary local_face_summary;
  sqmesh::mesh::Domain local_face_domain;
  if(!expect(
       sqmesh::mesh::mesh_summary(local_face_mesh, local_face_summary, context) ==
         sqmesh::base::StatusCode::ok &&
       sqmesh::mesh::domain_snapshot(local_face_mesh, local_face_domain, context) ==
         sqmesh::base::StatusCode::ok,
       "mesh_summary and domain_snapshot should resolve the locally controlled surface mesh"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       local_face_summary.node_count > coarse_summary.node_count &&
         local_face_summary.edge_count > coarse_summary.edge_count &&
         local_face_summary.face_count > coarse_summary.face_count,
       "a local face size control should refine the surface mesh without changing the coarse global settings"
     )) {
    return EXIT_FAILURE;
  }

  // EntityGroup count: 3 proxy + 3 surface = 6
  if(!expect(
       fine_domain.entity_group_count(sqmesh::mesh::EntityGroupRole::computational) == 3U,
       "surface meshing should produce exactly 3 computational entity_groups (node, edge, face)"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       coarse_summary.source_topology_revision == box_view.snapshot.topology_revision &&
         fine_summary.source_topology_revision == box_view.snapshot.topology_revision &&
         local_face_summary.source_topology_revision == box_view.snapshot.topology_revision &&
         coarse_domain.source_topology_revision() == box_view.snapshot.topology_revision &&
         fine_domain.source_topology_revision() == box_view.snapshot.topology_revision &&
         local_face_domain.source_topology_revision() == box_view.snapshot.topology_revision,
       "surface meshing should stamp mesh summaries and snapshots with the source topology revision"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       all_surface_faces_have_topology_owner(local_face_domain),
       "surface meshing should preserve a valid geometry-face owner for every generated triangle"
     )) {
    return EXIT_FAILURE;
  }

  const auto coarse_face_counts =
    count_surface_faces_per_geometry_face(model_handle, coarse_domain, context);
  const auto local_face_counts =
    count_surface_faces_per_geometry_face(model_handle, local_face_domain, context);
  if(!expect(
       coarse_face_counts.size() == box_view.faces.size() &&
         local_face_counts.size() == box_view.faces.size(),
       "surface local sizing checks should classify generated triangles back onto geometry faces"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       local_face_counts[target_face.index] > coarse_face_counts[target_face.index] &&
         local_face_counts[opposite_face.index] == coarse_face_counts[opposite_face.index],
       "the local face size control should densify the targeted box face while leaving the disjoint opposite face coarse"
     )) {
    return EXIT_FAILURE;
  }



  double x_min = std::numeric_limits<double>::infinity();
  double y_min = std::numeric_limits<double>::infinity();
  double z_min = std::numeric_limits<double>::infinity();
  double x_max = -std::numeric_limits<double>::infinity();
  double y_max = -std::numeric_limits<double>::infinity();
  double z_max = -std::numeric_limits<double>::infinity();

  for(const auto &entity_group : fine_domain.entity_groups()) {
    if(entity_group.order() != sqmesh::mesh::EntityOrder::node) {
      continue;
    }

    for(const auto &node : entity_group.nodes()) {
      x_min = std::min(x_min, node.coordinates[0]);
      y_min = std::min(y_min, node.coordinates[1]);
      z_min = std::min(z_min, node.coordinates[2]);
      x_max = std::max(x_max, node.coordinates[0]);
      y_max = std::max(y_max, node.coordinates[1]);
      z_max = std::max(z_max, node.coordinates[2]);
    }
  }

  if(!expect(
       std::abs(axis_span(x_min, x_max) - 10.0) < 1.0e-7 &&
         std::abs(axis_span(y_min, y_max) - 20.0) < 1.0e-7 &&
         std::abs(axis_span(z_min, z_max) - 30.0) < 1.0e-7,
       "surface mesh node coordinates should span the source box dimensions"
     )) {
    return EXIT_FAILURE;
  }

  if(!expect(
       all_nodes_project_to_model_faces(model_handle, fine_domain, context),
       "every surface mesh node should lie on or project back to the source geometry"
     )) {
    return EXIT_FAILURE;
  }

  // Trimmed face / cylinder / multi-hole tests removed — these were
  // Simple Surface Mesher specific geometries not suited for Auto CFD.

  if(!expect(
       sqmesh::base::shutdown_all() == sqmesh::base::StatusCode::ok,
       "shutdown_all should release the runtime after the surface mesher smoke test"
     )) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
#endif
}
// ---- removed Simple Surface Mesher tests below ----
#if 0
  sqmesh::geo::ModelHandle holed_face_handle = sqmesh::invalid_handle;
  const TopoDS_Face holed_face = make_rectangular_face_with_hole();
  if(!expect(
       sqmesh::cad::occ::store_shape(holed_face, holed_face_handle, context) ==
         sqmesh::base::StatusCode::ok,
       "store_shape should publish a trimmed OCC face with one inner loop"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::ParameterDictionary holed_coarse_parameters;
  holed_coarse_parameters.set_number("minimum_length", 3.0);
  holed_coarse_parameters.set_number("maximum_length", 6.0);
  holed_coarse_parameters.set_number("distortion_angle", 18.0);
  holed_coarse_parameters.set_number("growth_rate", 1.2);
  holed_coarse_parameters.set_text("element_type", "tri");

  sqmesh::mesh::ParameterDictionary holed_fine_parameters;
  holed_fine_parameters.set_number("minimum_length", 1.0);
  holed_fine_parameters.set_number("maximum_length", 3.0);
  holed_fine_parameters.set_number("distortion_angle", 18.0);
  holed_fine_parameters.set_number("growth_rate", 1.2);
  holed_fine_parameters.set_text("element_type", "triangle");

  sqmesh::mesh::MeshHandle holed_coarse_mesh = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::create_surface_mesh(
         holed_face_handle, "Auto CFD Surface Mesher",
         holed_coarse_parameters, holed_coarse_mesh, context
       ) == sqmesh::base::StatusCode::ok,
       "create_surface_mesh should mesh a trimmed face with one inner loop"
     )) {
    return EXIT_FAILURE;
  }

  // Snapshot coarse before fine overwrites
  sqmesh::mesh::MeshSummary holed_coarse_summary;
  if(!expect(
       sqmesh::mesh::mesh_summary(holed_coarse_mesh, holed_coarse_summary, context) ==
         sqmesh::base::StatusCode::ok,
       "mesh_summary should resolve the coarse trimmed-face mesh"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle holed_fine_mesh = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::create_surface_mesh(
         holed_face_handle, "Auto CFD Surface Mesher",
         holed_fine_parameters, holed_fine_mesh, context
       ) == sqmesh::base::StatusCode::ok,
       "create_surface_mesh should refine the one-hole trimmed-face mesh when minimum_length decreases"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshSummary holed_fine_summary;
  sqmesh::mesh::Domain holed_fine_domain;
  if(!expect(
       sqmesh::mesh::mesh_summary(holed_fine_mesh, holed_fine_summary, context) ==
         sqmesh::base::StatusCode::ok,
       "mesh_summary should resolve the fine trimmed-face mesh"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::mesh::domain_snapshot(holed_fine_mesh, holed_fine_domain, context) ==
         sqmesh::base::StatusCode::ok,
       "domain_snapshot should resolve the trimmed-face mesh domain"
     )) {
    return EXIT_FAILURE;
  }

  if(!expect(
       holed_coarse_summary.node_count > 0U &&
         holed_coarse_summary.edge_count > 0U &&
         holed_coarse_summary.face_count > 0U &&
         holed_coarse_summary.cell_count == 0U,
       "the coarse trimmed-face mesh should produce a non-empty triangle mesh"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       holed_fine_summary.node_count > holed_coarse_summary.node_count &&
         holed_fine_summary.edge_count > holed_coarse_summary.edge_count &&
         holed_fine_summary.face_count > holed_coarse_summary.face_count &&
         holed_fine_summary.cell_count == 0U,
       "the trimmed-face mesh should remain sensitive to the requested target size"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       holed_fine_domain.entity_group_count(sqmesh::mesh::EntityGroupRole::computational) == 3U,
       "trimmed-face meshing should produce exactly 3 computational entity_groups"
     )) {
    return EXIT_FAILURE;
  }

  const auto holed_stats = holed_fine_domain.statistics();
  if(!expect(
       holed_stats.boundary_edge_count > 0U &&
         boundary_edge_component_count(holed_fine_domain) == 2U,
       "the one-hole trimmed-face mesh should preserve separate outer and inner boundary loops"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       all_nodes_project_to_model_faces(holed_face_handle, holed_fine_domain, context),
       "every trimmed-face mesh node should lie on or project back to the trimmed source face"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::ModelHandle two_hole_handle = sqmesh::invalid_handle;
  const TopoDS_Face two_hole_face = make_rectangular_face_with_two_holes();
  if(!expect(
       sqmesh::cad::occ::store_shape(two_hole_face, two_hole_handle, context) ==
         sqmesh::base::StatusCode::ok,
       "store_shape should publish a trimmed OCC face with two inner loops for the multi-hole support check"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle two_hole_coarse_mesh = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::create_surface_mesh(
         two_hole_handle, "Auto CFD Surface Mesher",
         holed_coarse_parameters, two_hole_coarse_mesh, context
       ) == sqmesh::base::StatusCode::ok,
       "create_surface_mesh should mesh a trimmed face with multiple inner loops"
     )) {
    return EXIT_FAILURE;
  }

  // Snapshot coarse before fine overwrites
  sqmesh::mesh::MeshSummary two_hole_coarse_summary;
  if(!expect(
       sqmesh::mesh::mesh_summary(two_hole_coarse_mesh, two_hole_coarse_summary, context) ==
         sqmesh::base::StatusCode::ok,
       "mesh_summary should resolve the coarse multi-hole trimmed-face mesh"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle two_hole_fine_mesh = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::create_surface_mesh(
         two_hole_handle, "Auto CFD Surface Mesher",
         holed_fine_parameters, two_hole_fine_mesh, context
       ) == sqmesh::base::StatusCode::ok,
       "create_surface_mesh should refine the multi-hole trimmed-face mesh when minimum_length decreases"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshSummary two_hole_fine_summary;
  sqmesh::mesh::Domain two_hole_fine_domain;
  if(!expect(
       sqmesh::mesh::mesh_summary(two_hole_fine_mesh, two_hole_fine_summary, context) ==
         sqmesh::base::StatusCode::ok,
       "mesh_summary should resolve the fine multi-hole trimmed-face mesh"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::mesh::domain_snapshot(two_hole_fine_mesh, two_hole_fine_domain, context) ==
         sqmesh::base::StatusCode::ok,
       "domain_snapshot should resolve the multi-hole trimmed-face mesh domain"
     )) {
    return EXIT_FAILURE;
  }

  if(!expect(
       two_hole_coarse_summary.node_count > 0U &&
         two_hole_coarse_summary.edge_count > 0U &&
         two_hole_coarse_summary.face_count > 0U &&
         two_hole_coarse_summary.cell_count == 0U,
       "the coarse multi-hole trimmed-face mesh should produce a non-empty triangle mesh"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       two_hole_fine_summary.node_count > two_hole_coarse_summary.node_count &&
         two_hole_fine_summary.edge_count > two_hole_coarse_summary.edge_count &&
         two_hole_fine_summary.face_count > two_hole_coarse_summary.face_count &&
         two_hole_fine_summary.cell_count == 0U,
       "the multi-hole trimmed-face mesh should remain sensitive to the requested target size"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       two_hole_fine_domain.entity_group_count(sqmesh::mesh::EntityGroupRole::computational) == 3U,
       "multi-hole trimmed-face meshing should produce exactly 3 computational entity_groups"
     )) {
    return EXIT_FAILURE;
  }

  const auto two_hole_stats = two_hole_fine_domain.statistics();
  if(!expect(
       two_hole_stats.boundary_edge_count > 0U &&
         boundary_edge_component_count(two_hole_fine_domain) == 3U,
       "the multi-hole trimmed-face mesh should preserve separate outer and inner boundary components"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       all_nodes_project_to_model_faces(two_hole_handle, two_hole_fine_domain, context),
       "every multi-hole trimmed-face mesh node should lie on or project back to the trimmed source face"
     )) {
    return EXIT_FAILURE;
  }

  const TopoDS_Face cylinder_face =
    find_seam_face(BRepPrimAPI_MakeCylinder(2.0, 5.0).Shape(), context);
  if(!expect(
       !cylinder_face.IsNull(),
       "the seam-sensitive broadening check should be able to isolate a standalone cylinder side face"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::ModelHandle cylinder_handle = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::cad::occ::store_shape(cylinder_face, cylinder_handle, context) ==
         sqmesh::base::StatusCode::ok,
       "store_shape should publish the seam-bearing cylinder side face for the broadened surface-mesher check"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::ModelView cylinder_view;
  if(!expect(
       sqmesh::geo::model_view(cylinder_handle, cylinder_view, context) ==
         sqmesh::base::StatusCode::ok &&
         cylinder_view.faces.size() == 1U,
       "model_view should expose the standalone seam-bearing cylinder face"
     )) {
    return EXIT_FAILURE;
  }

  const auto &cylinder_boundary =
    sqmesh::geo::ordered_boundary_loops(cylinder_view.faces.front());
  if(!expect(
       cylinder_boundary.outer_loop_count == 1U &&
         cylinder_boundary.inner_loop_count == 0U &&
         cylinder_boundary.seam_edge_use_count == 2U &&
         cylinder_boundary.degenerate_edge_use_count == 0U &&
         cylinder_boundary.repeated_edge_use_count == 1U &&
         cylinder_boundary.has_seams,
       "the broadened surface-mesher fixture should remain explicitly seam-bearing and non-degenerate"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::ParameterDictionary cylinder_coarse_parameters;
  cylinder_coarse_parameters.set_number("minimum_length", 2.0);
  cylinder_coarse_parameters.set_number("maximum_length", 4.0);
  cylinder_coarse_parameters.set_number("distortion_angle", 18.0);
  cylinder_coarse_parameters.set_number("growth_rate", 1.2);
  cylinder_coarse_parameters.set_text("element_type", "tri");

  sqmesh::mesh::ParameterDictionary cylinder_fine_parameters;
  cylinder_fine_parameters.set_number("minimum_length", 0.75);
  cylinder_fine_parameters.set_number("maximum_length", 2.0);
  cylinder_fine_parameters.set_number("distortion_angle", 18.0);
  cylinder_fine_parameters.set_number("growth_rate", 1.2);
  cylinder_fine_parameters.set_text("element_type", "triangle");

  sqmesh::mesh::MeshHandle cylinder_coarse_mesh = sqmesh::invalid_handle;
  if(!expect_with_last_error(
       sqmesh::mesh::create_surface_mesh(
         cylinder_handle, "Auto CFD Surface Mesher",
         cylinder_coarse_parameters, cylinder_coarse_mesh, context
       ) == sqmesh::base::StatusCode::ok,
       "create_surface_mesh should mesh the seam-bearing cylinder side face"
     )) {
    return EXIT_FAILURE;
  }

  // Snapshot coarse before fine overwrites
  sqmesh::mesh::MeshSummary cylinder_coarse_summary;
  if(!expect(
       sqmesh::mesh::mesh_summary(cylinder_coarse_mesh, cylinder_coarse_summary, context) ==
         sqmesh::base::StatusCode::ok,
       "mesh_summary should resolve the coarse seam-bearing cylinder mesh"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle cylinder_fine_mesh = sqmesh::invalid_handle;
  if(!expect_with_last_error(
       sqmesh::mesh::create_surface_mesh(
         cylinder_handle, "Auto CFD Surface Mesher",
         cylinder_fine_parameters, cylinder_fine_mesh, context
       ) == sqmesh::base::StatusCode::ok,
       "create_surface_mesh should refine the seam-bearing cylinder side face when minimum_length decreases"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshSummary cylinder_fine_summary;
  sqmesh::mesh::Domain cylinder_fine_domain;
  if(!expect(
       sqmesh::mesh::mesh_summary(cylinder_fine_mesh, cylinder_fine_summary, context) ==
         sqmesh::base::StatusCode::ok,
       "mesh_summary should resolve the fine seam-bearing cylinder mesh"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::mesh::domain_snapshot(cylinder_fine_mesh, cylinder_fine_domain, context) ==
         sqmesh::base::StatusCode::ok,
       "domain_snapshot should resolve the fine seam-bearing cylinder mesh domain"
     )) {
    return EXIT_FAILURE;
  }

  if(!expect(
       cylinder_coarse_summary.node_count > 0U &&
         cylinder_coarse_summary.edge_count > 0U &&
         cylinder_coarse_summary.face_count > 0U &&
         cylinder_coarse_summary.cell_count == 0U,
       "the coarse seam-bearing cylinder mesh should produce a non-empty triangle mesh"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       cylinder_fine_summary.node_count > cylinder_coarse_summary.node_count &&
         cylinder_fine_summary.edge_count > cylinder_coarse_summary.edge_count &&
         cylinder_fine_summary.face_count > cylinder_coarse_summary.face_count &&
         cylinder_fine_summary.cell_count == 0U,
       "the seam-bearing cylinder mesh should remain sensitive to the requested target size"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       cylinder_fine_domain.entity_group_count(sqmesh::mesh::EntityGroupRole::computational) == 3U,
       "seam-bearing cylinder meshing should produce exactly 3 computational entity_groups"
     )) {
    return EXIT_FAILURE;
  }

  const auto cylinder_stats = cylinder_fine_domain.statistics();
  if(!expect(
       cylinder_stats.boundary_edge_count > 0U &&
         boundary_edge_component_count(cylinder_fine_domain) == 2U,
       "the seam-bearing cylinder mesh should keep only the two physical circular boundaries instead of exposing the seam as an extra boundary component"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       all_nodes_project_to_model_faces(cylinder_handle, cylinder_fine_domain, context),
       "every seam-bearing cylinder mesh node should lie on or project back to the non-planar source face"
     )) {
    return EXIT_FAILURE;
  }

  const TopoDS_Face cone_face =
    find_seam_face(BRepPrimAPI_MakeCone(2.0, 0.0, 5.0).Shape(), context);
  if(!expect(
       !cone_face.IsNull(),
       "the negative-path seam check should be able to isolate a cone side face with a degenerate trim edge"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::ModelHandle cone_handle = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::cad::occ::store_shape(cone_face, cone_handle, context) ==
         sqmesh::base::StatusCode::ok,
       "store_shape should publish the seam-bearing cone side face for the explicit unsupported check"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::ModelView cone_view;
  if(!expect(
       sqmesh::geo::model_view(cone_handle, cone_view, context) ==
         sqmesh::base::StatusCode::ok &&
         cone_view.faces.size() == 1U,
       "model_view should expose the standalone seam-bearing cone face"
     )) {
    return EXIT_FAILURE;
  }

  const auto &cone_boundary = sqmesh::geo::ordered_boundary_loops(cone_view.faces.front());
  if(!expect(
       cone_boundary.has_seams && cone_boundary.degenerate_edge_use_count > 0U,
       "the explicit unsupported seam fixture should remain degenerate-trim-bearing instead of collapsing into the supported cylinder subset"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle cone_mesh = sqmesh::invalid_handle;
  if(!expect_with_last_error(
       sqmesh::mesh::create_surface_mesh(
         cone_handle,
         "Auto CFD Surface Mesher",
         cylinder_coarse_parameters,
         cone_mesh,
         context
       ) == sqmesh::base::StatusCode::unsupported,
       "harder degenerate seam-bearing faces should remain explicitly unsupported"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       cone_mesh == sqmesh::invalid_handle,
       "unsupported seam-bearing faces should not allocate a mesh handle"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::base::last_error_message().find("degenerate trim edges") !=
         std::string_view::npos,
       "unsupported harder seam-bearing faces should fail with an explicit degenerate-trim diagnostic"
     )) {
    return EXIT_FAILURE;
  }

  if(!expect(
       sqmesh::base::shutdown_all() == sqmesh::base::StatusCode::ok,
       "shutdown_all should release the runtime after the surface mesher smoke test"
     )) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
#endif
