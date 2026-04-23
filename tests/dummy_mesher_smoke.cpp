// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "sqmesh/sqmesh.hpp"

#include <cstdio>
#include <cstdlib>

namespace {

bool expect(bool condition, const char *message)
{
  if(condition) {
    return true;
  }

  std::fprintf(stderr, "dummy_mesher_smoke: %s\n", message);
  return false;
}

} // namespace

int main()
{
  if(!expect(
       sqmesh::mesh::is_algorithm_registered("Dummy Mesher"),
       "Dummy Mesher should be registered in the algorithm registry"
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

  sqmesh::geo::ModelHandle model = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::geo::create_placeholder_model(model, context) ==
         sqmesh::base::StatusCode::ok,
       "create_placeholder_model should return a model handle"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       model != sqmesh::invalid_handle,
       "model handle should be valid"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::ModelSummary model_summary;
  if(!expect(
       sqmesh::geo::model_summary(model, model_summary, context) ==
         sqmesh::base::StatusCode::ok,
       "model_summary should resolve the placeholder model handle"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       model_summary.entity_count == 0U,
       "placeholder model should remain empty in Milestone 1"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::ParameterDictionary parameters;
  parameters.set_number("minimum_length", 0.25);
  parameters.set_number("growth_rate", 1.2);
  parameters.set_text("element_type", "tetra");

  sqmesh::mesh::MeshHandle mesh_handle = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::create_volume_mesh(
         model,
         "Dummy Mesher",
         parameters,
         mesh_handle,
         context
       ) == sqmesh::base::StatusCode::ok,
       "create_volume_mesh should run the registered Dummy Mesher"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       mesh_handle != sqmesh::invalid_handle,
       "create_volume_mesh should return a mesh handle"
     )) {
    return EXIT_FAILURE;
  }

  std::size_t node_count = 0U;
  std::size_t cell_count = 0U;
  sqmesh::mesh::MeshSummary mesh_summary;
  sqmesh::mesh::Domain mesh_domain;

  if(!expect(
       sqmesh::mesh::nodes_count(mesh_handle, node_count, context) ==
         sqmesh::base::StatusCode::ok,
       "nodes_count should resolve the generated mesh handle"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::mesh::cells_count(mesh_handle, cell_count, context) ==
         sqmesh::base::StatusCode::ok,
       "cells_count should resolve the generated mesh handle"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::mesh::mesh_summary(mesh_handle, mesh_summary, context) ==
         sqmesh::base::StatusCode::ok,
       "mesh_summary should resolve the generated mesh handle"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::mesh::domain_snapshot(mesh_handle, mesh_domain, context) ==
         sqmesh::base::StatusCode::ok,
       "domain_snapshot should resolve the generated mesh handle to a real mesh domain"
     )) {
    return EXIT_FAILURE;
  }

  if(!expect(node_count == 4U, "Dummy Mesher should generate four nodes")) {
    return EXIT_FAILURE;
  }
  if(!expect(cell_count == 1U, "Dummy Mesher should generate one tetrahedral cell")) {
    return EXIT_FAILURE;
  }
  if(!expect(
       mesh_summary.node_count == 4U &&
         mesh_summary.edge_count == 0U &&
         mesh_summary.face_count == 4U &&
         mesh_summary.cell_count == 1U,
       "mesh_summary should report the dummy tetra counts"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       mesh_domain.entity_group_count() == 3U,
       "the runtime-backed mesh should retain its domain entity_group table"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       mesh_domain.summary().node_count == 4U &&
         mesh_domain.summary().edge_count == 0U &&
         mesh_domain.summary().face_count == 4U &&
         mesh_domain.summary().cell_count == 1U,
       "domain_snapshot should preserve the generated dummy tetra counts"
     )) {
    return EXIT_FAILURE;
  }

  const auto face_entity_group = mesh_domain.entity_group(sqmesh::mesh::EntityGroupIndex(1)).id();
  const auto cell_ref = sqmesh::mesh::EntityRef {
    sqmesh::mesh::EntityGroupIndex(2),
    0U,
  };
  const auto face_ref = sqmesh::mesh::EntityRef {
    face_entity_group,
    0U,
  };
  const auto face_nodes = mesh_domain.face_nodes(face_ref);
  const auto cell_faces = mesh_domain.cell_faces(cell_ref);
  if(!expect(face_nodes.size == 3U, "runtime-backed face connectivity should stay intact")) {
    return EXIT_FAILURE;
  }
  if(!expect(cell_faces.size == 4U, "runtime-backed cell connectivity should stay intact")) {
    return EXIT_FAILURE;
  }
  if(!expect(
       mesh_domain.adjacent_cell(face_ref, sqmesh::mesh::FaceSide::left) == cell_ref,
       "runtime-backed face adjacency should preserve left-cell ownership"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       !sqmesh::mesh::is_valid(
         mesh_domain.adjacent_cell(face_ref, sqmesh::mesh::FaceSide::right)
       ),
       "runtime-backed boundary faces should still have no right cell"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle invalid_mesh = sqmesh::invalid_handle;
  sqmesh::mesh::ParameterDictionary invalid_parameters;
  invalid_parameters.set_number("minimum_length", -1.0);
  if(!expect(
       sqmesh::mesh::create_volume_mesh(
         model,
         "Dummy Mesher",
         invalid_parameters,
         invalid_mesh,
         context
       ) == sqmesh::base::StatusCode::invalid_argument,
       "invalid Dummy Mesher parameters should be rejected through the SDK"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       invalid_mesh == sqmesh::invalid_handle,
       "failed meshing should not allocate a mesh handle"
     )) {
    return EXIT_FAILURE;
  }

  if(!expect(
       sqmesh::base::shutdown_all() == sqmesh::base::StatusCode::ok,
       "shutdown_all should release the runtime after the smoke test"
     )) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
