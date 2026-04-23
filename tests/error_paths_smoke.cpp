// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "sqmesh/sqmesh.hpp"

#include <array>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>
#include <string_view>

namespace {

bool expect(bool condition, const char *message)
{
  if(condition) {
    return true;
  }

  std::fprintf(stderr, "error_paths_smoke: %s\n", message);
  return false;
}

bool expect_error(
  sqmesh::base::StatusCode actual,
  sqmesh::base::StatusCode expected,
  const char *message
)
{
  if(!expect(actual == expected, message)) {
    return false;
  }

  if(expected == sqmesh::base::StatusCode::ok) {
    return expect(
      sqmesh::base::last_error_message().empty(),
      "success paths should leave an empty diagnostic message"
    );
  }

  return expect(
    !sqmesh::base::last_error_message().empty(),
    "error paths should publish a diagnostic message"
  );
}

bool write_text_file(
  const std::filesystem::path &path,
  const char *contents
)
{
  std::ofstream output(path, std::ios::trunc);
  if(!output) {
    std::fprintf(stderr, "error_paths_smoke: failed to open %s for writing\n", path.string().c_str());
    return false;
  }

  output << contents;
  return static_cast<bool>(output);
}

bool write_binary_file(
  const std::filesystem::path &path,
  const std::string &contents
)
{
  std::ofstream output(path, std::ios::binary | std::ios::trunc);
  if(!output) {
    std::fprintf(stderr, "error_paths_smoke: failed to open %s for writing\n", path.string().c_str());
    return false;
  }

  output.write(contents.data(), static_cast<std::streamsize>(contents.size()));
  return static_cast<bool>(output);
}

template <typename T>
void append_binary_value(
  std::string &contents,
  const T &value
)
{
  const auto *bytes = reinterpret_cast<const char *>(&value);
  contents.append(bytes, sizeof(T));
}

template <typename T, std::size_t N>
void append_binary_values(
  std::string &contents,
  const std::array<T, N> &values
)
{
  const auto *bytes = reinterpret_cast<const char *>(values.data());
  contents.append(bytes, sizeof(T) * values.size());
}

bool write_non_finite_binary_msh_fixture(const std::filesystem::path &path)
{
  std::string contents;
  contents += "$MeshFormat\n2.2 1 8\n";
  const int one = 1;
  append_binary_value(contents, one);
  contents += "\n$EndMeshFormat\n";
  contents += "$Nodes\n1\n";
  const int node_id = 1;
  append_binary_value(contents, node_id);
  const std::array<double, 3> node = {
    0.0,
    std::numeric_limits<double>::infinity(),
    0.0,
  };
  append_binary_values(contents, node);
  contents += "\n$EndNodes\n";
  contents += "$Elements\n0\n$EndElements\n";
  return write_binary_file(path, contents);
}

bool expect_last_error_contains(
  std::string_view token,
  const char *message
)
{
  return expect(
    sqmesh::base::last_error_message().find(token) != std::string::npos,
    message
  );
}

} // namespace

int main()
{
  constexpr sqmesh::Handle bogus_handle = 999999U;

  sqmesh::base::ContextHandle context = sqmesh::invalid_handle;
  if(!expect_error(
       sqmesh::base::initialize(context),
       sqmesh::base::StatusCode::ok,
       "initialize should create a runtime context"
     )) {
    return EXIT_FAILURE;
  }

  if(!expect(
       context != sqmesh::invalid_handle,
       "initialize should return a valid context handle"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::ModelHandle model = sqmesh::invalid_handle;
  if(!expect_error(
       sqmesh::geo::create_placeholder_model(model, context),
       sqmesh::base::StatusCode::ok,
       "create_placeholder_model should seed a valid model handle"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::ModelSummary model_summary;
  if(!expect_error(
       sqmesh::geo::model_summary(
         static_cast<sqmesh::geo::ModelHandle>(bogus_handle),
         model_summary,
         context
       ),
       sqmesh::base::StatusCode::invalid_handle,
       "invalid model handles should return invalid_handle"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::TopologySnapshot topology_snapshot;
  if(!expect_error(
       sqmesh::geo::topology_snapshot(
         static_cast<sqmesh::geo::ModelHandle>(bogus_handle),
         topology_snapshot,
         context
       ),
       sqmesh::base::StatusCode::invalid_handle,
       "topology_snapshot should reject invalid model handles"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::FaceUvBounds face_bounds;
  if(!expect_error(
       sqmesh::geo::face_uv_bounds(
         static_cast<sqmesh::geo::ModelHandle>(bogus_handle),
         {sqmesh::geo::TopologyDimension::face, 0U},
         face_bounds,
         context
       ),
       sqmesh::base::StatusCode::invalid_handle,
       "face_uv_bounds should reject invalid model handles"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::FaceSample face_sample;
  if(!expect_error(
       sqmesh::geo::sample_face(
         static_cast<sqmesh::geo::ModelHandle>(bogus_handle),
         {sqmesh::geo::TopologyDimension::face, 0U},
         0.0,
         0.0,
         face_sample,
         context
       ),
       sqmesh::base::StatusCode::invalid_handle,
       "sample_face should reject invalid model handles"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::FaceCurvatureSample face_curvature;
  if(!expect_error(
       sqmesh::geo::sample_face_curvature(
         static_cast<sqmesh::geo::ModelHandle>(bogus_handle),
         {sqmesh::geo::TopologyDimension::face, 0U},
         0.0,
         0.0,
         face_curvature,
         context
       ),
       sqmesh::base::StatusCode::invalid_handle,
       "sample_face_curvature should reject invalid model handles"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::FaceDerivatives face_derivatives;
  if(!expect_error(
       sqmesh::geo::sample_face_derivatives(
         static_cast<sqmesh::geo::ModelHandle>(bogus_handle),
         {sqmesh::geo::TopologyDimension::face, 0U},
         0.0,
         0.0,
         face_derivatives,
         context
       ),
       sqmesh::base::StatusCode::invalid_handle,
       "sample_face_derivatives should reject invalid model handles"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::FaceProjection face_projection;
  if(!expect_error(
       sqmesh::geo::project_point_to_face(
         static_cast<sqmesh::geo::ModelHandle>(bogus_handle),
         {sqmesh::geo::TopologyDimension::face, 0U},
         {0.0, 0.0, 0.0},
         face_projection,
         context
       ),
       sqmesh::base::StatusCode::invalid_handle,
       "project_point_to_face should reject invalid model handles"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::FaceUvMapping face_uv_mapping;
  if(!expect_error(
       sqmesh::geo::recover_face_uv(
         static_cast<sqmesh::geo::ModelHandle>(bogus_handle),
         {sqmesh::geo::TopologyDimension::face, 0U},
         {0.0, 0.0, 0.0},
         face_uv_mapping,
         context
       ),
       sqmesh::base::StatusCode::invalid_handle,
       "recover_face_uv should reject invalid model handles"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::EdgeCurveInfo edge_info;
  if(!expect_error(
       sqmesh::geo::edge_curve_info(
         static_cast<sqmesh::geo::ModelHandle>(bogus_handle),
         {sqmesh::geo::TopologyDimension::edge, 0U},
         edge_info,
         context
       ),
       sqmesh::base::StatusCode::invalid_handle,
       "edge_curve_info should reject invalid model handles"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::EdgeTangentSample edge_tangent;
  if(!expect_error(
       sqmesh::geo::sample_edge_tangent(
         static_cast<sqmesh::geo::ModelHandle>(bogus_handle),
         {sqmesh::geo::TopologyDimension::edge, 0U},
         0.0,
         edge_tangent,
         context
       ),
       sqmesh::base::StatusCode::invalid_handle,
       "sample_edge_tangent should reject invalid model handles"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::EdgeCurveSamples edge_curve_samples;
  if(!expect_error(
       sqmesh::geo::sample_edge_curve(
         static_cast<sqmesh::geo::ModelHandle>(bogus_handle),
         {sqmesh::geo::TopologyDimension::edge, 0U},
         {},
         edge_curve_samples,
         context
       ),
       sqmesh::base::StatusCode::invalid_handle,
       "sample_edge_curve should reject invalid model handles"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::FaceBoundaryLoops face_boundary;
  if(!expect_error(
       sqmesh::geo::face_boundary_loops(
         static_cast<sqmesh::geo::ModelHandle>(bogus_handle),
         {sqmesh::geo::TopologyDimension::face, 0U},
         face_boundary,
         context
       ),
       sqmesh::base::StatusCode::invalid_handle,
       "face_boundary_loops should reject invalid model handles"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::FeatureEdgeReport feature_report;
  if(!expect_error(
       sqmesh::geo::feature_edges(
         static_cast<sqmesh::geo::ModelHandle>(bogus_handle),
         feature_report,
         {},
         context
       ),
       sqmesh::base::StatusCode::invalid_handle,
       "feature_edges should reject invalid model handles"
     )) {
    return EXIT_FAILURE;
  }

  if(!expect_error(
       sqmesh::base::shutdown(static_cast<sqmesh::base::ContextHandle>(bogus_handle)),
       sqmesh::base::StatusCode::invalid_handle,
       "invalid context handles should return invalid_handle"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle mesh_handle = sqmesh::invalid_handle;
  sqmesh::mesh::ParameterDictionary valid_parameters;
  valid_parameters.set_number("minimum_length", 0.25);
  valid_parameters.set_number("growth_rate", 1.25);
  valid_parameters.set_text("element_type", "tetra");

  sqmesh::mesh::ParameterDictionary valid_auto_cfd_parameters;
  valid_auto_cfd_parameters.set_number("minimum_length", 0.1);
  valid_auto_cfd_parameters.set_number("maximum_length", 0.5);
  valid_auto_cfd_parameters.set_number("distortion_angle", 20.0);
  valid_auto_cfd_parameters.set_number("growth_rate", 1.2);
  valid_auto_cfd_parameters.set_boolean("proximity", false);

  if(!expect_error(
       sqmesh::mesh::create_volume_mesh(
         static_cast<sqmesh::geo::ModelHandle>(bogus_handle),
         "Dummy Mesher",
         valid_parameters,
         mesh_handle,
         context
       ),
       sqmesh::base::StatusCode::invalid_handle,
       "meshing should reject invalid source model handles"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       mesh_handle == sqmesh::invalid_handle,
       "failed meshing should not allocate a mesh handle"
     )) {
    return EXIT_FAILURE;
  }

  if(!expect_error(
       sqmesh::mesh::create_volume_mesh(
         model,
         "",
         valid_parameters,
         mesh_handle,
         context
       ),
       sqmesh::base::StatusCode::invalid_argument,
       "meshing should reject an empty algorithm name"
     )) {
    return EXIT_FAILURE;
  }

  if(!expect_error(
       sqmesh::mesh::create_surface_mesh(
         static_cast<sqmesh::geo::ModelHandle>(bogus_handle),
         "Auto CFD Surface Mesher",
         valid_auto_cfd_parameters,
         mesh_handle,
         context
       ),
       sqmesh::base::StatusCode::invalid_handle,
       "surface meshing should reject invalid source model handles"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       mesh_handle == sqmesh::invalid_handle,
       "failed surface meshing should not allocate a mesh handle"
     )) {
    return EXIT_FAILURE;
  }

  if(!expect_error(
       sqmesh::mesh::create_surface_mesh(
         model,
         "",
         valid_auto_cfd_parameters,
         mesh_handle,
         context
       ),
       sqmesh::base::StatusCode::invalid_argument,
       "surface meshing should reject an empty algorithm name"
     )) {
    return EXIT_FAILURE;
  }

  if(!expect_error(
       sqmesh::mesh::create_surface_mesh(
         model,
         "Auto CFD Surface Mesher",
         valid_auto_cfd_parameters,
         mesh_handle,
         context
       ),
       sqmesh::base::StatusCode::invalid_argument,
       "Auto CFD surface meshing should reject placeholder models without real geometry faces"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       mesh_handle == sqmesh::invalid_handle,
       "failed Auto CFD surface meshing should not allocate a mesh handle"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::ParameterDictionary invalid_parameters;
  invalid_parameters.set_number("minimum_length", -1.0);
  if(!expect_error(
       sqmesh::mesh::create_volume_mesh(
         model,
         "Dummy Mesher",
         invalid_parameters,
         mesh_handle,
         context
       ),
       sqmesh::base::StatusCode::invalid_argument,
       "meshing should reject negative size parameters"
     )) {
    return EXIT_FAILURE;
  }

  if(!expect_error(
       sqmesh::mesh::create_surface_mesh(
         model,
         "Auto CFD Surface Mesher",
         invalid_parameters,
         mesh_handle,
         context
       ),
       sqmesh::base::StatusCode::invalid_argument,
       "surface meshing should reject negative size parameters"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::ParameterDictionary invalid_auto_cfd_range_parameters =
    valid_auto_cfd_parameters;
  invalid_auto_cfd_range_parameters.set_number("maximum_length", 0.05);
  if(!expect_error(
       sqmesh::mesh::create_surface_mesh(
         model,
         "Auto CFD Surface Mesher",
         invalid_auto_cfd_range_parameters,
         mesh_handle,
         context
       ),
       sqmesh::base::StatusCode::invalid_argument,
       "Auto CFD surface meshing should reject maximum_length values below minimum_length"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::ParameterDictionary invalid_auto_cfd_growth_parameters =
    valid_auto_cfd_parameters;
  invalid_auto_cfd_growth_parameters.set_number("growth_rate", 0.9);
  if(!expect_error(
       sqmesh::mesh::create_surface_mesh(
         model,
         "Auto CFD Surface Mesher",
         invalid_auto_cfd_growth_parameters,
         mesh_handle,
         context
       ),
       sqmesh::base::StatusCode::invalid_argument,
       "Auto CFD surface meshing should reject growth_rate values below 1.0"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::ParameterDictionary invalid_auto_cfd_proximity_parameters =
    valid_auto_cfd_parameters;
  invalid_auto_cfd_proximity_parameters.set_text("proximity", "yes");
  if(!expect_error(
       sqmesh::mesh::create_surface_mesh(
         model,
         "Auto CFD Surface Mesher",
         invalid_auto_cfd_proximity_parameters,
         mesh_handle,
         context
       ),
       sqmesh::base::StatusCode::invalid_argument,
       "Auto CFD surface meshing should reject non-boolean proximity parameters"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshingOptions invalid_local_surface_options;
  invalid_local_surface_options.parameters = valid_auto_cfd_parameters;
  invalid_local_surface_options.size_controls.add_local_size(
    {sqmesh::geo::TopologyDimension::vertex, 0U},
    0.25
  );
  if(!expect_error(
       sqmesh::mesh::create_surface_mesh(
         model,
         "Auto CFD Surface Mesher",
         invalid_local_surface_options,
         mesh_handle,
         context
       ),
       sqmesh::base::StatusCode::invalid_argument,
       "surface meshing should reject local size controls on unsupported topology dimensions"
     )) {
    return EXIT_FAILURE;
  }

  // Volume meshing local size control test removed — Simple Box Volume Mesher deleted.

  sqmesh::mesh::MeshSummary mesh_summary;
  if(!expect_error(
       sqmesh::mesh::mesh_summary(
         static_cast<sqmesh::mesh::MeshHandle>(bogus_handle),
         mesh_summary,
         context
       ),
       sqmesh::base::StatusCode::invalid_handle,
       "mesh queries should reject invalid mesh handles"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::Domain bogus_domain;
  if(!expect_error(
       sqmesh::mesh::domain_snapshot(
         static_cast<sqmesh::mesh::MeshHandle>(bogus_handle),
         bogus_domain,
         context
       ),
       sqmesh::base::StatusCode::invalid_handle,
       "domain_snapshot should reject invalid mesh handles"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       bogus_domain.entity_group_count() == 0U,
       "domain_snapshot should leave failed outputs in an empty state"
     )) {
    return EXIT_FAILURE;
  }

  std::size_t bogus_node_count = 123U;
  if(!expect_error(
       sqmesh::mesh::nodes_count(
         static_cast<sqmesh::mesh::MeshHandle>(bogus_handle),
         bogus_node_count,
         context
       ),
       sqmesh::base::StatusCode::invalid_handle,
       "nodes_count should reject invalid mesh handles"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       bogus_node_count == 0U,
       "nodes_count should leave failed count outputs in a stable zero state"
     )) {
    return EXIT_FAILURE;
  }

  std::size_t bogus_cell_count = 456U;
  if(!expect_error(
       sqmesh::mesh::cells_count(
         static_cast<sqmesh::mesh::MeshHandle>(bogus_handle),
         bogus_cell_count,
         context
       ),
       sqmesh::base::StatusCode::invalid_handle,
       "cells_count should reject invalid mesh handles"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       bogus_cell_count == 0U,
       "cells_count should leave failed count outputs in a stable zero state"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle io_mesh = sqmesh::invalid_handle;
  if(!expect_error(
       sqmesh::mesh::import_msh("", io_mesh, {}, context),
       sqmesh::base::StatusCode::invalid_argument,
       "MSH import should reject an empty path"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect_error(
       sqmesh::mesh::import_msh("missing.msh", io_mesh, {}, context),
       sqmesh::base::StatusCode::io_error,
       "MSH import should report io_error for a missing file"
     )) {
    return EXIT_FAILURE;
  }
  const auto malformed_binary_msh = std::filesystem::path("sqmesh_malformed_binary_input.msh");
  const auto binary_size_mismatch_msh = std::filesystem::path("sqmesh_binary_size_mismatch_input.msh");
  const auto unsupported_quad_msh = std::filesystem::path("sqmesh_quad_input.msh");
  const auto non_finite_binary_msh = std::filesystem::path("sqmesh_non_finite_node_input.msh");
  std::filesystem::remove(malformed_binary_msh);
  std::filesystem::remove(binary_size_mismatch_msh);
  std::filesystem::remove(unsupported_quad_msh);
  std::filesystem::remove(non_finite_binary_msh);
  if(!expect(
       write_text_file(
         malformed_binary_msh,
         "$MeshFormat\n2.2 1 8\n"
       ),
       "error_paths_smoke should create a malformed binary MSH fixture"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect_error(
       sqmesh::mesh::import_msh(malformed_binary_msh.string(), io_mesh, {}, context),
       sqmesh::base::StatusCode::io_error,
       "MSH import should reject malformed binary MeshFormat markers"
     )) {
    std::filesystem::remove(malformed_binary_msh);
    std::filesystem::remove(binary_size_mismatch_msh);
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::base::last_error_message().find("binary MeshFormat marker") !=
         std::string::npos,
       "malformed binary MSH rejection should explain the missing binary marker"
     )) {
    std::filesystem::remove(malformed_binary_msh);
    std::filesystem::remove(binary_size_mismatch_msh);
    return EXIT_FAILURE;
  }
  if(!expect(
       write_text_file(
         binary_size_mismatch_msh,
         "$MeshFormat\n4.1 1 16\n"
       ),
       "error_paths_smoke should create a binary size-mismatch MSH fixture"
     )) {
    std::filesystem::remove(malformed_binary_msh);
    return EXIT_FAILURE;
  }
  if(!expect_error(
       sqmesh::mesh::import_msh(binary_size_mismatch_msh.string(), io_mesh, {}, context),
       sqmesh::base::StatusCode::unsupported,
       "MSH import should reject binary Gmsh 4.1 files with mismatched MeshFormat data sizes"
     )) {
    std::filesystem::remove(malformed_binary_msh);
    std::filesystem::remove(binary_size_mismatch_msh);
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::base::last_error_message().find("4-byte or 8-byte") != std::string::npos,
       "binary Gmsh 4.1 rejection should explain the supported MeshFormat data-size boundary"
     )) {
    std::filesystem::remove(malformed_binary_msh);
    std::filesystem::remove(binary_size_mismatch_msh);
    std::filesystem::remove(unsupported_quad_msh);
    std::filesystem::remove(non_finite_binary_msh);
    return EXIT_FAILURE;
  }

  if(!expect(
       write_text_file(
         unsupported_quad_msh,
         "$MeshFormat\n2.2 0 8\n"
         "$EndMeshFormat\n"
         "$Nodes\n"
         "4\n"
         "1 0 0 0\n"
         "2 1 0 0\n"
         "3 1 1 0\n"
         "4 0 1 0\n"
         "$EndNodes\n"
         "$Elements\n"
         "1\n"
         "1 3 0 1 2 3 4\n"
         "$EndElements\n"
       ),
       "error_paths_smoke should create an MSH fixture with an unsupported quadrangle element"
     )) {
    std::filesystem::remove(malformed_binary_msh);
    std::filesystem::remove(binary_size_mismatch_msh);
    return EXIT_FAILURE;
  }
  if(!expect_error(
       sqmesh::mesh::import_msh(unsupported_quad_msh.string(), io_mesh, {}, context),
       sqmesh::base::StatusCode::unsupported,
       "MSH import should reject unsupported wider element families with a truthful diagnostic"
     )) {
    std::filesystem::remove(malformed_binary_msh);
    std::filesystem::remove(binary_size_mismatch_msh);
    std::filesystem::remove(unsupported_quad_msh);
    std::filesystem::remove(non_finite_binary_msh);
    return EXIT_FAILURE;
  }
  if(!expect_last_error_contains(
       "quadrangle (3)",
       "unsupported MSH element diagnostics should name the encountered element family"
     )) {
    std::filesystem::remove(malformed_binary_msh);
    std::filesystem::remove(binary_size_mismatch_msh);
    std::filesystem::remove(unsupported_quad_msh);
    std::filesystem::remove(non_finite_binary_msh);
    return EXIT_FAILURE;
  }

  if(!expect(
       write_non_finite_binary_msh_fixture(non_finite_binary_msh),
       "error_paths_smoke should create a binary MSH fixture with a non-finite node coordinate"
     )) {
    std::filesystem::remove(malformed_binary_msh);
    std::filesystem::remove(binary_size_mismatch_msh);
    std::filesystem::remove(unsupported_quad_msh);
    return EXIT_FAILURE;
  }
  if(!expect_error(
       sqmesh::mesh::import_msh(non_finite_binary_msh.string(), io_mesh, {}, context),
       sqmesh::base::StatusCode::io_error,
       "MSH import should reject non-finite node coordinates on the supported subset"
     )) {
    std::filesystem::remove(malformed_binary_msh);
    std::filesystem::remove(binary_size_mismatch_msh);
    std::filesystem::remove(unsupported_quad_msh);
    std::filesystem::remove(non_finite_binary_msh);
    return EXIT_FAILURE;
  }
  if(!expect_last_error_contains(
       "non-finite coordinate",
       "non-finite MSH node diagnostics should explain the coordinate validity failure"
     )) {
    std::filesystem::remove(malformed_binary_msh);
    std::filesystem::remove(binary_size_mismatch_msh);
    std::filesystem::remove(unsupported_quad_msh);
    std::filesystem::remove(non_finite_binary_msh);
    return EXIT_FAILURE;
  }

  std::filesystem::remove(malformed_binary_msh);
  std::filesystem::remove(binary_size_mismatch_msh);
  std::filesystem::remove(unsupported_quad_msh);
  std::filesystem::remove(non_finite_binary_msh);
  if(!expect_error(
       sqmesh::mesh::import_obj("", io_mesh, {}, context),
       sqmesh::base::StatusCode::invalid_argument,
       "OBJ import should reject an empty path"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect_error(
       sqmesh::mesh::import_obj("missing.obj", io_mesh, {}, context),
       sqmesh::base::StatusCode::io_error,
       "OBJ import should report io_error for a missing file"
     )) {
    return EXIT_FAILURE;
  }

  const auto malformed_obj = std::filesystem::path("sqmesh_obj_missing_vertex.obj");
  const auto unsupported_obj = std::filesystem::path("sqmesh_obj_unsupported_record.obj");
  const auto self_intersecting_obj =
    std::filesystem::path("sqmesh_obj_self_intersecting_face.obj");
  const auto degenerate_obj = std::filesystem::path("sqmesh_obj_degenerate_face.obj");
  const auto non_planar_obj = std::filesystem::path("sqmesh_obj_non_planar_face.obj");
  const auto repeated_vertex_obj = std::filesystem::path("sqmesh_obj_repeated_vertex_face.obj");
  std::filesystem::remove(malformed_obj);
  std::filesystem::remove(unsupported_obj);
  std::filesystem::remove(self_intersecting_obj);
  std::filesystem::remove(degenerate_obj);
  std::filesystem::remove(non_planar_obj);
  std::filesystem::remove(repeated_vertex_obj);
  if(!expect(
       write_text_file(
         malformed_obj,
         "v 0 0 0\n"
         "v 1 0 0\n"
         "f 1 2 3\n"
       ),
       "error_paths_smoke should create an OBJ fixture with an undefined face vertex"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect_error(
       sqmesh::mesh::import_obj(malformed_obj.string(), io_mesh, {}, context),
       sqmesh::base::StatusCode::io_error,
       "OBJ import should reject face records that reference undefined vertices"
     )) {
    std::filesystem::remove(malformed_obj);
    std::filesystem::remove(unsupported_obj);
    std::filesystem::remove(non_planar_obj);
    std::filesystem::remove(repeated_vertex_obj);
    return EXIT_FAILURE;
  }
  if(!expect_last_error_contains(
       "references a vertex that has not been defined",
       "OBJ face-reference failures should explain that the referenced vertex is undefined"
     )) {
    std::filesystem::remove(malformed_obj);
    std::filesystem::remove(unsupported_obj);
    std::filesystem::remove(non_planar_obj);
    std::filesystem::remove(repeated_vertex_obj);
    return EXIT_FAILURE;
  }
  if(!expect(
       write_text_file(
         unsupported_obj,
         "v 0 0 0\n"
         "curv 0 1 1\n"
       ),
       "error_paths_smoke should create an OBJ fixture with an unsupported record"
     )) {
    std::filesystem::remove(malformed_obj);
    std::filesystem::remove(self_intersecting_obj);
    std::filesystem::remove(degenerate_obj);
    std::filesystem::remove(non_planar_obj);
    std::filesystem::remove(repeated_vertex_obj);
    return EXIT_FAILURE;
  }
  if(!expect_error(
       sqmesh::mesh::import_obj(unsupported_obj.string(), io_mesh, {}, context),
       sqmesh::base::StatusCode::unsupported,
       "OBJ import should reject unsupported records cleanly"
     )) {
    std::filesystem::remove(malformed_obj);
    std::filesystem::remove(unsupported_obj);
    std::filesystem::remove(self_intersecting_obj);
    std::filesystem::remove(degenerate_obj);
    std::filesystem::remove(non_planar_obj);
    std::filesystem::remove(repeated_vertex_obj);
    return EXIT_FAILURE;
  }
  if(!expect_last_error_contains(
       "`curv`",
       "unsupported OBJ-record diagnostics should name the offending keyword"
     )) {
    std::filesystem::remove(malformed_obj);
    std::filesystem::remove(unsupported_obj);
    std::filesystem::remove(self_intersecting_obj);
    std::filesystem::remove(degenerate_obj);
    std::filesystem::remove(non_planar_obj);
    std::filesystem::remove(repeated_vertex_obj);
    return EXIT_FAILURE;
  }

  if(!expect(
       write_text_file(
         degenerate_obj,
         "v 0 0 0\n"
         "v 1 0 0\n"
         "v 2 0 0\n"
         "f 1 2 3\n"
       ),
       "error_paths_smoke should create an OBJ fixture with a degenerate polygon face"
     )) {
    std::filesystem::remove(malformed_obj);
    std::filesystem::remove(unsupported_obj);
    std::filesystem::remove(self_intersecting_obj);
    std::filesystem::remove(degenerate_obj);
    std::filesystem::remove(non_planar_obj);
    std::filesystem::remove(repeated_vertex_obj);
    return EXIT_FAILURE;
  }
  if(!expect_error(
       sqmesh::mesh::import_obj(degenerate_obj.string(), io_mesh, {}, context),
       sqmesh::base::StatusCode::unsupported,
       "OBJ import should reject degenerate polygon faces on the bounded simple-planar subset"
     )) {
    std::filesystem::remove(malformed_obj);
    std::filesystem::remove(unsupported_obj);
    std::filesystem::remove(self_intersecting_obj);
    std::filesystem::remove(degenerate_obj);
    std::filesystem::remove(non_planar_obj);
    std::filesystem::remove(repeated_vertex_obj);
    return EXIT_FAILURE;
  }
  if(!expect_last_error_contains(
       "degenerate",
       "degenerate OBJ polygon diagnostics should explain the non-degenerate boundary"
     )) {
    std::filesystem::remove(malformed_obj);
    std::filesystem::remove(unsupported_obj);
    std::filesystem::remove(self_intersecting_obj);
    std::filesystem::remove(degenerate_obj);
    std::filesystem::remove(non_planar_obj);
    std::filesystem::remove(repeated_vertex_obj);
    return EXIT_FAILURE;
  }

  if(!expect(
       write_text_file(
         non_planar_obj,
         "v 0 0 0\n"
         "v 1 0 0\n"
         "v 1 1 1\n"
         "v 0 1 0\n"
         "f 1 2 3 4\n"
       ),
       "error_paths_smoke should create an OBJ fixture with a non-planar polygon face"
     )) {
    std::filesystem::remove(malformed_obj);
    std::filesystem::remove(unsupported_obj);
    std::filesystem::remove(self_intersecting_obj);
    std::filesystem::remove(degenerate_obj);
    std::filesystem::remove(repeated_vertex_obj);
    return EXIT_FAILURE;
  }
  if(!expect_error(
       sqmesh::mesh::import_obj(non_planar_obj.string(), io_mesh, {}, context),
       sqmesh::base::StatusCode::unsupported,
       "OBJ import should reject non-planar polygon faces on the bounded simple-planar subset"
     )) {
    std::filesystem::remove(malformed_obj);
    std::filesystem::remove(unsupported_obj);
    std::filesystem::remove(self_intersecting_obj);
    std::filesystem::remove(degenerate_obj);
    std::filesystem::remove(non_planar_obj);
    std::filesystem::remove(repeated_vertex_obj);
    return EXIT_FAILURE;
  }
  if(!expect_last_error_contains(
       "not coplanar",
       "non-planar OBJ polygon diagnostics should explain the planar-face boundary"
     )) {
    std::filesystem::remove(malformed_obj);
    std::filesystem::remove(unsupported_obj);
    std::filesystem::remove(self_intersecting_obj);
    std::filesystem::remove(degenerate_obj);
    std::filesystem::remove(non_planar_obj);
    std::filesystem::remove(repeated_vertex_obj);
    return EXIT_FAILURE;
  }

  if(!expect(
       write_text_file(
         repeated_vertex_obj,
         "v 0 0 0\n"
         "v 1 0 0\n"
         "v 1 1 0\n"
         "v 0 1 0\n"
         "f 1 2 3 2\n"
       ),
       "error_paths_smoke should create an OBJ fixture with repeated polygon vertex references"
     )) {
    std::filesystem::remove(malformed_obj);
    std::filesystem::remove(unsupported_obj);
    std::filesystem::remove(self_intersecting_obj);
    std::filesystem::remove(degenerate_obj);
    std::filesystem::remove(non_planar_obj);
    return EXIT_FAILURE;
  }
  if(!expect_error(
       sqmesh::mesh::import_obj(repeated_vertex_obj.string(), io_mesh, {}, context),
       sqmesh::base::StatusCode::unsupported,
       "OBJ import should reject repeated polygon vertex references instead of silently triangulating them"
     )) {
    std::filesystem::remove(malformed_obj);
    std::filesystem::remove(unsupported_obj);
    std::filesystem::remove(self_intersecting_obj);
    std::filesystem::remove(degenerate_obj);
    std::filesystem::remove(non_planar_obj);
    std::filesystem::remove(repeated_vertex_obj);
    return EXIT_FAILURE;
  }
  if(!expect_last_error_contains(
       "repeated vertex references",
       "repeated-vertex OBJ diagnostics should explain the simple-polygon boundary"
     )) {
    std::filesystem::remove(malformed_obj);
    std::filesystem::remove(unsupported_obj);
    std::filesystem::remove(self_intersecting_obj);
    std::filesystem::remove(degenerate_obj);
    std::filesystem::remove(non_planar_obj);
    std::filesystem::remove(repeated_vertex_obj);
    return EXIT_FAILURE;
  }

  if(!expect(
       write_text_file(
         self_intersecting_obj,
         "v 0 0 0\n"
         "v 2 0 0\n"
         "v 0 2 0\n"
         "v 2 1 0\n"
         "v 0 1 0\n"
         "f 1 2 3 4 5\n"
       ),
       "error_paths_smoke should create an OBJ fixture with a self-intersecting polygon face"
     )) {
    std::filesystem::remove(malformed_obj);
    std::filesystem::remove(unsupported_obj);
    std::filesystem::remove(degenerate_obj);
    return EXIT_FAILURE;
  }
  if(!expect_error(
       sqmesh::mesh::import_obj(self_intersecting_obj.string(), io_mesh, {}, context),
       sqmesh::base::StatusCode::unsupported,
       "OBJ import should reject self-intersecting polygon faces instead of over-claiming polygon coverage"
     )) {
    std::filesystem::remove(malformed_obj);
    std::filesystem::remove(unsupported_obj);
    std::filesystem::remove(self_intersecting_obj);
    std::filesystem::remove(degenerate_obj);
    return EXIT_FAILURE;
  }
  if(!expect_last_error_contains(
       "self-intersecting",
       "self-intersecting OBJ polygon diagnostics should explain the simple-planar boundary"
     )) {
    std::filesystem::remove(malformed_obj);
    std::filesystem::remove(unsupported_obj);
    std::filesystem::remove(self_intersecting_obj);
    std::filesystem::remove(degenerate_obj);
    return EXIT_FAILURE;
  }
  std::filesystem::remove(malformed_obj);
  std::filesystem::remove(unsupported_obj);
  std::filesystem::remove(self_intersecting_obj);
  std::filesystem::remove(degenerate_obj);
  std::filesystem::remove(non_planar_obj);
  std::filesystem::remove(repeated_vertex_obj);

  const auto invalid_msh_export = std::filesystem::path("sqmesh_invalid_export.msh");
  const auto invalid_obj_export = std::filesystem::path("sqmesh_invalid_export.obj");
  std::filesystem::remove(invalid_msh_export);
  std::filesystem::remove(invalid_obj_export);
  if(!expect_error(
       sqmesh::mesh::export_msh(
         static_cast<sqmesh::mesh::MeshHandle>(bogus_handle),
         invalid_msh_export.string(),
         {},
         context
       ),
       sqmesh::base::StatusCode::invalid_handle,
       "MSH export should reject invalid mesh handles"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect_error(
       sqmesh::mesh::export_obj(
         static_cast<sqmesh::mesh::MeshHandle>(bogus_handle),
         invalid_obj_export.string(),
         {},
         context
       ),
       sqmesh::base::StatusCode::invalid_handle,
       "OBJ export should reject invalid mesh handles"
     )) {
    return EXIT_FAILURE;
  }
  sqmesh::mesh::MeshQualityReport invalid_quality_report;
  if(!expect_error(
       sqmesh::mesh::mesh_quality_report(
         static_cast<sqmesh::mesh::MeshHandle>(bogus_handle),
         invalid_quality_report,
         context
       ),
       sqmesh::base::StatusCode::invalid_handle,
       "mesh_quality_report should reject invalid mesh handles"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::ModelHandle io_model = sqmesh::invalid_handle;
  if(!expect_error(
       sqmesh::geo::import_stl("", io_model, {}, context),
       sqmesh::base::StatusCode::invalid_argument,
       "STL import should reject an empty path"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect_error(
       sqmesh::geo::import_stl("missing.stl", io_model, {}, context),
       sqmesh::base::StatusCode::io_error,
       "STL import should report io_error for a missing file"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::StlImportOptions invalid_stl_options;
  invalid_stl_options.relative_merge_tolerance = -1.0;
  if(!expect_error(
       sqmesh::geo::import_stl("missing.stl", io_model, invalid_stl_options, context),
       sqmesh::base::StatusCode::invalid_argument,
       "STL import should reject a negative relative merge tolerance"
     )) {
    return EXIT_FAILURE;
  }

  const auto malformed_ascii_stl = std::filesystem::path("sqmesh_malformed_ascii_input.stl");
  const auto truncated_binary_stl = std::filesystem::path("sqmesh_truncated_binary_input.stl");
  std::filesystem::remove(malformed_ascii_stl);
  std::filesystem::remove(truncated_binary_stl);
  if(!expect(
       write_text_file(
         malformed_ascii_stl,
         "solid broken\n"
         "facet normal 0 0 1\n"
         "outer loop\n"
         "vertex 0 0 0\n"
         "vertex 1 0 0\n"
         "vertex 0 1\n"
         "endloop\n"
         "endfacet\n"
         "endsolid broken\n"
       ),
       "error_paths_smoke should create a malformed ASCII STL fixture"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect_error(
       sqmesh::geo::import_stl(malformed_ascii_stl.string(), io_model, {}, context),
       sqmesh::base::StatusCode::io_error,
       "STL import should reject malformed ASCII facet vertex records"
     )) {
    std::filesystem::remove(malformed_ascii_stl);
    std::filesystem::remove(truncated_binary_stl);
    return EXIT_FAILURE;
  }
  if(!expect_last_error_contains(
       "vertex x y z",
       "malformed ASCII STL diagnostics should explain the expected vertex record shape"
     )) {
    std::filesystem::remove(malformed_ascii_stl);
    std::filesystem::remove(truncated_binary_stl);
    return EXIT_FAILURE;
  }

  std::string truncated_binary_bytes(84U, '\0');
  truncated_binary_bytes[0] = 'B';
  truncated_binary_bytes[1] = 'I';
  truncated_binary_bytes[2] = 'N';
  truncated_binary_bytes[3] = 'A';
  truncated_binary_bytes[4] = 'R';
  truncated_binary_bytes[5] = 'Y';
  truncated_binary_bytes[80] = 1;
  if(!expect(
       write_binary_file(truncated_binary_stl, truncated_binary_bytes),
       "error_paths_smoke should create a truncated binary STL fixture"
     )) {
    std::filesystem::remove(malformed_ascii_stl);
    return EXIT_FAILURE;
  }
  if(!expect_error(
       sqmesh::geo::import_stl(truncated_binary_stl.string(), io_model, {}, context),
       sqmesh::base::StatusCode::io_error,
       "STL import should reject binary facet tables that are shorter than their declared facet count"
     )) {
    std::filesystem::remove(malformed_ascii_stl);
    std::filesystem::remove(truncated_binary_stl);
    return EXIT_FAILURE;
  }
  if(!expect_last_error_contains(
       "expected at least 134 bytes",
       "truncated binary STL diagnostics should report the declared size mismatch"
     )) {
    std::filesystem::remove(malformed_ascii_stl);
    std::filesystem::remove(truncated_binary_stl);
    return EXIT_FAILURE;
  }
  std::filesystem::remove(malformed_ascii_stl);
  std::filesystem::remove(truncated_binary_stl);

  if(sqmesh::geo::cad_io_available()) {
    if(!expect_error(
         sqmesh::geo::import_step("", io_model, {}, context),
         sqmesh::base::StatusCode::invalid_argument,
         "STEP import should reject an empty path"
       )) {
      return EXIT_FAILURE;
    }
    if(!expect_error(
         sqmesh::geo::import_step("missing.step", io_model, {}, context),
         sqmesh::base::StatusCode::io_error,
         "STEP import should report io_error for a missing file"
       )) {
      return EXIT_FAILURE;
    }
    if(!expect_error(
         sqmesh::geo::import_iges("missing.igs", io_model, {}, context),
         sqmesh::base::StatusCode::io_error,
         "IGES import should report io_error for a missing file"
       )) {
      return EXIT_FAILURE;
    }

    const auto invalid_export_path = std::filesystem::path("sqmesh_invalid_export.step");
    std::filesystem::remove(invalid_export_path);
    if(!expect_error(
         sqmesh::geo::export_step(
           static_cast<sqmesh::geo::ModelHandle>(bogus_handle),
           invalid_export_path.string(),
           {},
           context
         ),
         sqmesh::base::StatusCode::invalid_handle,
         "STEP export should reject invalid model handles"
       )) {
      return EXIT_FAILURE;
    }
  } else {
    if(!expect_error(
         sqmesh::geo::import_step("missing.step", io_model, {}, context),
         sqmesh::base::StatusCode::unsupported,
         "STEP import should remain unsupported without OCC"
       )) {
      return EXIT_FAILURE;
    }
  }

  if(!expect_error(
       sqmesh::base::shutdown_all(),
       sqmesh::base::StatusCode::ok,
       "shutdown_all should clean up the runtime"
     )) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
