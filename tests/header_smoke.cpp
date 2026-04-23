// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "sqmesh/sqmesh.hpp"

#include <cstdlib>
#include <string_view>

int main()
{
  const auto model_summary = sqmesh::geo::placeholder_model_summary();
  const sqmesh::geo::StepImportOptions step_import_options {};
  const sqmesh::geo::StlImportOptions stl_import_options {};
  const sqmesh::geo::IgesImportOptions iges_import_options {};
  const sqmesh::geo::StepExportOptions step_export_options {};
  const sqmesh::geo::IgesExportOptions iges_export_options {};
  const sqmesh::geo::TopologyCheckReport topology_report {};
  const sqmesh::geo::TopoOptions topo_options {};
  const sqmesh::geo::TopoReport topo_report {};
  const sqmesh::geo::TopologyEntityId topology_entity {};
  const sqmesh::geo::TopologySnapshot topology_snapshot {};
  const sqmesh::geo::FaceUvBounds face_bounds {};
  const sqmesh::geo::FaceSample face_sample {};
  const sqmesh::geo::FaceCurvatureSample face_curvature {};
  const sqmesh::geo::FaceDerivatives face_derivatives {};
  const sqmesh::geo::FaceProjection face_projection {};
  const sqmesh::geo::FaceUvMapping face_uv_mapping {};
  const sqmesh::geo::EdgeCurveInfo edge_info {};
  const sqmesh::geo::EdgeTangentSample edge_tangent {};
  const sqmesh::geo::EdgeCurveSamplingOptions edge_curve_sampling_options {};
  const sqmesh::geo::EdgeCurveSamples edge_curve_samples {};
  const sqmesh::geo::FaceBoundaryEdgeUse face_boundary_edge {};
  const sqmesh::geo::FaceBoundaryLoop face_boundary_loop {};
  const sqmesh::geo::FaceBoundaryLoops face_boundary {};
  const sqmesh::geo::FeatureEdgeOptions feature_options {};
  const sqmesh::geo::FeatureEdgeReport feature_report {};
  const auto mesh_domain = sqmesh::mesh::make_dummy_tetra_domain();
  const sqmesh::mesh::MeshSizeControl mesh_size_control {};
  const sqmesh::mesh::MeshSizeControls mesh_size_controls {};
  const sqmesh::mesh::MeshingOptions meshing_options {};
  const sqmesh::mesh::Edge edge {};
  const auto mesh_summary = mesh_domain.summary();
  const auto mesh_stats = mesh_domain.statistics();

#if defined(_TopoDS_Shape_HeaderFile) || defined(_STEPControl_Reader_HeaderFile) || \
  defined(_IGESControl_Reader_HeaderFile)
  return EXIT_FAILURE;
#endif

  if(sqmesh::version().major_v != 0) {
    return EXIT_FAILURE;
  }
  if(sqmesh::version().minor_v != 1) {
    return EXIT_FAILURE;
  }
  if(sqmesh::version().patch_v != 0) {
    return EXIT_FAILURE;
  }
  if(sqmesh::version_string()[0] == '\0') {
    return EXIT_FAILURE;
  }
  if(sqmesh::base::module_name() != std::string_view("sqmesh::base")) {
    return EXIT_FAILURE;
  }
  if(sqmesh::geo::module_name() != std::string_view("sqmesh::geo")) {
    return EXIT_FAILURE;
  }
  if(step_import_options.system_length_unit != 0.0) {
    return EXIT_FAILURE;
  }
  if(stl_import_options.relative_merge_tolerance != 1.0e-9) {
    return EXIT_FAILURE;
  }
  if(iges_import_options.read_visible_only) {
    return EXIT_FAILURE;
  }
  if(step_export_options.linear_tolerance != 0.0) {
    return EXIT_FAILURE;
  }
  if(!step_export_options.unit_name.empty() || !step_export_options.schema.empty()) {
    return EXIT_FAILURE;
  }
  if(!iges_export_options.unit_name.empty()) {
    return EXIT_FAILURE;
  }
  if(iges_export_options.write_mode != sqmesh::geo::IgesWriteMode::brep) {
    return EXIT_FAILURE;
  }
  if(topology_report.is_valid || topology_report.free_edge_count != 0U ||
     topology_report.contiguous_edge_count != 0U ||
     topology_report.multiple_edge_count != 0U) {
    return EXIT_FAILURE;
  }
  if(topology_report.model_summary.entity_count != 0U) {
    return EXIT_FAILURE;
  }
  if(topo_options.tolerance != 1.0e-6 || topo_options.min_tolerance != 0.0 ||
     topo_options.max_tolerance != 0.0) {
    return EXIT_FAILURE;
  }
  if(!topo_options.fix_degenerated || !topo_options.fix_small_edges ||
     !topo_options.fix_small_faces || !topo_options.sew_faces ||
     !topo_options.make_solids) {
    return EXIT_FAILURE;
  }
  if(topo_report.topology_revision_before != 0U ||
     topo_report.topology_revision_after != 0U ||
     topo_report.modified || topo_report.topology_identity_changed ||
     topo_report.free_edges_reduced ||
     topo_report.sewing_performed || topo_report.sewing_modified ||
     topo_report.shape_fix_performed || topo_report.shape_fix_modified) {
    return EXIT_FAILURE;
  }
  if(topo_report.before.free_edge_count != 0U || topo_report.after.free_edge_count != 0U) {
    return EXIT_FAILURE;
  }
  if(sqmesh::geo::is_valid(topology_entity)) {
    return EXIT_FAILURE;
  }
  if(topology_snapshot.topology_revision != 0U ||
     topology_snapshot.entity_count(sqmesh::geo::TopologyDimension::region) != 0U ||
     topology_snapshot.entity_count(sqmesh::geo::TopologyDimension::vertex) != 0U) {
    return EXIT_FAILURE;
  }
  if(face_bounds.u_max != 0.0 || face_sample.normal_defined ||
     face_curvature.curvature_defined || face_derivatives.first_derivatives_defined ||
     face_derivatives.second_derivatives_defined || face_derivatives.normal_defined ||
     face_projection.distance != 0.0 || face_projection.normal_defined ||
     face_uv_mapping.distance != 0.0 || edge_info.parameter_max != 0.0 ||
     edge_tangent.speed != 0.0 || edge_tangent.tangent_defined ||
     edge_curve_sampling_options.target_segment_length != 0.0 ||
     edge_curve_sampling_options.min_segment_count != 1U ||
     !edge_curve_samples.samples.empty()) {
    return EXIT_FAILURE;
  }
  if(sqmesh::geo::is_valid(face_boundary_edge.edge) ||
     face_boundary_edge.same_orientation_as_edge || face_boundary_edge.is_seam ||
     face_boundary_edge.is_degenerate || face_boundary_loop.closed ||
     face_boundary_loop.continuous ||
     face_boundary_loop.seam_edge_use_count != 0U ||
     face_boundary_loop.degenerate_edge_use_count != 0U ||
     face_boundary_loop.repeated_edge_use_count != 0U ||
     !face_boundary_loop.edge_uses.empty() || !face_boundary_loop.vertex_ids.empty() ||
     sqmesh::geo::is_valid(face_boundary.face) || !face_boundary.loops.empty() ||
     face_boundary.primary_outer_loop_index != sqmesh::geo::invalid_boundary_loop_index ||
     face_boundary.outer_loop_count != 0U || face_boundary.inner_loop_count != 0U ||
     face_boundary.unknown_loop_count != 0U || face_boundary.closed_loop_count != 0U ||
     face_boundary.open_loop_count != 0U ||
     face_boundary.non_continuous_loop_count != 0U ||
     face_boundary.seam_loop_count != 0U ||
     face_boundary.seam_edge_use_count != 0U ||
     face_boundary.degenerate_loop_count != 0U ||
     face_boundary.degenerate_edge_use_count != 0U ||
     face_boundary.repeated_edge_loop_count != 0U ||
     face_boundary.repeated_edge_use_count != 0U ||
     face_boundary.has_holes || face_boundary.has_seams ||
     sqmesh::geo::primary_outer_boundary_loop(face_boundary) != nullptr) {
    return EXIT_FAILURE;
  }
  if(feature_options.feature_angle_degrees != 30.0 ||
     feature_report.sharp_edge_count != 0U) {
    return EXIT_FAILURE;
  }
  if(sqmesh::mesh::module_name() != std::string_view("sqmesh::mesh")) {
    return EXIT_FAILURE;
  }
  if(sqmesh::geo::is_valid(mesh_size_control.entity) ||
     mesh_size_control.target_size != 0.0 ||
     mesh_size_controls.topology_revision != 0U ||
     !mesh_size_controls.empty() ||
     meshing_options.parameters.size() != 0U ||
     !meshing_options.size_controls.empty()) {
    return EXIT_FAILURE;
  }
  if(model_summary.entity_count != 0U || edge.header.id != 0U) {
    return EXIT_FAILURE;
  }
  if(mesh_summary.node_count != 4U || mesh_summary.edge_count != 0U ||
     mesh_summary.face_count != 4U ||
     mesh_summary.cell_count != 1U) {
    return EXIT_FAILURE;
  }
  if(mesh_stats.entity_group_count != 3U || mesh_stats.boundary_edge_count != 0U ||
     mesh_stats.boundary_face_count != 4U || mesh_stats.line_edge_count != 0U ||
     mesh_stats.tetra_cell_count != 1U) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
