// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "sqmesh/geo/api.hpp"

#include <cstdlib>

int main()
{
  const sqmesh::geo::StepImportOptions step_import_options {};
  const sqmesh::geo::StlImportOptions stl_import_options {};
  const sqmesh::geo::IgesExportOptions iges_export_options {};
  const sqmesh::geo::ModelSummary model_summary {};
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
  const sqmesh::geo::VertexView vertex_view {};
  const sqmesh::geo::EdgeView edge_view {};
  const sqmesh::geo::FaceView face_view {};
  const sqmesh::geo::RegionView region_view {};
  const sqmesh::geo::ModelView geometry_view {};
  const sqmesh::geo::FeatureEdgeOptions feature_options {};
  const sqmesh::geo::FeatureEdgeReport feature_report {};

  if(step_import_options.system_length_unit != 0.0) {
    return EXIT_FAILURE;
  }
  if(stl_import_options.relative_merge_tolerance != 1.0e-9) {
    return EXIT_FAILURE;
  }
  if(!iges_export_options.unit_name.empty()) {
    return EXIT_FAILURE;
  }
  if(model_summary.entity_count != 0U) {
    return EXIT_FAILURE;
  }
  if(sqmesh::geo::is_valid(topology_entity)) {
    return EXIT_FAILURE;
  }
  if(topology_snapshot.topology_revision != 0U ||
     topology_snapshot.entity_count(sqmesh::geo::TopologyDimension::face) != 0U) {
    return EXIT_FAILURE;
  }
  if(face_bounds.u_min != 0.0 || face_bounds.v_max != 0.0) {
    return EXIT_FAILURE;
  }
  if(face_sample.normal_defined || face_curvature.curvature_defined ||
     face_derivatives.first_derivatives_defined ||
     face_derivatives.second_derivatives_defined || face_derivatives.normal_defined ||
     face_projection.distance != 0.0 || face_projection.normal_defined ||
     face_uv_mapping.distance != 0.0 || edge_info.approximate_length != 0.0 ||
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
  if(sqmesh::geo::is_valid(vertex_view.entity) || !vertex_view.edge_ids.empty() ||
     sqmesh::geo::is_valid(edge_view.entity) || !edge_view.face_ids.empty() ||
     !edge_view.vertex_ids.empty() || sqmesh::geo::is_valid(face_view.entity) ||
     !face_view.region_ids.empty() || !face_view.edge_ids.empty() ||
     !sqmesh::geo::ordered_boundary_loops(face_view).loops.empty() ||
     sqmesh::geo::is_valid(region_view.entity) || !region_view.face_ids.empty() ||
     geometry_view.entity_count(sqmesh::geo::TopologyDimension::face) != 0U ||
     geometry_view.find_face({sqmesh::geo::TopologyDimension::face, 0U}) != nullptr ||
     !geometry_view.root_face_ids.empty() || !geometry_view.root_edge_ids.empty() ||
     !geometry_view.root_vertex_ids.empty()) {
    return EXIT_FAILURE;
  }
  if(feature_options.feature_angle_degrees != 30.0 ||
     !feature_options.include_boundary_edges ||
     !feature_options.include_non_manifold_edges ||
     !feature_report.edges.empty()) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
