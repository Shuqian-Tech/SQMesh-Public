// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "sqmesh/sqmesh.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#if defined(SQMESH_TEST_OCC_ENABLED)
#include "../src/cad/occ/adapter.hpp"

#include <BRepPrimAPI_MakeBox.hxx>
#include <TopoDS_Shape.hxx>
#endif

namespace {

bool expect(bool condition, const char *message)
{
  if(condition) {
    return true;
  }

  std::fprintf(stderr, "geo_topology_query_smoke: %s\n", message);
  return false;
}

double vector_norm(const sqmesh::geo::Vector3 &vector)
{
  return std::sqrt(
    vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]
  );
}

double point_distance(
  const sqmesh::geo::Point3 &lhs,
  const sqmesh::geo::Point3 &rhs
)
{
  return std::sqrt(
    (lhs[0] - rhs[0]) * (lhs[0] - rhs[0]) +
    (lhs[1] - rhs[1]) * (lhs[1] - rhs[1]) +
    (lhs[2] - rhs[2]) * (lhs[2] - rhs[2])
  );
}

double dot_product(
  const sqmesh::geo::Vector3 &lhs,
  const sqmesh::geo::Vector3 &rhs
)
{
  return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2];
}

sqmesh::geo::Vector3 normalized(const sqmesh::geo::Vector3 &vector)
{
  const double length = vector_norm(vector);
  if(length <= 0.0) {
    return {0.0, 0.0, 0.0};
  }

  return {vector[0] / length, vector[1] / length, vector[2] / length};
}

sqmesh::geo::Point3 translated_point(
  const sqmesh::geo::Point3 &point,
  const sqmesh::geo::Vector3 &direction,
  double distance
)
{
  return {
    point[0] + direction[0] * distance,
    point[1] + direction[1] * distance,
    point[2] + direction[2] * distance,
  };
}

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

    sqmesh::geo::TopologySnapshot snapshot;
    if(!expect(
         sqmesh::geo::topology_snapshot(model_handle, snapshot, context) ==
           sqmesh::base::StatusCode::ok,
         "topology_snapshot should succeed for placeholder models"
       )) {
      return EXIT_FAILURE;
    }
    if(!expect(
         snapshot.entity_count(sqmesh::geo::TopologyDimension::region) == 0U &&
           snapshot.entity_count(sqmesh::geo::TopologyDimension::face) == 0U &&
           snapshot.entity_count(sqmesh::geo::TopologyDimension::edge) == 0U &&
           snapshot.entity_count(sqmesh::geo::TopologyDimension::vertex) == 0U,
         "placeholder models should expose an empty topology snapshot"
       )) {
      return EXIT_FAILURE;
    }

    std::vector<sqmesh::geo::TopologyEntityId> children;
    if(!expect(
         sqmesh::geo::topology_children(
           model_handle,
           {sqmesh::geo::TopologyDimension::face, 0U},
           children,
           context
         ) == sqmesh::base::StatusCode::invalid_argument,
         "placeholder models should reject entity traversal requests"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::geo::FaceUvBounds face_bounds;
    if(!expect(
         sqmesh::geo::face_uv_bounds(
           model_handle,
           {sqmesh::geo::TopologyDimension::face, 0U},
           face_bounds,
           context
         ) == sqmesh::base::StatusCode::invalid_argument,
         "placeholder models should reject face geometry requests"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::geo::FaceSample face_sample;
    if(!expect(
         sqmesh::geo::sample_face(
           model_handle,
           {sqmesh::geo::TopologyDimension::face, 0U},
           0.0,
           0.0,
           face_sample,
           context
         ) == sqmesh::base::StatusCode::invalid_argument,
         "placeholder models should reject face sampling requests"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::geo::FaceCurvatureSample face_curvature;
    if(!expect(
         sqmesh::geo::sample_face_curvature(
           model_handle,
           {sqmesh::geo::TopologyDimension::face, 0U},
           0.0,
           0.0,
           face_curvature,
           context
         ) == sqmesh::base::StatusCode::invalid_argument,
         "placeholder models should reject face curvature requests"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::geo::FaceDerivatives face_derivatives;
    if(!expect(
         sqmesh::geo::sample_face_derivatives(
           model_handle,
           {sqmesh::geo::TopologyDimension::face, 0U},
           0.0,
           0.0,
           face_derivatives,
           context
         ) == sqmesh::base::StatusCode::invalid_argument,
         "placeholder models should reject face derivative requests"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::geo::FaceProjection face_projection;
    if(!expect(
         sqmesh::geo::project_point_to_face(
           model_handle,
           {sqmesh::geo::TopologyDimension::face, 0U},
           {0.0, 0.0, 0.0},
           face_projection,
           context
         ) == sqmesh::base::StatusCode::invalid_argument,
         "placeholder models should reject face projection requests"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::geo::FaceUvMapping face_uv_mapping;
    if(!expect(
         sqmesh::geo::recover_face_uv(
           model_handle,
           {sqmesh::geo::TopologyDimension::face, 0U},
           {0.0, 0.0, 0.0},
           face_uv_mapping,
           context
         ) == sqmesh::base::StatusCode::invalid_argument,
         "placeholder models should reject face inverse-mapping requests"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::geo::EdgeCurveInfo edge_info;
    if(!expect(
         sqmesh::geo::edge_curve_info(
           model_handle,
           {sqmesh::geo::TopologyDimension::edge, 0U},
           edge_info,
           context
         ) == sqmesh::base::StatusCode::invalid_argument,
         "placeholder models should reject edge geometry requests"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::geo::EdgeTangentSample edge_tangent;
    if(!expect(
         sqmesh::geo::sample_edge_tangent(
           model_handle,
           {sqmesh::geo::TopologyDimension::edge, 0U},
           0.0,
           edge_tangent,
           context
         ) == sqmesh::base::StatusCode::invalid_argument,
         "placeholder models should reject edge tangent requests"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::geo::FaceBoundaryLoops face_boundary;
    if(!expect(
         sqmesh::geo::face_boundary_loops(
           model_handle,
           {sqmesh::geo::TopologyDimension::face, 0U},
           face_boundary,
           context
         ) == sqmesh::base::StatusCode::invalid_argument,
         "placeholder models should reject face boundary loop requests"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::geo::FeatureEdgeReport feature_report;
    if(!expect(
         sqmesh::geo::feature_edges(model_handle, feature_report, {}, context) ==
           sqmesh::base::StatusCode::invalid_argument,
         "placeholder models should reject feature edge extraction requests"
       )) {
      return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
  }

#if !defined(SQMESH_TEST_OCC_ENABLED)
  std::fprintf(
    stderr,
    "geo_topology_query_smoke: OCC runtime is enabled but the smoke target was not built with OCC.\n"
  );
  return EXIT_FAILURE;
#else
  sqmesh::geo::ModelHandle model_handle = sqmesh::invalid_handle;
  const TopoDS_Shape box = BRepPrimAPI_MakeBox(10.0, 20.0, 30.0).Shape();
  if(!expect(
       sqmesh::cad::occ::store_shape(box, model_handle, context) ==
         sqmesh::base::StatusCode::ok,
       "the test should be able to seed an OCC-backed solid model"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::TopologySnapshot snapshot;
  if(!expect(
       sqmesh::geo::topology_snapshot(model_handle, snapshot, context) ==
         sqmesh::base::StatusCode::ok,
       "topology_snapshot should inspect an OCC-backed model"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       snapshot.entity_count(sqmesh::geo::TopologyDimension::region) == 1U &&
         snapshot.entity_count(sqmesh::geo::TopologyDimension::face) == 6U &&
         snapshot.entity_count(sqmesh::geo::TopologyDimension::edge) == 12U &&
         snapshot.entity_count(sqmesh::geo::TopologyDimension::vertex) == 8U,
       "a solid OCC box should expose 1 region, 6 faces, 12 edges, and 8 vertices"
     )) {
    return EXIT_FAILURE;
  }

  for(const auto &face_info : snapshot.faces) {
    sqmesh::geo::FaceBoundaryLoops face_boundary;
    if(!expect(
         sqmesh::geo::face_boundary_loops(
           model_handle,
           face_info.entity,
           face_boundary,
           context
         ) == sqmesh::base::StatusCode::ok,
         "face_boundary_loops should succeed for each OCC box face"
       )) {
      return EXIT_FAILURE;
    }
    if(!expect(
         face_boundary.face == face_info.entity &&
           face_boundary.loops.size() == 1U &&
           face_boundary.primary_outer_loop_index == 0U &&
           face_boundary.outer_loop_count == 1U &&
           face_boundary.inner_loop_count == 0U &&
           face_boundary.unknown_loop_count == 0U &&
           face_boundary.closed_loop_count == 1U &&
           face_boundary.open_loop_count == 0U &&
           face_boundary.non_continuous_loop_count == 0U &&
           face_boundary.seam_loop_count == 0U &&
           face_boundary.seam_edge_use_count == 0U &&
           face_boundary.degenerate_loop_count == 0U &&
           face_boundary.degenerate_edge_use_count == 0U &&
           face_boundary.repeated_edge_loop_count == 0U &&
           face_boundary.repeated_edge_use_count == 0U &&
           !face_boundary.has_holes &&
           !face_boundary.has_seams &&
           face_boundary.loops[0].kind == sqmesh::geo::FaceBoundaryLoopKind::outer &&
           face_boundary.loops[0].closed &&
           face_boundary.loops[0].continuous &&
           face_boundary.loops[0].seam_edge_use_count == 0U &&
           face_boundary.loops[0].degenerate_edge_use_count == 0U &&
           face_boundary.loops[0].repeated_edge_use_count == 0U &&
           face_boundary.loops[0].vertex_ids.size() == 4U &&
           face_boundary.loops[0].edge_uses.size() == 4U,
         "each OCC box face should expose one classified closed outer loop with four ordered edges and four ordered corner vertices"
       )) {
      return EXIT_FAILURE;
    }

    std::vector<sqmesh::geo::TopologyEntityId> unordered_edges;
    if(!expect(
         sqmesh::geo::topology_children(
           model_handle,
           face_info.entity,
           unordered_edges,
           context
         ) == sqmesh::base::StatusCode::ok,
         "face child lookup should still succeed when validating ordered face loops"
       )) {
      return EXIT_FAILURE;
    }

    const auto &loop = face_boundary.loops[0];
    if(!expect(
         sqmesh::geo::primary_outer_boundary_loop(face_boundary) == &loop,
         "the richer boundary semantics should expose the primary outer loop directly"
       )) {
      return EXIT_FAILURE;
    }
    for(std::size_t edge_index = 0; edge_index < loop.edge_uses.size(); ++edge_index) {
      const auto &edge_use = loop.edge_uses[edge_index];
      if(!expect(
           sqmesh::geo::is_valid(edge_use.edge) &&
             sqmesh::geo::is_valid(edge_use.start_vertex) &&
             sqmesh::geo::is_valid(edge_use.end_vertex) &&
             !edge_use.is_seam && !edge_use.is_degenerate,
           "face boundary loops should return valid edge and vertex ids"
         )) {
        return EXIT_FAILURE;
      }

      bool found_in_face_children = false;
      for(const auto child_edge : unordered_edges) {
        if(child_edge == edge_use.edge) {
          found_in_face_children = true;
          break;
        }
      }
      if(!expect(
           found_in_face_children,
           "ordered face loops should reference only edges that belong to the face"
         )) {
        return EXIT_FAILURE;
      }

      const auto &next_edge_use = loop.edge_uses[(edge_index + 1U) % loop.edge_uses.size()];
      if(!expect(
           edge_use.end_vertex == next_edge_use.start_vertex,
           "ordered face loops should connect edge uses head-to-tail"
         )) {
        return EXIT_FAILURE;
      }
    }

    for(const auto child_edge : unordered_edges) {
      const auto use_count = static_cast<std::size_t>(std::count_if(
        loop.edge_uses.begin(),
        loop.edge_uses.end(),
        [&](const sqmesh::geo::FaceBoundaryEdgeUse &edge_use) {
          return edge_use.edge == child_edge;
        }
      ));
      if(!expect(
           use_count == 1U,
           "each face child edge should appear exactly once in the ordered face loop"
         )) {
        return EXIT_FAILURE;
      }
    }
  }

  std::vector<sqmesh::geo::TopologyEntityId> region_children;
  if(!expect(
       sqmesh::geo::topology_children(
         model_handle,
         snapshot.regions[0].entity,
         region_children,
         context
       ) == sqmesh::base::StatusCode::ok,
       "region child lookup should succeed"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       region_children.size() == 6U,
       "the single region should own six faces"
     )) {
    return EXIT_FAILURE;
  }

  std::vector<sqmesh::geo::TopologyEntityId> face_children;
  if(!expect(
       sqmesh::geo::topology_children(
         model_handle,
         snapshot.faces[0].entity,
         face_children,
         context
       ) == sqmesh::base::StatusCode::ok,
       "face child lookup should succeed"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       face_children.size() == 4U,
       "each box face should reference four edges"
     )) {
    return EXIT_FAILURE;
  }

  std::vector<sqmesh::geo::TopologyEntityId> face_parents;
  if(!expect(
       sqmesh::geo::topology_parents(
         model_handle,
         snapshot.faces[0].entity,
         face_parents,
         context
       ) == sqmesh::base::StatusCode::ok,
       "face parent lookup should succeed"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       face_parents.size() == 1U &&
         face_parents[0] == sqmesh::geo::TopologyEntityId {
           sqmesh::geo::TopologyDimension::region,
           0U,
         },
       "each box face should belong to the single solid region"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::FaceUvBounds face_bounds;
  if(!expect(
       sqmesh::geo::face_uv_bounds(
         model_handle,
         snapshot.faces[0].entity,
         face_bounds,
         context
       ) == sqmesh::base::StatusCode::ok,
       "face_uv_bounds should succeed"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       face_bounds.face == snapshot.faces[0].entity &&
         face_bounds.u_min < face_bounds.u_max &&
         face_bounds.v_min < face_bounds.v_max,
       "a box face should expose a valid trimmed UV box"
     )) {
    return EXIT_FAILURE;
  }

  const double u_mid = 0.5 * (face_bounds.u_min + face_bounds.u_max);
  const double v_mid = 0.5 * (face_bounds.v_min + face_bounds.v_max);
  sqmesh::geo::FaceSample face_sample;
  if(!expect(
       sqmesh::geo::sample_face(
         model_handle,
         snapshot.faces[0].entity,
         u_mid,
         v_mid,
         face_sample,
         context
       ) == sqmesh::base::StatusCode::ok,
       "sample_face should succeed inside the face UV box"
     )) {
    return EXIT_FAILURE;
  }
  const double normal_length = vector_norm(face_sample.normal);
  if(!expect(
       face_sample.face == snapshot.faces[0].entity &&
         face_sample.normal_defined &&
         normal_length > 0.99 && normal_length < 1.01,
       "sample_face should return a unit-length surface normal on the box face"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::FaceCurvatureSample face_curvature;
  if(!expect(
       sqmesh::geo::sample_face_curvature(
         model_handle,
         snapshot.faces[0].entity,
         u_mid,
         v_mid,
         face_curvature,
         context
       ) == sqmesh::base::StatusCode::ok,
       "sample_face_curvature should succeed inside the face UV box"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       face_curvature.face == snapshot.faces[0].entity &&
         face_curvature.curvature_defined &&
         std::abs(face_curvature.min_curvature) < 1.0e-12 &&
         std::abs(face_curvature.max_curvature) < 1.0e-12 &&
         std::abs(face_curvature.mean_curvature) < 1.0e-12 &&
         std::abs(face_curvature.gaussian_curvature) < 1.0e-12,
       "sample_face_curvature should return zero principal, mean, and Gaussian curvature on a box face"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::FaceDerivatives face_derivatives;
  if(!expect(
       sqmesh::geo::sample_face_derivatives(
         model_handle,
         snapshot.faces[0].entity,
         u_mid,
         v_mid,
         face_derivatives,
         context
       ) == sqmesh::base::StatusCode::ok,
       "sample_face_derivatives should succeed inside the face UV box"
     )) {
    return EXIT_FAILURE;
  }
  const double du_length = vector_norm(face_derivatives.du);
  const double dv_length = vector_norm(face_derivatives.dv);
  const sqmesh::geo::Vector3 derivative_normal {
    face_derivatives.du[1] * face_derivatives.dv[2] -
      face_derivatives.du[2] * face_derivatives.dv[1],
    face_derivatives.du[2] * face_derivatives.dv[0] -
      face_derivatives.du[0] * face_derivatives.dv[2],
    face_derivatives.du[0] * face_derivatives.dv[1] -
      face_derivatives.du[1] * face_derivatives.dv[0],
  };
  const double derivative_normal_length = vector_norm(derivative_normal);
  const double derivative_normal_alignment =
    derivative_normal_length > 0.0
      ? std::abs(dot_product(normalized(derivative_normal), face_sample.normal))
      : 0.0;
  if(!expect(
       face_derivatives.face == snapshot.faces[0].entity &&
         face_derivatives.first_derivatives_defined &&
         face_derivatives.second_derivatives_defined &&
         face_derivatives.normal_defined && du_length > 0.0 && dv_length > 0.0 &&
         derivative_normal_alignment > 0.99 &&
         std::abs(face_derivatives.duu[0]) < 1.0e-12 &&
         std::abs(face_derivatives.duu[1]) < 1.0e-12 &&
         std::abs(face_derivatives.duu[2]) < 1.0e-12 &&
         std::abs(face_derivatives.duv[0]) < 1.0e-12 &&
         std::abs(face_derivatives.duv[1]) < 1.0e-12 &&
         std::abs(face_derivatives.duv[2]) < 1.0e-12 &&
         std::abs(face_derivatives.dvv[0]) < 1.0e-12 &&
         std::abs(face_derivatives.dvv[1]) < 1.0e-12 &&
         std::abs(face_derivatives.dvv[2]) < 1.0e-12,
       "sample_face_derivatives should return a usable tangent basis and zero second derivatives on a box face"
     )) {
    return EXIT_FAILURE;
  }

  const auto projection_input = translated_point(face_sample.position, face_sample.normal, 7.5);
  sqmesh::geo::FaceProjection face_projection;
  if(!expect(
       sqmesh::geo::project_point_to_face(
         model_handle,
         snapshot.faces[0].entity,
         projection_input,
         face_projection,
         context
       ) == sqmesh::base::StatusCode::ok,
       "project_point_to_face should succeed for a point offset along the face normal"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       face_projection.face == snapshot.faces[0].entity &&
         face_projection.normal_defined &&
         point_distance(face_projection.projected_point, face_sample.position) < 1.0e-7 &&
         std::abs(face_projection.u - u_mid) < 1.0e-7 &&
         std::abs(face_projection.v - v_mid) < 1.0e-7 &&
         std::abs(face_projection.distance - 7.5) < 1.0e-7,
       "project_point_to_face should recover the sampled point, UV coordinates, and projection distance on a box face"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::FaceUvMapping face_uv_mapping;
  if(!expect(
       sqmesh::geo::recover_face_uv(
         model_handle,
         snapshot.faces[0].entity,
         face_sample.position,
         face_uv_mapping,
         context
       ) == sqmesh::base::StatusCode::ok,
       "recover_face_uv should succeed for a point already on the face"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       face_uv_mapping.face == snapshot.faces[0].entity &&
         point_distance(face_uv_mapping.mapped_point, face_sample.position) < 1.0e-7 &&
         std::abs(face_uv_mapping.u - u_mid) < 1.0e-7 &&
         std::abs(face_uv_mapping.v - v_mid) < 1.0e-7 &&
         face_uv_mapping.distance < 1.0e-7,
       "recover_face_uv should invert sample_face for a simple OCC box face"
     )) {
    return EXIT_FAILURE;
  }

  std::vector<sqmesh::geo::TopologyEntityId> edge_children;
  if(!expect(
       sqmesh::geo::topology_children(
         model_handle,
         snapshot.edges[0].entity,
         edge_children,
         context
       ) == sqmesh::base::StatusCode::ok,
       "edge child lookup should succeed"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       edge_children.size() == 2U,
       "each box edge should reference two vertices"
     )) {
    return EXIT_FAILURE;
  }

  std::vector<sqmesh::geo::TopologyEntityId> edge_parents;
  if(!expect(
       sqmesh::geo::topology_parents(
         model_handle,
         snapshot.edges[0].entity,
         edge_parents,
         context
       ) == sqmesh::base::StatusCode::ok,
       "edge parent lookup should succeed"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       edge_parents.size() == 2U,
       "each box edge should belong to two faces"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::EdgeCurveInfo edge_info;
  if(!expect(
       sqmesh::geo::edge_curve_info(
         model_handle,
         snapshot.edges[0].entity,
         edge_info,
         context
       ) == sqmesh::base::StatusCode::ok,
       "edge_curve_info should succeed"
     )) {
    return EXIT_FAILURE;
  }
  const double edge_dx = edge_info.end_point[0] - edge_info.start_point[0];
  const double edge_dy = edge_info.end_point[1] - edge_info.start_point[1];
  const double edge_dz = edge_info.end_point[2] - edge_info.start_point[2];
  const double endpoint_distance = std::sqrt(edge_dx * edge_dx + edge_dy * edge_dy + edge_dz * edge_dz);
  if(!expect(
       edge_info.edge == snapshot.edges[0].entity &&
         edge_info.parameter_min < edge_info.parameter_max &&
         edge_info.approximate_length > 0.0 &&
         endpoint_distance > 0.0,
       "edge_curve_info should return endpoints and a positive approximate length"
     )) {
    return EXIT_FAILURE;
  }

  const double parameter_mid = 0.5 * (edge_info.parameter_min + edge_info.parameter_max);
  sqmesh::geo::EdgeTangentSample edge_tangent;
  if(!expect(
       sqmesh::geo::sample_edge_tangent(
         model_handle,
         snapshot.edges[0].entity,
         parameter_mid,
         edge_tangent,
         context
       ) == sqmesh::base::StatusCode::ok,
       "sample_edge_tangent should succeed inside the edge parameter range"
     )) {
    return EXIT_FAILURE;
  }
  const sqmesh::geo::Vector3 edge_direction {
    edge_dx / endpoint_distance,
    edge_dy / endpoint_distance,
    edge_dz / endpoint_distance,
  };
  if(!expect(
       edge_tangent.edge == snapshot.edges[0].entity &&
         edge_tangent.tangent_defined && edge_tangent.speed > 0.0 &&
         dot_product(edge_tangent.tangent, edge_direction) > 0.99,
       "sample_edge_tangent should return a tangent aligned with the box edge direction"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::EdgeCurveSamples sampled_edge_curve;
  if(!expect(
       sqmesh::geo::sample_edge_curve(
         model_handle,
         snapshot.edges[0].entity,
         {
           0.0,
           4U,
         },
         sampled_edge_curve,
         context
       ) == sqmesh::base::StatusCode::ok,
       "sample_edge_curve should build reusable endpoint-inclusive edge samples"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sampled_edge_curve.curve.edge == snapshot.edges[0].entity &&
         sampled_edge_curve.samples.size() == 5U &&
         sampled_edge_curve.samples.front().edge == snapshot.edges[0].entity &&
         sampled_edge_curve.samples.back().edge == snapshot.edges[0].entity &&
         point_distance(
           sampled_edge_curve.samples.front().position,
           edge_info.start_point
         ) <= 1.0e-9 &&
         point_distance(
           sampled_edge_curve.samples.back().position,
           edge_info.end_point
         ) <= 1.0e-9 &&
         sampled_edge_curve.samples[2].tangent_defined &&
         dot_product(sampled_edge_curve.samples[2].tangent, edge_direction) > 0.99,
       "sample_edge_curve should expose ordered endpoint and midpoint tangent data that later meshers can reuse directly"
     )) {
    return EXIT_FAILURE;
  }

  const auto face_u_direction = normalized(face_derivatives.du);
  const double face_u_span = du_length * (face_bounds.u_max - face_bounds.u_min);
  const auto outside_projection_input = translated_point(
    translated_point(face_sample.position, face_u_direction, 2.0 * face_u_span),
    face_sample.normal,
    7.5
  );
  sqmesh::geo::FaceProjection outside_face_projection;
  if(!expect(
       sqmesh::geo::project_point_to_face(
         model_handle,
         snapshot.faces[0].entity,
         outside_projection_input,
         outside_face_projection,
         context
       ) == sqmesh::base::StatusCode::unsupported,
       "project_point_to_face should fail cleanly when the orthogonal projection leaves the trimmed face material"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::FeatureEdgeReport feature_report;
  sqmesh::geo::FeatureEdgeOptions feature_options;
  feature_options.feature_angle_degrees = 45.0;
  feature_options.include_boundary_edges = false;
  if(!expect(
       sqmesh::geo::feature_edges(
         model_handle,
         feature_report,
         feature_options,
         context
       ) == sqmesh::base::StatusCode::ok,
       "feature_edges should succeed for an OCC solid"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       feature_report.edges.size() == 12U &&
         feature_report.sharp_edge_count == 12U &&
         feature_report.boundary_edge_count == 0U &&
         feature_report.non_manifold_edge_count == 0U,
       "a box should classify all 12 solid edges as sharp feature edges"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::FeatureEdgeReport smooth_report;
  sqmesh::geo::FeatureEdgeOptions smooth_options;
  smooth_options.feature_angle_degrees = 120.0;
  smooth_options.include_boundary_edges = false;
  if(!expect(
       sqmesh::geo::feature_edges(
         model_handle,
         smooth_report,
         smooth_options,
         context
       ) == sqmesh::base::StatusCode::ok,
       "feature_edges should support higher threshold filters"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       smooth_report.edges.empty() &&
         smooth_report.sharp_edge_count == 0U,
       "a 120-degree threshold should exclude all box edges"
     )) {
    return EXIT_FAILURE;
  }

  std::vector<sqmesh::geo::TopologyEntityId> vertex_parents;
  if(!expect(
       sqmesh::geo::topology_parents(
         model_handle,
         snapshot.vertices[0].entity,
         vertex_parents,
         context
       ) == sqmesh::base::StatusCode::ok,
       "vertex parent lookup should succeed"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       vertex_parents.size() == 3U,
       "each box vertex should belong to three edges"
     )) {
    return EXIT_FAILURE;
  }

  std::vector<sqmesh::geo::TopologyEntityId> invalid_children;
  if(!expect(
       sqmesh::geo::topology_children(
         model_handle,
         {sqmesh::geo::TopologyDimension::face, 99U},
         invalid_children,
         context
       ) == sqmesh::base::StatusCode::invalid_argument,
       "topology traversal should reject out-of-range entity ids"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::FaceBoundaryLoops invalid_face_boundary;
  if(!expect(
       sqmesh::geo::face_boundary_loops(
         model_handle,
         {sqmesh::geo::TopologyDimension::face, 99U},
         invalid_face_boundary,
         context
       ) == sqmesh::base::StatusCode::invalid_argument,
       "face boundary loop lookup should reject out-of-range face entity ids"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::FaceDerivatives invalid_face_derivatives;
  if(!expect(
       sqmesh::geo::sample_face_derivatives(
         model_handle,
         {sqmesh::geo::TopologyDimension::face, 99U},
         u_mid,
         v_mid,
         invalid_face_derivatives,
         context
       ) == sqmesh::base::StatusCode::invalid_argument,
       "face derivative sampling should reject out-of-range face entity ids"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::FaceProjection invalid_face_projection;
  if(!expect(
       sqmesh::geo::project_point_to_face(
         model_handle,
         {sqmesh::geo::TopologyDimension::face, 99U},
         face_sample.position,
         invalid_face_projection,
         context
       ) == sqmesh::base::StatusCode::invalid_argument,
       "face projection should reject out-of-range face entity ids"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::FaceUvMapping invalid_face_uv_mapping;
  if(!expect(
       sqmesh::geo::recover_face_uv(
         model_handle,
         {sqmesh::geo::TopologyDimension::face, 99U},
         face_sample.position,
         invalid_face_uv_mapping,
         context
       ) == sqmesh::base::StatusCode::invalid_argument,
       "face inverse mapping should reject out-of-range face entity ids"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::EdgeTangentSample invalid_edge_tangent;
  if(!expect(
       sqmesh::geo::sample_edge_tangent(
         model_handle,
         {sqmesh::geo::TopologyDimension::edge, 99U},
         parameter_mid,
         invalid_edge_tangent,
         context
       ) == sqmesh::base::StatusCode::invalid_argument,
       "edge tangent sampling should reject out-of-range edge entity ids"
     )) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
#endif
}
