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

#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <gp_Dir.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#endif

namespace {

bool expect(bool condition, const char *message)
{
  if(condition) {
    return true;
  }

  std::fprintf(stderr, "geo_model_view_smoke: %s\n", message);
  return false;
}

double vector_norm(const sqmesh::geo::Vector3 &vector)
{
  return std::sqrt(
    vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]
  );
}

template <typename T>
bool all_marked(const std::vector<T> &visited)
{
  return std::all_of(visited.begin(), visited.end(), [](T value) { return value; });
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

    sqmesh::geo::ModelView view;
    if(!expect(
         sqmesh::geo::model_view(model_handle, view, context) ==
           sqmesh::base::StatusCode::ok,
         "model_view should succeed for placeholder models"
       )) {
      return EXIT_FAILURE;
    }

    if(!expect(
         view.model_handle == model_handle &&
           view.entity_count(sqmesh::geo::TopologyDimension::region) == 0U &&
           view.entity_count(sqmesh::geo::TopologyDimension::face) == 0U &&
           view.entity_count(sqmesh::geo::TopologyDimension::edge) == 0U &&
           view.entity_count(sqmesh::geo::TopologyDimension::vertex) == 0U &&
           view.regions.empty() &&
           view.faces.empty() &&
           view.edges.empty() &&
           view.vertices.empty() &&
           view.root_face_ids.empty() &&
           view.root_edge_ids.empty() &&
           view.root_vertex_ids.empty() &&
           view.find_face({sqmesh::geo::TopologyDimension::face, 0U}) == nullptr,
         "placeholder model views should be empty and OCC-free"
       )) {
      return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
  }

#if !defined(SQMESH_TEST_OCC_ENABLED)
  std::fprintf(
    stderr,
    "geo_model_view_smoke: OCC runtime is enabled but the smoke target was not built with OCC.\n"
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

  sqmesh::geo::ModelView view;
  if(!expect(
       sqmesh::geo::model_view(model_handle, view, context) ==
         sqmesh::base::StatusCode::ok,
       "model_view should build a stable geometry traversal surface"
     )) {
    return EXIT_FAILURE;
  }

  if(!expect(
       view.entity_count(sqmesh::geo::TopologyDimension::region) == 1U &&
         view.entity_count(sqmesh::geo::TopologyDimension::face) == 6U &&
         view.entity_count(sqmesh::geo::TopologyDimension::edge) == 12U &&
         view.entity_count(sqmesh::geo::TopologyDimension::vertex) == 8U &&
         view.regions.size() == 1U &&
         view.faces.size() == 6U &&
         view.edges.size() == 12U &&
         view.vertices.size() == 8U &&
         view.root_face_ids.empty() &&
         view.root_edge_ids.empty() &&
         view.root_vertex_ids.empty(),
       "a solid OCC box should traverse through region, face, edge, and vertex views"
     )) {
    return EXIT_FAILURE;
  }

  std::vector<bool> visited_faces(view.faces.size(), false);
  std::vector<bool> visited_edges(view.edges.size(), false);
  std::vector<bool> visited_vertices(view.vertices.size(), false);

  for(const auto &region : view.regions) {
    if(!expect(
         region.face_ids.size() == 6U,
         "the box region should expose six child faces from the model view"
       )) {
      return EXIT_FAILURE;
    }

    for(const auto face_id : region.face_ids) {
      const auto *face = view.find_face(face_id);
      if(!expect(face != nullptr, "region face ids should resolve to face views")) {
        return EXIT_FAILURE;
      }
      visited_faces[face->entity.index] = true;

      if(!expect(
           face->region_ids.size() == 1U &&
             face->region_ids[0] == region.entity &&
             face->edge_ids.size() == 4U,
           "face views should expose parent region ids and child edge ids"
         )) {
        return EXIT_FAILURE;
      }

      const auto &boundary = sqmesh::geo::ordered_boundary_loops(*face);
      const auto *outer_loop = sqmesh::geo::primary_outer_boundary_loop(boundary);
      if(!expect(
           boundary.face == face->entity &&
             boundary.loops.size() == 1U &&
             boundary.primary_outer_loop_index == 0U &&
             boundary.outer_loop_count == 1U &&
             boundary.inner_loop_count == 0U &&
             boundary.unknown_loop_count == 0U &&
             boundary.closed_loop_count == 1U &&
             boundary.open_loop_count == 0U &&
             boundary.non_continuous_loop_count == 0U &&
             boundary.seam_loop_count == 0U &&
             boundary.seam_edge_use_count == 0U &&
             boundary.degenerate_loop_count == 0U &&
             boundary.degenerate_edge_use_count == 0U &&
             boundary.repeated_edge_loop_count == 0U &&
             boundary.repeated_edge_use_count == 0U &&
             !boundary.has_holes &&
             !boundary.has_seams &&
             outer_loop != nullptr &&
             boundary.loops[0].kind == sqmesh::geo::FaceBoundaryLoopKind::outer &&
             boundary.loops[0].closed &&
             boundary.loops[0].continuous &&
             boundary.loops[0].seam_edge_use_count == 0U &&
             boundary.loops[0].degenerate_edge_use_count == 0U &&
             boundary.loops[0].repeated_edge_use_count == 0U &&
             boundary.loops[0].vertex_ids.size() == 4U &&
             boundary.loops[0].edge_uses.size() == 4U,
           "face views should carry richer ordered boundary semantics directly"
         )) {
        return EXIT_FAILURE;
      }

      sqmesh::geo::FaceBoundaryLoops boundary_copy;
      if(!expect(
           sqmesh::geo::face_boundary_loops(*face, boundary_copy) ==
             sqmesh::base::StatusCode::ok &&
             boundary_copy.face == face->entity &&
             boundary_copy.loops.size() == boundary.loops.size() &&
             boundary_copy.outer_loop_count == boundary.outer_loop_count &&
             boundary_copy.closed_loop_count == boundary.closed_loop_count &&
             sqmesh::geo::primary_outer_boundary_loop(boundary_copy) != nullptr,
           "face boundary access should continue directly from the face view"
         )) {
        return EXIT_FAILURE;
      }

      sqmesh::geo::FaceUvBounds bounds;
      if(!expect(
           sqmesh::geo::face_uv_bounds(*face, bounds) ==
             sqmesh::base::StatusCode::ok &&
             bounds.face == face->entity &&
             bounds.u_min < bounds.u_max &&
             bounds.v_min < bounds.v_max,
           "face queries should accept the face view directly"
         )) {
        return EXIT_FAILURE;
      }

      const double u_mid = 0.5 * (bounds.u_min + bounds.u_max);
      const double v_mid = 0.5 * (bounds.v_min + bounds.v_max);
      sqmesh::geo::FaceSample sample;
      if(!expect(
           sqmesh::geo::sample_face(*face, u_mid, v_mid, sample) ==
             sqmesh::base::StatusCode::ok &&
             sample.face == face->entity &&
             sample.normal_defined &&
             vector_norm(sample.normal) > 0.99 &&
             vector_norm(sample.normal) < 1.01,
           "sample_face should work from the face view without replaying ids"
         )) {
        return EXIT_FAILURE;
      }

      sqmesh::geo::FaceProjection projection;
      if(!expect(
           sqmesh::geo::project_point_to_face(*face, sample.position, projection) ==
             sqmesh::base::StatusCode::ok &&
             projection.face == face->entity &&
             projection.distance <= 1.0e-9,
           "project_point_to_face should continue directly from the face view"
         )) {
        return EXIT_FAILURE;
      }

      for(const auto &edge_use : boundary.loops.front().edge_uses) {
        const auto *edge = view.find_edge(edge_use.edge);
        if(!expect(edge != nullptr, "boundary edge ids should resolve to edge views")) {
          return EXIT_FAILURE;
        }
        visited_edges[edge->entity.index] = true;

        if(!expect(
             edge->vertex_ids.size() == 2U &&
               std::find(face->edge_ids.begin(), face->edge_ids.end(), edge->entity) !=
                 face->edge_ids.end() &&
               std::find(edge->vertex_ids.begin(), edge->vertex_ids.end(), edge_use.start_vertex) !=
                 edge->vertex_ids.end() &&
               std::find(edge->vertex_ids.begin(), edge->vertex_ids.end(), edge_use.end_vertex) !=
                 edge->vertex_ids.end(),
             "edge views should expose child vertices and match ordered face boundary uses"
           )) {
          return EXIT_FAILURE;
        }

        sqmesh::geo::EdgeCurveInfo edge_info;
        if(!expect(
             sqmesh::geo::edge_curve_info(*edge, edge_info) ==
               sqmesh::base::StatusCode::ok &&
               edge_info.edge == edge->entity &&
               edge_info.parameter_min < edge_info.parameter_max &&
               edge_info.approximate_length > 0.0,
             "edge queries should accept the edge view directly"
           )) {
          return EXIT_FAILURE;
        }

        sqmesh::geo::EdgeTangentSample tangent;
        if(!expect(
             sqmesh::geo::sample_edge_tangent(*edge, edge_info.parameter_min, tangent) ==
               sqmesh::base::StatusCode::ok &&
               tangent.edge == edge->entity &&
               tangent.speed > 0.0,
             "sample_edge_tangent should continue directly from the edge view"
           )) {
          return EXIT_FAILURE;
        }

        sqmesh::geo::EdgeCurveSamples sampled_curve;
        if(!expect(
             sqmesh::geo::sample_edge_curve(
               *edge,
               {
                 0.0,
                 2U,
               },
               sampled_curve
             ) == sqmesh::base::StatusCode::ok &&
               sampled_curve.curve.edge == edge->entity &&
               sampled_curve.samples.size() == 3U &&
               sampled_curve.samples.front().edge == edge->entity &&
               sampled_curve.samples.back().edge == edge->entity,
             "sample_edge_curve should continue directly from the edge view"
           )) {
          return EXIT_FAILURE;
        }

        for(const auto vertex_id : edge->vertex_ids) {
          const auto *vertex = view.find_vertex(vertex_id);
          if(!expect(vertex != nullptr, "edge vertex ids should resolve to vertex views")) {
            return EXIT_FAILURE;
          }
          visited_vertices[vertex->entity.index] = true;

          if(!expect(
               !vertex->edge_ids.empty(),
               "vertex views should carry their incident edge ids"
             )) {
            return EXIT_FAILURE;
          }
        }
      }
    }
  }

  if(!expect(
       all_marked(visited_faces) &&
         all_marked(visited_edges) &&
         all_marked(visited_vertices),
       "top-down region-face-edge-vertex traversal should cover the full OCC box topology"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::ModelHandle holed_face_handle = sqmesh::invalid_handle;
  const TopoDS_Face holed_face = make_rectangular_face_with_hole();
  if(!expect(
       sqmesh::cad::occ::store_shape(holed_face, holed_face_handle, context) ==
         sqmesh::base::StatusCode::ok,
       "the test should be able to seed an OCC-backed face with an inner trim loop"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::ModelView holed_face_view;
  if(!expect(
       sqmesh::geo::model_view(holed_face_handle, holed_face_view, context) ==
         sqmesh::base::StatusCode::ok,
       "model_view should build trim-aware semantics for an OCC face with a hole"
     )) {
    return EXIT_FAILURE;
  }

  if(!expect(
       holed_face_view.regions.empty() &&
         holed_face_view.faces.size() == 1U &&
         holed_face_view.edges.size() == 8U &&
         holed_face_view.vertices.size() == 8U &&
         holed_face_view.root_face_ids.size() == 1U &&
         holed_face_view.root_face_ids[0] == holed_face_view.faces[0].entity,
       "a standalone trimmed face should remain reachable through the stable OCC-free model view"
     )) {
    return EXIT_FAILURE;
  }

  const auto &trimmed_face = holed_face_view.faces[0];
  const auto &trimmed_boundary = sqmesh::geo::ordered_boundary_loops(trimmed_face);
  const auto *trimmed_outer_loop =
    sqmesh::geo::primary_outer_boundary_loop(trimmed_boundary);
  const auto inner_loop_it = std::find_if(
    trimmed_boundary.loops.begin(),
    trimmed_boundary.loops.end(),
    [](const sqmesh::geo::FaceBoundaryLoop &loop) {
      return loop.kind == sqmesh::geo::FaceBoundaryLoopKind::inner;
    }
  );
  if(!expect(
       trimmed_face.edge_ids.size() == 8U &&
         trimmed_boundary.face == trimmed_face.entity &&
         trimmed_boundary.loops.size() == 2U &&
         trimmed_boundary.primary_outer_loop_index != sqmesh::geo::invalid_boundary_loop_index &&
         trimmed_boundary.outer_loop_count == 1U &&
         trimmed_boundary.inner_loop_count == 1U &&
         trimmed_boundary.unknown_loop_count == 0U &&
         trimmed_boundary.closed_loop_count == 2U &&
         trimmed_boundary.open_loop_count == 0U &&
         trimmed_boundary.non_continuous_loop_count == 0U &&
         trimmed_boundary.seam_loop_count == 0U &&
         trimmed_boundary.seam_edge_use_count == 0U &&
         trimmed_boundary.degenerate_loop_count == 0U &&
         trimmed_boundary.degenerate_edge_use_count == 0U &&
         trimmed_boundary.repeated_edge_loop_count == 0U &&
         trimmed_boundary.repeated_edge_use_count == 0U &&
         trimmed_boundary.has_holes &&
         !trimmed_boundary.has_seams &&
         trimmed_outer_loop != nullptr &&
         trimmed_outer_loop->closed &&
         trimmed_outer_loop->continuous &&
         trimmed_outer_loop->seam_edge_use_count == 0U &&
         trimmed_outer_loop->degenerate_edge_use_count == 0U &&
         trimmed_outer_loop->repeated_edge_use_count == 0U &&
         trimmed_outer_loop->vertex_ids.size() == 4U &&
         inner_loop_it != trimmed_boundary.loops.end() &&
         inner_loop_it->closed &&
         inner_loop_it->continuous &&
         inner_loop_it->seam_edge_use_count == 0U &&
         inner_loop_it->degenerate_edge_use_count == 0U &&
         inner_loop_it->repeated_edge_use_count == 0U &&
         inner_loop_it->vertex_ids.size() == 4U,
       "trimmed OCC faces should expose explicit outer/inner loop semantics and ordered loop vertices through the stable view API"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::FaceBoundaryLoops trimmed_boundary_copy;
  if(!expect(
       sqmesh::geo::face_boundary_loops(trimmed_face, trimmed_boundary_copy) ==
         sqmesh::base::StatusCode::ok &&
         trimmed_boundary_copy.face == trimmed_face.entity &&
         trimmed_boundary_copy.outer_loop_count == 1U &&
         trimmed_boundary_copy.inner_loop_count == 1U &&
         trimmed_boundary_copy.has_holes &&
         sqmesh::geo::primary_outer_boundary_loop(trimmed_boundary_copy) != nullptr,
       "trim-aware face semantics should remain consumable from the face view without replaying backend traversal"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::ModelHandle cylinder_handle = sqmesh::invalid_handle;
  const TopoDS_Shape cylinder = BRepPrimAPI_MakeCylinder(2.0, 5.0).Shape();
  if(!expect(
       sqmesh::cad::occ::store_shape(cylinder, cylinder_handle, context) ==
         sqmesh::base::StatusCode::ok,
       "the test should be able to seed an OCC-backed cylinder with an implicit seam face"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::ModelView cylinder_view;
  if(!expect(
       sqmesh::geo::model_view(cylinder_handle, cylinder_view, context) ==
         sqmesh::base::StatusCode::ok,
       "model_view should preserve seam-aware ordered boundary semantics"
     )) {
    return EXIT_FAILURE;
  }

  const auto seam_face_count = static_cast<std::size_t>(std::count_if(
    cylinder_view.faces.begin(),
    cylinder_view.faces.end(),
    [](const sqmesh::geo::FaceView &face) {
      return sqmesh::geo::ordered_boundary_loops(face).has_seams;
    }
  ));
  const auto seam_face_it = std::find_if(
    cylinder_view.faces.begin(),
    cylinder_view.faces.end(),
    [](const sqmesh::geo::FaceView &face) {
      return sqmesh::geo::ordered_boundary_loops(face).has_seams;
    }
  );
  if(!expect(
       seam_face_count == 1U && seam_face_it != cylinder_view.faces.end(),
       "the cylinder should expose exactly one seam-bearing face through the stable geometry view"
     )) {
    return EXIT_FAILURE;
  }

  const auto &seam_boundary = sqmesh::geo::ordered_boundary_loops(*seam_face_it);
  if(!expect(
       !seam_boundary.loops.empty(),
       "the seam-bearing cylinder face should still expose at least one ordered boundary loop"
     )) {
    return EXIT_FAILURE;
  }
  const auto &seam_loop = seam_boundary.loops[0];
  const auto seam_use_count = static_cast<std::size_t>(std::count_if(
    seam_loop.edge_uses.begin(),
    seam_loop.edge_uses.end(),
    [](const sqmesh::geo::FaceBoundaryEdgeUse &edge_use) {
      return edge_use.is_seam;
    }
  ));
  if(!expect(
       seam_boundary.loops.size() == 1U &&
         seam_boundary.outer_loop_count == 1U &&
         seam_boundary.inner_loop_count == 0U &&
         seam_boundary.unknown_loop_count == 0U &&
         seam_boundary.closed_loop_count == 1U &&
         seam_boundary.open_loop_count == 0U &&
         seam_boundary.non_continuous_loop_count == 0U &&
         seam_boundary.seam_loop_count == 1U &&
         seam_boundary.seam_edge_use_count == 2U &&
         seam_boundary.degenerate_loop_count == 0U &&
         seam_boundary.degenerate_edge_use_count == 0U &&
         seam_boundary.repeated_edge_loop_count == 1U &&
         seam_boundary.repeated_edge_use_count == 1U &&
         !seam_boundary.has_holes &&
         seam_boundary.has_seams &&
         seam_face_it->edge_ids.size() == 3U &&
         seam_loop.closed &&
         seam_loop.continuous &&
         seam_loop.seam_edge_use_count == 2U &&
         seam_loop.degenerate_edge_use_count == 0U &&
         seam_loop.repeated_edge_use_count == 1U &&
         seam_loop.edge_uses.size() == 4U &&
         seam_loop.vertex_ids.size() == 4U &&
         seam_use_count == seam_loop.seam_edge_use_count,
       "seam-bearing OCC faces should expose explicit seam and repeated-edge diagnostics instead of leaving them implicit in ordered loops"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::FaceBoundaryLoops seam_boundary_copy;
  if(!expect(
       sqmesh::geo::face_boundary_loops(*seam_face_it, seam_boundary_copy) ==
         sqmesh::base::StatusCode::ok &&
         seam_boundary_copy.has_seams &&
         seam_boundary_copy.seam_loop_count == 1U &&
         seam_boundary_copy.seam_edge_use_count == 2U &&
         seam_boundary_copy.repeated_edge_loop_count == 1U &&
         seam_boundary_copy.repeated_edge_use_count == 1U,
       "seam-aware boundary diagnostics should remain consumable directly from the face view"
     )) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
#endif
}
