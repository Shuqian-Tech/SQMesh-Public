// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "sqmesh/sqmesh.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <vector>

namespace {

bool expect(bool condition, const char *message)
{
  if(condition) {
    return true;
  }

  std::fprintf(stderr, "geo_stl_smoke: %s\n", message);
  return false;
}

double vector_norm(const sqmesh::geo::Vector3 &vector)
{
  return std::sqrt(
    vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]
  );
}

bool write_patch_stl(const std::filesystem::path &path)
{
  std::ofstream stream(path);
  if(!stream) {
    return false;
  }

  stream <<
    "solid square_patch\n"
    "  facet normal 0 0 1\n"
    "    outer loop\n"
    "      vertex 0 0 0\n"
    "      vertex 1 0 0\n"
    "      vertex 1 1 0\n"
    "    endloop\n"
    "  endfacet\n"
    "  facet normal 0 0 1\n"
    "    outer loop\n"
    "      vertex 0 0 0\n"
    "      vertex 1 1 0\n"
    "      vertex 0 1 0\n"
    "    endloop\n"
    "  endfacet\n"
    "endsolid square_patch\n";

  return static_cast<bool>(stream);
}

bool contains_edge(
  const std::vector<sqmesh::geo::TopologyEntityId> &entities,
  sqmesh::geo::TopologyEntityId edge
)
{
  return std::find(entities.begin(), entities.end(), edge) != entities.end();
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

  const auto stl_path = std::filesystem::path("sqmesh_geo_stl_patch_ascii.stl");
  std::filesystem::remove(stl_path);
  if(!expect(
       write_patch_stl(stl_path),
       "the STL smoke test should be able to write its input file"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::ModelHandle model_handle = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::geo::import_stl(stl_path.string(), model_handle, {}, context) ==
         sqmesh::base::StatusCode::ok,
       "import_stl should create a usable geometry model handle"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }

  sqmesh::geo::ModelSummary summary;
  if(!expect(
       sqmesh::geo::model_summary(model_handle, summary, context) ==
         sqmesh::base::StatusCode::ok,
       "model_summary should succeed for STL-backed models"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }
  if(!expect(
       summary.face_count == 1U &&
         summary.edge_count == 5U &&
         summary.vertex_count == 4U &&
         summary.shell_count == 0U,
       "the square patch STL should expose one discrete face, five edges, and four vertices"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }

  sqmesh::geo::TopologySnapshot snapshot;
  if(!expect(
       sqmesh::geo::topology_snapshot(model_handle, snapshot, context) ==
         sqmesh::base::StatusCode::ok,
       "topology_snapshot should succeed for STL-backed models"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }
  if(!expect(
       snapshot.entity_count(sqmesh::geo::TopologyDimension::region) == 0U &&
         snapshot.entity_count(sqmesh::geo::TopologyDimension::face) == 1U &&
         snapshot.entity_count(sqmesh::geo::TopologyDimension::edge) == 5U &&
         snapshot.entity_count(sqmesh::geo::TopologyDimension::vertex) == 4U,
       "the STL patch should expose only discrete face, edge, and vertex topology"
     )) {
    std::filesystem::remove(stl_path);
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
       "topology_children should succeed for a discrete STL face"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }
  if(!expect(
       face_children.size() == 5U,
       "discrete STL face children should enumerate all incident triangle edges"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }

  sqmesh::geo::FaceBoundaryLoops boundary;
  if(!expect(
       sqmesh::geo::face_boundary_loops(
         model_handle,
         snapshot.faces[0].entity,
         boundary,
         context
       ) == sqmesh::base::StatusCode::ok,
       "face_boundary_loops should succeed for a discrete STL patch"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }
  if(!expect(
       boundary.face == snapshot.faces[0].entity &&
         boundary.loops.size() == 1U &&
         boundary.primary_outer_loop_index == sqmesh::geo::invalid_boundary_loop_index &&
         boundary.outer_loop_count == 0U &&
         boundary.inner_loop_count == 0U &&
         boundary.unknown_loop_count == 1U &&
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
         boundary.loops[0].kind == sqmesh::geo::FaceBoundaryLoopKind::unknown &&
         boundary.loops[0].closed &&
         boundary.loops[0].continuous &&
         boundary.loops[0].seam_edge_use_count == 0U &&
         boundary.loops[0].degenerate_edge_use_count == 0U &&
         boundary.loops[0].repeated_edge_use_count == 0U &&
         boundary.loops[0].vertex_ids.size() == 4U &&
         boundary.loops[0].edge_uses.size() == 4U,
       "the STL patch should expose one closed boundary loop with explicit unknown kind and direct ordered vertices"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::geo::primary_outer_boundary_loop(boundary) == nullptr,
       "discrete STL boundary semantics should stay OCC-free without inventing outer-loop trim classification"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }

  std::vector<sqmesh::geo::TopologyEntityId> boundary_edges;
  boundary_edges.reserve(boundary.loops[0].edge_uses.size());
  for(const auto &edge_use : boundary.loops[0].edge_uses) {
    if(!expect(
         !edge_use.is_seam && !edge_use.is_degenerate,
         "discrete STL boundary edge uses should keep seam and degenerate trim metadata explicitly false"
       )) {
      std::filesystem::remove(stl_path);
      return EXIT_FAILURE;
    }
    boundary_edges.push_back(edge_use.edge);
  }

  std::size_t face_interior_edge_count = 0U;
  for(const auto child_edge : face_children) {
    if(!contains_edge(boundary_edges, child_edge)) {
      ++face_interior_edge_count;
    }
  }
  if(!expect(
       face_interior_edge_count == 1U,
       "the discrete STL face should distinguish incident edges from boundary-loop edges"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }

  std::vector<sqmesh::geo::TopologyEntityId> edge_children;
  if(!expect(
       sqmesh::geo::topology_children(
         model_handle,
         boundary_edges[0],
         edge_children,
         context
       ) == sqmesh::base::StatusCode::ok,
       "topology_children should succeed for a discrete STL edge"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }
  if(!expect(
       edge_children.size() == 2U,
       "each discrete STL edge should expose two endpoint vertices"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }

  std::vector<sqmesh::geo::TopologyEntityId> edge_parents;
  if(!expect(
       sqmesh::geo::topology_parents(
         model_handle,
         boundary_edges[0],
         edge_parents,
         context
       ) == sqmesh::base::StatusCode::ok,
       "topology_parents should succeed for a discrete STL edge"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }
  if(!expect(
       edge_parents.size() == 1U && edge_parents[0] == snapshot.faces[0].entity,
       "each discrete STL edge in the patch should belong to the single discrete face"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }

  sqmesh::geo::EdgeCurveInfo edge_info;
  if(!expect(
       sqmesh::geo::edge_curve_info(model_handle, boundary_edges[0], edge_info, context) ==
         sqmesh::base::StatusCode::ok,
       "edge_curve_info should succeed for discrete STL edges"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }
  if(!expect(
       edge_info.edge == boundary_edges[0] &&
         edge_info.parameter_min == 0.0 &&
         edge_info.parameter_max == 1.0 &&
         edge_info.approximate_length > 0.0,
       "discrete STL edges should expose line-segment geometry information"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }

  sqmesh::geo::EdgeTangentSample edge_tangent;
  if(!expect(
       sqmesh::geo::sample_edge_tangent(
         model_handle,
         boundary_edges[0],
         0.25,
         edge_tangent,
         context
       ) == sqmesh::base::StatusCode::ok,
       "sample_edge_tangent should succeed for discrete STL edges"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }
  if(!expect(
       edge_tangent.edge == boundary_edges[0] &&
         edge_tangent.tangent_defined &&
         edge_tangent.speed > 0.0 &&
         vector_norm(edge_tangent.tangent) > 0.99 &&
         vector_norm(edge_tangent.tangent) < 1.01,
       "discrete STL edge tangent samples should expose a usable unit tangent"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }

  sqmesh::geo::EdgeCurveSamples sampled_edge_curve;
  if(!expect(
       sqmesh::geo::sample_edge_curve(
         model_handle,
         boundary_edges[0],
         {
           0.0,
           3U,
         },
         sampled_edge_curve,
         context
       ) == sqmesh::base::StatusCode::ok,
       "sample_edge_curve should succeed for discrete STL edges"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }
  if(!expect(
       sampled_edge_curve.curve.edge == boundary_edges[0] &&
         sampled_edge_curve.samples.size() == 4U &&
         sampled_edge_curve.samples.front().position == edge_info.start_point &&
         sampled_edge_curve.samples.back().position == edge_info.end_point &&
         sampled_edge_curve.samples[1].tangent_defined &&
         vector_norm(sampled_edge_curve.samples[1].tangent) > 0.99 &&
         vector_norm(sampled_edge_curve.samples[1].tangent) < 1.01,
       "discrete STL edges should expose reusable endpoint-inclusive edge samples through the shared geometry helper"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }

  sqmesh::geo::FeatureEdgeReport feature_report;
  sqmesh::geo::FeatureEdgeOptions feature_options;
  feature_options.feature_angle_degrees = 45.0;
  feature_options.include_boundary_edges = true;
  if(!expect(
       sqmesh::geo::feature_edges(
         model_handle,
         feature_report,
         feature_options,
         context
       ) == sqmesh::base::StatusCode::ok,
       "feature_edges should succeed for discrete STL models"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }
  if(!expect(
       feature_report.edges.size() == 4U &&
         feature_report.boundary_edge_count == 4U &&
         feature_report.sharp_edge_count == 0U &&
         feature_report.non_manifold_edge_count == 0U,
       "the square patch should classify only its four open-boundary edges as feature edges"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }

  sqmesh::geo::FaceUvBounds face_bounds;
  if(!expect(
       sqmesh::geo::face_uv_bounds(
         model_handle,
         snapshot.faces[0].entity,
         face_bounds,
         context
       ) == sqmesh::base::StatusCode::unsupported,
       "discrete STL models should reject UV-bound queries explicitly"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }

  sqmesh::geo::FaceSample face_sample;
  if(!expect(
       sqmesh::geo::sample_face(
         model_handle,
         snapshot.faces[0].entity,
         0.0,
         0.0,
         face_sample,
         context
       ) == sqmesh::base::StatusCode::unsupported,
       "discrete STL models should reject parametric face sampling explicitly"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }

  sqmesh::geo::FaceProjection face_projection;
  if(!expect(
       sqmesh::geo::project_point_to_face(
         model_handle,
         snapshot.faces[0].entity,
         {0.25, 0.25, 1.0},
         face_projection,
         context
       ) == sqmesh::base::StatusCode::unsupported,
       "discrete STL models should reject point-to-face projection explicitly"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }

  sqmesh::geo::FaceUvMapping face_uv_mapping;
  if(!expect(
       sqmesh::geo::recover_face_uv(
         model_handle,
         snapshot.faces[0].entity,
         {0.25, 0.25, 0.0},
         face_uv_mapping,
         context
       ) == sqmesh::base::StatusCode::unsupported,
       "discrete STL models should reject inverse-mapping queries explicitly"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }

  sqmesh::geo::TopologyCheckReport topology_report;
  if(!expect(
       sqmesh::geo::check_topology(model_handle, topology_report, context) ==
         sqmesh::base::StatusCode::unsupported,
       "discrete STL models should reject OCC-only topology inspection explicitly"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }

  std::filesystem::remove(stl_path);
  return EXIT_SUCCESS;
}
