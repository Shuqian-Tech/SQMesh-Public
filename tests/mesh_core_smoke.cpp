// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "sqmesh/sqmesh.hpp"

#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <type_traits>
#include <utility>

namespace {

bool expect(bool condition, const char *message)
{
  if(condition) {
    return true;
  }

  std::fprintf(stderr, "mesh_core_smoke: %s\n", message);
  return false;
}

} // namespace

constexpr auto kMeshCoreLayout = sqmesh::mesh::mesh_core_layout();
constexpr auto kExpectedMeshCoreLayout = sqmesh::mesh::expected_mesh_core_layout();

static_assert(sqmesh::mesh::matches_mesh_core_layout_baseline(kMeshCoreLayout),
              "Mesh-core layout baseline drifted.");
static_assert(kMeshCoreLayout.entity_ref_size == kExpectedMeshCoreLayout.entity_ref_size,
              "EntityRef size drifted.");
static_assert(kMeshCoreLayout.entity_ref_alignment == kExpectedMeshCoreLayout.entity_ref_alignment,
              "EntityRef alignment drifted.");
static_assert(kMeshCoreLayout.connectivity_span_size == kExpectedMeshCoreLayout.connectivity_span_size,
              "ConnectivitySpan size drifted.");
static_assert(
  kMeshCoreLayout.connectivity_span_alignment ==
    kExpectedMeshCoreLayout.connectivity_span_alignment,
  "ConnectivitySpan alignment drifted."
);
static_assert(kMeshCoreLayout.entity_header_size == kExpectedMeshCoreLayout.entity_header_size,
              "EntityHeader size drifted.");
static_assert(kMeshCoreLayout.entity_header_alignment ==
                kExpectedMeshCoreLayout.entity_header_alignment,
              "EntityHeader alignment drifted.");
static_assert(kMeshCoreLayout.node_size == kExpectedMeshCoreLayout.node_size,
              "Node size drifted.");
static_assert(kMeshCoreLayout.node_alignment == kExpectedMeshCoreLayout.node_alignment,
              "Node alignment drifted.");
static_assert(kMeshCoreLayout.edge_size == kExpectedMeshCoreLayout.edge_size,
              "Edge size drifted.");
static_assert(kMeshCoreLayout.edge_alignment == kExpectedMeshCoreLayout.edge_alignment,
              "Edge alignment drifted.");
static_assert(kMeshCoreLayout.face_size == kExpectedMeshCoreLayout.face_size,
              "Face size drifted.");
static_assert(kMeshCoreLayout.face_alignment == kExpectedMeshCoreLayout.face_alignment,
              "Face alignment drifted.");
static_assert(kMeshCoreLayout.cell_size == kExpectedMeshCoreLayout.cell_size,
              "Cell size drifted.");
static_assert(kMeshCoreLayout.cell_alignment == kExpectedMeshCoreLayout.cell_alignment,
              "Cell alignment drifted.");
static_assert(kMeshCoreLayout.node_header_offset == kExpectedMeshCoreLayout.node_header_offset,
              "Node header offset drifted.");
static_assert(kMeshCoreLayout.edge_header_offset == kExpectedMeshCoreLayout.edge_header_offset,
              "Edge header offset drifted.");
static_assert(kMeshCoreLayout.edge_node_span_offset ==
                kExpectedMeshCoreLayout.edge_node_span_offset,
              "Edge node-span offset drifted.");
static_assert(kMeshCoreLayout.edge_left_face_offset ==
                kExpectedMeshCoreLayout.edge_left_face_offset,
              "Edge left-face offset drifted.");
static_assert(kMeshCoreLayout.edge_right_face_offset ==
                kExpectedMeshCoreLayout.edge_right_face_offset,
              "Edge right-face offset drifted.");
static_assert(kMeshCoreLayout.face_header_offset == kExpectedMeshCoreLayout.face_header_offset,
              "Face header offset drifted.");
static_assert(kMeshCoreLayout.face_node_span_offset ==
                kExpectedMeshCoreLayout.face_node_span_offset,
              "Face node-span offset drifted.");
static_assert(kMeshCoreLayout.face_left_cell_offset ==
                kExpectedMeshCoreLayout.face_left_cell_offset,
              "Face left-cell offset drifted.");
static_assert(kMeshCoreLayout.face_right_cell_offset ==
                kExpectedMeshCoreLayout.face_right_cell_offset,
              "Face right-cell offset drifted.");
static_assert(kMeshCoreLayout.cell_header_offset == kExpectedMeshCoreLayout.cell_header_offset,
              "Cell header offset drifted.");
static_assert(kMeshCoreLayout.cell_node_span_offset ==
                kExpectedMeshCoreLayout.cell_node_span_offset,
              "Cell node-span offset drifted.");
static_assert(kMeshCoreLayout.cell_face_span_offset ==
                kExpectedMeshCoreLayout.cell_face_span_offset,
              "Cell face-span offset drifted.");
static_assert(sqmesh::mesh::expected_layout_benchmark_bytes_per_cell() == 584U,
              "Layout benchmark bytes-per-cell baseline drifted.");
static_assert(std::is_standard_layout_v<sqmesh::mesh::Node>, "Node must stay standard-layout.");
static_assert(std::is_standard_layout_v<sqmesh::mesh::Edge>, "Edge must stay standard-layout.");
static_assert(std::is_standard_layout_v<sqmesh::mesh::Face>, "Face must stay standard-layout.");
static_assert(std::is_standard_layout_v<sqmesh::mesh::Cell>, "Cell must stay standard-layout.");

int main()
{
  using sqmesh::mesh::Domain;
  using sqmesh::mesh::EntityKind;
  using sqmesh::mesh::EntityOrder;
  using sqmesh::mesh::FaceSide;
  using sqmesh::mesh::NastranEntityGroupSourceCard;
  using sqmesh::mesh::EntityGroupSemantic;
  using sqmesh::mesh::EntityGroupRole;
  using sqmesh::mesh::EntityGroupImportFormat;
  using sqmesh::mesh::invalid_index;
  using sqmesh::mesh::is_interface_entity;
  using sqmesh::mesh::is_valid;

  Domain domain("mesh_core_smoke");

  const auto node_entity_group = domain.create_entity_group({EntityOrder::node, "nodes"});
  sqmesh::mesh::EntityGroupDefinition face_definition;
  face_definition.order = EntityOrder::face;
  face_definition.name = "boundary_faces";
  face_definition.boundary = true;
  face_definition.default_kind = EntityKind::face_triangle;
  face_definition.import_info.format = EntityGroupImportFormat::nastran;
  face_definition.import_info.nastran.source_card = NastranEntityGroupSourceCard::ctria3;
  const auto face_entity_group = domain.create_entity_group(std::move(face_definition));
  const auto cell_entity_group = domain.create_entity_group(
    {EntityOrder::cell, "cells", invalid_index, false, EntityKind::cell_tetra}
  );

  domain.reserve_entity_group_storage(node_entity_group, 4U);
  domain.reserve_entity_group_storage(face_entity_group, 4U, 12U);
  domain.reserve_entity_group_storage(cell_entity_group, 1U, 4U, 4U);

  const auto n0 = domain.add_node(node_entity_group, {0.0, 0.0, 0.0});
  const auto n1 = domain.add_node(node_entity_group, {1.0, 0.0, 0.0});
  const auto n2 = domain.add_node(node_entity_group, {0.0, 1.0, 0.0});
  const auto n3 = domain.add_node(node_entity_group, {0.0, 0.0, 1.0});

  const auto f0 = domain.add_triangle_face(face_entity_group, {n0, n2, n1});
  const auto f1 = domain.add_triangle_face(face_entity_group, {n0, n1, n3});
  const auto f2 = domain.add_triangle_face(face_entity_group, {n1, n2, n3});
  const auto f3 = domain.add_triangle_face(face_entity_group, {n2, n0, n3});

  const auto c0 = domain.add_tetra_cell(cell_entity_group, {n0, n1, n2, n3}, {f0, f1, f2, f3});

  domain.set_face_cells(f0, c0);
  domain.set_face_cells(f1, c0);
  domain.set_face_cells(f2, c0);
  domain.set_face_cells(f3, c0);
  domain.set_face_topology_owner(
    f0,
    {
      sqmesh::geo::TopologyDimension::face,
      7U,
    }
  );
  domain.set_face_source_entity_tag(f0, 101U);
  domain.set_source_topology_revision(33U);

  const auto summary = domain.summary();
  const auto stats = domain.statistics();

  if(!expect(summary.node_count == 4U, "summary should report four nodes")) {
    return EXIT_FAILURE;
  }
  if(!expect(summary.edge_count == 0U, "summary should report zero edges when none were inserted")) {
    return EXIT_FAILURE;
  }
  if(!expect(summary.face_count == 4U, "summary should report four faces")) {
    return EXIT_FAILURE;
  }
  if(!expect(summary.cell_count == 1U, "summary should report one cell")) {
    return EXIT_FAILURE;
  }
  if(!expect(summary.source_topology_revision == 33U,
             "summary should preserve the source topology revision")) {
    return EXIT_FAILURE;
  }
  if(!expect(stats.entity_group_count == 3U, "statistics should report three entity_groups")) {
    return EXIT_FAILURE;
  }
  if(!expect(stats.boundary_edge_count == 0U, "statistics should report zero boundary edges")) {
    return EXIT_FAILURE;
  }
  if(!expect(stats.interior_edge_count == 0U, "statistics should report zero interior edges")) {
    return EXIT_FAILURE;
  }
  if(!expect(stats.boundary_face_count == 4U, "single tetra should expose four boundary faces")) {
    return EXIT_FAILURE;
  }
  if(!expect(stats.interior_face_count == 0U, "single tetra should expose zero interior faces")) {
    return EXIT_FAILURE;
  }
  if(!expect(stats.line_edge_count == 0U, "line edge count should stay zero without edge entities")) {
    return EXIT_FAILURE;
  }
  if(!expect(stats.triangle_face_count == 4U, "triangle face count should match inserted faces")) {
    return EXIT_FAILURE;
  }
  if(!expect(stats.tetra_cell_count == 1U, "tetra cell count should match inserted cells")) {
    return EXIT_FAILURE;
  }
  if(!expect(domain.entity_group(node_entity_group).semantic() == EntityGroupSemantic::node,
             "node entity_groups should expose explicit node semantics")) {
    return EXIT_FAILURE;
  }
  if(!expect(domain.entity_group(face_entity_group).semantic() == EntityGroupSemantic::boundary,
             "boundary face entity_groups should expose explicit boundary semantics")) {
    return EXIT_FAILURE;
  }
  if(!expect(
       domain.entity_group(face_entity_group).import_info().format == EntityGroupImportFormat::nastran &&
         domain.entity_group(face_entity_group).import_info().nastran.source_card ==
           NastranEntityGroupSourceCard::ctria3,
       "entity_group import metadata should preserve narrow entity_group-level import provenance"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(domain.entity_group(cell_entity_group).semantic() == EntityGroupSemantic::region,
             "cell entity_groups should expose explicit region semantics")) {
    return EXIT_FAILURE;
  }
  if(!expect(domain.entity_group(cell_entity_group).primary_region_zone_id() == cell_entity_group,
             "default region semantics should bind the cell entity_group zone id as the region id")) {
    return EXIT_FAILURE;
  }

  const auto &node_storage = domain.entity_group(node_entity_group).nodes();
  const auto &face_storage = domain.entity_group(face_entity_group).faces();
  const auto &cell_storage = domain.entity_group(cell_entity_group).cells();

  if(!expect(node_storage.size() == 4U, "node entity_group should own four nodes")) {
    return EXIT_FAILURE;
  }
  if(!expect(face_storage.size() == 4U, "face entity_group should own four faces")) {
    return EXIT_FAILURE;
  }
  if(!expect(cell_storage.size() == 1U, "cell entity_group should own one cell")) {
    return EXIT_FAILURE;
  }
  if(!expect(&node_storage[0] + 1 == &node_storage[1], "node pool should stay contiguous")) {
    return EXIT_FAILURE;
  }
  if(!expect(&face_storage[0] + 1 == &face_storage[1], "face pool should stay contiguous")) {
    return EXIT_FAILURE;
  }

  const auto face_nodes = domain.face_nodes(f0);
  const auto cell_nodes = domain.cell_nodes(c0);
  const auto cell_faces = domain.cell_faces(c0);

  if(!expect(face_nodes.size == 3U, "tri face should reference three nodes")) {
    return EXIT_FAILURE;
  }
  if(!expect(cell_nodes.size == 4U, "tet cell should reference four nodes")) {
    return EXIT_FAILURE;
  }
  if(!expect(cell_faces.size == 4U, "tet cell should reference four faces")) {
    return EXIT_FAILURE;
  }
  if(!expect(face_nodes[0] == n0 && face_nodes[1] == n2 && face_nodes[2] == n1,
             "face connectivity should preserve insertion order")) {
    return EXIT_FAILURE;
  }
  if(!expect(cell_nodes[0] == n0 && cell_nodes[3] == n3,
             "cell node channel should preserve references")) {
    return EXIT_FAILURE;
  }
  if(!expect(cell_faces[0] == f0 && cell_faces[3] == f3,
             "cell face channel should preserve references")) {
    return EXIT_FAILURE;
  }
  if(!expect(&face_nodes[0] + 1 == &face_nodes[1], "face node side-array should stay contiguous")) {
    return EXIT_FAILURE;
  }
  if(!expect(&cell_faces[0] + 1 == &cell_faces[1], "cell face side-array should stay contiguous")) {
    return EXIT_FAILURE;
  }

  if(!expect(domain.adjacent_cell(f0, FaceSide::left) == c0,
             "left face adjacency should resolve in O(1)")) {
    return EXIT_FAILURE;
  }
  if(!expect(
       domain.face_topology_owner(f0) ==
         sqmesh::geo::TopologyEntityId {sqmesh::geo::TopologyDimension::face, 7U},
       "face topology ownership should preserve per-face source metadata"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       !sqmesh::geo::is_valid(domain.face_topology_owner(f1)),
       "faces without an assigned topology owner should remain unset"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(domain.face_source_entity_tag(f0) == 101U,
             "face source-entity tags should preserve per-face import provenance")) {
    return EXIT_FAILURE;
  }
  if(!expect(domain.face_source_entity_tag(f1) == sqmesh::mesh::invalid_index,
             "faces without an assigned source-entity tag should remain unset")) {
    return EXIT_FAILURE;
  }
  if(!expect(!is_valid(domain.adjacent_cell(f0, FaceSide::right)),
             "single tetra boundary faces should have no right cell")) {
    return EXIT_FAILURE;
  }
  if(!expect(domain.node(n3).coordinates[2] == 1.0, "node coordinates should be preserved")) {
    return EXIT_FAILURE;
  }
  if(!expect(domain.face(f2).header.entity_group == face_entity_group, "face header should point back to owner entity_group")) {
    return EXIT_FAILURE;
  }
  if(!expect(domain.cell(c0).header.index == 0U, "first cell should retain local index zero")) {
    return EXIT_FAILURE;
  }

  Domain edge_domain("edge_domain");
  const auto edge_node_entity_group = edge_domain.create_entity_group({EntityOrder::node, "nodes"});
  const auto edge_entity_group = edge_domain.create_entity_group(
    {EntityOrder::edge, "feature_edges", invalid_index, true, EntityKind::edge_line}
  );
  const auto edge_face_entity_group = edge_domain.create_entity_group(
    {EntityOrder::face, "surface_faces", invalid_index, true, EntityKind::face_triangle}
  );

  edge_domain.reserve_entity_group_storage(edge_node_entity_group, 3U);
  edge_domain.reserve_entity_group_storage(edge_entity_group, 3U, 6U);
  edge_domain.reserve_entity_group_storage(edge_face_entity_group, 1U, 3U);

  const auto e0n = edge_domain.add_node(edge_node_entity_group, {0.0, 0.0, 0.0});
  const auto e1n = edge_domain.add_node(edge_node_entity_group, {1.0, 0.0, 0.0});
  const auto e2n = edge_domain.add_node(edge_node_entity_group, {0.0, 1.0, 0.0});

  const auto e0 = edge_domain.add_edge(edge_entity_group, {e0n, e1n});
  const auto e1 = edge_domain.add_edge(edge_entity_group, {e1n, e2n});
  const auto e2 = edge_domain.add_edge(edge_entity_group, {e2n, e0n});
  const auto ef0 = edge_domain.add_triangle_face(edge_face_entity_group, {e0n, e1n, e2n});

  edge_domain.set_edge_faces(e0, ef0);
  edge_domain.set_edge_faces(e1, ef0);
  edge_domain.set_edge_faces(e2, ef0);
  edge_domain.set_edge_topology_owner(
    e0,
    {
      sqmesh::geo::TopologyDimension::edge,
      17U,
    }
  );

  const auto edge_summary = edge_domain.summary();
  const auto edge_stats = edge_domain.statistics();
  const auto &edge_storage = edge_domain.entity_group(edge_entity_group).edges();
  const auto edge_nodes = edge_domain.edge_nodes(e0);

  if(!expect(edge_summary.node_count == 3U, "edge domain should report three nodes")) {
    return EXIT_FAILURE;
  }
  if(!expect(edge_summary.edge_count == 3U, "edge domain should report three edges")) {
    return EXIT_FAILURE;
  }
  if(!expect(edge_summary.face_count == 1U, "edge domain should report one face")) {
    return EXIT_FAILURE;
  }
  if(!expect(edge_domain.entity_group_count(EntityOrder::edge) == 1U,
             "edge domain should expose one edge entity_group")) {
    return EXIT_FAILURE;
  }
  if(!expect(edge_storage.size() == 3U, "edge entity_group should own three edges")) {
    return EXIT_FAILURE;
  }
  if(!expect(&edge_storage[0] + 1 == &edge_storage[1], "edge pool should stay contiguous")) {
    return EXIT_FAILURE;
  }
  if(!expect(edge_nodes.size == 2U, "line edge should reference two nodes")) {
    return EXIT_FAILURE;
  }
  if(!expect(edge_nodes[0] == e0n && edge_nodes[1] == e1n,
             "edge connectivity should preserve insertion order")) {
    return EXIT_FAILURE;
  }
  if(!expect(
       edge_domain.adjacent_face(e0, FaceSide::left) == ef0,
       "left edge adjacency should preserve the incident face"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       edge_domain.edge_topology_owner(e0) ==
         sqmesh::geo::TopologyEntityId {sqmesh::geo::TopologyDimension::edge, 17U},
       "edge topology ownership should preserve per-edge source metadata"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       !sqmesh::geo::is_valid(edge_domain.edge_topology_owner(e1)),
       "edges without an assigned topology owner should remain unset"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       !is_valid(edge_domain.adjacent_face(e0, FaceSide::right)),
       "boundary edges should have no right-adjacent face"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(edge_domain.edge(e2).header.entity_group == edge_entity_group,
             "edge header should point back to owner entity_group")) {
    return EXIT_FAILURE;
  }
  if(!expect(edge_stats.boundary_edge_count == 3U,
             "boundary edge count should match inserted boundary edges")) {
    return EXIT_FAILURE;
  }
  if(!expect(edge_stats.interior_edge_count == 0U,
             "single-face edges should not count as interior edges")) {
    return EXIT_FAILURE;
  }
  if(!expect(edge_stats.line_edge_count == 3U,
             "line edge count should match inserted edge entities")) {
    return EXIT_FAILURE;
  }

  Domain contiguous_cells_domain("contiguous_cells");
  const auto contiguous_node_entity_group =
    contiguous_cells_domain.create_entity_group({EntityOrder::node, "nodes"});
  const auto contiguous_face_entity_group = contiguous_cells_domain.create_entity_group(
    {EntityOrder::face, "faces", invalid_index, true, EntityKind::face_triangle}
  );
  const auto contiguous_cell_entity_group = contiguous_cells_domain.create_entity_group(
    {EntityOrder::cell, "cells", invalid_index, false, EntityKind::cell_tetra}
  );

  contiguous_cells_domain.reserve_entity_group_storage(contiguous_node_entity_group, 8U);
  contiguous_cells_domain.reserve_entity_group_storage(contiguous_face_entity_group, 8U, 24U);
  contiguous_cells_domain.reserve_entity_group_storage(contiguous_cell_entity_group, 2U, 8U, 8U);

  const auto a0 = contiguous_cells_domain.add_node(contiguous_node_entity_group, {0.0, 0.0, 0.0});
  const auto a1 = contiguous_cells_domain.add_node(contiguous_node_entity_group, {1.0, 0.0, 0.0});
  const auto a2 = contiguous_cells_domain.add_node(contiguous_node_entity_group, {0.0, 1.0, 0.0});
  const auto a3 = contiguous_cells_domain.add_node(contiguous_node_entity_group, {0.0, 0.0, 1.0});
  const auto a4 = contiguous_cells_domain.add_node(contiguous_node_entity_group, {2.0, 0.0, 0.0});
  const auto a5 = contiguous_cells_domain.add_node(contiguous_node_entity_group, {3.0, 0.0, 0.0});
  const auto a6 = contiguous_cells_domain.add_node(contiguous_node_entity_group, {2.0, 1.0, 0.0});
  const auto a7 = contiguous_cells_domain.add_node(contiguous_node_entity_group, {2.0, 0.0, 1.0});

  const auto g0 = contiguous_cells_domain.add_triangle_face(contiguous_face_entity_group, {a0, a2, a1});
  const auto g1 = contiguous_cells_domain.add_triangle_face(contiguous_face_entity_group, {a0, a1, a3});
  const auto g2 = contiguous_cells_domain.add_triangle_face(contiguous_face_entity_group, {a1, a2, a3});
  const auto g3 = contiguous_cells_domain.add_triangle_face(contiguous_face_entity_group, {a2, a0, a3});
  const auto g4 = contiguous_cells_domain.add_triangle_face(contiguous_face_entity_group, {a4, a6, a5});
  const auto g5 = contiguous_cells_domain.add_triangle_face(contiguous_face_entity_group, {a4, a5, a7});
  const auto g6 = contiguous_cells_domain.add_triangle_face(contiguous_face_entity_group, {a5, a6, a7});
  const auto g7 = contiguous_cells_domain.add_triangle_face(contiguous_face_entity_group, {a6, a4, a7});

  const auto d0 = contiguous_cells_domain.add_tetra_cell(
    contiguous_cell_entity_group,
    {a0, a1, a2, a3},
    {g0, g1, g2, g3}
  );
  const auto d1 = contiguous_cells_domain.add_tetra_cell(
    contiguous_cell_entity_group,
    {a4, a5, a6, a7},
    {g4, g5, g6, g7}
  );

  contiguous_cells_domain.set_face_cells(g0, d0);
  contiguous_cells_domain.set_face_cells(g1, d0);
  contiguous_cells_domain.set_face_cells(g2, d0);
  contiguous_cells_domain.set_face_cells(g3, d0);
  contiguous_cells_domain.set_face_cells(g4, d1);
  contiguous_cells_domain.set_face_cells(g5, d1);
  contiguous_cells_domain.set_face_cells(g6, d1);
  contiguous_cells_domain.set_face_cells(g7, d1);

  const auto &contiguous_cell_storage =
    contiguous_cells_domain.entity_group(contiguous_cell_entity_group).cells();
  if(!expect(
       contiguous_cell_storage.size() == 2U,
       "contiguity domain should own two tetrahedral cells"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       &contiguous_cell_storage[0] + 1 == &contiguous_cell_storage[1],
       "cell pool should stay contiguous across adjacent tetra allocations"
     )) {
    return EXIT_FAILURE;
  }

  Domain remap_domain("entity_group_remap_domain");
  const auto stale_node_entity_group =
    remap_domain.create_entity_group({EntityOrder::node, "surface_old_nodes"});
  const auto remap_node_entity_group =
    remap_domain.create_entity_group({EntityOrder::node, "shared_nodes"});
  const auto stale_face_entity_group = remap_domain.create_entity_group(
    {EntityOrder::face, "surface_old_faces", invalid_index, true, EntityKind::face_triangle}
  );
  const auto remap_face_entity_group = remap_domain.create_entity_group(
    {EntityOrder::face, "volume_faces", invalid_index, true, EntityKind::face_triangle}
  );
  const auto remap_cell_entity_group = remap_domain.create_entity_group(
    {EntityOrder::cell, "volume_cells", invalid_index, false, EntityKind::cell_tetra}
  );

  remap_domain.reserve_entity_group_storage(stale_node_entity_group, 1U);
  remap_domain.reserve_entity_group_storage(remap_node_entity_group, 4U);
  remap_domain.reserve_entity_group_storage(stale_face_entity_group, 1U, 3U);
  remap_domain.reserve_entity_group_storage(remap_face_entity_group, 4U, 12U);
  remap_domain.reserve_entity_group_storage(remap_cell_entity_group, 1U, 4U, 4U);

  static_cast<void>(remap_domain.add_node(stale_node_entity_group, {-1.0, -1.0, -1.0}));
  const auto r0 = remap_domain.add_node(remap_node_entity_group, {0.0, 0.0, 0.0});
  const auto r1 = remap_domain.add_node(remap_node_entity_group, {1.0, 0.0, 0.0});
  const auto r2 = remap_domain.add_node(remap_node_entity_group, {0.0, 1.0, 0.0});
  const auto r3 = remap_domain.add_node(remap_node_entity_group, {0.0, 0.0, 1.0});

  static_cast<void>(remap_domain.add_triangle_face(stale_face_entity_group, {r0, r2, r1}));
  const auto rf0 = remap_domain.add_triangle_face(remap_face_entity_group, {r0, r2, r1});
  const auto rf1 = remap_domain.add_triangle_face(remap_face_entity_group, {r0, r1, r3});
  const auto rf2 = remap_domain.add_triangle_face(remap_face_entity_group, {r1, r2, r3});
  const auto rf3 = remap_domain.add_triangle_face(remap_face_entity_group, {r2, r0, r3});
  const auto rc0 = remap_domain.add_tetra_cell(
    remap_cell_entity_group,
    {r0, r1, r2, r3},
    {rf0, rf1, rf2, rf3}
  );

  remap_domain.set_face_cells(rf0, rc0);
  remap_domain.set_face_cells(rf1, rc0);
  remap_domain.set_face_cells(rf2, rc0);
  remap_domain.set_face_cells(rf3, rc0);

  remap_domain.remove_entity_groups_with_prefix("surface_");

  const sqmesh::mesh::EntityRef remapped_face {1U, 0U};
  const sqmesh::mesh::EntityRef remapped_cell {2U, 0U};
  const auto remapped_face_nodes = remap_domain.face_nodes(remapped_face);
  const auto remapped_cell_nodes = remap_domain.cell_nodes(remapped_cell);
  const auto remapped_cell_faces = remap_domain.cell_faces(remapped_cell);

  if(!expect(remap_domain.entity_group_count() == 3U,
             "prefix entity_group removal should compact surviving entity_groups")) {
    return EXIT_FAILURE;
  }
  if(!expect(remap_domain.entity_group(0U).name() == "shared_nodes" &&
               remap_domain.entity_group(1U).name() == "volume_faces" &&
               remap_domain.entity_group(2U).name() == "volume_cells",
             "surviving entity_groups should be reindexed in-order after prefix removal")) {
    return EXIT_FAILURE;
  }
  if(!expect(remap_domain.node({0U, 0U}).header.entity_group == 0U,
             "surviving node headers should point at the remapped entity_group id")) {
    return EXIT_FAILURE;
  }
  if(!expect(remap_domain.face(remapped_face).header.entity_group == remapped_face.entity_group,
             "surviving face headers should point at the remapped entity_group id")) {
    return EXIT_FAILURE;
  }
  if(!expect(remap_domain.cell(remapped_cell).header.entity_group == remapped_cell.entity_group,
             "surviving cell headers should point at the remapped entity_group id")) {
    return EXIT_FAILURE;
  }
  if(!expect(remapped_face_nodes.size == 3U &&
               remapped_face_nodes[0].entity_group == 0U &&
               remapped_face_nodes[1].entity_group == 0U &&
               remapped_face_nodes[2].entity_group == 0U,
             "face node-channel references should be remapped after prefix removal")) {
    return EXIT_FAILURE;
  }
  if(!expect(remapped_cell_nodes.size == 4U &&
               remapped_cell_nodes[0].entity_group == 0U &&
               remapped_cell_nodes[3].entity_group == 0U,
             "cell node-channel references should preserve valid remapped node entity_groups")) {
    return EXIT_FAILURE;
  }
  if(!expect(remapped_cell_faces.size == 4U &&
               remapped_cell_faces[0].entity_group == 1U &&
               remapped_cell_faces[3].entity_group == 1U,
             "cell face-channel references should be remapped after prefix removal")) {
    return EXIT_FAILURE;
  }
  if(!expect(remap_domain.adjacent_cell(remapped_face, FaceSide::left) == remapped_cell &&
               !is_valid(remap_domain.adjacent_cell(remapped_face, FaceSide::right)),
             "surviving adjacency references should keep pointing at the remapped cell entity_group")) {
    return EXIT_FAILURE;
  }

  Domain interface_domain("interface_domain");
  const auto interface_node_entity_group =
    interface_domain.create_entity_group({EntityOrder::node, "nodes"});
  const auto fluid_boundary_entity_group = interface_domain.create_entity_group(
    {
      EntityOrder::face,
      "fluid_wall",
      301U,
      true,
      EntityKind::face_triangle,
      EntityGroupSemantic::boundary,
      EntityGroupRole::computational,
      101U,
    }
  );
  const auto solid_boundary_entity_group = interface_domain.create_entity_group(
    {
      EntityOrder::face,
      "solid_wall",
      302U,
      true,
      EntityKind::face_triangle,
      EntityGroupSemantic::boundary,
      EntityGroupRole::computational,
      202U,
    }
  );
  const auto interface_entity_group = interface_domain.create_entity_group(
    {
      EntityOrder::face,
      "fluid_solid_interface",
      401U,
      false,
      EntityKind::face_triangle,
      EntityGroupSemantic::interface,
      EntityGroupRole::computational,
      101U,
      202U,
    }
  );
  const auto fluid_cell_entity_group = interface_domain.create_entity_group(
    {
      EntityOrder::cell,
      "fluid_region",
      101U,
      false,
      EntityKind::cell_tetra,
      EntityGroupSemantic::region,
      EntityGroupRole::computational,
      101U,
    }
  );
  const auto solid_cell_entity_group = interface_domain.create_entity_group(
    {
      EntityOrder::cell,
      "solid_region",
      202U,
      false,
      EntityKind::cell_tetra,
      EntityGroupSemantic::region,
      EntityGroupRole::computational,
      202U,
    }
  );

  interface_domain.reserve_entity_group_storage(interface_node_entity_group, 5U);
  interface_domain.reserve_entity_group_storage(fluid_boundary_entity_group, 3U, 9U);
  interface_domain.reserve_entity_group_storage(solid_boundary_entity_group, 3U, 9U);
  interface_domain.reserve_entity_group_storage(interface_entity_group, 1U, 3U);
  interface_domain.reserve_entity_group_storage(fluid_cell_entity_group, 1U, 4U, 4U);
  interface_domain.reserve_entity_group_storage(solid_cell_entity_group, 1U, 4U, 4U);

  const auto i0 = interface_domain.add_node(interface_node_entity_group, {0.0, 0.0, 0.0});
  const auto i1 = interface_domain.add_node(interface_node_entity_group, {1.0, 0.0, 0.0});
  const auto i2 = interface_domain.add_node(interface_node_entity_group, {0.0, 1.0, 0.0});
  const auto i3 = interface_domain.add_node(interface_node_entity_group, {0.0, 0.0, 1.0});
  const auto i4 = interface_domain.add_node(interface_node_entity_group, {1.0, 1.0, 1.0});

  const auto fluid_f0 = interface_domain.add_triangle_face(fluid_boundary_entity_group, {i0, i2, i1});
  const auto fluid_f1 = interface_domain.add_triangle_face(fluid_boundary_entity_group, {i0, i1, i3});
  const auto fluid_f2 = interface_domain.add_triangle_face(fluid_boundary_entity_group, {i2, i0, i3});
  const auto shared_f = interface_domain.add_triangle_face(interface_entity_group, {i1, i2, i3});
  const auto solid_f0 = interface_domain.add_triangle_face(solid_boundary_entity_group, {i4, i1, i2});
  const auto solid_f1 = interface_domain.add_triangle_face(solid_boundary_entity_group, {i4, i3, i1});
  const auto solid_f2 = interface_domain.add_triangle_face(solid_boundary_entity_group, {i4, i2, i3});

  const auto fluid_c = interface_domain.add_tetra_cell(
    fluid_cell_entity_group,
    {i0, i1, i2, i3},
    {fluid_f0, fluid_f1, shared_f, fluid_f2}
  );
  const auto solid_c = interface_domain.add_tetra_cell(
    solid_cell_entity_group,
    {i4, i1, i2, i3},
    {solid_f0, solid_f1, shared_f, solid_f2}
  );

  interface_domain.set_face_cells(fluid_f0, fluid_c);
  interface_domain.set_face_cells(fluid_f1, fluid_c);
  interface_domain.set_face_cells(fluid_f2, fluid_c);
  interface_domain.set_face_cells(shared_f, fluid_c, solid_c);
  interface_domain.set_face_cells(solid_f0, solid_c);
  interface_domain.set_face_cells(solid_f1, solid_c);
  interface_domain.set_face_cells(solid_f2, solid_c);

  if(!expect(interface_domain.entity_group(interface_entity_group).semantic() == EntityGroupSemantic::interface,
             "shared cross-region faces should expose explicit interface semantics")) {
    return EXIT_FAILURE;
  }
  if(!expect(!interface_domain.entity_group(interface_entity_group).is_boundary(),
             "interface entity_groups should stay distinct from boundary-only entity_groups")) {
    return EXIT_FAILURE;
  }
  if(!expect(interface_domain.entity_group(interface_entity_group).primary_region_zone_id() == 101U &&
               interface_domain.entity_group(interface_entity_group).secondary_region_zone_id() == 202U,
             "interface entity_groups should record both adjacent region zone ids")) {
    return EXIT_FAILURE;
  }
  if(!expect(interface_domain.entity_group(fluid_boundary_entity_group).semantic() == EntityGroupSemantic::boundary &&
               interface_domain.entity_group(fluid_boundary_entity_group).primary_region_zone_id() == 101U,
             "boundary entity_groups should preserve their owning region zone id")) {
    return EXIT_FAILURE;
  }
  if(!expect(interface_domain.entity_group(solid_cell_entity_group).semantic() == EntityGroupSemantic::region &&
               interface_domain.entity_group(solid_cell_entity_group).primary_region_zone_id() == 202U,
             "region entity_groups should preserve their explicit region zone id")) {
    return EXIT_FAILURE;
  }
  if(!expect(is_interface_entity(interface_domain.face(shared_f).header.flags),
             "interface faces should carry the explicit interface entity flag")) {
    return EXIT_FAILURE;
  }
  if(!expect(!is_interface_entity(interface_domain.face(fluid_f0).header.flags),
             "boundary faces should not carry the interface entity flag")) {
    return EXIT_FAILURE;
  }
  if(!expect(interface_domain.adjacent_cell(shared_f, FaceSide::left) == fluid_c &&
               interface_domain.adjacent_cell(shared_f, FaceSide::right) == solid_c,
             "interface faces should preserve left/right region adjacency")) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
