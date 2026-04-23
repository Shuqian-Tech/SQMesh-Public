// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "sqmesh/sqmesh.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>

namespace {

constexpr double kSqrt3Over2 = 0.86602540378443864676372317075293618;
constexpr double kSqrt2Over3 = 0.8164965809277260327324280249019638;
constexpr double kIdealTriangleAngleDegrees = 60.0;
constexpr double kIdealTetraDihedralDegrees = 70.52877936550931;
constexpr double kTolerance = 1.0e-9;

bool expect(bool condition, const char *message)
{
  if(condition) {
    return true;
  }

  std::fprintf(stderr, "mesh_quality_smoke: %s\n", message);
  return false;
}

bool expect_near(double actual, double expected, double tolerance, const char *message)
{
  return expect(std::abs(actual - expected) <= tolerance, message);
}

const sqmesh::mesh::KindQualitySummary *find_kind_summary(
  const sqmesh::mesh::MeshQualityReport &report,
  sqmesh::mesh::EntityKind kind
)
{
  for(const auto &summary : report.kinds) {
    if(summary.kind == kind) {
      return &summary;
    }
  }
  return nullptr;
}

bool write_text_file(
  const std::filesystem::path &path,
  const char *contents
)
{
  std::ofstream output(path, std::ios::trunc);
  if(!output) {
    std::fprintf(
      stderr,
      "mesh_quality_smoke: failed to open %s for writing\n",
      path.string().c_str()
    );
    return false;
  }

  output << contents;
  return static_cast<bool>(output);
}

} // namespace

int main()
{
  using sqmesh::mesh::Domain;
  using sqmesh::mesh::ElementQuality;
  using sqmesh::mesh::EntityKind;
  using sqmesh::mesh::EntityOrder;
  using sqmesh::mesh::MeshHandle;
  using sqmesh::mesh::MeshQualityReport;
  using sqmesh::mesh::QualityStatus;
  using sqmesh::mesh::invalid_index;

  {
    Domain domain("surface_quality");
    const auto node_entity_group = domain.create_entity_group({EntityOrder::node, "nodes"});
    const auto face_entity_group = domain.create_entity_group(
      {EntityOrder::face, "surface_faces", invalid_index, true, EntityKind::face_triangle}
    );

    domain.reserve_entity_group_storage(node_entity_group, 3U);
    domain.reserve_entity_group_storage(face_entity_group, 1U, 3U);

    const auto n0 = domain.add_node(node_entity_group, {0.0, 0.0, 0.0});
    const auto n1 = domain.add_node(node_entity_group, {1.0, 0.0, 0.0});
    const auto n2 = domain.add_node(node_entity_group, {0.5, kSqrt3Over2, 0.0});
    const auto f0 = domain.add_triangle_face(face_entity_group, {n0, n1, n2});

    const ElementQuality quality = domain.element_quality(f0);
    if(!expect(quality.supported, "triangle quality should be supported")) {
      return EXIT_FAILURE;
    }
    if(!expect(quality.status == QualityStatus::valid, "equilateral triangle should be valid")) {
      return EXIT_FAILURE;
    }
    if(!expect(!quality.degenerate && !quality.inverted,
               "equilateral triangle should not be degenerate or inverted")) {
      return EXIT_FAILURE;
    }
    if(!expect_near(quality.jacobian, 1.0, kTolerance,
                    "equilateral triangle should have unit Jacobian quality")) {
      return EXIT_FAILURE;
    }
    if(!expect_near(quality.skewness, 0.0, kTolerance,
                    "equilateral triangle should have zero skewness")) {
      return EXIT_FAILURE;
    }
    if(!expect_near(quality.aspect_ratio, 1.0, kTolerance,
                    "equilateral triangle should have unit aspect ratio")) {
      return EXIT_FAILURE;
    }
    if(!expect_near(quality.radius_ratio, 1.0, kTolerance,
                    "equilateral triangle should have unit radius ratio")) {
      return EXIT_FAILURE;
    }
    if(!expect_near(quality.min_angle, kIdealTriangleAngleDegrees, kTolerance,
                    "equilateral triangle should have a 60 degree minimum angle")) {
      return EXIT_FAILURE;
    }
    if(!expect_near(quality.max_angle, kIdealTriangleAngleDegrees, kTolerance,
                    "equilateral triangle should have a 60 degree maximum angle")) {
      return EXIT_FAILURE;
    }

    const auto report = domain.quality_report();
    if(!expect(report.supported_element_count == 1U,
               "surface report should contain one supported element")) {
      return EXIT_FAILURE;
    }
    if(!expect(report.valid_element_count == 1U,
               "surface report should count one valid element")) {
      return EXIT_FAILURE;
    }

    const auto *triangle_summary = find_kind_summary(report, EntityKind::face_triangle);
    if(!expect(triangle_summary != nullptr,
               "surface report should include the triangle kind summary")) {
      return EXIT_FAILURE;
    }
    if(!expect_near(triangle_summary->jacobian.average, 1.0, kTolerance,
                    "triangle summary should average to unit Jacobian quality")) {
      return EXIT_FAILURE;
    }
    if(!expect_near(triangle_summary->radius_ratio.average, 1.0, kTolerance,
                    "triangle summary should average to unit radius ratio")) {
      return EXIT_FAILURE;
    }
    if(!expect_near(triangle_summary->min_angle.average, kIdealTriangleAngleDegrees, kTolerance,
                    "triangle summary should preserve the 60 degree minimum angle")) {
      return EXIT_FAILURE;
    }
    if(!expect_near(triangle_summary->max_angle.average, kIdealTriangleAngleDegrees, kTolerance,
                    "triangle summary should preserve the 60 degree maximum angle")) {
      return EXIT_FAILURE;
    }
  }

  {
    Domain domain("skewed_surface_angle_quality");
    const auto node_entity_group = domain.create_entity_group({EntityOrder::node, "nodes"});
    const auto face_entity_group = domain.create_entity_group(
      {EntityOrder::face, "surface_faces", invalid_index, true, EntityKind::face_triangle}
    );

    domain.reserve_entity_group_storage(node_entity_group, 3U);
    domain.reserve_entity_group_storage(face_entity_group, 1U, 3U);

    const auto n0 = domain.add_node(node_entity_group, {0.0, 0.0, 0.0});
    const auto n1 = domain.add_node(node_entity_group, {1.0, 0.0, 0.0});
    const auto n2 = domain.add_node(node_entity_group, {0.99, 0.01, 0.0});
    const auto f0 = domain.add_triangle_face(face_entity_group, {n0, n1, n2});

    const ElementQuality quality = domain.element_quality(f0);
    if(!expect(quality.status == QualityStatus::valid,
               "near-collinear triangle should remain valid when the area is non-zero")) {
      return EXIT_FAILURE;
    }
    if(!expect(quality.min_angle < 1.0,
               "near-collinear triangle should expose a very small minimum angle")) {
      return EXIT_FAILURE;
    }
    if(!expect(quality.max_angle > 130.0,
               "near-collinear triangle should expose an obtuse maximum angle")) {
      return EXIT_FAILURE;
    }
    if(!expect(quality.radius_ratio < 0.05,
               "near-collinear triangle should expose a poor radius ratio")) {
      return EXIT_FAILURE;
    }

    const auto report = domain.quality_report();
    const auto *triangle_summary = find_kind_summary(report, EntityKind::face_triangle);
    if(!expect(triangle_summary != nullptr,
               "skewed triangle report should include the triangle kind summary")) {
      return EXIT_FAILURE;
    }
    if(!expect(triangle_summary->min_angle.minimum < 1.0,
               "triangle angle summary should preserve the small minimum angle")) {
      return EXIT_FAILURE;
    }
    if(!expect(triangle_summary->max_angle.maximum > 130.0,
               "triangle angle summary should preserve the large maximum angle")) {
      return EXIT_FAILURE;
    }
    if(!expect(triangle_summary->radius_ratio.minimum < 0.05,
               "triangle radius-ratio summary should preserve the poor shape")) {
      return EXIT_FAILURE;
    }
  }

  {
    Domain domain("degenerate_surface_quality");
    const auto node_entity_group = domain.create_entity_group({EntityOrder::node, "nodes"});
    const auto face_entity_group = domain.create_entity_group(
      {EntityOrder::face, "surface_faces", invalid_index, true, EntityKind::face_triangle}
    );

    domain.reserve_entity_group_storage(node_entity_group, 3U);
    domain.reserve_entity_group_storage(face_entity_group, 1U, 3U);

    const auto n0 = domain.add_node(node_entity_group, {0.0, 0.0, 0.0});
    const auto n1 = domain.add_node(node_entity_group, {1.0, 0.0, 0.0});
    const auto n2 = domain.add_node(node_entity_group, {2.0, 0.0, 0.0});
    const auto f0 = domain.add_triangle_face(face_entity_group, {n0, n1, n2});

    const ElementQuality quality = domain.element_quality(f0);
    if(!expect(quality.status == QualityStatus::degenerate,
               "collinear triangle should be reported as degenerate")) {
      return EXIT_FAILURE;
    }
    if(!expect(quality.degenerate && !quality.inverted,
               "degenerate triangle should not be reported as inverted")) {
      return EXIT_FAILURE;
    }
    if(!expect(quality.jacobian == 0.0,
               "degenerate triangle should have zero Jacobian quality")) {
      return EXIT_FAILURE;
    }
    if(!expect(quality.skewness == 1.0,
               "degenerate triangle should clamp skewness to one")) {
      return EXIT_FAILURE;
    }
    if(!expect(std::isinf(quality.aspect_ratio),
               "degenerate triangle should have infinite aspect ratio")) {
      return EXIT_FAILURE;
    }
    if(!expect(quality.radius_ratio == 0.0,
               "degenerate triangle should have zero radius ratio")) {
      return EXIT_FAILURE;
    }
    if(!expect(std::isnan(quality.min_angle) && std::isnan(quality.max_angle),
               "degenerate triangle should leave angle metrics undefined")) {
      return EXIT_FAILURE;
    }

    const auto report = domain.quality_report();
    if(!expect(report.degenerate_element_count == 1U,
               "degenerate triangle report should still count one bad element")) {
      return EXIT_FAILURE;
    }
    const auto *triangle_summary = find_kind_summary(report, EntityKind::face_triangle);
    if(!expect(triangle_summary != nullptr,
               "degenerate triangle report should still include the triangle kind summary")) {
      return EXIT_FAILURE;
    }
    if(!expect(triangle_summary->min_angle.count == 0U &&
                 triangle_summary->max_angle.count == 0U,
               "degenerate triangle angle summaries should not hide the bad-element count with fake values")) {
      return EXIT_FAILURE;
    }
    if(!expect(triangle_summary->radius_ratio.count == 1U &&
                 triangle_summary->radius_ratio.maximum == 0.0,
               "degenerate triangle radius-ratio summaries should preserve the bad element with zero shape quality")) {
      return EXIT_FAILURE;
    }
  }

  {
    Domain domain("volume_quality");
    const auto node_entity_group = domain.create_entity_group({EntityOrder::node, "nodes"});
    const auto face_entity_group = domain.create_entity_group(
      {EntityOrder::face, "boundary_faces", invalid_index, true, EntityKind::face_triangle}
    );
    const auto cell_entity_group = domain.create_entity_group(
      {EntityOrder::cell, "cells", invalid_index, false, EntityKind::cell_tetra}
    );

    domain.reserve_entity_group_storage(node_entity_group, 4U);
    domain.reserve_entity_group_storage(face_entity_group, 4U, 12U);
    domain.reserve_entity_group_storage(cell_entity_group, 1U, 4U, 4U);

    const auto n0 = domain.add_node(node_entity_group, {0.0, 0.0, 0.0});
    const auto n1 = domain.add_node(node_entity_group, {1.0, 0.0, 0.0});
    const auto n2 = domain.add_node(node_entity_group, {0.5, kSqrt3Over2, 0.0});
    const auto n3 = domain.add_node(node_entity_group, {0.5, kSqrt3Over2 / 3.0, kSqrt2Over3});

    const auto f0 = domain.add_triangle_face(face_entity_group, {n0, n2, n1});
    const auto f1 = domain.add_triangle_face(face_entity_group, {n0, n1, n3});
    const auto f2 = domain.add_triangle_face(face_entity_group, {n1, n2, n3});
    const auto f3 = domain.add_triangle_face(face_entity_group, {n2, n0, n3});
    const auto c0 = domain.add_tetra_cell(cell_entity_group, {n0, n1, n2, n3}, {f0, f1, f2, f3});

    domain.set_face_cells(f0, c0);
    domain.set_face_cells(f1, c0);
    domain.set_face_cells(f2, c0);
    domain.set_face_cells(f3, c0);

    const ElementQuality quality = domain.element_quality(c0);
    if(!expect(quality.status == QualityStatus::valid,
               "regular tetrahedron should be valid")) {
      return EXIT_FAILURE;
    }
    if(!expect(!quality.degenerate && !quality.inverted,
               "regular tetrahedron should not be degenerate or inverted")) {
      return EXIT_FAILURE;
    }
    if(!expect_near(quality.jacobian, 1.0, kTolerance,
                    "regular tetrahedron should have unit Jacobian quality")) {
      return EXIT_FAILURE;
    }
    if(!expect_near(quality.skewness, 0.0, kTolerance,
                    "regular tetrahedron should have zero skewness")) {
      return EXIT_FAILURE;
    }
    if(!expect_near(quality.aspect_ratio, 1.0, kTolerance,
                    "regular tetrahedron should have unit aspect ratio")) {
      return EXIT_FAILURE;
    }
    if(!expect_near(quality.radius_ratio, 1.0, kTolerance,
                    "regular tetrahedron should have unit radius ratio")) {
      return EXIT_FAILURE;
    }
    if(!expect_near(quality.min_angle, kIdealTetraDihedralDegrees, kTolerance,
                    "regular tetrahedron should have the ideal minimum dihedral angle")) {
      return EXIT_FAILURE;
    }
    if(!expect_near(quality.max_angle, kIdealTetraDihedralDegrees, kTolerance,
                    "regular tetrahedron should have the ideal maximum dihedral angle")) {
      return EXIT_FAILURE;
    }

    const auto report = domain.quality_report();
    if(!expect(report.supported_element_count == 5U,
               "tetra quality report should include four faces and one cell")) {
      return EXIT_FAILURE;
    }
    if(!expect(report.valid_element_count == 5U,
               "regular tetra quality report should count five valid elements")) {
      return EXIT_FAILURE;
    }

    const auto *tetra_summary = find_kind_summary(report, EntityKind::cell_tetra);
    if(!expect(tetra_summary != nullptr,
               "tetra quality report should include a tetra kind summary")) {
      return EXIT_FAILURE;
    }
    if(!expect_near(tetra_summary->jacobian.minimum, 1.0, kTolerance,
                    "tetra summary should preserve unit Jacobian quality")) {
      return EXIT_FAILURE;
    }
    if(!expect_near(tetra_summary->radius_ratio.average, 1.0, kTolerance,
                    "tetra summary should preserve unit radius ratio")) {
      return EXIT_FAILURE;
    }
    if(!expect_near(tetra_summary->min_angle.average, kIdealTetraDihedralDegrees, kTolerance,
                    "tetra summary should preserve the ideal minimum dihedral angle")) {
      return EXIT_FAILURE;
    }
    if(!expect_near(tetra_summary->max_angle.average, kIdealTetraDihedralDegrees, kTolerance,
                    "tetra summary should preserve the ideal maximum dihedral angle")) {
      return EXIT_FAILURE;
    }
  }

  {
    Domain domain("sliver_volume_angle_quality");
    const auto node_entity_group = domain.create_entity_group({EntityOrder::node, "nodes"});
    const auto face_entity_group = domain.create_entity_group(
      {EntityOrder::face, "boundary_faces", invalid_index, true, EntityKind::face_triangle}
    );
    const auto cell_entity_group = domain.create_entity_group(
      {EntityOrder::cell, "cells", invalid_index, false, EntityKind::cell_tetra}
    );

    domain.reserve_entity_group_storage(node_entity_group, 4U);
    domain.reserve_entity_group_storage(face_entity_group, 4U, 12U);
    domain.reserve_entity_group_storage(cell_entity_group, 1U, 4U, 4U);

    const auto n0 = domain.add_node(node_entity_group, {0.0, 0.0, 0.0});
    const auto n1 = domain.add_node(node_entity_group, {1.0, 0.0, 0.0});
    const auto n2 = domain.add_node(node_entity_group, {0.5, kSqrt3Over2, 0.0});
    const auto n3 = domain.add_node(node_entity_group, {0.5, kSqrt3Over2 / 3.0, 0.05});

    const auto f0 = domain.add_triangle_face(face_entity_group, {n0, n2, n1});
    const auto f1 = domain.add_triangle_face(face_entity_group, {n0, n1, n3});
    const auto f2 = domain.add_triangle_face(face_entity_group, {n1, n2, n3});
    const auto f3 = domain.add_triangle_face(face_entity_group, {n2, n0, n3});
    const auto c0 = domain.add_tetra_cell(cell_entity_group, {n0, n1, n2, n3}, {f0, f1, f2, f3});

    domain.set_face_cells(f0, c0);
    domain.set_face_cells(f1, c0);
    domain.set_face_cells(f2, c0);
    domain.set_face_cells(f3, c0);

    const ElementQuality quality = domain.element_quality(c0);
    if(!expect(quality.status == QualityStatus::valid,
               "flattened tetrahedron should remain valid while exposing poor dihedral angles")) {
      return EXIT_FAILURE;
    }
    if(!expect(quality.min_angle < 10.0,
               "flattened tetrahedron should expose a small minimum dihedral angle")) {
      return EXIT_FAILURE;
    }
    if(!expect(quality.max_angle > 160.0,
               "flattened tetrahedron should expose a large maximum dihedral angle")) {
      return EXIT_FAILURE;
    }
    if(!expect(quality.radius_ratio < 0.2,
               "flattened tetrahedron should expose a poor radius ratio")) {
      return EXIT_FAILURE;
    }

    const auto report = domain.quality_report();
    const auto *tetra_summary = find_kind_summary(report, EntityKind::cell_tetra);
    if(!expect(tetra_summary != nullptr,
               "sliver tetrahedron report should include the tetra kind summary")) {
      return EXIT_FAILURE;
    }
    if(!expect(tetra_summary->min_angle.minimum < 10.0,
               "tetra angle summary should preserve the small minimum dihedral angle")) {
      return EXIT_FAILURE;
    }
    if(!expect(tetra_summary->max_angle.maximum > 160.0,
               "tetra angle summary should preserve the large maximum dihedral angle")) {
      return EXIT_FAILURE;
    }
    if(!expect(tetra_summary->radius_ratio.minimum < 0.2,
               "tetra radius-ratio summary should preserve the poor sliver quality")) {
      return EXIT_FAILURE;
    }
  }

  {
    Domain domain("inverted_volume_quality");
    const auto node_entity_group = domain.create_entity_group({EntityOrder::node, "nodes"});
    const auto face_entity_group = domain.create_entity_group(
      {EntityOrder::face, "boundary_faces", invalid_index, true, EntityKind::face_triangle}
    );
    const auto cell_entity_group = domain.create_entity_group(
      {EntityOrder::cell, "cells", invalid_index, false, EntityKind::cell_tetra}
    );

    domain.reserve_entity_group_storage(node_entity_group, 4U);
    domain.reserve_entity_group_storage(face_entity_group, 4U, 12U);
    domain.reserve_entity_group_storage(cell_entity_group, 1U, 4U, 4U);

    const auto n0 = domain.add_node(node_entity_group, {0.0, 0.0, 0.0});
    const auto n1 = domain.add_node(node_entity_group, {1.0, 0.0, 0.0});
    const auto n2 = domain.add_node(node_entity_group, {0.5, kSqrt3Over2, 0.0});
    const auto n3 = domain.add_node(node_entity_group, {0.5, kSqrt3Over2 / 3.0, kSqrt2Over3});

    const auto f0 = domain.add_triangle_face(face_entity_group, {n0, n2, n1});
    const auto f1 = domain.add_triangle_face(face_entity_group, {n0, n1, n3});
    const auto f2 = domain.add_triangle_face(face_entity_group, {n1, n2, n3});
    const auto f3 = domain.add_triangle_face(face_entity_group, {n2, n0, n3});
    const auto c0 = domain.add_tetra_cell(cell_entity_group, {n0, n2, n1, n3}, {f0, f1, f2, f3});

    domain.set_face_cells(f0, c0);
    domain.set_face_cells(f1, c0);
    domain.set_face_cells(f2, c0);
    domain.set_face_cells(f3, c0);

    const ElementQuality quality = domain.element_quality(c0);
    if(!expect(quality.status == QualityStatus::inverted,
               "swapped tetra winding should be reported as inverted")) {
      return EXIT_FAILURE;
    }
    if(!expect(!quality.degenerate && quality.inverted,
               "inverted tetrahedron should not be reported as degenerate")) {
      return EXIT_FAILURE;
    }
    if(!expect(quality.jacobian < -0.99,
               "inverted regular tetrahedron should keep a strongly negative Jacobian")) {
      return EXIT_FAILURE;
    }
    if(!expect_near(quality.skewness, 0.0, kTolerance,
                    "tetra skewness should be geometric and stay ideal after inversion")) {
      return EXIT_FAILURE;
    }
    if(!expect_near(quality.aspect_ratio, 1.0, kTolerance,
                    "tetra aspect ratio should stay ideal after inversion")) {
      return EXIT_FAILURE;
    }
    if(!expect_near(quality.radius_ratio, 1.0, kTolerance,
                    "tetra radius ratio should stay geometric after inversion")) {
      return EXIT_FAILURE;
    }
    if(!expect_near(quality.min_angle, kIdealTetraDihedralDegrees, kTolerance,
                    "tetra minimum dihedral angle should stay geometric after inversion")) {
      return EXIT_FAILURE;
    }
    if(!expect_near(quality.max_angle, kIdealTetraDihedralDegrees, kTolerance,
                    "tetra maximum dihedral angle should stay geometric after inversion")) {
      return EXIT_FAILURE;
    }

    const auto report = domain.quality_report();
    if(!expect(report.inverted_element_count == 1U,
               "inverted tetra report should count one inverted element")) {
      return EXIT_FAILURE;
    }
    if(!expect(report.valid_element_count == 4U,
               "inverted tetra report should still count the four valid boundary faces")) {
      return EXIT_FAILURE;
    }
    const auto *tetra_summary = find_kind_summary(report, EntityKind::cell_tetra);
    if(!expect(tetra_summary != nullptr &&
                 tetra_summary->inverted_element_count == 1U &&
                 std::abs(tetra_summary->radius_ratio.average - 1.0) <= kTolerance,
               "inverted tetra summaries should preserve the geometric radius ratio without hiding inversion")) {
      return EXIT_FAILURE;
    }
  }

  {
    sqmesh::base::ContextHandle context = sqmesh::invalid_handle;
    if(!expect(
         sqmesh::base::initialize(context) == sqmesh::base::StatusCode::ok,
         "initialize should create a runtime context for handle-based quality queries"
       )) {
      return EXIT_FAILURE;
    }

    const auto temp_path =
      std::filesystem::temp_directory_path() / "sqmesh_mesh_quality_smoke.obj";
    constexpr const char *kObjContents =
      "o tri_patch\n"
      "v 0.0 0.0 0.0\n"
      "v 1.0 0.0 0.0\n"
      "v 0.5 0.8660254037844386 0.0\n"
      "g feature_loop\n"
      "l 1 2\n"
      "l 2 3\n"
      "l 3 1\n"
      "g surface_patch\n"
      "f 1 3 2\n";

    if(!expect(write_text_file(temp_path, kObjContents),
               "quality smoke should write a temporary OBJ fixture")) {
      return EXIT_FAILURE;
    }

    MeshHandle mesh_handle = sqmesh::invalid_handle;
    if(!expect(
         sqmesh::mesh::import_obj(temp_path.string(), mesh_handle, {}, context) ==
           sqmesh::base::StatusCode::ok,
         "import_obj should create a mesh handle for quality reporting"
       )) {
      std::filesystem::remove(temp_path);
      return EXIT_FAILURE;
    }

    MeshQualityReport report;
    if(!expect(
         sqmesh::mesh::mesh_quality_report(mesh_handle, report, context) ==
           sqmesh::base::StatusCode::ok,
         "mesh_quality_report should evaluate the imported OBJ mesh"
       )) {
      std::filesystem::remove(temp_path);
      return EXIT_FAILURE;
    }

    if(!expect(report.supported_element_count == 1U,
               "handle-based quality reporting should evaluate the imported face")) {
      std::filesystem::remove(temp_path);
      return EXIT_FAILURE;
    }
    if(!expect(report.valid_element_count == 1U,
               "handle-based quality reporting should preserve the valid triangle")) {
      std::filesystem::remove(temp_path);
      return EXIT_FAILURE;
    }
    const auto *triangle_summary = find_kind_summary(report, EntityKind::face_triangle);
    if(!expect(triangle_summary != nullptr &&
                 std::abs(triangle_summary->radius_ratio.average - 1.0) <=
                   kTolerance &&
                 std::abs(triangle_summary->min_angle.average - kIdealTriangleAngleDegrees) <=
                   kTolerance &&
                 std::abs(triangle_summary->max_angle.average - kIdealTriangleAngleDegrees) <=
                   kTolerance,
               "handle-based quality reporting should preserve the public triangle quality semantics")) {
      std::filesystem::remove(temp_path);
      return EXIT_FAILURE;
    }

    std::filesystem::remove(temp_path);
    if(!expect(
         sqmesh::base::shutdown_all() == sqmesh::base::StatusCode::ok,
         "shutdown_all should release the handle-based quality test context"
       )) {
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}
