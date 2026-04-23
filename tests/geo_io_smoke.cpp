// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "sqmesh/sqmesh.hpp"

#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <filesystem>
#include <string>

namespace {

bool expect(bool condition, const char *message)
{
  if(condition) {
    return true;
  }

  std::fprintf(stderr, "geo_io_smoke: %s\n", message);
  return false;
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
         sqmesh::geo::import_step("missing.step", model_handle, {}, context) ==
           sqmesh::base::StatusCode::unsupported,
         "STEP import should report unsupported when CAD IO is disabled"
       )) {
      return EXIT_FAILURE;
    }
    if(!expect(
         sqmesh::geo::export_iges(
           sqmesh::invalid_handle,
           "missing.igs",
           {},
           context
         ) == sqmesh::base::StatusCode::unsupported,
         "IGES export should report unsupported when CAD IO is disabled"
       )) {
      return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
  }

#ifndef sqmesh_TEST_OCC_STEP_SAMPLE_PATH
  std::fprintf(
    stderr,
    "geo_io_smoke: OCC is enabled but no STEP sample path was configured.\n"
  );
  return EXIT_FAILURE;
#else
  namespace fs = std::filesystem;

  const fs::path step_input = sqmesh_TEST_OCC_STEP_SAMPLE_PATH;
  if(!expect(
       fs::exists(step_input),
       "configured STEP sample path should exist"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::ModelHandle source_model = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::geo::import_step(step_input.string(), source_model, {}, context) ==
         sqmesh::base::StatusCode::ok,
       "STEP import should produce a geometry model handle"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::ModelSummary source_summary;
  if(!expect(
       sqmesh::geo::model_summary(source_model, source_summary, context) ==
         sqmesh::base::StatusCode::ok,
       "model_summary should resolve the imported STEP model"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       source_summary.face_count > 0U && source_summary.edge_count > 0U,
       "imported STEP sample should expose non-zero face and edge counts"
     )) {
    return EXIT_FAILURE;
  }

  const fs::path temp_dir =
    fs::temp_directory_path() /
    ("sqmesh_occ_smoke_" +
     std::to_string(
       std::chrono::steady_clock::now().time_since_epoch().count()
     ));
  fs::create_directories(temp_dir);

  const fs::path iges_output = temp_dir / "roundtrip.igs";
  const fs::path step_output = temp_dir / "roundtrip.step";

  sqmesh::geo::IgesExportOptions iges_options;
  iges_options.unit_name = "MM";
  iges_options.write_mode = sqmesh::geo::IgesWriteMode::brep;

  if(!expect(
       sqmesh::geo::export_iges(source_model, iges_output.string(), iges_options, context) ==
         sqmesh::base::StatusCode::ok,
       "IGES export should write the OCC-backed model"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::ModelHandle iges_model = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::geo::import_iges(iges_output.string(), iges_model, {}, context) ==
         sqmesh::base::StatusCode::ok,
       "IGES import should reload the exported geometry"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::ModelSummary iges_summary;
  if(!expect(
       sqmesh::geo::model_summary(iges_model, iges_summary, context) ==
         sqmesh::base::StatusCode::ok,
       "model_summary should resolve the imported IGES model"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       iges_summary.face_count == source_summary.face_count &&
         iges_summary.edge_count == source_summary.edge_count,
       "STEP -> IGES round-trip should preserve face and edge counts"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::StepExportOptions step_options;
  step_options.linear_tolerance = 1.0e-6;
  step_options.unit_name = "MM";
  step_options.schema = "AP214IS";

  if(!expect(
       sqmesh::geo::export_step(iges_model, step_output.string(), step_options, context) ==
         sqmesh::base::StatusCode::ok,
       "STEP export should write the IGES-imported OCC model"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::ModelHandle step_roundtrip_model = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::geo::import_step(step_output.string(), step_roundtrip_model, {}, context) ==
         sqmesh::base::StatusCode::ok,
       "STEP import should reload the exported STEP geometry"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::ModelSummary step_roundtrip_summary;
  if(!expect(
       sqmesh::geo::model_summary(step_roundtrip_model, step_roundtrip_summary, context) ==
         sqmesh::base::StatusCode::ok,
       "model_summary should resolve the STEP round-trip model"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       step_roundtrip_summary.face_count == source_summary.face_count &&
         step_roundtrip_summary.edge_count == source_summary.edge_count,
       "STEP -> IGES -> STEP round-trip should preserve face and edge counts"
     )) {
    return EXIT_FAILURE;
  }

  fs::remove_all(temp_dir);
  return EXIT_SUCCESS;
#endif
}
