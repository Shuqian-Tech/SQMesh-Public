// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
// Minimal surface-only example:
//   import STEP/IGES -> Auto CFD Surface Mesher -> export OBJ
//
// Usage:
//   surface_mesh_example <model.step|model.iges> [min_len=1.0] [max_len=20.0]

#include "sqmesh/sqmesh.hpp"
#include "../src/core/log.hpp"

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>

namespace {

bool ends_with_ci(const std::string& s, const std::string& suffix) {
    if (s.size() < suffix.size()) return false;
    for (std::size_t i = 0; i < suffix.size(); ++i) {
        const char a = s[s.size() - suffix.size() + i];
        const char b = suffix[i];
        const char al = (a >= 'A' && a <= 'Z') ? static_cast<char>(a + 32) : a;
        const char bl = (b >= 'A' && b <= 'Z') ? static_cast<char>(b + 32) : b;
        if (al != bl) return false;
    }
    return true;
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0]
                  << " <model.step|model.iges> [min_len=1.0] [max_len=20.0]"
                  << " [dist_angle=15] [growth=1.2] [proximity:0|1=1]\n";
        return 1;
    }

    const char* file_path = argv[1];
    const double min_len = argc > 2 ? std::atof(argv[2]) : 1.0;
    const double max_len = argc > 3 ? std::atof(argv[3]) : 20.0;
    const double dist_angle = argc > 4 ? std::atof(argv[4]) : 18;
    const double growth_rate = argc > 5 ? std::atof(argv[5]) : 1.2;
    const bool   proximity = argc > 6 ? (std::atoi(argv[6]) != 0) : false;

    const std::string log_path =
        std::filesystem::path(file_path).stem().string() + ".log";
    sqmesh::log::init(log_path, "info");

    sqmesh::base::ContextHandle context;
    if (sqmesh::base::initialize(context) != sqmesh::base::StatusCode::ok) {
        std::cerr << "[ERROR] SQMesh initialization failed.\n";
        return 1;
    }

    const std::string path_str(file_path);
    const bool is_iges =
        ends_with_ci(path_str, ".iges") || ends_with_ci(path_str, ".igs");

    sqmesh::geo::ModelHandle model;
    std::cout << "Loading " << (is_iges ? "IGES" : "STEP")
              << " model: " << file_path << " ...\n";
    const auto import_status = is_iges
        ? sqmesh::geo::import_iges(file_path, model, {}, context)
        : sqmesh::geo::import_step(file_path, model, {}, context);
    if (import_status != sqmesh::base::StatusCode::ok) {
        std::cerr << "[ERROR] Failed to import geometry. status="
                  << static_cast<int>(import_status) << "\n";
        return 1;
    }

    sqmesh::mesh::MeshingOptions meshing_options;
    meshing_options.parameters.set_number("minimum_length", min_len);
    meshing_options.parameters.set_number("maximum_length", max_len);
    meshing_options.parameters.set_number("distortion_angle", dist_angle);
    meshing_options.parameters.set_number("growth_rate", growth_rate);
    meshing_options.parameters.set_boolean("proximity", proximity);
    meshing_options.parameters.set_boolean("allow_quality_gate_failure", true);
    meshing_options.parameters.set_boolean("allow_final_screen_failure", true);

    std::cout << "\n======== Auto CFD Surface Mesher ========\n";
    std::cout << " min_length = " << min_len << "\n";
    std::cout << " max_length = " << max_len << "\n";
    std::cout << "=========================================\n";

    sqmesh::mesh::MeshHandle surface_mesh;
    const auto t0 = std::chrono::high_resolution_clock::now();
    const auto status = sqmesh::mesh::create_surface_mesh(
        model, "Auto CFD Surface Mesher", meshing_options, surface_mesh, context);
    const auto t1 = std::chrono::high_resolution_clock::now();
    const double elapsed_ms =
        std::chrono::duration<double, std::milli>(t1 - t0).count();

    if (status != sqmesh::base::StatusCode::ok) {
        std::cerr << "[ERROR] Surface meshing failed. status="
                  << static_cast<int>(status) << " ("
                  << sqmesh::base::status_code_name(status) << ")\n";
        const auto err = sqmesh::base::last_error_message();
        if (!err.empty()) std::cerr << "  " << err << "\n";
        return 1;
    }
    std::cout << "[OK] Surface meshing succeeded! Elapsed: " << elapsed_ms
              << " ms\n";

    sqmesh::mesh::MeshSummary summary;
    if (sqmesh::mesh::mesh_summary(surface_mesh, summary, context) ==
        sqmesh::base::StatusCode::ok) {
        std::cout << "  Nodes: " << summary.node_count << "\n";
        std::cout << "  Edges: " << summary.edge_count << "\n";
        std::cout << "  Faces: " << summary.face_count << "\n";
    }

    const std::string out_path = "surface_mesh.obj";
    sqmesh::mesh::ObjExportOptions export_opts;
    const auto export_status =
        sqmesh::mesh::export_obj(surface_mesh, out_path, export_opts, context);
    if (export_status == sqmesh::base::StatusCode::ok) {
        std::cout << "[OK] Surface mesh exported to " << out_path << "\n";
    } else {
        std::cerr << "[WARN] OBJ export failed. status="
                  << static_cast<int>(export_status) << "\n";
    }

    return 0;
}
