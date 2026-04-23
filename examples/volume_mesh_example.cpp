#include "sqmesh/sqmesh.hpp"
#include "../src/mesh/auto_cfd/auto_cfd_surface_test_hook.hpp"
#include "../src/core/log.hpp"
#include "../src/core/runtime_registry.hpp"
#include "../src/mesh/region/region_detector.hpp"

#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

std::uint64_t pack_entity_ref(sqmesh::mesh::EntityRef ref) {
    return (static_cast<std::uint64_t>(ref.entity_group) << 32U) |
           static_cast<std::uint64_t>(ref.index);
}

std::array<double, 3> parse_point3(const std::string &text,
                                   std::array<double, 3> fallback) {
    std::array<double, 3> out = fallback;
    std::size_t start = 0, comma_idx = 0;
    for (int i = 0; i < 3 && start <= text.size(); ++i) {
        comma_idx = text.find(',', start);
        const std::string token = text.substr(
            start,
            comma_idx == std::string::npos ? std::string::npos : comma_idx - start);
        if (!token.empty()) {
            try { out[i] = std::stod(token); } catch (...) { out[i] = fallback[i]; }
        }
        if (comma_idx == std::string::npos) break;
        start = comma_idx + 1;
    }
    return out;
}

bool ray_intersects_triangle(
    const std::array<double, 3> &origin,
    const std::array<double, 3> &dir,
    const std::array<double, 3> &v0,
    const std::array<double, 3> &v1,
    const std::array<double, 3> &v2) {
    const std::array<double, 3> e1{v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]};
    const std::array<double, 3> e2{v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2]};
    const std::array<double, 3> h{
        dir[1]*e2[2]-dir[2]*e2[1],
        dir[2]*e2[0]-dir[0]*e2[2],
        dir[0]*e2[1]-dir[1]*e2[0]};
    const double a = e1[0]*h[0]+e1[1]*h[1]+e1[2]*h[2];
    if (std::abs(a) < 1e-20) return false;
    const double f = 1.0 / a;
    const std::array<double, 3> s{origin[0]-v0[0], origin[1]-v0[1], origin[2]-v0[2]};
    const double u = f * (s[0]*h[0]+s[1]*h[1]+s[2]*h[2]);
    if (u < 0.0 || u > 1.0) return false;
    const std::array<double, 3> q{
        s[1]*e1[2]-s[2]*e1[1],
        s[2]*e1[0]-s[0]*e1[2],
        s[0]*e1[1]-s[1]*e1[0]};
    const double v = f * (dir[0]*q[0]+dir[1]*q[1]+dir[2]*q[2]);
    if (v < 0.0 || u + v > 1.0) return false;
    const double t = f * (e2[0]*q[0]+e2[1]*q[1]+e2[2]*q[2]);
    return t > 1e-10;
}

struct VolumeRegionSelection {
    int target_region_id = -1;
    std::set<std::uint32_t> target_region_face_ids;   // farfield side
    std::set<std::uint32_t> boundary_layer_face_ids;  // object side
};

VolumeRegionSelection select_volume_regions(
    const sqmesh::mesh::Domain &domain,
    const sqmesh::mesh::RegionDetectionResult &region_result,
    const std::array<double, 3> &material_point) {
    VolumeRegionSelection out;
    if (region_result.regions.empty()) return out;

    const std::array<double, 3> ray_dir{1.0, 0.3713906763, 0.2189474511};
    auto region_contains_point = [&](const sqmesh::mesh::MeshRegion &r) {
        std::size_t hits = 0;
        for (const auto boundary_face_id : r.boundary_face_ids) {
            if (boundary_face_id >= region_result.face_refs.size()) continue;
            const auto face_ref = region_result.face_refs[boundary_face_id];
            const auto fn = domain.face_nodes(face_ref);
            if (fn.size < 3) continue;
            std::array<std::array<double, 3>, 4> vs{};
            for (std::size_t i = 0; i < fn.size && i < 4; ++i) {
                const auto &nd = domain.node(fn[i]);
                vs[i] = {nd.coordinates[0], nd.coordinates[1], nd.coordinates[2]};
            }
            if (ray_intersects_triangle(material_point, ray_dir, vs[0], vs[1], vs[2])) ++hits;
            if (fn.size == 4 &&
                ray_intersects_triangle(material_point, ray_dir, vs[0], vs[2], vs[3])) ++hits;
        }
        return (hits % 2U) == 1U;
    };

    for (const auto &r : region_result.regions) {
        if (region_contains_point(r)) { out.target_region_id = r.id; break; }
    }
    if (out.target_region_id < 0) {
        double best = std::numeric_limits<double>::infinity();
        for (const auto &r : region_result.regions) {
            const double dx = r.interior_point[0] - material_point[0];
            const double dy = r.interior_point[1] - material_point[1];
            const double dz = r.interior_point[2] - material_point[2];
            const double d2 = dx*dx + dy*dy + dz*dz;
            if (d2 < best) { best = d2; out.target_region_id = r.id; }
        }
    }

    auto collect = [&](const sqmesh::mesh::MeshRegion &r,
                       std::set<std::uint32_t> &sink) {
        for (const auto boundary_face_id : r.boundary_face_ids) {
            if (boundary_face_id >= region_result.face_refs.size()) continue;
            const auto face_ref = region_result.face_refs[boundary_face_id];
            auto topo = domain.face_topology_owner(face_ref);
            if (!sqmesh::geo::is_valid(topo)) continue;
            sink.insert(topo.index);
        }
    };
    for (const auto &r : region_result.regions) {
        if (r.id == out.target_region_id) collect(r, out.target_region_face_ids);
        else collect(r, out.boundary_layer_face_ids);
    }
    return out;
}

std::string join_topology_ids(const std::set<std::uint32_t> &ids) {
    std::string s;
    bool first = true;
    for (auto id : ids) {
        if (!first) s.push_back(',');
        s.append(std::to_string(id));
        first = false;
    }
    return s;
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0]
                  << " <model.step|.iges> [min_len] [max_len] [dist_angle]"
                  << " [growth] [proximity:0|1] [extra_topo:0|1]"
                  << " [bl_first_height] [bl_growth_rate] [bl_num_layers]"
                  << " [tet_max_length] [tet_growth_rate] [material_point_x,y,z]"
                  << " [bl_first_height_mode:absolute|aspect] [bl_first_height_aspect]\n";
        std::cerr << "  Surface defaults: min=1.0, max=20.0, dist=15.0, growth=1.2, proximity=0, extra_topo=1\n";
        std::cerr << "  Volume defaults:  bl_first_height=min_len, bl_growth_rate=1.2, bl_num_layers=5,\n";
        std::cerr << "                    tet_max_length=max_len, tet_growth_rate=growth_rate,\n";
        std::cerr << "                    material_point=-500,0,0 (farfield seed for multi-region geometries),\n";
        std::cerr << "                    bl_first_height_mode=absolute, bl_first_height_aspect=0.5\n";
        std::cerr << "  BL face IDs are auto-classified from the material point and detected regions.\n";
        return 1;
    }

    const char* file_path = argv[1];
    const std::string log_path = std::filesystem::path(file_path).stem().string() + ".log";
    sqmesh::log::init(log_path, "info");

    const double min_len = argc > 2 ? std::atof(argv[2]) : 1;
    const double max_len = argc > 3 ? std::atof(argv[3]) : 20;
    const double distortion_angle = argc > 4 ? std::atof(argv[4]) : 15.0;
    const double growth_rate = argc > 5 ? std::atof(argv[5]) : 1.2;
    const bool   proximity = argc > 6 ? (std::atoi(argv[6]) != 0) : true;
    const bool   auto_topo = argc > 7 ? (std::atoi(argv[7]) != 0) : true;
    const double bl_first_height = argc > 8 ? std::atof(argv[8]) : min_len;
    const double bl_growth_rate = argc > 9 ? std::atof(argv[9]) : 1.2;
    const int    bl_num_layers = argc > 10 ? std::atoi(argv[10]) : 5;
    const double tet_max_length = argc > 11 ? std::atof(argv[11]) : max_len;
    const double tet_growth_rate = argc > 12 ? std::atof(argv[12]) : growth_rate;
    const std::string material_point_text = argc > 13 ? argv[13] : "-500,0,0";
    const std::string bl_first_height_mode = argc > 14 ? argv[14] : "absolute";
    const double bl_first_height_aspect = argc > 15 ? std::atof(argv[15]) : 0.5;
    const std::array<double, 3> material_point =
        parse_point3(material_point_text, {-500.0, 0.0, 0.0});

    // ── 1. Initialize ──────────────────────────────────────────────────
    sqmesh::base::ContextHandle context;
    if (sqmesh::base::initialize(context) != sqmesh::base::StatusCode::ok) {
        std::cerr << "[ERROR] SQMesh Initialization failed.\n";
        return 1;
    }

    // ── 2. Import geometry (STEP or IGES auto-detected by extension) ──
    sqmesh::geo::ModelHandle model;

    auto ends_with_ci = [](const std::string& s, const std::string& suffix) {
        if (s.size() < suffix.size()) return false;
        for (std::size_t i = 0; i < suffix.size(); ++i) {
            const char a = s[s.size() - suffix.size() + i];
            const char b = suffix[i];
            const char al = (a >= 'A' && a <= 'Z') ? (a + 32) : a;
            const char bl = (b >= 'A' && b <= 'Z') ? (b + 32) : b;
            if (al != bl) return false;
        }
        return true;
    };

    const std::string path_str(file_path);
    const bool is_iges = ends_with_ci(path_str, ".iges") || ends_with_ci(path_str, ".igs");

    auto import_geometry = [&](const char* import_path,
                               sqmesh::geo::ModelHandle& target_model,
                               bool verbose) -> sqmesh::base::StatusCode {
        if (verbose) {
            std::cout << "Loading " << (is_iges ? "IGES" : "STEP")
                      << " model: " << import_path << " ...\n";
        }

        return is_iges
            ? sqmesh::geo::import_iges(import_path, target_model, {}, context)
            : sqmesh::geo::import_step(import_path, target_model, {}, context);
    };

    auto import_status = import_geometry(file_path, model, true);
    if (import_status != sqmesh::base::StatusCode::ok) {
        std::cerr << "[ERROR] Failed to load "
                  << (is_iges ? "IGES" : "STEP")
                  << " model. status="
                  << static_cast<int>(import_status) << "\n";
        return 1;
    }

    auto print_model_summary = [&](const char* label) {
        sqmesh::geo::ModelSummary summary;
        if (sqmesh::geo::model_summary(model, summary, context) !=
            sqmesh::base::StatusCode::ok) {
            return;
        }
        std::cout << "=====================================\n";
        if (label && label[0] != '\0') {
            std::cout << " [Geom] " << label << "\n";
        }
        std::cout << " [Geom] Faces (surfaces): " << summary.face_count << "\n";
        std::cout << " [Geom] Edges (curves):   " << summary.edge_count << "\n";
        std::cout << " [Geom] Shells:           " << summary.shell_count << "\n";
        std::cout << " [Geom] Solids:           " << summary.solid_count << "\n";
        std::cout << "=====================================\n";
    };

    auto print_topology_check = [&](const char* label) {
        sqmesh::geo::TopologyCheckReport report;
        if (sqmesh::geo::check_topology(model, report, context) !=
            sqmesh::base::StatusCode::ok) {
            return;
        }
        std::cout << "\n── " << label << " ─────────────────────\n";
        std::cout << "  Valid:            " << (report.is_valid ? "YES" : "NO") << "\n";
        std::cout << "  Free edges:       " << report.free_edge_count << "\n";
        std::cout << "  Contiguous edges: " << report.contiguous_edge_count << "\n";
        std::cout << "  Multiple edges:   " << report.multiple_edge_count << "\n";
        std::cout << "────────────────────────────────────────\n";
    };

    print_model_summary("Imported");
    print_topology_check("Topology After Import");

    if (auto_topo) {
        sqmesh::geo::TopoReport topo_report;
        sqmesh::geo::TopoOptions topo_options;
        const auto topo_status = sqmesh::geo::topo(
            model, topo_report, topo_options, context);
        if (topo_status == sqmesh::base::StatusCode::ok) {
            std::cout << "\n── Extra Topo Repair ──────────────────\n";
            std::cout << "  Modified:                " << (topo_report.modified ? "YES" : "NO") << "\n";
            std::cout << "  Topology ID changed:     "
                << (topo_report.topology_identity_changed ? "YES" : "NO") << "\n";
            std::cout << "  Sewing performed:        " << (topo_report.sewing_performed ? "YES" : "NO") << "\n";
            std::cout << "  Sewing modified:         " << (topo_report.sewing_modified ? "YES" : "NO") << "\n";
            std::cout << "  Shape-fix performed:     " << (topo_report.shape_fix_performed ? "YES" : "NO") << "\n";
            std::cout << "  Shape-fix modified:      " << (topo_report.shape_fix_modified ? "YES" : "NO") << "\n";
            std::cout << "  Free edges reduced:      " << (topo_report.free_edges_reduced ? "YES" : "NO") << "\n";
            std::cout << "  Free edges:              " << topo_report.before.free_edge_count
                << " -> " << topo_report.after.free_edge_count << "\n";
            std::cout << "  Multiple edges:          " << topo_report.before.multiple_edge_count
                << " -> " << topo_report.after.multiple_edge_count << "\n";
            std::cout << "  Faces:                   " << topo_report.before.model_summary.face_count
                << " -> " << topo_report.after.model_summary.face_count << "\n";
            std::cout << "  Edges:                   " << topo_report.before.model_summary.edge_count
                << " -> " << topo_report.after.model_summary.edge_count << "\n";
            std::cout << "────────────────────────────────────────\n";
            print_model_summary("After Auto Topo");
            print_topology_check("Topology After Auto Topo");
        }
        else {
            std::cerr << "[WARN] Auto topo repair failed.\n";
        }
    }

    // ── 3. Run Auto CFD Surface Mesher ──────────────────────────────────
    std::cout << "\n======== Auto CFD Surface Mesher ========\n";
    std::cout << " min_length       = " << min_len << "\n";
    std::cout << " max_length       = " << max_len << "\n";
    std::cout << " distortion_angle = " << distortion_angle << "\n";
    std::cout << " growth_rate      = " << growth_rate << "\n";
    std::cout << " proximity        = " << (proximity ? "ON" : "OFF") << "\n";
    std::cout << " validation_mode  = diagnostic"
              << " (allow_quality_gate_failure=ON, allow_final_screen_failure=ON)\n";
    std::cout << "=========================================\n";

    sqmesh::mesh::MeshingOptions meshing_options;
    meshing_options.parameters.set_number("minimum_length", min_len);
    meshing_options.parameters.set_number("maximum_length", max_len);
    meshing_options.parameters.set_number("distortion_angle", distortion_angle);
    meshing_options.parameters.set_number("growth_rate", growth_rate);
    meshing_options.parameters.set_boolean("proximity", proximity);
    meshing_options.parameters.set_boolean("allow_quality_gate_failure", true);
    meshing_options.parameters.set_boolean("allow_final_screen_failure", true);

    sqmesh::mesh::MeshHandle cfd_mesh;

    auto t_start = std::chrono::high_resolution_clock::now();

    auto mesh_status = sqmesh::mesh::create_surface_mesh(
        model, "Auto CFD Surface Mesher", meshing_options, cfd_mesh, context);

    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed_ms =
        std::chrono::duration<double, std::milli>(t_end - t_start).count();

    if (mesh_status != sqmesh::base::StatusCode::ok) {
        std::cerr << "[ERROR] CFD Surface Meshing FAILED (StatusCode = "
            << static_cast<int>(mesh_status) << ")\n";
        std::cerr << "  Status name: "
            << sqmesh::base::status_code_name(mesh_status) << "\n";
        auto err_msg = sqmesh::base::last_error_message();
        if (!err_msg.empty()) {
            std::cerr << "  Error message: " << err_msg << "\n";
        }
        std::cerr << "  Elapsed: " << elapsed_ms << " ms\n";
        return 1;
    }
    std::cout << "\n[OK] CFD Surface Meshing succeeded!\n";
    std::cout << "  Elapsed: " << elapsed_ms << " ms\n";

    // ── 5. Print mesh summary ──────────────────────────────────────────
    sqmesh::mesh::MeshSummary mesh_summary;
    if (sqmesh::mesh::mesh_summary(cfd_mesh, mesh_summary, context) ==
        sqmesh::base::StatusCode::ok) {
        std::cout << "\n── Mesh Summary ────────────────────────\n";
        std::cout << "  Nodes: " << mesh_summary.node_count << "\n";
        std::cout << "  Edges: " << mesh_summary.edge_count << "\n";
        std::cout << "  Faces: " << mesh_summary.face_count << "\n";
        std::cout << "────────────────────────────────────────\n";
    }

    // ── 6. Quality report ──────────────────────────────────────────────
    sqmesh::mesh::MeshQualityReport quality_report;
    if (sqmesh::mesh::mesh_quality_report(cfd_mesh, quality_report, context) ==
        sqmesh::base::StatusCode::ok) {
        std::cout << "\n── Quality Report ──────────────────────\n";
        std::cout << "  Supported elements:   " << quality_report.supported_element_count << "\n";
        std::cout << "  Valid elements:       " << quality_report.valid_element_count << "\n";
        std::cout << "  Degenerate elements:  " << quality_report.degenerate_element_count << "\n";
        std::cout << "  Inverted elements:    " << quality_report.inverted_element_count << "\n";

        for (const auto& kind_summary : quality_report.kinds) {
            if (kind_summary.kind != sqmesh::mesh::EntityKind::face_triangle)
                continue;
            std::cout << "\n  [Triangle Quality]\n";
            std::cout << "    Count:        " << kind_summary.supported_element_count << "\n";
            std::cout << "    Min angle:    " << kind_summary.min_angle.minimum
                << " ~ " << kind_summary.min_angle.maximum
                << "  (avg " << kind_summary.min_angle.average << ")\n";
            std::cout << "    Max angle:    " << kind_summary.max_angle.minimum
                << " ~ " << kind_summary.max_angle.maximum
                << "  (avg " << kind_summary.max_angle.average << ")\n";
            std::cout << "    Aspect ratio: " << kind_summary.aspect_ratio.minimum
                << " ~ " << kind_summary.aspect_ratio.maximum
                << "  (avg " << kind_summary.aspect_ratio.average << ")\n";
            std::cout << "    Skewness:     " << kind_summary.skewness.minimum
                << " ~ " << kind_summary.skewness.maximum
                << "  (avg " << kind_summary.skewness.average << ")\n";
            std::cout << "    Radius ratio: " << kind_summary.radius_ratio.minimum
                << " ~ " << kind_summary.radius_ratio.maximum
                << "  (avg " << kind_summary.radius_ratio.average << ")\n";
        }
        std::cout << "────────────────────────────────────────\n";
    }

    // ── 7. Export surface mesh to OBJ (before volume meshing) ────────
    {
      const std::string surface_path = "cfd_surface_mesh.obj";
      sqmesh::mesh::ObjExportOptions surface_export_options;
      auto surface_export_status =
        sqmesh::mesh::export_obj(cfd_mesh, surface_path, surface_export_options, context);
      if(surface_export_status == sqmesh::base::StatusCode::ok) {
        std::cout << "\n[OK] CFD surface mesh exported to " << surface_path << "\n";
      }
    }

    // ── 7b. Print face topology owners ─────────────────────────────
    (void)sqmesh::core::detail::with_mesh_domain(
      cfd_mesh, [](const sqmesh::mesh::Domain &domain) {
        struct FaceStats {
          std::uint32_t count = 0;
          double sx = 0.0, sy = 0.0, sz = 0.0;   // centroid accumulator
          double minx = 1e30, miny = 1e30, minz = 1e30;
          double maxx = -1e30, maxy = -1e30, maxz = -1e30;
        };
        std::map<std::uint32_t, FaceStats> topo_face_stats;
        for(const auto &entity_group : domain.entity_groups()) {
          if(entity_group.order() != sqmesh::mesh::EntityOrder::face) continue;
          if(entity_group.role() != sqmesh::mesh::EntityGroupRole::computational) continue;
          for(std::uint32_t fi = 0; fi < entity_group.faces().size(); ++fi) {
            sqmesh::mesh::EntityRef fref{entity_group.id(), fi};
            auto topo = domain.face_topology_owner(fref);
            auto &s = topo_face_stats[topo.index];
            auto nodes = domain.face_nodes(fref);
            for(std::size_t ni = 0; ni < nodes.size; ++ni) {
              const auto &p = domain.node(nodes[ni]).coordinates;
              s.sx += p[0]; s.sy += p[1]; s.sz += p[2];
              s.minx = std::min(s.minx, p[0]); s.maxx = std::max(s.maxx, p[0]);
              s.miny = std::min(s.miny, p[1]); s.maxy = std::max(s.maxy, p[1]);
              s.minz = std::min(s.minz, p[2]); s.maxz = std::max(s.maxz, p[2]);
            }
            s.count += 1;
          }
        }
        std::cout << "\n[Topology] Face topology owners (0-indexed):\n";
        std::cout << std::fixed << std::setprecision(1);
        for(auto &[idx, s] : topo_face_stats) {
          const double n_nodes = static_cast<double>(s.count) * 3.0;
          std::cout << "  Face " << idx << ": " << s.count << " tri, "
                    << "centroid=(" << s.sx/n_nodes << ", " << s.sy/n_nodes
                    << ", " << s.sz/n_nodes << "), "
                    << "x[" << s.minx << "," << s.maxx << "] "
                    << "y[" << s.miny << "," << s.maxy << "] "
                    << "z[" << s.minz << "," << s.maxz << "]\n";
        }
        std::cout.unsetf(std::ios_base::floatfield);
        return sqmesh::base::StatusCode::ok;
      }, context
    );

    // ── 7c. Detect regions + auto-classify BL / target faces ────────
    sqmesh::mesh::RegionDetectionResult region_result;
    (void)sqmesh::core::detail::with_mesh_domain(
      cfd_mesh, [&region_result](const sqmesh::mesh::Domain &domain) {
        region_result = sqmesh::mesh::detect_regions(domain);
        std::cout << "\n[Regions] Detected " << region_result.regions.size() << " enclosed regions:\n";
        for(const auto &r : region_result.regions) {
          std::cout << "  Region " << r.id << ": " << r.boundary_face_ids.size()
                    << " boundary faces, interior=("
                    << r.interior_point[0] << ", "
                    << r.interior_point[1] << ", "
                    << r.interior_point[2] << ")\n";
        }
        return sqmesh::base::StatusCode::ok;
      }, context
    );

    VolumeRegionSelection volume_selection;
    (void)sqmesh::core::detail::with_mesh_domain(
      cfd_mesh,
      [&](const sqmesh::mesh::Domain &domain) {
        volume_selection = select_volume_regions(domain, region_result, material_point);
        return sqmesh::base::StatusCode::ok;
      }, context
    );
    const std::string bl_face_ids_text = join_topology_ids(volume_selection.boundary_layer_face_ids);
    std::cout << "\n[Volume] Material point=(" << material_point[0]
              << ", " << material_point[1] << ", " << material_point[2]
              << ") -> target region=" << volume_selection.target_region_id << "\n";
    std::cout << "  Target (farfield) face IDs: "
              << join_topology_ids(volume_selection.target_region_face_ids) << "\n";
    std::cout << "  Boundary layer face IDs:    " << bl_face_ids_text << "\n";

    // ── 8. Boundary Layer mesh ──────────────────────────────────────
    sqmesh::mesh::MeshHandle bl_mesh = sqmesh::invalid_handle;
    {
      std::cout << "\n======== Boundary Layer Mesher ========\n";
      sqmesh::mesh::ParameterDictionary bl_params;
      bl_params.set_number("bl_num_layers", bl_num_layers);
      bl_params.set_number("bl_first_height", bl_first_height);
      bl_params.set_number("bl_growth_rate", bl_growth_rate);
      bl_params.set_text("bl_first_height_mode", bl_first_height_mode);
      bl_params.set_number("bl_first_height_aspect", bl_first_height_aspect);
      bl_params.set_number("bl_max_skewness", 0.98);
      bl_params.set_number("bl_smooth_iterations", 5);
      bl_params.set_number("bl_max_normal_deviation", 40.0);
      bl_params.set_boolean("bl_smooth_normals", true);
      bl_params.set_boolean("bl_allow_proximity", true);
      bl_params.set_number("bl_proximity_factor", 0.5);
      bl_params.set_text("bl_treatment", "collapse");
      bl_params.set_boolean("bl_allow_squeeze", false);
      if (!bl_face_ids_text.empty()) {
        bl_params.set_text("bl_face_ids", bl_face_ids_text);
      }
      bl_params.set_text("bl_target_point", material_point_text);

      auto bl_t0 = std::chrono::high_resolution_clock::now();
      auto bl_status = sqmesh::mesh::create_volume_mesh(
        model, "Boundary Layer Mesher", bl_params, bl_mesh, context
      );
      auto bl_t1 = std::chrono::high_resolution_clock::now();
      double bl_elapsed = std::chrono::duration<double, std::milli>(bl_t1 - bl_t0).count();

      if(bl_status != sqmesh::base::StatusCode::ok) {
        std::cerr << "[WARN] Boundary layer meshing failed (status="
                  << static_cast<int>(bl_status) << "). Tet will fill the full volume.\n";
        auto err_msg = sqmesh::base::last_error_message();
        if(!err_msg.empty()) {
          std::cerr << "  Error: " << err_msg << "\n";
        }
      }
      else {
        std::cout << "[OK] Boundary layer meshing succeeded!\n";
        std::cout << "  Elapsed: " << bl_elapsed << " ms\n";

        sqmesh::mesh::MeshSummary bl_summary;
        if(sqmesh::mesh::mesh_summary(bl_mesh, bl_summary, context) ==
           sqmesh::base::StatusCode::ok) {
          std::cout << "  Nodes: " << bl_summary.node_count << "\n";
          std::cout << "  Faces: " << bl_summary.face_count << "\n";
          std::cout << "  Cells: " << bl_summary.cell_count << "\n";
        }
      }
      std::cout << "========================================\n";
    }

    // ── 9. Tetrahedral volume mesh ───────────────────────────────────
    sqmesh::mesh::MeshHandle tet_mesh = sqmesh::invalid_handle;
    {
      std::cout << "\n======== Tetrahedral Volume Mesher ========\n";
      sqmesh::mesh::ParameterDictionary tet_params;
      tet_params.set_number("minimum_length", min_len);
      tet_params.set_number("maximum_length", tet_max_length);
      tet_params.set_number("growth_rate", tet_growth_rate);
      tet_params.set_number("quality_ratio", 1.414);
      tet_params.set_text("material_points", material_point_text);
      auto tet_t0 = std::chrono::high_resolution_clock::now();
      auto tet_status = sqmesh::mesh::create_volume_mesh(
        model, "Tetrahedral Volume Mesher", tet_params, tet_mesh, context
      );
      auto tet_t1 = std::chrono::high_resolution_clock::now();
      double tet_elapsed = std::chrono::duration<double, std::milli>(tet_t1 - tet_t0).count();

      if(tet_status != sqmesh::base::StatusCode::ok) {
        std::cerr << "[WARN] Tetrahedral volume meshing failed (status="
                  << static_cast<int>(tet_status) << "). Skipping.\n";
      }
      else {
        std::cout << "[OK] Tetrahedral volume meshing succeeded!\n";
        std::cout << "  Elapsed: " << tet_elapsed << " ms\n";

        sqmesh::mesh::MeshSummary tet_summary;
        if(sqmesh::mesh::mesh_summary(tet_mesh, tet_summary, context) ==
           sqmesh::base::StatusCode::ok) {
          std::cout << "  Nodes: " << tet_summary.node_count << "\n";
          std::cout << "  Faces: " << tet_summary.face_count << "\n";
          std::cout << "  Cells: " << tet_summary.cell_count << "\n";
        }

        // Volume quality report (per-kind: tetra / prism / hexa / pyramid).
        sqmesh::mesh::MeshQualityReport tet_quality;
        if (sqmesh::mesh::mesh_quality_report(tet_mesh, tet_quality, context) ==
            sqmesh::base::StatusCode::ok) {
          std::cout << "\n── Volume Quality Report ───────────────\n";
          std::cout << "  Supported elements:   " << tet_quality.supported_element_count << "\n";
          std::cout << "  Valid elements:       " << tet_quality.valid_element_count << "\n";
          std::cout << "  Degenerate elements:  " << tet_quality.degenerate_element_count << "\n";
          std::cout << "  Inverted elements:    " << tet_quality.inverted_element_count << "\n";
          for (const auto &ks : tet_quality.kinds) {
            const char *kind_name = "";
            switch (ks.kind) {
              case sqmesh::mesh::EntityKind::cell_tetra:   kind_name = "tetra";   break;
              case sqmesh::mesh::EntityKind::cell_prism:   kind_name = "prism";   break;
              case sqmesh::mesh::EntityKind::cell_hexa:    kind_name = "hexa";    break;
              case sqmesh::mesh::EntityKind::cell_pyramid: kind_name = "pyramid"; break;
              default: continue;
            }
            std::cout << "\n  [" << kind_name << "]  count=" << ks.supported_element_count
                      << "  valid=" << ks.valid_element_count
                      << "  degen=" << ks.degenerate_element_count
                      << "  inverted=" << ks.inverted_element_count << "\n";
            std::cout << "    Min dihedral: "   << ks.min_angle.minimum
                      << " ~ " << ks.min_angle.maximum
                      << " (avg " << ks.min_angle.average << ")\n";
            std::cout << "    Max dihedral: "   << ks.max_angle.minimum
                      << " ~ " << ks.max_angle.maximum
                      << " (avg " << ks.max_angle.average << ")\n";
            std::cout << "    Aspect ratio: "  << ks.aspect_ratio.minimum
                      << " ~ " << ks.aspect_ratio.maximum
                      << " (avg " << ks.aspect_ratio.average << ")\n";
            std::cout << "    Skewness:     "  << ks.skewness.minimum
                      << " ~ " << ks.skewness.maximum
                      << " (avg " << ks.skewness.average << ")\n";
            std::cout << "    Radius ratio: "  << ks.radius_ratio.minimum
                      << " ~ " << ks.radius_ratio.maximum
                      << " (avg " << ks.radius_ratio.average << ")\n";
          }
          std::cout << "────────────────────────────────────────\n";
        }

        // ── 9. Export volume mesh (BL prisms + tets) to VTK ──────────
        sqmesh::mesh::Domain tet_domain;
        if (sqmesh::mesh::domain_snapshot(tet_mesh, tet_domain, context) == sqmesh::base::StatusCode::ok) {
            const std::string vtk_path = "cfd_volume_mesh.vtk";
            std::ofstream out(vtk_path, std::ios::out | std::ios::trunc);
            if (out) {
                out << "# vtk DataFile Version 3.0\n";
                out << "SQMesh Volume Mesh\n";
                out << "ASCII\n";
                out << "DATASET UNSTRUCTURED_GRID\n";

                std::size_t n_count = 0;
                for (const auto& entity_group : tet_domain.entity_groups()) {
                    if (entity_group.role() == sqmesh::mesh::EntityGroupRole::computational &&
                        entity_group.order() == sqmesh::mesh::EntityOrder::node) {
                        n_count += entity_group.nodes().size();
                    }
                }
                out << "POINTS " << n_count << " double\n";

                std::unordered_map<std::uint64_t, std::size_t> vtk_node_indices;
                vtk_node_indices.reserve(n_count);
                std::size_t next_v_idx = 0;
                for (const auto& entity_group : tet_domain.entity_groups()) {
                    if (entity_group.role() == sqmesh::mesh::EntityGroupRole::computational &&
                        entity_group.order() == sqmesh::mesh::EntityOrder::node) {
                        for (std::uint32_t i = 0; i < entity_group.nodes().size(); ++i) {
                            sqmesh::mesh::EntityRef ref{entity_group.id(), i};
                            const auto& pos = tet_domain.node(ref).coordinates;
                            out << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
                            vtk_node_indices[pack_entity_ref(ref)] = next_v_idx++;
                        }
                    }
                }

                // Collect all volume cells (tet=10, wedge=13, hex=12).
                struct VtkCell { int type; std::vector<std::size_t> nodes; };
                std::vector<VtkCell> cells;
                std::size_t tet_count = 0, prism_count = 0, hex_count = 0;
                for (const auto& entity_group : tet_domain.entity_groups()) {
                    if (entity_group.role() != sqmesh::mesh::EntityGroupRole::computational ||
                        entity_group.order() != sqmesh::mesh::EntityOrder::cell) {
                        continue;
                    }
                    for (std::uint32_t i = 0; i < entity_group.cells().size(); ++i) {
                        sqmesh::mesh::EntityRef ref{entity_group.id(), i};
                        auto nodes = tet_domain.cell_nodes(ref);
                        int type = 0;
                        if (nodes.size == 4) { type = 10; ++tet_count; }
                        else if (nodes.size == 6) { type = 13; ++prism_count; }
                        else if (nodes.size == 8) { type = 12; ++hex_count; }
                        else { continue; }
                        VtkCell c;
                        c.type = type;
                        c.nodes.reserve(nodes.size);
                        for (std::size_t ni = 0; ni < nodes.size; ++ni) {
                            c.nodes.push_back(vtk_node_indices[pack_entity_ref(nodes[ni])]);
                        }
                        cells.push_back(std::move(c));
                    }
                }

                std::size_t total_conn = 0;
                for (const auto& c : cells) total_conn += c.nodes.size() + 1;

                out << "\nCELLS " << cells.size() << " " << total_conn << "\n";
                for (const auto& c : cells) {
                    out << c.nodes.size();
                    for (auto n : c.nodes) out << " " << n;
                    out << "\n";
                }

                out << "\nCELL_TYPES " << cells.size() << "\n";
                for (const auto& c : cells) out << c.type << "\n";

                out.close();
                std::cout << "[OK] Volume mesh exported to " << vtk_path
                          << " (" << tet_count << " tets, "
                          << prism_count << " prisms, "
                          << hex_count << " hexas)\n";

                // Also export only the tet cells to a separate VTK so the
                // tet-only region can be inspected in isolation from the BL
                // prism shell.  Points are re-indexed so the file is
                // standalone (no need to reference the full mesh points).
                const std::string tet_vtk_path = "cfd_tet_only.vtk";
                std::ofstream out_tet(tet_vtk_path, std::ios::out | std::ios::trunc);
                if (out_tet) {
                    std::unordered_map<std::size_t, std::size_t> tet_idx_remap;
                    tet_idx_remap.reserve(tet_count * 4);
                    std::vector<std::size_t> tet_pts;
                    tet_pts.reserve(tet_count * 4);
                    std::vector<std::array<std::size_t, 4>> tets_only;
                    tets_only.reserve(tet_count);
                    for (const auto &c : cells) {
                        if (c.type != 10) continue;
                        std::array<std::size_t, 4> t{};
                        for (int k = 0; k < 4; ++k) {
                            auto it = tet_idx_remap.find(c.nodes[k]);
                            if (it == tet_idx_remap.end()) {
                                const std::size_t new_idx = tet_pts.size();
                                tet_idx_remap.emplace(c.nodes[k], new_idx);
                                tet_pts.push_back(c.nodes[k]);
                                t[k] = new_idx;
                            } else {
                                t[k] = it->second;
                            }
                        }
                        tets_only.push_back(t);
                    }

                    // Rebuild coord table for the subset of points used.
                    std::vector<std::array<double, 3>> all_coords;
                    all_coords.reserve(n_count);
                    for (const auto &entity_group : tet_domain.entity_groups()) {
                        if (entity_group.role() == sqmesh::mesh::EntityGroupRole::computational &&
                            entity_group.order() == sqmesh::mesh::EntityOrder::node) {
                            for (std::uint32_t i = 0; i < entity_group.nodes().size(); ++i) {
                                sqmesh::mesh::EntityRef ref{entity_group.id(), i};
                                const auto &pos = tet_domain.node(ref).coordinates;
                                all_coords.push_back({pos[0], pos[1], pos[2]});
                            }
                        }
                    }

                    out_tet << "# vtk DataFile Version 3.0\n";
                    out_tet << "SQMesh Tet-Only Mesh\n";
                    out_tet << "ASCII\n";
                    out_tet << "DATASET UNSTRUCTURED_GRID\n";
                    out_tet << "POINTS " << tet_pts.size() << " double\n";
                    for (std::size_t global_idx : tet_pts) {
                        const auto &p = all_coords[global_idx];
                        out_tet << p[0] << " " << p[1] << " " << p[2] << "\n";
                    }
                    out_tet << "\nCELLS " << tets_only.size() << " "
                            << (tets_only.size() * 5) << "\n";
                    for (const auto &t : tets_only) {
                        out_tet << "4 " << t[0] << " " << t[1] << " " << t[2]
                                << " " << t[3] << "\n";
                    }
                    out_tet << "\nCELL_TYPES " << tets_only.size() << "\n";
                    for (std::size_t i = 0; i < tets_only.size(); ++i) {
                        out_tet << "10\n";
                    }
                    out_tet.close();
                    std::cout << "[OK] Tet-only mesh exported to " << tet_vtk_path
                              << " (" << tet_pts.size() << " nodes, "
                              << tets_only.size() << " tets)\n";
                }
            }
        }
      }
      std::cout << "===========================================\n";
    }

    return 0;
}
