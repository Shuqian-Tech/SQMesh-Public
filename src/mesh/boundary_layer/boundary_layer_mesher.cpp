// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "boundary_layer_mesher.hpp"
#include "boundary_layer_solver.hpp"
#include "../framework/meshing_framework.hpp"
#include "../region/region_detector.hpp"

#include "core/log.hpp"
#include "core/runtime_registry.hpp"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdarg>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>

namespace sqmesh::mesh::detail {

namespace {

constexpr std::string_view kBLMesherName = "Boundary Layer Mesher";
constexpr bool kEnableBoundaryLayerDebugFileOutput = false;

void bl_log_vprint(const char *level, const char *fmt, std::va_list args)
{
  std::fprintf(stderr, "[%s] ", level);
  std::vfprintf(stderr, fmt, args);
  std::fputc('\n', stderr);
}

void bl_log_info(const char *fmt, ...)
{
  std::va_list args;
  va_start(args, fmt);
  bl_log_vprint("info", fmt, args);
  va_end(args);
}

void bl_log_warn(const char *fmt, ...)
{
  std::va_list args;
  va_start(args, fmt);
  bl_log_vprint("warn", fmt, args);
  va_end(args);
}

// ===========================================================================
// BoundaryLayerAlgorithm: MeshingAlgorithm implementation
// ===========================================================================
class BoundaryLayerAlgorithm final : public MeshingAlgorithm
{
public:
  [[nodiscard]] std::string_view name() const noexcept override
  {
    return kBLMesherName;
  }

  [[nodiscard]] std::string_view entity_group_prefix() const noexcept override
  {
    return "bl_";
  }

protected:
  void reset() noexcept override
  {
    request_ = {};
    params_ = {};
  }

  [[nodiscard]] base::StatusCode on_initialize(
    const MeshingRequest &request
  ) override
  {
    if(request.target_dimension != MeshingDimension::volume) {
      return core::detail::publish_error(
        base::StatusCode::unsupported,
        "Boundary Layer Mesher only supports volume meshing requests."
      );
    }
    if(false) {
      // Session/model handles are optional (not needed for OBJ import pipeline)
      return core::detail::publish_error(
        base::StatusCode::invalid_handle,
        "Boundary Layer Mesher requires valid context handle."
      );
    }
    request_ = request;
    return core::detail::clear_error_state();
  }

  [[nodiscard]] base::StatusCode on_configure(
    const ParameterDictionary &parameters
  ) override
  {
    // Read BL-specific parameters from dictionary
    double val = 0;
    std::string_view text;
    bool bval = false;
    if(parameters.try_get_number("bl_num_layers", val) && val > 0)
      params_.num_layers = static_cast<int>(val);

    if(parameters.try_get_number("bl_first_height", val) && val > 0)
      params_.first_height = val;

    if(parameters.try_get_number("bl_growth_rate", val) && val > 0)
      params_.growth_rate = val;

    if(parameters.try_get_number("bl_first_height_aspect", val) && val > 0)
      params_.first_height_aspect = val;

    if(parameters.try_get_number("bl_aspect", val) && val > 0)
      params_.first_height_aspect = val;

    if(parameters.try_get_number("bl_max_skewness", val))
      params_.max_skewness = val;

    if(parameters.try_get_number("bl_proximity_factor", val))
      params_.proximity_factor = val;

    if(parameters.try_get_number("bl_max_normal_deviation", val))
      params_.max_normal_deviation_deg = val;

    if(parameters.try_get_number("bl_smooth_iterations", val))
      params_.smooth_iterations = static_cast<int>(val);

    if(parameters.try_get_number("bl_min_first_height", val) && val > 0)
      params_.min_first_height = val;

    if(parameters.try_get_number("bl_max_warping", val))
      params_.max_warping_deg = val;

    // y+ auto-calculation
    if(parameters.try_get_number("bl_target_yplus", val) && val > 0)
      params_.target_y_plus = val;

    if(parameters.try_get_number("bl_velocity", val) && val > 0)
      params_.freestream_velocity = val;

    if(parameters.try_get_number("bl_density", val) && val > 0)
      params_.density = val;

    if(parameters.try_get_number("bl_viscosity", val) && val > 0)
      params_.dynamic_viscosity = val;

    if(parameters.try_get_number("bl_ref_length", val) && val > 0)
      params_.reference_length = val;

    // Treatment mode
    if(parameters.try_get_text("bl_treatment", text)) {
      if(text == "collapse")
        params_.treatment = BoundaryLayerParams::ProblematicAreaTreatment::Collapse;
      else if(text == "exclude")
        params_.treatment = BoundaryLayerParams::ProblematicAreaTreatment::Exclude;
      else if(text == "stop")
        params_.treatment = BoundaryLayerParams::ProblematicAreaTreatment::Stop;
    }

    if(parameters.try_get_text("bl_first_height_mode", text)) {
      std::string mode(text);
      std::transform(mode.begin(), mode.end(), mode.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
      });
      if(mode == "aspect" || mode == "proportional") {
        params_.first_height_mode = BoundaryLayerParams::FirstHeightMode::Aspect;
      } else if(mode == "absolute") {
        params_.first_height_mode = BoundaryLayerParams::FirstHeightMode::Absolute;
      }
    }

    if(parameters.try_get_boolean("bl_proportional_height", bval)) {
      params_.first_height_mode = bval
        ? BoundaryLayerParams::FirstHeightMode::Aspect
        : BoundaryLayerParams::FirstHeightMode::Absolute;
    }

    if(parameters.try_get_boolean("bl_smooth_normals", bval))
      params_.smooth_normals = bval;

    if(parameters.try_get_boolean("bl_allow_squeeze", bval))
      params_.allow_squeeze = bval;

    if(parameters.try_get_boolean("bl_allow_proximity", bval))
      params_.allow_proximity = bval;

    // bl_face_ids: comma-separated list of topology face IDs to grow BL on
    // e.g. "0" or "0,1,3" — if empty/missing, grow on ALL faces
    std::string_view face_ids_text;
    if(parameters.try_get_text("bl_face_ids", face_ids_text)) {
      bl_face_ids_.clear();
      std::string s(face_ids_text);
      std::istringstream iss(s);
      std::string token;
      while(std::getline(iss, token, ',')) {
        if(!token.empty()) {
          bl_face_ids_.insert(static_cast<std::uint32_t>(std::stoul(token)));
        }
      }
      bl_log_info("BL face IDs filter: %zu faces specified", bl_face_ids_.size());
    }

    // bl_target_region: which region to push BL into (-1 = auto centroid)
    if(parameters.try_get_number("bl_target_region", val))
      bl_target_region_ = static_cast<int>(val);

    // bl_target_point: "x,y,z" — BL pushes toward this point (per-face direction)
    std::string_view tp_text;
    if(parameters.try_get_text("bl_target_point", tp_text)) {
      std::string s(tp_text);
      std::istringstream iss(s);
      std::string coord;
      int ci = 0;
      while(std::getline(iss, coord, ',') && ci < 3) {
        bl_target_point_[ci++] = std::stod(coord);
      }
      if(ci == 3) {
        has_target_point_ = true;
      }
    }

    return core::detail::clear_error_state();
  }

  [[nodiscard]] base::StatusCode on_generate(Domain &output) override
  {
    auto total_start = std::chrono::high_resolution_clock::now();

    // -- 1. Extract surface from the existing Domain ------------------
    // The output Domain already contains the surface mesh (from the
    // proxy path or from a prior surface mesher run). We read
    // surface face entity_groups as the BL input.
    BLWorkData work;
    auto status = extract_surface_to_work_data(output, work);
    if(status != base::StatusCode::ok) {
      return status;
    }

    // -- 2b. Determine normal orientation from region detection ------
    bool force_flip = false;
    if(bl_target_region_ >= 0) {
      auto regions = detect_regions(output);
      if(bl_target_region_ < static_cast<int>(regions.regions.size())) {
        // Find a BL face and check which side faces the target region
        // Use the first extracted BL face to determine orientation
        // The region detector's face_to_regions uses the same global face indexing
        // We need to map our work faces back to the detector's face list

        // Build a mapping: (entity_group_id, face_index) -> detector face index
        std::unordered_map<std::uint64_t, std::uint32_t> ref_to_detector_idx;
        for(std::uint32_t di = 0; di < regions.face_refs.size(); ++di) {
          auto key = (static_cast<std::uint64_t>(regions.face_refs[di].entity_group) << 32U) |
                     regions.face_refs[di].index;
          ref_to_detector_idx[key] = di;
        }

        // Check the first BL face: which side faces the target region?
        int pos_votes = 0, neg_votes = 0;
        for(std::size_t wi = 0; wi < work.surface_node_refs.size() && wi < 1; ++wi) {
          // We stored face_refs during extraction — but we don't have per-face refs in work
          // Instead, check a few faces from the work data
        }

        // Use the NON-TARGET region's interior point to determine direction.
        // BL faces sit on the boundary between two regions. Normals should
        // point AWAY from the non-target region (= toward the target region).
        // Using the non-target region's interior (close to the surface) is
        // more reliable than the target region's interior (which may be far away).

        // Find which non-target region shares the BL faces.
        // The BL faces form the boundary of the target region AND one other region.
        // Use the centroid of BL faces to find the nearest non-target region.
        Vec3 bl_centroid = {0, 0, 0};
        const auto &nodes = work.surface_nodes;
        for(const auto &nd : nodes) {
          bl_centroid[0] += nd[0]; bl_centroid[1] += nd[1]; bl_centroid[2] += nd[2];
        }
        if(!nodes.empty()) {
          bl_centroid[0] /= nodes.size();
          bl_centroid[1] /= nodes.size();
          bl_centroid[2] /= nodes.size();
        }

        // Find the non-target region whose interior point is closest to BL centroid
        int other_region = -1;
        double min_dist = 1e30;
        for(int ri = 0; ri < static_cast<int>(regions.regions.size()); ++ri) {
          if(ri == bl_target_region_) continue;
          const auto &r = regions.regions[ri];
          double dx = r.interior_point[0] - bl_centroid[0];
          double dy = r.interior_point[1] - bl_centroid[1];
          double dz = r.interior_point[2] - bl_centroid[2];
          double d = dx*dx + dy*dy + dz*dz;
          if(d < min_dist) { min_dist = d; other_region = ri; }
        }

        if(other_region >= 0) {
          const auto &other = regions.regions[other_region];
          Vec3 other_pt = {other.interior_point[0], other.interior_point[1], other.interior_point[2]};

          // Check if normals point toward the other (non-target) region
          // If they do -> flip (we want them pointing AWAY from other, toward target)
          int toward_other = 0, away_from_other = 0;
          int sample = std::min(static_cast<int>(work.surface_faces.size()), 100);
          int step_sz = std::max(1, static_cast<int>(work.surface_faces.size()) / sample);
          for(int i = 0; i < static_cast<int>(work.surface_faces.size()); i += step_sz) {
            const auto &face = work.surface_faces[i];
            Vec3 center = {
              (nodes[face.vertices[0]][0] + nodes[face.vertices[1]][0] + nodes[face.vertices[2]][0]) / 3.0,
              (nodes[face.vertices[0]][1] + nodes[face.vertices[1]][1] + nodes[face.vertices[2]][1]) / 3.0,
              (nodes[face.vertices[0]][2] + nodes[face.vertices[1]][2] + nodes[face.vertices[2]][2]) / 3.0
            };
            auto e1 = vec3_sub(nodes[face.vertices[1]], nodes[face.vertices[0]]);
            auto e2 = vec3_sub(nodes[face.vertices[2]], nodes[face.vertices[0]]);
            Vec3 n = vec3_normalized(vec3_cross(e1, e2));
            Vec3 to_other = vec3_sub(other_pt, center);
            double dot = vec3_dot(n, to_other);
            if(dot > 0) toward_other++;
            else away_from_other++;
          }

          // If normals point toward other region -> flip to point toward target
          force_flip = (toward_other > away_from_other);
        }
      }
    }

    // -- 3. Run the BL solver ----------------------------------------─
    BoundaryLayerSolver solver;
    solver.set_params(params_);
    if(has_target_point_) {
      solver.set_target_point(bl_target_point_[0], bl_target_point_[1], bl_target_point_[2]);
    } else {
      solver.set_force_flip(force_flip);
    }
    solver.initialize(work);
    solver.generate();

    // -- 3b. Optional debug exports ----------------------------------
    if constexpr (kEnableBoundaryLayerDebugFileOutput) {
      export_bl_vtk(work, "bl_output.vtk");
      export_top_cap_vtk(work, "bl_top_cap.vtk");
      export_top_cap_obj(work, "bl_top_cap.obj");
    }

    // -- 4. Write results to output Domain ----------------------------
    status = write_bl_to_domain(work, output, output);
    if(status != base::StatusCode::ok) {
      return status;
    }

    // -- 5. Rebuild surface boundary for tet mesher --------------------─
    // BL'd faces are replaced by bl_bottom (wall) + bl_top (cavity boundary).
    // Non-BL faces (e.g. farfield) must be preserved as computational for
    // the tet mesher to form a closed cavity.
    auto is_original_boundary_face_entity_group = [](const EntityGroup &entity_group) {
      if(entity_group.order() != EntityOrder::face ||
         entity_group.role() != EntityGroupRole::computational ||
         !entity_group.is_boundary()) {
        return false;
      }
      const std::string name(entity_group.name());
      if(name.rfind("bl_", 0) == 0 || name.rfind("kept_", 0) == 0) {
        return false;
      }
      return true;
    };

    if(!bl_face_ids_.empty()) {
      // Selective: collect non-BL faces first, then batch-add to new entity_group.
      // NOTE: Must not modify Domain while iterating entity_groups (iterator invalidation).
      struct FaceCopy {
        std::array<EntityRef, 4> nodes;
        int num_nodes;
        geo::TopologyEntityId topo;
      };
      std::vector<FaceCopy> tri_faces_to_keep;
      std::vector<FaceCopy> quad_faces_to_keep;
      std::set<std::string> entity_groups_to_remove;

      for(const auto &entity_group : output.entity_groups()) {
        if(!is_original_boundary_face_entity_group(entity_group)) continue;
        entity_groups_to_remove.insert(std::string(entity_group.name()));

        for(std::uint32_t fi = 0; fi < entity_group.faces().size(); ++fi) {
          EntityRef face_ref = {entity_group.id(), fi};
          auto topo = output.face_topology_owner(face_ref);
          if(!geo::is_valid(topo)) {
            topo = {geo::TopologyDimension::face, entity_group.zone_id()};
          }

          // Skip BL'd faces — they're now in bl_bottom/bl_top entity_groups
          if(bl_face_ids_.count(topo.index)) continue;

          auto nodes = output.face_nodes(face_ref);
          if(nodes.size < 3 || nodes.size > 4) continue;
          FaceCopy fc;
          fc.num_nodes = static_cast<int>(nodes.size);
          fc.topo = topo;
          for(std::size_t ni = 0; ni < nodes.size && ni < 4; ++ni) {
            fc.nodes[ni] = nodes[ni];
          }
          if(fc.num_nodes == 3) {
            tri_faces_to_keep.push_back(fc);
          } else if(fc.num_nodes == 4) {
            quad_faces_to_keep.push_back(fc);
          }
        }
      }

      for(const auto &name : entity_groups_to_remove) {
        output.remove_entity_groups_with_prefix(name);
      }

      EntityGroupIndex kept_tri_entity_group = invalid_index;
      if(!tri_faces_to_keep.empty()) {
        EntityGroupDefinition kept_tri_def;
        kept_tri_def.order = EntityOrder::face;
        kept_tri_def.name = "kept_surface_tri_faces";
        kept_tri_def.boundary = true;
        kept_tri_def.role = EntityGroupRole::computational;
        kept_tri_def.default_kind = EntityKind::face_triangle;
        kept_tri_entity_group = output.create_entity_group(kept_tri_def);
      }

      EntityGroupIndex kept_quad_entity_group = invalid_index;
      if(!quad_faces_to_keep.empty()) {
        EntityGroupDefinition kept_quad_def;
        kept_quad_def.order = EntityOrder::face;
        kept_quad_def.name = "kept_surface_quad_faces";
        kept_quad_def.boundary = true;
        kept_quad_def.role = EntityGroupRole::computational;
        kept_quad_def.default_kind = EntityKind::face_quad;
        kept_quad_entity_group = output.create_entity_group(kept_quad_def);
      }

      for(const auto &fc : tri_faces_to_keep) {
        std::array<EntityRef, 3> refs = {fc.nodes[0], fc.nodes[1], fc.nodes[2]};
        auto new_ref = output.add_triangle_face(kept_tri_entity_group, refs);
        if(geo::is_valid(fc.topo)) {
          output.set_face_topology_owner(new_ref, fc.topo);
        }
      }

      for(const auto &fc : quad_faces_to_keep) {
        std::array<EntityRef, 4> refs = {fc.nodes[0], fc.nodes[1], fc.nodes[2], fc.nodes[3]};
        auto new_ref = output.add_quad_face(kept_quad_entity_group, refs);
        if(geo::is_valid(fc.topo)) {
          output.set_face_topology_owner(new_ref, fc.topo);
        }
      }

      bl_log_info(
        "Kept %zu non-BL surface faces for tet boundary (%zu tris, %zu quads)",
        tri_faces_to_keep.size() + quad_faces_to_keep.size(),
        tri_faces_to_keep.size(),
        quad_faces_to_keep.size());
    } else {
      // All faces are BL'd: remove all original computational face entity_groups.
      // bl_top_cap replaces the original surface as the new boundary.
      // For deactivated faces (no BL generated), the original face is already
      // included in bl_top by write_bl_to_domain (as an un-pushed top cap face).
      //
      // Collect original face entity_group names to remove (avoid iterator invalidation)
      std::set<std::string> entity_groups_to_remove;
      for(const auto &entity_group : output.entity_groups()) {
        if(is_original_boundary_face_entity_group(entity_group)) {
          entity_groups_to_remove.insert(std::string(entity_group.name()));
        }
      }
      for(const auto &name : entity_groups_to_remove) {
        output.remove_entity_groups_with_prefix(name);
      }
      bl_log_info(
        "Removed %zu original face entity groups, bl_top is now the boundary",
        entity_groups_to_remove.size());
    }

    auto total_end = std::chrono::high_resolution_clock::now();
    double total_ms = std::chrono::duration<double, std::milli>(total_end - total_start).count();
    bl_log_info("Boundary Layer Mesher complete in %.1f ms", total_ms);

    return core::detail::clear_error_state();
  }

private:
  // Extract surface mesh from Domain into BLWorkData
  [[nodiscard]] base::StatusCode extract_surface_to_work_data(
    const Domain &domain,
    BLWorkData &work
  )
  {
    work.clear();

    std::unordered_map<std::uint64_t, std::uint32_t> node_map;
    std::unordered_map<std::uint64_t, std::uint32_t> non_bl_node_map;
    auto pack_ref = [](EntityRef ref) -> std::uint64_t {
      return (static_cast<std::uint64_t>(ref.entity_group) << 32U) |
             static_cast<std::uint64_t>(ref.index);
    };

    auto ensure_node = [&](EntityRef node_ref) -> std::uint32_t {
      auto key = pack_ref(node_ref);
      auto it = node_map.find(key);
      if(it != node_map.end()) return it->second;

      const auto &node = domain.node(node_ref);
      auto idx = static_cast<std::uint32_t>(work.surface_nodes.size());
      work.surface_nodes.push_back({node.coordinates[0], node.coordinates[1], node.coordinates[2]});
      work.surface_node_refs.push_back(node_ref);
      node_map[key] = idx;
      return idx;
    };

    auto ensure_non_bl_node = [&](EntityRef node_ref) -> std::uint32_t {
      auto key = pack_ref(node_ref);
      auto it = non_bl_node_map.find(key);
      if(it != non_bl_node_map.end()) return it->second;

      const auto &node = domain.node(node_ref);
      auto idx = static_cast<std::uint32_t>(work.non_bl_nodes.size());
      work.non_bl_nodes.push_back({node.coordinates[0], node.coordinates[1], node.coordinates[2]});
      work.non_bl_node_refs.push_back(node_ref);
      non_bl_node_map[key] = idx;
      return idx;
    };

    auto resolve_topology_owner = [&](const EntityGroup &entity_group, EntityRef face_ref) {
      auto topo = domain.face_topology_owner(face_ref);
      if(!geo::is_valid(topo)) {
        topo = {geo::TopologyDimension::face, entity_group.zone_id()};
      }
      return topo;
    };

    // Collect face topology IDs for filtering
    const bool filter_by_face_id = !bl_face_ids_.empty();

    // All boundary (computational) face entity_groups are BL candidates
    for(const auto &entity_group : domain.entity_groups()) {
      if(entity_group.order() != EntityOrder::face ||
         entity_group.role() != EntityGroupRole::computational ||
         !entity_group.is_boundary()) {
        continue;
      }

      for(std::uint32_t fi = 0; fi < entity_group.faces().size(); ++fi) {
        EntityRef face_ref = {entity_group.id(), fi};
        const auto topo_id = resolve_topology_owner(entity_group, face_ref);
        const bool selected_for_bl =
          !filter_by_face_id ||
          bl_face_ids_.find(topo_id.index) != bl_face_ids_.end();

        const auto face_nodes = domain.face_nodes(face_ref);
        if(face_nodes.size < 3 || face_nodes.size > 4) continue;

        if(selected_for_bl) {
          SurfaceFace sf;
          sf.num_vertices = static_cast<int>(face_nodes.size);
          sf.topology_owner = topo_id.index;

          for(std::size_t ni = 0; ni < face_nodes.size; ++ni) {
            sf.vertices[ni] = ensure_node(face_nodes[ni]);
          }
          work.surface_faces.push_back(sf);
          continue;
        }

        if(!filter_by_face_id) {
          continue;
        }

        std::array<std::uint32_t, 4> verts = {0U, 0U, 0U, 0U};
        for(std::size_t ni = 0; ni < face_nodes.size; ++ni) {
          verts[ni] = ensure_non_bl_node(face_nodes[ni]);
        }
        work.non_bl_face_indices.push_back(static_cast<std::uint32_t>(work.non_bl_face_verts.size()));
        work.non_bl_face_verts.push_back(verts);
        work.non_bl_face_nv.push_back(static_cast<int>(face_nodes.size));
      }
    }

    bl_log_info(
      "Extracted %zu BL surface nodes, %zu BL surface faces, %zu non-BL boundary faces",
      work.surface_nodes.size(),
      work.surface_faces.size(),
      work.non_bl_face_verts.size());

    if(work.surface_faces.empty()) {
      return core::detail::publish_error(
        base::StatusCode::internal_error,
        "No surface faces found for boundary layer generation."
      );
    }

    return core::detail::clear_error_state();
  }

  // Write BL layers + top cap to output Domain
  [[nodiscard]] base::StatusCode write_bl_to_domain(
    const BLWorkData &work,
    const Domain &surface_domain,
    Domain &output
  )
  {
    static_cast<void>(surface_domain);
    const auto num_layers = work.layer_positions.size() - 1; // Excluding base
    if(num_layers == 0) {
      bl_log_warn("No boundary layers generated");
      return core::detail::clear_error_state();
    }

    bl_log_info("Writing %zu BL layers to Domain", num_layers);

    // We add BL entity_groups to the existing output Domain (which contains the surface mesh)

    // -- Create entity_groups --
    EntityGroupDefinition node_def;
    node_def.order = EntityOrder::node;
    node_def.name = "bl_nodes";
    node_def.role = EntityGroupRole::computational;
    auto node_entity_group = output.create_entity_group(node_def);

    EntityGroupDefinition tri_face_def;
    tri_face_def.order = EntityOrder::face;
    tri_face_def.name = "bl_tri_faces";
    tri_face_def.default_kind = EntityKind::face_triangle;
    tri_face_def.role = EntityGroupRole::annotation; // Internal BL faces, not for tet mesher
    auto tri_face_entity_group = output.create_entity_group(tri_face_def);

    EntityGroupDefinition quad_face_def;
    quad_face_def.order = EntityOrder::face;
    quad_face_def.name = "bl_quad_faces";
    quad_face_def.default_kind = EntityKind::face_quad;
    quad_face_def.role = EntityGroupRole::annotation; // Internal BL faces, not for tet mesher
    auto quad_face_entity_group = output.create_entity_group(quad_face_def);

    EntityGroupDefinition prism_cell_def;
    prism_cell_def.order = EntityOrder::cell;
    prism_cell_def.name = "bl_prism_cells";
    prism_cell_def.default_kind = EntityKind::cell_prism;
    prism_cell_def.role = EntityGroupRole::computational;
    auto prism_cell_entity_group = output.create_entity_group(prism_cell_def);

    EntityGroupDefinition hexa_cell_def;
    hexa_cell_def.order = EntityOrder::cell;
    hexa_cell_def.name = "bl_hexa_cells";
    hexa_cell_def.default_kind = EntityKind::cell_hexa;
    hexa_cell_def.role = EntityGroupRole::computational;
    auto hexa_cell_entity_group = output.create_entity_group(hexa_cell_def);

    EntityGroupDefinition top_tri_face_def;
    top_tri_face_def.order = EntityOrder::face;
    top_tri_face_def.name = "bl_top_tri_faces";
    top_tri_face_def.default_kind = EntityKind::face_triangle;
    top_tri_face_def.role = EntityGroupRole::computational;
    top_tri_face_def.boundary = true;
    auto top_tri_face_entity_group = output.create_entity_group(top_tri_face_def);

    EntityGroupDefinition top_quad_face_def;
    top_quad_face_def.order = EntityOrder::face;
    top_quad_face_def.name = "bl_top_quad_faces";
    top_quad_face_def.default_kind = EntityKind::face_quad;
    top_quad_face_def.role = EntityGroupRole::computational;
    top_quad_face_def.boundary = true;
    auto top_quad_face_entity_group = output.create_entity_group(top_quad_face_def);

    EntityGroupDefinition transition_quad_face_def;
    transition_quad_face_def.order = EntityOrder::face;
    transition_quad_face_def.name = "bl_transition_quad_faces";
    transition_quad_face_def.default_kind = EntityKind::face_quad;
    transition_quad_face_def.role = EntityGroupRole::computational;
    transition_quad_face_def.boundary = true;
    auto transition_quad_face_entity_group = output.create_entity_group(transition_quad_face_def);

    EntityGroupDefinition pyramid_base_face_def;
    pyramid_base_face_def.order = EntityOrder::face;
    pyramid_base_face_def.name = "bl_pyramid_base_quad_faces";
    pyramid_base_face_def.default_kind = EntityKind::face_quad;
    pyramid_base_face_def.role = EntityGroupRole::annotation;
    auto pyramid_base_face_entity_group = output.create_entity_group(pyramid_base_face_def);

    EntityGroupDefinition pyramid_side_face_def;
    pyramid_side_face_def.order = EntityOrder::face;
    pyramid_side_face_def.name = "bl_pyramid_side_tri_faces";
    pyramid_side_face_def.default_kind = EntityKind::face_triangle;
    pyramid_side_face_def.role = EntityGroupRole::computational;
    pyramid_side_face_def.boundary = true;
    auto pyramid_side_face_entity_group = output.create_entity_group(pyramid_side_face_def);

    EntityGroupDefinition pyramid_cell_def;
    pyramid_cell_def.order = EntityOrder::cell;
    pyramid_cell_def.name = "bl_pyramid_cells";
    pyramid_cell_def.default_kind = EntityKind::cell_pyramid;
    pyramid_cell_def.role = EntityGroupRole::computational;
    auto pyramid_cell_entity_group = output.create_entity_group(pyramid_cell_def);

    EntityGroupDefinition bot_tri_face_def;
    bot_tri_face_def.order = EntityOrder::face;
    bot_tri_face_def.name = "bl_bottom_tri_faces";
    bot_tri_face_def.default_kind = EntityKind::face_triangle;
    bot_tri_face_def.role = EntityGroupRole::annotation; // Wall boundary, not for tet mesher
    bot_tri_face_def.boundary = true;
    auto bot_tri_face_entity_group = output.create_entity_group(bot_tri_face_def);

    EntityGroupDefinition bot_quad_face_def;
    bot_quad_face_def.order = EntityOrder::face;
    bot_quad_face_def.name = "bl_bottom_quad_faces";
    bot_quad_face_def.default_kind = EntityKind::face_quad;
    bot_quad_face_def.role = EntityGroupRole::annotation; // Wall boundary, not for tet mesher
    bot_quad_face_def.boundary = true;
    auto bot_quad_face_entity_group = output.create_entity_group(bot_quad_face_def);

    // -- Add all layer nodes --
    // node_refs[layer][surface_node_index] = EntityRef in output Domain
    const auto nn = work.surface_nodes.size();
    std::vector<std::vector<EntityRef>> node_refs(num_layers + 1);
    node_refs[0] = work.surface_node_refs;

    for(std::size_t layer = 1; layer <= num_layers; ++layer) {
      node_refs[layer].resize(nn);
      for(std::uint32_t ni = 0; ni < nn; ++ni) {
        const auto &pos = work.layer_positions[layer][ni];
        node_refs[layer][ni] = output.add_node(
          node_entity_group, {pos[0], pos[1], pos[2]});
      }
    }

    bl_log_info(
      "Added %zu nodes (%zu layers x %zu surface nodes)",
      (num_layers + 1) * nn,
      num_layers + 1,
      nn);

    std::vector<int> face_levels(work.surface_faces.size(), 0);
    for(std::size_t fi = 0; fi < work.surface_faces.size(); ++fi) {
      face_levels[fi] = std::max(0, std::min(work.max_face_level[fi], static_cast<int>(num_layers)));
    }

    auto face_topology = [&](std::uint32_t face_index) {
      return geo::TopologyEntityId {
        geo::TopologyDimension::face,
        work.surface_faces[face_index].topology_owner
      };
    };

    auto node_position = [&](EntityRef ref) -> Vec3 {
      const auto &node = output.node(ref);
      return {node.coordinates[0], node.coordinates[1], node.coordinates[2]};
    };

    const Vec3 target_point = {
      bl_target_point_[0],
      bl_target_point_[1],
      bl_target_point_[2]
    };

    auto reverse_quad = [](const std::array<EntityRef, 4> &refs) {
      return std::array<EntityRef, 4> {refs[0], refs[3], refs[2], refs[1]};
    };

    auto reverse_tri = [](const std::array<EntityRef, 3> &refs) {
      return std::array<EntityRef, 3> {refs[0], refs[2], refs[1]};
    };

    auto quad_center = [&](const std::array<EntityRef, 4> &refs) -> Vec3 {
      Vec3 sum {0.0, 0.0, 0.0};
      for(int i = 0; i < 4; ++i) {
        const auto p = node_position(refs[static_cast<std::size_t>(i)]);
        sum[0] += p[0];
        sum[1] += p[1];
        sum[2] += p[2];
      }
      return {sum[0] * 0.25, sum[1] * 0.25, sum[2] * 0.25};
    };

    auto tri_center = [&](const std::array<EntityRef, 3> &refs) -> Vec3 {
      Vec3 sum {0.0, 0.0, 0.0};
      for(int i = 0; i < 3; ++i) {
        const auto p = node_position(refs[static_cast<std::size_t>(i)]);
        sum[0] += p[0];
        sum[1] += p[1];
        sum[2] += p[2];
      }
      return {sum[0] / 3.0, sum[1] / 3.0, sum[2] / 3.0};
    };

    auto quad_normal = [&](const std::array<EntityRef, 4> &refs) -> Vec3 {
      const auto p0 = node_position(refs[0]);
      const auto p1 = node_position(refs[1]);
      const auto p2 = node_position(refs[2]);
      const auto p3 = node_position(refs[3]);
      return vec3_normalized(vec3_add(
        vec3_cross(vec3_sub(p1, p0), vec3_sub(p2, p0)),
        vec3_cross(vec3_sub(p2, p0), vec3_sub(p3, p0))
      ));
    };

    auto tri_normal = [&](const std::array<EntityRef, 3> &refs) -> Vec3 {
      const auto p0 = node_position(refs[0]);
      const auto p1 = node_position(refs[1]);
      const auto p2 = node_position(refs[2]);
      return vec3_normalized(vec3_cross(vec3_sub(p1, p0), vec3_sub(p2, p0)));
    };

    auto orient_quad_outward = [&](std::array<EntityRef, 4> refs, const Vec3 &inside_point) {
      const auto normal = quad_normal(refs);
      const auto center = quad_center(refs);
      if(vec3_dot(normal, vec3_sub(center, inside_point)) < 0.0) {
        refs = reverse_quad(refs);
      }
      return refs;
    };

    auto orient_tri_outward = [&](std::array<EntityRef, 3> refs, const Vec3 &inside_point) {
      const auto normal = tri_normal(refs);
      const auto center = tri_center(refs);
      if(vec3_dot(normal, vec3_sub(center, inside_point)) < 0.0) {
        refs = reverse_tri(refs);
      }
      return refs;
    };

    auto polygon_scale = [&](const std::array<EntityRef, 4> &refs, int nv) -> double {
      double min_len = std::numeric_limits<double>::max();
      for(int i = 0; i < nv; ++i) {
        const auto a = node_position(refs[static_cast<std::size_t>(i)]);
        const auto b = node_position(refs[static_cast<std::size_t>((i + 1) % nv)]);
        const double len = vec3_length(vec3_sub(b, a));
        if(len > 1e-12) {
          min_len = std::min(min_len, len);
        }
      }
      if(!std::isfinite(min_len) || min_len == std::numeric_limits<double>::max()) {
        return 1e-6;
      }
      return std::max(min_len, 1e-6);
    };

    auto pyramid_volume_measure =
      [&](const std::array<EntityRef, 4> &base_refs,
          const Vec3 &apex) -> double {
        const auto p0 = node_position(base_refs[0]);
        const auto p1 = node_position(base_refs[1]);
        const auto p3 = node_position(base_refs[3]);
        return vec3_dot(
          vec3_cross(vec3_sub(p1, p0), vec3_sub(p3, p0)),
          vec3_sub(apex, p0));
      };

    constexpr double kShellIntersectionEps = 1e-9;
    constexpr double kMaximumSharedEdgeDihedralDegrees = 150.0;
    constexpr std::uint64_t kInvalidVertexKey = std::numeric_limits<std::uint64_t>::max();

    struct Triangle3 {
      Vec3 a;
      Vec3 b;
      Vec3 c;
    };

    struct ShellTriangle {
      Triangle3 tri;
      std::array<std::uint64_t, 3> vertex_keys {kInvalidVertexKey, kInvalidVertexKey, kInvalidVertexKey};
      AABB box;
    };

    struct ShellEdge {
      Vec3 a;
      Vec3 b;
      std::array<std::uint64_t, 2> vertex_keys {kInvalidVertexKey, kInvalidVertexKey};
      AABB box;
    };

    struct ValidationFace {
      std::array<Vec3, 4> points {};
      std::array<std::uint64_t, 4> vertex_keys {
        kInvalidVertexKey, kInvalidVertexKey, kInvalidVertexKey, kInvalidVertexKey};
      EntityRef output_face_ref {};
      int num_vertices = 0;
      bool active = true;
      AABB box;
    };

    auto pack_ref = [](EntityRef ref) -> std::uint64_t {
      return (static_cast<std::uint64_t>(ref.entity_group) << 32U) |
             static_cast<std::uint64_t>(ref.index);
    };

    auto unpack_ref = [](std::uint64_t key) -> EntityRef {
      if(key == kInvalidVertexKey || (key & (std::uint64_t {1} << 63)) != 0) {
        return {};
      }
      return EntityRef {
        static_cast<EntityGroupIndex>(key >> 32U),
        static_cast<std::uint32_t>(key & 0xffffffffULL)
      };
    };

    auto make_triangle_aabb = [&](const Triangle3 &tri) {
      AABB box;
      box.expand(tri.a);
      box.expand(tri.b);
      box.expand(tri.c);
      box.pad(kShellIntersectionEps);
      return box;
    };

    auto make_edge_aabb = [&](const Vec3 &a, const Vec3 &b) {
      AABB box;
      box.expand(a);
      box.expand(b);
      box.pad(kShellIntersectionEps);
      return box;
    };

    auto triangle_area = [&](const Triangle3 &tri) -> double {
      return 0.5 * vec3_length(vec3_cross(
        vec3_sub(tri.b, tri.a),
        vec3_sub(tri.c, tri.a)));
    };

    auto make_triangle_from_refs = [&](const std::array<EntityRef, 3> &refs) {
      return Triangle3 {
        node_position(refs[0]),
        node_position(refs[1]),
        node_position(refs[2])
      };
    };

    auto make_shell_triangle = [&](const std::array<EntityRef, 3> &refs) {
      ShellTriangle tri;
      tri.tri = make_triangle_from_refs(refs);
      tri.vertex_keys = {pack_ref(refs[0]), pack_ref(refs[1]), pack_ref(refs[2])};
      tri.box = make_triangle_aabb(tri.tri);
      return tri;
    };

    auto make_candidate_triangle =
      [&](EntityRef a, EntityRef b, const Vec3 &apex) {
        ShellTriangle tri;
        tri.tri = {
          node_position(a),
          node_position(b),
          apex
        };
        tri.vertex_keys = {pack_ref(a), pack_ref(b), kInvalidVertexKey};
        tri.box = make_triangle_aabb(tri.tri);
        return tri;
      };

    auto make_shell_edge = [&](EntityRef a, EntityRef b) {
      ShellEdge edge;
      edge.a = node_position(a);
      edge.b = node_position(b);
      edge.vertex_keys = {pack_ref(a), pack_ref(b)};
      edge.box = make_edge_aabb(edge.a, edge.b);
      return edge;
    };

    auto make_candidate_edge =
      [&](std::uint64_t key_a, const Vec3 &a, std::uint64_t key_b, const Vec3 &b) {
        ShellEdge edge;
        edge.a = a;
        edge.b = b;
        edge.vertex_keys = {key_a, key_b};
        edge.box = make_edge_aabb(edge.a, edge.b);
        return edge;
      };

    auto make_validation_face_from_refs =
      [&](const std::array<EntityRef, 4> &refs, int nv) {
        ValidationFace face;
        face.num_vertices = nv;
        face.active = true;
        for(int i = 0; i < nv; ++i) {
          face.points[static_cast<std::size_t>(i)] = node_position(refs[static_cast<std::size_t>(i)]);
          face.vertex_keys[static_cast<std::size_t>(i)] = pack_ref(refs[static_cast<std::size_t>(i)]);
          face.box.expand(face.points[static_cast<std::size_t>(i)]);
        }
        face.box.pad(kShellIntersectionEps);
        return face;
      };

    auto make_validation_face_from_shell_triangle =
      [&](const ShellTriangle &tri) {
        ValidationFace face;
        face.num_vertices = 3;
        face.active = true;
        face.points[0] = tri.tri.a;
        face.points[1] = tri.tri.b;
        face.points[2] = tri.tri.c;
        face.vertex_keys[0] = tri.vertex_keys[0];
        face.vertex_keys[1] = tri.vertex_keys[1];
        face.vertex_keys[2] = tri.vertex_keys[2];
        face.box = tri.box;
        return face;
      };

    auto validation_face_triangle_count = [](const ValidationFace &face) noexcept {
      return face.num_vertices == 4 ? 2 : 1;
    };

    auto validation_face_triangle =
      [&](const ValidationFace &face, int tri_index) noexcept {
        if(face.num_vertices == 4 && tri_index != 0) {
          return Triangle3 {
            face.points[0],
            face.points[2],
            face.points[3]
          };
        }

        return Triangle3 {
          face.points[0],
          face.points[1],
          face.points[2]
        };
      };

    auto validation_face_triangle_keys =
      [&](const ValidationFace &face, int tri_index) noexcept {
        if(face.num_vertices == 4 && tri_index != 0) {
          return std::array<std::uint64_t, 3> {
            face.vertex_keys[0],
            face.vertex_keys[2],
            face.vertex_keys[3]
          };
        }

        return std::array<std::uint64_t, 3> {
          face.vertex_keys[0],
          face.vertex_keys[1],
          face.vertex_keys[2]
        };
      };

    auto validation_face_edge =
      [&](const ValidationFace &face, int edge_index) {
        const int next = (edge_index + 1) % face.num_vertices;
        return make_candidate_edge(
          kInvalidVertexKey,
          face.points[static_cast<std::size_t>(edge_index)],
          kInvalidVertexKey,
          face.points[static_cast<std::size_t>(next)]);
      };

    auto validation_face_normal =
      [&](const ValidationFace &face) -> Vec3 {
        if(face.num_vertices == 4) {
          return vec3_normalized(vec3_add(
            vec3_cross(vec3_sub(face.points[1], face.points[0]),
                       vec3_sub(face.points[2], face.points[0])),
            vec3_cross(vec3_sub(face.points[2], face.points[0]),
                       vec3_sub(face.points[3], face.points[0]))));
        }
        if(face.num_vertices == 3) {
          return vec3_normalized(vec3_cross(
            vec3_sub(face.points[1], face.points[0]),
            vec3_sub(face.points[2], face.points[0])));
        }
        return {0.0, 0.0, 0.0};
      };

    auto shared_edge_dihedral_degrees =
      [&](const ValidationFace &quad_face,
          int quad_edge_index,
          const ValidationFace &tri_face) -> double {
        if(quad_face.num_vertices != 4 ||
           quad_edge_index < 0 ||
           quad_edge_index >= 4 ||
           tri_face.num_vertices != 3) {
          return 360.0;
        }

        const auto a = quad_face.points[static_cast<std::size_t>(quad_edge_index)];
        const auto b =
          quad_face.points[static_cast<std::size_t>((quad_edge_index + 1) % 4)];
        const auto edge_dir = vec3_normalized(vec3_sub(b, a));
        if(vec3_length(edge_dir) <= 1e-12) {
          return 360.0;
        }

        const auto quad_n = validation_face_normal(quad_face);
        const auto tri_n = validation_face_normal(tri_face);
        if(vec3_length(quad_n) <= 1e-12 || vec3_length(tri_n) <= 1e-12) {
          return 360.0;
        }

        constexpr double kPi = 3.14159265358979323846;
        const double dot_n = std::clamp(vec3_dot(quad_n, tri_n), -1.0, 1.0);
        const double angle_between_normals = std::acos(dot_n);
        const auto cross_n = vec3_cross(quad_n, tri_n);
        const double cross_dot_edge = vec3_dot(cross_n, edge_dir);

        // Directed dihedral following the user's convention:
        //   angle a = acos(n1 . n2)
        //   if edge_dir and (n1 x n2) are opposite -> pi - a
        //   if edge_dir and (n1 x n2) are aligned  -> pi + a
        double dihedral = kPi;
        if(cross_dot_edge >= 0.0) {
          dihedral += angle_between_normals;
        } else {
          dihedral -= angle_between_normals;
        }
        return dihedral * 180.0 / kPi;
      };

    auto oriented_shared_edge_dihedral_degrees =
      [&](const std::array<EntityRef, 4> &quad_refs,
          int quad_edge_index,
          const ValidationFace &tri_face) -> double {
        return shared_edge_dihedral_degrees(
          make_validation_face_from_refs(quad_refs, 4),
          quad_edge_index,
          tri_face);
      };

    auto shell_triangles_share_vertex =
      [=](const std::array<std::uint64_t, 3> &lhs,
          const std::array<std::uint64_t, 3> &rhs) noexcept {
        for(const auto lhs_key : lhs) {
          if(lhs_key == kInvalidVertexKey) {
            continue;
          }
          for(const auto rhs_key : rhs) {
            if(rhs_key == kInvalidVertexKey) {
              continue;
            }
            if(lhs_key == rhs_key) {
              return true;
            }
          }
        }
        return false;
      };

    auto count_shared_triangle_vertices =
      [=](const std::array<std::uint64_t, 3> &lhs,
          const std::array<std::uint64_t, 3> &rhs) noexcept {
        int count = 0;
        for(const auto lhs_key : lhs) {
          if(lhs_key == kInvalidVertexKey) {
            continue;
          }
          for(const auto rhs_key : rhs) {
            if(rhs_key == kInvalidVertexKey) {
              continue;
            }
            if(lhs_key == rhs_key) {
              ++count;
              break;
            }
          }
        }
        return count;
      };

    auto face_shares_candidate_base_edge =
      [=](const ValidationFace &face,
          const ShellTriangle &candidate) noexcept {
        const auto a = candidate.vertex_keys[0];
        const auto b = candidate.vertex_keys[1];
        if(a == kInvalidVertexKey || b == kInvalidVertexKey) {
          return false;
        }
        for(int edge_index = 0; edge_index < face.num_vertices; ++edge_index) {
          const auto fa = face.vertex_keys[static_cast<std::size_t>(edge_index)];
          const auto fb =
            face.vertex_keys[static_cast<std::size_t>((edge_index + 1) % face.num_vertices)];
          if((fa == a && fb == b) || (fa == b && fb == a)) {
            return true;
          }
        }
        return false;
      };

    auto points_close = [&](const Vec3 &lhs, const Vec3 &rhs, double tol) noexcept {
      return vec3_length(vec3_sub(lhs, rhs)) <= tol;
    };

    auto edges_same_geometry =
      [&](const ShellEdge &lhs, const ShellEdge &rhs, double tol) noexcept {
        return (points_close(lhs.a, rhs.a, tol) && points_close(lhs.b, rhs.b, tol)) ||
               (points_close(lhs.a, rhs.b, tol) && points_close(lhs.b, rhs.a, tol));
      };

    auto edges_have_nonconforming_overlap =
      [&](const ShellEdge &lhs, const ShellEdge &rhs) noexcept {
        const Vec3 lhs_dir = vec3_sub(lhs.b, lhs.a);
        const Vec3 rhs_dir = vec3_sub(rhs.b, rhs.a);
        const double lhs_len = vec3_length(lhs_dir);
        const double rhs_len = vec3_length(rhs_dir);
        if(lhs_len <= kShellIntersectionEps || rhs_len <= kShellIntersectionEps) {
          return false;
        }

        const double geom_tol = std::max(1e-9, 1e-8 * std::max(lhs_len, rhs_len));
        const double parallel_measure =
          vec3_length(vec3_cross(lhs_dir, rhs_dir)) / (lhs_len * rhs_len);
        if(parallel_measure > 1e-8) {
          return false;
        }

        const auto point_line_distance =
          [&](const Vec3 &point) noexcept {
            return vec3_length(vec3_cross(lhs_dir, vec3_sub(point, lhs.a))) / lhs_len;
          };
        if(point_line_distance(rhs.a) > geom_tol || point_line_distance(rhs.b) > geom_tol) {
          return false;
        }

        if(edges_same_geometry(lhs, rhs, geom_tol)) {
          return false;
        }

        int axis = 0;
        double axis_span = std::abs(lhs_dir[0]);
        for(int i = 1; i < 3; ++i) {
          const double cur_span = std::abs(lhs_dir[i]);
          if(cur_span > axis_span) {
            axis = i;
            axis_span = cur_span;
          }
        }
        if(axis_span <= geom_tol) {
          return false;
        }

        const double lhs_min = std::min(lhs.a[axis], lhs.b[axis]);
        const double lhs_max = std::max(lhs.a[axis], lhs.b[axis]);
        const double rhs_min = std::min(rhs.a[axis], rhs.b[axis]);
        const double rhs_max = std::max(rhs.a[axis], rhs.b[axis]);
        const double overlap = std::min(lhs_max, rhs_max) - std::max(lhs_min, rhs_min);
        return overlap > geom_tol;
      };

    auto face_triangle_count = [](const SurfaceFace &face) noexcept {
      return face.num_vertices == 4 ? 2 : 1;
    };

    auto face_triangle =
      [&](const std::vector<Vec3> &positions,
          const SurfaceFace &face,
          int tri_index) noexcept {
        if(face.num_vertices == 4 && tri_index != 0) {
          return Triangle3 {
            positions[face.vertices[0]],
            positions[face.vertices[2]],
            positions[face.vertices[3]]
          };
        }

        return Triangle3 {
          positions[face.vertices[0]],
          positions[face.vertices[1]],
          positions[face.vertices[2]]
        };
      };

    auto project_point_2d =
      [](const Vec3 &point, int drop_axis) noexcept -> std::array<double, 2> {
        switch(drop_axis) {
        case 0: return {point[1], point[2]};
        case 1: return {point[0], point[2]};
        default: return {point[0], point[1]};
        }
      };

    auto orient_2d =
      [](const std::array<double, 2> &a,
         const std::array<double, 2> &b,
         const std::array<double, 2> &c) noexcept {
        return (b[0] - a[0]) * (c[1] - a[1]) -
               (b[1] - a[1]) * (c[0] - a[0]);
      };

    auto point_in_triangle_2d_strict =
      [&](const std::array<double, 2> &p,
          const std::array<double, 2> &a,
          const std::array<double, 2> &b,
          const std::array<double, 2> &c) noexcept {
        const double ab = orient_2d(a, b, p);
        const double bc = orient_2d(b, c, p);
        const double ca = orient_2d(c, a, p);

        const bool positive =
          ab > kShellIntersectionEps &&
          bc > kShellIntersectionEps &&
          ca > kShellIntersectionEps;
        const bool negative =
          ab < -kShellIntersectionEps &&
          bc < -kShellIntersectionEps &&
          ca < -kShellIntersectionEps;
        return positive || negative;
      };

    auto segments_intersect_2d_strict =
      [&](const std::array<double, 2> &a0,
          const std::array<double, 2> &a1,
          const std::array<double, 2> &b0,
          const std::array<double, 2> &b1) noexcept {
        const double o1 = orient_2d(a0, a1, b0);
        const double o2 = orient_2d(a0, a1, b1);
        const double o3 = orient_2d(b0, b1, a0);
        const double o4 = orient_2d(b0, b1, a1);

        const bool separated =
          (o1 > kShellIntersectionEps && o2 > kShellIntersectionEps) ||
          (o1 < -kShellIntersectionEps && o2 < -kShellIntersectionEps) ||
          (o3 > kShellIntersectionEps && o4 > kShellIntersectionEps) ||
          (o3 < -kShellIntersectionEps && o4 < -kShellIntersectionEps);
        if(separated) {
          return false;
        }

        const bool straddles_a =
          (o1 > kShellIntersectionEps && o2 < -kShellIntersectionEps) ||
          (o1 < -kShellIntersectionEps && o2 > kShellIntersectionEps);
        const bool straddles_b =
          (o3 > kShellIntersectionEps && o4 < -kShellIntersectionEps) ||
          (o3 < -kShellIntersectionEps && o4 > kShellIntersectionEps);
        return straddles_a && straddles_b;
      };

    auto coplanar_triangles_intersect_strict =
      [&](const Triangle3 &lhs,
          const Triangle3 &rhs,
          const Vec3 &plane_normal) noexcept {
        int drop_axis = 0;
        double max_abs = std::abs(plane_normal[0]);
        for(int axis = 1; axis < 3; ++axis) {
          const double cur_abs = std::abs(plane_normal[axis]);
          if(cur_abs > max_abs) {
            max_abs = cur_abs;
            drop_axis = axis;
          }
        }

        const auto la = project_point_2d(lhs.a, drop_axis);
        const auto lb = project_point_2d(lhs.b, drop_axis);
        const auto lc = project_point_2d(lhs.c, drop_axis);
        const auto ra = project_point_2d(rhs.a, drop_axis);
        const auto rb = project_point_2d(rhs.b, drop_axis);
        const auto rc = project_point_2d(rhs.c, drop_axis);

        const std::array<std::pair<std::array<double, 2>, std::array<double, 2>>, 3> lhs_edges = {
          std::pair {la, lb},
          std::pair {lb, lc},
          std::pair {lc, la}
        };
        const std::array<std::pair<std::array<double, 2>, std::array<double, 2>>, 3> rhs_edges = {
          std::pair {ra, rb},
          std::pair {rb, rc},
          std::pair {rc, ra}
        };

        for(const auto &lhs_edge : lhs_edges) {
          for(const auto &rhs_edge : rhs_edges) {
            if(segments_intersect_2d_strict(
                 lhs_edge.first,
                 lhs_edge.second,
                 rhs_edge.first,
                 rhs_edge.second)) {
              return true;
            }
          }
        }

        return point_in_triangle_2d_strict(la, ra, rb, rc) ||
               point_in_triangle_2d_strict(ra, la, lb, lc);
      };

    auto segment_triangle_intersect_strict =
      [&](const Vec3 &p0,
          const Vec3 &p1,
          const Triangle3 &tri) noexcept {
        const Vec3 dir = vec3_sub(p1, p0);
        const Vec3 e1 = vec3_sub(tri.b, tri.a);
        const Vec3 e2 = vec3_sub(tri.c, tri.a);
        const Vec3 h = vec3_cross(dir, e2);
        const double det = vec3_dot(e1, h);
        if(std::abs(det) < kShellIntersectionEps) {
          return false;
        }

        const double inv_det = 1.0 / det;
        const Vec3 s = vec3_sub(p0, tri.a);
        const double u = inv_det * vec3_dot(s, h);
        if(u <= kShellIntersectionEps || u >= 1.0 - kShellIntersectionEps) {
          return false;
        }

        const Vec3 q = vec3_cross(s, e1);
        const double v = inv_det * vec3_dot(dir, q);
        if(v <= kShellIntersectionEps || u + v >= 1.0 - kShellIntersectionEps) {
          return false;
        }

        const double t = inv_det * vec3_dot(e2, q);
        return t > kShellIntersectionEps && t < 1.0 - kShellIntersectionEps;
      };

    auto triangles_intersect_strict =
      [&](const Triangle3 &lhs,
          const Triangle3 &rhs) noexcept {
        const Vec3 lhs_normal = vec3_cross(vec3_sub(lhs.b, lhs.a), vec3_sub(lhs.c, lhs.a));
        const Vec3 rhs_normal = vec3_cross(vec3_sub(rhs.b, rhs.a), vec3_sub(rhs.c, rhs.a));

        if(vec3_length(lhs_normal) < kShellIntersectionEps ||
           vec3_length(rhs_normal) < kShellIntersectionEps) {
          return false;
        }

        const std::array<double, 3> rhs_to_lhs = {
          vec3_dot(lhs_normal, vec3_sub(rhs.a, lhs.a)),
          vec3_dot(lhs_normal, vec3_sub(rhs.b, lhs.a)),
          vec3_dot(lhs_normal, vec3_sub(rhs.c, lhs.a))
        };
        const std::array<double, 3> lhs_to_rhs = {
          vec3_dot(rhs_normal, vec3_sub(lhs.a, rhs.a)),
          vec3_dot(rhs_normal, vec3_sub(lhs.b, rhs.a)),
          vec3_dot(rhs_normal, vec3_sub(lhs.c, rhs.a))
        };

        const auto same_side = [=](const std::array<double, 3> &distances) noexcept {
          bool has_pos = false;
          bool has_neg = false;
          for(double d : distances) {
            if(d > kShellIntersectionEps) has_pos = true;
            if(d < -kShellIntersectionEps) has_neg = true;
          }
          return !(has_pos && has_neg);
        };

        const bool coplanar =
          std::abs(rhs_to_lhs[0]) <= kShellIntersectionEps &&
          std::abs(rhs_to_lhs[1]) <= kShellIntersectionEps &&
          std::abs(rhs_to_lhs[2]) <= kShellIntersectionEps &&
          std::abs(lhs_to_rhs[0]) <= kShellIntersectionEps &&
          std::abs(lhs_to_rhs[1]) <= kShellIntersectionEps &&
          std::abs(lhs_to_rhs[2]) <= kShellIntersectionEps;
        if(coplanar) {
          return coplanar_triangles_intersect_strict(lhs, rhs, lhs_normal);
        }

        if(same_side(rhs_to_lhs) || same_side(lhs_to_rhs)) {
          return false;
        }

        const std::array<std::pair<Vec3, Vec3>, 3> lhs_edges = {
          std::pair {lhs.a, lhs.b},
          std::pair {lhs.b, lhs.c},
          std::pair {lhs.c, lhs.a}
        };
        const std::array<std::pair<Vec3, Vec3>, 3> rhs_edges = {
          std::pair {rhs.a, rhs.b},
          std::pair {rhs.b, rhs.c},
          std::pair {rhs.c, rhs.a}
        };

        for(const auto &edge : lhs_edges) {
          if(segment_triangle_intersect_strict(edge.first, edge.second, rhs)) {
            return true;
          }
        }
        for(const auto &edge : rhs_edges) {
          if(segment_triangle_intersect_strict(edge.first, edge.second, lhs)) {
            return true;
          }
        }

        return false;
      };

    std::vector<ValidationFace> validation_faces;
    validation_faces.reserve(work.surface_faces.size() * 8);

    auto add_validation_face =
      [&](const ValidationFace &face) {
        validation_faces.push_back(face);
        return static_cast<int>(validation_faces.size() - 1);
      };

    auto set_validation_face_output_ref =
      [&](int face_index, EntityRef face_ref) {
        if(face_index >= 0 &&
           face_index < static_cast<int>(validation_faces.size())) {
          validation_faces[static_cast<std::size_t>(face_index)].output_face_ref = face_ref;
        }
      };

    auto replace_validation_face_with_pyramid =
      [&](int source_face_index,
          const std::array<ShellTriangle, 4> &side_triangles,
          const std::array<EntityRef, 4> &side_face_refs,
          const std::array<bool, 4> &add_side_validation) {
        if(source_face_index >= 0 &&
           source_face_index < static_cast<int>(validation_faces.size())) {
          validation_faces[static_cast<std::size_t>(source_face_index)].active = false;
        }
        for(int side = 0; side < 4; ++side) {
          if(!add_side_validation[static_cast<std::size_t>(side)]) {
            continue;
          }
          const auto &side_triangle = side_triangles[static_cast<std::size_t>(side)];
          if(triangle_area(side_triangle.tri) <= kShellIntersectionEps) {
            continue;
          }
          const auto face_index =
            add_validation_face(make_validation_face_from_shell_triangle(side_triangle));
          set_validation_face_output_ref(
            face_index,
            side_face_refs[static_cast<std::size_t>(side)]);
        }
      };

    std::vector<std::uint32_t> wall_candidates;

    auto shell_triangle_hits_non_bl_wall =
      [&](const ShellTriangle &candidate) {
        if(work.non_bl_face_verts.empty()) {
          return false;
        }

        work.non_bl_bvh.query_overlap(candidate.box, wall_candidates);
        for(const auto wall_local : wall_candidates) {
          SurfaceFace wall_face {};
          wall_face.vertices = work.non_bl_face_verts[wall_local];
          wall_face.num_vertices = work.non_bl_face_nv[wall_local];
          for(int tri_index = 0; tri_index < face_triangle_count(wall_face); ++tri_index) {
            const auto wall_tri = face_triangle(work.non_bl_nodes, wall_face, tri_index);
            if(triangles_intersect_strict(candidate.tri, wall_tri)) {
              return true;
            }
          }
          for(int edge_index = 0; edge_index < wall_face.num_vertices; ++edge_index) {
            const auto wall_a = work.non_bl_nodes[wall_face.vertices[static_cast<std::size_t>(edge_index)]];
            const auto wall_b =
              work.non_bl_nodes[wall_face.vertices[static_cast<std::size_t>((edge_index + 1) % wall_face.num_vertices)]];
            const auto wall_edge = make_candidate_edge(
              kInvalidVertexKey,
              wall_a,
              kInvalidVertexKey,
              wall_b);
            if(!candidate.box.overlaps(wall_edge.box)) {
              continue;
            }
            if(segment_triangle_intersect_strict(wall_edge.a, wall_edge.b, candidate.tri)) {
              return true;
            }
          }
        }

        return false;
      };

    auto shell_edge_hits_non_bl_wall =
      [&](const ShellEdge &candidate) {
        if(work.non_bl_face_verts.empty()) {
          return false;
        }

        work.non_bl_bvh.query_overlap(candidate.box, wall_candidates);
        for(const auto wall_local : wall_candidates) {
          SurfaceFace wall_face {};
          wall_face.vertices = work.non_bl_face_verts[wall_local];
          wall_face.num_vertices = work.non_bl_face_nv[wall_local];
          for(int tri_index = 0; tri_index < face_triangle_count(wall_face); ++tri_index) {
            const auto wall_tri = face_triangle(work.non_bl_nodes, wall_face, tri_index);
            if(segment_triangle_intersect_strict(candidate.a, candidate.b, wall_tri)) {
              return true;
            }
          }

          const auto &wall_face_indices = work.non_bl_face_verts[wall_local];
          const int wall_nv = work.non_bl_face_nv[wall_local];
          for(int edge_index = 0; edge_index < wall_nv; ++edge_index) {
            const auto wall_a = work.non_bl_nodes[wall_face_indices[static_cast<std::size_t>(edge_index)]];
            const auto wall_b = work.non_bl_nodes[wall_face_indices[static_cast<std::size_t>((edge_index + 1) % wall_nv)]];
            const auto wall_edge = make_candidate_edge(
              kInvalidVertexKey,
              wall_a,
              kInvalidVertexKey,
              wall_b);
            if(!candidate.box.overlaps(wall_edge.box)) {
              continue;
            }
            if(edges_have_nonconforming_overlap(candidate, wall_edge)) {
              return true;
            }
          }
        }

        return false;
      };

    auto shell_triangle_hits_shell =
      [&](int excluded_face_index,
          const ShellTriangle &candidate) {
        for(int face_index = 0; face_index < static_cast<int>(validation_faces.size()); ++face_index) {
          if(face_index == excluded_face_index) {
            continue;
          }
          const auto &other_face = validation_faces[static_cast<std::size_t>(face_index)];
          if(!other_face.active || !candidate.box.overlaps(other_face.box)) {
            continue;
          }
          const bool shares_boundary_edge =
            face_shares_candidate_base_edge(other_face, candidate);
          for(int tri_index = 0; tri_index < validation_face_triangle_count(other_face); ++tri_index) {
            if(shares_boundary_edge && other_face.num_vertices == 4) {
              continue;
            }
            const auto other_triangle = validation_face_triangle(other_face, tri_index);
            const auto other_triangle_keys =
              validation_face_triangle_keys(other_face, tri_index);
            const int shared_vertex_count =
              count_shared_triangle_vertices(candidate.vertex_keys, other_triangle_keys);
            if(shared_vertex_count >= 3) {
              // An exact reused-face match must be handled by face reuse, not by
              // silently emitting another shell triangle with the same three nodes.
              return true;
            }
            if(shared_vertex_count == 2) {
              // Legitimate adjacent shell triangles often share one edge and only
              // meet on that edge. However, if the two triangles also overlap away
              // from the shared edge (typically because a pyramid was emitted on the
              // wrong side or the same side triangle was created twice), we must
              // reject it here instead of skipping the check entirely.
              if(triangles_intersect_strict(candidate.tri, other_triangle)) {
                return true;
              }
              continue;
            }
            if(triangles_intersect_strict(candidate.tri, other_triangle)) {
              return true;
            }
          }
          for(int edge_index = 0; edge_index < other_face.num_vertices; ++edge_index) {
            const auto other_edge = validation_face_edge(other_face, edge_index);
            if(!candidate.box.overlaps(other_edge.box)) {
              continue;
            }
            if(segment_triangle_intersect_strict(other_edge.a, other_edge.b, candidate.tri)) {
              return true;
            }
          }
        }
        return false;
      };

    auto face_index_in_extra_exclusions =
      [&](int face_index, const std::array<int, 4> &extra_excluded_faces) {
        for(const int excluded : extra_excluded_faces) {
          if(face_index == excluded) {
            return true;
          }
        }
        return false;
      };

    auto shell_triangle_hits_shell_excluding =
      [&](int excluded_face_index,
          const std::array<int, 4> &extra_excluded_faces,
          const ShellTriangle &candidate) {
        for(int face_index = 0; face_index < static_cast<int>(validation_faces.size()); ++face_index) {
          if(face_index == excluded_face_index ||
             face_index_in_extra_exclusions(face_index, extra_excluded_faces)) {
            continue;
          }
          const auto &other_face = validation_faces[static_cast<std::size_t>(face_index)];
          if(!other_face.active || !candidate.box.overlaps(other_face.box)) {
            continue;
          }
          const bool shares_boundary_edge =
            face_shares_candidate_base_edge(other_face, candidate);
          for(int tri_index = 0; tri_index < validation_face_triangle_count(other_face); ++tri_index) {
            if(shares_boundary_edge && other_face.num_vertices == 4) {
              continue;
            }
            const auto other_triangle = validation_face_triangle(other_face, tri_index);
            const auto other_triangle_keys =
              validation_face_triangle_keys(other_face, tri_index);
            const int shared_vertex_count =
              count_shared_triangle_vertices(candidate.vertex_keys, other_triangle_keys);
            if(shared_vertex_count >= 3) {
              return true;
            }
            if(shared_vertex_count == 2) {
              if(triangles_intersect_strict(candidate.tri, other_triangle)) {
                return true;
              }
              continue;
            }
            if(triangles_intersect_strict(candidate.tri, other_triangle)) {
              return true;
            }
          }
          for(int edge_index = 0; edge_index < other_face.num_vertices; ++edge_index) {
            const auto other_edge = validation_face_edge(other_face, edge_index);
            if(!candidate.box.overlaps(other_edge.box)) {
              continue;
            }
            if(segment_triangle_intersect_strict(other_edge.a, other_edge.b, candidate.tri)) {
              return true;
            }
          }
        }
        return false;
      };

    auto shell_edge_hits_shell =
      [&](int excluded_face_index,
          const ShellEdge &candidate) {
        for(int face_index = 0; face_index < static_cast<int>(validation_faces.size()); ++face_index) {
          if(face_index == excluded_face_index) {
            continue;
          }
          const auto &other_face = validation_faces[static_cast<std::size_t>(face_index)];
          if(!other_face.active || !candidate.box.overlaps(other_face.box)) {
            continue;
          }
          for(int tri_index = 0; tri_index < validation_face_triangle_count(other_face); ++tri_index) {
            const auto other_triangle = validation_face_triangle(other_face, tri_index);
            if(segment_triangle_intersect_strict(candidate.a, candidate.b, other_triangle)) {
              return true;
            }
          }
          for(int edge_index = 0; edge_index < other_face.num_vertices; ++edge_index) {
            const auto other_edge = validation_face_edge(other_face, edge_index);
            if(!candidate.box.overlaps(other_edge.box)) {
              continue;
            }
            if(edges_have_nonconforming_overlap(candidate, other_edge)) {
              return true;
            }
          }
        }
        return false;
      };

    auto shell_edge_hits_shell_excluding =
      [&](int excluded_face_index,
          const std::array<int, 4> &extra_excluded_faces,
          const ShellEdge &candidate) {
        for(int face_index = 0; face_index < static_cast<int>(validation_faces.size()); ++face_index) {
          if(face_index == excluded_face_index ||
             face_index_in_extra_exclusions(face_index, extra_excluded_faces)) {
            continue;
          }
          const auto &other_face = validation_faces[static_cast<std::size_t>(face_index)];
          if(!other_face.active || !candidate.box.overlaps(other_face.box)) {
            continue;
          }
          for(int tri_index = 0; tri_index < validation_face_triangle_count(other_face); ++tri_index) {
            const auto other_triangle = validation_face_triangle(other_face, tri_index);
            if(segment_triangle_intersect_strict(candidate.a, candidate.b, other_triangle)) {
              return true;
            }
          }
          for(int edge_index = 0; edge_index < other_face.num_vertices; ++edge_index) {
            const auto other_edge = validation_face_edge(other_face, edge_index);
            if(!candidate.box.overlaps(other_edge.box)) {
              continue;
            }
            if(edges_have_nonconforming_overlap(candidate, other_edge)) {
              return true;
            }
          }
        }
        return false;
      };

    auto face_anchor_center = [&](std::uint32_t fi, std::size_t layer) -> Vec3 {
      const auto &face = work.surface_faces[fi];
      Vec3 sum {0.0, 0.0, 0.0};
      double count = 0.0;
      for(int vi = 0; vi < face.num_vertices; ++vi) {
        const auto bot = node_position(node_refs[layer][face.vertices[vi]]);
        const auto top = node_position(node_refs[layer + 1][face.vertices[vi]]);
        sum[0] += bot[0] + top[0];
        sum[1] += bot[1] + top[1];
        sum[2] += bot[2] + top[2];
        count += 2.0;
      }
      return {sum[0] / count, sum[1] / count, sum[2] / count};
    };

    auto collapsed_face_inside = [&](std::uint32_t fi) -> Vec3 {
      const auto &face = work.surface_faces[fi];
      std::array<EntityRef, 4> refs {{}};
      for(int vi = 0; vi < face.num_vertices; ++vi) {
        refs[static_cast<std::size_t>(vi)] = node_refs[0][face.vertices[vi]];
      }
      const auto center = face.num_vertices == 3
        ? tri_center({refs[0], refs[1], refs[2]})
        : quad_center(refs);
      Vec3 normal {0.0, 0.0, 0.0};
      if(fi < work.face_normal.size()) {
        normal = vec3_normalized(work.face_normal[fi]);
      }
      if(vec3_length(normal) < 1e-12) {
        normal = face.num_vertices == 3
          ? tri_normal({refs[0], refs[1], refs[2]})
          : quad_normal(refs);
      }
      return vec3_sub(center, vec3_scale(normal, 0.25 * polygon_scale(refs, face.num_vertices)));
    };

    struct EdgeUse {
      std::uint32_t face = 0;
      std::uint32_t v0 = 0;
      std::uint32_t v1 = 0;
    };

    std::map<std::pair<std::uint32_t, std::uint32_t>, std::vector<EdgeUse>> edge_uses;
    for(std::uint32_t fi = 0; fi < work.surface_faces.size(); ++fi) {
      const auto &face = work.surface_faces[fi];
      for(int ei = 0; ei < face.num_vertices; ++ei) {
        const auto v0 = face.vertices[static_cast<std::size_t>(ei)];
        const auto v1 = face.vertices[static_cast<std::size_t>((ei + 1) % face.num_vertices)];
        const auto key = std::make_pair(std::min(v0, v1), std::max(v0, v1));
        edge_uses[key].push_back({fi, v0, v1});
      }
    }

    std::vector<int> top_validation_face_indices(work.surface_faces.size(), -1);

    struct PendingTopQuadFace {
      std::array<EntityRef, 4> quad {};
      Vec3 inside_point {0.0, 0.0, 0.0};
      geo::TopologyEntityId topo {};
      int validation_face_index = -1;
    };

    struct PendingTransitionFace {
      std::array<EntityRef, 4> quad {};
      Vec3 inside_point {0.0, 0.0, 0.0};
      geo::TopologyEntityId topo {};
      int validation_face_index = -1;
      int preferred_shared_edge_index = -1;
    };
    std::vector<PendingTopQuadFace> pending_top_quad_faces;
    pending_top_quad_faces.reserve(work.surface_faces.size());
    std::vector<PendingTransitionFace> pending_transition_faces;
    pending_transition_faces.reserve(edge_uses.size() * 2);

    for(std::uint32_t fi = 0; fi < work.surface_faces.size(); ++fi) {
      const auto &face = work.surface_faces[fi];
      const auto actual_layer = static_cast<std::size_t>(face_levels[fi]);
      const auto inside_point = actual_layer > 0
        ? face_anchor_center(fi, actual_layer - 1)
        : collapsed_face_inside(fi);
      const auto &v = face.vertices;
      if(face.num_vertices == 3) {
        std::array<EntityRef, 3> tri = {
          node_refs[actual_layer][v[0]],
          node_refs[actual_layer][v[1]],
          node_refs[actual_layer][v[2]]
        };
        tri = orient_tri_outward(tri, inside_point);
        std::array<EntityRef, 4> tri_refs {tri[0], tri[1], tri[2], {}};
        top_validation_face_indices[fi] = add_validation_face(make_validation_face_from_refs(tri_refs, 3));
      } else if(face.num_vertices == 4) {
        std::array<EntityRef, 4> quad = {
          node_refs[actual_layer][v[0]],
          node_refs[actual_layer][v[1]],
          node_refs[actual_layer][v[2]],
          node_refs[actual_layer][v[3]]
        };
        quad = orient_quad_outward(quad, inside_point);
        top_validation_face_indices[fi] = add_validation_face(make_validation_face_from_refs(quad, 4));
      }
    }

    for(const auto &entry : edge_uses) {
      const auto &uses = entry.second;
      if(uses.empty()) continue;

      int min_level = std::numeric_limits<int>::max();
      for(const auto &use : uses) {
        min_level = std::min(min_level, face_levels[use.face]);
      }
      if(min_level == std::numeric_limits<int>::max()) {
        min_level = 0;
      }

      for(const auto &use : uses) {
        const int face_level = face_levels[use.face];
        const int begin_layer = (uses.size() == 1) ? 0 : min_level;
        if(face_level <= begin_layer) continue;

        const auto topo = face_topology(use.face);
        for(int layer = begin_layer; layer < face_level; ++layer) {
          const auto top_edge_a =
            node_refs[static_cast<std::size_t>(layer + 1)][use.v1];
          const auto top_edge_b =
            node_refs[static_cast<std::size_t>(layer + 1)][use.v0];
          std::array<EntityRef, 4> quad = {
            node_refs[static_cast<std::size_t>(layer)][use.v0],
            node_refs[static_cast<std::size_t>(layer)][use.v1],
            node_refs[static_cast<std::size_t>(layer + 1)][use.v1],
            node_refs[static_cast<std::size_t>(layer + 1)][use.v0]
          };
          const auto inside_point = face_anchor_center(use.face, static_cast<std::size_t>(layer));
          quad = orient_quad_outward(quad, inside_point);
          int preferred_shared_edge_index = -1;
          for(int edge_index = 0; edge_index < 4; ++edge_index) {
            const auto edge_a = quad[static_cast<std::size_t>(edge_index)];
            const auto edge_b = quad[static_cast<std::size_t>((edge_index + 1) % 4)];
            const bool same_edge =
              (edge_a == top_edge_a && edge_b == top_edge_b) ||
              (edge_a == top_edge_b && edge_b == top_edge_a);
            if(same_edge) {
              preferred_shared_edge_index = edge_index;
              break;
            }
          }
          pending_transition_faces.push_back({
            quad,
            inside_point,
            topo,
            add_validation_face(make_validation_face_from_refs(quad, 4)),
            preferred_shared_edge_index
          });
        }
      }
    }

    std::vector<int> preferred_shared_edge_by_validation_face(
      validation_faces.size(),
      -1);
    for(const auto &pending_transition_face : pending_transition_faces) {
      if(pending_transition_face.validation_face_index < 0 ||
         pending_transition_face.validation_face_index >=
           static_cast<int>(preferred_shared_edge_by_validation_face.size())) {
        continue;
      }
      preferred_shared_edge_by_validation_face[
        static_cast<std::size_t>(pending_transition_face.validation_face_index)] =
        pending_transition_face.preferred_shared_edge_index;
    }

    if constexpr (kEnableBoundaryLayerDebugFileOutput) {
      auto export_pre_pyramid_closed_shell_obj =
        [&](const std::string &path) {
        auto pack_entity = [](EntityRef ref) -> std::uint64_t {
          return (static_cast<std::uint64_t>(ref.entity_group) << 32U) |
                 static_cast<std::uint64_t>(ref.index);
        };
        auto pack_non_bl_node = [](std::uint32_t node_index) -> std::uint64_t {
          return (std::uint64_t {1} << 63) | static_cast<std::uint64_t>(node_index);
        };

        struct ObjFace {
          int nv = 0;
          std::array<std::size_t, 4> v {0, 0, 0, 0};
        };

        std::map<std::uint64_t, std::size_t> node_map;
        std::vector<Vec3> points;
        std::vector<ObjFace> faces;
        faces.reserve(work.surface_faces.size() * 2 + pending_transition_faces.size());

        auto get_node = [&](EntityRef ref) -> std::size_t {
          const auto key = pack_entity(ref);
          const auto it = node_map.find(key);
          if(it != node_map.end()) {
            return it->second;
          }
          const auto idx = points.size();
          points.push_back(node_position(ref));
          node_map.emplace(key, idx);
          return idx;
        };

        auto append_face = [&](const std::array<EntityRef, 4> &refs, int nv) {
          ObjFace face;
          face.nv = nv;
          for(int i = 0; i < nv; ++i) {
            face.v[static_cast<std::size_t>(i)] = get_node(refs[static_cast<std::size_t>(i)]);
          }
          faces.push_back(face);
        };

        // Export the cavity shell seen by the pyramid step:
        //   1. BL top cap faces at their reached layer
        //   2. transition quads between mismatched reached layers
        //   3. preserved non-BL boundary faces
        // Do not include the BL bottom wall here, otherwise the OBJ looks like
        // the whole BL volume shell instead of the surface handed to the
        // top/transition pyramid conversion + downstream Tet closure.
        for(std::uint32_t fi = 0; fi < work.surface_faces.size(); ++fi) {
          const auto &face = work.surface_faces[fi];
          const auto actual_layer = static_cast<std::size_t>(face_levels[fi]);
          const auto inside_point = actual_layer > 0
            ? face_anchor_center(fi, actual_layer - 1)
            : collapsed_face_inside(fi);
          const auto &v = face.vertices;

          if(face.num_vertices == 3) {
            std::array<EntityRef, 3> top_tri = {
              node_refs[actual_layer][v[0]],
              node_refs[actual_layer][v[1]],
              node_refs[actual_layer][v[2]]
            };
            top_tri = orient_tri_outward(top_tri, inside_point);
            append_face({top_tri[0], top_tri[1], top_tri[2], {}}, 3);
          } else if(face.num_vertices == 4) {
            std::array<EntityRef, 4> top_quad = {
              node_refs[actual_layer][v[0]],
              node_refs[actual_layer][v[1]],
              node_refs[actual_layer][v[2]],
              node_refs[actual_layer][v[3]]
            };
            top_quad = orient_quad_outward(top_quad, inside_point);
            append_face(top_quad, 4);
          }
        }

        for(const auto &transition_face : pending_transition_faces) {
          append_face(transition_face.quad, 4);
        }

        for(std::size_t face_index = 0; face_index < work.non_bl_face_verts.size(); ++face_index) {
          const auto &verts = work.non_bl_face_verts[face_index];
          const int nv = work.non_bl_face_nv[face_index];
          ObjFace face;
          face.nv = nv;
          for(int i = 0; i < nv; ++i) {
            const auto node_index = verts[static_cast<std::size_t>(i)];
            const auto key = pack_non_bl_node(node_index);
            const auto it = node_map.find(key);
            if(it != node_map.end()) {
              face.v[static_cast<std::size_t>(i)] = it->second;
              continue;
            }
            const auto idx = points.size();
            points.push_back(work.non_bl_nodes[node_index]);
            node_map.emplace(key, idx);
            face.v[static_cast<std::size_t>(i)] = idx;
          }
          faces.push_back(face);
        }

        std::ofstream out(path, std::ios::trunc);
        if(!out) {
          return;
        }

        out << "# Closed cavity shell before pyramid conversion\n";
        out << "# " << points.size() << " vertices, " << faces.size() << " faces\n";
        out << std::setprecision(15);
        for(const auto &p : points) {
          out << "v " << p[0] << " " << p[1] << " " << p[2] << "\n";
        }
        for(const auto &face : faces) {
          out << "f";
          for(int i = 0; i < face.nv; ++i) {
            out << " " << (face.v[static_cast<std::size_t>(i)] + 1);
          }
          out << "\n";
        }

        };

      export_pre_pyramid_closed_shell_obj("bl_closed_shell_pre_pyramid.obj");
    }

    struct PrecomputedQuadApexCandidate {
      EntityRef apex_ref {};
      std::array<int, 4> matched_face_indices {-1, -1, -1, -1};
    };

    struct PrecomputedQuadApexInfo {
      std::vector<PrecomputedQuadApexCandidate> candidates;
    };

    struct SharedEdgeQuadTriPair {
      int quad_face_index = -1;
      int tri_face_index = -1;
      int shared_edge_index = -1;
      double dihedral_degrees = 360.0;
    };

    struct AssignedSharedEdgePyramidPlan {
      bool assigned = false;
      EntityRef apex_ref {};
      std::array<int, 4> matched_face_indices {-1, -1, -1, -1};
      int selected_shared_edge_index = -1;
      double selected_dihedral_degrees = 360.0;
    };

    const auto pre_pyramid_validation_faces = validation_faces;
    std::vector<PrecomputedQuadApexInfo> precomputed_quad_apex_info(
      validation_faces.size());
    std::vector<SharedEdgeQuadTriPair> precomputed_shared_edge_pairs;
    {
      for(int face_index = 0;
          face_index < static_cast<int>(pre_pyramid_validation_faces.size());
          ++face_index) {
        const auto &quad_face =
          pre_pyramid_validation_faces[static_cast<std::size_t>(face_index)];
        if(!quad_face.active || quad_face.num_vertices != 4) {
          continue;
        }

        std::map<std::uint64_t, PrecomputedQuadApexCandidate> apex_candidates;
        for(int other_face_index = 0;
            other_face_index < static_cast<int>(pre_pyramid_validation_faces.size());
            ++other_face_index) {
          if(other_face_index == face_index) {
            continue;
          }

          const auto &other_face =
            pre_pyramid_validation_faces[static_cast<std::size_t>(other_face_index)];
          if(!other_face.active || other_face.num_vertices != 3) {
            continue;
          }

          const auto tri_keys = validation_face_triangle_keys(other_face, 0);
          int matched_edge_index = -1;
          std::uint64_t apex_key = kInvalidVertexKey;
          for(int edge_index = 0; edge_index < 4; ++edge_index) {
            const auto a =
              quad_face.vertex_keys[static_cast<std::size_t>(edge_index)];
            const auto b =
              quad_face.vertex_keys[static_cast<std::size_t>((edge_index + 1) % 4)];
            const bool has_a =
              tri_keys[0] == a || tri_keys[1] == a || tri_keys[2] == a;
            const bool has_b =
              tri_keys[0] == b || tri_keys[1] == b || tri_keys[2] == b;
            if(!has_a || !has_b) {
              continue;
            }

            matched_edge_index = edge_index;
            for(const auto tri_key : tri_keys) {
              if(tri_key != a && tri_key != b) {
                apex_key = tri_key;
                break;
              }
            }
            break;
          }

          if(matched_edge_index < 0 || apex_key == kInvalidVertexKey) {
            continue;
          }

          precomputed_shared_edge_pairs.push_back({
            face_index,
            other_face_index,
            matched_edge_index,
            shared_edge_dihedral_degrees(quad_face, matched_edge_index, other_face)
          });

          auto &candidate = apex_candidates[apex_key];
          candidate.apex_ref = unpack_ref(apex_key);
          candidate.matched_face_indices[static_cast<std::size_t>(matched_edge_index)] =
            other_face_index;
        }

        auto &info =
          precomputed_quad_apex_info[static_cast<std::size_t>(face_index)];
        info.candidates.reserve(apex_candidates.size());
        for(const auto &[apex_key, candidate] : apex_candidates) {
          static_cast<void>(apex_key);
          info.candidates.push_back(candidate);
        }
      }
    }

    std::vector<AssignedSharedEdgePyramidPlan> assigned_shared_edge_pyramid_plans(
      validation_faces.size());
    std::vector<int> claimed_shared_edge_triangle_owner(
      validation_faces.size(),
      -1);
    std::size_t assigned_shared_edge_pyramid_count = 0;
    std::size_t claimed_shared_edge_triangle_count = 0;
    {
      auto candidate_is_unclaimed =
        [&](const PrecomputedQuadApexCandidate &candidate) {
          for(const int matched_face_index : candidate.matched_face_indices) {
            if(matched_face_index < 0) {
              continue;
            }
            if(matched_face_index >=
               static_cast<int>(claimed_shared_edge_triangle_owner.size())) {
              return false;
            }
            if(claimed_shared_edge_triangle_owner[
                 static_cast<std::size_t>(matched_face_index)] >= 0) {
              return false;
            }
          }
          return true;
        };

      for(int face_index = 0;
          face_index < static_cast<int>(pre_pyramid_validation_faces.size());
          ++face_index) {
        const auto &quad_face =
          pre_pyramid_validation_faces[static_cast<std::size_t>(face_index)];
        if(!quad_face.active || quad_face.num_vertices != 4) {
          continue;
        }

        if(face_index >= static_cast<int>(precomputed_quad_apex_info.size())) {
          continue;
        }

        const auto &apex_info =
          precomputed_quad_apex_info[static_cast<std::size_t>(face_index)];
        if(apex_info.candidates.empty()) {
          continue;
        }

        std::array<int, 4> edge_order {0, 1, 2, 3};
        const int preferred_edge =
          face_index < static_cast<int>(preferred_shared_edge_by_validation_face.size())
            ? preferred_shared_edge_by_validation_face[static_cast<std::size_t>(face_index)]
            : -1;
        if(preferred_edge >= 0 && preferred_edge < 4) {
          edge_order = {preferred_edge, (preferred_edge + 1) % 4,
                        (preferred_edge + 2) % 4, (preferred_edge + 3) % 4};
        }

        bool assigned = false;
        for(const int edge_index : edge_order) {
          for(const auto &candidate : apex_info.candidates) {
            const int matched_face_index =
              candidate.matched_face_indices[static_cast<std::size_t>(edge_index)];
            if(matched_face_index < 0 ||
               matched_face_index >=
                 static_cast<int>(pre_pyramid_validation_faces.size())) {
              continue;
            }
            if(!is_valid(candidate.apex_ref) || !candidate_is_unclaimed(candidate)) {
              continue;
            }

            const auto &tri_face =
              pre_pyramid_validation_faces[static_cast<std::size_t>(matched_face_index)];
            const double dihedral_degrees =
              shared_edge_dihedral_degrees(quad_face, edge_index, tri_face);
            if(!(dihedral_degrees < kMaximumSharedEdgeDihedralDegrees)) {
              continue;
            }

            auto &plan =
              assigned_shared_edge_pyramid_plans[static_cast<std::size_t>(face_index)];
            plan.assigned = true;
            plan.apex_ref = candidate.apex_ref;
            plan.matched_face_indices = candidate.matched_face_indices;
            plan.selected_shared_edge_index = edge_index;
            plan.selected_dihedral_degrees = dihedral_degrees;
            ++assigned_shared_edge_pyramid_count;

            for(const int claimed_face_index : candidate.matched_face_indices) {
              if(claimed_face_index < 0) {
                continue;
              }
              claimed_shared_edge_triangle_owner[
                static_cast<std::size_t>(claimed_face_index)] = face_index;
              ++claimed_shared_edge_triangle_count;
            }
            assigned = true;
            break;
          }
          if(assigned) {
            break;
          }
        }
      }
    }

    if constexpr (kEnableBoundaryLayerDebugFileOutput) {
      auto export_shared_edge_quad_tri_pairs_obj =
        [&](const std::string &path,
            bool export_only_over_limit) {
        struct ObjRecord {
          std::string comment;
          int nv = 0;
          std::array<std::size_t, 4> v {0, 0, 0, 0};
        };

        std::vector<Vec3> points;
        std::vector<ObjRecord> records;
        points.reserve(precomputed_shared_edge_pairs.size() * 7);
        records.reserve(precomputed_shared_edge_pairs.size() * 2);

        auto append_face =
          [&](const ValidationFace &face,
              int nv,
              const std::string &comment) {
            ObjRecord record;
            record.comment = comment;
            record.nv = nv;
            for(int i = 0; i < nv; ++i) {
              record.v[static_cast<std::size_t>(i)] = points.size();
              points.push_back(face.points[static_cast<std::size_t>(i)]);
            }
            records.push_back(std::move(record));
          };

        std::size_t exported_pair_count = 0;
        for(const auto &pair : precomputed_shared_edge_pairs) {
          const bool is_over_limit =
            pair.dihedral_degrees > kMaximumSharedEdgeDihedralDegrees;
          if(export_only_over_limit && !is_over_limit) {
            continue;
          }
          if(!export_only_over_limit && pair.quad_face_index < 0) {
            continue;
          }

          const auto &quad_face =
            pre_pyramid_validation_faces[static_cast<std::size_t>(pair.quad_face_index)];
          const auto &tri_face =
            pre_pyramid_validation_faces[static_cast<std::size_t>(pair.tri_face_index)];

          std::ostringstream pair_header;
          pair_header << std::setprecision(15)
                      << "# quad_validation_face_index=" << pair.quad_face_index
                      << " tri_validation_face_index=" << pair.tri_face_index
                      << " shared_edge_index=" << pair.shared_edge_index
                      << " dihedral_degrees=" << pair.dihedral_degrees;
          const auto pair_header_text = pair_header.str();
          append_face(quad_face, 4, pair_header_text + " role=quad");
          append_face(tri_face, 3, pair_header_text + " role=tri");
          ++exported_pair_count;
        }

        std::ofstream out(path, std::ios::trunc);
        if(!out) {
          return;
        }

        out << "# Shared-edge quad-tri pairs before pyramid conversion\n";
        out << "# " << points.size() << " vertices, " << records.size()
            << " faces, " << exported_pair_count << " pairs\n";
        if(export_only_over_limit) {
          out << "# Only pairs with dihedral greater than "
              << kMaximumSharedEdgeDihedralDegrees << " degrees are included\n";
        } else {
          out << "# Includes all quad-tri pairs sharing an edge in the pre-pyramid shell\n";
        }
        out << std::setprecision(15);
        for(const auto &p : points) {
          out << "v " << p[0] << " " << p[1] << " " << p[2] << "\n";
        }
        for(const auto &record : records) {
          if(!record.comment.empty()) {
            out << record.comment << "\n";
          }
          out << "f";
          for(int i = 0; i < record.nv; ++i) {
            out << " " << (record.v[static_cast<std::size_t>(i)] + 1);
          }
          out << "\n";
        }

        };

      export_shared_edge_quad_tri_pairs_obj(
        "bl_shared_edge_quad_tri_pairs.obj",
        false);
      export_shared_edge_quad_tri_pairs_obj(
        "bl_shared_edge_quad_tri_pairs_over_150deg.obj",
        true);
    }

    std::size_t pyramid_count = 0;
    std::size_t direct_top_quad_count = 0;
    std::size_t direct_transition_quad_count = 0;
    std::size_t top_tri_count = 0;
    std::size_t transition_target_side_override_count = 0;
    std::size_t transition_target_same_side_count = 0;
    std::size_t rejected_pyramid_quality_count = 0;
    std::size_t rejected_pyramid_shell_intersection_count = 0;
    std::size_t rejected_pyramid_wall_intersection_count = 0;
    std::size_t preferred_edge_candidate_found_count = 0;
    std::size_t preferred_edge_candidate_missing_count = 0;
    std::size_t preferred_edge_candidate_success_count = 0;
    std::size_t preferred_edge_candidate_quality_count = 0;
    std::size_t preferred_edge_candidate_shell_collision_count = 0;
    std::size_t preferred_edge_candidate_wall_collision_count = 0;
    std::size_t preferred_edge_no_candidate_on_edge_count = 0;
    std::size_t preferred_edge_multiple_candidates_on_edge_count = 0;
    std::size_t preferred_edge_global_fallback_reject_count = 0;
    std::size_t preferred_edge_invalid_apex_ref_count = 0;
    std::size_t preferred_edge_dihedral_reject_count = 0;
    std::size_t preferred_edge_missing_output_face_ref_count = 0;
    std::size_t preferred_edge_negative_volume_reject_count = 0;

    enum class QuadEmitOutcome {
      pyramid,
      quality_or_height,
      shell_collision,
      wall_collision
    };

    auto record_rejected_pyramid =
      [&](QuadEmitOutcome outcome) {
        switch(outcome) {
        case QuadEmitOutcome::quality_or_height:
          ++rejected_pyramid_quality_count;
          break;
        case QuadEmitOutcome::shell_collision:
          ++rejected_pyramid_shell_intersection_count;
          break;
        case QuadEmitOutcome::wall_collision:
          ++rejected_pyramid_wall_intersection_count;
          break;
        case QuadEmitOutcome::pyramid:
          break;
        }
      };

    auto emit_boundary_quad =
      [&](const std::array<EntityRef, 4> &quad,
          geo::TopologyEntityId topo,
          EntityGroupIndex quad_entity_group,
          std::size_t &quad_count) {
        auto quad_ref = output.add_quad_face(quad_entity_group, quad);
        if(geo::is_valid(topo)) {
          output.set_face_topology_owner(quad_ref, topo);
        }
        ++quad_count;
        return quad_ref;
      };

    auto emit_pyramid =
      [&](std::array<EntityRef, 4> boundary_quad,
          const Vec3 &inside_point,
          geo::TopologyEntityId topo,
          int validation_face_index,
          bool prefer_target_point_direction,
          bool allow_preassigned_shared_edge_candidate,
          bool allow_inserted_apex_candidate,
          int preferred_shared_edge_index = -1) -> QuadEmitOutcome {
        boundary_quad = orient_quad_outward(boundary_quad, inside_point);

        const auto outward = quad_normal(boundary_quad);
        const auto base_center = quad_center(boundary_quad);
        const double scale = polygon_scale(boundary_quad, 4);
        auto pyramid_direction = outward;
        if(prefer_target_point_direction && has_target_point_) {
          const auto to_target = vec3_sub(target_point, base_center);
          if(vec3_length(to_target) > 1e-12) {
            const auto target_direction = vec3_normalized(to_target);
            const double target_alignment = vec3_dot(target_direction, outward);
            // Transition pyramids stay purely on the prism-exterior normal.
            // We keep the target-point relation only as diagnostics so we can
            // see how often the global target would have pointed to the wrong side.
            if(target_alignment > 1e-8) {
              ++transition_target_same_side_count;
            } else {
              ++transition_target_side_override_count;
            }
          }
        }

        double height = 0.35 * scale;
        if(vec3_length(pyramid_direction) > 1e-12) {
          const double query_dist = std::max(height * 6.0, scale * 2.0);
          const double non_bl_hit =
            work.non_bl_bvh.ray_nearest_intersection(base_center, pyramid_direction, query_dist, {});
          if(non_bl_hit > 0.0) {
            height = std::min(height, 0.45 * non_bl_hit);
          }
          const double surface_hit =
            work.face_bvh.ray_nearest_intersection(base_center, pyramid_direction, query_dist, {});
          if(surface_hit > 0.0) {
            height = std::min(height, 0.45 * surface_hit);
          }
        }

        if(height <= 1e-7 || vec3_length(pyramid_direction) <= 1e-12) {
          return QuadEmitOutcome::quality_or_height;
        }

        const double geom_scale = std::max(scale, 1e-6);
        const double min_height = std::max(1e-12, 1e-8 * geom_scale);
        const double min_area = std::max(1e-18, 1e-10 * geom_scale * geom_scale);
        const auto boundary_quad_keys = std::array<std::uint64_t, 4> {
          pack_ref(boundary_quad[0]),
          pack_ref(boundary_quad[1]),
          pack_ref(boundary_quad[2]),
          pack_ref(boundary_quad[3])
        };
        const auto base_positions = std::array<Vec3, 4> {
          node_position(boundary_quad[0]),
          node_position(boundary_quad[1]),
          node_position(boundary_quad[2]),
          node_position(boundary_quad[3])
        };
        const auto base_tri0 = Triangle3 {
          base_positions[0],
          base_positions[1],
          base_positions[2]
        };
        const auto base_tri1 = Triangle3 {
          base_positions[0],
          base_positions[2],
          base_positions[3]
        };
        const double base_area = triangle_area(base_tri0) + triangle_area(base_tri1);
        if(base_area <= min_area) {
          return QuadEmitOutcome::quality_or_height;
        }

        constexpr int kMaxHeightShrinkAttempts = 12;
        constexpr double kHeightShrinkFactor = 0.5;
        bool saw_shell_collision = false;
        bool saw_wall_collision = false;
        bool accepted_pyramid = false;
        EntityRef accepted_apex_ref {};
        Vec3 accepted_apex {0.0, 0.0, 0.0};
        std::array<EntityRef, 4> quad_base = boundary_quad;
        std::array<EntityRef, 4> reused_side_face_refs {};
        std::array<bool, 4> reuse_existing_side_face = {false, false, false, false};
        std::array<ShellTriangle, 4> candidate_side_tris = {
          make_candidate_triangle(boundary_quad[0], boundary_quad[1], base_center),
          make_candidate_triangle(boundary_quad[1], boundary_quad[2], base_center),
          make_candidate_triangle(boundary_quad[2], boundary_quad[3], base_center),
          make_candidate_triangle(boundary_quad[3], boundary_quad[0], base_center)
        };

        auto find_existing_apex_candidate =
          [&](std::array<EntityRef, 4> &candidate_base,
              EntityRef &candidate_apex_ref,
              Vec3 &candidate_apex,
              std::array<EntityRef, 4> &side_face_refs,
              std::array<bool, 4> &reuse_side_face_refs,
              std::array<int, 4> &side_validation_face_indices) -> bool {
            if(!allow_preassigned_shared_edge_candidate) {
              return false;
            }

            const bool track_preferred_edge_failure =
              prefer_target_point_direction && preferred_shared_edge_index >= 0;
            if(validation_face_index < 0 ||
               validation_face_index >=
                 static_cast<int>(assigned_shared_edge_pyramid_plans.size())) {
              return false;
            }

            const auto &plan =
              assigned_shared_edge_pyramid_plans[static_cast<std::size_t>(validation_face_index)];
            if(!plan.assigned) {
              if(track_preferred_edge_failure) {
                ++preferred_edge_no_candidate_on_edge_count;
              }
              return false;
            }

            if(!is_valid(plan.apex_ref)) {
              if(track_preferred_edge_failure) {
                ++preferred_edge_invalid_apex_ref_count;
              }
              return false;
            }

            if(plan.selected_shared_edge_index < 0 || plan.selected_shared_edge_index >= 4) {
              if(track_preferred_edge_failure) {
                ++preferred_edge_global_fallback_reject_count;
              }
              return false;
            }
            if(!(plan.selected_dihedral_degrees < kMaximumSharedEdgeDihedralDegrees)) {
              if(track_preferred_edge_failure) {
                ++preferred_edge_dihedral_reject_count;
              }
              return false;
            }

            candidate_apex_ref = plan.apex_ref;
            candidate_apex = node_position(plan.apex_ref);
            for(int side = 0; side < 4; ++side) {
              const int matched_face_index =
                plan.matched_face_indices[static_cast<std::size_t>(side)];
              side_validation_face_indices[static_cast<std::size_t>(side)] =
                matched_face_index;
              if(matched_face_index < 0 ||
                 matched_face_index >=
                   static_cast<int>(validation_faces.size())) {
                continue;
              }

              const auto &matched_face =
                validation_faces[static_cast<std::size_t>(matched_face_index)];
              if(!matched_face.active || !is_valid(matched_face.output_face_ref)) {
                if(track_preferred_edge_failure) {
                  ++preferred_edge_missing_output_face_ref_count;
                }
                return false;
              }

              side_face_refs[static_cast<std::size_t>(side)] =
                matched_face.output_face_ref;
              reuse_side_face_refs[static_cast<std::size_t>(side)] = true;
            }

            if(pyramid_volume_measure(candidate_base, candidate_apex) < 0.0) {
              if(track_preferred_edge_failure) {
                ++preferred_edge_negative_volume_reject_count;
              }
              return false;
            }

            return true;
          };

        bool have_existing_apex_candidate = false;
        EntityRef existing_apex_ref {};
        Vec3 existing_apex_candidate {0.0, 0.0, 0.0};
        std::array<EntityRef, 4> existing_side_face_refs {};
        std::array<bool, 4> existing_side_face_valid {false, false, false, false};
        std::array<int, 4> existing_side_validation_indices {-1, -1, -1, -1};
        auto existing_base = boundary_quad;
        have_existing_apex_candidate = find_existing_apex_candidate(
          existing_base,
          existing_apex_ref,
          existing_apex_candidate,
          existing_side_face_refs,
          existing_side_face_valid,
          existing_side_validation_indices);
        if(allow_preassigned_shared_edge_candidate &&
           prefer_target_point_direction &&
           preferred_shared_edge_index >= 0) {
          if(have_existing_apex_candidate) {
            ++preferred_edge_candidate_found_count;
          } else {
            ++preferred_edge_candidate_missing_count;
          }
        }

        auto try_specific_apex =
          [&](EntityRef apex_ref,
              const Vec3 &apex,
              std::array<EntityRef, 4> attempt_base,
              const std::array<EntityRef, 4> &existing_side_refs,
              const std::array<bool, 4> &existing_side_valid,
              const std::array<int, 4> &existing_side_indices) -> QuadEmitOutcome {
            const double projected_height =
              std::abs(vec3_dot(vec3_sub(apex, base_center), outward));
            if(projected_height <= min_height) {
              return QuadEmitOutcome::quality_or_height;
            }
            if(std::abs(pyramid_volume_measure(attempt_base, apex)) <= min_height * min_area) {
              return QuadEmitOutcome::quality_or_height;
            }

            const auto candidate_edges = std::array<ShellEdge, 8> {
              make_candidate_edge(pack_ref(attempt_base[0]), node_position(attempt_base[0]), pack_ref(attempt_base[1]), node_position(attempt_base[1])),
              make_candidate_edge(pack_ref(attempt_base[1]), node_position(attempt_base[1]), pack_ref(attempt_base[2]), node_position(attempt_base[2])),
              make_candidate_edge(pack_ref(attempt_base[2]), node_position(attempt_base[2]), pack_ref(attempt_base[3]), node_position(attempt_base[3])),
              make_candidate_edge(pack_ref(attempt_base[3]), node_position(attempt_base[3]), pack_ref(attempt_base[0]), node_position(attempt_base[0])),
              make_candidate_edge(pack_ref(attempt_base[0]), node_position(attempt_base[0]), kInvalidVertexKey, apex),
              make_candidate_edge(pack_ref(attempt_base[1]), node_position(attempt_base[1]), kInvalidVertexKey, apex),
              make_candidate_edge(pack_ref(attempt_base[2]), node_position(attempt_base[2]), kInvalidVertexKey, apex),
              make_candidate_edge(pack_ref(attempt_base[3]), node_position(attempt_base[3]), kInvalidVertexKey, apex)
            };

            bool shell_collision = false;
            bool wall_collision = false;
            for(const auto &candidate_edge : candidate_edges) {
              if(shell_edge_hits_shell_excluding(
                   validation_face_index,
                   existing_side_indices,
                   candidate_edge)) {
                shell_collision = true;
                break;
              }
              if(shell_edge_hits_non_bl_wall(candidate_edge)) {
                wall_collision = true;
                break;
              }
            }
            if(shell_collision || wall_collision) {
              saw_shell_collision = saw_shell_collision || shell_collision;
              saw_wall_collision = saw_wall_collision || wall_collision;
              return wall_collision ? QuadEmitOutcome::wall_collision : QuadEmitOutcome::shell_collision;
            }

            auto attempt_side_tris = std::array<ShellTriangle, 4> {
              make_candidate_triangle(attempt_base[0], attempt_base[1], apex),
              make_candidate_triangle(attempt_base[1], attempt_base[2], apex),
              make_candidate_triangle(attempt_base[2], attempt_base[3], apex),
              make_candidate_triangle(attempt_base[3], attempt_base[0], apex)
            };
            const auto apex_key = is_valid(apex_ref)
              ? pack_ref(apex_ref)
              : kInvalidVertexKey;
            for(auto &candidate_tri : attempt_side_tris) {
              candidate_tri.vertex_keys[2] = apex_key;
            }
            bool quality_failure = false;
            for(int side = 0; side < 4; ++side) {
              if(existing_side_valid[static_cast<std::size_t>(side)]) {
                continue;
              }
              const auto &candidate_tri = attempt_side_tris[static_cast<std::size_t>(side)];
              if(triangle_area(candidate_tri.tri) <= min_area) {
                quality_failure = true;
                break;
              }
              if(shell_triangle_hits_shell_excluding(
                   validation_face_index,
                   existing_side_indices,
                   candidate_tri)) {
                shell_collision = true;
                break;
              }
              if(shell_triangle_hits_non_bl_wall(candidate_tri)) {
                wall_collision = true;
                break;
              }
            }
            if(quality_failure) {
              return QuadEmitOutcome::quality_or_height;
            }
            if(shell_collision || wall_collision) {
              saw_shell_collision = saw_shell_collision || shell_collision;
              saw_wall_collision = saw_wall_collision || wall_collision;
              return wall_collision ? QuadEmitOutcome::wall_collision : QuadEmitOutcome::shell_collision;
            }

            accepted_pyramid = true;
            accepted_apex_ref = apex_ref;
            accepted_apex = apex;
            quad_base = attempt_base;
            candidate_side_tris = attempt_side_tris;
            reused_side_face_refs = existing_side_refs;
            reuse_existing_side_face = existing_side_valid;
            return QuadEmitOutcome::pyramid;
          };

        if(have_existing_apex_candidate) {
          const auto outcome = try_specific_apex(
            existing_apex_ref,
            existing_apex_candidate,
            existing_base,
            existing_side_face_refs,
            existing_side_face_valid,
            existing_side_validation_indices);
          if(allow_preassigned_shared_edge_candidate &&
             prefer_target_point_direction &&
             preferred_shared_edge_index >= 0) {
            switch(outcome) {
            case QuadEmitOutcome::pyramid:
              ++preferred_edge_candidate_success_count;
              break;
            case QuadEmitOutcome::quality_or_height:
              ++preferred_edge_candidate_quality_count;
              break;
            case QuadEmitOutcome::shell_collision:
              ++preferred_edge_candidate_shell_collision_count;
              break;
            case QuadEmitOutcome::wall_collision:
              ++preferred_edge_candidate_wall_collision_count;
              break;
            }
          }
          if(outcome == QuadEmitOutcome::pyramid) {
            // accepted_pyramid populated through try_specific_apex
          }
        }

        if(allow_inserted_apex_candidate) {
          for(int attempt = accepted_pyramid ? kMaxHeightShrinkAttempts : 0;
              attempt < kMaxHeightShrinkAttempts;
              ++attempt) {
          const double attempt_height = height * std::pow(kHeightShrinkFactor, attempt);
          if(attempt_height <= min_height) {
            break;
          }

          const auto apex = vec3_add(base_center, vec3_scale(pyramid_direction, attempt_height));
          const double projected_height =
            std::abs(vec3_dot(vec3_sub(apex, base_center), outward));
          if(projected_height <= min_height) {
            return QuadEmitOutcome::quality_or_height;
          }

          auto attempt_base = boundary_quad;
          if(pyramid_volume_measure(attempt_base, apex) < 0.0) {
            attempt_base = reverse_quad(attempt_base);
          }
          if(std::abs(pyramid_volume_measure(attempt_base, apex)) <= min_height * min_area) {
            return QuadEmitOutcome::quality_or_height;
          }

          const auto candidate_edges = std::array<ShellEdge, 8> {
            make_candidate_edge(pack_ref(attempt_base[0]), node_position(attempt_base[0]), pack_ref(attempt_base[1]), node_position(attempt_base[1])),
            make_candidate_edge(pack_ref(attempt_base[1]), node_position(attempt_base[1]), pack_ref(attempt_base[2]), node_position(attempt_base[2])),
            make_candidate_edge(pack_ref(attempt_base[2]), node_position(attempt_base[2]), pack_ref(attempt_base[3]), node_position(attempt_base[3])),
            make_candidate_edge(pack_ref(attempt_base[3]), node_position(attempt_base[3]), pack_ref(attempt_base[0]), node_position(attempt_base[0])),
            make_candidate_edge(pack_ref(attempt_base[0]), node_position(attempt_base[0]), kInvalidVertexKey, apex),
            make_candidate_edge(pack_ref(attempt_base[1]), node_position(attempt_base[1]), kInvalidVertexKey, apex),
            make_candidate_edge(pack_ref(attempt_base[2]), node_position(attempt_base[2]), kInvalidVertexKey, apex),
            make_candidate_edge(pack_ref(attempt_base[3]), node_position(attempt_base[3]), kInvalidVertexKey, apex)
          };

          bool shell_collision = false;
          bool wall_collision = false;
          for(const auto &candidate_edge : candidate_edges) {
            if(shell_edge_hits_shell(validation_face_index, candidate_edge)) {
              shell_collision = true;
              break;
            }
            if(shell_edge_hits_non_bl_wall(candidate_edge)) {
              wall_collision = true;
              break;
            }
          }
          if(shell_collision || wall_collision) {
            saw_shell_collision = saw_shell_collision || shell_collision;
            saw_wall_collision = saw_wall_collision || wall_collision;
            continue;
          }

          auto attempt_side_tris = std::array<ShellTriangle, 4> {
            make_candidate_triangle(attempt_base[0], attempt_base[1], apex),
            make_candidate_triangle(attempt_base[1], attempt_base[2], apex),
            make_candidate_triangle(attempt_base[2], attempt_base[3], apex),
            make_candidate_triangle(attempt_base[3], attempt_base[0], apex)
          };
          bool quality_failure = false;
          for(const auto &candidate_tri : attempt_side_tris) {
            if(triangle_area(candidate_tri.tri) <= min_area) {
              quality_failure = true;
              break;
            }
            if(shell_triangle_hits_shell(validation_face_index, candidate_tri)) {
              shell_collision = true;
              break;
            }
            if(shell_triangle_hits_non_bl_wall(candidate_tri)) {
              wall_collision = true;
              break;
            }
          }
          if(quality_failure) {
            return QuadEmitOutcome::quality_or_height;
          }
          if(shell_collision || wall_collision) {
            saw_shell_collision = saw_shell_collision || shell_collision;
            saw_wall_collision = saw_wall_collision || wall_collision;
            continue;
          }

          accepted_pyramid = true;
          accepted_apex_ref = {};
          accepted_apex = apex;
          quad_base = attempt_base;
          candidate_side_tris = attempt_side_tris;
          break;
        }
        }

        if(!accepted_pyramid) {
          if(saw_wall_collision) {
            return QuadEmitOutcome::wall_collision;
          }
          if(saw_shell_collision) {
            return QuadEmitOutcome::shell_collision;
          }
          return QuadEmitOutcome::quality_or_height;
        }

        const auto pyramid_base_face = reverse_quad(quad_base);

        const auto apex_ref = is_valid(accepted_apex_ref)
          ? accepted_apex_ref
          : output.add_node(node_entity_group, {accepted_apex[0], accepted_apex[1], accepted_apex[2]});
        Vec3 pyramid_center = accepted_apex;
        for(int i = 0; i < 4; ++i) {
          const auto p = node_position(quad_base[static_cast<std::size_t>(i)]);
          pyramid_center[0] += p[0];
          pyramid_center[1] += p[1];
          pyramid_center[2] += p[2];
        }
        pyramid_center = {pyramid_center[0] / 5.0, pyramid_center[1] / 5.0, pyramid_center[2] / 5.0};

        auto base_ref = output.add_quad_face(pyramid_base_face_entity_group, pyramid_base_face);
        if(geo::is_valid(topo)) {
          output.set_face_topology_owner(base_ref, topo);
        }

        std::array<EntityRef, 5> pyramid_faces {base_ref, EntityRef {}, EntityRef {}, EntityRef {}, EntityRef {}};
        std::array<bool, 4> add_side_validation = {false, false, false, false};
        std::array<EntityRef, 4> side_face_refs {};
        for(int side = 0; side < 4; ++side) {
          if(reuse_existing_side_face[static_cast<std::size_t>(side)] &&
             is_valid(reused_side_face_refs[static_cast<std::size_t>(side)])) {
            side_face_refs[static_cast<std::size_t>(side)] =
              reused_side_face_refs[static_cast<std::size_t>(side)];
          } else {
            std::array<EntityRef, 3> tri = {
              quad_base[static_cast<std::size_t>(side)],
              quad_base[static_cast<std::size_t>((side + 1) % 4)],
              apex_ref
            };
            tri = orient_tri_outward(tri, pyramid_center);
            auto tri_ref = output.add_triangle_face(pyramid_side_face_entity_group, tri);
            if(geo::is_valid(topo)) {
              output.set_face_topology_owner(tri_ref, topo);
            }
            side_face_refs[static_cast<std::size_t>(side)] = tri_ref;
            add_side_validation[static_cast<std::size_t>(side)] = true;
          }
          pyramid_faces[static_cast<std::size_t>(side + 1)] =
            side_face_refs[static_cast<std::size_t>(side)];
        }

        std::array<EntityRef, 5> pyramid_nodes = {
          quad_base[0], quad_base[1], quad_base[2], quad_base[3], apex_ref
        };
        output.add_pyramid_cell(pyramid_cell_entity_group, pyramid_nodes, pyramid_faces);
        replace_validation_face_with_pyramid(
          validation_face_index,
          candidate_side_tris,
          side_face_refs,
          add_side_validation);
        ++pyramid_count;
        return QuadEmitOutcome::pyramid;
      };

    // -- Create prism/hex cells layer by layer --
    std::size_t prism_count = 0;
    std::size_t hexa_count = 0;

    for(std::size_t layer = 0; layer < num_layers; ++layer) {
      for(std::uint32_t fi = 0; fi < work.surface_faces.size(); ++fi) {
        const auto &face = work.surface_faces[fi];

        // Check if this face is active at this layer
        if(static_cast<int>(layer) >= face_levels[fi]) continue;

        if(face.num_vertices == 3) {
          // Prism: bottom tri + top tri
          const auto &v = face.vertices;
          auto bot0 = node_refs[layer][v[0]];
          auto bot1 = node_refs[layer][v[1]];
          auto bot2 = node_refs[layer][v[2]];
          auto top0 = node_refs[layer + 1][v[0]];
          auto top1 = node_refs[layer + 1][v[1]];
          auto top2 = node_refs[layer + 1][v[2]];

          // Prism nodes: {bot0, bot1, bot2, top0, top1, top2}
          std::array<EntityRef, 6> prism_nodes = {
            bot0, bot1, bot2, top0, top1, top2
          };

          // Prism faces (5): bottom tri, top tri, 3 quad sides
          auto f_bot = output.add_triangle_face(tri_face_entity_group, {bot0, bot2, bot1}); // Reversed for outward
          auto f_top = output.add_triangle_face(tri_face_entity_group, {top0, top1, top2});
          auto f_s0 = output.add_quad_face(quad_face_entity_group, {bot0, bot1, top1, top0});
          auto f_s1 = output.add_quad_face(quad_face_entity_group, {bot1, bot2, top2, top1});
          auto f_s2 = output.add_quad_face(quad_face_entity_group, {bot2, bot0, top0, top2});

          std::array<EntityRef, 5> prism_faces = {f_bot, f_top, f_s0, f_s1, f_s2};
          output.add_prism_cell(prism_cell_entity_group, prism_nodes, prism_faces);
          ++prism_count;

        } else if(face.num_vertices == 4) {
          // Hexa: bottom quad + top quad
          const auto &v = face.vertices;
          auto bot0 = node_refs[layer][v[0]];
          auto bot1 = node_refs[layer][v[1]];
          auto bot2 = node_refs[layer][v[2]];
          auto bot3 = node_refs[layer][v[3]];
          auto top0 = node_refs[layer + 1][v[0]];
          auto top1 = node_refs[layer + 1][v[1]];
          auto top2 = node_refs[layer + 1][v[2]];
          auto top3 = node_refs[layer + 1][v[3]];

          std::array<EntityRef, 8> hexa_nodes = {
            bot0, bot1, bot2, bot3, top0, top1, top2, top3
          };

          auto f_bot = output.add_quad_face(quad_face_entity_group, {bot0, bot3, bot2, bot1});
          auto f_top = output.add_quad_face(quad_face_entity_group, {top0, top1, top2, top3});
          auto f_s0 = output.add_quad_face(quad_face_entity_group, {bot0, bot1, top1, top0});
          auto f_s1 = output.add_quad_face(quad_face_entity_group, {bot1, bot2, top2, top1});
          auto f_s2 = output.add_quad_face(quad_face_entity_group, {bot2, bot3, top3, top2});
          auto f_s3 = output.add_quad_face(quad_face_entity_group, {bot3, bot0, top0, top3});

          std::array<EntityRef, 6> hexa_faces = {f_bot, f_top, f_s0, f_s1, f_s2, f_s3};
          output.add_hexa_cell(hexa_cell_entity_group, hexa_nodes, hexa_faces);
          ++hexa_count;
        }
      }
    }

    bl_log_info("Created %zu prism cells, %zu hexa cells", prism_count, hexa_count);

    // -- Extract top cap faces (for subsequent tet fill) --
    std::size_t top_face_count = 0;
    for(std::uint32_t fi = 0; fi < work.surface_faces.size(); ++fi) {
      const auto &face = work.surface_faces[fi];
      const auto actual_layer = static_cast<std::size_t>(face_levels[fi]);
      const auto topo = face_topology(fi);
      const auto inside_point = actual_layer > 0
        ? face_anchor_center(fi, actual_layer - 1)
        : collapsed_face_inside(fi);

      const auto &v = face.vertices;
      if(face.num_vertices == 3) {
        std::array<EntityRef, 3> tri = {
          node_refs[actual_layer][v[0]],
          node_refs[actual_layer][v[1]],
          node_refs[actual_layer][v[2]]
        };
        tri = orient_tri_outward(tri, inside_point);
        auto tri_ref = output.add_triangle_face(top_tri_face_entity_group, tri);
        if(geo::is_valid(topo)) {
          output.set_face_topology_owner(tri_ref, topo);
        }
        set_validation_face_output_ref(top_validation_face_indices[fi], tri_ref);
        ++top_tri_count;
      } else if(face.num_vertices == 4) {
        std::array<EntityRef, 4> quad = {
          node_refs[actual_layer][v[0]],
          node_refs[actual_layer][v[1]],
          node_refs[actual_layer][v[2]],
          node_refs[actual_layer][v[3]]
        };
        quad = orient_quad_outward(quad, inside_point);
        pending_top_quad_faces.push_back({
          quad,
          inside_point,
          topo,
          top_validation_face_indices[fi]
        });
      }
      ++top_face_count;
    }

    const std::size_t transition_face_count = pending_transition_faces.size();
    std::vector<bool> top_quad_finished(pending_top_quad_faces.size(), false);
    std::vector<bool> transition_quad_finished(pending_transition_faces.size(), false);

    auto has_assigned_shared_edge_plan = [&](int validation_face_index) {
      return validation_face_index >= 0 &&
             validation_face_index < static_cast<int>(assigned_shared_edge_pyramid_plans.size()) &&
             assigned_shared_edge_pyramid_plans[static_cast<std::size_t>(validation_face_index)].assigned;
    };

    std::vector<std::size_t> deferred_shared_edge_top_indices;
    std::vector<std::size_t> deferred_shared_edge_transition_indices;
    deferred_shared_edge_top_indices.reserve(pending_top_quad_faces.size());
    deferred_shared_edge_transition_indices.reserve(pending_transition_faces.size());
    for(std::size_t top_index = 0; top_index < pending_top_quad_faces.size(); ++top_index) {
      if(has_assigned_shared_edge_plan(
           pending_top_quad_faces[top_index].validation_face_index)) {
        deferred_shared_edge_top_indices.push_back(top_index);
      }
    }
    for(std::size_t transition_index = 0;
        transition_index < pending_transition_faces.size();
        ++transition_index) {
      if(has_assigned_shared_edge_plan(
           pending_transition_faces[transition_index].validation_face_index)) {
        deferred_shared_edge_transition_indices.push_back(transition_index);
      }
    }

    constexpr int kMaxSharedEdgeRetryPasses = 16;
    for(int pass = 0;
        pass < kMaxSharedEdgeRetryPasses &&
        (!deferred_shared_edge_top_indices.empty() ||
         !deferred_shared_edge_transition_indices.empty());
        ++pass) {
      std::vector<std::size_t> next_shared_edge_top_indices;
      std::vector<std::size_t> next_shared_edge_transition_indices;
      next_shared_edge_top_indices.reserve(deferred_shared_edge_top_indices.size());
      next_shared_edge_transition_indices.reserve(deferred_shared_edge_transition_indices.size());
      std::size_t accepted_this_pass = 0;

      for(const auto top_index : deferred_shared_edge_top_indices) {
        if(top_quad_finished[top_index]) {
          continue;
        }
        const auto &pending_top_quad = pending_top_quad_faces[top_index];
        const auto outcome = emit_pyramid(
          pending_top_quad.quad,
          pending_top_quad.inside_point,
          pending_top_quad.topo,
          pending_top_quad.validation_face_index,
          false,
          true,
          false,
          -1);
        if(outcome == QuadEmitOutcome::pyramid) {
          top_quad_finished[top_index] = true;
          ++accepted_this_pass;
          continue;
        }
        if(outcome == QuadEmitOutcome::shell_collision) {
          next_shared_edge_top_indices.push_back(top_index);
        }
      }

      for(const auto transition_index : deferred_shared_edge_transition_indices) {
        if(transition_quad_finished[transition_index]) {
          continue;
        }
        const auto &pending_transition_face =
          pending_transition_faces[transition_index];
        const auto outcome = emit_pyramid(
          pending_transition_face.quad,
          pending_transition_face.inside_point,
          pending_transition_face.topo,
          pending_transition_face.validation_face_index,
          true,
          true,
          false,
          pending_transition_face.preferred_shared_edge_index);
        if(outcome == QuadEmitOutcome::pyramid) {
          transition_quad_finished[transition_index] = true;
          ++accepted_this_pass;
          continue;
        }
        if(outcome == QuadEmitOutcome::shell_collision) {
          next_shared_edge_transition_indices.push_back(transition_index);
        }
      }

      deferred_shared_edge_top_indices.swap(next_shared_edge_top_indices);
      deferred_shared_edge_transition_indices.swap(next_shared_edge_transition_indices);
      if(accepted_this_pass == 0) {
        break;
      }
    }

    for(std::size_t top_index = 0; top_index < pending_top_quad_faces.size(); ++top_index) {
      if(top_quad_finished[top_index]) {
        continue;
      }
      const auto &pending_top_quad = pending_top_quad_faces[top_index];
      const auto outcome = emit_pyramid(
        pending_top_quad.quad,
        pending_top_quad.inside_point,
        pending_top_quad.topo,
        pending_top_quad.validation_face_index,
        false,
        false,
        true,
        -1);
      if(outcome != QuadEmitOutcome::pyramid) {
        record_rejected_pyramid(outcome);
        const auto quad_ref = emit_boundary_quad(
          pending_top_quad.quad,
          pending_top_quad.topo,
          top_quad_face_entity_group,
          direct_top_quad_count);
        set_validation_face_output_ref(pending_top_quad.validation_face_index, quad_ref);
      }
      top_quad_finished[top_index] = true;
    }

    std::vector<std::size_t> deferred_transition_indices;
    deferred_transition_indices.reserve(pending_transition_faces.size());
    for(std::size_t transition_index = 0;
        transition_index < pending_transition_faces.size();
        ++transition_index) {
      if(transition_quad_finished[transition_index]) {
        continue;
      }
      const auto &pending_transition_face = pending_transition_faces[transition_index];
      const auto outcome = emit_pyramid(
        pending_transition_face.quad,
        pending_transition_face.inside_point,
        pending_transition_face.topo,
        pending_transition_face.validation_face_index,
        true,
        false,
        true,
        pending_transition_face.preferred_shared_edge_index);
      if(outcome == QuadEmitOutcome::pyramid) {
        transition_quad_finished[transition_index] = true;
        continue;
      }
      if(outcome == QuadEmitOutcome::shell_collision) {
        deferred_transition_indices.push_back(transition_index);
        continue;
      }

      record_rejected_pyramid(outcome);
      const auto quad_ref = emit_boundary_quad(
        pending_transition_face.quad,
        pending_transition_face.topo,
        transition_quad_face_entity_group,
        direct_transition_quad_count);
      set_validation_face_output_ref(pending_transition_face.validation_face_index, quad_ref);
      transition_quad_finished[transition_index] = true;
    }

    constexpr int kMaxTransitionRetryPasses = 16;
    for(int pass = 0;
        pass < kMaxTransitionRetryPasses && !deferred_transition_indices.empty();
        ++pass) {
      std::vector<std::size_t> next_deferred_indices;
      next_deferred_indices.reserve(deferred_transition_indices.size());
      std::size_t accepted_this_pass = 0;

      for(const auto transition_index : deferred_transition_indices) {
        const auto &pending_transition_face =
          pending_transition_faces[transition_index];
        const auto outcome = emit_pyramid(
          pending_transition_face.quad,
          pending_transition_face.inside_point,
          pending_transition_face.topo,
          pending_transition_face.validation_face_index,
          true,
          false,
          true,
          pending_transition_face.preferred_shared_edge_index);
        if(outcome == QuadEmitOutcome::pyramid) {
          transition_quad_finished[transition_index] = true;
          ++accepted_this_pass;
          continue;
        }
        if(outcome == QuadEmitOutcome::shell_collision) {
          next_deferred_indices.push_back(transition_index);
          continue;
        }

        record_rejected_pyramid(outcome);
        const auto quad_ref = emit_boundary_quad(
          pending_transition_face.quad,
          pending_transition_face.topo,
          transition_quad_face_entity_group,
          direct_transition_quad_count);
        set_validation_face_output_ref(pending_transition_face.validation_face_index, quad_ref);
        transition_quad_finished[transition_index] = true;
      }

      deferred_transition_indices.swap(next_deferred_indices);
      if(accepted_this_pass == 0) {
        break;
      }
    }

    for(const auto transition_index : deferred_transition_indices) {
      const auto &pending_transition_face =
        pending_transition_faces[transition_index];
      record_rejected_pyramid(QuadEmitOutcome::shell_collision);
      const auto quad_ref = emit_boundary_quad(
        pending_transition_face.quad,
        pending_transition_face.topo,
        transition_quad_face_entity_group,
        direct_transition_quad_count);
      set_validation_face_output_ref(pending_transition_face.validation_face_index, quad_ref);
      transition_quad_finished[transition_index] = true;
    }

    // -- Extract bottom faces (original surface = wall boundary) --
    std::size_t bot_face_count = 0;
    for(std::uint32_t fi = 0; fi < work.surface_faces.size(); ++fi) {
      const auto &face = work.surface_faces[fi];
      const auto &v = face.vertices;
      const auto topo = face_topology(fi);

      if(face.num_vertices == 3) {
        auto bottom_ref = output.add_triangle_face(bot_tri_face_entity_group,
          {node_refs[0][v[0]], node_refs[0][v[2]], node_refs[0][v[1]]});
        if(geo::is_valid(topo)) {
          output.set_face_topology_owner(bottom_ref, topo);
        }
        ++bot_face_count;
      } else if(face.num_vertices == 4) {
        auto bottom_ref = output.add_quad_face(bot_quad_face_entity_group,
          {node_refs[0][v[0]], node_refs[0][v[3]], node_refs[0][v[2]], node_refs[0][v[1]]});
        if(geo::is_valid(topo)) {
          output.set_face_topology_owner(bottom_ref, topo);
        }
        ++bot_face_count;
      }
    }

    bl_log_info(
      "Boundary shell: %zu top faces, %zu top triangles, %zu transition faces, %zu pyramids, %zu direct top quads, %zu direct transition quads, %zu pyramid quality collapses, %zu pyramid shell collisions, %zu pyramid wall collisions, %zu bottom wall faces",
      top_face_count,
      top_tri_count,
      transition_face_count,
      pyramid_count,
      direct_top_quad_count,
      direct_transition_quad_count,
      rejected_pyramid_quality_count,
      rejected_pyramid_shell_intersection_count,
      rejected_pyramid_wall_intersection_count,
      bot_face_count);
    std::fprintf(
      stderr,
      "[info] Boundary shell summary: top=%zu, top_tri=%zu, transition=%zu, pyramids=%zu, direct_top_quads=%zu, direct_transition_quads=%zu, transition_target_side_overrides=%zu, transition_target_same_side=%zu, pyramid_quality_rejects=%zu, pyramid_shell_collisions=%zu, pyramid_wall_collisions=%zu, bottom=%zu\n",
      top_face_count,
      top_tri_count,
      transition_face_count,
      pyramid_count,
      direct_top_quad_count,
      direct_transition_quad_count,
      transition_target_side_override_count,
      transition_target_same_side_count,
      rejected_pyramid_quality_count,
      rejected_pyramid_shell_intersection_count,
      rejected_pyramid_wall_intersection_count,
      bot_face_count);
    std::fprintf(
      stderr,
      "[info] Preferred shared-edge apex: found=%zu, missing=%zu, success=%zu, quality_failures=%zu, shell_collisions=%zu, wall_collisions=%zu\n",
      preferred_edge_candidate_found_count,
      preferred_edge_candidate_missing_count,
      preferred_edge_candidate_success_count,
      preferred_edge_candidate_quality_count,
      preferred_edge_candidate_shell_collision_count,
      preferred_edge_candidate_wall_collision_count);
    std::fprintf(
      stderr,
      "[info] Preferred shared-edge reject breakdown: no_edge_match=%zu, multi_edge_match=%zu, global_fallback_reject=%zu, invalid_apex=%zu, dihedral=%zu, missing_output_face=%zu, negative_volume=%zu\n",
      preferred_edge_no_candidate_on_edge_count,
      preferred_edge_multiple_candidates_on_edge_count,
      preferred_edge_global_fallback_reject_count,
      preferred_edge_invalid_apex_ref_count,
      preferred_edge_dihedral_reject_count,
      preferred_edge_missing_output_face_ref_count,
      preferred_edge_negative_volume_reject_count);

    return core::detail::clear_error_state();
  }

  void export_bl_vtk(const BLWorkData &work, const std::string &path)
  {
    const auto num_layers = work.layer_positions.size() - 1;
    if(num_layers == 0) return;

    const auto nn = work.surface_nodes.size();
    const auto total_nodes = nn * (num_layers + 1);

    // Count cells
    std::size_t prism_count = 0, hexa_count = 0;
    for(std::uint32_t fi = 0; fi < work.surface_faces.size(); ++fi) {
      const auto &face = work.surface_faces[fi];
      int face_layers = std::min(work.max_face_level[fi], static_cast<int>(num_layers));
      if(face.num_vertices == 3) prism_count += face_layers;
      else if(face.num_vertices == 4) hexa_count += face_layers;
    }

    std::ofstream out(path, std::ios::trunc);
    if(!out) return;

    out << "# vtk DataFile Version 3.0\n";
    out << "Boundary Layer Mesh\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";
    out << std::setprecision(15);

    // Write points: layer 0 .. layer N
    out << "POINTS " << total_nodes << " double\n";
    for(std::size_t layer = 0; layer <= num_layers; ++layer) {
      for(std::uint32_t ni = 0; ni < nn; ++ni) {
        const auto &p = work.layer_positions[layer][ni];
        out << p[0] << " " << p[1] << " " << p[2] << "\n";
      }
    }

    auto node_id = [&](std::size_t layer, std::uint32_t ni) -> std::size_t {
      return layer * nn + ni;
    };

    // Connectivity size
    const std::size_t total_cells = prism_count + hexa_count;
    const std::size_t conn_size = prism_count * 7 + hexa_count * 9; // 1+6 per prism, 1+8 per hexa

    out << "CELLS " << total_cells << " " << conn_size << "\n";
    for(std::size_t layer = 0; layer < num_layers; ++layer) {
      for(std::uint32_t fi = 0; fi < work.surface_faces.size(); ++fi) {
        const auto &face = work.surface_faces[fi];
        if(static_cast<int>(layer) >= work.max_face_level[fi]) continue;

        if(face.num_vertices == 3) {
          // VTK_WEDGE (13): bot0 bot1 bot2 top0 top1 top2
          out << "6"
              << " " << node_id(layer, face.vertices[0])
              << " " << node_id(layer, face.vertices[1])
              << " " << node_id(layer, face.vertices[2])
              << " " << node_id(layer+1, face.vertices[0])
              << " " << node_id(layer+1, face.vertices[1])
              << " " << node_id(layer+1, face.vertices[2])
              << "\n";
        } else if(face.num_vertices == 4) {
          // VTK_HEXAHEDRON (12): bot0 bot1 bot2 bot3 top0 top1 top2 top3
          out << "8"
              << " " << node_id(layer, face.vertices[0])
              << " " << node_id(layer, face.vertices[1])
              << " " << node_id(layer, face.vertices[2])
              << " " << node_id(layer, face.vertices[3])
              << " " << node_id(layer+1, face.vertices[0])
              << " " << node_id(layer+1, face.vertices[1])
              << " " << node_id(layer+1, face.vertices[2])
              << " " << node_id(layer+1, face.vertices[3])
              << "\n";
        }
      }
    }

    out << "CELL_TYPES " << total_cells << "\n";
    for(std::size_t layer = 0; layer < num_layers; ++layer) {
      for(std::uint32_t fi = 0; fi < work.surface_faces.size(); ++fi) {
        if(static_cast<int>(layer) >= work.max_face_level[fi]) continue;
        const auto &face = work.surface_faces[fi];
        if(face.num_vertices == 3) out << "13\n";      // VTK_WEDGE
        else if(face.num_vertices == 4) out << "12\n";  // VTK_HEXAHEDRON
      }
    }

    bl_log_info(
      "Exported BL mesh to %s: %zu prisms, %zu hexas, %zu nodes",
      path.c_str(),
      prism_count,
      hexa_count,
      total_nodes);
  }

  void export_top_cap_vtk(const BLWorkData &work, const std::string &path)
  {
    const auto num_layers = work.layer_positions.size() - 1;
    if(num_layers == 0) return;

    const auto nn = work.surface_nodes.size();

    // Collect top cap face vertices and faces
    // Each face's top cap is at its actual_layer = min(max_face_level, num_layers)
    struct TopFace {
      int num_verts;
      std::uint32_t verts[4];
      std::size_t layer;
    };
    std::vector<TopFace> top_faces;

    for(std::uint32_t fi = 0; fi < work.surface_faces.size(); ++fi) {
      const auto &face = work.surface_faces[fi];
      int face_max = work.max_face_level[fi];
      if(face_max <= 0) continue;
      auto actual_layer = static_cast<std::size_t>(
        std::min(face_max, static_cast<int>(num_layers)));

      TopFace tf;
      tf.num_verts = face.num_vertices;
      tf.layer = actual_layer;
      for(int v = 0; v < face.num_vertices; ++v)
        tf.verts[v] = face.vertices[v];
      top_faces.push_back(tf);
    }

    if(top_faces.empty()) return;

    // Collect unique node indices (layer, node_id) -> output index
    std::map<std::pair<std::size_t, std::uint32_t>, std::size_t> node_map;
    std::vector<std::array<double, 3>> points;

    auto get_node = [&](std::size_t layer, std::uint32_t ni) -> std::size_t {
      auto key = std::make_pair(layer, ni);
      auto it = node_map.find(key);
      if(it != node_map.end()) return it->second;
      auto idx = points.size();
      points.push_back({work.layer_positions[layer][ni][0],
                        work.layer_positions[layer][ni][1],
                        work.layer_positions[layer][ni][2]});
      node_map[key] = idx;
      return idx;
    };

    // Resolve all face node indices
    std::vector<std::vector<std::size_t>> face_indices(top_faces.size());
    for(std::size_t i = 0; i < top_faces.size(); ++i) {
      const auto &tf = top_faces[i];
      face_indices[i].resize(tf.num_verts);
      for(int v = 0; v < tf.num_verts; ++v)
        face_indices[i][v] = get_node(tf.layer, tf.verts[v]);
    }

    std::ofstream out(path, std::ios::trunc);
    if(!out) return;

    out << "# vtk DataFile Version 3.0\n";
    out << "BL Top Cap Surface\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";
    out << std::setprecision(15);

    out << "POINTS " << points.size() << " double\n";
    for(const auto &p : points)
      out << p[0] << " " << p[1] << " " << p[2] << "\n";

    // Connectivity
    std::size_t conn_size = 0;
    for(const auto &tf : top_faces)
      conn_size += 1 + tf.num_verts;

    out << "CELLS " << top_faces.size() << " " << conn_size << "\n";
    for(std::size_t i = 0; i < top_faces.size(); ++i) {
      out << top_faces[i].num_verts;
      for(int v = 0; v < top_faces[i].num_verts; ++v)
        out << " " << face_indices[i][v];
      out << "\n";
    }

    out << "CELL_TYPES " << top_faces.size() << "\n";
    for(const auto &tf : top_faces) {
      if(tf.num_verts == 3) out << "5\n";       // VTK_TRIANGLE
      else if(tf.num_verts == 4) out << "9\n";   // VTK_QUAD
    }

    bl_log_info(
      "Exported BL top cap to %s: %zu faces, %zu nodes",
      path.c_str(),
      top_faces.size(),
      points.size());
  }

  void export_top_cap_obj(const BLWorkData &work, const std::string &path)
  {
    const auto num_layers = work.layer_positions.size() - 1;
    if(num_layers == 0) return;

    // Collect unique nodes and faces
    std::map<std::pair<std::size_t, std::uint32_t>, std::size_t> node_map;
    std::vector<std::array<double, 3>> points;

    auto get_node = [&](std::size_t layer, std::uint32_t ni) -> std::size_t {
      auto key = std::make_pair(layer, ni);
      auto it = node_map.find(key);
      if(it != node_map.end()) return it->second;
      auto idx = points.size();
      points.push_back({work.layer_positions[layer][ni][0],
                        work.layer_positions[layer][ni][1],
                        work.layer_positions[layer][ni][2]});
      node_map[key] = idx;
      return idx;
    };

    struct ObjFace { int nv; std::size_t v[4]; };
    std::vector<ObjFace> faces;

    for(std::uint32_t fi = 0; fi < work.surface_faces.size(); ++fi) {
      const auto &face = work.surface_faces[fi];
      int face_max = work.max_face_level[fi];
      if(face_max <= 0) continue;
      auto actual_layer = static_cast<std::size_t>(
        std::min(face_max, static_cast<int>(num_layers)));

      ObjFace of;
      of.nv = face.num_vertices;
      for(int v = 0; v < face.num_vertices; ++v)
        of.v[v] = get_node(actual_layer, face.vertices[v]);
      faces.push_back(of);
    }

    std::ofstream out(path, std::ios::trunc);
    if(!out) return;

    out << "# BL Top Cap Surface\n";
    out << "# " << points.size() << " vertices, " << faces.size() << " faces\n";
    out << std::setprecision(15);

    for(const auto &p : points)
      out << "v " << p[0] << " " << p[1] << " " << p[2] << "\n";

    for(const auto &f : faces) {
      out << "f";
      for(int v = 0; v < f.nv; ++v)
        out << " " << (f.v[v] + 1); // OBJ is 1-indexed
      out << "\n";
    }

    bl_log_info(
      "Exported BL top cap OBJ to %s: %zu faces, %zu nodes",
      path.c_str(),
      faces.size(),
      points.size());
  }

  MeshingRequest request_ {};
  BoundaryLayerParams params_ {};
  std::set<std::uint32_t> bl_face_ids_; // Empty = all faces
  int bl_target_region_ = -1;           // -1 = auto (centroid), >=0 = target region
  std::array<double, 3> bl_target_point_ = {0, 0, 0};
  bool has_target_point_ = false;
};

} // namespace

MeshingAlgorithmPtr create_boundary_layer_mesher()
{
  return std::make_unique<BoundaryLayerAlgorithm>();
}

} // namespace sqmesh::mesh::detail
