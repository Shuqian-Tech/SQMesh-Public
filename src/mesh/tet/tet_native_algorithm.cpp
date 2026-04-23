#include "mesh/framework/meshing_framework.hpp"

#include "core/log.hpp"
#include "core/runtime_registry.hpp"
#include "mesh/tet/core/tet_core.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace sqmesh::mesh::detail {

using sqmesh::mesh::tet::detail::TetMeshCore;
using sqmesh::mesh::tet::detail::TetMeshData;
using sqmesh::mesh::tet::detail::TetMeshBehavior;
using sqmesh::mesh::tet::detail::run_tet_mesh_core;

namespace {

[[nodiscard]] base::StatusCode fail_tet_core(
  base::StatusCode code,
  std::string_view message
)
{
  return core::detail::publish_error(code, message);
}

} // namespace

class NativeTetVolumeAlgorithm final : public MeshingAlgorithm {
public:
  [[nodiscard]] std::string_view name() const noexcept override {
    return "Tetrahedral Volume Mesher";
  }

  [[nodiscard]] std::string_view entity_group_prefix() const noexcept override {
    return "volume_";
  }

protected:
  double min_length_ = 0.0;
  double max_length_ = 0.0;
  double growth_rate_ = 1.5;
  double quality_ratio_ = 2.0;
  double max_volume_ = 0.0;
  std::vector<std::array<double, 3>> material_points_;
  std::vector<std::array<double, 3>> hole_points_;

  void reset() noexcept override {
    MeshingAlgorithm::reset();
    material_points_.clear();
    hole_points_.clear();
  }

  [[nodiscard]] base::StatusCode on_initialize(
    const MeshingRequest &request
  ) override {
    if(request.target_dimension != MeshingDimension::volume) {
      return base::StatusCode::invalid_argument;
    }
    return base::StatusCode::ok;
  }

  [[nodiscard]] base::StatusCode on_configure(
    const ParameterDictionary &parameters
  ) override {
    static_cast<void>(parameters.try_get_number("minimum_length", min_length_));
    static_cast<void>(parameters.try_get_number("maximum_length", max_length_));
    static_cast<void>(parameters.try_get_number("growth_rate", growth_rate_));
    static_cast<void>(parameters.try_get_number("quality_ratio", quality_ratio_));
    static_cast<void>(parameters.try_get_number("maximum_volume", max_volume_));

    if(max_volume_ <= 0.0 && max_length_ > 0.0) {
      max_volume_ = (max_length_ * max_length_ * max_length_) * 0.11785113;
    }

    // Shared "x1,y1,z1;x2,y2,z2;..." parser for material_points and hole_points.
    auto parse_points = [](std::string_view text) {
      std::vector<std::array<double, 3>> out;
      std::istringstream iss((std::string(text)));
      std::string point_str;
      while(std::getline(iss, point_str, ';')) {
        std::istringstream pss(point_str);
        std::string coord;
        std::array<double, 3> pt{{0.0, 0.0, 0.0}};
        int ci = 0;
        while(std::getline(pss, coord, ',') && ci < 3) {
          try {
            pt[ci++] = std::stod(coord);
          } catch(...) {
            ci = -1;
            break;
          }
        }
        if(ci == 3) {
          out.push_back(pt);
        }
      }
      return out;
    };

    // "material_points" = "x1,y1,z1;x2,y2,z2;..." — regions to keep.
    // Each material point tags one enclosed region; tets outside any
    // tagged region are discarded after meshing. When unspecified, the
    // full enclosed volume is filled.
    material_points_.clear();
    std::string_view mp_text;
    if(parameters.try_get_text("material_points", mp_text)) {
      material_points_ = parse_points(mp_text);
      SQMESH_LOG_INFO(
        "Tet material_points: {} region point(s) parsed",
        material_points_.size()
      );
    }

    // "hole_points" = "x1,y1,z1;x2,y2,z2;..." — regions to skip.
    // Useful when a region is not fully sealed by the surface and cannot
    // be isolated via material_points alone (e.g. BL top with pinholes).
    hole_points_.clear();
    std::string_view hp_text;
    if(parameters.try_get_text("hole_points", hp_text)) {
      hole_points_ = parse_points(hp_text);
      SQMESH_LOG_INFO(
        "Tet hole_points: {} skip point(s) parsed",
        hole_points_.size()
      );
    }
    return base::StatusCode::ok;
  }

  [[nodiscard]] base::StatusCode on_generate(Domain &output) override {
#ifdef SQMESH_HAS_SPDLOG
    SQMESH_LOG_INFO("Tetrahedral Volume Mesher running...");
#else
    SQMESH_LOG_INFO("%s", "Tetrahedral Volume Mesher running...");
#endif

    std::unordered_map<std::uint64_t, int> node_ref_to_index;
    std::vector<std::array<double, 3>> input_points;
    struct InputFace {
      int v[4];            // up to 4 verts (triangle or quad)
      int num_vertices;
      int marker;
    };
    std::vector<InputFace> input_faces;

    auto get_or_add_point = [&](sqmesh::mesh::EntityRef node_ref) -> int {
      const auto key =
        (static_cast<std::uint64_t>(node_ref.entity_group) << 32U) |
        static_cast<std::uint64_t>(node_ref.index);
      auto it = node_ref_to_index.find(key);
      if(it != node_ref_to_index.end()) return it->second;
      const int new_idx = static_cast<int>(input_points.size());
      node_ref_to_index[key] = new_idx;
      const auto &node = output.node(node_ref);
      input_points.push_back(
        {node.coordinates[0], node.coordinates[1], node.coordinates[2]}
      );
      return new_idx;
    };

    for(const auto &entity_group : output.entity_groups()) {
      if(entity_group.order() != EntityOrder::face ||
         entity_group.role() != EntityGroupRole::computational) {
        continue;
      }
      const int marker = entity_group.zone_id();
      for(std::uint32_t fi = 0; fi < entity_group.faces().size(); ++fi) {
        const auto face_nodes = output.face_nodes({entity_group.id(), fi});
        if(face_nodes.size != 3U && face_nodes.size != 4U) {
          continue;
        }
        // Pass quads through as 4-vertex polygons rather than pre-splitting
        // them into two triangles. Forcing a diagonal on the SQMesh side
        // can introduce edges that the tet boundary-recovery cannot recover
        // when neighbouring facets disagree about the diagonal choice.
        InputFace face;
        face.marker = marker;
        face.num_vertices = static_cast<int>(face_nodes.size);
        for(std::size_t ni = 0; ni < face_nodes.size; ++ni) {
          face.v[ni] = get_or_add_point(face_nodes[ni]);
        }
        input_faces.push_back(face);
      }
    }

    if(input_points.empty() || input_faces.empty()) {
#ifdef SQMESH_HAS_SPDLOG
      SQMESH_LOG_ERROR("No surface mesh found in Domain — cannot generate volume mesh");
#else
      SQMESH_LOG_ERROR("%s", "No surface mesh found in Domain — cannot generate volume mesh");
#endif
      return fail_tet_core(
        base::StatusCode::invalid_argument,
        "Tet volume mesher requires a seeded surface mesh domain"
      );
    }

#ifdef SQMESH_HAS_SPDLOG
    SQMESH_LOG_INFO(
      "Extracted surface: {} vertices, {} triangles",
      input_points.size(),
      input_faces.size()
    );
#else
    SQMESH_LOG_INFO(
      "Extracted surface: %zu vertices, %zu triangles",
      input_points.size(),
      input_faces.size()
    );
#endif

    TetMeshData in, out_tet;
    in.firstnumber = 0;

    in.numberofpoints = static_cast<int>(input_points.size());
    in.pointlist = new double[in.numberofpoints * 3];
    for(int i = 0; i < in.numberofpoints; ++i) {
      in.pointlist[i * 3 + 0] = input_points[i][0];
      in.pointlist[i * 3 + 1] = input_points[i][1];
      in.pointlist[i * 3 + 2] = input_points[i][2];
    }

    in.numberoffacets = static_cast<int>(input_faces.size());
    in.facetlist = new TetMeshData::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];

    for(int i = 0; i < in.numberoffacets; ++i) {
      TetMeshData::facet *f = &in.facetlist[i];
      f->polygonlist = new TetMeshData::polygon[1];
      f->numberofpolygons = 1;
      f->holelist = nullptr;
      f->numberofholes = 0;

      TetMeshData::polygon *p = &f->polygonlist[0];
      const int nv = input_faces[i].num_vertices;
      p->numberofvertices = nv;
      p->vertexlist = new int[nv];
      for(int vi = 0; vi < nv; ++vi) {
        p->vertexlist[vi] = input_faces[i].v[vi];
      }

      in.facetmarkerlist[i] = input_faces[i].marker;
    }

    // Per-point target size: average incident surface edge length,
    // clamped to [min_length, max_length].
    {
      const int np = in.numberofpoints;
      const int nf = static_cast<int>(input_faces.size());
      std::vector<double> edge_len_sum(static_cast<std::size_t>(np), 0.0);
      std::vector<int> edge_count(static_cast<std::size_t>(np), 0);
      for(int fi = 0; fi < nf; ++fi) {
        const int nvf = input_faces[fi].num_vertices;
        for(int ei = 0; ei < nvf; ++ei) {
          const int v0 = input_faces[fi].v[ei];
          const int v1 = input_faces[fi].v[(ei + 1) % nvf];
          const double dx = in.pointlist[v0 * 3]     - in.pointlist[v1 * 3];
          const double dy = in.pointlist[v0 * 3 + 1] - in.pointlist[v1 * 3 + 1];
          const double dz = in.pointlist[v0 * 3 + 2] - in.pointlist[v1 * 3 + 2];
          const double len = std::sqrt(dx * dx + dy * dy + dz * dz);
          edge_len_sum[v0] += len; edge_count[v0] += 1;
          edge_len_sum[v1] += len; edge_count[v1] += 1;
        }
      }
      const double min_size = min_length_ > 0.0 ? min_length_ : 1e-6;
      const double max_size = max_length_ > 0.0 ? max_length_ : 1e10;
      in.numberofpointmtrs = 1;
      in.pointmtrlist = new double[np];
      for(int i = 0; i < np; ++i) {
        double local_size = edge_count[i] > 0
          ? edge_len_sum[i] / static_cast<double>(edge_count[i])
          : max_size;
        local_size = std::max(local_size, min_size);
        local_size = std::min(local_size, max_size);
        in.pointmtrlist[i] = local_size;
      }
    }

    // Configure the tet core via the structured TetMeshBehavior carrier.
    TetMeshBehavior b;
    b.plc = 1;
    b.nobisect = 1;
    b.quality = 1;
    b.metric = 1;
    b.minratio = quality_ratio_;
    b.verbose = 1;

    // Region-attribute path: each material point tags one region with a
    // unique positive attribute; attributes propagate to every tet inside
    // the enclosing closed surface, so the ingestion loop below can drop
    // tets from untagged regions. Volume constraints are carried per-region
    // when regions are active; a global fixedvolume switch is used only in
    // the single-region fallback.
    const bool use_regions = !material_points_.empty();
    if(use_regions) {
      b.regionattrib = 1;
      in.numberofregions = static_cast<int>(material_points_.size());
      in.regionlist = new double[in.numberofregions * 5];
      for(int ri = 0; ri < in.numberofregions; ++ri) {
        in.regionlist[ri * 5 + 0] = material_points_[ri][0];
        in.regionlist[ri * 5 + 1] = material_points_[ri][1];
        in.regionlist[ri * 5 + 2] = material_points_[ri][2];
        in.regionlist[ri * 5 + 3] = static_cast<double>(ri + 1);
        in.regionlist[ri * 5 + 4] = max_volume_ > 0.0 ? max_volume_ : -1.0;
      }
    } else if(max_volume_ > 0.0) {
      b.fixedvolume = 1;
      b.maxvolume = max_volume_;
      b.maxvolume_length = std::pow(max_volume_, 1.0 / 3.0) / 3.0;
    }

    if(!hole_points_.empty()) {
      in.numberofholes = static_cast<int>(hole_points_.size());
      in.holelist = new double[in.numberofholes * 3];
      for(int hi = 0; hi < in.numberofholes; ++hi) {
        in.holelist[hi * 3 + 0] = hole_points_[hi][0];
        in.holelist[hi * 3 + 1] = hole_points_[hi][1];
        in.holelist[hi * 3 + 2] = hole_points_[hi][2];
      }
    }

#ifdef SQMESH_HAS_SPDLOG
    SQMESH_LOG_INFO(
      "Running tet core with behavior: plc={} nobisect={} quality={} minratio={} verbose={} fixedvolume={} maxvolume={}",
      b.plc, b.nobisect, b.quality, b.minratio, b.verbose, b.fixedvolume, b.maxvolume);
#else
    SQMESH_LOG_INFO(
      "Running tet core with behavior: plc=%d nobisect=%d quality=%d minratio=%g verbose=%d fixedvolume=%d maxvolume=%g",
      b.plc, b.nobisect, b.quality, b.minratio, b.verbose, b.fixedvolume, b.maxvolume);
#endif
    try {
      run_tet_mesh_core(&b, &in, &out_tet);
    }
    catch(...) {
#ifdef SQMESH_HAS_SPDLOG
      SQMESH_LOG_ERROR("{}", "Exception caught during run_tet_mesh_core()");
#else
      SQMESH_LOG_ERROR("%s", "Exception caught during run_tet_mesh_core()");
#endif
      return fail_tet_core(
        base::StatusCode::internal_error,
        "run_tet_mesh_core() threw an exception"
      );
    }

#ifdef SQMESH_HAS_SPDLOG
    SQMESH_LOG_INFO(
      "Tet core produced: {} points, {} tetrahedra",
      out_tet.numberofpoints,
      out_tet.numberoftetrahedra
    );
#else
    SQMESH_LOG_INFO(
      "Tet core produced: %d points, %d tetrahedra",
      out_tet.numberofpoints,
      out_tet.numberoftetrahedra
    );
#endif

    if(out_tet.numberoftetrahedra == 0) {
#ifdef SQMESH_HAS_SPDLOG
      SQMESH_LOG_ERROR("{}", "Tet core failed to produce tetrahedra.");
#else
      SQMESH_LOG_ERROR("%s", "Tet core failed to produce tetrahedra.");
#endif
      return fail_tet_core(
        base::StatusCode::internal_error,
        "Tet core produced zero tetrahedra"
      );
    }

    EntityGroupDefinition node_def;
    node_def.order = EntityOrder::node;
    node_def.name = "volume_nodes";
    node_def.role = EntityGroupRole::computational;
    auto node_entity_group = output.create_entity_group(node_def);

    std::vector<EntityRef> vertex_map(out_tet.numberofpoints);
    for(int i = 0; i < out_tet.numberofpoints; ++i) {
      std::array<double, 3> pos = {
        out_tet.pointlist[i * 3 + 0],
        out_tet.pointlist[i * 3 + 1],
        out_tet.pointlist[i * 3 + 2]
      };
      vertex_map[i] = output.add_node(node_entity_group, pos);
    }

    struct Array3Less {
      bool operator()(
        const std::array<int, 3> &a,
        const std::array<int, 3> &b
      ) const {
        if(a[0] != b[0]) return a[0] < b[0];
        if(a[1] != b[1]) return a[1] < b[1];
        return a[2] < b[2];
      }
    };

    std::map<std::array<int, 3>, int, Array3Less> boundary_markers;
    if(out_tet.trifacelist != nullptr && out_tet.trifacemarkerlist != nullptr) {
      for(int i = 0; i < out_tet.numberoftrifaces; ++i) {
        std::array<int, 3> key = {
          out_tet.trifacelist[i * 3 + 0],
          out_tet.trifacelist[i * 3 + 1],
          out_tet.trifacelist[i * 3 + 2]
        };
        std::sort(key.begin(), key.end());
        const int marker = out_tet.trifacemarkerlist[i];
        if(marker != 0) {
          boundary_markers[key] = marker;
        }
      }
    } else {
      for(const auto &f : input_faces) {
        std::array<int, 3> key = {f.v[0], f.v[1], f.v[2]};
        std::sort(key.begin(), key.end());
        boundary_markers[key] = f.marker;
      }
    }

    EntityGroupDefinition int_face_def;
    int_face_def.order = EntityOrder::face;
    int_face_def.name = "interior_faces";
    int_face_def.default_kind = EntityKind::face_triangle;
    int_face_def.semantic = EntityGroupSemantic::interior;
    int_face_def.role = EntityGroupRole::computational;
    auto interior_face_entity_group = output.create_entity_group(int_face_def);

    std::map<int, EntityGroupIndex> pid_face_entity_groups;
    for(const auto &pair : boundary_markers) {
      const int marker = pair.second;
      if(pid_face_entity_groups.find(marker) == pid_face_entity_groups.end()) {
        EntityGroupDefinition def;
        def.order = EntityOrder::face;
        def.name = "boundary_face_" + std::to_string(marker);
        def.zone_id = marker;
        def.source_entity_tag = marker;
        def.boundary = true;
        def.default_kind = EntityKind::face_triangle;
        def.semantic = EntityGroupSemantic::boundary;
        def.role = EntityGroupRole::computational;
        pid_face_entity_groups[marker] = output.create_entity_group(def);
      }
    }

    EntityGroupDefinition cell_def;
    cell_def.order = EntityOrder::cell;
    cell_def.name = "volume_cells";
    cell_def.default_kind = EntityKind::cell_tetra;
    cell_def.role = EntityGroupRole::computational;
    auto cell_entity_group = output.create_entity_group(cell_def);

    std::map<std::array<int, 3>, EntityRef, Array3Less> created_faces;
    constexpr int face_table[4][3] = {
      {1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}
    };

    int kept_tets = 0;
    int skipped_tets = 0;
    for(int i = 0; i < out_tet.numberoftetrahedra; ++i) {
      if(use_regions && out_tet.tetrahedronattributelist != nullptr &&
         out_tet.numberoftetrahedronattributes > 0) {
        const int attr = static_cast<int>(
          out_tet.tetrahedronattributelist[
            i * out_tet.numberoftetrahedronattributes
          ]
        );
        if(attr < 1 || attr > static_cast<int>(material_points_.size())) {
          ++skipped_tets;
          continue;
        }
      }
      ++kept_tets;

      int *t = &out_tet.tetrahedronlist[i * 4];
      std::array<EntityRef, 4> nrefs = {
        vertex_map[t[0]], vertex_map[t[1]], vertex_map[t[2]], vertex_map[t[3]]
      };

      std::array<EntityRef, 4> frefs;
      for(int fi = 0; fi < 4; ++fi) {
        std::array<int, 3> key = {
          t[face_table[fi][0]], t[face_table[fi][1]], t[face_table[fi][2]]
        };
        std::array<EntityRef, 3> fn = {
          vertex_map[key[0]], vertex_map[key[1]], vertex_map[key[2]]
        };
        std::sort(key.begin(), key.end());

        auto it = created_faces.find(key);
        if(it != created_faces.end()) {
          frefs[fi] = it->second;
        } else {
          auto b_it = boundary_markers.find(key);
          if(b_it != boundary_markers.end()) {
            frefs[fi] = output.add_triangle_face(
              pid_face_entity_groups[b_it->second], fn
            );
          } else {
            frefs[fi] = output.add_triangle_face(interior_face_entity_group, fn);
          }
          created_faces[key] = frefs[fi];
        }
      }

      static_cast<void>(output.add_tetra_cell(cell_entity_group, nrefs, frefs));
    }

    if(use_regions) {
      SQMESH_LOG_INFO(
        "Tet region filter: kept {} tets, dropped {} tets outside material regions",
        kept_tets, skipped_tets
      );
    }

#ifdef SQMESH_HAS_SPDLOG
    SQMESH_LOG_INFO(
      "Successfully injected {} nodes and {} cells into Domain",
      out_tet.numberofpoints,
      use_regions ? kept_tets : out_tet.numberoftetrahedra
    );
#else
    SQMESH_LOG_INFO(
      "Successfully injected %d nodes and %d cells into Domain",
      out_tet.numberofpoints,
      use_regions ? kept_tets : out_tet.numberoftetrahedra
    );
#endif

    return base::StatusCode::ok;
  }
};

MeshingAlgorithmPtr create_native_tet_volume_mesher()
{
  return std::make_unique<NativeTetVolumeAlgorithm>();
}

} // namespace sqmesh::mesh::detail
