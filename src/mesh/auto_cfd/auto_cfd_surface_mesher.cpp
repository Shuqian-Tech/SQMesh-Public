// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "auto_cfd_surface_mesher.hpp"

#include "auto_cfd_surface_pipeline.hpp"
#include "auto_cfd_surface_test_hook.hpp"

#include "core/log.hpp"
#include "core/predicates.hpp"
#include "core/runtime_registry.hpp"

#include <algorithm>
#include <chrono>
#include <functional>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <set>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace sqmesh::mesh::detail {
namespace {

constexpr std::string_view kAutoCfdSurfaceMesherName = "Auto CFD Surface Mesher";
constexpr double kPi = 3.14159265358979323846264338327950288;
constexpr double kProjectionToleranceFloor = 1.0e-7;
constexpr double kUvAreaTolerance = 1.0e-12;
constexpr double kMetricTolerance = 1.0e-9;
constexpr double kMetricAcceptLength = 1.35;
constexpr double kMetricAcceptAspectRatio = 1.8;
constexpr double kMetricAcceptMinAngleDegrees = 32.0;
constexpr double kMetricAcceptMaxAngleDegrees = 118.0;
constexpr double kInternalQualityGuardrailMinimumAngleDegrees = 1.5;
constexpr double kInternalQualityGuardrailMaximumAngleDegrees = 175.0;
constexpr double kInternalQualityGuardrailMaximumAspectRatio = 40.0;
constexpr double kInternalQualityGuardrailMinimumRadiusRatio = 0.005;
constexpr double kInternalQualityGuardrailMaximumSkewness = 0.98;
constexpr double kDeliveredQualityMinimumAngleDegrees = 2.6;
constexpr double kDeliveredQualityMaximumAngleDegrees = 170.0;
constexpr double kDeliveredQualityMaximumAspectRatio = 20.5;
constexpr double kDeliveredQualityMinimumRadiusRatio = 0.018;
constexpr double kDeliveredQualityMaximumSkewness = 0.96;
constexpr double kFrontIdealHeight = 0.86602540378443864676372317075293618;
constexpr double kFrontRadiusLimit = 0.7071067811865476;
constexpr std::size_t kMaximumInsertedVerticesPerFace = 100000U;
constexpr std::size_t kMaximumFrontIterationsPerFace = 500000U;
constexpr std::size_t kMaximumRepairPasses = 48U;


struct LocalEdgeKey final {
  std::uint32_t first = invalid_index;
  std::uint32_t second = invalid_index;
};

[[nodiscard]] bool operator==(LocalEdgeKey lhs, LocalEdgeKey rhs) noexcept
{
  return lhs.first == rhs.first && lhs.second == rhs.second;
}

struct LocalEdgeKeyHash final {
  [[nodiscard]] std::size_t operator()(LocalEdgeKey key) const noexcept
  {
    return (static_cast<std::size_t>(key.first) << 32U) ^ key.second;
  }
};

struct MeshEdgeKey final {
  std::uint64_t first = 0U;
  std::uint64_t second = 0U;
};

[[nodiscard]] bool operator==(MeshEdgeKey lhs, MeshEdgeKey rhs) noexcept
{
  return lhs.first == rhs.first && lhs.second == rhs.second;
}

struct MeshEdgeKeyHash final {
  [[nodiscard]] std::size_t operator()(MeshEdgeKey key) const noexcept
  {
    return static_cast<std::size_t>(key.first ^ (key.second << 1U));
  }
};

struct MeshEdgeIncidence final {
  std::array<EntityRef, 2> nodes {};
  EntityRef left_face {};
  EntityRef right_face {};
  geo::TopologyEntityId topology_owner {};
  bool topology_owner_conflicted = false;
};

struct MeshFaceKey final {
  std::array<std::uint64_t, 3> nodes {0U, 0U, 0U};
};

[[nodiscard]] bool operator==(const MeshFaceKey &lhs, const MeshFaceKey &rhs) noexcept
{
  return lhs.nodes == rhs.nodes;
}

struct MeshFaceKeyHash final {
  [[nodiscard]] std::size_t operator()(const MeshFaceKey &key) const noexcept
  {
    return static_cast<std::size_t>(
      key.nodes[0] ^ (key.nodes[1] << 1U) ^ (key.nodes[2] << 2U)
    );
  }
};

struct LocalVertex final {
  std::uint32_t boundary_node_index = invalid_index;
  geo::Point3 position {0.0, 0.0, 0.0};
  std::array<double, 2> uv {0.0, 0.0};
  bool boundary = false;
  double local_size = 0.0;
  double background_size = 0.0;
  geo::TopologyEntityId boundary_edge_owner {};
  std::uint32_t periodic_counterpart = invalid_index;
};

struct LocalTriangle final {
  std::array<std::uint32_t, 3> vertices {invalid_index, invalid_index, invalid_index};
  bool valid = false;
  bool accepted = false;
  double front_radius = 0.0;
  double front_radius_override = std::numeric_limits<double>::quiet_NaN();
};

struct SeedLoopData final {
  geo::FaceBoundaryLoopKind kind = geo::FaceBoundaryLoopKind::unknown;
  std::vector<std::uint32_t> vertices {};
  std::vector<std::array<double, 2>> uvs {};
  std::vector<geo::TopologyEntityId> edge_owners {};
};

// Reverse the winding of a closed loop while keeping the starting vertex.
// vertices[i]/uvs[i] mirror; edge_owners[i] owns segment i → (i+1)%N, so the
// new segment i runs old-vertex N-1-i → N-2-i = old segment N-2-i; that is,
// reverse the first N-1 owners while leaving the wrap-around owner in place.
void reverse_seed_loop(SeedLoopData &loop)
{
  std::reverse(loop.vertices.begin(), loop.vertices.end());
  std::reverse(loop.uvs.begin(), loop.uvs.end());
  if(loop.edge_owners.size() >= 2U) {
    std::reverse(loop.edge_owners.begin(), loop.edge_owners.end() - 1);
  }
}

// Forward declaration: defined near other math helpers.
[[nodiscard]] double uv_loop_winding_sum(
  const std::vector<std::array<double, 2>> &polygon
) noexcept;
void enforce_outer_inner_loop_winding(
  SeedLoopData &outer_loop,
  std::vector<SeedLoopData> &inner_loops
)
{
  if(uv_loop_winding_sum(outer_loop.uvs) < 0.0) {
    reverse_seed_loop(outer_loop);
  }
  for(auto &inner : inner_loops) {
    if(uv_loop_winding_sum(inner.uvs) > 0.0) {
      reverse_seed_loop(inner);
    }
  }
}

void normalize_seed_loop_cycle(SeedLoopData &loop)
{
  if(loop.vertices.empty()) {
    loop.uvs.clear();
    loop.edge_owners.clear();
    return;
  }

  std::vector<std::uint32_t> normalized_vertices;
  std::vector<std::array<double, 2>> normalized_uvs;
  normalized_vertices.reserve(loop.vertices.size());
  normalized_uvs.reserve(loop.uvs.size());

  for(std::size_t index = 0U; index < loop.vertices.size(); ++index) {
    const auto vertex = loop.vertices[index];
    if(!normalized_vertices.empty() &&
       normalized_vertices.back() == vertex) {
      continue;
    }

    normalized_vertices.push_back(vertex);
    if(index < loop.uvs.size()) {
      normalized_uvs.push_back(loop.uvs[index]);
    }
  }

  if(normalized_vertices.size() >= 2U &&
     normalized_vertices.front() == normalized_vertices.back()) {
    normalized_vertices.pop_back();
    if(!normalized_uvs.empty()) {
      normalized_uvs.pop_back();
    }
  }

  loop.vertices = std::move(normalized_vertices);
  loop.uvs = std::move(normalized_uvs);
  loop.edge_owners.clear();
}

[[nodiscard]] bool boundary_loop_samples_are_coincident(
  const geo::Point3 &lhs_position,
  const std::array<double, 2> &lhs_uv,
  const geo::Point3 &rhs_position,
  const std::array<double, 2> &rhs_uv
) noexcept
{
  constexpr double kPositionTolerance = 1.0e-6;
  constexpr double kUvTolerance = 1.0e-8;
  for(std::size_t axis = 0U; axis < 3U; ++axis) {
    if(std::abs(lhs_position[axis] - rhs_position[axis]) > kPositionTolerance) {
      return false;
    }
  }
  for(std::size_t axis = 0U; axis < 2U; ++axis) {
    if(std::abs(lhs_uv[axis] - rhs_uv[axis]) > kUvTolerance) {
      return false;
    }
  }
  return true;
}

[[nodiscard]] bool boundary_loop_samples_are_coincident(
  const LocalVertex &vertex,
  const AutoCfdSurfaceFaceLoopPoint &point
) noexcept
{
  return boundary_loop_samples_are_coincident(
    vertex.position,
    vertex.uv,
    point.position,
    point.uv
  );
}

using BoundaryNodeRepresentativeMap =
  std::unordered_map<std::uint32_t, std::uint32_t>;

[[nodiscard]] BoundaryNodeRepresentativeMap build_boundary_node_representatives(
  const std::vector<AutoCfdSurfaceFaceLoopSegment> &segments
)
{
  struct RepresentativeSample final {
    std::uint32_t node_index = invalid_index;
    geo::Point3 position {0.0, 0.0, 0.0};
    std::array<double, 2> uv {0.0, 0.0};
  };

  BoundaryNodeRepresentativeMap representative_by_node;
  std::vector<RepresentativeSample> representatives;
  for(const auto &segment : segments) {
    for(const auto &point : segment.points) {
      if(point.node_index == invalid_index || !point.uv_defined) {
        continue;
      }

      auto rep_it = representative_by_node.find(point.node_index);
      if(rep_it != representative_by_node.end()) {
        continue;
      }

      auto existing_rep = std::find_if(
        representatives.begin(),
        representatives.end(),
        [&](const RepresentativeSample &sample) noexcept {
          return boundary_loop_samples_are_coincident(
            sample.position,
            sample.uv,
            point.position,
            point.uv
          );
        }
      );
      if(existing_rep != representatives.end()) {
        representative_by_node.emplace(point.node_index, existing_rep->node_index);
        continue;
      }

      representative_by_node.emplace(point.node_index, point.node_index);
      representatives.push_back({point.node_index, point.position, point.uv});
    }
  }

  return representative_by_node;
}

using LocalEdgeOwnerMap =
  std::unordered_map<LocalEdgeKey, geo::TopologyEntityId, LocalEdgeKeyHash>;

struct MetricSampleTrace final {
  bool used_reference_metric = false;
  bool used_planar_fallback = false;
  base::StatusCode direct_status = base::StatusCode::ok;
  AutoCfdSurfaceMetricFallbackKind direct_fallback_kind =
    AutoCfdSurfaceMetricFallbackKind::none;
};

struct MetricSamplingDiagnostics final {
  std::size_t total_calls = 0U;
  std::size_t planar_fallback_calls = 0U;
  std::size_t direct_metric_calls = 0U;
  std::size_t direct_internal_fallback_calls = 0U;
  std::size_t reference_metric_fallback_calls = 0U;
  std::size_t reference_metric_status_failure_calls = 0U;
  std::size_t reference_metric_unusable_calls = 0U;
};

struct TriangleMetricQuality final {
  bool valid = false;
  std::array<double, 3> edge_lengths {0.0, 0.0, 0.0};
  double min_angle_degrees = 0.0;
  double max_angle_degrees = 180.0;
  double aspect_ratio = std::numeric_limits<double>::infinity();
};

struct QualityGateFailures final {
  bool minimum_angle = false;
  bool maximum_angle = false;
  bool aspect_ratio = false;
  bool radius_ratio = false;
  bool skewness = false;

  [[nodiscard]] bool any() const noexcept
  {
    return minimum_angle || maximum_angle || aspect_ratio || radius_ratio || skewness;
  }
};

struct FrontEdgeCandidate final {
  std::uint32_t triangle_index = invalid_index;
  std::uint8_t edge_slot = 0U;
  double priority = std::numeric_limits<double>::infinity();
};

[[nodiscard]] double effective_triangle_front_radius(
  const LocalTriangle &triangle
) noexcept
{
  if(!triangle.valid) {
    return -std::numeric_limits<double>::infinity();
  }
  if(std::isfinite(triangle.front_radius_override)) {
    return triangle.front_radius_override;
  }
  return triangle.front_radius;
}

struct FrontTrianglePriorityComparator final {
  const std::vector<LocalTriangle> *triangles = nullptr;

  [[nodiscard]] bool operator()(
    std::uint32_t lhs,
    std::uint32_t rhs
  ) const noexcept
  {
    if(lhs == rhs) {
      return false;
    }

    const double lhs_radius =
      triangles != nullptr && lhs < triangles->size()
        ? effective_triangle_front_radius((*triangles)[lhs])
        : -std::numeric_limits<double>::infinity();
    const double rhs_radius =
      triangles != nullptr && rhs < triangles->size()
        ? effective_triangle_front_radius((*triangles)[rhs])
        : -std::numeric_limits<double>::infinity();

    if(lhs_radius > rhs_radius) {
      return true;
    }
    if(lhs_radius < rhs_radius) {
      return false;
    }
    return lhs < rhs;
  }
};

using FrontTriangleSet = std::set<std::uint32_t, FrontTrianglePriorityComparator>;

struct FaceMeshPlan final {
  geo::TopologyEntityId face {};
  std::vector<LocalVertex> vertices {};
  std::vector<LocalTriangle> triangles {};
  std::unordered_set<LocalEdgeKey, LocalEdgeKeyHash> constrained_edges {};
  LocalEdgeOwnerMap constrained_edge_owners {};
  std::array<double, 2> reference_uv {0.0, 0.0};
  geo::Vector3 reference_normal {0.0, 0.0, 0.0};
  bool reference_normal_defined = false;
  bool planar_fallback = false;
  geo::Point3 planar_origin {0.0, 0.0, 0.0};
  geo::Vector3 planar_u {1.0, 0.0, 0.0};
  geo::Vector3 planar_v {0.0, 1.0, 0.0};
  double target_size = 0.0;
  std::size_t seed_vertex_count = 0U;
  std::size_t seed_triangle_count = 0U;
  std::size_t inserted_vertex_count = 0U;
  std::size_t front_iteration_count = 0U;
  bool inserted_vertex_limit_reached = false;
  bool front_iteration_limit_reached = false;

  
  std::vector<std::array<double, 2>> periodic_true_boundary_segments {};
  std::array<double, 2> periodic_true_boundary_far {0.0, 0.0};
  bool periodic_true_boundary_defined = false;

  mutable MetricSamplingDiagnostics metric_sampling {};
  testing::AutoCfdSurfaceMeshingInstabilityDiagnostics instability {};
};

enum class LocalTriangleSourceFaceFailureKind : std::uint8_t {
  none = 0,
  node_reprojection_failure = 1,
  source_face_containment_failure = 2,
  orientation_flip = 3,
};

void record_source_face_failure(
  testing::AutoCfdSurfaceSourceFaceFailureCounts &counts,
  LocalTriangleSourceFaceFailureKind failure
) noexcept
{
  switch(failure) {
  case LocalTriangleSourceFaceFailureKind::none:
    break;
  case LocalTriangleSourceFaceFailureKind::node_reprojection_failure:
    ++counts.node_reprojection_failure_count;
    break;
  case LocalTriangleSourceFaceFailureKind::source_face_containment_failure:
    ++counts.source_face_containment_failure_count;
    break;
  case LocalTriangleSourceFaceFailureKind::orientation_flip:
    ++counts.orientation_flip_count;
    break;
  }
}

[[nodiscard]] bool triangle_passes_source_face_screen(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  const FaceMeshPlan &plan,
  const std::array<std::uint32_t, 3> &triangle,
  LocalTriangleSourceFaceFailureKind *failure
);

[[nodiscard]] base::StatusCode inject_face_meshes(
  const geo::ModelView &model_view,
  const AutoCfdSurfacePipelineState &pipeline_state,
  const std::vector<FaceMeshPlan> &face_meshes,
  Domain &output
);

[[nodiscard]] bool should_flip_triangle_winding(
  const geo::FaceView &face_view,
  const FaceMeshPlan &plan,
  const std::array<std::uint32_t, 3> &triangle
) noexcept;

[[nodiscard]] std::uint8_t find_active_edge(
  const FaceMeshPlan &plan,
  const std::vector<std::array<std::uint32_t, 3>> &neighbors,
  std::uint32_t triangle_index,
  const std::unordered_set<LocalEdgeKey, LocalEdgeKeyHash> *front_edges
) noexcept;

void enqueue_front_layer_candidates(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceSizingFieldState &field,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  const FaceMeshPlan &plan,
  const std::vector<std::array<std::uint32_t, 3>> &neighbors,
  const std::unordered_set<LocalEdgeKey, LocalEdgeKeyHash> &front_edges,
  std::uint32_t begin_triangle_index,
  std::uint32_t end_triangle_index,
  std::unordered_set<std::uint32_t> &queued_triangles,
  std::vector<FrontEdgeCandidate> &candidates
);

void accept_all_valid_triangles(FaceMeshPlan &plan) noexcept;

void align_uv_to_reference_periodically(
  const geo::FaceUvBounds &bounds,
  const std::array<double, 2> &reference,
  std::array<double, 2> &uv
) noexcept;

void wrap_uv_to_face_bounds(
  const geo::FaceUvBounds &bounds,
  std::array<double, 2> &uv
) noexcept;

[[nodiscard]] bool uses_supported_seam_material_screen(
  const AutoCfdSurfaceFacePreprocessState &face_preprocess
) noexcept;


[[nodiscard]] bool build_boundary_cdt_with_holes(
  const std::vector<LocalVertex> &vertices,
  const SeedLoopData &outer_loop,
  const std::vector<SeedLoopData> &inner_loops,
  std::vector<std::array<std::uint32_t, 3>> &out_triangles,
  std::unordered_set<std::uint64_t> *out_not_recovered = nullptr
);

[[nodiscard]] double projection_tolerance(double target_size) noexcept
{
  return std::max(kProjectionToleranceFloor, target_size * 1.0e-6);
}

[[nodiscard]] LocalEdgeKey canonical_edge(
  std::uint32_t first,
  std::uint32_t second
) noexcept
{
  if(second < first) {
    std::swap(first, second);
  }
  return {first, second};
}

[[nodiscard]] bool edge_keys_share_vertex(
  LocalEdgeKey lhs,
  LocalEdgeKey rhs
) noexcept
{
  return lhs.first == rhs.first ||
         lhs.first == rhs.second ||
         lhs.second == rhs.first ||
         lhs.second == rhs.second;
}

[[nodiscard]] bool triangle_contains_edge(
  const LocalTriangle &triangle,
  LocalEdgeKey edge
) noexcept
{
  if(!triangle.valid) {
    return false;
  }

  for(std::uint8_t edge_slot = 0U; edge_slot < 3U; ++edge_slot) {
    if(canonical_edge(
         triangle.vertices[edge_slot],
         triangle.vertices[(edge_slot + 1U) % 3U]
       ) == edge) {
      return true;
    }
  }

  return false;
}

[[nodiscard]] std::uint8_t find_triangle_edge_slot(
  const LocalTriangle &triangle,
  LocalEdgeKey edge
) noexcept
{
  if(!triangle.valid) {
    return 3U;
  }

  for(std::uint8_t edge_slot = 0U; edge_slot < 3U; ++edge_slot) {
    if(canonical_edge(
         triangle.vertices[edge_slot],
         triangle.vertices[(edge_slot + 1U) % 3U]
       ) == edge) {
      return edge_slot;
    }
  }

  return 3U;
}

[[nodiscard]] std::uint64_t pack_entity_ref(EntityRef ref) noexcept
{
  return (static_cast<std::uint64_t>(ref.entity_group) << 32U) |
         static_cast<std::uint64_t>(ref.index);
}

[[nodiscard]] std::uint64_t pack_topology_entity_id(
  geo::TopologyEntityId entity
) noexcept
{
  return (static_cast<std::uint64_t>(entity.dimension) << 32U) |
         static_cast<std::uint64_t>(entity.index);
}

[[nodiscard]] std::string topology_entity_debug_label(
  geo::TopologyEntityId entity
)
{
  if(!geo::is_valid(entity)) {
    return "invalid";
  }

  std::string label;
  switch(entity.dimension) {
  case geo::TopologyDimension::vertex:
    label = "vertex";
    break;
  case geo::TopologyDimension::edge:
    label = "edge";
    break;
  case geo::TopologyDimension::face:
    label = "face";
    break;
  case geo::TopologyDimension::region:
    label = "region";
    break;
  default:
    label = "unknown";
    break;
  }
  label += ":";
  label += std::to_string(entity.index);
  return label;
}

void log_constrained_edge_owner_breakdown(
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  const FaceMeshPlan &plan,
  const std::vector<LocalEdgeKey> &edges,
  std::string_view stage
)
{
  if(edges.empty()) {
    return;
  }

  struct OwnerBucket final {
    geo::TopologyEntityId owner {};
    std::size_t count = 0U;
    LocalEdgeKey sample {invalid_index, invalid_index};
  };

  std::unordered_map<std::uint64_t, OwnerBucket> buckets;
  buckets.reserve(edges.size());
  constexpr std::uint64_t kInvalidOwnerKey = std::numeric_limits<std::uint64_t>::max();

  for(const auto edge : edges) {
    const auto owner_it = plan.constrained_edge_owners.find(edge);
    const auto owner =
      owner_it != plan.constrained_edge_owners.end()
        ? owner_it->second
        : geo::TopologyEntityId {};
    const auto bucket_key =
      geo::is_valid(owner) ? pack_topology_entity_id(owner) : kInvalidOwnerKey;
    auto [it, inserted] = buckets.emplace(bucket_key, OwnerBucket {});
    if(inserted) {
      it->second.owner = owner;
      it->second.sample = edge;
    }
    ++it->second.count;
  }

  std::vector<OwnerBucket> ordered_buckets;
  ordered_buckets.reserve(buckets.size());
  for(const auto &[key, bucket] : buckets) {
    static_cast<void>(key);
    ordered_buckets.push_back(bucket);
  }
  std::sort(
    ordered_buckets.begin(),
    ordered_buckets.end(),
    [](const OwnerBucket &lhs, const OwnerBucket &rhs) noexcept {
      if(lhs.count != rhs.count) {
        return lhs.count > rhs.count;
      }
      const auto lhs_owner =
        geo::is_valid(lhs.owner) ? pack_topology_entity_id(lhs.owner) : kInvalidOwnerKey;
      const auto rhs_owner =
        geo::is_valid(rhs.owner) ? pack_topology_entity_id(rhs.owner) : kInvalidOwnerKey;
      return lhs_owner < rhs_owner;
    }
  );

  std::string summary;
  for(std::size_t i = 0U; i < ordered_buckets.size(); ++i) {
    if(i > 0U) {
      summary += ", ";
    }
    summary += topology_entity_debug_label(ordered_buckets[i].owner);
    summary += "(x";
    summary += std::to_string(ordered_buckets[i].count);
    summary += ")";
  }
  SQMESH_LOG_WARN(
    "{} on {}: {} unrecovered constrained edge(s) across {} owner(s) — {}",
    stage,
    topology_entity_debug_label(face_preprocess.face),
    edges.size(),
    ordered_buckets.size(),
    summary);
}

void log_open_mesh_edge_owner_breakdown(
  const Domain &domain
)
{
  struct OwnerBucket final {
    geo::TopologyEntityId owner {};
    std::size_t count = 0U;
    EntityRef sample_edge {};
    geo::TopologyEntityId sample_face_owner {};
  };

  std::unordered_map<std::uint64_t, OwnerBucket> buckets;
  buckets.reserve(128U);
  constexpr std::uint64_t kInvalidOwnerKey = std::numeric_limits<std::uint64_t>::max();
  std::size_t total_open_edge_count = 0U;

  for(const auto &entity_group : domain.entity_groups()) {
    if(entity_group.order() != EntityOrder::edge ||
       entity_group.role() != EntityGroupRole::computational) {
      continue;
    }

    for(std::uint32_t edge_index = 0U; edge_index < entity_group.edges().size(); ++edge_index) {
      const EntityRef edge_ref {entity_group.id(), edge_index};
      const auto left_face = domain.adjacent_face(edge_ref, FaceSide::left);
      const auto right_face = domain.adjacent_face(edge_ref, FaceSide::right);
      if(is_valid(left_face) && is_valid(right_face)) {
        continue;
      }

      ++total_open_edge_count;
      const auto owner = domain.edge_topology_owner(edge_ref);
      const auto bucket_key =
        geo::is_valid(owner) ? pack_topology_entity_id(owner) : kInvalidOwnerKey;
      auto [it, inserted] = buckets.emplace(bucket_key, OwnerBucket {});
      if(inserted) {
        it->second.owner = owner;
        it->second.sample_edge = edge_ref;
        const auto face_ref = is_valid(left_face) ? left_face : right_face;
        it->second.sample_face_owner =
          is_valid(face_ref) ? domain.face_topology_owner(face_ref) : geo::TopologyEntityId {};
      }
      ++it->second.count;
    }
  }

  if(total_open_edge_count == 0U) {
    return;
  }

  std::vector<OwnerBucket> ordered_buckets;
  ordered_buckets.reserve(buckets.size());
  for(const auto &[key, bucket] : buckets) {
    static_cast<void>(key);
    ordered_buckets.push_back(bucket);
  }
  std::sort(
    ordered_buckets.begin(),
    ordered_buckets.end(),
    [](const OwnerBucket &lhs, const OwnerBucket &rhs) noexcept {
      if(lhs.count != rhs.count) {
        return lhs.count > rhs.count;
      }
      const auto lhs_owner =
        geo::is_valid(lhs.owner) ? pack_topology_entity_id(lhs.owner) : kInvalidOwnerKey;
      const auto rhs_owner =
        geo::is_valid(rhs.owner) ? pack_topology_entity_id(rhs.owner) : kInvalidOwnerKey;
      return lhs_owner < rhs_owner;
    }
  );

}

void log_plan_open_edge_breakdown(
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  const FaceMeshPlan &plan,
  std::string_view stage
)
{
  struct Bucket final {
    bool constrained = false;
    geo::TopologyEntityId owner {};
    std::size_t count = 0U;
    LocalEdgeKey sample {invalid_index, invalid_index};
  };

  std::unordered_map<LocalEdgeKey, std::size_t, LocalEdgeKeyHash> edge_counts;
  edge_counts.reserve(plan.triangles.size() * 3U);
  for(const auto &triangle : plan.triangles) {
    if(!triangle.valid || !triangle.accepted) {
      continue;
    }
    for(std::uint8_t edge_slot = 0U; edge_slot < 3U; ++edge_slot) {
      ++edge_counts[canonical_edge(
        triangle.vertices[edge_slot],
        triangle.vertices[(edge_slot + 1U) % 3U]
      )];
    }
  }

  std::size_t total_open_edge_count = 0U;
  std::vector<Bucket> buckets;
  for(const auto &[edge, count] : edge_counts) {
    if(count != 1U) {
      continue;
    }
    ++total_open_edge_count;

    const bool constrained =
      plan.constrained_edges.find(edge) != plan.constrained_edges.end();
    geo::TopologyEntityId owner {};
    const auto owner_it = plan.constrained_edge_owners.find(edge);
    if(owner_it != plan.constrained_edge_owners.end()) {
      owner = owner_it->second;
    }

    auto bucket_it = std::find_if(
      buckets.begin(),
      buckets.end(),
      [&](const Bucket &bucket) noexcept {
        return bucket.constrained == constrained && bucket.owner == owner;
      }
    );
    if(bucket_it == buckets.end()) {
      buckets.push_back({constrained, owner, 1U, edge});
    }
    else {
      ++bucket_it->count;
    }
  }

  if(total_open_edge_count == 0U) {
    return;
  }

  std::sort(
    buckets.begin(),
    buckets.end(),
    [](const Bucket &lhs, const Bucket &rhs) noexcept {
      if(lhs.constrained != rhs.constrained) {
        return lhs.constrained > rhs.constrained;
      }
      if(lhs.count != rhs.count) {
        return lhs.count > rhs.count;
      }
      constexpr std::uint64_t kInvalidOwnerKey = std::numeric_limits<std::uint64_t>::max();
      const auto lhs_owner =
        geo::is_valid(lhs.owner) ? pack_topology_entity_id(lhs.owner) : kInvalidOwnerKey;
      const auto rhs_owner =
        geo::is_valid(rhs.owner) ? pack_topology_entity_id(rhs.owner) : kInvalidOwnerKey;
      return lhs_owner < rhs_owner;
    }
  );

}

[[nodiscard]] MeshFaceKey canonical_face(
  EntityRef first,
  EntityRef second,
  EntityRef third
) noexcept
{
  MeshFaceKey key {
    {
      pack_entity_ref(first),
      pack_entity_ref(second),
      pack_entity_ref(third),
    }
  };
  std::sort(key.nodes.begin(), key.nodes.end());
  return key;
}

[[nodiscard]] testing::AutoCfdSurfaceQualityGate
internal_quality_guardrail() noexcept
{
  return {
    kInternalQualityGuardrailMinimumAngleDegrees,
    kInternalQualityGuardrailMaximumAngleDegrees,
    kInternalQualityGuardrailMaximumAspectRatio,
    kInternalQualityGuardrailMinimumRadiusRatio,
    kInternalQualityGuardrailMaximumSkewness,
  };
}

[[nodiscard]] testing::AutoCfdSurfaceQualityGate delivered_quality_gate() noexcept
{
  return {
    kDeliveredQualityMinimumAngleDegrees,
    kDeliveredQualityMaximumAngleDegrees,
    kDeliveredQualityMaximumAspectRatio,
    kDeliveredQualityMinimumRadiusRatio,
    kDeliveredQualityMaximumSkewness,
  };
}

[[nodiscard]] QualityGateFailures evaluate_quality_gate_failures(
  const ElementQuality &quality,
  const testing::AutoCfdSurfaceQualityGate &gate
) noexcept
{
  return {
    !std::isfinite(quality.min_angle) ||
      quality.min_angle < gate.minimum_min_angle_degrees,
    !std::isfinite(quality.max_angle) ||
      quality.max_angle > gate.maximum_max_angle_degrees,
    !std::isfinite(quality.aspect_ratio) ||
      quality.aspect_ratio > gate.maximum_aspect_ratio,
    !std::isfinite(quality.radius_ratio) ||
      quality.radius_ratio < gate.minimum_radius_ratio,
    !std::isfinite(quality.skewness) ||
      quality.skewness > gate.maximum_skewness,
  };
}

[[nodiscard]] std::array<double, 2> uv_subtract(
  const std::array<double, 2> &lhs,
  const std::array<double, 2> &rhs
) noexcept
{
  return {lhs[0] - rhs[0], lhs[1] - rhs[1]};
}

[[nodiscard]] std::array<double, 2> uv_add(
  const std::array<double, 2> &lhs,
  const std::array<double, 2> &rhs
) noexcept
{
  return {lhs[0] + rhs[0], lhs[1] + rhs[1]};
}

[[nodiscard]] std::array<double, 2> uv_scale(
  const std::array<double, 2> &value,
  double scale
) noexcept
{
  return {value[0] * scale, value[1] * scale};
}

[[nodiscard]] double uv_dot(
  const std::array<double, 2> &lhs,
  const std::array<double, 2> &rhs
) noexcept
{
  return lhs[0] * rhs[0] + lhs[1] * rhs[1];
}

[[nodiscard]] double uv_cross(
  const std::array<double, 2> &lhs,
  const std::array<double, 2> &rhs
) noexcept
{
  return lhs[0] * rhs[1] - lhs[1] * rhs[0];
}

[[nodiscard]] double uv_norm(const std::array<double, 2> &value) noexcept
{
  return std::sqrt(uv_dot(value, value));
}

[[nodiscard]] std::array<double, 2> uv_left_normal(
  const std::array<double, 2> &value
) noexcept
{
  return {-value[1], value[0]};
}

[[nodiscard]] std::array<double, 2> uv_normalized(
  const std::array<double, 2> &value
) noexcept
{
  const double length = uv_norm(value);
  if(length <= 0.0 || !std::isfinite(length)) {
    return {0.0, 0.0};
  }
  return {value[0] / length, value[1] / length};
}

[[nodiscard]] geo::Point3 point_subtract(
  const geo::Point3 &lhs,
  const geo::Point3 &rhs
) noexcept
{
  return {
    lhs[0] - rhs[0],
    lhs[1] - rhs[1],
    lhs[2] - rhs[2],
  };
}

[[nodiscard]] geo::Vector3 cross_product(
  const geo::Point3 &lhs,
  const geo::Point3 &rhs
) noexcept
{
  return {
    lhs[1] * rhs[2] - lhs[2] * rhs[1],
    lhs[2] * rhs[0] - lhs[0] * rhs[2],
    lhs[0] * rhs[1] - lhs[1] * rhs[0],
  };
}

[[nodiscard]] double dot_product(
  const geo::Vector3 &lhs,
  const geo::Vector3 &rhs
) noexcept
{
  return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2];
}

[[nodiscard]] geo::Vector3 normalized_vector(
  const geo::Vector3 &value
) noexcept
{
  const double length = std::sqrt(dot_product(value, value));
  if(length <= 0.0 || !std::isfinite(length)) {
    return {0.0, 0.0, 0.0};
  }

  return {
    value[0] / length,
    value[1] / length,
    value[2] / length,
  };
}

[[nodiscard]] geo::Point3 point_add_scaled(
  const geo::Point3 &origin,
  const geo::Vector3 &direction,
  double scale
) noexcept
{
  return {
    origin[0] + direction[0] * scale,
    origin[1] + direction[1] * scale,
    origin[2] + direction[2] * scale,
  };
}

[[nodiscard]] double signed_triangle_area_twice(
  const std::array<double, 2> &a,
  const std::array<double, 2> &b,
  const std::array<double, 2> &c
) noexcept
{
  return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
}

[[nodiscard]] double uv_loop_winding_sum(
  const std::vector<std::array<double, 2>> &polygon
) noexcept
{
  double winding_sum = 0.0;
  for(std::size_t index = 0U; index < polygon.size(); ++index) {
    const auto &current = polygon[index];
    const auto &previous = polygon[(index + polygon.size() - 1U) % polygon.size()];
    winding_sum += (previous[0] - current[0]) * (current[1] + previous[1]);
  }
  return winding_sum;
}

[[nodiscard]] bool orient_triangle_ccw(
  const std::vector<LocalVertex> &vertices,
  std::array<std::uint32_t, 3> &triangle
) noexcept
{
  if(triangle[0] >= vertices.size() ||
     triangle[1] >= vertices.size() ||
     triangle[2] >= vertices.size()) {
    return false;
  }

  if(std::abs(signed_triangle_area_twice(
       vertices[triangle[0]].uv,
       vertices[triangle[1]].uv,
       vertices[triangle[2]].uv
     )) <= kUvAreaTolerance) {
    return false;
  }

  if(signed_triangle_area_twice(
       vertices[triangle[0]].uv,
       vertices[triangle[1]].uv,
       vertices[triangle[2]].uv
     ) < 0.0) {
    std::swap(triangle[1], triangle[2]);
  }
  return true;
}

void record_mesh_edge(
  std::unordered_map<MeshEdgeKey, MeshEdgeIncidence, MeshEdgeKeyHash> &mesh_edges,
  EntityRef a,
  EntityRef b,
  EntityRef face_ref,
  geo::TopologyEntityId topology_owner = {}
)
{
  const std::uint64_t packed_a =
    (static_cast<std::uint64_t>(a.entity_group) << 32U) | a.index;
  const std::uint64_t packed_b =
    (static_cast<std::uint64_t>(b.entity_group) << 32U) | b.index;
  const MeshEdgeKey key {
    std::min(packed_a, packed_b),
    std::max(packed_a, packed_b),
  };

  auto [it, inserted] = mesh_edges.emplace(
    key,
    MeshEdgeIncidence {{a, b}, face_ref, {}, topology_owner, false}
  );
  if(inserted) {
    return;
  }

  if(!is_valid(it->second.right_face)) {
    it->second.right_face = face_ref;
  }

  if(geo::is_valid(topology_owner)) {
    if(!geo::is_valid(it->second.topology_owner)) {
      it->second.topology_owner = topology_owner;
    }
    else if(it->second.topology_owner != topology_owner) {
      it->second.topology_owner = {};
      it->second.topology_owner_conflicted = true;
    }
  }
}

[[nodiscard]] std::size_t dominant_uv_dimension(
  const std::vector<std::array<double, 2>> &lhs,
  const std::vector<std::array<double, 2>> &rhs
) noexcept
{
  std::array<double, 2> minima {
    std::numeric_limits<double>::infinity(),
    std::numeric_limits<double>::infinity(),
  };
  std::array<double, 2> maxima {
    -std::numeric_limits<double>::infinity(),
    -std::numeric_limits<double>::infinity(),
  };

  const auto accumulate = [&minima, &maxima](const std::vector<std::array<double, 2>> &uvs) {
    for(const auto &uv : uvs) {
      minima[0] = std::min(minima[0], uv[0]);
      minima[1] = std::min(minima[1], uv[1]);
      maxima[0] = std::max(maxima[0], uv[0]);
      maxima[1] = std::max(maxima[1], uv[1]);
    }
  };

  accumulate(lhs);
  accumulate(rhs);

  const double span0 = maxima[0] - minima[0];
  const double span1 = maxima[1] - minima[1];
  return span1 > span0 ? 1U : 0U;
}

[[nodiscard]] bool build_planar_fallback_basis(
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  FaceMeshPlan &plan
) noexcept
{
  const AutoCfdSurfaceFaceLoopPoint *origin_point = nullptr;
  const AutoCfdSurfaceFaceLoopPoint *u_point = nullptr;
  const AutoCfdSurfaceFaceLoopPoint *v_point = nullptr;

  for(const auto &loop : face_preprocess.loops) {
    if(loop.points.size() < 3U) {
      continue;
    }

    origin_point = &loop.points[0];
    for(std::size_t point_index = 1U; point_index < loop.points.size(); ++point_index) {
      const auto edge = point_subtract(loop.points[point_index].position, origin_point->position);
      if(std::sqrt(dot_product(edge, edge)) > kMetricTolerance) {
        u_point = &loop.points[point_index];
        break;
      }
    }
    if(u_point == nullptr) {
      continue;
    }

    for(std::size_t point_index = 1U; point_index < loop.points.size(); ++point_index) {
      const auto edge = point_subtract(loop.points[point_index].position, origin_point->position);
      const auto normal = cross_product(
        point_subtract(u_point->position, origin_point->position),
        edge
      );
      if(std::sqrt(dot_product(normal, normal)) > kMetricTolerance) {
        v_point = &loop.points[point_index];
        break;
      }
    }
    if(v_point != nullptr) {
      break;
    }
  }

  if(origin_point == nullptr || u_point == nullptr || v_point == nullptr) {
    return false;
  }

  const auto raw_u = point_subtract(u_point->position, origin_point->position);
  const auto raw_normal = cross_product(
    raw_u,
    point_subtract(v_point->position, origin_point->position)
  );
  const auto normal = normalized_vector(raw_normal);
  if(std::sqrt(dot_product(normal, normal)) <= kMetricTolerance) {
    return false;
  }

  const auto u_axis = normalized_vector(raw_u);
  const auto v_axis = normalized_vector(cross_product(normal, u_axis));
  if(std::sqrt(dot_product(v_axis, v_axis)) <= kMetricTolerance) {
    return false;
  }

  plan.planar_fallback = true;
  plan.planar_origin = origin_point->position;
  plan.planar_u = u_axis;
  plan.planar_v = v_axis;
  plan.reference_normal = normal;
  plan.reference_normal_defined = true;
  return true;
}

[[nodiscard]] std::array<double, 2> planar_local_uv(
  const FaceMeshPlan &plan,
  const geo::Point3 &position
) noexcept
{
  const auto delta = point_subtract(position, plan.planar_origin);
  return {
    dot_product(delta, plan.planar_u),
    dot_product(delta, plan.planar_v),
  };
}

[[nodiscard]] geo::Point3 planar_world_point(
  const FaceMeshPlan &plan,
  const std::array<double, 2> &uv
) noexcept
{
  auto point = point_add_scaled(plan.planar_origin, plan.planar_u, uv[0]);
  point = point_add_scaled(point, plan.planar_v, uv[1]);
  return point;
}

[[nodiscard]] double angle_from_edges(
  double lhs,
  double rhs,
  double opposite
) noexcept
{
  const double denominator = 2.0 * lhs * rhs;
  if(denominator <= 0.0) {
    return 0.0;
  }

  const double cosine =
    std::clamp((lhs * lhs + rhs * rhs - opposite * opposite) / denominator, -1.0, 1.0);
  return std::acos(cosine);
}

[[nodiscard]] base::StatusCode sample_metric_or_reference(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceSizingFieldState &field,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  const FaceMeshPlan *plan,
  double u,
  double v,
  AutoCfdSurfaceFaceMetricTensor &metric,
  MetricSampleTrace *trace = nullptr
)
{
  if(trace != nullptr) {
    *trace = {};
  }
  if(plan != nullptr) {
    ++plan->metric_sampling.total_calls;
  }

  if(plan != nullptr && plan->planar_fallback) {
    metric = {};
    metric.face = face_view.entity;
    metric.u = u;
    metric.v = v;
    const auto position = planar_world_point(*plan, {u, v});
    metric.target_size =
      query_auto_cfd_surface_sizing_field(field, face_view.entity, position);
    const double inverse_target_size_squared =
      1.0 / std::max(metric.target_size * metric.target_size, kMetricTolerance);
    metric.tensor = {
      inverse_target_size_squared,
      0.0,
      inverse_target_size_squared,
    };
    metric.clamped_eigenvalues = {
      inverse_target_size_squared,
      inverse_target_size_squared,
    };
    metric.determinant = inverse_target_size_squared * inverse_target_size_squared;
    metric.usable = true;
    metric.fallback_kind = AutoCfdSurfaceMetricFallbackKind::isotropic_from_uv_identity;
    if(trace != nullptr) {
      trace->used_planar_fallback = true;
      trace->direct_fallback_kind = metric.fallback_kind;
    }
    if(plan != nullptr) {
      ++plan->metric_sampling.planar_fallback_calls;
    }
    return core::detail::clear_error_state();
  }

  auto status = sample_auto_cfd_surface_face_metric(face_view, field, u, v, metric);
  if(uses_supported_seam_material_screen(face_preprocess)) {
    auto wrapped_uv = std::array<double, 2> {u, v};
    wrap_uv_to_face_bounds(face_preprocess.uv_bounds, wrapped_uv);
    status = sample_auto_cfd_surface_face_metric(
      face_view,
      field,
      wrapped_uv[0],
      wrapped_uv[1],
      metric
    );
  }
  if(status == base::StatusCode::ok && metric.usable) {
    metric.u = u;
    metric.v = v;
    if(trace != nullptr) {
      trace->direct_status = status;
      trace->direct_fallback_kind = metric.fallback_kind;
    }
    if(plan != nullptr) {
      ++plan->metric_sampling.direct_metric_calls;
      if(metric.fallback_kind != AutoCfdSurfaceMetricFallbackKind::none) {
        ++plan->metric_sampling.direct_internal_fallback_calls;
      }
    }
    return status;
  }

  if(trace != nullptr) {
    trace->used_reference_metric = true;
    trace->direct_status = status;
    trace->direct_fallback_kind = metric.fallback_kind;
  }
  if(plan != nullptr) {
    ++plan->metric_sampling.reference_metric_fallback_calls;
    if(status != base::StatusCode::ok) {
      ++plan->metric_sampling.reference_metric_status_failure_calls;
    }
    else {
      ++plan->metric_sampling.reference_metric_unusable_calls;
    }
  }

  metric = face_preprocess.reference_metric;
  metric.u = u;
  metric.v = v;
  if(metric.usable) {
    return core::detail::clear_error_state();
  }

  metric = {};
  metric.face = face_view.entity;
  metric.u = u;
  metric.v = v;
  metric.usable = true;
  metric.tensor = {1.0, 0.0, 1.0};
  metric.clamped_eigenvalues = {1.0, 1.0};
  metric.determinant = 1.0;
  metric.target_size = field.maximum_length;
  return core::detail::clear_error_state();
}

void log_metric_sampling_summary(
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  const FaceMeshPlan &plan,
  std::string_view stage
)
{
  const auto &diag = plan.metric_sampling;
}


[[nodiscard]] bool circumcenter_metric(
  const std::array<double, 2> &pa,
  const std::array<double, 2> &pb,
  const std::array<double, 2> &pc,
  const std::array<double, 3> &metric,
  std::array<double, 2> &center,
  double &radius_squared
) noexcept;

[[nodiscard]] bool point_in_circumcircle_aniso(
  const std::array<double, 2> &a,
  const std::array<double, 2> &b,
  const std::array<double, 2> &c,
  const std::array<double, 2> &p,
  const std::array<double, 3> &metric
) noexcept;

[[nodiscard]] std::array<double, 3> select_surface_metric_tensor(
  const AutoCfdSurfaceFaceMetricTensor &metric
) noexcept
{
  if(metric.first_fundamental_form_defined) {
    return metric.first_fundamental_form;
  }
  if(std::isfinite(metric.target_size) &&
     metric.target_size > kMetricTolerance) {
    const double scale = metric.target_size * metric.target_size;
    return {
      metric.tensor[0] * scale,
      metric.tensor[1] * scale,
      metric.tensor[2] * scale,
    };
  }
  return {1.0, 0.0, 1.0};
}

[[nodiscard]] double raw_metric_distance_squared(
  const std::array<double, 3> &metric,
  const std::array<double, 2> &delta
) noexcept
{
  return metric[0] * delta[0] * delta[0] +
         metric[2] * delta[1] * delta[1] +
         2.0 * metric[1] * delta[0] * delta[1];
}

[[nodiscard]] TriangleMetricQuality evaluate_triangle_metric_quality(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceSizingFieldState &field,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  const FaceMeshPlan &plan,
  const std::array<std::uint32_t, 3> &triangle
)
{
  TriangleMetricQuality quality;
  if(triangle[0] >= plan.vertices.size() ||
     triangle[1] >= plan.vertices.size() ||
     triangle[2] >= plan.vertices.size()) {
    return quality;
  }

  const auto &a = plan.vertices[triangle[0]].uv;
  const auto &b = plan.vertices[triangle[1]].uv;
  const auto &c = plan.vertices[triangle[2]].uv;
  const std::array<double, 2> centroid {
    (a[0] + b[0] + c[0]) / 3.0,
    (a[1] + b[1] + c[1]) / 3.0,
  };

  AutoCfdSurfaceFaceMetricTensor metric;
  static_cast<void>(sample_metric_or_reference(
    face_view,
    field,
    face_preprocess,
    &plan,
    centroid[0],
    centroid[1],
    metric
  ));

  
  const auto metric_tensor = select_surface_metric_tensor(metric);
  if(!std::isfinite(metric_tensor[0]) ||
     !std::isfinite(metric_tensor[1]) ||
     !std::isfinite(metric_tensor[2]) ||
     metric_tensor[0] <= kMetricTolerance ||
     metric_tensor[2] <= kMetricTolerance) {
    return quality;
  }
  const double metric_det =
    metric_tensor[0] * metric_tensor[2] - metric_tensor[1] * metric_tensor[1];
  if(!std::isfinite(metric_det) || metric_det <= kMetricTolerance) {
    return quality;
  }

  const auto ab_uv = uv_subtract(b, a);
  const auto ac_uv = uv_subtract(c, a);
  const auto bc_uv = uv_subtract(c, b);

  const double l01 = std::sqrt(raw_metric_distance_squared(metric_tensor, ab_uv));
  const double l12 = std::sqrt(raw_metric_distance_squared(metric_tensor, bc_uv));
  const double l20 = std::sqrt(raw_metric_distance_squared(metric_tensor, ac_uv));
  const double longest = std::max({l01, l12, l20});
  // Metric area: dA_M = sqrt(det(M)) du dv, so the metric triangle area is
  // 0.5 * sqrt(det(M)) * |Euclidean UV cross of (b-a, c-a)|.
  const double area = 0.5 * std::sqrt(metric_det) * std::abs(uv_cross(ab_uv, ac_uv));
  if(longest <= kMetricTolerance || area <= kMetricTolerance) {
    return quality;
  }

  const double shortest_altitude = (2.0 * area) / longest;
  if(shortest_altitude <= kMetricTolerance) {
    return quality;
  }

  const double angle0 = angle_from_edges(l01, l20, l12);
  const double angle1 = angle_from_edges(l01, l12, l20);
  const double angle2 = angle_from_edges(l12, l20, l01);

  quality.valid = true;
  quality.edge_lengths = {l01, l12, l20};
  quality.min_angle_degrees =
    std::min({angle0, angle1, angle2}) * (180.0 / kPi);
  quality.max_angle_degrees =
    std::max({angle0, angle1, angle2}) * (180.0 / kPi);
  quality.aspect_ratio =
    0.86602540378443864676372317075293618 * longest / shortest_altitude;
  return quality;
}

[[nodiscard]] bool triangle_is_good_enough(
  const TriangleMetricQuality &quality
) noexcept
{
  if(!quality.valid) {
    return false;
  }

  const double longest = std::max({
    quality.edge_lengths[0],
    quality.edge_lengths[1],
    quality.edge_lengths[2],
  });
  return longest <= kMetricAcceptLength &&
         quality.aspect_ratio <= kMetricAcceptAspectRatio &&
         quality.min_angle_degrees >= kMetricAcceptMinAngleDegrees &&
         quality.max_angle_degrees <= kMetricAcceptMaxAngleDegrees;
}

[[nodiscard]] double triangle_metric_badness(
  const TriangleMetricQuality &quality
) noexcept
{
  if(!quality.valid) {
    return std::numeric_limits<double>::infinity();
  }

  const double longest = std::max({
    quality.edge_lengths[0],
    quality.edge_lengths[1],
    quality.edge_lengths[2],
  });
  const double longest_penalty = std::max(0.0, longest - kMetricAcceptLength);
  const double aspect_penalty =
    std::max(0.0, quality.aspect_ratio - kMetricAcceptAspectRatio);
  const double min_angle_penalty =
    std::max(0.0, kMetricAcceptMinAngleDegrees - quality.min_angle_degrees);
  const double max_angle_penalty =
    std::max(0.0, quality.max_angle_degrees - kMetricAcceptMaxAngleDegrees);
  return longest_penalty +
         aspect_penalty +
         (min_angle_penalty / kMetricAcceptMinAngleDegrees) +
         (max_angle_penalty / kMetricAcceptMaxAngleDegrees);
}

[[nodiscard]] bool point_in_triangle(
  const std::array<double, 2> &point,
  const std::array<double, 2> &a,
  const std::array<double, 2> &b,
  const std::array<double, 2> &c
) noexcept
{
  const double area0 = signed_triangle_area_twice(a, b, point);
  const double area1 = signed_triangle_area_twice(b, c, point);
  const double area2 = signed_triangle_area_twice(c, a, point);
  const bool has_negative =
    area0 < -kUvAreaTolerance || area1 < -kUvAreaTolerance || area2 < -kUvAreaTolerance;
  const bool has_positive =
    area0 > kUvAreaTolerance || area1 > kUvAreaTolerance || area2 > kUvAreaTolerance;
  return !(has_negative && has_positive);
}

enum class UvLoopContainmentResult : std::uint8_t {
  outside = 0,
  inside = 1,
  on_boundary = 2,
};

struct SupportedSeamStripBoundaryData final {
  std::array<std::vector<std::array<double, 2>>, 2> boundaries {};
  std::size_t periodic_dimension = 0U;
  double periodic_span = 0.0;
  bool valid = false;
};

[[nodiscard]] bool uv_point_on_segment(
  const std::array<double, 2> &point,
  const std::array<double, 2> &a,
  const std::array<double, 2> &b
) noexcept
{
  const auto ab = uv_subtract(b, a);
  const auto ap = uv_subtract(point, a);
  if(std::abs(uv_cross(ab, ap)) > kUvAreaTolerance) {
    return false;
  }

  const double dot = uv_dot(ap, ab);
  if(dot < -kUvAreaTolerance) {
    return false;
  }

  const double length_squared = uv_dot(ab, ab);
  return dot <= length_squared + kUvAreaTolerance;
}

[[nodiscard]] bool uv_segments_intersect(
  const std::array<double, 2> &a,
  const std::array<double, 2> &b,
  const std::array<double, 2> &c,
  const std::array<double, 2> &d
) noexcept
{
  const auto orientation_sign =
    [](double value) noexcept {
      if(value > kUvAreaTolerance) {
        return 1;
      }
      if(value < -kUvAreaTolerance) {
        return -1;
      }
      return 0;
    };

  const auto o1 = orientation_sign(signed_triangle_area_twice(a, b, c));
  const auto o2 = orientation_sign(signed_triangle_area_twice(a, b, d));
  const auto o3 = orientation_sign(signed_triangle_area_twice(c, d, a));
  const auto o4 = orientation_sign(signed_triangle_area_twice(c, d, b));

  if(o1 != o2 && o3 != o4) {
    return true;
  }
  if(o1 == 0 && uv_point_on_segment(c, a, b)) {
    return true;
  }
  if(o2 == 0 && uv_point_on_segment(d, a, b)) {
    return true;
  }
  if(o3 == 0 && uv_point_on_segment(a, c, d)) {
    return true;
  }
  if(o4 == 0 && uv_point_on_segment(b, c, d)) {
    return true;
  }

  return false;
}

void align_uv_to_reference_periodically(
  const geo::FaceUvBounds &bounds,
  const std::array<double, 2> &reference,
  std::array<double, 2> &uv
) noexcept
{
  const std::array<double, 2> spans {
    bounds.u_max - bounds.u_min,
    bounds.v_max - bounds.v_min,
  };

  for(std::size_t dimension = 0U; dimension < spans.size(); ++dimension) {
    const double span = spans[dimension];
    if(!std::isfinite(span) || span <= kUvAreaTolerance) {
      continue;
    }

    uv[dimension] -=
      std::round((uv[dimension] - reference[dimension]) / span) * span;
  }
}

void wrap_uv_to_face_bounds(
  const geo::FaceUvBounds &bounds,
  std::array<double, 2> &uv
) noexcept
{
  const std::array<double, 2> minima {
    bounds.u_min,
    bounds.v_min,
  };
  const std::array<double, 2> maxima {
    bounds.u_max,
    bounds.v_max,
  };

  for(std::size_t dimension = 0U; dimension < minima.size(); ++dimension) {
    const double span = maxima[dimension] - minima[dimension];
    if(!std::isfinite(span) || span <= kUvAreaTolerance) {
      continue;
    }

    if(uv[dimension] < minima[dimension] - kUvAreaTolerance ||
       uv[dimension] > maxima[dimension] + kUvAreaTolerance) {
      uv[dimension] -=
        std::floor((uv[dimension] - minima[dimension]) / span) * span;
    }
    if(uv[dimension] > maxima[dimension] &&
       uv[dimension] - maxima[dimension] <= kUvAreaTolerance) {
      uv[dimension] = maxima[dimension];
    }
  }
}

[[nodiscard]] UvLoopContainmentResult classify_uv_point_in_loop(
  const std::vector<AutoCfdSurfaceFaceLoopPoint> &points,
  const std::array<double, 2> &uv
) noexcept
{
  if(points.size() < 3U) {
    return UvLoopContainmentResult::outside;
  }

  bool inside = false;
  for(std::size_t index = 0U; index < points.size(); ++index) {
    const auto &current = points[index];
    const auto &next = points[(index + 1U) % points.size()];
    if(!current.uv_defined || !next.uv_defined) {
      continue;
    }

    if(uv_point_on_segment(uv, current.uv, next.uv)) {
      return UvLoopContainmentResult::on_boundary;
    }

    const bool intersects =
      ((current.uv[1] > uv[1]) != (next.uv[1] > uv[1])) &&
      (uv[0] <
       ((next.uv[0] - current.uv[0]) * (uv[1] - current.uv[1]) /
          (next.uv[1] - current.uv[1])) +
         current.uv[0]);
    if(intersects) {
      inside = !inside;
    }
  }

  return inside ? UvLoopContainmentResult::inside : UvLoopContainmentResult::outside;
}

[[nodiscard]] bool supports_uv_loop_containment(
  const AutoCfdSurfaceFacePreprocessState &face_preprocess
) noexcept
{
  return face_preprocess.uv_reconstruction_available &&
         !face_preprocess.boundary.has_seams &&
         face_preprocess.boundary.outer_loop_count == 1U &&
         face_preprocess.boundary.unknown_loop_count == 0U;
}

[[nodiscard]] bool uses_supported_seam_material_screen(
  const AutoCfdSurfaceFacePreprocessState &face_preprocess
) noexcept
{
  return face_preprocess.boundary.has_seams &&
         face_preprocess.seam_support ==
           AutoCfdSurfaceSeamSupportKind::supported_subset &&
         face_preprocess.seam_unwrap_applied &&
         face_preprocess.uv_bounds_defined &&
         face_preprocess.loops.size() == 1U;
}

[[nodiscard]] bool uses_uv_material_face_screen(
  const AutoCfdSurfaceFacePreprocessState &face_preprocess
) noexcept
{
  return uses_supported_seam_material_screen(face_preprocess) ||
         (supports_uv_loop_containment(face_preprocess) &&
          face_preprocess.boundary.inner_loop_count == 0U);
}

[[nodiscard]] bool face_has_nontrivial_curvature(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess
)
{
  std::array<std::array<double, 2>, 2> sample_uvs {{
    {
      face_preprocess.reference_metric.u,
      face_preprocess.reference_metric.v,
    },
    {
      0.5 * (face_preprocess.recovered_uv_min[0] + face_preprocess.recovered_uv_max[0]),
      0.5 * (face_preprocess.recovered_uv_min[1] + face_preprocess.recovered_uv_max[1]),
    },
  }};

  for(auto uv : sample_uvs) {
    if(!std::isfinite(uv[0]) || !std::isfinite(uv[1])) {
      continue;
    }

    geo::FaceCurvatureSample sample;
    if(geo::sample_face_curvature(face_view, uv[0], uv[1], sample) !=
         base::StatusCode::ok ||
       !sample.curvature_defined) {
      continue;
    }

    return std::max(
             std::abs(sample.min_curvature),
             std::abs(sample.max_curvature)
           ) > 1.0e-6;
  }

  return false;
}

[[nodiscard]] bool build_supported_seam_strip_boundary_data(
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  SupportedSeamStripBoundaryData &strip
) noexcept
{
  strip = {};
  if(!uses_supported_seam_material_screen(face_preprocess)) {
    return false;
  }

  std::array<std::vector<std::array<double, 2>>, 2> boundaries;
  std::size_t boundary_index = 0U;
  for(const auto &segment : face_preprocess.loops.front().segments) {
    if(segment.edge_use.is_seam) {
      continue;
    }
    if(boundary_index >= boundaries.size()) {
      return false;
    }

    auto &boundary = boundaries[boundary_index++];
    boundary.reserve(segment.points.size());
    for(const auto &point : segment.points) {
      if(!point.uv_defined) {
        return false;
      }
      boundary.push_back(point.uv);
    }
  }

  if(boundary_index != boundaries.size() ||
     boundaries[0].size() < 2U ||
     boundaries[0].size() != boundaries[1].size()) {
    return false;
  }

  strip.periodic_dimension =
    dominant_uv_dimension(boundaries[0], boundaries[1]);
  strip.periodic_span =
    strip.periodic_dimension == 0U
      ? face_preprocess.uv_bounds.u_max - face_preprocess.uv_bounds.u_min
      : face_preprocess.uv_bounds.v_max - face_preprocess.uv_bounds.v_min;
  if(!std::isfinite(strip.periodic_span) ||
     strip.periodic_span <= kUvAreaTolerance) {
    return false;
  }

  const double first_direction =
    boundaries[0].back()[strip.periodic_dimension] -
    boundaries[0].front()[strip.periodic_dimension];
  const double second_direction =
    boundaries[1].back()[strip.periodic_dimension] -
    boundaries[1].front()[strip.periodic_dimension];
  if(first_direction * second_direction < 0.0) {
    std::reverse(boundaries[1].begin(), boundaries[1].end());
  }

  const auto append_periodic_closure =
    [periodic_dimension = strip.periodic_dimension, periodic_span = strip.periodic_span](
      std::vector<std::array<double, 2>> &boundary
    ) noexcept {
      auto closing = boundary.front();
      double direction =
        boundary.back()[periodic_dimension] - boundary.front()[periodic_dimension];
      if(std::abs(direction) <= kUvAreaTolerance) {
        direction = periodic_span;
      }
      closing[periodic_dimension] += direction >= 0.0 ? periodic_span : -periodic_span;
      boundary.push_back(closing);
    };

  append_periodic_closure(boundaries[0]);
  append_periodic_closure(boundaries[1]);

  strip.boundaries = std::move(boundaries);
  strip.valid = true;
  return true;
}

[[nodiscard]] bool uv_point_inside_supported_seam_strip(
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  const std::array<double, 2> &uv
) noexcept
{
  SupportedSeamStripBoundaryData strip;
  if(!build_supported_seam_strip_boundary_data(face_preprocess, strip)) {
    return true;
  }

  for(int shift = -1; shift <= 1; ++shift) {
    auto shifted_uv = uv;
    shifted_uv[strip.periodic_dimension] +=
      static_cast<double>(shift) * strip.periodic_span;
    for(std::size_t index = 0U;
        index + 1U < strip.boundaries[0].size();
        ++index) {
      const auto &first_current = strip.boundaries[0][index];
      const auto &first_next = strip.boundaries[0][index + 1U];
      const auto &second_current = strip.boundaries[1][index];
      const auto &second_next = strip.boundaries[1][index + 1U];

      if(point_in_triangle(
           shifted_uv,
           first_current,
           second_current,
           first_next
         ) ||
         point_in_triangle(
           shifted_uv,
           first_next,
           second_current,
           second_next
         )) {
        return true;
      }
    }
  }

  return false;
}

[[nodiscard]] bool uv_point_inside_face_material(
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  const std::array<double, 2> &uv
) noexcept
{
  if(uses_supported_seam_material_screen(face_preprocess)) {
    return uv_point_inside_supported_seam_strip(face_preprocess, uv);
  }

  if(!supports_uv_loop_containment(face_preprocess)) {
    return true;
  }

  bool inside_outer = false;
  for(const auto &loop : face_preprocess.loops) {
    const auto containment = classify_uv_point_in_loop(loop.points, uv);
    switch(loop.kind) {
    case geo::FaceBoundaryLoopKind::outer:
      if(containment == UvLoopContainmentResult::inside ||
         containment == UvLoopContainmentResult::on_boundary) {
        inside_outer = true;
      }
      break;
    case geo::FaceBoundaryLoopKind::inner:
      if(containment == UvLoopContainmentResult::inside) {
        return false;
      }
      break;
    case geo::FaceBoundaryLoopKind::unknown:
      break;
    }
  }

  return inside_outer;
}


[[nodiscard]] bool point_inside_periodic_true_boundary(
  const std::vector<std::array<double, 2>> &bnd,
  const std::array<double, 2> &far,
  const std::array<double, 2> &p
) noexcept
{
  if(bnd.size() < 2U || (bnd.size() % 2U) != 0U) {
    return true;
  }

  int count = 0;
  for(std::size_t i = 0U; i + 1U < bnd.size(); i += 2U) {
    const auto &p1 = bnd[i];
    const auto &p2 = bnd[i + 1U];

    // orient2d(p1, p2, p) vs orient2d(p1, p2, far): opposite sign means
    // p and far are on opposite sides of the segment line.
    const double a1 = signed_triangle_area_twice(p1, p2, p);
    const double b1 = signed_triangle_area_twice(p1, p2, far);
    if(a1 * b1 >= 0.0) {
      continue;
    }

    // orient2d(p, far, p1) vs orient2d(p, far, p2): opposite sign means
    // p1 and p2 are on opposite sides of the ray from p to far. Combined
    // with the first test, the two line segments (p1,p2) and (p,far)
    // actually cross.
    const double a2 = signed_triangle_area_twice(p, far, p1);
    const double b2 = signed_triangle_area_twice(p, far, p2);
    if(a2 * b2 < 0.0) {
      ++count;
    }
  }

  return (count & 1) == 1;
}

[[nodiscard]] bool triangle_uv_samples_inside_face_material(
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  const std::array<std::array<double, 2>, 3> &triangle_uvs
) noexcept
{
  if(!uses_uv_material_face_screen(face_preprocess)) {
    return true;
  }

  for(const auto &triangle_uv : triangle_uvs) {
    if(!uv_point_inside_face_material(face_preprocess, triangle_uv)) {
      return false;
    }
  }

  if(!uses_supported_seam_material_screen(face_preprocess)) {
    return true;
  }

  for(std::size_t edge_index = 0U; edge_index < triangle_uvs.size(); ++edge_index) {
    const auto &edge_start = triangle_uvs[edge_index];
    const auto &edge_end = triangle_uvs[(edge_index + 1U) % triangle_uvs.size()];
    for(const double edge_fraction : {0.25, 0.5, 0.75}) {
      const std::array<double, 2> edge_sample {
        edge_start[0] + (edge_end[0] - edge_start[0]) * edge_fraction,
        edge_start[1] + (edge_end[1] - edge_start[1]) * edge_fraction,
      };
      if(!uv_point_inside_face_material(face_preprocess, edge_sample)) {
        return false;
      }
    }
  }

  return true;
}

[[nodiscard]] bool sample_face_at_material_uv(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  std::array<double, 2> uv,
  geo::FaceSample &sample
)
{
  if(uses_supported_seam_material_screen(face_preprocess)) {
    wrap_uv_to_face_bounds(face_preprocess.uv_bounds, uv);
  }

  return geo::sample_face(face_view, uv[0], uv[1], sample) ==
         base::StatusCode::ok;
}

[[nodiscard]] bool triangle_passes_uv_material_face_screen(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  const std::array<geo::Point3, 3> &triangle_points,
  const std::array<std::array<double, 2>, 3> &triangle_uvs,
  LocalTriangleSourceFaceFailureKind *failure
)
{
  const std::array<double, 2> centroid_uv {
    (triangle_uvs[0][0] + triangle_uvs[1][0] + triangle_uvs[2][0]) / 3.0,
    (triangle_uvs[0][1] + triangle_uvs[1][1] + triangle_uvs[2][1]) / 3.0,
  };
  if(!uv_point_inside_face_material(face_preprocess, centroid_uv) ||
     !triangle_uv_samples_inside_face_material(face_preprocess, triangle_uvs)) {
    if(failure != nullptr) {
      *failure = LocalTriangleSourceFaceFailureKind::source_face_containment_failure;
    }
    return false;
  }

  geo::FaceSample centroid_sample;
  if(!sample_face_at_material_uv(
       face_view,
       face_preprocess,
       centroid_uv,
       centroid_sample
     )) {
    if(failure != nullptr) {
      *failure = LocalTriangleSourceFaceFailureKind::source_face_containment_failure;
    }
    return false;
  }

  if(centroid_sample.normal_defined) {
    const auto ab = point_subtract(triangle_points[1], triangle_points[0]);
    const auto ac = point_subtract(triangle_points[2], triangle_points[0]);
    const auto normal = normalized_vector(cross_product(ab, ac));
    if(dot_product(normal, centroid_sample.normal) <= 0.0) {
      if(failure != nullptr) {
        *failure = LocalTriangleSourceFaceFailureKind::orientation_flip;
      }
      return false;
    }
  }

  return true;
}

[[nodiscard]] bool uv_point_strictly_inside_triangle(
  const std::array<double, 2> &point,
  const std::array<double, 2> &a,
  const std::array<double, 2> &b,
  const std::array<double, 2> &c
) noexcept
{
  return point_in_triangle(point, a, b, c) &&
         !uv_point_on_segment(point, a, b) &&
         !uv_point_on_segment(point, b, c) &&
         !uv_point_on_segment(point, c, a);
}

[[nodiscard]] bool triangle_overlaps_inner_loops(
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  const FaceMeshPlan &plan,
  const std::array<std::uint32_t, 3> &triangle,
  const std::array<std::array<double, 2>, 3> &triangle_uvs
) noexcept
{
  if(face_preprocess.boundary.inner_loop_count == 0U ||
     !supports_uv_loop_containment(face_preprocess)) {
    return false;
  }

  const auto edge_is_inner_loop_chord =
    [&](std::uint32_t first_vertex, std::uint32_t second_vertex) noexcept {
      if(first_vertex >= plan.vertices.size() || second_vertex >= plan.vertices.size()) {
        return false;
      }

      const auto &first = plan.vertices[first_vertex];
      const auto &second = plan.vertices[second_vertex];
      if(!first.boundary || !second.boundary ||
         first.boundary_node_index == invalid_index ||
         second.boundary_node_index == invalid_index ||
         first.boundary_node_index == second.boundary_node_index) {
        return false;
      }

      for(const auto &loop : face_preprocess.loops) {
        if(loop.kind != geo::FaceBoundaryLoopKind::inner || loop.points.size() < 3U) {
          continue;
        }

        std::vector<std::size_t> first_positions;
        std::vector<std::size_t> second_positions;
        first_positions.reserve(2U);
        second_positions.reserve(2U);
        for(std::size_t point_index = 0U; point_index < loop.points.size(); ++point_index) {
          if(loop.points[point_index].node_index == first.boundary_node_index) {
            first_positions.push_back(point_index);
          }
          if(loop.points[point_index].node_index == second.boundary_node_index) {
            second_positions.push_back(point_index);
          }
        }

        if(first_positions.empty() || second_positions.empty()) {
          continue;
        }

        const auto loop_size = loop.points.size();
        for(const auto first_position : first_positions) {
          for(const auto second_position : second_positions) {
            const bool adjacent =
              ((first_position + 1U) % loop_size) == second_position ||
              ((second_position + 1U) % loop_size) == first_position;
            if(!adjacent) {
              return true;
            }
          }
        }
      }

      return false;
    };

  for(const auto &loop : face_preprocess.loops) {
    if(loop.kind != geo::FaceBoundaryLoopKind::inner || loop.points.size() < 3U) {
      continue;
    }

    for(const auto &triangle_uv : triangle_uvs) {
      if(classify_uv_point_in_loop(loop.points, triangle_uv) ==
         UvLoopContainmentResult::inside) {
        return true;
      }
    }

    for(std::size_t edge_index = 0U; edge_index < triangle.size(); ++edge_index) {
      if(edge_is_inner_loop_chord(
           triangle[edge_index],
           triangle[(edge_index + 1U) % triangle.size()]
         )) {
        return true;
      }
    }

    for(std::size_t edge_index = 0U; edge_index < triangle_uvs.size(); ++edge_index) {
      const auto &edge_start = triangle_uvs[edge_index];
      const auto &edge_end = triangle_uvs[(edge_index + 1U) % triangle_uvs.size()];
      for(const double edge_fraction : {0.25, 0.5, 0.75}) {
        const std::array<double, 2> edge_sample {
          edge_start[0] + (edge_end[0] - edge_start[0]) * edge_fraction,
          edge_start[1] + (edge_end[1] - edge_start[1]) * edge_fraction,
        };
        if(classify_uv_point_in_loop(loop.points, edge_sample) ==
           UvLoopContainmentResult::inside) {
          return true;
        }
      }
    }

    for(const auto &loop_point : loop.points) {
      if(!loop_point.uv_defined) {
        continue;
      }
      if(uv_point_strictly_inside_triangle(
           loop_point.uv,
           triangle_uvs[0],
           triangle_uvs[1],
           triangle_uvs[2]
          )) {
        return true;
      }
    }
  }

  return false;
}

[[nodiscard]] bool resolve_face_reference_normal(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  geo::Vector3 &normal
) noexcept
{
  geo::FaceSample sample;
  if(geo::sample_face(
       face_view,
       face_preprocess.reference_metric.u,
       face_preprocess.reference_metric.v,
       sample
     ) == base::StatusCode::ok &&
     sample.normal_defined) {
    normal = sample.normal;
    return true;
  }

  for(const auto &loop : face_preprocess.loops) {
    if(loop.points.size() < 3U) {
      continue;
    }

    const auto ab = point_subtract(loop.points[1].position, loop.points[0].position);
    const auto ac = point_subtract(loop.points[2].position, loop.points[0].position);
    normal = cross_product(ab, ac);
    if(std::abs(dot_product(normal, normal)) > kMetricTolerance) {
      return true;
    }
  }

  normal = {0.0, 0.0, 0.0};
  return false;
}

void build_boundary_node_edge_owner_map(
  const std::vector<AutoCfdSurfaceFaceLoopSegment> &segments,
  LocalEdgeOwnerMap &owner_by_boundary_edge
)
{
  owner_by_boundary_edge.clear();
  const auto representative_by_node =
    build_boundary_node_representatives(segments);
  std::size_t estimated_edge_count = 0U;
  for(const auto &segment : segments) {
    if(segment.points.size() >= 2U) {
      estimated_edge_count += segment.points.size() - 1U;
    }
  }
  owner_by_boundary_edge.reserve(estimated_edge_count);

  for(const auto &segment : segments) {
    if(!geo::is_valid(segment.edge_use.edge) || segment.points.size() < 2U) {
      continue;
    }

    for(std::size_t point_index = 0U; point_index + 1U < segment.points.size();
        ++point_index) {
      auto first = segment.points[point_index].node_index;
      auto second = segment.points[point_index + 1U].node_index;
      const auto first_rep = representative_by_node.find(first);
      const auto second_rep = representative_by_node.find(second);
      if(first_rep != representative_by_node.end()) {
        first = first_rep->second;
      }
      if(second_rep != representative_by_node.end()) {
        second = second_rep->second;
      }
      if(first == invalid_index || second == invalid_index || first == second) {
        continue;
      }

      owner_by_boundary_edge[canonical_edge(first, second)] =
        segment.edge_use.edge;
    }
  }

  for(std::size_t segment_index = 0U; segment_index < segments.size(); ++segment_index) {
    const auto &segment = segments[segment_index];
    if(!geo::is_valid(segment.edge_use.edge) || segment.points.empty()) {
      continue;
    }

    const auto next_index = (segment_index + 1U) % segments.size();
    if(next_index >= segments.size()) {
      continue;
    }
    const auto &next_segment = segments[next_index];
    if(next_segment.points.empty()) {
      continue;
    }

    
    if(segment.edge_use.end_vertex != next_segment.edge_use.start_vertex) {
      continue;
    }

    auto first = segment.points.back().node_index;
    auto second = next_segment.points.front().node_index;
    const auto first_rep = representative_by_node.find(first);
    const auto second_rep = representative_by_node.find(second);
    if(first_rep != representative_by_node.end()) {
      first = first_rep->second;
    }
    if(second_rep != representative_by_node.end()) {
      second = second_rep->second;
    }
    if(first == invalid_index || second == invalid_index || first == second) {
      continue;
    }

    owner_by_boundary_edge[canonical_edge(first, second)] =
      segment.edge_use.edge;
  }
}

[[nodiscard]] geo::TopologyEntityId lookup_boundary_edge_owner(
  const LocalEdgeOwnerMap *owner_by_boundary_edge,
  std::uint32_t first_boundary_node,
  std::uint32_t second_boundary_node
) noexcept
{
  if(owner_by_boundary_edge == nullptr ||
     first_boundary_node == invalid_index ||
     second_boundary_node == invalid_index ||
     first_boundary_node == second_boundary_node) {
    return {};
  }

  const auto it = owner_by_boundary_edge->find(
    canonical_edge(first_boundary_node, second_boundary_node)
  );
  if(it == owner_by_boundary_edge->end()) {
    return {};
  }
  return it->second;
}

void populate_seed_loop_edge_owners(
  const std::vector<LocalVertex> &vertices,
  const LocalEdgeOwnerMap *owner_by_boundary_edge,
  SeedLoopData &loop
)
{
  loop.edge_owners.clear();
  loop.edge_owners.reserve(loop.vertices.size());
  if(loop.vertices.empty()) {
    return;
  }

  for(std::size_t index = 0U; index < loop.vertices.size(); ++index) {
    const auto current = loop.vertices[index];
    const auto next = loop.vertices[(index + 1U) % loop.vertices.size()];
    geo::TopologyEntityId owner {};
    if(current < vertices.size() && next < vertices.size()) {
      const auto &current_vertex = vertices[current];
      const auto &next_vertex = vertices[next];
      owner = lookup_boundary_edge_owner(
        owner_by_boundary_edge,
        current_vertex.boundary_node_index,
        next_vertex.boundary_node_index
      );
      if(!geo::is_valid(owner) && !current_vertex.boundary &&
         geo::is_valid(current_vertex.boundary_edge_owner)) {
        owner = current_vertex.boundary_edge_owner;
      }
      if(!geo::is_valid(owner) && !next_vertex.boundary &&
         geo::is_valid(next_vertex.boundary_edge_owner)) {
        owner = next_vertex.boundary_edge_owner;
      }
    }
    loop.edge_owners.push_back(owner);
  }
}

[[nodiscard]] base::StatusCode append_seed_loop_from_points(
  const std::vector<AutoCfdSurfaceFaceLoopPoint> &points,
  geo::FaceBoundaryLoopKind kind,
  std::unordered_map<std::uint32_t, std::uint32_t> &boundary_to_local,
  std::vector<LocalVertex> &vertices,
  SeedLoopData &loop,
  bool densify_segments = false,
  const LocalEdgeOwnerMap *owner_by_boundary_edge = nullptr
)
{
  loop = {};
  loop.kind = kind;

  loop.vertices.reserve(densify_segments ? points.size() * 2U : points.size());
  loop.uvs.reserve(densify_segments ? points.size() * 2U : points.size());
  for(std::size_t point_index = 0U; point_index < points.size(); ++point_index) {
    const auto &point = points[point_index];
    if(!point.uv_defined || point.node_index == invalid_index) {
      return core::detail::publish_error(
        base::StatusCode::unsupported,
        "Auto CFD Surface Mesher requires UV-defined recovered boundary points for face seeding."
      );
    }

    auto local_it = boundary_to_local.find(point.node_index);
    if(local_it == boundary_to_local.end()) {
      if(!loop.vertices.empty()) {
        const auto previous_local = loop.vertices.back();
        if(previous_local < vertices.size() &&
           boundary_loop_samples_are_coincident(vertices[previous_local], point)) {
          local_it = boundary_to_local.emplace(point.node_index, previous_local).first;
        }
      }

      if(local_it == boundary_to_local.end() &&
         point_index + 1U == points.size() &&
         !loop.vertices.empty()) {
        const auto first_local = loop.vertices.front();
        if(first_local < vertices.size() &&
           boundary_loop_samples_are_coincident(vertices[first_local], point)) {
          local_it = boundary_to_local.emplace(point.node_index, first_local).first;
        }
      }

      if(local_it == boundary_to_local.end()) {
        const auto local_index = static_cast<std::uint32_t>(vertices.size());
        vertices.push_back(
          {
            point.node_index,
            point.position,
            point.uv,
            true,
          }
        );
        local_it = boundary_to_local.emplace(point.node_index, local_index).first;
      }
    }

    if(!loop.vertices.empty() && loop.vertices.back() == local_it->second) {
      continue;
    }

    loop.vertices.push_back(local_it->second);
    loop.uvs.push_back(point.uv);

    if(!densify_segments) {
      continue;
    }

    const auto &next_point = points[(point_index + 1U) % points.size()];
    const auto segment_owner = lookup_boundary_edge_owner(
      owner_by_boundary_edge,
      point.node_index,
      next_point.node_index
    );
    if(!next_point.uv_defined ||
       point.position == next_point.position ||
       point.uv == next_point.uv) {
      continue;
    }

    const auto midpoint_index = static_cast<std::uint32_t>(vertices.size());
    vertices.push_back(
      {
        invalid_index,
        {
          0.5 * (point.position[0] + next_point.position[0]),
          0.5 * (point.position[1] + next_point.position[1]),
          0.5 * (point.position[2] + next_point.position[2]),
        },
        {
          0.5 * (point.uv[0] + next_point.uv[0]),
          0.5 * (point.uv[1] + next_point.uv[1]),
        },
        false,
        0.0,
        0.0,
        segment_owner,
      }
    );
    loop.vertices.push_back(midpoint_index);
    loop.uvs.push_back(vertices[midpoint_index].uv);
  }

  if(loop.vertices.size() >= 2U && loop.vertices.front() == loop.vertices.back()) {
    loop.vertices.pop_back();
    loop.uvs.pop_back();
  }

  if(loop.vertices.size() < 3U) {
    return core::detail::publish_error(
      base::StatusCode::unsupported,
      "Auto CFD Surface Mesher requires each supported face loop to resolve to at least three unique boundary samples."
    );
  }

  populate_seed_loop_edge_owners(vertices, owner_by_boundary_edge, loop);

  return core::detail::clear_error_state();
}

void add_constrained_loop_edges(
  const SeedLoopData &loop,
  std::unordered_set<LocalEdgeKey, LocalEdgeKeyHash> &constrained_edges,
  LocalEdgeOwnerMap *constrained_edge_owners = nullptr
)
{
  for(std::size_t index = 0U; index < loop.vertices.size(); ++index) {
    const auto edge = canonical_edge(
      loop.vertices[index],
      loop.vertices[(index + 1U) % loop.vertices.size()]
    );
    constrained_edges.insert(edge);
    if(constrained_edge_owners != nullptr &&
       index < loop.edge_owners.size() &&
       geo::is_valid(loop.edge_owners[index])) {
      auto [it, inserted] = constrained_edge_owners->emplace(
        edge,
        loop.edge_owners[index]
      );
      if(!inserted && it->second != loop.edge_owners[index]) {
        it->second = {};
      }
    }
  }
}

[[nodiscard]] base::StatusCode build_planar_seed_mesh(
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  FaceMeshPlan &plan
)
{
  if(face_preprocess.boundary.outer_loop_count != 1U ||
     face_preprocess.boundary.unknown_loop_count != 0U) {
    return core::detail::publish_error(
      base::StatusCode::unsupported,
      "Auto CFD Surface Mesher currently supports only faces with one classified outer loop and zero unknown loops."
    );
  }

  std::unordered_map<std::uint32_t, std::uint32_t> boundary_to_local;
  boundary_to_local.reserve(face_preprocess.loops.size() * 8U);

  std::vector<SeedLoopData> seed_loops;
  seed_loops.reserve(face_preprocess.loops.size());
  for(std::size_t loop_index = 0U; loop_index < face_preprocess.loops.size(); ++loop_index) {
    const auto &loop = face_preprocess.loops[loop_index];
    if(!loop.closed || !loop.continuous) {
      return core::detail::publish_error(
        base::StatusCode::unsupported,
        "Auto CFD Surface Mesher requires closed, continuous face trim loops."
      );
    }

    SeedLoopData seed_loop;
    LocalEdgeOwnerMap boundary_edge_owners;
    build_boundary_node_edge_owner_map(loop.segments, boundary_edge_owners);
    auto status = append_seed_loop_from_points(
      loop.points,
      face_preprocess.boundary.loops[loop_index].kind,
      boundary_to_local,
      plan.vertices,
      seed_loop,
      false,  // densify_segments: OFF for CDT path — midpoints break shared boundary node matching
      &boundary_edge_owners
    );
    if(status != base::StatusCode::ok) {
      return status;
    }
    add_constrained_loop_edges(
      seed_loop,
      plan.constrained_edges,
      &plan.constrained_edge_owners
    );
    seed_loops.push_back(std::move(seed_loop));
  }

  if(face_preprocess.boundary.primary_outer_loop_index >= seed_loops.size()) {
    return core::detail::publish_error(
      base::StatusCode::internal_error,
      "Auto CFD Surface Mesher could not resolve the primary outer loop from face preprocessing state."
    );
  }

  auto outer_loop = seed_loops[face_preprocess.boundary.primary_outer_loop_index];
  std::vector<SeedLoopData> inner_loops;
  inner_loops.reserve(seed_loops.size());
  for(std::size_t loop_index = 0U; loop_index < seed_loops.size(); ++loop_index) {
    if(loop_index == face_preprocess.boundary.primary_outer_loop_index) {
      continue;
    }
    if(seed_loops[loop_index].kind != geo::FaceBoundaryLoopKind::inner) {
      return core::detail::publish_error(
        base::StatusCode::unsupported,
        "Auto CFD Surface Mesher currently supports trimmed faces only through classified outer and inner loops."
      );
    }
    inner_loops.push_back(seed_loops[loop_index]);
  }

  enforce_outer_inner_loop_winding(outer_loop, inner_loops);

  std::vector<std::array<std::uint32_t, 3>> triangles;
  if(!build_boundary_cdt_with_holes(plan.vertices, outer_loop, inner_loops, triangles)) {
    return core::detail::publish_error(
      base::StatusCode::unsupported,
      "Auto CFD Surface Mesher planar CDT seeding failed."
    );
  }

  plan.triangles.reserve(triangles.size());
  for(auto triangle : triangles) {
    const double area = signed_triangle_area_twice(
      plan.vertices[triangle[0]].uv,
      plan.vertices[triangle[1]].uv,
      plan.vertices[triangle[2]].uv
    );
    if(area < 0.0) {
      std::swap(triangle[1], triangle[2]);
    }
    plan.triangles.push_back({triangle, true, false});
  }

  if(plan.triangles.empty()) {
    return core::detail::publish_error(
      base::StatusCode::unsupported,
      "Auto CFD Surface Mesher seeded only degenerate planar triangles."
    );
  }

  return core::detail::clear_error_state();
}

struct PatchSeedCenterCandidate final {
  std::array<double, 2> uv {0.0, 0.0};
  geo::Point3 position {0.0, 0.0, 0.0};
};

void append_single_loop_patch_seed_center_candidate(
  std::vector<PatchSeedCenterCandidate> &candidates,
  const std::array<double, 2> &uv,
  const geo::Point3 &position
) noexcept
{
  const auto duplicate_it = std::find_if(
    candidates.begin(),
    candidates.end(),
    [&uv, &position](const PatchSeedCenterCandidate &candidate) {
      const auto uv_delta = uv_subtract(candidate.uv, uv);
      const auto point_delta = point_subtract(candidate.position, position);
      return uv_dot(uv_delta, uv_delta) <= kUvAreaTolerance &&
             dot_product(point_delta, point_delta) <= kProjectionToleranceFloor;
    }
  );
  if(duplicate_it == candidates.end()) {
    candidates.push_back({uv, position});
  }
}

void collect_single_loop_patch_seed_centers(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  std::vector<PatchSeedCenterCandidate> &candidates
)
{
  candidates.clear();
  std::array<double, 2> averaged_uv {0.0, 0.0};
  std::size_t averaged_uv_count = 0U;
  geo::Point3 averaged_point {0.0, 0.0, 0.0};
  std::size_t averaged_point_count = 0U;
  if(!face_preprocess.loops.empty()) {
    for(const auto &point : face_preprocess.loops.front().points) {
      if(point.uv_defined) {
        averaged_uv[0] += point.uv[0];
        averaged_uv[1] += point.uv[1];
        ++averaged_uv_count;
      }
      averaged_point[0] += point.position[0];
      averaged_point[1] += point.position[1];
      averaged_point[2] += point.position[2];
      ++averaged_point_count;
    }
  }
  if(averaged_uv_count > 0U) {
    const double scale = 1.0 / static_cast<double>(averaged_uv_count);
    averaged_uv[0] *= scale;
    averaged_uv[1] *= scale;
  }

  const auto try_material_uv =
    [&](const std::array<double, 2> &candidate_uv) {
      if(!std::isfinite(candidate_uv[0]) ||
         !std::isfinite(candidate_uv[1]) ||
         !uv_point_inside_face_material(face_preprocess, candidate_uv)) {
        return;
      }

      geo::FaceSample sample;
      if(sample_face_at_material_uv(face_view, face_preprocess, candidate_uv, sample)) {
        append_single_loop_patch_seed_center_candidate(
          candidates,
          candidate_uv,
          sample.position
        );
      }
    };

  const std::array<std::array<double, 2>, 3> material_uv_candidates {{
    {
      face_preprocess.reference_metric.u,
      face_preprocess.reference_metric.v,
    },
    averaged_uv,
    {
      0.5 * (face_preprocess.recovered_uv_min[0] + face_preprocess.recovered_uv_max[0]),
      0.5 * (face_preprocess.recovered_uv_min[1] + face_preprocess.recovered_uv_max[1]),
    },
  }};
  for(const auto candidate_uv : material_uv_candidates) {
    try_material_uv(candidate_uv);
  }

  if(averaged_point_count > 0U) {
    const double scale = 1.0 / static_cast<double>(averaged_point_count);
    averaged_point[0] *= scale;
    averaged_point[1] *= scale;
    averaged_point[2] *= scale;

    geo::FaceProjection projection;
    if(geo::project_point_to_face(face_view, averaged_point, projection) ==
         base::StatusCode::ok) {
      {
        const std::array<double, 2> candidate_uv {projection.u, projection.v};
        if(uv_point_inside_face_material(face_preprocess, candidate_uv)) {
          append_single_loop_patch_seed_center_candidate(
            candidates,
            candidate_uv,
            projection.projected_point
          );
        }
      }
    }
  }
}

// ── Constrained Delaunay Triangulation (CDT) for seed mesh ──────────────
// Replaces the fan-from-center seed with a proper CDT of boundary points
// in UV space. This produces much higher quality initial triangles, greatly
// reducing the work required by the frontal Delaunay refinement step.
//
// Algorithm:
//   1. Create a super-triangle enclosing all boundary UV points.
//   2. Incrementally insert each boundary point using Bowyer-Watson.
//   3. Remove super-triangle vertices and incident triangles.
//   4. Recover constrained boundary edges by edge-flipping.
//   5. Remove triangles outside the boundary polygon.

namespace cdt_detail {

namespace robust_pred {

// splitter = 2^ceil(53/2) + 1 for IEEE 754 double.
// epsilon = 2^(-53).
// These are computed at first use via init().
static double s_splitter = 0.0;
static double s_epsilon = 0.0;
static double s_ccwerrboundA = 0.0;
static double s_ccwerrboundB = 0.0;
static double s_ccwerrboundC = 0.0;
static double s_resulterrbound = 0.0;
static bool s_initialized = false;

inline void init() noexcept
{
  if(s_initialized) return;
  double half = 0.5;
  double eps = 1.0;
  double spl = 1.0;
  double check = 1.0;
  double lastcheck;
  int every_other = 1;
  do {
    lastcheck = check;
    eps *= half;
    if(every_other) spl *= 2.0;
    every_other = !every_other;
    check = 1.0 + eps;
  } while(check != 1.0 && check != lastcheck);
  spl += 1.0;
  s_splitter = spl;
  s_epsilon = eps;
  s_resulterrbound = (3.0 + 8.0 * eps) * eps;
  s_ccwerrboundA = (3.0 + 16.0 * eps) * eps;
  s_ccwerrboundB = (2.0 + 12.0 * eps) * eps;
  s_ccwerrboundC = (9.0 + 64.0 * eps) * eps * eps;
  s_initialized = true;
}

// Error-free arithmetic macros (inlined as lambdas for safety).
// These match Shewchuk's Two_Sum, Two_Diff, Split, Two_Product, etc.

#define RP_Fast_Two_Sum(a, b, x, y) \
  do { (x) = (a) + (b); double rp_bv = (x) - (a); (y) = (b) - rp_bv; } while(0)

#define RP_Two_Sum(a, b, x, y) \
  do { \
    (x) = (a) + (b); \
    double rp_bv = (x) - (a); \
    double rp_av = (x) - rp_bv; \
    double rp_br = (b) - rp_bv; \
    double rp_ar = (a) - rp_av; \
    (y) = rp_ar + rp_br; \
  } while(0)

#define RP_Two_Diff(a, b, x, y) \
  do { \
    (x) = (a) - (b); \
    double rp_bv = (a) - (x); \
    double rp_av = (x) + rp_bv; \
    double rp_br = rp_bv - (b); \
    double rp_ar = (a) - rp_av; \
    (y) = rp_ar + rp_br; \
  } while(0)

#define RP_Split(a, ahi, alo) \
  do { \
    double rp_c = s_splitter * (a); \
    double rp_ab = rp_c - (a); \
    (ahi) = rp_c - rp_ab; \
    (alo) = (a) - (ahi); \
  } while(0)

#define RP_Two_Product(a, b, x, y) \
  do { \
    (x) = (a) * (b); \
    double rp_ahi, rp_alo, rp_bhi, rp_blo; \
    RP_Split(a, rp_ahi, rp_alo); \
    RP_Split(b, rp_bhi, rp_blo); \
    double rp_e1 = (x) - rp_ahi * rp_bhi; \
    double rp_e2 = rp_e1 - rp_alo * rp_bhi; \
    double rp_e3 = rp_e2 - rp_ahi * rp_blo; \
    (y) = rp_alo * rp_blo - rp_e3; \
  } while(0)

inline int fast_expansion_sum_zeroelim(
  int elen, const double *e, int flen, const double *f, double *h
) noexcept
{
  if(elen == 0) { for(int i = 0; i < flen; ++i) h[i] = f[i]; return flen; }
  if(flen == 0) { for(int i = 0; i < elen; ++i) h[i] = e[i]; return elen; }

  double Q, Qnew, hh;
  int eindex = 0, findex = 0, hindex = 0;
  double enow = e[0], fnow = f[0];
  if((fnow > enow) == (fnow > -enow)) { Q = enow; ++eindex; }
  else { Q = fnow; ++findex; }

  if(eindex < elen && findex < flen) {
    enow = e[eindex]; fnow = f[findex];
    if((fnow > enow) == (fnow > -enow)) { RP_Fast_Two_Sum(enow, Q, Qnew, hh); ++eindex; }
    else { RP_Fast_Two_Sum(fnow, Q, Qnew, hh); ++findex; }
    Q = Qnew;
    if(hh != 0.0) h[hindex++] = hh;
    while(eindex < elen && findex < flen) {
      enow = e[eindex]; fnow = f[findex];
      if((fnow > enow) == (fnow > -enow)) { RP_Two_Sum(Q, enow, Qnew, hh); ++eindex; }
      else { RP_Two_Sum(Q, fnow, Qnew, hh); ++findex; }
      Q = Qnew;
      if(hh != 0.0) h[hindex++] = hh;
    }
  }
  while(eindex < elen) { enow = e[eindex]; RP_Two_Sum(Q, enow, Qnew, hh); ++eindex; Q = Qnew; if(hh != 0.0) h[hindex++] = hh; }
  while(findex < flen) { fnow = f[findex]; RP_Two_Sum(Q, fnow, Qnew, hh); ++findex; Q = Qnew; if(hh != 0.0) h[hindex++] = hh; }
  if(Q != 0.0 || hindex == 0) h[hindex++] = Q;
  return hindex;
}

inline double estimate(int elen, const double *e) noexcept
{
  double Q = e[0];
  for(int i = 1; i < elen; ++i) Q += e[i];
  return Q;
}

double orient2dadapt(const double *pa, const double *pb, const double *pc, double detsum) noexcept
{
  double acx = pa[0] - pc[0], bcx = pb[0] - pc[0];
  double acy = pa[1] - pc[1], bcy = pb[1] - pc[1];
  double detleft, detlefttail, detright, detrighttail;
  RP_Two_Product(acx, bcy, detleft, detlefttail);
  RP_Two_Product(acy, bcx, detright, detrighttail);

  double B[4];
  // Two_Two_Diff(detleft, detlefttail, detright, detrighttail, B3, B[2], B[1], B[0])
  double _i, _j, _0, B3;
  RP_Two_Diff(detlefttail, detrighttail, _i, B[0]);
  RP_Two_Sum(detleft, _i, _j, _0);
  // Wait, Two_Two_Diff is: Two_One_Diff(a1,a0,b0,...) then Two_One_Diff(...)
  // Let me inline it properly:
  // Two_One_Diff(detleft, detlefttail, detright, _j, _0, B[0]):
  //   Two_Diff(detlefttail, detright, _i, B[0]);
  //   Two_Sum(detleft, _i, _j, _0 /*=B[1]*/);
  // Two_One_Diff(_j, _0, detrighttail, B3, B[2], B[1]):
  //   Two_Diff(_0, detrighttail, _i2, B[1]);
  //   Two_Sum(_j, _i2, B3, B[2]);
  // Redo:
  RP_Two_Diff(detlefttail, detrighttail, _i, B[0]);
  RP_Two_Sum(detleft, _i, _j, _0);
  // Now second part:
  double _i2;
  RP_Two_Diff(_0, 0.0, _i2, B[1]); // no, this is wrong
  
  (void)_i; (void)_j; (void)_0; (void)B3; (void)_i2;
  double t_i, t_j, t_0;
  RP_Two_Diff(detlefttail, detrighttail, t_i, B[0]);
  RP_Two_Sum(detleft, t_i, t_j, t_0);
  double t_i2;
  RP_Two_Diff(t_0, detright, t_i2, B[1]);
  RP_Two_Sum(t_j, t_i2, B3, B[2]);
  B[3] = B3;

  double det = estimate(4, B);
  double errbound = s_ccwerrboundB * detsum;
  if(det >= errbound || -det >= errbound) return det;

  double acxtail, acytail, bcxtail, bcytail;
  RP_Two_Diff(pa[0], pc[0], acx, acxtail);
  RP_Two_Diff(pb[0], pc[0], bcx, bcxtail);
  RP_Two_Diff(pa[1], pc[1], acy, acytail);
  RP_Two_Diff(pb[1], pc[1], bcy, bcytail);

  if(acxtail == 0.0 && acytail == 0.0 && bcxtail == 0.0 && bcytail == 0.0) {
    return det;
  }

  errbound = s_ccwerrboundC * detsum + s_resulterrbound * std::abs(det);
  det += (acx * bcytail + bcy * acxtail) - (acy * bcxtail + bcx * acytail);
  if(det >= errbound || -det >= errbound) return det;

  double s1, s0, t1, t0;
  double u[4];
  RP_Two_Product(acxtail, bcy, s1, s0);
  RP_Two_Product(acytail, bcx, t1, t0);
  // Two_Two_Diff(s1, s0, t1, t0, u[3], u[2], u[1], u[0])
  double u3;
  {
    double ti, tj, t0_;
    RP_Two_Diff(s0, t0, ti, u[0]);
    RP_Two_Sum(s1, ti, tj, t0_);
    double ti2;
    RP_Two_Diff(t0_, t1, ti2, u[1]);
    RP_Two_Sum(tj, ti2, u3, u[2]);
  }
  u[3] = u3;
  double C1[8];
  int C1length = fast_expansion_sum_zeroelim(4, B, 4, u, C1);

  RP_Two_Product(acx, bcytail, s1, s0);
  RP_Two_Product(acy, bcxtail, t1, t0);
  {
    double ti, tj, t0_;
    RP_Two_Diff(s0, t0, ti, u[0]);
    RP_Two_Sum(s1, ti, tj, t0_);
    double ti2;
    RP_Two_Diff(t0_, t1, ti2, u[1]);
    RP_Two_Sum(tj, ti2, u3, u[2]);
  }
  u[3] = u3;
  double C2[12];
  int C2length = fast_expansion_sum_zeroelim(C1length, C1, 4, u, C2);

  RP_Two_Product(acxtail, bcytail, s1, s0);
  RP_Two_Product(acytail, bcxtail, t1, t0);
  {
    double ti, tj, t0_;
    RP_Two_Diff(s0, t0, ti, u[0]);
    RP_Two_Sum(s1, ti, tj, t0_);
    double ti2;
    RP_Two_Diff(t0_, t1, ti2, u[1]);
    RP_Two_Sum(tj, ti2, u3, u[2]);
  }
  u[3] = u3;
  double D[16];
  int Dlength = fast_expansion_sum_zeroelim(C2length, C2, 4, u, D);

  return D[Dlength - 1];
}

// Adaptive exact 2D orientation test (Shewchuk, public domain).
// Returns positive if pa, pb, pc are in CCW order; negative if CW; zero if collinear.
double orient2d(const double *pa, const double *pb, const double *pc) noexcept
{
  init();
  double detleft = (pa[0] - pc[0]) * (pb[1] - pc[1]);
  double detright = (pa[1] - pc[1]) * (pb[0] - pc[0]);
  double det = detleft - detright;

  if(detleft > 0.0) {
    if(detright <= 0.0) return det;
  } else if(detleft < 0.0) {
    if(detright >= 0.0) return det;
  } else {
    return det;
  }

  double detsum = (detleft > 0.0) ? (detleft + detright) : (-detleft - detright);
  double errbound = s_ccwerrboundA * detsum;
  if(det >= errbound || -det >= errbound) return det;

  return orient2dadapt(pa, pb, pc, detsum);
}

#undef RP_Fast_Two_Sum
#undef RP_Two_Sum
#undef RP_Two_Diff
#undef RP_Split
#undef RP_Two_Product

} // namespace robust_pred

struct CdtTriangle {
  std::array<std::uint32_t, 3> vertices {invalid_index, invalid_index, invalid_index};
  std::array<std::uint32_t, 3> neighbors {invalid_index, invalid_index, invalid_index};
  bool valid = true;
};

[[nodiscard]] std::uint32_t find_neighbor_edge_slot(
  const CdtTriangle &tri,
  std::uint32_t neighbor_index
) noexcept
{
  for(std::uint32_t i = 0; i < 3; ++i) {
    if(tri.neighbors[i] == neighbor_index) return i;
  }
  return invalid_index;
}

// Find which triangle contains the point p (linear scan, fine for seed meshes).
[[nodiscard]] std::uint32_t locate_triangle(
  const std::vector<std::array<double, 2>> &uvs,
  const std::vector<CdtTriangle> &triangles,
  const std::array<double, 2> &p
) noexcept
{
  for(std::uint32_t i = 0; i < static_cast<std::uint32_t>(triangles.size()); ++i) {
    if(!triangles[i].valid) continue;
    const auto &a = uvs[triangles[i].vertices[0]];
    const auto &b = uvs[triangles[i].vertices[1]];
    const auto &c = uvs[triangles[i].vertices[2]];
    const double d0 = (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0]);
    const double d1 = (c[0] - b[0]) * (p[1] - b[1]) - (c[1] - b[1]) * (p[0] - b[0]);
    const double d2 = (a[0] - c[0]) * (p[1] - c[1]) - (a[1] - c[1]) * (p[0] - c[0]);
    if(d0 >= -1.0e-15 && d1 >= -1.0e-15 && d2 >= -1.0e-15) return i;
  }
  return invalid_index;
}

// Rebuild adjacency for all valid triangles (simple O(n²) but n is small for seed).
void rebuild_adjacency(std::vector<CdtTriangle> &triangles) noexcept
{
  for(auto &tri : triangles) {
    tri.neighbors = {invalid_index, invalid_index, invalid_index};
  }
  const auto n = static_cast<std::uint32_t>(triangles.size());
  for(std::uint32_t i = 0; i < n; ++i) {
    if(!triangles[i].valid) continue;
    for(std::uint32_t j = i + 1; j < n; ++j) {
      if(!triangles[j].valid) continue;
      for(std::uint32_t ei = 0; ei < 3; ++ei) {
        const auto ia = triangles[i].vertices[ei];
        const auto ib = triangles[i].vertices[(ei + 1) % 3];
        for(std::uint32_t ej = 0; ej < 3; ++ej) {
          const auto ja = triangles[j].vertices[ej];
          const auto jb = triangles[j].vertices[(ej + 1) % 3];
          if((ia == jb && ib == ja) || (ia == ja && ib == jb)) {
            triangles[i].neighbors[ei] = j;
            triangles[j].neighbors[ej] = i;
          }
        }
      }
    }
  }
}

// Isotropic in-circumcircle test (positive if p is inside circumcircle of abc).
[[nodiscard]] bool in_circumcircle(
  const std::vector<std::array<double, 2>> &uvs,
  std::uint32_t ai, std::uint32_t bi, std::uint32_t ci,
  const std::array<double, 2> &p
) noexcept
{
  const auto &a = uvs[ai];
  const auto &b = uvs[bi];
  const auto &c = uvs[ci];
  const double ax = a[0] - p[0], ay = a[1] - p[1];
  const double bx = b[0] - p[0], by = b[1] - p[1];
  const double cx = c[0] - p[0], cy = c[1] - p[1];
  const double det =
    (ax*ax + ay*ay) * (bx*cy - by*cx) -
    (bx*bx + by*by) * (ax*cy - ay*cx) +
    (cx*cx + cy*cy) * (ax*by - ay*bx);
  return det > 1.0e-15;
}

// Bowyer-Watson insertion of a single point into the triangulation.
void insert_point(
  const std::vector<std::array<double, 2>> &uvs,
  std::vector<CdtTriangle> &triangles,
  std::uint32_t point_index
)
{
  const auto &p = uvs[point_index];
  const auto containing = locate_triangle(uvs, triangles, p);
  if(containing == invalid_index) return;

  // Find all triangles whose circumcircle contains the new point (cavity).
  std::vector<std::uint32_t> cavity;
  std::vector<bool> in_cavity(triangles.size(), false);
  std::vector<std::uint32_t> stack;
  stack.push_back(containing);
  in_cavity[containing] = true;
  while(!stack.empty()) {
    const auto ti = stack.back();
    stack.pop_back();
    cavity.push_back(ti);
    for(std::uint32_t e = 0; e < 3; ++e) {
      const auto ni = triangles[ti].neighbors[e];
      if(ni == invalid_index || in_cavity[ni] || !triangles[ni].valid) continue;
      const auto &nt = triangles[ni];
      if(in_circumcircle(uvs, nt.vertices[0], nt.vertices[1], nt.vertices[2], p)) {
        in_cavity[ni] = true;
        stack.push_back(ni);
      }
    }
  }

  // Extract the boundary shell of the cavity.
  struct ShellEdge { std::uint32_t v0, v1; std::uint32_t outside_tri; };
  std::vector<ShellEdge> shell;
  for(auto ti : cavity) {
    for(std::uint32_t e = 0; e < 3; ++e) {
      const auto ni = triangles[ti].neighbors[e];
      if(ni == invalid_index || !in_cavity[ni]) {
        shell.push_back({
          triangles[ti].vertices[e],
          triangles[ti].vertices[(e + 1) % 3],
          ni
        });
      }
    }
  }

  // Invalidate cavity triangles.
  for(auto ti : cavity) {
    triangles[ti].valid = false;
  }

  // Create new triangles from shell edges to the new point.
  std::vector<std::uint32_t> new_tri_indices;
  new_tri_indices.reserve(shell.size());
  for(const auto &edge : shell) {
    CdtTriangle nt;
    nt.vertices = {edge.v0, edge.v1, point_index};
    // Check orientation
    const auto &a = uvs[edge.v0];
    const auto &b = uvs[edge.v1];
    const double area = (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0]);
    if(area < 0.0) {
      std::swap(nt.vertices[0], nt.vertices[1]);
    }
    const auto idx = static_cast<std::uint32_t>(triangles.size());
    triangles.push_back(nt);
    new_tri_indices.push_back(idx);
  }

  // Rebuild adjacency (simple approach for seed meshes).
  rebuild_adjacency(triangles);
}

// Check if edge (a,b) exists in the triangulation.
[[nodiscard]] bool edge_exists(
  const std::vector<CdtTriangle> &triangles,
  std::uint32_t a, std::uint32_t b
) noexcept
{
  for(const auto &tri : triangles) {
    if(!tri.valid) continue;
    for(std::uint32_t e = 0; e < 3; ++e) {
      if((tri.vertices[e] == a && tri.vertices[(e+1)%3] == b) ||
         (tri.vertices[e] == b && tri.vertices[(e+1)%3] == a)) {
        return true;
      }
    }
  }
  return false;
}


[[nodiscard]] bool segments_cross_strict(
  const std::array<double, 2> &p1, const std::array<double, 2> &p2,
  const std::array<double, 2> &q1, const std::array<double, 2> &q2
) noexcept
{
  // Solve: p1 + t*(p2-p1) = q1 + s*(q2-q1)
  // => [p2x-p1x, -(q2x-q1x)] [t]   [q1x-p1x]
  //    [p2y-p1y, -(q2y-q1y)] [s] = [q1y-p1y]
  const double a00 = p2[0] - p1[0];
  const double a01 = -(q2[0] - q1[0]);
  const double a10 = p2[1] - p1[1];
  const double a11 = -(q2[1] - q1[1]);
  const double r0 = q1[0] - p1[0];
  const double r1 = q1[1] - p1[1];

  const double det = a00 * a11 - a01 * a10;
  if(std::abs(det) < 1.0e-30) return false;  // parallel / collinear

  const double inv_det = 1.0 / det;
  const double t = (a11 * r0 - a01 * r1) * inv_det;
  const double s = (-a10 * r0 + a00 * r1) * inv_det;

  
  constexpr double eps = 1.0e-12;
  return (t > eps && t < 1.0 - eps && s > eps && s < 1.0 - eps);
}

bool try_cdt_swap(
  const std::vector<std::array<double, 2>> &uvs,
  std::vector<CdtTriangle> &triangles,
  std::uint32_t ti,
  std::uint32_t edge_slot
)
{
  auto &tri = triangles[ti];
  if(!tri.valid) return false;
  const auto ni = tri.neighbors[edge_slot];
  if(ni == invalid_index) return false;
  auto &neigh = triangles[ni];
  if(!neigh.valid) return false;

  const auto e0 = tri.vertices[edge_slot];
  const auto e1 = tri.vertices[(edge_slot + 1) % 3];
  const auto e2 = tri.vertices[(edge_slot + 2) % 3];  // opposite in tri

  // Find opposite vertex in neighbor.
  std::uint32_t opp_slot = invalid_index;
  for(std::uint32_t ne = 0; ne < 3; ++ne) {
    if(neigh.vertices[ne] != e0 && neigh.vertices[ne] != e1) {
      opp_slot = ne;
      break;
    }
  }
  if(opp_slot == invalid_index) return false;
  const auto opp_v = neigh.vertices[opp_slot];
  if(opp_v == e2) return false;  // degenerate

  const auto area = [&](std::uint32_t a, std::uint32_t b, std::uint32_t c) {
    return (uvs[b][0]-uvs[a][0])*(uvs[c][1]-uvs[a][1])
         - (uvs[b][1]-uvs[a][1])*(uvs[c][0]-uvs[a][0]);
  };

  // Convexity check: orient2d(e2, e0, opp_v) and orient2d(e2, opp_v, e1)
  // must have the same sign — the quadrilateral (e0, e2, e1, opp_v) is
  // strictly convex around the new diagonal (e2, opp_v).
  const double ori_t1 = area(opp_v, e0, e2);
  const double ori_t2 = area(opp_v, e2, e1);
  if(ori_t1 * ori_t2 <= 0.0) return false;

  // Determine orientation of each new triangle for CCW output.
  const double a1_raw = area(e2, opp_v, e0);
  const double a2_raw = area(opp_v, e2, e1);

  // Perform the swap in-place, ensuring CCW orientation.
  tri.vertices = {e2, opp_v, e0};
  if(a1_raw < 0.0) std::swap(tri.vertices[1], tri.vertices[2]);
  neigh.vertices = {opp_v, e2, e1};
  if(a2_raw < 0.0) std::swap(neigh.vertices[1], neigh.vertices[2]);

  
  rebuild_adjacency(triangles);
  return true;
}

void divide_and_conquer_triangulate(
  const std::vector<std::array<double, 2>> &points,
  std::vector<std::array<std::uint32_t, 3>> &out_triangles
)
{
  // Initialise Shewchuk's adaptive predicates on first entry.
  // Thread-safe per C++11 magic-statics; cheap on subsequent calls.
  static const bool s_robust_predicates_ready = []{
    ::exactinit(0, 0, 0, 1.0, 1.0, 1.0);
    return true;
  }();
  (void)s_robust_predicates_ready;

  using PointNum = std::int32_t;

  const auto n = static_cast<PointNum>(points.size());
  if(n < 3) return;

  
  struct AdjNode {
    PointNum point_num = -1;
    AdjNode *next = nullptr;
    AdjNode *prev = nullptr;
  };

  
  std::vector<PointNum> sorted_to_orig(static_cast<std::size_t>(n));
  for(PointNum i = 0; i < n; ++i) sorted_to_orig[static_cast<std::size_t>(i)] = i;
  std::sort(sorted_to_orig.begin(), sorted_to_orig.end(),
    [&](PointNum a, PointNum b) {
      const double dx = points[static_cast<std::size_t>(a)][0]
                      - points[static_cast<std::size_t>(b)][0];
      if(dx != 0.0) return dx < 0.0;
      return points[static_cast<std::size_t>(a)][1]
           < points[static_cast<std::size_t>(b)][1];
    }
  );
  // orig_to_sorted[original_index] = sorted position.
  std::vector<PointNum> orig_to_sorted(static_cast<std::size_t>(n));
  for(PointNum i = 0; i < n; ++i) {
    orig_to_sorted[static_cast<std::size_t>(sorted_to_orig[static_cast<std::size_t>(i)])] = i;
  }

  // Accessor: get UV coords for sorted index s.
  const auto px = [&](PointNum s) -> double {
    return points[static_cast<std::size_t>(sorted_to_orig[static_cast<std::size_t>(s)])][0];
  };
  const auto py = [&](PointNum s) -> double {
    return points[static_cast<std::size_t>(sorted_to_orig[static_cast<std::size_t>(s)])][1];
  };

  // Per-point adjacency list heads (sorted index -> list).
  std::vector<AdjNode *> adj(static_cast<std::size_t>(n), nullptr);

  // Track all allocated AdjNodes for cleanup.
  std::vector<AdjNode *> all_nodes;
  all_nodes.reserve(static_cast<std::size_t>(n * 6));

  
  const auto is_left_of = [&](PointNum x, PointNum y, PointNum check) -> bool {
    const double pa[2] = {px(x), py(x)};
    const double pb[2] = {px(y), py(y)};
    const double pc[2] = {px(check), py(check)};
    return robust_pred::orient2d(pa, pb, pc) > 0.0;
  };

  const auto is_right_of = [&](PointNum x, PointNum y, PointNum check) -> bool {
    return is_left_of(y, x, check);
  };

  
  const auto qtest = [&](PointNum h, PointNum i, PointNum j, PointNum k) -> bool {
    // Non-const arrays: ::incircle() takes double* (Shewchuk's PD signature).
    double ph[2] = {px(h), py(h)};
    double pi[2] = {px(i), py(i)};
    double pj[2] = {px(j), py(j)};
    double pk[2] = {px(k), py(k)};

    const double orient = robust_pred::orient2d(ph, pi, pj);
    const double incircle = ::incircle(ph, pi, pj, pk);

    return (incircle * orient < 0.0);
  };

  
  const auto predecessor = [&](PointNum a, PointNum b) -> PointNum {
    AdjNode *p = adj[static_cast<std::size_t>(a)];
    if(!p) return -1;
    AdjNode *start = p;
    do {
      if(p->point_num == b) return p->prev->point_num;
      p = p->prev;
    } while(p != start);
    return -1;
  };

  
  const auto successor = [&](PointNum a, PointNum b) -> PointNum {
    AdjNode *p = adj[static_cast<std::size_t>(a)];
    if(!p) return -1;
    AdjNode *start = p;
    do {
      if(p->point_num == b) return p->next->point_num;
      p = p->next;
    } while(p != start);
    return -1;
  };

  // First: first neighbor of x.
  const auto first = [&](PointNum x) -> PointNum {
    return adj[static_cast<std::size_t>(x)]->point_num;
  };

  
  const auto fix_first = [&](PointNum x, PointNum f) -> bool {
    AdjNode *p = adj[static_cast<std::size_t>(x)];
    if(!p) return false;
    AdjNode *copy = p;
    do {
      if(p->point_num == f) {
        adj[static_cast<std::size_t>(x)] = p;
        return true;
      }
      p = p->next;
    } while(p != copy);
    return false;
  };

  
  const auto dlist_insert = [&](PointNum centerPoint, PointNum newPoint) -> bool {
    auto *newp = new AdjNode;
    newp->point_num = newPoint;
    all_nodes.push_back(newp);

    AdjNode *&dlist = adj[static_cast<std::size_t>(centerPoint)];
    if(dlist == nullptr) {
      dlist = newp;
      newp->prev = newp;
      newp->next = newp;
      return true;
    }
    if(dlist->next == dlist) {
      dlist->prev = newp;
      dlist->next = newp;
      newp->prev = dlist;
      newp->next = dlist;
      return true;
    }

    // 2 or more elements: insert by polar angle.
    AdjNode *p = dlist;
    const PointNum first_pn = p->point_num;

    const double center_x = px(centerPoint);
    const double center_y = py(centerPoint);

    const double yy0 = py(first_pn) - center_y;
    const double xx0 = px(first_pn) - center_x;
    const double alpha1 = std::atan2(yy0, xx0);

    const double yyn = py(newPoint) - center_y;
    const double xxn = px(newPoint) - center_x;
    double beta = std::atan2(yyn, xxn) - alpha1;
    if(beta <= 0.0) beta += 2.0 * kPi;

    do {
      const double yys = py(p->next->point_num) - center_y;
      const double xxs = px(p->next->point_num) - center_x;
      double alpha = std::atan2(yys, xxs) - alpha1;
      if(alpha <= -1.0e-15 || p->next->point_num == first_pn)
        alpha += 2.0 * kPi;
      else if(std::abs(alpha) <= 1e-15 &&
              is_right_of(centerPoint, first_pn, p->next->point_num))
        alpha += 2.0 * kPi;

      if(alpha >= beta + 1e-15) {
        newp->next = p->next;
        p->next = newp;
        newp->prev = p;
        newp->next->prev = newp;
        return true;
      }
      else if(alpha >= beta - 1e-15) {
        // Nearly same angle — use orientation to break tie.
        if(is_right_of(centerPoint, p->next->point_num, newPoint)) {
          newp->next = p->next;
          p->next = newp;
          newp->prev = p;
          newp->next->prev = newp;
        }
        else {
          newp->next = p->next->next;
          p->next->next = newp;
          newp->prev = p->next;
          newp->next->prev = newp;
        }
        return true;
      }
      p = p->next;
    } while(p != dlist);

    return false;
  };

 
  const auto insert = [&](PointNum a, PointNum b) -> bool {
    bool r = dlist_insert(a, b);
    r &= dlist_insert(b, a);
    return r;
  };

  
  const auto dlist_delete = [&](AdjNode *&dlist, PointNum oldPoint) -> bool {
    if(dlist == nullptr) return false;
    if(dlist->next == dlist) {
      if(dlist->point_num == oldPoint) {
        // Don't delete — tracked in all_nodes for bulk cleanup.
        dlist = nullptr;
        return true;
      }
      return false;
    }
    AdjNode *p = dlist;
    do {
      if(p->point_num == oldPoint) {
        p->next->prev = p->prev;
        p->prev->next = p->next;
        if(p == dlist) dlist = p->next;
        // Don't delete — tracked in all_nodes for bulk cleanup.
        return true;
      }
      p = p->next;
    } while(p != dlist);
    return false;
  };

 
  const auto delete_edge = [&](PointNum a, PointNum b) -> bool {
    bool r = dlist_delete(adj[static_cast<std::size_t>(a)], b);
    r &= dlist_delete(adj[static_cast<std::size_t>(b)], a);
    return r;
  };

  // ── DT range struct (begin, end of sorted indices) ───────────────
  struct DT { PointNum begin; PointNum end; };

 
  const auto lower_common_tangent = [&](DT vl, DT vr) -> std::pair<PointNum, PointNum> {
    PointNum x = vl.end;
    PointNum y = vr.begin;
    PointNum z = first(y);
    PointNum z1 = first(x);
    PointNum z2 = predecessor(x, z1);
    for(;;) {
      if(is_right_of(x, y, z)) {
        PointNum temp = z;
        z = successor(z, y);
        y = temp;
      }
      else if(is_right_of(x, y, z2)) {
        PointNum temp = z2;
        z2 = predecessor(z2, x);
        x = temp;
      }
      else {
        return {x, y};
      }
    }
  };

  
  const auto upper_common_tangent = [&](DT vl, DT vr) -> std::pair<PointNum, PointNum> {
    PointNum x = vl.end;
    PointNum y = vr.begin;
    PointNum z = first(y);
    PointNum z1 = first(x);
    PointNum z2 = predecessor(y, z);
    for(;;) {
      if(is_left_of(x, y, z2)) {
        PointNum temp = z2;
        z2 = predecessor(z2, y);
        y = temp;
      }
      else if(is_left_of(x, y, z1)) {
        PointNum temp = z1;
        z1 = successor(z1, x);
        x = temp;
      }
      else {
        return {x, y};
      }
    }
  };

  
  const auto merge = [&](DT vl, DT vr) -> bool {
    auto [bt_l, bt_r] = lower_common_tangent(vl, vr);
    auto [ut_l, ut_r] = upper_common_tangent(vl, vr);
    PointNum l = bt_l;
    PointNum r = bt_r;

    while(l != ut_l || r != ut_r) {
      int a = 0, b = 0;
      if(!insert(l, r)) return false;

      PointNum r1 = predecessor(r, l);
      if(r1 == -1) return false;
      if(is_right_of(l, r, r1)) {
        a = 1;
      }
      else {
        bool out = false;
        while(!out) {
          PointNum r2 = predecessor(r, r1);
          if(r2 == -1) return false;
          if(r2 < vr.begin)
            out = true;
          else if(qtest(l, r, r1, r2))
            out = true;
          else {
            if(!delete_edge(r, r1)) return false;
            r1 = r2;
            if(is_right_of(l, r, r1)) { out = true; a = 1; }
          }
        }
      }

      PointNum l1 = successor(l, r);
      if(l1 == -1) return false;
      if(is_left_of(r, l, l1)) {
        b = 1;
      }
      else {
        bool out = false;
        while(!out) {
          PointNum l2 = successor(l, l1);
          if(l2 == -1) return false;
          if(l2 > vl.end)
            out = true;
          else if(qtest(r, l, l1, l2))
            out = true;
          else {
            if(!delete_edge(l, l1)) return false;
            l1 = l2;
            if(is_left_of(r, l, l1)) { out = true; b = 1; }
          }
        }
      }

      if(a)
        l = l1;
      else if(b)
        r = r1;
      else {
        if(qtest(l, r, r1, l1))
          r = r1;
        else
          l = l1;
      }
    }
    if(!insert(l, r)) return false;

    fix_first(ut_r, ut_l);
    fix_first(bt_l, bt_r);
    return true;
  };

  
  std::function<DT(PointNum, PointNum)> recur_trig =
    [&](PointNum left, PointNum right) -> DT
  {
    DT dt {left, right};
    const PointNum count = right - left + 1;
    switch(count) {
    case 0:
    case 1:
      break;
    case 2:
      insert(left, right);
      fix_first(left, right);
      fix_first(right, left);
      break;
    case 3:
      insert(left, right);
      insert(left, left + 1);
      insert(left + 1, right);
      if(is_right_of(left, right, left + 1)) {
        fix_first(left, left + 1);
        fix_first(left + 1, right);
        fix_first(right, left);
      }
      else {
        fix_first(left, right);
        fix_first(left + 1, left);
        fix_first(right, left + 1);
      }
      break;
    default: {
      PointNum m = (left + right) >> 1;
      DT left_dt = recur_trig(left, m);
      DT right_dt = recur_trig(m + 1, right);
      merge(left_dt, right_dt);
      break;
    }
    }
    return dt;
  };

  // ── Run the divide-and-conquer algorithm ─────────────────────────
  recur_trig(0, n - 1);

  
  out_triangles.clear();
  for(PointNum i = 0; i < n; ++i) {
    AdjNode *p = adj[static_cast<std::size_t>(i)];
    if(!p) continue;
    std::vector<PointNum> neighbors;
    AdjNode *start = p;
    do {
      neighbors.push_back(p->point_num);
      p = p->prev;  
    } while(p != start);
    neighbors.push_back(neighbors[0]);  // close the ring

    for(std::size_t t = 0; t + 1 < neighbors.size(); ++t) {
      const PointNum j = neighbors[t];
      const PointNum k = neighbors[t + 1];
      if(j > i && k > i && is_right_of(i, j, k)) {
        // Map sorted indices back to original point indices.
        out_triangles.push_back({
          static_cast<std::uint32_t>(sorted_to_orig[static_cast<std::size_t>(i)]),
          static_cast<std::uint32_t>(sorted_to_orig[static_cast<std::size_t>(j)]),
          static_cast<std::uint32_t>(sorted_to_orig[static_cast<std::size_t>(k)])
        });
      }
    }
  }

  // ── Cleanup all allocated AdjNodes ───────────────────────────────
  for(auto *node : all_nodes) {
    delete node;
  }
}


[[nodiscard]] inline std::uint64_t make_edge_key(
  std::uint32_t a, std::uint32_t b
) noexcept
{
  const std::uint32_t lo = a < b ? a : b;
  const std::uint32_t hi = a < b ? b : a;
  return (static_cast<std::uint64_t>(lo) << 32) | static_cast<std::uint64_t>(hi);
}

bool recover_edge(
  const std::vector<std::array<double, 2>> &uvs,
  std::vector<CdtTriangle> &triangles,
  std::uint32_t va, std::uint32_t vb,
  const std::unordered_set<std::uint64_t> *e2r = nullptr,
  std::unordered_set<std::uint64_t> *not_recovered = nullptr
)
{
  if(edge_exists(triangles, va, vb)) return true;

  for(int ix = 0; ix < 300; ++ix) {
    if(edge_exists(triangles, va, vb)) return true;

    // Collect all edges that intersect segment (va, vb).
    struct IntersectingEdge {
      std::uint32_t tri_index;
      std::uint32_t edge_slot;
      std::uint32_t v0;
      std::uint32_t v1;
    };
    std::vector<IntersectingEdge> intersected;
    bool self_intersection = false;

    for(std::uint32_t ti = 0; ti < static_cast<std::uint32_t>(triangles.size()); ++ti) {
      if(!triangles[ti].valid) continue;
      const auto &tri = triangles[ti];
      for(std::uint32_t e = 0; e < 3; ++e) {
        const auto v0 = tri.vertices[e];
        const auto v1 = tri.vertices[(e + 1) % 3];
        if(v0 == va || v0 == vb || v1 == va || v1 == vb) continue;
        if(v0 > v1) continue;
        if(!segments_cross_strict(uvs[v0], uvs[v1], uvs[va], uvs[vb])) continue;

        
        if(e2r != nullptr &&
           e2r->find(make_edge_key(v0, v1)) != e2r->end()) {
          if(not_recovered != nullptr) {
            not_recovered->insert(make_edge_key(va, vb));
            not_recovered->insert(make_edge_key(v0, v1));
          }
          self_intersection = true;
        }
        intersected.push_back({ti, e, v0, v1});
      }
    }

    if(self_intersection) return false;

    if(intersected.empty()) break;  // no intersections — edge recovered or unrecoverable

    bool any_swapped = false;
    for(const auto &ie : intersected) {
      if(!triangles[ie.tri_index].valid) continue;
      const auto v0 = triangles[ie.tri_index].vertices[ie.edge_slot];
      const auto v1 = triangles[ie.tri_index].vertices[(ie.edge_slot + 1) % 3];
      if(v0 == va || v0 == vb || v1 == va || v1 == vb) continue;
      if(!segments_cross_strict(uvs[v0], uvs[v1], uvs[va], uvs[vb])) continue;

      if(try_cdt_swap(uvs, triangles, ie.tri_index, ie.edge_slot)) {
        any_swapped = true;
        break;  
      }
    }

    if(!any_swapped) break;  // all intersecting edges are non-swappable
  }

  const bool ok = edge_exists(triangles, va, vb);
  if(!ok && not_recovered != nullptr) {
    not_recovered->insert(make_edge_key(va, vb));
  }
  return ok;
}

// Check if a point is inside a polygon (ray-casting).
[[nodiscard]] bool point_in_polygon(
  const std::vector<std::array<double, 2>> &polygon,
  const std::array<double, 2> &p
) noexcept
{
  bool inside = false;
  const auto n = polygon.size();
  for(std::size_t i = 0, j = n - 1; i < n; j = i++) {
    const auto &vi = polygon[i];
    const auto &vj = polygon[j];
    if(((vi[1] > p[1]) != (vj[1] > p[1])) &&
       (p[0] < (vj[0] - vi[0]) * (p[1] - vi[1]) / (vj[1] - vi[1]) + vi[0])) {
      inside = !inside;
    }
  }
  return inside;
}

} // namespace cdt_detail

// Build a constrained Delaunay triangulation of boundary points in UV space.
// Returns the triangles as (vertex_index, vertex_index, vertex_index) tuples
// referring to indices in plan.vertices. Only produces triangles inside the
// boundary polygon.
[[nodiscard]] bool build_boundary_cdt(
  const std::vector<LocalVertex> &/*vertices*/,
  const std::vector<std::uint32_t> &boundary_vertex_indices,
  const std::vector<std::array<double, 2>> &boundary_uvs,
  std::vector<std::array<std::uint32_t, 3>> &out_triangles
)
{
  const auto n = boundary_vertex_indices.size();
  if(n < 3U) return false;

  
  double umin = boundary_uvs[0][0], umax = umin;
  double vmin = boundary_uvs[0][1], vmax = vmin;
  for(const auto &uv : boundary_uvs) {
    umin = std::min(umin, uv[0]); umax = std::max(umax, uv[0]);
    vmin = std::min(vmin, uv[1]); vmax = std::max(vmax, uv[1]);
  }
  const double du = umax - umin;
  const double dv = vmax - vmin;
  
  const double lc_2d = std::sqrt(du * du + dv * dv);
  constexpr double kRandFactor = 1.0e-9;
  std::srand(42U);  // deterministic seed for reproducibility

  // Create UV array with perturbation.
  std::vector<std::array<double, 2>> uvs;
  uvs.reserve(n + 4U);
  for(std::size_t i = 0; i < n; ++i) {
    auto uv = boundary_uvs[i];
    uv[0] += kRandFactor * lc_2d * static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);
    uv[1] += kRandFactor * lc_2d * static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);
    uvs.push_back(uv);
  }

  
  const auto real_n = static_cast<std::uint32_t>(uvs.size());
  {
    double bb_du = du, bb_dv = dv;
    if(du > 1e-30 && dv > 1e-30 && du / dv < 1200.0 && dv / du < 1200.0) {
      bb_du = bb_dv = std::max(du, dv);  // makeCube
    }
    const double hdu = bb_du * 3.5 * 0.5;
    const double hdv = bb_dv * 3.5 * 0.5;
    const double cu = (umin + umax) * 0.5;
    const double cv = (vmin + vmax) * 0.5;
    uvs.push_back({cu - hdu, cv - hdv});
    uvs.push_back({cu - hdu, cv + hdv});
    uvs.push_back({cu + hdu, cv - hdv});
    uvs.push_back({cu + hdu, cv + hdv});
  }

  // Divide-and-conquer Delaunay triangulation.
  std::vector<std::array<std::uint32_t, 3>> dc_triangles;
  cdt_detail::divide_and_conquer_triangulate(uvs, dc_triangles);

  if(dc_triangles.empty()) return false;

  std::vector<cdt_detail::CdtTriangle> triangles;
  triangles.reserve(dc_triangles.size());
  for(const auto &tri : dc_triangles) {
    cdt_detail::CdtTriangle ct;
    ct.vertices = {tri[0], tri[1], tri[2]};
    ct.valid = true;
    triangles.push_back(ct);
  }
  cdt_detail::rebuild_adjacency(triangles);

  // Recover constrained boundary edges in the FULL triangulation.
  for(std::uint32_t i = 0; i < static_cast<std::uint32_t>(n); ++i) {
    const std::uint32_t va = i;
    const std::uint32_t vb = (i + 1) % static_cast<std::uint32_t>(n);
    cdt_detail::recover_edge(uvs, triangles, va, vb);
  }

  // Remove bbox vertices and exterior triangles.
  uvs.resize(real_n);
  for(auto &tri : triangles) {
    if(!tri.valid) continue;
    if(tri.vertices[0] >= real_n || tri.vertices[1] >= real_n || tri.vertices[2] >= real_n) {
      tri.valid = false;
      continue;
    }
    const auto &a = uvs[tri.vertices[0]];
    const auto &b = uvs[tri.vertices[1]];
    const auto &c = uvs[tri.vertices[2]];
    const std::array<double, 2> centroid {
      (a[0] + b[0] + c[0]) / 3.0,
      (a[1] + b[1] + c[1]) / 3.0,
    };
    if(!cdt_detail::point_in_polygon(boundary_uvs, centroid)) {
      tri.valid = false;
    }
  }

  // Collect valid triangles, mapping CDT vertex indices back to
  // boundary_vertex_indices.
  out_triangles.clear();
  for(const auto &tri : triangles) {
    if(!tri.valid) continue;
    std::array<std::uint32_t, 3> mapped {
      boundary_vertex_indices[tri.vertices[0]],
      boundary_vertex_indices[tri.vertices[1]],
      boundary_vertex_indices[tri.vertices[2]],
    };
    out_triangles.push_back(mapped);
  }

  return !out_triangles.empty();
}

[[nodiscard]] bool run_pure_delaunay_cdt(
  const std::vector<std::array<double, 2>> &uvs,
  std::vector<cdt_detail::CdtTriangle> &out_triangles
)
{
  std::vector<std::array<std::uint32_t, 3>> dc_triangles;
  cdt_detail::divide_and_conquer_triangulate(uvs, dc_triangles);
  if(dc_triangles.empty()) return false;

  out_triangles.clear();
  out_triangles.reserve(dc_triangles.size());
  for(const auto &tri : dc_triangles) {
    cdt_detail::CdtTriangle ct;
    ct.vertices = {tri[0], tri[1], tri[2]};
    ct.valid = true;
    out_triangles.push_back(ct);
  }
  cdt_detail::rebuild_adjacency(out_triangles);
  return true;
}


void classify_interior_and_drop_exterior(
  std::vector<cdt_detail::CdtTriangle> &triangles,
  const std::unordered_set<std::uint64_t> &constrained_edge_keys,
  std::uint32_t real_point_count
)
{
  const auto tri_touches_bbox = [&](const cdt_detail::CdtTriangle &t) noexcept {
    return t.vertices[0] >= real_point_count ||
           t.vertices[1] >= real_point_count ||
           t.vertices[2] >= real_point_count;
  };

  constexpr std::uint8_t kTagNone = 0U;
  constexpr std::uint8_t kTagExterior = 1U;
  constexpr std::uint8_t kTagInterior = 2U;
  std::vector<std::uint8_t> tags(triangles.size(), kTagNone);

  const auto flood = [&](std::uint32_t seed, std::uint8_t target) {
    if(seed >= triangles.size() || !triangles[seed].valid ||
       tags[seed] != kTagNone) {
      return;
    }
    std::vector<std::uint32_t> stack;
    stack.push_back(seed);
    while(!stack.empty()) {
      const std::uint32_t cur = stack.back();
      stack.pop_back();
      if(tags[cur] != kTagNone) continue;
      tags[cur] = target;
      const auto &ct = triangles[cur];
      for(std::uint32_t e = 0U; e < 3U; ++e) {
        const std::uint32_t nb = ct.neighbors[e];
        if(nb == invalid_index || nb >= triangles.size()) continue;
        if(!triangles[nb].valid) continue;
        const std::uint32_t v0 = ct.vertices[e];
        const std::uint32_t v1 = ct.vertices[(e + 1U) % 3U];
        if(constrained_edge_keys.count(cdt_detail::make_edge_key(v0, v1)) != 0U)
          continue;
        if(tags[nb] == kTagNone) stack.push_back(nb);
      }
    }
  };

  // E1: seed EXTERIOR from the first bbox-touching triangle found.
  for(std::uint32_t ti = 0U; ti < triangles.size(); ++ti) {
    if(!triangles[ti].valid) continue;
    if(tri_touches_bbox(triangles[ti])) {
      flood(ti, kTagExterior);
      break;
    }
  }

  // E2: seed INTERIOR from the non-EXTERIOR side of any constrained
  // edge whose other side is EXTERIOR (i.e. jump across the outer loop).
  for(std::uint32_t ti = 0U; ti < triangles.size(); ++ti) {
    if(!triangles[ti].valid) continue;
    const auto &ct = triangles[ti];
    bool seeded_interior = false;
    for(std::uint32_t e = 0U; e < 3U; ++e) {
      const std::uint32_t nb = ct.neighbors[e];
      if(nb == invalid_index || nb >= triangles.size()) continue;
      if(!triangles[nb].valid) continue;
      const std::uint32_t v0 = ct.vertices[e];
      const std::uint32_t v1 = ct.vertices[(e + 1U) % 3U];
      if(constrained_edge_keys.count(cdt_detail::make_edge_key(v0, v1)) == 0U)
        continue;
      const std::uint8_t a = tags[ti];
      const std::uint8_t b = tags[nb];
      if(a == kTagExterior && b == kTagNone) {
        flood(nb, kTagInterior);
        seeded_interior = true;
        break;
      }
      if(b == kTagExterior && a == kTagNone) {
        flood(ti, kTagInterior);
        seeded_interior = true;
        break;
      }
    }
    if(seeded_interior) break;
  }

  // E4: drop every triangle not tagged INTERIOR (exterior + hole +
  // any stray bbox triangles the EXTERIOR seed didn't reach).
  for(std::uint32_t ti = 0U; ti < triangles.size(); ++ti) {
    if(!triangles[ti].valid) continue;
    if(tags[ti] != kTagInterior) {
      triangles[ti].valid = false;
    }
  }
}


[[nodiscard]] bool build_boundary_cdt_with_holes(
  const std::vector<LocalVertex> &/*vertices*/,
  const SeedLoopData &outer_loop,
  const std::vector<SeedLoopData> &inner_loops,
  std::vector<std::array<std::uint32_t, 3>> &out_triangles,
  std::unordered_set<std::uint64_t> *out_not_recovered
)
{
  if(outer_loop.vertices.size() < 3U || outer_loop.uvs.size() < 3U) {
    return false;
  }

  // Build combined UV array: NO super-triangle needed for divide-and-conquer.
  // Track which original vertex index each CDT vertex maps to.
  struct LoopRange {
    std::uint32_t cdt_start = 0U;
    std::uint32_t count = 0U;
  };

  std::vector<std::array<double, 2>> uvs;
  std::vector<std::uint32_t> cdt_to_original;  // CDT index -> original vertex index

  // Reserve space.
  std::size_t total_points = outer_loop.vertices.size();
  for(const auto &inner : inner_loops) {
    total_points += inner.vertices.size();
  }
  uvs.reserve(total_points);
  cdt_to_original.reserve(total_points);

  // Compute UV bounding box over ALL loops.
  double umin = outer_loop.uvs[0][0], umax = umin;
  double vmin = outer_loop.uvs[0][1], vmax = vmin;
  for(const auto &uv : outer_loop.uvs) {
    umin = std::min(umin, uv[0]); umax = std::max(umax, uv[0]);
    vmin = std::min(vmin, uv[1]); vmax = std::max(vmax, uv[1]);
  }
  for(const auto &inner : inner_loops) {
    for(const auto &uv : inner.uvs) {
      umin = std::min(umin, uv[0]); umax = std::max(umax, uv[0]);
      vmin = std::min(vmin, uv[1]); vmax = std::max(vmax, uv[1]);
    }
  }
  const double du = umax - umin;
  const double dv = vmax - vmin;

  // Normalize UV coordinates to unit square before CDT. Extreme aspect
  // ratios (e.g. the 56:1 cylinder on missile.step) are handled natively
  // by divide-and-conquer, but normalization still helps numerical
  // conditioning for the incircle test.
  const double inv_du = du > 1.0e-30 ? (1.0 / du) : 1.0;
  const double inv_dv = dv > 1.0e-30 ? (1.0 / dv) : 1.0;
  const auto normalize_uv = [&](std::array<double, 2> uv) noexcept {
    return std::array<double, 2> {(uv[0] - umin) * inv_du, (uv[1] - vmin) * inv_dv};
  };

  
  constexpr double kRandFactor = 1.0e-9;
  const double lc_2d = std::sqrt(2.0);  // diagonal of unit square
  std::srand(42U);  // deterministic seed for reproducibility

  // Outer loop points (starting at index 0, no super-triangle offset).
  LoopRange outer_range;
  outer_range.cdt_start = static_cast<std::uint32_t>(uvs.size());
  outer_range.count = static_cast<std::uint32_t>(outer_loop.vertices.size());
  for(std::size_t i = 0U; i < outer_loop.vertices.size(); ++i) {
    auto uv = normalize_uv(outer_loop.uvs[i]);
    uv[0] += kRandFactor * lc_2d * static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);
    uv[1] += kRandFactor * lc_2d * static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);
    uvs.push_back(uv);
    cdt_to_original.push_back(outer_loop.vertices[i]);
  }

  // Inner loop points.
  std::vector<LoopRange> inner_ranges;
  inner_ranges.reserve(inner_loops.size());
  for(const auto &inner : inner_loops) {
    LoopRange range;
    range.cdt_start = static_cast<std::uint32_t>(uvs.size());
    range.count = static_cast<std::uint32_t>(inner.vertices.size());
    for(std::size_t i = 0U; i < inner.vertices.size(); ++i) {
      auto uv = normalize_uv(inner.uvs[i]);
      uv[0] += kRandFactor * lc_2d * static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);
      uv[1] += kRandFactor * lc_2d * static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);
      uvs.push_back(uv);
      cdt_to_original.push_back(inner.vertices[i]);
    }
    inner_ranges.push_back(range);
  }

  
  const auto real_point_count = static_cast<std::uint32_t>(uvs.size());
  {
    
    const double bbox_lo = -2.5;
    const double bbox_hi =  3.5;
    uvs.push_back({bbox_lo, bbox_lo});
    uvs.push_back({bbox_lo, bbox_hi});
    uvs.push_back({bbox_hi, bbox_lo});
    uvs.push_back({bbox_hi, bbox_hi});
  }

  
  std::vector<cdt_detail::CdtTriangle> triangles;
  if(!run_pure_delaunay_cdt(uvs, triangles)) {
    return false;
  }

  std::unordered_set<std::uint64_t> e2r_set;
  const auto total_constraint_count =
    outer_range.count + [&] {
      std::uint32_t n = 0U;
      for(const auto &r : inner_ranges) n += r.count;
      return n;
    }();
  e2r_set.reserve(total_constraint_count);
  for(std::uint32_t i = 0U; i < outer_range.count; ++i) {
    const std::uint32_t va = outer_range.cdt_start + i;
    const std::uint32_t vb = outer_range.cdt_start + ((i + 1) % outer_range.count);
    e2r_set.insert(cdt_detail::make_edge_key(va, vb));
  }
  for(const auto &range : inner_ranges) {
    for(std::uint32_t i = 0U; i < range.count; ++i) {
      const std::uint32_t va = range.cdt_start + i;
      const std::uint32_t vb = range.cdt_start + ((i + 1) % range.count);
      e2r_set.insert(cdt_detail::make_edge_key(va, vb));
    }
  }

  std::unordered_set<std::uint64_t> not_recovered;
  for(std::uint32_t i = 0U; i < outer_range.count; ++i) {
    const std::uint32_t va = outer_range.cdt_start + i;
    const std::uint32_t vb = outer_range.cdt_start + ((i + 1) % outer_range.count);
    cdt_detail::recover_edge(uvs, triangles, va, vb, &e2r_set, &not_recovered);
  }
  for(const auto &range : inner_ranges) {
    for(std::uint32_t i = 0U; i < range.count; ++i) {
      const std::uint32_t va = range.cdt_start + i;
      const std::uint32_t vb = range.cdt_start + ((i + 1) % range.count);
      cdt_detail::recover_edge(uvs, triangles, va, vb, &e2r_set, &not_recovered);
    }
  }

  // Stage 3: topological flood-fill classification (reuses e2r_set as the
  // constrained-edge set — same set used as foresight above).
  classify_interior_and_drop_exterior(triangles, e2r_set, real_point_count);

  // Bbox vertex slots are no longer referenced; trim the UV array.
  uvs.resize(real_point_count);

  // Collect valid triangles, mapping CDT indices back to original vertex indices.
  out_triangles.clear();
  for(const auto &tri : triangles) {
    if(!tri.valid) continue;
    std::array<std::uint32_t, 3> mapped {
      cdt_to_original[tri.vertices[0]],
      cdt_to_original[tri.vertices[1]],
      cdt_to_original[tri.vertices[2]],
    };
    out_triangles.push_back(mapped);
  }

  // Unrecovered constraint edges, rewritten in terms of original vertex
  // indices. The retry loop splits these failed 1D edges at their
  // midpoint and restarts meshing.
  if(out_not_recovered != nullptr) {
    out_not_recovered->clear();
    out_not_recovered->reserve(not_recovered.size());
    for(const auto key : not_recovered) {
      const auto lo = static_cast<std::uint32_t>(key >> 32);
      const auto hi = static_cast<std::uint32_t>(key & 0xFFFFFFFFull);
      const std::uint32_t a = cdt_to_original[lo];
      const std::uint32_t b = cdt_to_original[hi];
      out_not_recovered->insert(cdt_detail::make_edge_key(a, b));
    }
  }

  return !out_triangles.empty();
}

[[nodiscard]] base::StatusCode build_nonplanar_single_loop_patch_seed_mesh(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  FaceMeshPlan &plan
)
{
  if(face_preprocess.boundary.has_seams ||
     face_preprocess.boundary.outer_loop_count != 1U ||
     face_preprocess.boundary.inner_loop_count != 0U ||
     face_preprocess.boundary.unknown_loop_count != 0U ||
     face_preprocess.loops.size() != 1U ||
     !face_preprocess.loops.front().closed ||
     !face_preprocess.loops.front().continuous) {
    return core::detail::publish_error(
      base::StatusCode::unsupported,
      "Auto CFD Surface Mesher curved single-loop patch seeding requires one closed, continuous non-seam outer loop."
    );
  }

  std::unordered_map<std::uint32_t, std::uint32_t> boundary_to_local;
  boundary_to_local.reserve(face_preprocess.loops.front().points.size() * 2U);

  SeedLoopData outer_loop;
  LocalEdgeOwnerMap boundary_edge_owners;
  build_boundary_node_edge_owner_map(
    face_preprocess.loops.front().segments,
    boundary_edge_owners
  );
  auto status = append_seed_loop_from_points(
    face_preprocess.loops.front().points,
    geo::FaceBoundaryLoopKind::outer,
    boundary_to_local,
    plan.vertices,
    outer_loop,
    false,
    &boundary_edge_owners
  );
  if(status != base::StatusCode::ok) {
    return status;
  }

  // Preserve the recovered boundary walk order.  The previous "remove any
  // duplicate vertex anywhere in the loop" pass could rewire non-adjacent
  // boundary samples into a new edge that never existed on the recovered trim,
  // which in turn drops constrained-edge owners and weakens diagnostics.
  normalize_seed_loop_cycle(outer_loop);
  if(outer_loop.vertices.size() < 3U) {
    return core::detail::publish_error(
      base::StatusCode::unsupported,
      "Auto CFD Surface Mesher curved single-loop patch seeding requires at least three unique boundary vertices after duplicate-loop compression."
    );
  }

  populate_seed_loop_edge_owners(plan.vertices, &boundary_edge_owners, outer_loop);
  add_constrained_loop_edges(
    outer_loop,
    plan.constrained_edges,
    &plan.constrained_edge_owners
  );

  
  std::vector<std::array<std::uint32_t, 3>> triangles;
  if(!build_boundary_cdt(
       plan.vertices,
       outer_loop.vertices,
       outer_loop.uvs,
       triangles)) {
    return core::detail::publish_error(
      base::StatusCode::unsupported,
      "Auto CFD Surface Mesher curved single-loop CDT seeding failed."
    );
  }

  plan.triangles.reserve(triangles.size());
  for(auto triangle : triangles) {
    const double area = signed_triangle_area_twice(
      plan.vertices[triangle[0]].uv,
      plan.vertices[triangle[1]].uv,
      plan.vertices[triangle[2]].uv
    );
    if(area < 0.0) {
      std::swap(triangle[1], triangle[2]);
    }
    plan.triangles.push_back({triangle, true, false});
  }

  if(plan.triangles.empty()) {
    return core::detail::publish_error(
      base::StatusCode::unsupported,
      "Auto CFD Surface Mesher curved single-loop patch seeding produced no non-degenerate triangles."
    );
  }

  return core::detail::clear_error_state();
}

void wire_periodic_counterparts_from_boundary_nodes(FaceMeshPlan &plan) noexcept
{
  std::unordered_map<std::uint32_t, std::uint32_t> first_by_boundary_node;
  first_by_boundary_node.reserve(plan.vertices.size());

  for(std::uint32_t vertex_index = 0U;
      vertex_index < plan.vertices.size();
      ++vertex_index) {
    auto &vertex = plan.vertices[vertex_index];
    vertex.periodic_counterpart = invalid_index;
    if(vertex.boundary_node_index == invalid_index) {
      continue;
    }
  }

  for(std::uint32_t vertex_index = 0U;
      vertex_index < plan.vertices.size();
      ++vertex_index) {
    auto &vertex = plan.vertices[vertex_index];
    if(vertex.boundary_node_index == invalid_index) {
      continue;
    }

    const auto it = first_by_boundary_node.find(vertex.boundary_node_index);
    if(it != first_by_boundary_node.end()) {
      const auto previous_vertex_index = it->second;
      if(previous_vertex_index < plan.vertices.size()) {
        plan.vertices[previous_vertex_index].periodic_counterpart =
          vertex_index;
        vertex.periodic_counterpart = previous_vertex_index;
      }
    }
    
    first_by_boundary_node[vertex.boundary_node_index] = vertex_index;
  }
}


[[nodiscard]] base::StatusCode append_consecutive_list_of_vertices(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceFaceLoopState &loop,
  FaceMeshPlan &plan,
  SeedLoopData &out_polygon,
  bool seam_the_first = false
)
{
  out_polygon = {};
  out_polygon.kind = loop.kind;

  const auto same_3d =
    [](const LocalVertex &lhs,
       const geo::Point3 &rhs_pos) noexcept {
      constexpr double kPositionTolSq = 1.0e-20;
      const double dx = lhs.position[0] - rhs_pos[0];
      const double dy = lhs.position[1] - rhs_pos[1];
      const double dz = lhs.position[2] - rhs_pos[2];
      return (dx*dx + dy*dy + dz*dz) <= kPositionTolSq;
    };

  const auto uv_dist_sq = [](const std::array<double, 2> &a,
                              const std::array<double, 2> &b) noexcept {
    const double du = a[0] - b[0], dv = a[1] - b[1];
    return du * du + dv * dv;
  };

  
  const auto get_alt_uv = [&](const AutoCfdSurfaceFaceLoopSegment &seg,
                               const AutoCfdSurfaceFaceLoopPoint &pt,
                               std::array<double, 2> &alt_uv) -> bool {
    geo::FaceBoundaryEdgeUse alt_eu = seg.edge_use;
    alt_eu.same_orientation_as_edge = !seg.edge_use.same_orientation_as_edge;
    geo::FaceUvMapping alt_map;
    if(geo::recover_face_uv_from_edge_use(face_view, alt_eu, pt.parameter, alt_map)
       == base::StatusCode::ok) {
      alt_uv = {alt_map.u, alt_map.v};
      return true;
    }
    return false;
  };

  // Helper: emit one point into the polygon.
  const auto emit_point = [&](const AutoCfdSurfaceFaceLoopPoint &pt,
                               const std::array<double, 2> &uv,
                               const geo::TopologyEntityId &owner) {
    
    if(!out_polygon.uvs.empty()) {
      constexpr double kJunctionUvTolSq = 1.0e-20;
      if(uv_dist_sq(out_polygon.uvs.back(), uv) <= kJunctionUvTolSq) {
        return;
      }
    }
    const auto idx = static_cast<std::uint32_t>(plan.vertices.size());
    plan.vertices.push_back({pt.node_index, pt.position, uv, true, 0.0, 0.0, {}});
    out_polygon.vertices.push_back(idx);
    out_polygon.uvs.push_back(uv);
    out_polygon.edge_owners.push_back(owner);
  };

  for(std::size_t seg_i = 0; seg_i < loop.segments.size(); ++seg_i) {
    const auto &segment = loop.segments[seg_i];
    if(segment.points.empty()) continue;

    for(const auto &pt : segment.points) {
      if(!pt.uv_defined || pt.node_index == invalid_index) {
        return core::detail::publish_error(
          base::StatusCode::unsupported,
          "Auto CFD Surface Mesher periodic face seeding requires every "
          "segment point to have a defined UV and boundary node index."
        );
      }
    }

    const bool is_seam = segment.edge_use.is_seam;
    const auto n_pts = segment.points.size();

    
    std::vector<std::array<double, 2>> p(n_pts);
    std::vector<std::array<double, 2>> p_alt;
    for(std::size_t k = 0; k < n_pts; ++k) {
      p[k] = segment.points[k].uv;
    }
    bool have_alt = false;
    if(is_seam && geo::is_valid(segment.edge_use.edge)) {
      p_alt.resize(n_pts);
      have_alt = true;
      for(std::size_t k = 0; k < n_pts; ++k) {
        if(!get_alt_uv(segment, segment.points[k], p_alt[k])) {
          have_alt = false;
          break;
        }
      }
    }

    
    const std::vector<std::array<double, 2>> *chosen_uvs = &p;
    bool chosen_reverse = false;

    if(seg_i == 0) {
      if(seam_the_first) {
        // Reverse the first segment's traversal order.
        chosen_reverse = true;
      }
    }
    else if(!out_polygon.uvs.empty()) {
      
      const auto &prev_uv = out_polygon.uvs.back();

      double best_dist = uv_dist_sq(prev_uv, p.front());
      chosen_uvs = &p;
      chosen_reverse = false;

      // Try reversed.
      double d = uv_dist_sq(prev_uv, p.back());
      if(d < best_dist) { best_dist = d; chosen_uvs = &p; chosen_reverse = true; }

      if(have_alt) {
        // Try alternate forward.
        d = uv_dist_sq(prev_uv, p_alt.front());
        if(d < best_dist) { best_dist = d; chosen_uvs = &p_alt; chosen_reverse = false; }
        // Try alternate reversed.
        d = uv_dist_sq(prev_uv, p_alt.back());
        if(d < best_dist) { best_dist = d; chosen_uvs = &p_alt; chosen_reverse = true; }
      }
    }

    // Emit points in chosen order.
    if(chosen_reverse) {
      for(std::size_t k = n_pts; k > 0; --k) {
        emit_point(segment.points[k - 1], (*chosen_uvs)[k - 1], segment.edge_use.edge);
      }
    }
    else {
      for(std::size_t k = 0; k < n_pts; ++k) {
        emit_point(segment.points[k], (*chosen_uvs)[k], segment.edge_use.edge);
      }
    }
  }

  
  if(out_polygon.uvs.size() >= 2U) {
    const double close_dist_sq =
      uv_dist_sq(out_polygon.uvs.front(), out_polygon.uvs.back());
    double u_lo = out_polygon.uvs[0][0], u_hi = u_lo;
    double v_lo = out_polygon.uvs[0][1], v_hi = v_lo;
    for(const auto &uv : out_polygon.uvs) {
      u_lo = std::min(u_lo, uv[0]); u_hi = std::max(u_hi, uv[0]);
      v_lo = std::min(v_lo, uv[1]); v_hi = std::max(v_hi, uv[1]);
    }
    const double lc2d = std::sqrt((u_hi-u_lo)*(u_hi-u_lo) + (v_hi-v_lo)*(v_hi-v_lo));
    const double tol = 1.0e-6 * lc2d;
    const double tol_sq = tol * tol;

    if(close_dist_sq < tol_sq) {
      out_polygon.vertices.pop_back();
      out_polygon.uvs.pop_back();
      if(!out_polygon.edge_owners.empty()) {
        out_polygon.edge_owners.pop_back();
      }
    }
    else {
      return core::detail::publish_error(
        base::StatusCode::unsupported,
        "Auto CFD Surface Mesher periodic polygon UV did not close "
        "(distance too large). Retrying with seam_the_first."
      );
    }
  }

  if(out_polygon.vertices.size() < 3U) {
    return core::detail::publish_error(
      base::StatusCode::unsupported,
      "Auto CFD Surface Mesher periodic face seeding produced fewer than three polygon vertices from the edge loop."
    );
  }

  return core::detail::clear_error_state();
}

[[nodiscard]] base::StatusCode build_periodic_face_seed_mesh(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  FaceMeshPlan &plan
)
{
  if(face_preprocess.loops.empty()) {
    return core::detail::publish_error(
      base::StatusCode::unsupported,
      "Auto CFD Surface Mesher periodic face seeding received no boundary loops."
    );
  }

  
  std::vector<SeedLoopData> all_polygons;
  all_polygons.reserve(face_preprocess.loops.size());
  std::size_t total_segment_count = 0U;

  for(const auto &loop : face_preprocess.loops) {
    if(!loop.closed) {
      return core::detail::publish_error(
        base::StatusCode::unsupported,
        "Auto CFD Surface Mesher periodic face seeding requires every boundary loop to be closed."
      );
    }
    if(loop.segments.empty()) {
      continue;
    }

   
    SeedLoopData polygon;
    auto saved_plan_vertices_size = plan.vertices.size();
    auto status = append_consecutive_list_of_vertices(
      face_view, loop, plan, polygon, /*seam_the_first=*/false
    );
    if(status != base::StatusCode::ok) {
      // Retry with seam_the_first=true.
      plan.vertices.resize(saved_plan_vertices_size);
      polygon = {};
      status = append_consecutive_list_of_vertices(
        face_view, loop, plan, polygon, /*seam_the_first=*/true
      );
      if(status != base::StatusCode::ok) {
        return status;
      }
    }
    total_segment_count += loop.segments.size();
    all_polygons.push_back(std::move(polygon));
  }

  if(all_polygons.empty()) {
    return core::detail::publish_error(
      base::StatusCode::unsupported,
      "Auto CFD Surface Mesher periodic face seeding produced no valid loop polygons."
    );
  }

  
  plan.periodic_true_boundary_segments.clear();
  for(const auto &polygon : all_polygons) {
    for(std::size_t i = 0U; i < polygon.uvs.size(); ++i) {
      const std::size_t next = (i + 1U) % polygon.uvs.size();
      plan.periodic_true_boundary_segments.push_back(polygon.uvs[i]);
      plan.periodic_true_boundary_segments.push_back(polygon.uvs[next]);
    }
  }
  const auto &uv_bounds = face_preprocess.uv_bounds;
  plan.periodic_true_boundary_far = {
    uv_bounds.u_max + (uv_bounds.u_max - uv_bounds.u_min),
    uv_bounds.v_max + (uv_bounds.v_max - uv_bounds.v_min),
  };
  plan.periodic_true_boundary_defined = true;

  
  wire_periodic_counterparts_from_boundary_nodes(plan);

  std::size_t counterpart_pair_count = 0U;
  for(std::uint32_t vi = 0U; vi < plan.vertices.size(); ++vi) {
    const auto cp = plan.vertices[vi].periodic_counterpart;
    if(cp != invalid_index && cp > vi) {
      ++counterpart_pair_count;
    }
  }

  std::size_t total_polygon_vertices = 0U;
  for(const auto &p : all_polygons) {
    total_polygon_vertices += p.vertices.size();
  }

  // Diagnostic: dump UV bounds of each polygon to verify seam walk.
  for(std::size_t pi = 0U; pi < all_polygons.size(); ++pi) {
    const auto &polygon = all_polygons[pi];
    if(polygon.uvs.empty()) continue;
    double u_lo = polygon.uvs[0][0], u_hi = u_lo;
    double v_lo = polygon.uvs[0][1], v_hi = v_lo;
    for(const auto &uv : polygon.uvs) {
      u_lo = std::min(u_lo, uv[0]); u_hi = std::max(u_hi, uv[0]);
      v_lo = std::min(v_lo, uv[1]); v_hi = std::max(v_hi, uv[1]);
    }
  }

  // Record boundary edges as constrained for ALL loops.
  for(const auto &polygon : all_polygons) {
    add_constrained_loop_edges(
      polygon,
      plan.constrained_edges,
      &plan.constrained_edge_owners
    );
  }

  
  SeedLoopData outer_loop;
  std::vector<SeedLoopData> inner_loops;

  // Find outer loop: prefer the loop marked as `outer`, otherwise
  // fall back to the primary_outer_loop_index from preprocessing,
  // otherwise use the first loop.
  std::size_t outer_index = 0U;
  for(std::size_t i = 0U; i < all_polygons.size(); ++i) {
    if(all_polygons[i].kind == geo::FaceBoundaryLoopKind::outer) {
      outer_index = i;
      break;
    }
  }
  if(outer_index == 0U &&
     face_preprocess.boundary.primary_outer_loop_index < all_polygons.size()) {
    outer_index = face_preprocess.boundary.primary_outer_loop_index;
  }

  outer_loop = all_polygons[outer_index];
  inner_loops.reserve(all_polygons.size());
  for(std::size_t i = 0U; i < all_polygons.size(); ++i) {
    if(i != outer_index) {
      inner_loops.push_back(all_polygons[i]);
    }
  }

  enforce_outer_inner_loop_winding(outer_loop, inner_loops);

  std::vector<std::array<std::uint32_t, 3>> triangles;
  if(inner_loops.empty()) {
    // Single-loop: use the simpler build_boundary_cdt.
    if(!build_boundary_cdt(
         plan.vertices,
         outer_loop.vertices,
         outer_loop.uvs,
         triangles
       )) {
      return core::detail::publish_error(
        base::StatusCode::unsupported,
        "Auto CFD Surface Mesher periodic face seeding could not build a boundary-constrained Delaunay triangulation from the loop polygon."
      );
    }
  }
  else {
    // Multi-loop: use build_boundary_cdt_with_holes which handles
    // outer + inner (hole) loops natively.
    if(!build_boundary_cdt_with_holes(
         plan.vertices,
         outer_loop,
         inner_loops,
         triangles
       )) {
      return core::detail::publish_error(
        base::StatusCode::unsupported,
        "Auto CFD Surface Mesher periodic face seeding could not build a multi-loop CDT."
      );
    }
  }

  plan.triangles.reserve(triangles.size());
  for(auto triangle : triangles) {
    const double area = signed_triangle_area_twice(
      plan.vertices[triangle[0]].uv,
      plan.vertices[triangle[1]].uv,
      plan.vertices[triangle[2]].uv
    );
    if(area < 0.0) {
      std::swap(triangle[1], triangle[2]);
    }
    plan.triangles.push_back({triangle, true, false});
  }

  if(plan.triangles.empty()) {
    return core::detail::publish_error(
      base::StatusCode::unsupported,
      "Auto CFD Surface Mesher periodic face seeding produced no non-degenerate CDT triangles."
    );
  }

  return core::detail::clear_error_state();
}

[[nodiscard]] base::StatusCode build_planar_fallback_seed_mesh(
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  FaceMeshPlan &plan
)
{
  if(face_preprocess.boundary.has_seams) {
    return core::detail::publish_error(
      base::StatusCode::unsupported,
      "Auto CFD Surface Mesher does not broaden planar fallback seeding onto seam-bearing faces."
    );
  }
  if(face_preprocess.boundary.outer_loop_count != 1U ||
     face_preprocess.boundary.unknown_loop_count != 0U) {
    return core::detail::publish_error(
      base::StatusCode::unsupported,
      "Auto CFD Surface Mesher planar fallback currently supports only one classified outer loop and zero unknown loops."
    );
  }
  if(!build_planar_fallback_basis(face_preprocess, plan)) {
    return core::detail::publish_error(
      base::StatusCode::unsupported,
      "Auto CFD Surface Mesher requires a bounded planar face to seed fallback UV coordinates."
    );
  }

  const double plane_tolerance = projection_tolerance(plan.target_size) * 100.0;
  for(const auto &loop : face_preprocess.loops) {
    for(const auto &point : loop.points) {
      const auto delta = point_subtract(point.position, plan.planar_origin);
      if(std::abs(dot_product(delta, plan.reference_normal)) > plane_tolerance) {
        return core::detail::publish_error(
          base::StatusCode::unsupported,
          "Auto CFD Surface Mesher planar fallback is currently limited to near-planar faces."
        );
      }
    }
  }

  std::unordered_map<std::uint32_t, std::uint32_t> boundary_to_local;
  boundary_to_local.reserve(face_preprocess.loops.size() * 8U);

  std::vector<SeedLoopData> seed_loops;
  seed_loops.reserve(face_preprocess.loops.size());
  for(std::size_t loop_index = 0U; loop_index < face_preprocess.loops.size(); ++loop_index) {
    const auto &loop = face_preprocess.loops[loop_index];
    if(!loop.closed || !loop.continuous) {
      return core::detail::publish_error(
        base::StatusCode::unsupported,
        "Auto CFD Surface Mesher planar fallback requires closed, continuous face trim loops."
      );
    }

    SeedLoopData seed_loop;
    LocalEdgeOwnerMap boundary_edge_owners;
    build_boundary_node_edge_owner_map(loop.segments, boundary_edge_owners);
    seed_loop.kind = face_preprocess.boundary.loops[loop_index].kind;
    seed_loop.vertices.reserve(loop.points.size());
    seed_loop.uvs.reserve(loop.points.size());
    for(const auto &point : loop.points) {
      if(point.node_index == invalid_index) {
        return core::detail::publish_error(
          base::StatusCode::unsupported,
          "Auto CFD Surface Mesher planar fallback requires recovered boundary-node ids."
        );
      }

      auto local_it = boundary_to_local.find(point.node_index);
      if(local_it == boundary_to_local.end()) {
        const auto local_index = static_cast<std::uint32_t>(plan.vertices.size());
        plan.vertices.push_back(
          {
            point.node_index,
            point.position,
            planar_local_uv(plan, point.position),
            true,
          }
        );
        local_it = boundary_to_local.emplace(point.node_index, local_index).first;
      }

      if(!seed_loop.vertices.empty() && seed_loop.vertices.back() == local_it->second) {
        continue;
      }

      seed_loop.vertices.push_back(local_it->second);
      seed_loop.uvs.push_back(plan.vertices[local_it->second].uv);
    }

    if(seed_loop.vertices.size() >= 2U &&
       seed_loop.vertices.front() == seed_loop.vertices.back()) {
      seed_loop.vertices.pop_back();
      seed_loop.uvs.pop_back();
    }
    if(seed_loop.vertices.size() < 3U) {
      return core::detail::publish_error(
        base::StatusCode::unsupported,
        "Auto CFD Surface Mesher planar fallback requires each loop to resolve to at least three unique boundary samples."
      );
    }

    populate_seed_loop_edge_owners(plan.vertices, &boundary_edge_owners, seed_loop);
    add_constrained_loop_edges(
      seed_loop,
      plan.constrained_edges,
      &plan.constrained_edge_owners
    );
    seed_loops.push_back(std::move(seed_loop));
  }

  if(face_preprocess.boundary.primary_outer_loop_index >= seed_loops.size()) {
    return core::detail::publish_error(
      base::StatusCode::internal_error,
      "Auto CFD Surface Mesher could not resolve the primary outer loop from planar fallback preprocessing state."
    );
  }

  auto outer_loop = seed_loops[face_preprocess.boundary.primary_outer_loop_index];
  std::vector<SeedLoopData> inner_loops;
  inner_loops.reserve(seed_loops.size());
  for(std::size_t loop_index = 0U; loop_index < seed_loops.size(); ++loop_index) {
    if(loop_index == face_preprocess.boundary.primary_outer_loop_index) {
      continue;
    }
    if(seed_loops[loop_index].kind != geo::FaceBoundaryLoopKind::inner) {
      return core::detail::publish_error(
        base::StatusCode::unsupported,
        "Auto CFD Surface Mesher planar fallback currently supports trimmed faces only through classified outer and inner loops."
      );
    }
    inner_loops.push_back(seed_loops[loop_index]);
  }

  enforce_outer_inner_loop_winding(outer_loop, inner_loops);

  std::vector<std::array<std::uint32_t, 3>> triangles;
  if(!build_boundary_cdt_with_holes(plan.vertices, outer_loop, inner_loops, triangles)) {
    return core::detail::publish_error(
      base::StatusCode::unsupported,
      "Auto CFD Surface Mesher planar fallback CDT seeding failed."
    );
  }

  plan.triangles.reserve(triangles.size());
  for(auto triangle : triangles) {
    const double area = signed_triangle_area_twice(
      plan.vertices[triangle[0]].uv,
      plan.vertices[triangle[1]].uv,
      plan.vertices[triangle[2]].uv
    );
    if(area < 0.0) {
      std::swap(triangle[1], triangle[2]);
    }
    plan.triangles.push_back({triangle, true, false});
  }

  if(plan.triangles.empty()) {
    return core::detail::publish_error(
      base::StatusCode::unsupported,
      "Auto CFD Surface Mesher planar fallback seeded only degenerate triangles."
    );
  }

  return core::detail::clear_error_state();
}

void initialize_face_mesh_plan(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  const AutoCfdSurfaceSizingFieldState &field,
  FaceMeshPlan &plan
) noexcept
{
  plan = {};
  plan.face = face_preprocess.face;
  plan.target_size =
    effective_auto_cfd_surface_owner_target_size(field, face_preprocess.face);
  plan.reference_uv = {
    face_preprocess.reference_metric.u,
    face_preprocess.reference_metric.v,
  };
  plan.reference_normal_defined =
    resolve_face_reference_normal(face_view, face_preprocess, plan.reference_normal);
}

[[nodiscard]] base::StatusCode build_face_seed_mesh(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  const AutoCfdSurfaceSizingFieldState &field,
  FaceMeshPlan &plan
)
{
  initialize_face_mesh_plan(face_view, face_preprocess, field, plan);

  if(face_preprocess.boundary.has_seams) {
    
    return build_periodic_face_seed_mesh(face_view, face_preprocess, plan);
  }

  if(face_preprocess.boundary.inner_loop_count == 0U &&
     supports_uv_loop_containment(face_preprocess) &&
     face_has_nontrivial_curvature(face_view, face_preprocess)) {
    return build_nonplanar_single_loop_patch_seed_mesh(
      face_view,
      face_preprocess,
      plan
    );
  }

  if(face_preprocess.disposition == AutoCfdSurfaceFacePreprocessDisposition::fallback_only) {
    return build_planar_fallback_seed_mesh(face_preprocess, plan);
  }

  
  return build_planar_seed_mesh(face_preprocess, plan);
}

void build_triangle_neighbors(
  const std::vector<LocalTriangle> &triangles,
  std::vector<std::array<std::uint32_t, 3>> &neighbors
)
{
  neighbors.assign(triangles.size(), {invalid_index, invalid_index, invalid_index});

  std::unordered_map<LocalEdgeKey, std::pair<std::uint32_t, std::uint8_t>, LocalEdgeKeyHash>
    edge_to_triangle;
  edge_to_triangle.reserve(triangles.size() * 3U);

  for(std::uint32_t triangle_index = 0U;
      triangle_index < triangles.size();
      ++triangle_index) {
    const auto &triangle = triangles[triangle_index];
    if(!triangle.valid) {
      continue;
    }

    for(std::uint8_t edge_slot = 0U; edge_slot < 3U; ++edge_slot) {
      const auto first = triangle.vertices[edge_slot];
      const auto second = triangle.vertices[(edge_slot + 1U) % 3U];
      const auto key = canonical_edge(first, second);
      auto [it, inserted] =
        edge_to_triangle.emplace(key, std::pair<std::uint32_t, std::uint8_t> {triangle_index, edge_slot});
      if(inserted) {
        continue;
      }

      // Preserve the first manifold pairing we observed. If a third triangle
      // hits the same edge, keep the earlier adjacency intact instead of
      // silently overwriting it with broken non-manifold connectivity.
      if(neighbors[triangle_index][edge_slot] != invalid_index ||
         neighbors[it->second.first][it->second.second] != invalid_index) {
        continue;
      }

      neighbors[triangle_index][edge_slot] = it->second.first;
      neighbors[it->second.first][it->second.second] = triangle_index;
    }
  }
}

void build_active_front_edges(
  const FaceMeshPlan &plan,
  const std::vector<std::array<std::uint32_t, 3>> &neighbors,
  std::unordered_set<LocalEdgeKey, LocalEdgeKeyHash> &front_edges
)
{
  front_edges.clear();
  front_edges.reserve(plan.triangles.size() * 2U);

  for(std::uint32_t triangle_index = 0U;
      triangle_index < plan.triangles.size();
      ++triangle_index) {
    const auto &triangle = plan.triangles[triangle_index];
    if(!triangle.valid) {
      continue;
    }

    for(std::uint8_t edge_slot = 0U; edge_slot < 3U; ++edge_slot) {
      const auto neighbor_index = neighbors[triangle_index][edge_slot];
      bool active_edge = false;
      if(neighbor_index == invalid_index ||
         neighbor_index >= plan.triangles.size()) {
        active_edge = true;
      }
      else {
        const auto &neighbor = plan.triangles[neighbor_index];
        const double neighbor_radius = effective_triangle_front_radius(neighbor);
        active_edge =
          !neighbor.valid ||
          (neighbor_radius < kFrontRadiusLimit && neighbor_radius > 0.0);
      }

      if(!active_edge) {
        continue;
      }

      front_edges.insert(canonical_edge(
        triangle.vertices[edge_slot],
        triangle.vertices[(edge_slot + 1U) % 3U]
      ));
    }
  }
}

// ── Bowyer-Watson frontal Delaunay helpers ──────────────────────────────────

[[nodiscard]] double triangle_circumradius_3d(
  const geo::Point3 &a,
  const geo::Point3 &b,
  const geo::Point3 &c
) noexcept
{
  const auto ab = point_subtract(b, a);
  const auto ac = point_subtract(c, a);
  const auto bc = point_subtract(c, b);
  const double ab_length = std::sqrt(dot_product(ab, ab));
  const double ac_length = std::sqrt(dot_product(ac, ac));
  const double bc_length = std::sqrt(dot_product(bc, bc));
  const auto cross = cross_product(ab, ac);
  const double doubled_area = std::sqrt(dot_product(cross, cross));
  if(!std::isfinite(ab_length) ||
     !std::isfinite(ac_length) ||
     !std::isfinite(bc_length) ||
     !std::isfinite(doubled_area) ||
     ab_length <= kMetricTolerance ||
     ac_length <= kMetricTolerance ||
     bc_length <= kMetricTolerance ||
     doubled_area <= kMetricTolerance) {
    return 0.0;
  }

  return (ab_length * ac_length * bc_length) / (2.0 * doubled_area);
}

[[nodiscard]] double triangle_front_radius(
  const FaceMeshPlan &plan,
  const std::array<std::uint32_t, 3> &vertices,
  double *circumradius_3d_out = nullptr,
  double *lc_out = nullptr,
  double *lc_bgm_out = nullptr
) noexcept
{
  if(vertices[0] >= plan.vertices.size() ||
     vertices[1] >= plan.vertices.size() ||
     vertices[2] >= plan.vertices.size()) {
    return 0.0;
  }

  const auto &va = plan.vertices[vertices[0]];
  const auto &vb = plan.vertices[vertices[1]];
  const auto &vc = plan.vertices[vertices[2]];
  const double circumradius_3d =
    triangle_circumradius_3d(va.position, vb.position, vc.position);
  const double lc =
    (va.local_size + vb.local_size + vc.local_size) / 3.0;
  const double lc_bgm =
    (va.background_size + vb.background_size + vc.background_size) / 3.0;
  const double ll = std::max(std::min(lc, lc_bgm), kMetricTolerance);

  if(circumradius_3d_out != nullptr) {
    *circumradius_3d_out = circumradius_3d;
  }
  if(lc_out != nullptr) {
    *lc_out = lc;
  }
  if(lc_bgm_out != nullptr) {
    *lc_bgm_out = lc_bgm;
  }

  if(!std::isfinite(circumradius_3d) ||
     circumradius_3d <= 0.0 ||
     !std::isfinite(ll) ||
     ll <= kMetricTolerance) {
    return 0.0;
  }

  return circumradius_3d / ll;
}


[[nodiscard]] double triangle_metric_circumradius(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceSizingFieldState &field,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  const FaceMeshPlan &plan,
  const std::array<std::uint32_t, 3> &vertices
) noexcept
{
  if(vertices[0] >= plan.vertices.size() ||
     vertices[1] >= plan.vertices.size() ||
     vertices[2] >= plan.vertices.size()) {
    return 0.0;
  }

  const auto &a = plan.vertices[vertices[0]].uv;
  const auto &b = plan.vertices[vertices[1]].uv;
  const auto &c = plan.vertices[vertices[2]].uv;
  const std::array<double, 2> centroid {
    (a[0] + b[0] + c[0]) / 3.0,
    (a[1] + b[1] + c[1]) / 3.0,
  };

  AutoCfdSurfaceFaceMetricTensor metric;
  static_cast<void>(sample_metric_or_reference(
    face_view, field, face_preprocess, &plan,
    centroid[0], centroid[1], metric
  ));

  const auto metric_tensor = select_surface_metric_tensor(metric);
  if(!std::isfinite(metric_tensor[0]) ||
     !std::isfinite(metric_tensor[1]) ||
     !std::isfinite(metric_tensor[2]) ||
     metric_tensor[0] <= kMetricTolerance ||
     metric_tensor[2] <= kMetricTolerance) {
    return 0.0;
  }
  const double metric_det =
    metric_tensor[0] * metric_tensor[2] - metric_tensor[1] * metric_tensor[1];
  if(!std::isfinite(metric_det) || metric_det <= kMetricTolerance) {
    return 0.0;
  }

  std::array<double, 2> center_uv;
  double r2 = 0.0;
  if(!circumcenter_metric(a, b, c, metric_tensor, center_uv, r2)) {
    return 0.0;
  }
  const double radius = std::sqrt(r2);

  // Diagnostic: log large circumradius triangles to understand over-refinement
  static thread_local std::size_t s_diag_count = 0U;
  if(radius > 1.5 && s_diag_count < 30U) {
    ++s_diag_count;
    const auto &va = plan.vertices[vertices[0]];
    const auto &vb = plan.vertices[vertices[1]];
    const auto &vc = plan.vertices[vertices[2]];
    // Compute physical edge lengths
    const double edge_ab = std::sqrt(
      (va.position[0] - vb.position[0]) * (va.position[0] - vb.position[0]) +
      (va.position[1] - vb.position[1]) * (va.position[1] - vb.position[1]) +
      (va.position[2] - vb.position[2]) * (va.position[2] - vb.position[2])
    );
    const double edge_bc = std::sqrt(
      (vb.position[0] - vc.position[0]) * (vb.position[0] - vc.position[0]) +
      (vb.position[1] - vc.position[1]) * (vb.position[1] - vc.position[1]) +
      (vb.position[2] - vc.position[2]) * (vb.position[2] - vc.position[2])
    );
    const double edge_ca = std::sqrt(
      (vc.position[0] - va.position[0]) * (vc.position[0] - va.position[0]) +
      (vc.position[1] - va.position[1]) * (vc.position[1] - va.position[1]) +
      (vc.position[2] - va.position[2]) * (vc.position[2] - va.position[2])
    );
  }

  return radius;
}

void refresh_triangle_front_radius(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceSizingFieldState &field,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  FaceMeshPlan &plan,
  std::uint32_t triangle_index
)
{
  if(triangle_index >= plan.triangles.size()) {
    return;
  }

  auto &triangle = plan.triangles[triangle_index];
  if(!triangle.valid) {
    triangle.front_radius = 0.0;
    return;
  }

  double circumradius_3d = 0.0;
  triangle.front_radius = triangle_front_radius(
    plan,
    triangle.vertices,
    &circumradius_3d
  );

  
  if((!std::isfinite(triangle.front_radius) ||
      triangle.front_radius <= kMetricTolerance ||
      circumradius_3d <= kMetricTolerance)) {
    const double metric_radius = triangle_metric_circumradius(
      face_view, field, face_preprocess, plan, triangle.vertices
    );
    if(std::isfinite(metric_radius) && metric_radius > kMetricTolerance) {
      triangle.front_radius = metric_radius;
    }
  }
}

void refresh_all_triangle_front_radii(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceSizingFieldState &field,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  FaceMeshPlan &plan
)
{
  for(std::uint32_t triangle_index = 0U;
      triangle_index < plan.triangles.size();
      ++triangle_index) {
    auto &triangle = plan.triangles[triangle_index];
    triangle.front_radius_override = std::numeric_limits<double>::quiet_NaN();
    refresh_triangle_front_radius(
      face_view,
      field,
      face_preprocess,
      plan,
      triangle_index
    );
  }
}

void refresh_all_vertex_front_sizes(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceSizingFieldState &field,
  FaceMeshPlan &plan
)
{
  for(auto &vertex : plan.vertices) {
    const double queried_size =
      query_auto_cfd_surface_sizing_field(
        field,
        face_view,
        vertex.uv[0],
        vertex.uv[1],
        vertex.position
      );
    vertex.local_size = queried_size;
    vertex.background_size = queried_size;
  }
}

// Find the triangle containing point p using edge-walk from a hint triangle.
// Returns invalid_index if not found.
[[nodiscard]] std::uint32_t find_containing_triangle(
  const FaceMeshPlan &plan,
  const std::vector<std::array<std::uint32_t, 3>> &neighbors,
  std::uint32_t hint_triangle,
  const std::array<double, 2> &p
) noexcept
{
  if(hint_triangle >= plan.triangles.size() ||
     !plan.triangles[hint_triangle].valid) {
    hint_triangle = invalid_index;
    // Find any valid triangle as starting point
    for(std::uint32_t i = 0U; i < plan.triangles.size(); ++i) {
      if(plan.triangles[i].valid) {
        hint_triangle = i;
        break;
      }
    }
    if(hint_triangle == invalid_index) {
      return invalid_index;
    }
  }

  auto current = hint_triangle;
  const std::size_t max_walk = plan.triangles.size() * 2U;

  for(std::size_t step = 0U; step < max_walk; ++step) {
    const auto &tri = plan.triangles[current];
    if(!tri.valid) {
      return invalid_index;
    }

    const auto &va = plan.vertices[tri.vertices[0]].uv;
    const auto &vb = plan.vertices[tri.vertices[1]].uv;
    const auto &vc = plan.vertices[tri.vertices[2]].uv;

    if(point_in_triangle(p, va, vb, vc)) {
      return current;
    }

    // Walk toward p: find the edge whose opposite side contains p
    // Use the centroid-to-p direction to pick the crossing edge
    const std::array<double, 2> centroid {
      (va[0] + vb[0] + vc[0]) / 3.0,
      (va[1] + vb[1] + vc[1]) / 3.0,
    };

    bool moved = false;
    for(std::uint8_t edge = 0U; edge < 3U; ++edge) {
      const auto &e0 = plan.vertices[tri.vertices[edge]].uv;
      const auto &e1 = plan.vertices[tri.vertices[(edge + 1U) % 3U]].uv;

      // Check if the segment centroid->p crosses edge (e0, e1)
      const double d1 = signed_triangle_area_twice(e0, e1, centroid);
      const double d2 = signed_triangle_area_twice(e0, e1, p);
      if(d1 * d2 < 0.0) {
        // Crosses this edge
        const auto neigh = neighbors[current][edge];
        if(neigh != invalid_index && plan.triangles[neigh].valid) {
          current = neigh;
          moved = true;
          break;
        }
      }
    }

    if(!moved) {
      // Fallback: brute force
      break;
    }
  }

  // Brute force fallback
  for(std::uint32_t i = 0U; i < plan.triangles.size(); ++i) {
    if(!plan.triangles[i].valid) {
      continue;
    }
    const auto &va = plan.vertices[plan.triangles[i].vertices[0]].uv;
    const auto &vb = plan.vertices[plan.triangles[i].vertices[1]].uv;
    const auto &vc = plan.vertices[plan.triangles[i].vertices[2]].uv;
    if(point_in_triangle(p, va, vb, vc)) {
      return i;
    }
  }
  return invalid_index;
}

[[nodiscard]] std::array<double, 2> triangle_barycenter_uv(
  const FaceMeshPlan &plan,
  std::uint32_t triangle_index
) noexcept
{
  const auto &triangle = plan.triangles[triangle_index];
  return {
    (plan.vertices[triangle.vertices[0]].uv[0] +
     plan.vertices[triangle.vertices[1]].uv[0] +
     plan.vertices[triangle.vertices[2]].uv[0]) / 3.0,
    (plan.vertices[triangle.vertices[0]].uv[1] +
     plan.vertices[triangle.vertices[1]].uv[1] +
     plan.vertices[triangle.vertices[2]].uv[1]) / 3.0,
  };
}

[[nodiscard]] bool triangle_barycentric_coordinates(
  const std::array<double, 2> &point,
  const std::array<double, 2> &a,
  const std::array<double, 2> &b,
  const std::array<double, 2> &c,
  std::array<double, 3> &barycentric,
  double tolerance
) noexcept
{
  const double determinant =
    (b[0] - a[0]) * (c[1] - a[1]) -
    (b[1] - a[1]) * (c[0] - a[0]);
  if(std::abs(determinant) <= tolerance) {
    barycentric = {0.0, 0.0, 0.0};
    return false;
  }

  const double inverse_determinant = 1.0 / determinant;
  const double u =
    ((point[0] - a[0]) * (c[1] - a[1]) -
     (point[1] - a[1]) * (c[0] - a[0])) * inverse_determinant;
  const double v =
    ((b[0] - a[0]) * (point[1] - a[1]) -
     (b[1] - a[1]) * (point[0] - a[0])) * inverse_determinant;
  const double w = 1.0 - u - v;

  barycentric = {w, u, v};
  return barycentric[0] >= -tolerance &&
         barycentric[1] >= -tolerance &&
         barycentric[2] >= -tolerance &&
         barycentric[0] <= 1.0 + tolerance &&
         barycentric[1] <= 1.0 + tolerance &&
         barycentric[2] <= 1.0 + tolerance;
}

[[nodiscard]] bool triangle_contains_point_barycentric(
  const FaceMeshPlan &plan,
  std::uint32_t triangle_index,
  const std::array<double, 2> &point,
  std::array<double, 3> &barycentric,
  double tolerance
) noexcept
{
  if(triangle_index >= plan.triangles.size() ||
     !plan.triangles[triangle_index].valid) {
    barycentric = {0.0, 0.0, 0.0};
    return false;
  }

  const auto &triangle = plan.triangles[triangle_index];
  return triangle_barycentric_coordinates(
    point,
    plan.vertices[triangle.vertices[0]].uv,
    plan.vertices[triangle.vertices[1]].uv,
    plan.vertices[triangle.vertices[2]].uv,
    barycentric,
    tolerance
  );
}

[[nodiscard]] std::uint32_t search_for_triangle(
  const FaceMeshPlan &plan,
  const std::vector<std::array<std::uint32_t, 3>> &neighbors,
  const FrontTriangleSet &all_triangles,
  std::uint32_t start_triangle,
  const std::array<double, 2> &point,
  std::array<double, 3> &barycentric,
  bool force
) noexcept
{
  constexpr double kSearchTolerance = 1.0e-8;

  if(start_triangle >= plan.triangles.size() ||
     !plan.triangles[start_triangle].valid) {
    if(!all_triangles.empty()) {
      start_triangle = *all_triangles.begin();
    }
    else {
      start_triangle = invalid_index;
    }
  }

  if(start_triangle == invalid_index) {
    barycentric = {0.0, 0.0, 0.0};
    return invalid_index;
  }

  if(triangle_contains_point_barycentric(
       plan,
       start_triangle,
       point,
       barycentric,
       kSearchTolerance
     )) {
    return start_triangle;
  }

  std::uint32_t current_triangle = start_triangle;
  std::size_t iteration = 0U;
  const std::size_t max_iteration =
    std::max<std::size_t>(all_triangles.size(), plan.triangles.size());

  while(current_triangle != invalid_index &&
        current_triangle < plan.triangles.size() &&
        plan.triangles[current_triangle].valid) {
    const auto q2 = triangle_barycenter_uv(plan, current_triangle);
    bool moved = false;

    for(std::uint8_t edge_slot = 0U; edge_slot < 3U; ++edge_slot) {
      const auto &triangle = plan.triangles[current_triangle];
      const auto &p0 = plan.vertices[triangle.vertices[edge_slot]].uv;
      const auto &p1 = plan.vertices[triangle.vertices[(edge_slot + 1U) % 3U]].uv;
      if(!uv_segments_intersect(p0, p1, point, q2)) {
        continue;
      }

      current_triangle = neighbors[current_triangle][edge_slot];
      moved = true;
      break;
    }

    if(!moved ||
       current_triangle == invalid_index ||
       current_triangle >= plan.triangles.size() ||
       !plan.triangles[current_triangle].valid) {
      break;
    }

    if(triangle_contains_point_barycentric(
         plan,
         current_triangle,
         point,
         barycentric,
         kSearchTolerance
       )) {
      return current_triangle;
    }

    if(iteration++ > max_iteration) {
      break;
    }
  }

  if(!force) {
    barycentric = {0.0, 0.0, 0.0};
    return invalid_index;
  }

  for(const auto triangle_index : all_triangles) {
    if(triangle_contains_point_barycentric(
         plan,
         triangle_index,
         point,
         barycentric,
         kSearchTolerance
       )) {
      return triangle_index;
    }
  }

  barycentric = {0.0, 0.0, 0.0};
  return invalid_index;
}

// Edge of the cavity shell: the edge of a cavity triangle that borders
// a non-cavity triangle (or the boundary).
struct CavityShellEdge final {
  std::uint32_t v0 = invalid_index;
  std::uint32_t v1 = invalid_index;
  // The neighbor triangle on the other side (invalid_index if boundary)
  std::uint32_t neighbor_triangle = invalid_index;
};

// Iteratively find the Bowyer-Watson cavity: all triangles whose metric-space
// circumcircle contains the new point. Respects constrained edges.
// Uses an explicit stack to avoid stack overflow on large cavities.
void find_cavity_aniso(
  const FaceMeshPlan &plan,
  const std::vector<std::array<std::uint32_t, 3>> &neighbors,
  const std::array<double, 3> &metric_tensor,
  const std::array<double, 2> &new_point,
  std::uint32_t start_triangle,
  std::vector<bool> &deleted,
  std::vector<std::uint32_t> &cavity,
  std::vector<CavityShellEdge> &shell
)
{
  // Use an explicit stack instead of recursion to avoid stack overflow
  std::vector<std::uint32_t> stack;
  stack.reserve(32U);
  stack.push_back(start_triangle);
  deleted[start_triangle] = true;

  // Safety limit: cavity should never exceed total triangle count
  const std::size_t max_cavity_size = plan.triangles.size();

  while(!stack.empty()) {
    const auto triangle_index = stack.back();
    stack.pop_back();
    cavity.push_back(triangle_index);

    if(cavity.size() > max_cavity_size) {
      // Something went wrong — bail out
      break;
    }

    const auto &tri = plan.triangles[triangle_index];
    for(std::uint8_t edge = 0U; edge < 3U; ++edge) {
      const auto v0 = tri.vertices[edge];
      const auto v1 = tri.vertices[(edge + 1U) % 3U];

      // Check if this edge is constrained
      if(plan.constrained_edges.find(canonical_edge(v0, v1)) !=
         plan.constrained_edges.end()) {
        shell.push_back({v0, v1, neighbors[triangle_index][edge]});
        continue;
      }

      const auto neigh = neighbors[triangle_index][edge];
      if(neigh == invalid_index) {
        shell.push_back({v0, v1, invalid_index});
        continue;
      }

      if(deleted[neigh]) {
        continue;
      }

      if(!plan.triangles[neigh].valid) {
        shell.push_back({v0, v1, neigh});
        continue;
      }

      const auto &ntri = plan.triangles[neigh];
      const auto &na = plan.vertices[ntri.vertices[0]].uv;
      const auto &nb = plan.vertices[ntri.vertices[1]].uv;
      const auto &nc = plan.vertices[ntri.vertices[2]].uv;

      if(point_in_circumcircle_aniso(na, nb, nc, new_point, metric_tensor)) {
        deleted[neigh] = true;
        stack.push_back(neigh);
      }
      else {
        shell.push_back({v0, v1, neigh});
      }
    }
  }
}

// Insert a vertex into the Bowyer-Watson cavity by re-triangulating the shell.
// Returns true on success. On failure, restores the deleted flags.
[[nodiscard]] bool insert_vertex_into_cavity(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceSizingFieldState &field,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  FaceMeshPlan &plan,
  std::vector<bool> &deleted,
  const std::vector<std::uint32_t> &cavity,
  const std::vector<CavityShellEdge> &shell,
  std::uint32_t new_vertex_index
)
{
  static_cast<void>(face_view);
  static_cast<void>(field);
  static_cast<void>(face_preprocess);

  const auto restore_cavity =
    [&deleted, &cavity]() noexcept {
      for(const auto cavity_triangle : cavity) {
        deleted[cavity_triangle] = false;
      }
    };

  
  if(cavity.size() < 2U || shell.size() < 3U) {
    restore_cavity();
    return false;
  }

  if(shell.size() != cavity.size() + 2U) {
   
    restore_cavity();
    return false;
  }

  if(new_vertex_index >= plan.vertices.size()) {
    restore_cavity();
    return false;
  }

  double old_area = 0.0;
  for(auto ci : cavity) {
    const auto &tri = plan.triangles[ci];
    old_area += std::abs(signed_triangle_area_twice(
      plan.vertices[tri.vertices[0]].uv,
      plan.vertices[tri.vertices[1]].uv,
      plan.vertices[tri.vertices[2]].uv
    ));
  }

  // Build new triangles from shell edges to new vertex
  double new_area = 0.0;
  bool one_point_is_too_close = false;
  std::vector<std::array<std::uint32_t, 3>> new_triangles;
  new_triangles.reserve(shell.size());

  for(const auto &edge : shell) {
    std::array<std::uint32_t, 3> tri_verts {edge.v0, edge.v1, new_vertex_index};
    if(!orient_triangle_ccw(plan.vertices, tri_verts)) {
      restore_cavity();
      return false;
    }
    const double signed_area = signed_triangle_area_twice(
      plan.vertices[tri_verts[0]].uv,
      plan.vertices[tri_verts[1]].uv,
      plan.vertices[tri_verts[2]].uv
    );
    if(!std::isfinite(signed_area) || std::abs(signed_area) < 1.0e-25) {
      restore_cavity();
      return false;
    }

    new_area += std::abs(signed_area);
    new_triangles.push_back(tri_verts);

    
    {
      const auto &p0 = plan.vertices[edge.v0].position;
      const auto &p1 = plan.vertices[edge.v1].position;
      const auto &pv = plan.vertices[new_vertex_index].position;

      const double dx0 = pv[0] - p0[0], dy0 = pv[1] - p0[1], dz0 = pv[2] - p0[2];
      const double dx1 = pv[0] - p1[0], dy1 = pv[1] - p1[1], dz1 = pv[2] - p1[2];
      const double d1 = std::sqrt(dx0*dx0 + dy0*dy0 + dz0*dz0);
      const double d2 = std::sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);

      const double ex = p1[0] - p0[0], ey = p1[1] - p0[1], ez = p1[2] - p0[2];
      const double d3 = std::sqrt(ex*ex + ey*ey + ez*ez);

      
      const double lc =
        (plan.vertices[edge.v0].local_size +
         plan.vertices[edge.v1].local_size +
         plan.vertices[new_vertex_index].local_size) / 3.0;
      const double lc_bgm =
        (plan.vertices[edge.v0].background_size +
         plan.vertices[edge.v1].background_size +
         plan.vertices[new_vertex_index].background_size) / 3.0;
      const double LL = std::max(
        std::min(lc, lc_bgm),
        kMetricTolerance
      );

      
      double d4 = 1.0e22;
      if(plan.vertices[edge.v0].boundary && plan.vertices[edge.v1].boundary) {
        // cross product of edge vector and v0-to-new vector
        const double cx = ey * dz0 - ez * dy0;
        const double cy = ez * dx0 - ex * dz0;
        const double cz = ex * dy0 - ey * dx0;
        if(d3 > 1.0e-30) {
          d4 = std::sqrt(cx*cx + cy*cy + cz*cz) / d3;
        }
      }

      // Angle check at new vertex: cos(angle) = (d1² + d2² - d3²) / (2·d1·d2)
      double cosv = 0.0;
      if(d1 > 1.0e-30 && d2 > 1.0e-30) {
        cosv = (d1*d1 + d2*d2 - d3*d3) / (2.0 * d1 * d2);
      }

      if(d1 < LL * 0.5 || d2 < LL * 0.5 || d4 < LL * 0.4 || cosv < -0.9999) {
        one_point_is_too_close = true;
      }
    }
  }

  const double EPS = 1.0e-12;
  if(old_area > kUvAreaTolerance &&
     std::abs(old_area - new_area) > EPS * old_area) {
    restore_cavity();
    return false;
  }
  if(one_point_is_too_close) {
    restore_cavity();
    return false;
  }

  // Commit: invalidate cavity triangles, add new triangles
  for(auto ci : cavity) {
    plan.triangles[ci].valid = false;
    plan.triangles[ci].accepted = false;
  }

  for(const auto &tri_verts : new_triangles) {
    plan.triangles.push_back({tri_verts, true, false});
  }

  return true;
}


[[nodiscard]] bool solve_linear_3x3(
  const double A[3][3],
  const double b[3],
  double x[3]
) noexcept
{
  const double det =
    A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
    A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
    A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
  if(!std::isfinite(det) || std::abs(det) < 1.0e-30) {
    return false;
  }
  const double inv_det = 1.0 / det;
  x[0] = inv_det * (
    b[0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
    A[0][1] * (b[1] * A[2][2] - A[1][2] * b[2]) +
    A[0][2] * (b[1] * A[2][1] - A[1][1] * b[2]));
  x[1] = inv_det * (
    A[0][0] * (b[1] * A[2][2] - A[1][2] * b[2]) -
    b[0] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
    A[0][2] * (A[1][0] * b[2] - b[1] * A[2][0]));
  x[2] = inv_det * (
    A[0][0] * (A[1][1] * b[2] - b[1] * A[2][1]) -
    A[0][1] * (A[1][0] * b[2] - b[1] * A[2][0]) +
    b[0] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]));
  return std::isfinite(x[0]) && std::isfinite(x[1]) && std::isfinite(x[2]);
}


[[nodiscard]] bool refine_frontal_point_3d(
  const geo::FaceView &face_view,
  const std::array<double, 2> &initial_uv,
  const geo::Point3 &pos_a,
  const geo::Point3 &pos_b,
  const geo::Point3 &pos_opposite,
  double ideal_height_3d,
  std::array<double, 2> &refined_uv,
  geo::Point3 &refined_position
) noexcept
{
  // 3D midpoint of the active edge.
  const geo::Point3 middle {
    0.5 * (pos_a[0] + pos_b[0]),
    0.5 * (pos_a[1] + pos_b[1]),
    0.5 * (pos_a[2] + pos_b[2]),
  };

  
  const auto v1v2 = point_subtract(pos_b, pos_a);
  const auto tmp = point_subtract(pos_opposite, middle);
  auto n1_normal = cross_product(v1v2, tmp);
  const double n1_len = std::sqrt(dot_product(n1_normal, n1_normal));
  if(!std::isfinite(n1_len) || n1_len < 1.0e-12) {
    return false;
  }
  n1_normal = {n1_normal[0] / n1_len, n1_normal[1] / n1_len, n1_normal[2] / n1_len};

  auto n2_inplane = cross_product(n1_normal, v1v2);
  const double n2_len = std::sqrt(dot_product(n2_inplane, n2_inplane));
  if(!std::isfinite(n2_len) || n2_len < 1.0e-12) {
    return false;
  }
  n2_inplane = {n2_inplane[0] / n2_len, n2_inplane[1] / n2_len, n2_inplane[2] / n2_len};

  const double d = ideal_height_3d;
  if(!std::isfinite(d) || d <= 0.0) {
    return false;
  }

  
  const auto circle_point = [&](double t, double out[3]) noexcept {
    const double ct = std::cos(t);
    const double st = std::sin(t);
    out[0] = middle[0] + d * (n2_inplane[0] * ct + n1_normal[0] * st);
    out[1] = middle[1] + d * (n2_inplane[1] * ct + n1_normal[1] * st);
    out[2] = middle[2] + d * (n2_inplane[2] * ct + n1_normal[2] * st);
  };

  // residual r(u, v, t) = S(u, v) − C(t)
  const auto evaluate_residual =
    [&](double u, double v, double t, double res[3]) noexcept -> bool {
      geo::FaceSample sample;
      if(geo::sample_face(face_view, u, v, sample) != base::StatusCode::ok) {
        return false;
      }
      double cp[3];
      circle_point(t, cp);
      res[0] = sample.position[0] - cp[0];
      res[1] = sample.position[1] - cp[1];
      res[2] = sample.position[2] - cp[2];
      return true;
    };

  // Initial guess: (u_2D, v_2D, 0) — t=0 places the circle point at
  // middle + d·n2_inplane (in the tangent plane on the opposite-vertex side),
  // which is what the 2D guess from optimalPointFrontal already approximates.
  double uvt[3] = {initial_uv[0], initial_uv[1], 0.0};

  // Pre-Newton early exit: if the residual at the initial guess is
  // already below epsilon, return immediately without iterating.
  const double kEpsilon = d * 1.0e-8;
  double f[3];
  if(!evaluate_residual(uvt[0], uvt[1], uvt[2], f)) {
    return false;
  }
  double f_norm = std::sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
  if(f_norm < kEpsilon) {
    geo::FaceSample sample;
    if(geo::sample_face(face_view, uvt[0], uvt[1], sample) != base::StatusCode::ok) {
      return false;
    }
    refined_uv = {uvt[0], uvt[1]};
    refined_position = sample.position;
    return true;
  }

  
  constexpr int kMaxIterations = 100;
  constexpr double kFiniteDifferenceStep = 1.0e-4;
  constexpr double kStepTolerance = 1.0e-6;

  for(int iter = 0; iter < kMaxIterations; ++iter) {
    if(std::abs(uvt[0]) > 1.0e6 ||
       std::abs(uvt[1]) > 1.0e6 ||
       std::abs(uvt[2]) > 1.0e6) {
      return false;
    }

    
    double J[3][3];
    for(int j = 0; j < 3; ++j) {
      double h = kFiniteDifferenceStep * std::abs(uvt[j]);
      if(h == 0.0) {
        h = kFiniteDifferenceStep;
      }
      uvt[j] += h;
      double feps[3];
      if(!evaluate_residual(uvt[0], uvt[1], uvt[2], feps)) {
        uvt[j] -= h;
        return false;
      }
      uvt[j] -= h;
      for(int i = 0; i < 3; ++i) {
        J[i][j] = (feps[i] - f[i]) / h;
      }
    }

    // Solve J · dx = f → dx
    double dx[3];
    if(!solve_linear_3x3(J, f, dx)) {
      return false;
    }

    uvt[0] -= dx[0];
    uvt[1] -= dx[1];
    uvt[2] -= dx[2];

    const double dx_norm = std::sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
    if(dx_norm < kStepTolerance) {
      geo::FaceSample sample;
      if(geo::sample_face(face_view, uvt[0], uvt[1], sample) != base::StatusCode::ok) {
        return false;
      }
      refined_uv = {uvt[0], uvt[1]};
      refined_position = sample.position;
      return true;
    }

    // Re-evaluate residual at the updated x for the next iteration.
    if(!evaluate_residual(uvt[0], uvt[1], uvt[2], f)) {
      return false;
    }
  }

  return false;
}


[[nodiscard]] bool compute_frontal_optimal_point(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceSizingFieldState &field,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  const FaceMeshPlan &plan,
  const std::array<std::uint32_t, 3> &tri_vertices,
  std::uint8_t active_edge,
  std::array<double, 2> &candidate_uv,
  geo::Point3 &candidate_position,
  std::array<double, 3> *metric_tensor_out = nullptr
)
{
  // active_edge is the edge slot from find_active_edge:
  //   slot e = edge from vertices[e] to vertices[(e+1)%3]
  // The "opposite" vertex (the one NOT on the active edge) is vertices[(e+2)%3].
  // a and b are the two endpoints of the active edge; c is the opposite vertex.
  const auto a_index = tri_vertices[active_edge];
  const auto b_index = tri_vertices[(active_edge + 1U) % 3U];
  if(a_index >= plan.vertices.size() || b_index >= plan.vertices.size()) {
    return false;
  }

  const auto &a = plan.vertices[a_index].uv;
  const auto &b = plan.vertices[b_index].uv;
  const std::array<double, 2> midpoint {
    0.5 * (a[0] + b[0]),
    0.5 * (a[1] + b[1]),
  };
  // Sample metric at active edge midpoint — the candidate point will be
  // placed near this location, so the metric here is more representative
  // than the triangle centroid (which includes the far opposite vertex).
  AutoCfdSurfaceFaceMetricTensor metric;
  static_cast<void>(sample_metric_or_reference(
    face_view, field, face_preprocess, &plan,
    midpoint[0], midpoint[1], metric
  ));

  
  const auto metric_tensor = select_surface_metric_tensor(metric);
  if(!std::isfinite(metric_tensor[0]) ||
     !std::isfinite(metric_tensor[1]) ||
     !std::isfinite(metric_tensor[2]) ||
     metric_tensor[0] <= kMetricTolerance ||
     metric_tensor[2] <= kMetricTolerance) {
    return false;
  }
  const double metric_det =
    metric_tensor[0] * metric_tensor[2] - metric_tensor[1] * metric_tensor[1];
  if(!std::isfinite(metric_det) || metric_det <= kMetricTolerance) {
    return false;
  }

  
  std::array<double, 2> circumcenter_uv;
  double r2 = 0.0;
  if(!circumcenter_metric(
       plan.vertices[tri_vertices[0]].uv,
       plan.vertices[tri_vertices[1]].uv,
       plan.vertices[tri_vertices[2]].uv,
       metric_tensor,
       circumcenter_uv,
       r2
     )) {
    return false;
  }

  // Direction from midpoint of active edge toward circumcenter
  auto dir = uv_subtract(circumcenter_uv, midpoint);
  const double dir_norm = uv_norm(dir);
  if(dir_norm <= kMetricTolerance) {
    return false;
  }
  dir = uv_scale(dir, 1.0 / dir_norm);

  
  const double ratio_squared = raw_metric_distance_squared(metric_tensor, dir);
  if(!std::isfinite(ratio_squared) || ratio_squared <= kMetricTolerance) {
    return false;
  }
  const double ratio = std::sqrt(ratio_squared);

  const double rhoM1 =
    0.5 * (plan.vertices[a_index].local_size + plan.vertices[b_index].local_size);
  const double rhoM2 = 0.5 * (
    plan.vertices[a_index].background_size +
    plan.vertices[b_index].background_size
  );
  const double rhoM = std::min(rhoM1, rhoM2);
  if(!std::isfinite(rhoM) || rhoM <= kMetricTolerance) {
    return false;
  }

  
  const double d = rhoM * kFrontIdealHeight;

  
  const double q = std::sqrt(raw_metric_distance_squared(
    metric_tensor,
    uv_subtract(circumcenter_uv, midpoint)
  ));

  const double L = std::min(d, q);

  
  candidate_uv = {
    midpoint[0] + L * dir[0] / ratio,
    midpoint[1] + L * dir[1] / ratio,
  };

  if(!plan.planar_fallback) {
    const auto c_index = tri_vertices[(active_edge + 2U) % 3U];
    if(c_index >= plan.vertices.size()) {
      return false;
    }

    const auto &pos_a = plan.vertices[a_index].position;
    const auto &pos_b = plan.vertices[b_index].position;
    const auto &pos_c = plan.vertices[c_index].position;

    std::array<double, 2> refined_uv {};
    geo::Point3 refined_position {};
    if(refine_frontal_point_3d(
         face_view,
         candidate_uv,
         pos_a, pos_b, pos_c,
         d,
         refined_uv, refined_position
       )) {
      // Seam alignment
      if(uses_supported_seam_material_screen(face_preprocess)) {
        align_uv_to_reference_periodically(
          face_preprocess.uv_bounds,
          candidate_uv,
          refined_uv
        );
      }
      if(uv_point_inside_face_material(face_preprocess, refined_uv)) {
        candidate_uv = refined_uv;
        candidate_position = refined_position;
        if(metric_tensor_out != nullptr) {
          *metric_tensor_out = metric_tensor;
        }
        return true;
      }
    }

    // 3D refinement failed — fall back to UV-space result with simple projection
  }

  // Planar fallback or 3D refinement failure: project UV candidate to surface
  if(plan.planar_fallback) {
    // Planar: compute world point from UV, project back to get snapped position.
    const auto candidate_point = planar_world_point(plan, candidate_uv);
    geo::FaceProjection projection;
    if(!std::isfinite(candidate_point[0]) ||
       geo::project_point_to_face(face_view, candidate_point, projection) !=
         base::StatusCode::ok ||
       projection.distance > projection_tolerance(metric.target_size)) {
      return false;
    }
    candidate_uv = planar_local_uv(plan, projection.projected_point);
    candidate_position = projection.projected_point;
  }
  else {
    // Non-planar: sample_face already gives us the exact 3D position for the
    // UV coordinates — no need to project back.  This eliminates one expensive
    // OCC call per iteration.
    geo::FaceSample sample;
    if(geo::sample_face(face_view, candidate_uv[0], candidate_uv[1], sample) !=
       base::StatusCode::ok) {
      return false;
    }
    if(uses_supported_seam_material_screen(face_preprocess)) {
      auto recovered_uv = candidate_uv;
      align_uv_to_reference_periodically(
        face_preprocess.uv_bounds,
        candidate_uv,
        recovered_uv
      );
      candidate_uv = recovered_uv;
    }
    if(!uv_point_inside_face_material(face_preprocess, candidate_uv)) {
      return false;
    }
    candidate_position = sample.position;
  }
  if(metric_tensor_out != nullptr) {
    *metric_tensor_out = metric_tensor;
  }
  return true;
}

void update_triangle_front_override(
  FaceMeshPlan &plan,
  FrontTriangleSet &all_triangles,
  std::uint32_t triangle_index,
  double radius_override
)
{
  if(triangle_index >= plan.triangles.size() ||
     !plan.triangles[triangle_index].valid) {
    return;
  }

  all_triangles.erase(triangle_index);
  plan.triangles[triangle_index].front_radius_override = radius_override;
  all_triangles.insert(triangle_index);
}

void reconnect_local_triangle_neighbors(
  const std::vector<LocalTriangle> &triangles,
  const std::vector<CavityShellEdge> &shell,
  std::uint32_t new_triangles_begin,
  std::vector<std::array<std::uint32_t, 3>> &neighbors
)
{
  if(neighbors.size() < triangles.size()) {
    neighbors.resize(
      triangles.size(),
      {invalid_index, invalid_index, invalid_index}
    );
  }
  // Freshly-created triangles start with all neighbors = invalid; the
  // sort-and-pair below will fill them in.
  for(std::uint32_t t = new_triangles_begin; t < triangles.size(); ++t) {
    neighbors[t] = {invalid_index, invalid_index, invalid_index};
  }

  // Collect the unique set of shell-outside neighbor triangles. Multiple
  // shell edges can in principle share the same outside triangle when a
  // cavity bites deep into a neighbor ring, so we dedup with sort+unique.
  std::vector<std::uint32_t> outside_triangles;
  outside_triangles.reserve(shell.size());
  for(const auto &shell_edge : shell) {
    if(shell_edge.neighbor_triangle == invalid_index ||
       shell_edge.neighbor_triangle >= triangles.size() ||
       !triangles[shell_edge.neighbor_triangle].valid) {
      continue;
    }
    outside_triangles.push_back(shell_edge.neighbor_triangle);
  }
  std::sort(outside_triangles.begin(), outside_triangles.end());
  outside_triangles.erase(
    std::unique(outside_triangles.begin(), outside_triangles.end()),
    outside_triangles.end()
  );

  
  struct EdgeSlot final {
    LocalEdgeKey key;
    std::uint32_t triangle;
    std::uint8_t slot;
  };

  std::vector<EdgeSlot> records;
  records.reserve(
    (triangles.size() - new_triangles_begin + outside_triangles.size()) * 3U
  );

  const auto push_triangle_edges = [&](std::uint32_t t) noexcept {
    for(std::uint8_t slot = 0U; slot < 3U; ++slot) {
      const auto first = triangles[t].vertices[slot];
      const auto second = triangles[t].vertices[(slot + 1U) % 3U];
      records.push_back({canonical_edge(first, second), t, slot});
    }
  };

  for(std::uint32_t t = new_triangles_begin; t < triangles.size(); ++t) {
    if(!triangles[t].valid) {
      continue;
    }
    push_triangle_edges(t);
  }
  for(const auto t : outside_triangles) {
    push_triangle_edges(t);
  }

  std::sort(
    records.begin(),
    records.end(),
    [](const EdgeSlot &lhs, const EdgeSlot &rhs) noexcept {
      if(lhs.key.first != rhs.key.first) {
        return lhs.key.first < rhs.key.first;
      }
      return lhs.key.second < rhs.key.second;
    }
  );

  for(std::size_t i = 0U; i + 1U < records.size(); ++i) {
    auto &f1 = records[i];
    auto &f2 = records[i + 1U];
    if(f1.key == f2.key && f1.triangle != f2.triangle) {
      neighbors[f1.triangle][f1.slot] = f2.triangle;
      neighbors[f2.triangle][f2.slot] = f1.triangle;
      ++i;
    }
  }
}

[[nodiscard]] bool insert_a_point(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceSizingFieldState &field,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  FaceMeshPlan &plan,
  std::vector<std::array<std::uint32_t, 3>> &neighbors,
  FrontTriangleSet &all_triangles,
  FrontTriangleSet &active_triangles,
  std::uint32_t worst_triangle_index,
  const std::array<double, 3> &metric_tensor,
  const std::array<double, 2> &candidate_uv,
  const geo::Point3 &candidate_position,
  std::vector<bool> &deleted,
  std::vector<std::uint32_t> &cavity,
  std::vector<CavityShellEdge> &shell
)
{
  if(worst_triangle_index >= plan.triangles.size() ||
     !plan.triangles[worst_triangle_index].valid) {
    return false;
  }

  const auto &worst_triangle = plan.triangles[worst_triangle_index];

  deleted.assign(plan.triangles.size(), false);
  cavity.clear();
  shell.clear();
  cavity.reserve(16U);
  shell.reserve(24U);

  std::uint32_t point_triangle_index = invalid_index;
  std::array<double, 3> barycentric {0.0, 0.0, 0.0};

  const auto &wa = plan.vertices[worst_triangle.vertices[0]].uv;
  const auto &wb = plan.vertices[worst_triangle.vertices[1]].uv;
  const auto &wc = plan.vertices[worst_triangle.vertices[2]].uv;
  if(point_in_circumcircle_aniso(wa, wb, wc, candidate_uv, metric_tensor)) {
    find_cavity_aniso(
      plan,
      neighbors,
      metric_tensor,
      candidate_uv,
      worst_triangle_index,
      deleted,
      cavity,
      shell
    );
    for(const auto cavity_triangle : cavity) {
      if(triangle_contains_point_barycentric(
           plan,
           cavity_triangle,
           candidate_uv,
           barycentric,
           1.0e-8
         )) {
        point_triangle_index = cavity_triangle;
        break;
      }
    }
  }
  else {
    point_triangle_index = search_for_triangle(
      plan,
      neighbors,
      all_triangles,
      worst_triangle_index,
      candidate_uv,
      barycentric,
      false
    );
    if(point_triangle_index != invalid_index) {
      deleted.assign(plan.triangles.size(), false);
      cavity.clear();
      shell.clear();
      find_cavity_aniso(
        plan,
        neighbors,
        metric_tensor,
        candidate_uv,
        point_triangle_index,
        deleted,
        cavity,
        shell
      );
    }
  }

  if(point_triangle_index == invalid_index) {
    update_triangle_front_override(plan, all_triangles, worst_triangle_index, 0.0);
    return false;
  }

  if(barycentric == std::array<double, 3> {0.0, 0.0, 0.0}) {
    if(!triangle_contains_point_barycentric(
         plan,
         point_triangle_index,
         candidate_uv,
         barycentric,
         1.0e-8
       )) {
      update_triangle_front_override(plan, all_triangles, worst_triangle_index, 0.0);
      return false;
    }
  }

  const auto &point_triangle = plan.triangles[point_triangle_index];
  const double local_size =
    barycentric[0] * plan.vertices[point_triangle.vertices[0]].local_size +
    barycentric[1] * plan.vertices[point_triangle.vertices[1]].local_size +
    barycentric[2] * plan.vertices[point_triangle.vertices[2]].local_size;
  const double background_size =
    query_auto_cfd_surface_sizing_field(
      field,
      face_view,
      candidate_uv[0],
      candidate_uv[1],
      candidate_position
    );

  const auto new_index = static_cast<std::uint32_t>(plan.vertices.size());
  const auto triangle_count_before_insertion =
    static_cast<std::uint32_t>(plan.triangles.size());
  plan.vertices.push_back({
    invalid_index,
    candidate_position,
    candidate_uv,
    false,
    local_size,
    background_size,
  });

  if(!insert_vertex_into_cavity(
       face_view,
       field,
       face_preprocess,
       plan,
       deleted,
       cavity,
       shell,
       new_index
     )) {
    plan.vertices.pop_back();
    update_triangle_front_override(plan, all_triangles, worst_triangle_index, -1.0);
    return false;
  }

  for(const auto cavity_triangle : cavity) {
    active_triangles.erase(cavity_triangle);
    all_triangles.erase(cavity_triangle);
  }

  for(std::uint32_t triangle_index = triangle_count_before_insertion;
      triangle_index < plan.triangles.size();
      ++triangle_index) {
    refresh_triangle_front_radius(
      face_view,
      field,
      face_preprocess,
      plan,
      triangle_index
    );
    all_triangles.insert(triangle_index);
  }

  reconnect_local_triangle_neighbors(
    plan.triangles,
    shell,
    triangle_count_before_insertion,
    neighbors
  );

  for(std::uint32_t triangle_index = triangle_count_before_insertion;
      triangle_index < plan.triangles.size();
      ++triangle_index) {
    if(!plan.triangles[triangle_index].valid) {
      continue;
    }

    const auto active_edge =
      find_active_edge(plan, neighbors, triangle_index, nullptr);
    if(active_edge < 3U &&
       effective_triangle_front_radius(plan.triangles[triangle_index]) >
         kFrontRadiusLimit) {
      active_triangles.insert(triangle_index);
    }
  }

  for(const auto &shell_edge : shell) {
    const auto neighbor_triangle = shell_edge.neighbor_triangle;
    if(neighbor_triangle == invalid_index ||
       neighbor_triangle >= plan.triangles.size() ||
       !plan.triangles[neighbor_triangle].valid) {
      continue;
    }

    const auto active_edge =
      find_active_edge(plan, neighbors, neighbor_triangle, nullptr);
    if(active_edge < 3U &&
       effective_triangle_front_radius(plan.triangles[neighbor_triangle]) >
         kFrontRadiusLimit) {
      active_triangles.insert(neighbor_triangle);
    }
  }

  return true;
}

// Determine the active (front) edge of a triangle: an edge where the neighbor
// is either missing, invalid, or already accepted. Returns the edge slot
// (0,1,2) or 3 if no active edge.
[[nodiscard]] std::uint8_t find_active_edge(
  const FaceMeshPlan &plan,
  const std::vector<std::array<std::uint32_t, 3>> &neighbors,
  std::uint32_t triangle_index,
  const std::unordered_set<LocalEdgeKey, LocalEdgeKeyHash> *front_edges = nullptr
) noexcept
{
  const auto &tri = plan.triangles[triangle_index];
  for(std::uint8_t edge = 0U; edge < 3U; ++edge) {
    const auto neigh = neighbors[triangle_index][edge];
    bool active_edge = false;
    if(neigh == invalid_index || neigh >= plan.triangles.size()) {
      active_edge = true;
    }
    else {
      const auto &neighbor = plan.triangles[neigh];
      const double neighbor_radius = effective_triangle_front_radius(neighbor);
      active_edge =
        !neighbor.valid ||
        (neighbor_radius < kFrontRadiusLimit && neighbor_radius > 0.0);
    }

    if(!active_edge) {
      continue;
    }

    if(front_edges != nullptr &&
       front_edges->find(canonical_edge(
         tri.vertices[edge],
         tri.vertices[(edge + 1U) % 3U]
       )) == front_edges->end()) {
      continue;
    }

    return edge;
  }
  return 3U;
}

// ── End Bowyer-Watson helpers ───────────────────────────────────────────────

[[nodiscard]] FrontEdgeCandidate select_front_edge(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceSizingFieldState &field,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  const FaceMeshPlan &plan,
  const std::vector<std::array<std::uint32_t, 3>> &neighbors
)
{
  FrontEdgeCandidate selected;
  std::unordered_set<LocalEdgeKey, LocalEdgeKeyHash> front_edges;
  build_active_front_edges(plan, neighbors, front_edges);

  for(std::uint32_t triangle_index = 0U;
      triangle_index < plan.triangles.size();
      ++triangle_index) {
    const auto &triangle = plan.triangles[triangle_index];
    if(!triangle.valid) {
      continue;
    }

    const auto active_edge =
      find_active_edge(plan, neighbors, triangle_index, &front_edges);
    if(active_edge >= 3U) {
      continue;
    }

    const double radius = effective_triangle_front_radius(triangle);
    if(!std::isfinite(radius)) {
      continue;
    }

    if(selected.triangle_index != invalid_index &&
       radius <= selected.priority) {
      continue;
    }

    selected.triangle_index = triangle_index;
    selected.edge_slot = active_edge;
    selected.priority = radius;
  }

  return selected;
}

void add_local_triangle(
  FaceMeshPlan &plan,
  std::array<std::uint32_t, 3> triangle_vertices,
  bool accepted
)
{
  if(!orient_triangle_ccw(plan.vertices, triangle_vertices)) {
    return;
  }

  plan.triangles.push_back({triangle_vertices, true, accepted});
}

[[nodiscard]] bool triangle_passes_source_face_screen(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  const FaceMeshPlan &plan,
  const std::array<std::uint32_t, 3> &triangle,
  LocalTriangleSourceFaceFailureKind *failure = nullptr
)
{
  if(failure != nullptr) {
    *failure = LocalTriangleSourceFaceFailureKind::none;
  }
  if(triangle[0] >= plan.vertices.size() ||
     triangle[1] >= plan.vertices.size() ||
     triangle[2] >= plan.vertices.size()) {
    if(failure != nullptr) {
      *failure = LocalTriangleSourceFaceFailureKind::node_reprojection_failure;
    }
    return false;
  }

  const auto &a = plan.vertices[triangle[0]].position;
  const auto &b = plan.vertices[triangle[1]].position;
  const auto &c = plan.vertices[triangle[2]].position;
  const std::array<geo::Point3, 3> triangle_points {
    a,
    b,
    c,
  };
  const auto &a_uv = plan.vertices[triangle[0]].uv;
  const auto &b_uv = plan.vertices[triangle[1]].uv;
  const auto &c_uv = plan.vertices[triangle[2]].uv;
  const auto ab = point_subtract(b, a);
  const auto bc = point_subtract(c, b);
  const auto ca = point_subtract(a, c);
  const double tolerance = projection_tolerance(std::max({
    std::sqrt(dot_product(ab, ab)),
    std::sqrt(dot_product(bc, bc)),
    std::sqrt(dot_product(ca, ca)),
  }));

  const auto point_projects_to_face =
    [&](const geo::Point3 &point, geo::FaceProjection *projection = nullptr) {
      geo::FaceProjection local_projection;
      auto &target_projection = projection == nullptr ? local_projection : *projection;
      return geo::project_point_to_face(face_view, point, target_projection) ==
               base::StatusCode::ok &&
             target_projection.distance <= tolerance;
    };

  if(!point_projects_to_face(a) ||
     !point_projects_to_face(b) ||
     !point_projects_to_face(c)) {
    if(failure != nullptr) {
      *failure = LocalTriangleSourceFaceFailureKind::node_reprojection_failure;
    }
    return false;
  }

  const geo::Point3 centroid {
    (a[0] + b[0] + c[0]) / 3.0,
    (a[1] + b[1] + c[1]) / 3.0,
    (a[2] + b[2] + c[2]) / 3.0,
  };
  const std::array<double, 2> centroid_uv {
    (a_uv[0] + b_uv[0] + c_uv[0]) / 3.0,
    (a_uv[1] + b_uv[1] + c_uv[1]) / 3.0,
  };
  const std::array<std::array<double, 2>, 3> triangle_uvs {
    a_uv,
    b_uv,
    c_uv,
  };
  if(uses_uv_material_face_screen(face_preprocess)) {
    if(!triangle_passes_uv_material_face_screen(
         face_view,
         face_preprocess,
         triangle_points,
         triangle_uvs,
         failure
       )) {
      return false;
    }
    if(triangle_overlaps_inner_loops(
         face_preprocess,
         plan,
         triangle,
         triangle_uvs
       )) {
      if(failure != nullptr) {
        *failure = LocalTriangleSourceFaceFailureKind::source_face_containment_failure;
      }
      return false;
    }
    return true;
  }

  geo::FaceProjection centroid_projection;
  const bool centroid_inside_material =
    uv_point_inside_face_material(face_preprocess, centroid_uv);
  if(!point_projects_to_face(centroid, &centroid_projection) ||
     !centroid_inside_material) {
    if(failure != nullptr) {
      *failure = LocalTriangleSourceFaceFailureKind::source_face_containment_failure;
    }
    return false;
  }

  if(triangle_overlaps_inner_loops(
       face_preprocess,
       plan,
       triangle,
       triangle_uvs
     )) {
    if(failure != nullptr) {
      *failure = LocalTriangleSourceFaceFailureKind::source_face_containment_failure;
    }
    return false;
  }

  if(centroid_projection.normal_defined) {
    const auto normal = normalized_vector(cross_product(ab, point_subtract(c, a)));
    if(dot_product(normal, centroid_projection.normal) <= 0.0) {
      if(failure != nullptr) {
        *failure = LocalTriangleSourceFaceFailureKind::orientation_flip;
      }
      return false;
    }
  }

  return true;
}

void record_plan_triangle_source_face_failures(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  const FaceMeshPlan &plan,
  testing::AutoCfdSurfaceSourceFaceFailureCounts &counts
) noexcept
{
  for(const auto &triangle : plan.triangles) {
    if(!triangle.valid) {
      continue;
    }

    LocalTriangleSourceFaceFailureKind failure =
      LocalTriangleSourceFaceFailureKind::none;
    if(!triangle_passes_source_face_screen(
         face_view,
         face_preprocess,
         plan,
         triangle.vertices,
         &failure
       )) {
      record_source_face_failure(counts, failure);
    }
  }
}

[[nodiscard]] double triangle_min_angle_3d(
  const geo::Point3 &pa,
  const geo::Point3 &pb,
  const geo::Point3 &pc
) noexcept
{
  const double abx = pb[0] - pa[0], aby = pb[1] - pa[1], abz = pb[2] - pa[2];
  const double acx = pc[0] - pa[0], acy = pc[1] - pa[1], acz = pc[2] - pa[2];
  const double bcx = pc[0] - pb[0], bcy = pc[1] - pb[1], bcz = pc[2] - pb[2];
  const double lab = std::sqrt(abx * abx + aby * aby + abz * abz);
  const double lac = std::sqrt(acx * acx + acy * acy + acz * acz);
  const double lbc = std::sqrt(bcx * bcx + bcy * bcy + bcz * bcz);
  if(lab <= kMetricTolerance || lac <= kMetricTolerance || lbc <= kMetricTolerance) {
    return 0.0;
  }

  const double a0 = angle_from_edges(lab, lac, lbc);
  const double a1 = angle_from_edges(lab, lbc, lac);
  const double a2 = angle_from_edges(lac, lbc, lab);
  return std::min({a0, a1, a2}) * (180.0 / kPi);
}

// Compute the minimum angle across all valid triangles containing a vertex.
[[nodiscard]] double vertex_star_min_angle(
  const FaceMeshPlan &plan,
  std::uint32_t vertex_index,
  const std::vector<std::vector<std::uint32_t>> &vertex_to_triangles
) noexcept
{
  double worst = 180.0;
  if(vertex_index >= vertex_to_triangles.size()) {
    return 0.0;
  }

  for(const auto triangle_index : vertex_to_triangles[vertex_index]) {
    if(triangle_index >= plan.triangles.size() || !plan.triangles[triangle_index].valid) {
      continue;
    }
    const auto &tri = plan.triangles[triangle_index];
    const double min_angle = triangle_min_angle_3d(
      plan.vertices[tri.vertices[0]].position,
      plan.vertices[tri.vertices[1]].position,
      plan.vertices[tri.vertices[2]].position
    );
    worst = std::min(worst, min_angle);
  }
  return worst;
}

// Build a vertex-to-triangle adjacency map.
void build_vertex_to_triangles(
  const FaceMeshPlan &plan,
  std::vector<std::vector<std::uint32_t>> &vertex_to_triangles
)
{
  vertex_to_triangles.assign(plan.vertices.size(), {});
  for(std::uint32_t tri_idx = 0U; tri_idx < plan.triangles.size(); ++tri_idx) {
    const auto &tri = plan.triangles[tri_idx];
    if(!tri.valid) {
      continue;
    }
    for(const auto v : tri.vertices) {
      if(v < vertex_to_triangles.size()) {
        vertex_to_triangles[v].push_back(tri_idx);
      }
    }
  }
}

// Laplacian smoothing: move each interior vertex toward the centroid of its
// neighbors, project onto the surface, and accept the move only if it
// improves the worst triangle quality in the vertex star.
void smooth_interior_vertices(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceSizingFieldState &field,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  FaceMeshPlan &plan,
  std::vector<std::array<std::uint32_t, 3>> &neighbors
)
{
  constexpr std::size_t kSmoothingPasses = 10U;
  constexpr double kRelaxation = 0.5; // blend factor: 0=stay, 1=full move
  constexpr double kEquilateralHeight = 0.86602540378443864676; // sqrt(3)/2

  std::vector<std::vector<std::uint32_t>> vertex_to_triangles;
  build_vertex_to_triangles(plan, vertex_to_triangles);

  std::size_t total_moved = 0U;
  for(std::size_t pass = 0U; pass < kSmoothingPasses; ++pass) {
    std::size_t moved_this_pass = 0U;

    for(std::uint32_t vi = 0U; vi < plan.vertices.size(); ++vi) {
      auto &vertex = plan.vertices[vi];
      if(vertex.boundary) {
        continue;
      }
      if(vi >= vertex_to_triangles.size() || vertex_to_triangles[vi].empty()) {
        continue;
      }

      // Optimization-based smoothing: for each triangle in the vertex star,
      // compute the "ideal" position that would make it equilateral, then
      // average weighted by inverse quality (worse triangles pull harder).
      std::array<double, 2> target_uv {0.0, 0.0};
      double total_weight = 0.0;

      for(const auto tri_idx : vertex_to_triangles[vi]) {
        if(tri_idx >= plan.triangles.size() || !plan.triangles[tri_idx].valid) {
          continue;
        }
        const auto &tri = plan.triangles[tri_idx];

        // Find the two other vertices (a, b) of this triangle
        std::uint32_t ai = invalid_index, bi = invalid_index;
        for(const auto tv : tri.vertices) {
          if(tv == vi) continue;
          if(ai == invalid_index) ai = tv;
          else bi = tv;
        }
        if(ai == invalid_index || bi == invalid_index ||
           ai >= plan.vertices.size() || bi >= plan.vertices.size()) {
          continue;
        }

        // Compute the ideal position for v to make triangle (v, a, b) equilateral.
        // In UV space: midpoint of edge (a, b), then offset by sqrt(3)/2 * edge_length
        // in the perpendicular direction (toward v's current side).
        const auto &a_uv = plan.vertices[ai].uv;
        const auto &b_uv = plan.vertices[bi].uv;
        const std::array<double, 2> mid {
          0.5 * (a_uv[0] + b_uv[0]),
          0.5 * (a_uv[1] + b_uv[1]),
        };
        const std::array<double, 2> edge {
          b_uv[0] - a_uv[0],
          b_uv[1] - a_uv[1],
        };
        // Perpendicular direction (90° CCW rotation)
        std::array<double, 2> perp {-edge[1], edge[0]};
        const double perp_len = std::sqrt(perp[0] * perp[0] + perp[1] * perp[1]);
        if(perp_len <= kMetricTolerance) {
          continue;
        }
        perp[0] /= perp_len;
        perp[1] /= perp_len;

        // Choose the side where v currently is
        const std::array<double, 2> to_v {
          vertex.uv[0] - mid[0],
          vertex.uv[1] - mid[1],
        };
        const double side = to_v[0] * perp[0] + to_v[1] * perp[1];
        if(side < 0.0) {
          perp[0] = -perp[0];
          perp[1] = -perp[1];
        }

        const double ideal_height = perp_len * kEquilateralHeight;
        const std::array<double, 2> ideal {
          mid[0] + ideal_height * perp[0],
          mid[1] + ideal_height * perp[1],
        };

        // Weight by inverse of current min angle — worse triangles pull harder
        const double min_angle = triangle_min_angle_3d(
          vertex.position,
          plan.vertices[ai].position,
          plan.vertices[bi].position
        );
        const double w = 1.0 / std::max(min_angle, 0.1);

        target_uv[0] += w * ideal[0];
        target_uv[1] += w * ideal[1];
        total_weight += w;
      }

      if(total_weight <= 1.0e-30) {
        continue;
      }
      target_uv[0] /= total_weight;
      target_uv[1] /= total_weight;

      // Blend between current and target position (relaxation)
      std::array<double, 2> new_uv {
        vertex.uv[0] + kRelaxation * (target_uv[0] - vertex.uv[0]),
        vertex.uv[1] + kRelaxation * (target_uv[1] - vertex.uv[1]),
      };

      if(!uv_point_inside_face_material(face_preprocess, new_uv)) {
        continue;
      }

      geo::FaceSample sample;
      if(geo::sample_face(face_view, new_uv[0], new_uv[1], sample) !=
         base::StatusCode::ok) {
        continue;
      }

      // Evaluate quality before move
      const double quality_before = vertex_star_min_angle(plan, vi, vertex_to_triangles);

      // Tentatively move
      const auto old_position = vertex.position;
      const auto old_uv = vertex.uv;
      vertex.position = sample.position;
      vertex.uv = new_uv;

      // Check for degenerate/inverted triangles
      bool degenerate = false;
      for(const auto tri_idx : vertex_to_triangles[vi]) {
        if(tri_idx >= plan.triangles.size() || !plan.triangles[tri_idx].valid) {
          continue;
        }
        const auto &tri = plan.triangles[tri_idx];
        const auto d01 = uv_subtract(plan.vertices[tri.vertices[1]].uv, plan.vertices[tri.vertices[0]].uv);
        const auto d02 = uv_subtract(plan.vertices[tri.vertices[2]].uv, plan.vertices[tri.vertices[0]].uv);
        if(uv_cross(d01, d02) <= kUvAreaTolerance) {
          degenerate = true;
          break;
        }
      }

      if(degenerate) {
        vertex.position = old_position;
        vertex.uv = old_uv;
        continue;
      }

      // Accept only if quality improves
      const double quality_after = vertex_star_min_angle(plan, vi, vertex_to_triangles);
      if(quality_after <= quality_before) {
        vertex.position = old_position;
        vertex.uv = old_uv;
      }
      else {
        vertex.local_size = query_auto_cfd_surface_sizing_field(
          field, face_view, new_uv[0], new_uv[1], sample.position
        );
        vertex.background_size = vertex.local_size;
        ++moved_this_pass;
      }
    }

    total_moved += moved_this_pass;
    if(moved_this_pass == 0U) {
      break;
    }
  }

  if(total_moved > 0U) {
    build_triangle_neighbors(plan.triangles, neighbors);
  }
}

// Quality-based edge swapping: for each non-constrained, non-boundary internal
// edge, swap the diagonal if it improves the minimum angle of the
// quadrilateral formed by two adjacent triangles.
void quality_edge_swaps(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  FaceMeshPlan &plan,
  std::vector<std::array<std::uint32_t, 3>> &neighbors
)
{
  constexpr std::size_t kMaxSwapPasses = 5U;
  constexpr double kMinImprovement = 0.5; // only swap if >= 0.5° improvement

  std::size_t total_swaps = 0U;
  std::vector<bool> dirty; // marks triangles modified in the current pass

  for(std::size_t pass = 0U; pass < kMaxSwapPasses; ++pass) {
    std::size_t swaps_this_pass = 0U;
    dirty.assign(plan.triangles.size(), false);

    for(std::uint32_t ti = 0U; ti < plan.triangles.size(); ++ti) {
      if(!plan.triangles[ti].valid || dirty[ti]) {
        continue;
      }

      for(std::uint8_t edge = 0U; edge < 3U; ++edge) {
        const auto ni = neighbors[ti][edge];
        if(ni == invalid_index || ni >= plan.triangles.size() ||
           !plan.triangles[ni].valid || dirty[ni]) {
          continue;
        }
        // Only process each edge once (lower index first)
        if(ni < ti) {
          continue;
        }

        const auto &tri = plan.triangles[ti];
        const auto v0 = tri.vertices[edge];
        const auto v1 = tri.vertices[(edge + 1U) % 3U];
        const auto v2 = tri.vertices[(edge + 2U) % 3U]; // opposite in tri

        // Check constrained
        if(plan.constrained_edges.find(canonical_edge(v0, v1)) !=
           plan.constrained_edges.end()) {
          continue;
        }

        // Find the opposite vertex in the neighbor triangle
        const auto &ntri = plan.triangles[ni];
        std::uint32_t v3 = invalid_index;
        for(const auto nv : ntri.vertices) {
          if(nv != v0 && nv != v1) {
            v3 = nv;
            break;
          }
        }
        if(v3 == invalid_index || v3 >= plan.vertices.size()) {
          continue;
        }

        
        if(plan.vertices[v2].periodic_counterpart != invalid_index &&
           plan.vertices[v3].periodic_counterpart != invalid_index) {
          continue;
        }

        // Current min angle
        const double ma_t1 = triangle_min_angle_3d(
          plan.vertices[v0].position,
          plan.vertices[v1].position,
          plan.vertices[v2].position
        );
        const double ma_t2 = triangle_min_angle_3d(
          plan.vertices[v0].position,
          plan.vertices[v1].position,
          plan.vertices[v3].position
        );
        const double current_min = std::min(ma_t1, ma_t2);

        // Convexity check: v0 and v1 must be on OPPOSITE sides of the
        // v2→v3 diagonal. cross(v3-v2, v0-v2) and cross(v3-v2, v1-v2) must
        // have opposite signs.
        const auto d23 = uv_subtract(plan.vertices[v3].uv, plan.vertices[v2].uv);
        const auto d20 = uv_subtract(plan.vertices[v0].uv, plan.vertices[v2].uv);
        const auto d21 = uv_subtract(plan.vertices[v1].uv, plan.vertices[v2].uv);
        const double cross_v0 = uv_cross(d23, d20); // v0 side
        const double cross_v1 = uv_cross(d23, d21); // v1 side
        if(cross_v0 * cross_v1 >= 0.0) {
          continue; // not convex — both on same side of diagonal
        }

        // Swapped min angle
        const double ma_s1 = triangle_min_angle_3d(
          plan.vertices[v2].position,
          plan.vertices[v3].position,
          plan.vertices[v0].position
        );
        const double ma_s2 = triangle_min_angle_3d(
          plan.vertices[v2].position,
          plan.vertices[v3].position,
          plan.vertices[v1].position
        );
        const double swapped_min = std::min(ma_s1, ma_s2);

        if(swapped_min <= current_min + kMinImprovement) {
          continue;
        }

        // Save old vertices for potential revert
        const auto old_ti_verts = tri.vertices;
        const auto old_ni_verts = ntri.vertices;

        // Perform the swap with correct CCW orientation:
        // cross_v0 > 0 means v0 is "left" of v2→v3, so (v2,v3,v0) is CCW.
        // cross_v1 < 0 means v1 is "right" of v2→v3, so (v2,v1,v3) is CCW.
        auto &t1 = plan.triangles[ti];
        auto &t2 = plan.triangles[ni];
        if(cross_v0 > 0.0) {
          t1.vertices = {v2, v3, v0};
        } else {
          t1.vertices = {v2, v0, v3};
        }
        if(cross_v1 > 0.0) {
          t2.vertices = {v2, v3, v1};
        } else {
          t2.vertices = {v2, v1, v3};
        }

        // Verify both triangles have proper (positive) area after swap
        {
          const auto c1a = uv_subtract(plan.vertices[t1.vertices[1]].uv, plan.vertices[t1.vertices[0]].uv);
          const auto c1b = uv_subtract(plan.vertices[t1.vertices[2]].uv, plan.vertices[t1.vertices[0]].uv);
          const auto c2a = uv_subtract(plan.vertices[t2.vertices[1]].uv, plan.vertices[t2.vertices[0]].uv);
          const auto c2b = uv_subtract(plan.vertices[t2.vertices[2]].uv, plan.vertices[t2.vertices[0]].uv);
          if(uv_cross(c1a, c1b) <= kUvAreaTolerance ||
             uv_cross(c2a, c2b) <= kUvAreaTolerance) {
            // Revert using saved copies
            t1.vertices = old_ti_verts;
            t2.vertices = old_ni_verts;
            continue;
          }
        }

        // Mark both triangles as dirty to avoid re-processing in this pass
        dirty[ti] = true;
        dirty[ni] = true;
        ++swaps_this_pass;
        break; // move to next triangle (edges of ti are now invalid)
      }
    }

    total_swaps += swaps_this_pass;
    if(swaps_this_pass == 0U) {
      break;
    }

    // Rebuild neighbors after swaps
    build_triangle_neighbors(plan.triangles, neighbors);
  }
}

// Edge collapse: remove short edges by merging endpoints. For each interior
// edge shorter than collapse_ratio * local_size, merge the two endpoints
// into their midpoint (projected onto the surface), invalidating the two
// adjacent triangles and rewiring all references to the removed vertex.
std::size_t edge_collapse_pass(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceSizingFieldState &field,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  FaceMeshPlan &plan,
  std::vector<std::array<std::uint32_t, 3>> &neighbors
)
{
  constexpr double kCollapseRatio = 0.5;

  std::vector<std::vector<std::uint32_t>> vertex_to_triangles;
  build_vertex_to_triangles(plan, vertex_to_triangles);

  std::size_t collapse_count = 0U;

  for(std::uint32_t ti = 0U; ti < plan.triangles.size(); ++ti) {
    if(!plan.triangles[ti].valid) {
      continue;
    }

    for(std::uint8_t edge = 0U; edge < 3U; ++edge) {
      const auto v0 = plan.triangles[ti].vertices[edge];
      const auto v1 = plan.triangles[ti].vertices[(edge + 1U) % 3U];
      if(v0 >= plan.vertices.size() || v1 >= plan.vertices.size()) {
        continue;
      }

      // Skip constrained (boundary) edges
      if(plan.constrained_edges.find(canonical_edge(v0, v1)) !=
         plan.constrained_edges.end()) {
        continue;
      }

      // Skip if either vertex is boundary
      if(plan.vertices[v0].boundary || plan.vertices[v1].boundary) {
        continue;
      }

      // Skip if either vertex has periodic counterpart
      if(plan.vertices[v0].periodic_counterpart != invalid_index ||
         plan.vertices[v1].periodic_counterpart != invalid_index) {
        continue;
      }

      // Compute 3D edge length
      const auto &p0 = plan.vertices[v0].position;
      const auto &p1 = plan.vertices[v1].position;
      const double dx = p1[0] - p0[0];
      const double dy = p1[1] - p0[1];
      const double dz = p1[2] - p0[2];
      const double edge_length = std::sqrt(dx * dx + dy * dy + dz * dz);

      // Compare against local size threshold
      const double target = 0.5 * (plan.vertices[v0].local_size +
                                   plan.vertices[v1].local_size);
      if(target <= kMetricTolerance || edge_length >= kCollapseRatio * target) {
        continue;
      }

      // Compute midpoint UV and sample surface
      const std::array<double, 2> mid_uv {
        0.5 * (plan.vertices[v0].uv[0] + plan.vertices[v1].uv[0]),
        0.5 * (plan.vertices[v0].uv[1] + plan.vertices[v1].uv[1]),
      };
      if(!uv_point_inside_face_material(face_preprocess, mid_uv)) {
        continue;
      }
      geo::FaceSample sample;
      if(geo::sample_face(face_view, mid_uv[0], mid_uv[1], sample) !=
         base::StatusCode::ok) {
        continue;
      }

      // Check quality: compute min angle of ALL triangles in the stars of
      // v0 and v1 BEFORE the collapse.
      double quality_before = 180.0;
      std::unordered_set<std::uint32_t> affected_tris;
      for(const auto idx : vertex_to_triangles[v0]) {
        if(idx < plan.triangles.size() && plan.triangles[idx].valid) {
          affected_tris.insert(idx);
        }
      }
      for(const auto idx : vertex_to_triangles[v1]) {
        if(idx < plan.triangles.size() && plan.triangles[idx].valid) {
          affected_tris.insert(idx);
        }
      }
      for(const auto idx : affected_tris) {
        const auto &tri = plan.triangles[idx];
        quality_before = std::min(quality_before, triangle_min_angle_3d(
          plan.vertices[tri.vertices[0]].position,
          plan.vertices[tri.vertices[1]].position,
          plan.vertices[tri.vertices[2]].position
        ));
      }

      // Save original state for potential revert
      const auto old_v0 = plan.vertices[v0];
      std::vector<std::pair<std::uint32_t, std::array<std::uint32_t, 3>>> saved_tris;
      saved_tris.reserve(affected_tris.size());
      for(const auto idx : affected_tris) {
        saved_tris.push_back({idx, plan.triangles[idx].vertices});
      }

      // Tentatively move v0 to midpoint
      plan.vertices[v0].position = sample.position;
      plan.vertices[v0].uv = mid_uv;

      // Identify triangles to remove (those containing BOTH v0 and v1)
      std::vector<std::uint32_t> tris_to_remove;
      for(const auto idx : affected_tris) {
        const auto &tri = plan.triangles[idx];
        bool has_v0 = false, has_v1 = false;
        for(const auto v : tri.vertices) {
          if(v == v0) has_v0 = true;
          if(v == v1) has_v1 = true;
        }
        if(has_v0 && has_v1) {
          tris_to_remove.push_back(idx);
        }
      }

      // Rewire v1 → v0 in all affected triangles
      for(const auto idx : affected_tris) {
        for(auto &v : plan.triangles[idx].vertices) {
          if(v == v1) v = v0;
        }
      }

      // Check quality after: skip removed triangles, check for degenerate
      bool degenerate = false;
      double quality_after = 180.0;
      for(const auto idx : affected_tris) {
        if(std::find(tris_to_remove.begin(), tris_to_remove.end(), idx) !=
           tris_to_remove.end()) {
          continue;
        }
        const auto &tri = plan.triangles[idx];
        if(tri.vertices[0] == tri.vertices[1] ||
           tri.vertices[0] == tri.vertices[2] ||
           tri.vertices[1] == tri.vertices[2]) {
          degenerate = true;
          break;
        }
        const auto d01 = uv_subtract(plan.vertices[tri.vertices[1]].uv,
                                      plan.vertices[tri.vertices[0]].uv);
        const auto d02 = uv_subtract(plan.vertices[tri.vertices[2]].uv,
                                      plan.vertices[tri.vertices[0]].uv);
        if(uv_cross(d01, d02) <= kUvAreaTolerance) {
          degenerate = true;
          break;
        }
        quality_after = std::min(quality_after, triangle_min_angle_3d(
          plan.vertices[tri.vertices[0]].position,
          plan.vertices[tri.vertices[1]].position,
          plan.vertices[tri.vertices[2]].position
        ));
      }

      if(degenerate || quality_after < quality_before) {
        // Revert: restore vertex and all triangle vertices
        plan.vertices[v0] = old_v0;
        for(const auto &[idx, verts] : saved_tris) {
          plan.triangles[idx].vertices = verts;
        }
        continue;
      }

      // Accept collapse: invalidate removed triangles
      for(const auto idx : tris_to_remove) {
        plan.triangles[idx].valid = false;
      }

      // Update sizing at collapsed vertex
      plan.vertices[v0].local_size = query_auto_cfd_surface_sizing_field(
        field, face_view, mid_uv[0], mid_uv[1], sample.position
      );
      plan.vertices[v0].background_size = plan.vertices[v0].local_size;

      ++collapse_count;
      // Rebuild adjacency and restart
      build_triangle_neighbors(plan.triangles, neighbors);
      build_vertex_to_triangles(plan, vertex_to_triangles);
      break; // restart outer loop since topology changed
    }

    // If a collapse happened, we broke out of the edge loop; restart from ti=0
    if(collapse_count > 0U &&
       ti + 1U < plan.triangles.size()) {
      // Continue scanning from the next triangle
    }
  }

  if(collapse_count > 0U) {
    build_triangle_neighbors(plan.triangles, neighbors);
  }
  return collapse_count;
}

// Edge split: split long edges by inserting a midpoint. For each interior
// edge longer than split_ratio * local_size, insert a new vertex at the
// edge midpoint (projected onto the surface) and replace the two adjacent
// triangles with four.
std::size_t edge_split_pass(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceSizingFieldState &field,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  FaceMeshPlan &plan,
  std::vector<std::array<std::uint32_t, 3>> &neighbors
)
{
  constexpr double kSplitRatio = 1.5;

  std::size_t split_count = 0U;

  for(std::uint32_t ti = 0U; ti < plan.triangles.size(); ++ti) {
    if(!plan.triangles[ti].valid) {
      continue;
    }

    for(std::uint8_t edge = 0U; edge < 3U; ++edge) {
      const auto ni = neighbors[ti][edge];
      if(ni == invalid_index || ni >= plan.triangles.size() ||
         !plan.triangles[ni].valid) {
        continue;
      }
      // Process each edge once
      if(ni < ti) {
        continue;
      }

      const auto v0 = plan.triangles[ti].vertices[edge];
      const auto v1 = plan.triangles[ti].vertices[(edge + 1U) % 3U];
      const auto v2 = plan.triangles[ti].vertices[(edge + 2U) % 3U];
      if(v0 >= plan.vertices.size() || v1 >= plan.vertices.size()) {
        continue;
      }

      // Skip constrained (boundary) edges
      if(plan.constrained_edges.find(canonical_edge(v0, v1)) !=
         plan.constrained_edges.end()) {
        continue;
      }

      // Compute 3D edge length
      const auto &p0 = plan.vertices[v0].position;
      const auto &p1 = plan.vertices[v1].position;
      const double dx = p1[0] - p0[0];
      const double dy = p1[1] - p0[1];
      const double dz = p1[2] - p0[2];
      const double edge_length = std::sqrt(dx * dx + dy * dy + dz * dz);

      // Compare against local size threshold
      const double target = 0.5 * (plan.vertices[v0].local_size +
                                   plan.vertices[v1].local_size);
      if(target <= kMetricTolerance || edge_length <= kSplitRatio * target) {
        continue;
      }

      // Find opposite vertex in neighbor triangle
      const auto &ntri = plan.triangles[ni];
      std::uint32_t v3 = invalid_index;
      for(const auto nv : ntri.vertices) {
        if(nv != v0 && nv != v1) {
          v3 = nv;
          break;
        }
      }
      if(v3 == invalid_index || v3 >= plan.vertices.size()) {
        continue;
      }

      // Compute midpoint UV and sample surface
      const std::array<double, 2> mid_uv {
        0.5 * (plan.vertices[v0].uv[0] + plan.vertices[v1].uv[0]),
        0.5 * (plan.vertices[v0].uv[1] + plan.vertices[v1].uv[1]),
      };
      if(!uv_point_inside_face_material(face_preprocess, mid_uv)) {
        continue;
      }
      geo::FaceSample sample;
      if(geo::sample_face(face_view, mid_uv[0], mid_uv[1], sample) !=
         base::StatusCode::ok) {
        continue;
      }

      // Create new vertex at midpoint
      const auto vm = static_cast<std::uint32_t>(plan.vertices.size());
      const double mid_size = query_auto_cfd_surface_sizing_field(
        field, face_view, mid_uv[0], mid_uv[1], sample.position
      );
      plan.vertices.push_back({
        invalid_index,
        sample.position,
        mid_uv,
        false,
        mid_size,
        mid_size,
      });

      // Replace 2 triangles with 4:
      //   tri (v0, v1, v2) → (v0, vm, v2) + (vm, v1, v2)
      //   ntri(v0, v1, v3) → (v0, vm, v3) + (vm, v1, v3)
      // (preserving winding by checking which triangle has which vertex order)
      auto &t1 = plan.triangles[ti];
      auto &t2 = plan.triangles[ni];

      // Determine correct winding for new triangles
      // t1 originally has (v0, v1, v2) in some CCW order
      // t2 originally has (v0, v1, v3) in some CCW order (v0,v1 reversed)
      const auto old_t1 = t1.vertices;
      const auto old_t2 = t2.vertices;

      // Find the edge position in t1
      std::uint32_t t1_e = 0U;
      for(std::uint32_t i = 0U; i < 3U; ++i) {
        if((old_t1[i] == v0 && old_t1[(i+1)%3] == v1) ||
           (old_t1[i] == v1 && old_t1[(i+1)%3] == v0)) {
          t1_e = i;
          break;
        }
      }
      const auto t1_a = old_t1[t1_e];
      const auto t1_b = old_t1[(t1_e+1)%3];
      const auto t1_c = old_t1[(t1_e+2)%3];

      std::uint32_t t2_e = 0U;
      for(std::uint32_t i = 0U; i < 3U; ++i) {
        if((old_t2[i] == v0 && old_t2[(i+1)%3] == v1) ||
           (old_t2[i] == v1 && old_t2[(i+1)%3] == v0)) {
          t2_e = i;
          break;
        }
      }
      const auto t2_a = old_t2[t2_e];
      const auto t2_b = old_t2[(t2_e+1)%3];
      const auto t2_c = old_t2[(t2_e+2)%3];

      // Quality check: the 4 new triangles must not be worse than the 2 originals
      const double quality_before = std::min(
        triangle_min_angle_3d(
          plan.vertices[t1_a].position, plan.vertices[t1_b].position, plan.vertices[t1_c].position),
        triangle_min_angle_3d(
          plan.vertices[t2_a].position, plan.vertices[t2_b].position, plan.vertices[t2_c].position)
      );

      // t1: (t1_a, t1_b, t1_c) → (t1_a, vm, t1_c) + new (vm, t1_b, t1_c)
      // t2: (t2_a, t2_b, t2_c) → (t2_a, vm, t2_c) + new (vm, t2_b, t2_c)
      const double quality_after = std::min({
        triangle_min_angle_3d(
          plan.vertices[t1_a].position, sample.position, plan.vertices[t1_c].position),
        triangle_min_angle_3d(
          sample.position, plan.vertices[t1_b].position, plan.vertices[t1_c].position),
        triangle_min_angle_3d(
          plan.vertices[t2_a].position, sample.position, plan.vertices[t2_c].position),
        triangle_min_angle_3d(
          sample.position, plan.vertices[t2_b].position, plan.vertices[t2_c].position),
      });

      if(quality_after < quality_before) {
        continue;
      }

      t1.vertices = {t1_a, vm, t1_c};
      t2.vertices = {t2_a, vm, t2_c};
      plan.triangles.push_back({{vm, t1_b, t1_c}, true, false});
      plan.triangles.push_back({{vm, t2_b, t2_c}, true, false});

      ++split_count;
      build_triangle_neighbors(plan.triangles, neighbors);
      break;
    }
  }

  if(split_count > 0U) {
    build_triangle_neighbors(plan.triangles, neighbors);
  }
  return split_count;
}

// ── End Post-BW quality optimization ────────────────────────────────────────

void run_front_refinement(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceSizingFieldState &field,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  FaceMeshPlan &plan
)
{
  plan.inserted_vertex_count = 0U;
  plan.front_iteration_count = 0U;
  plan.inserted_vertex_limit_reached = false;
  plan.front_iteration_limit_reached = false;
  bool front_converged = false;

  std::vector<std::array<std::uint32_t, 3>> neighbors;
  // Reusable buffers to avoid per-iteration allocation
  std::vector<bool> deleted;
  std::vector<std::uint32_t> cavity;
  std::vector<CavityShellEdge> shell;
  FrontTriangleSet all_triangles(
    FrontTrianglePriorityComparator {&plan.triangles}
  );
  FrontTriangleSet active_triangles(
    FrontTrianglePriorityComparator {&plan.triangles}
  );

  // std::cerr << "[BW] Face " << face_preprocess.face.index
  //           << ": starting frontal refinement, seed_tris="
  //           << plan.triangles.size() << " seed_verts=" << plan.vertices.size()
  //           << " target_size=" << plan.target_size
  //           << " planar_fallback=" << plan.planar_fallback << "\n";
  // std::cerr << "[BW] Face " << face_preprocess.face.index
  //           << ": uv_extent recovered=["
  //           << face_preprocess.recovered_uv_min[0] << ","
  //           << face_preprocess.recovered_uv_max[0] << "]x["
  //           << face_preprocess.recovered_uv_min[1] << ","
  //           << face_preprocess.recovered_uv_max[1] << "]"
  //           << " bounds=["
  //           << face_preprocess.uv_bounds.u_min << ","
  //           << face_preprocess.uv_bounds.u_max << "]x["
  //           << face_preprocess.uv_bounds.v_min << ","
  //           << face_preprocess.uv_bounds.v_max << "]"
  //           << " ref_uv=(" << face_preprocess.reference_metric.u
  //           << "," << face_preprocess.reference_metric.v << ")"
  //           << " ref_target=" << face_preprocess.reference_metric.target_size
  //           << "\n";

  refresh_all_vertex_front_sizes(face_view, field, plan);
  refresh_all_triangle_front_radii(face_view, field, face_preprocess, plan);
  build_triangle_neighbors(plan.triangles, neighbors);

  for(std::uint32_t triangle_index = 0U;
      triangle_index < plan.triangles.size();
      ++triangle_index) {
    if(plan.triangles[triangle_index].valid) {
      all_triangles.insert(triangle_index);
    }
  }

  for(auto triangle_it = all_triangles.begin();
      triangle_it != all_triangles.end();
      ++triangle_it) {
    const auto triangle_index = *triangle_it;
    const auto &triangle = plan.triangles[triangle_index];
    if(!triangle.valid) {
      continue;
    }

    const auto active_edge =
      find_active_edge(plan, neighbors, triangle_index, nullptr);
    if(active_edge >= 3U ||
       effective_triangle_front_radius(triangle) <= kFrontRadiusLimit) {
      if(effective_triangle_front_radius(triangle) < kFrontRadiusLimit) {
        break;
      }
      continue;
    }
    active_triangles.insert(triangle_index);
  }

  std::size_t iteration = 0U;
  std::size_t skip_count = 0U;
  std::size_t insert_ok_count = 0U;
  std::size_t insert_fail_count = 0U;
  std::size_t true_boundary_reject_count = 0U;
  long long ms_find_edge = 0, ms_optimal = 0, ms_insert = 0, ms_diag = 0;
  while(iteration < kMaximumFrontIterationsPerFace) {
    if(active_triangles.empty()) {
      front_converged = true;
      break;
    }

    const auto active_it = active_triangles.begin();
    const auto triangle_index = *active_it;
    active_triangles.erase(active_it);

    if(triangle_index >= plan.triangles.size() ||
       !plan.triangles[triangle_index].valid) {
      continue;
    }

    auto te0 = std::chrono::steady_clock::now();
    const auto front_edge_slot =
      find_active_edge(plan, neighbors, triangle_index, nullptr);
    const double front_priority =
      effective_triangle_front_radius(plan.triangles[triangle_index]);
    auto te1 = std::chrono::steady_clock::now();
    ms_find_edge += std::chrono::duration_cast<std::chrono::microseconds>(te1 - te0).count();
    if(front_edge_slot >= 3U || front_priority <= kFrontRadiusLimit) {
      ++skip_count;
      continue;
    }

    const auto current_iteration = iteration++;
    ++plan.front_iteration_count;

    if(current_iteration < 20U || current_iteration % 500U == 0U) {

      const auto &triangle = plan.triangles[triangle_index];
      const std::array<double, 2> centroid_uv {
        (plan.vertices[triangle.vertices[0]].uv[0] +
         plan.vertices[triangle.vertices[1]].uv[0] +
         plan.vertices[triangle.vertices[2]].uv[0]) / 3.0,
        (plan.vertices[triangle.vertices[0]].uv[1] +
         plan.vertices[triangle.vertices[1]].uv[1] +
         plan.vertices[triangle.vertices[2]].uv[1]) / 3.0,
      };
      geo::Point3 centroid_position {
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
      };
      geo::FaceSample centroid_sample;
      if(sample_face_at_material_uv(
           face_view,
           face_preprocess,
           centroid_uv,
           centroid_sample
         )) {
        centroid_position = centroid_sample.position;
      }
      MetricSampleTrace metric_trace;
      AutoCfdSurfaceFaceMetricTensor centroid_metric;
      static_cast<void>(sample_metric_or_reference(
        face_view,
        field,
        face_preprocess,
        &plan,
        centroid_uv[0],
        centroid_uv[1],
        centroid_metric,
        &metric_trace
      ));
      double owner_size = std::numeric_limits<double>::quiet_NaN();
      double uv_size = std::numeric_limits<double>::quiet_NaN();
      if(std::isfinite(centroid_position[0])) {
        owner_size = query_auto_cfd_surface_sizing_field(
          field,
          face_preprocess.face,
          centroid_position
        );
        uv_size = query_auto_cfd_surface_sizing_field(
          field,
          face_view,
          centroid_uv[0],
          centroid_uv[1],
          centroid_position
        );
      }
    }

    if(plan.inserted_vertex_count >= kMaximumInsertedVerticesPerFace) {
      plan.inserted_vertex_limit_reached = true;
      break;
    }

    const auto worst_vertices = plan.triangles[triangle_index].vertices;

    auto to0 = std::chrono::steady_clock::now();
    std::array<double, 2> candidate_uv {};
    geo::Point3 candidate_position {0.0, 0.0, 0.0};
    std::array<double, 3> metric_tensor {1.0, 0.0, 1.0};
    if(!compute_frontal_optimal_point(
         face_view, field, face_preprocess, plan,
         worst_vertices, front_edge_slot,
         candidate_uv, candidate_position,
         &metric_tensor
       )) {
      continue;
    }

    auto to1 = std::chrono::steady_clock::now();
    ms_optimal += std::chrono::duration_cast<std::chrono::microseconds>(to1 - to0).count();

    if(!uv_point_inside_face_material(face_preprocess, candidate_uv)) {
      continue;
    }

    if(plan.periodic_true_boundary_defined &&
       !point_inside_periodic_true_boundary(
         plan.periodic_true_boundary_segments,
         plan.periodic_true_boundary_far,
         candidate_uv
       )) {
      ++true_boundary_reject_count;
      continue;
    }

    const auto vertex_count_before_insertion =
      static_cast<std::uint32_t>(plan.vertices.size());
    const auto triangle_count_before_insertion =
      static_cast<std::uint32_t>(plan.triangles.size());

    auto ti0 = std::chrono::steady_clock::now();
    const auto insertion_succeeded = insert_a_point(
      face_view,
      field,
      face_preprocess,
      plan,
      neighbors,
      all_triangles,
      active_triangles,
      triangle_index,
      metric_tensor,
      candidate_uv,
      candidate_position,
      deleted,
      cavity,
      shell
    );

    auto ti1 = std::chrono::steady_clock::now();
    ms_insert += std::chrono::duration_cast<std::chrono::microseconds>(ti1 - ti0).count();

    if(!insertion_succeeded) {
      ++insert_fail_count;
      continue;
    }

    ++insert_ok_count;
    ++plan.inserted_vertex_count;

    if(current_iteration < 20U &&
       plan.vertices.size() > vertex_count_before_insertion) {
    }
  }

  if(!front_converged &&
     plan.front_iteration_count >= kMaximumFrontIterationsPerFace) {
    plan.front_iteration_limit_reached = true;
  }

  log_metric_sampling_summary(face_preprocess, plan, "post_bw");

 
#if 0
  {
    constexpr std::size_t kMaxOptPasses = 5U;
    auto t_opt_start = std::chrono::steady_clock::now();
    std::size_t total_collapses = 0U;
    std::size_t total_splits = 0U;

    for(std::size_t opt_pass = 0U; opt_pass < kMaxOptPasses; ++opt_pass) {
      smooth_interior_vertices(face_view, field, face_preprocess, plan, neighbors);
      quality_edge_swaps(face_view, face_preprocess, plan, neighbors);

      const auto collapses = edge_collapse_pass(
        face_view, field, face_preprocess, plan, neighbors
      );
      const auto splits = edge_split_pass(
        face_view, field, face_preprocess, plan, neighbors
      );
      total_collapses += collapses;
      total_splits += splits;

      if(collapses == 0U && splits == 0U) {
        break;
      }
    }

    auto t_opt_end = std::chrono::steady_clock::now();
  }
#endif

  accept_all_valid_triangles(plan);
}

[[nodiscard]] double orientation2d(
  const std::array<double, 2> &a,
  const std::array<double, 2> &b,
  const std::array<double, 2> &c
) noexcept
{
  return signed_triangle_area_twice(a, b, c);
}


[[nodiscard]] bool circumcenter_metric(
  const std::array<double, 2> &pa,
  const std::array<double, 2> &pb,
  const std::array<double, 2> &pc,
  const std::array<double, 3> &metric, // {a, b, d}  (E, F, G of first fundamental form)
  std::array<double, 2> &center,
  double &radius_squared
) noexcept
{
  const double a = metric[0];
  const double b = metric[1];
  const double d = metric[2];

  // 2×2 linear system: sys · x = rhs
  const double sys00 = 2.0 * a * (pa[0] - pb[0]) + 2.0 * b * (pa[1] - pb[1]);
  const double sys01 = 2.0 * d * (pa[1] - pb[1]) + 2.0 * b * (pa[0] - pb[0]);
  const double sys10 = 2.0 * a * (pa[0] - pc[0]) + 2.0 * b * (pa[1] - pc[1]);
  const double sys11 = 2.0 * d * (pa[1] - pc[1]) + 2.0 * b * (pa[0] - pc[0]);

  const double rhs0 =
    a * (pa[0]*pa[0] - pb[0]*pb[0]) +
    d * (pa[1]*pa[1] - pb[1]*pb[1]) +
    2.0 * b * (pa[0]*pa[1] - pb[0]*pb[1]);
  const double rhs1 =
    a * (pa[0]*pa[0] - pc[0]*pc[0]) +
    d * (pa[1]*pa[1] - pc[1]*pc[1]) +
    2.0 * b * (pa[0]*pa[1] - pc[0]*pc[1]);

  const double det = sys00 * sys11 - sys01 * sys10;
  if(std::abs(det) < 1.0e-30) {
    center = {0.0, 0.0};
    radius_squared = 1.0e22;
    return false;
  }

  center[0] = (rhs0 * sys11 - rhs1 * sys01) / det;
  center[1] = (sys00 * rhs1 - sys10 * rhs0) / det;

  const double dx = center[0] - pa[0];
  const double dy = center[1] - pa[1];
  radius_squared = dx * dx * a + dy * dy * d + 2.0 * dx * dy * b;
  return true;
}

[[nodiscard]] bool point_in_circumcircle_aniso(
  const std::array<double, 2> &a,
  const std::array<double, 2> &b,
  const std::array<double, 2> &c,
  const std::array<double, 2> &p,
  const std::array<double, 3> &metric // {E, F, G}
) noexcept
{
  std::array<double, 2> center;
  double radius_squared;
  if(!circumcenter_metric(a, b, c, metric, center, radius_squared)) {
    return false;
  }

  const double m_a = metric[0];
  const double m_b = metric[1];
  const double m_d = metric[2];
  const double d0 = center[0] - p[0];
  const double d1 = center[1] - p[1];
  const double dist_sq = d0 * d0 * m_a + d1 * d1 * m_d + 2.0 * d0 * d1 * m_b;

  
  const double tolerance =
    radius_squared <= 1.0e3 ? 1.0e-12 :
    radius_squared <= 1.0e5 ? 1.0e-11 : 1.0e-9;
  return dist_sq < radius_squared - tolerance;
}

// Isotropic in-circumcircle test (original).
[[nodiscard]] bool point_in_circumcircle(
  const std::array<double, 2> &a,
  const std::array<double, 2> &b,
  const std::array<double, 2> &c,
  const std::array<double, 2> &p
) noexcept
{
  const double ax = a[0] - p[0];
  const double ay = a[1] - p[1];
  const double bx = b[0] - p[0];
  const double by = b[1] - p[1];
  const double cx = c[0] - p[0];
  const double cy = c[1] - p[1];

  const double determinant =
    (ax * ax + ay * ay) * (bx * cy - by * cx) -
    (bx * bx + by * by) * (ax * cy - ay * cx) +
    (cx * cx + cy * cy) * (ax * by - ay * bx);
  return determinant > kUvAreaTolerance;
}

[[nodiscard]] unsigned triangle_source_face_penalty(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  const FaceMeshPlan &plan,
  const std::array<std::uint32_t, 3> &candidate,
  testing::AutoCfdSurfaceSourceFaceFailureCounts *counts
) noexcept
{
  LocalTriangleSourceFaceFailureKind failure =
    LocalTriangleSourceFaceFailureKind::none;
  if(triangle_passes_source_face_screen(
       face_view,
       face_preprocess,
       plan,
       candidate,
       &failure
     )) {
    return 0U;
  }
  if(counts != nullptr) {
    record_source_face_failure(*counts, failure);
  }

  switch(failure) {
  case LocalTriangleSourceFaceFailureKind::source_face_containment_failure:
    return 1U;
  case LocalTriangleSourceFaceFailureKind::orientation_flip:
    return 2U;
  case LocalTriangleSourceFaceFailureKind::node_reprojection_failure:
  default:
    return 3U;
  }
}

void enqueue_front_layer_candidates(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceSizingFieldState &field,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  const FaceMeshPlan &plan,
  const std::vector<std::array<std::uint32_t, 3>> &neighbors,
  const std::unordered_set<LocalEdgeKey, LocalEdgeKeyHash> &front_edges,
  std::uint32_t begin_triangle_index,
  std::uint32_t end_triangle_index,
  std::unordered_set<std::uint32_t> &queued_triangles,
  std::vector<FrontEdgeCandidate> &candidates
)
{
  const auto triangle_count =
    static_cast<std::uint32_t>(plan.triangles.size());
  const auto clamped_end = std::min(end_triangle_index, triangle_count);
  for(std::uint32_t triangle_index = begin_triangle_index;
      triangle_index < clamped_end;
      ++triangle_index) {
    const auto &triangle = plan.triangles[triangle_index];
    if(!triangle.valid ||
       queued_triangles.find(triangle_index) != queued_triangles.end()) {
      continue;
    }

    const auto active_edge =
      find_active_edge(plan, neighbors, triangle_index, &front_edges);
    if(active_edge >= 3U) {
      continue;
    }

    const double radius = effective_triangle_front_radius(triangle);
    if(!std::isfinite(radius)) {
      continue;
    }

    queued_triangles.insert(triangle_index);
    candidates.push_back({triangle_index, active_edge, radius});
  }
}

[[nodiscard]] bool edge_is_present_in_valid_mesh(
  const FaceMeshPlan &plan,
  LocalEdgeKey edge
) noexcept
{
  for(const auto &triangle : plan.triangles) {
    if(!triangle.valid) {
      continue;
    }
    if(triangle_contains_edge(triangle, edge)) {
      return true;
    }
  }

  return false;
}

void collect_missing_constrained_edges(
  const FaceMeshPlan &plan,
  std::vector<LocalEdgeKey> &missing_edges
)
{
  missing_edges.clear();
  missing_edges.reserve(plan.constrained_edges.size());

  for(const auto edge : plan.constrained_edges) {
    if(edge.first >= plan.vertices.size() ||
       edge.second >= plan.vertices.size()) {
      continue;
    }
    if(edge_is_present_in_valid_mesh(plan, edge)) {
      continue;
    }
    missing_edges.push_back(edge);
  }

  std::sort(
    missing_edges.begin(),
    missing_edges.end(),
    [](LocalEdgeKey lhs, LocalEdgeKey rhs) noexcept {
      if(lhs.first != rhs.first) {
        return lhs.first < rhs.first;
      }
      return lhs.second < rhs.second;
    }
  );
}

[[nodiscard]] bool extract_seed_loops_from_boundary_graph(
  const FaceMeshPlan &plan,
  std::vector<SeedLoopData> &loops
)
{
  loops.clear();
  if(plan.constrained_edges.empty()) {
    return false;
  }

  std::unordered_map<std::uint32_t, std::vector<std::uint32_t>> adjacency;
  adjacency.reserve(plan.constrained_edges.size() * 2U);
  for(const auto edge : plan.constrained_edges) {
    if(edge.first >= plan.vertices.size() || edge.second >= plan.vertices.size()) {
      return false;
    }
    adjacency[edge.first].push_back(edge.second);
    adjacency[edge.second].push_back(edge.first);
  }

  for(const auto &[vertex, neighbors] : adjacency) {
    static_cast<void>(vertex);
    if(neighbors.size() != 2U) {
      return false;
    }
  }

  std::unordered_set<LocalEdgeKey, LocalEdgeKeyHash> visited_edges;
  visited_edges.reserve(plan.constrained_edges.size() * 2U);
  for(const auto edge : plan.constrained_edges) {
    if(visited_edges.find(edge) != visited_edges.end()) {
      continue;
    }

    SeedLoopData loop;
    std::uint32_t start = edge.first;
    std::uint32_t previous = invalid_index;
    std::uint32_t current = start;
    while(true) {
      if(current >= plan.vertices.size()) {
        return false;
      }

      loop.vertices.push_back(current);
      loop.uvs.push_back(plan.vertices[current].uv);

      const auto adj_it = adjacency.find(current);
      if(adj_it == adjacency.end() || adj_it->second.size() != 2U) {
        return false;
      }

      const auto &neighbors = adj_it->second;
      const std::uint32_t next =
        neighbors[0] == previous ? neighbors[1] : neighbors[0];
      const auto traversed = canonical_edge(current, next);
      visited_edges.insert(traversed);
      const auto owner_it = plan.constrained_edge_owners.find(traversed);
      loop.edge_owners.push_back(
        owner_it != plan.constrained_edge_owners.end()
          ? owner_it->second
          : geo::TopologyEntityId {}
      );
      previous = current;
      current = next;
      if(current == start) {
        break;
      }
      if(loop.vertices.size() > plan.constrained_edges.size()) {
        return false;
      }
    }

    if(loop.vertices.size() < 3U) {
      return false;
    }
    loops.push_back(std::move(loop));
  }

  if(loops.empty()) {
    return false;
  }

  std::size_t outer_index = 0U;
  double outer_area = -1.0;
  for(std::size_t index = 0U; index < loops.size(); ++index) {
    const auto &loop = loops[index];
    double area = 0.0;
    for(std::size_t i = 0U; i < loop.uvs.size(); ++i) {
      const auto &a = loop.uvs[i];
      const auto &b = loop.uvs[(i + 1U) % loop.uvs.size()];
      area += a[0] * b[1] - b[0] * a[1];
    }
    area = std::abs(area);
    if(area > outer_area) {
      outer_area = area;
      outer_index = index;
    }
  }

  for(std::size_t index = 0U; index < loops.size(); ++index) {
    loops[index].kind =
      index == outer_index
        ? geo::FaceBoundaryLoopKind::outer
        : geo::FaceBoundaryLoopKind::inner;
  }
  return true;
}

[[nodiscard]] bool densify_unrecovered_boundary_edges_and_reseed(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  FaceMeshPlan &plan,
  const std::vector<LocalEdgeKey> &missing_edges
)
{
  std::vector<SeedLoopData> loops;
  if(!extract_seed_loops_from_boundary_graph(plan, loops)) {
    return false;
  }

  std::unordered_set<LocalEdgeKey, LocalEdgeKeyHash> edges_to_split(
    missing_edges.begin(),
    missing_edges.end()
  );
  bool split_any_edge = false;
  for(auto &loop : loops) {
    SeedLoopData densified_loop;
    densified_loop.kind = loop.kind;
    densified_loop.vertices.reserve(loop.vertices.size() + edges_to_split.size());
    densified_loop.uvs.reserve(loop.uvs.size() + edges_to_split.size());
    densified_loop.edge_owners.reserve(
      loop.edge_owners.size() + edges_to_split.size() * 2U
    );
    for(std::size_t index = 0U; index < loop.vertices.size(); ++index) {
      const auto current = loop.vertices[index];
      const auto next = loop.vertices[(index + 1U) % loop.vertices.size()];
      const auto owner =
        index < loop.edge_owners.size() ? loop.edge_owners[index] : geo::TopologyEntityId {};
      densified_loop.vertices.push_back(current);
      densified_loop.uvs.push_back(loop.uvs[index]);

      if(edges_to_split.find(canonical_edge(current, next)) == edges_to_split.end()) {
        densified_loop.edge_owners.push_back(owner);
        continue;
      }

      const auto &current_uv = plan.vertices[current].uv;
      const auto &next_uv = plan.vertices[next].uv;
      const std::array<double, 2> midpoint_uv {
        0.5 * (current_uv[0] + next_uv[0]),
        0.5 * (current_uv[1] + next_uv[1]),
      };

      geo::Point3 midpoint_position {
        0.5 * (plan.vertices[current].position[0] + plan.vertices[next].position[0]),
        0.5 * (plan.vertices[current].position[1] + plan.vertices[next].position[1]),
        0.5 * (plan.vertices[current].position[2] + plan.vertices[next].position[2]),
      };
      if(plan.planar_fallback) {
        midpoint_position = planar_world_point(plan, midpoint_uv);
      }
      else {
        geo::FaceSample sample;
        if(!sample_face_at_material_uv(
             face_view,
             face_preprocess,
             midpoint_uv,
             sample
           )) {
          return false;
        }
        midpoint_position = sample.position;
      }

      const auto midpoint_index =
        static_cast<std::uint32_t>(plan.vertices.size());
      plan.vertices.push_back(
        {
          invalid_index,
          midpoint_position,
          midpoint_uv,
          false,
          0.0,
          0.0,
          owner,
        }
      );
      densified_loop.edge_owners.push_back(owner);
      densified_loop.vertices.push_back(midpoint_index);
      densified_loop.uvs.push_back(midpoint_uv);
      densified_loop.edge_owners.push_back(owner);
      split_any_edge = true;
    }
    loop = std::move(densified_loop);
  }

  if(!split_any_edge) {
    return false;
  }

  std::size_t outer_loop_index = invalid_index;
  for(std::size_t index = 0U; index < loops.size(); ++index) {
    if(loops[index].kind == geo::FaceBoundaryLoopKind::outer) {
      outer_loop_index = index;
      break;
    }
  }
  if(outer_loop_index == invalid_index) {
    return false;
  }

  std::vector<std::array<std::uint32_t, 3>> triangles;
  SeedLoopData outer_for_cdt = loops[outer_loop_index];
  std::vector<SeedLoopData> inner_loops;
  inner_loops.reserve(loops.size());
  for(std::size_t index = 0U; index < loops.size(); ++index) {
    if(index == outer_loop_index) {
      continue;
    }
    inner_loops.push_back(loops[index]);
  }

  enforce_outer_inner_loop_winding(outer_for_cdt, inner_loops);

  if(!build_boundary_cdt_with_holes(
       plan.vertices,
       outer_for_cdt,
       inner_loops,
       triangles)) {
    return false;
  }

  plan.constrained_edges.clear();
  plan.constrained_edge_owners.clear();
  for(const auto &loop : loops) {
    add_constrained_loop_edges(
      loop,
      plan.constrained_edges,
      &plan.constrained_edge_owners
    );
  }

  plan.triangles.clear();
  plan.triangles.reserve(triangles.size());
  for(auto triangle : triangles) {
    const double area = signed_triangle_area_twice(
      plan.vertices[triangle[0]].uv,
      plan.vertices[triangle[1]].uv,
      plan.vertices[triangle[2]].uv
    );
    if(area < 0.0) {
      std::swap(triangle[1], triangle[2]);
    }
    plan.triangles.push_back({triangle, true, false});
  }
  return !plan.triangles.empty();
}

struct ConstrainedRecoverStats final {
  std::size_t calls = 0U;
  std::size_t target_invalid = 0U;
  std::size_t target_already_present = 0U;
  std::size_t empty_intersected = 0U;
  std::size_t no_swap_choice_succeeded = 0U;
  std::size_t swap_iter_limit_hit = 0U;
  std::size_t intersected_total = 0U;
  std::size_t constrained_blocking_total = 0U;
  std::size_t reject_re_intersect = 0U;
  std::size_t reject_neighbor_slot = 0U;
  std::size_t reject_same_opposite = 0U;
  std::size_t reject_new_diag_intersects = 0U;
  std::size_t reject_non_convex = 0U;
  std::size_t reject_orient = 0U;
  std::size_t swapped = 0U;
};

static thread_local ConstrainedRecoverStats g_recover_stats {};

// Recover one constrained edge by repeatedly swapping mesh edges that cross
// the target segment in UV space. Standard CDT edge-recovery flip pass
// (Lawson-style edge swaps applied along the missing constraint; see
// Shewchuk, "Constrained Delaunay Triangulations and Robust Geometric
// Predicates", and Chew 1989). Each outer iteration recollects the
// intersected list (topology shifts after every swap), tries swap
// candidates until one succeeds, and exits when the target edge is present
// (intersected list empty) or no candidate can be swapped (fatal). The
// previous version collected `intersected` once and returned after a
// single swap, leaving outer callers to retry without state continuity —
// that wasted swaps without ever recovering CRM's face:0.
[[nodiscard]] bool try_recover_constrained_edge_by_swap(
  FaceMeshPlan &plan,
  LocalEdgeKey target_edge
)
{
  ++g_recover_stats.calls;
  if(target_edge.first >= plan.vertices.size() ||
     target_edge.second >= plan.vertices.size()) {
    ++g_recover_stats.target_invalid;
    return false;
  }

  if(edge_is_present_in_valid_mesh(plan, target_edge)) {
    ++g_recover_stats.target_already_present;
    return true;
  }

  const auto &target_a = plan.vertices[target_edge.first].uv;
  const auto &target_b = plan.vertices[target_edge.second].uv;

  constexpr std::size_t kMaxSwapIterations = 300U;

  struct IntersectingEdge {
    std::uint32_t triangle_index;
    std::uint8_t edge_slot;
  };

  std::vector<std::array<std::uint32_t, 3>> neighbors;
  std::vector<IntersectingEdge> intersected;

  for(std::size_t iter = 0U; iter < kMaxSwapIterations; ++iter) {
    neighbors.clear();
    build_triangle_neighbors(plan.triangles, neighbors);

    intersected.clear();
    std::size_t constrained_intersect_count = 0U;

    for(std::uint32_t triangle_index = 0U;
        triangle_index < plan.triangles.size();
        ++triangle_index) {
      const auto &triangle = plan.triangles[triangle_index];
      if(!triangle.valid) continue;

      for(std::uint8_t edge_slot = 0U; edge_slot < 3U; ++edge_slot) {
        const auto neighbor_index = neighbors[triangle_index][edge_slot];
        if(neighbor_index == invalid_index) continue;
        if(!plan.triangles[neighbor_index].valid) continue;

        const auto shared_a = triangle.vertices[edge_slot];
        const auto shared_b = triangle.vertices[(edge_slot + 1U) % 3U];

        // Deduplicate: only process each undirected edge once.
        if(shared_a > shared_b) continue;

        if(shared_a == target_edge.first || shared_a == target_edge.second ||
           shared_b == target_edge.first || shared_b == target_edge.second) {
          continue;
        }

        if(!uv_segments_intersect(
             plan.vertices[shared_a].uv, plan.vertices[shared_b].uv,
             target_a, target_b)) {
          continue;
        }

        const auto shared_edge = canonical_edge(shared_a, shared_b);
        if(plan.constrained_edges.find(shared_edge) != plan.constrained_edges.end()) {
          ++constrained_intersect_count;
          continue;
        }

        intersected.push_back({triangle_index, edge_slot});
      }
    }

    g_recover_stats.constrained_blocking_total += constrained_intersect_count;
    g_recover_stats.intersected_total += intersected.size();

    if(intersected.empty()) {
      // Either target is now in the mesh, or it shares a face with no
      // crossing edge — verify and report accordingly.
      if(edge_is_present_in_valid_mesh(plan, target_edge)) {
        return true;
      }
      ++g_recover_stats.empty_intersected;
      return false;
    }

    bool swapped = false;
    for(const auto &ie : intersected) {
      auto &triangle = plan.triangles[ie.triangle_index];
      if(!triangle.valid) continue;

      const auto neighbor_index = neighbors[ie.triangle_index][ie.edge_slot];
      if(neighbor_index == invalid_index) continue;
      auto &neighbor = plan.triangles[neighbor_index];
      if(!neighbor.valid) continue;

      const auto shared_a = triangle.vertices[ie.edge_slot];
      const auto shared_b = triangle.vertices[(ie.edge_slot + 1U) % 3U];
      const auto shared_edge = canonical_edge(shared_a, shared_b);

      if(!uv_segments_intersect(
           plan.vertices[shared_a].uv, plan.vertices[shared_b].uv,
           target_a, target_b)) {
        ++g_recover_stats.reject_re_intersect;
        continue;
      }

      const auto neighbor_edge_slot = find_triangle_edge_slot(neighbor, shared_edge);
      if(neighbor_edge_slot >= 3U) {
        ++g_recover_stats.reject_neighbor_slot;
        continue;
      }

      const auto opposite_c = triangle.vertices[(ie.edge_slot + 2U) % 3U];
      const auto opposite_d = neighbor.vertices[(neighbor_edge_slot + 2U) % 3U];
      if(opposite_c == opposite_d) {
        ++g_recover_stats.reject_same_opposite;
        continue;
      }

      const auto &opposite_c_uv = plan.vertices[opposite_c].uv;
      const auto &opposite_d_uv = plan.vertices[opposite_d].uv;
      const bool new_diag_touches_target =
        opposite_c == target_edge.first || opposite_c == target_edge.second ||
        opposite_d == target_edge.first || opposite_d == target_edge.second;
      if(!new_diag_touches_target &&
         uv_segments_intersect(opposite_c_uv, opposite_d_uv, target_a, target_b)) {
        ++g_recover_stats.reject_new_diag_intersects;
        continue;
      }

      const auto &p1_uv = plan.vertices[shared_a].uv;
      const auto &p2_uv = plan.vertices[shared_b].uv;
      const auto &op1_uv = plan.vertices[opposite_c].uv;
      const auto &op2_uv = plan.vertices[opposite_d].uv;

      const double ori_t1 = signed_triangle_area_twice(op1_uv, p1_uv, op2_uv);
      const double ori_t2 = signed_triangle_area_twice(op1_uv, op2_uv, p2_uv);
      if(ori_t1 * ori_t2 <= 0.0) {
        ++g_recover_stats.reject_non_convex;
        continue;
      }

      std::array<std::uint32_t, 3> first_triangle {shared_b, opposite_c, opposite_d};
      std::array<std::uint32_t, 3> second_triangle {opposite_d, opposite_c, shared_a};
      if(!orient_triangle_ccw(plan.vertices, first_triangle) ||
         !orient_triangle_ccw(plan.vertices, second_triangle)) {
        ++g_recover_stats.reject_orient;
        continue;
      }

      triangle.vertices = first_triangle;
      neighbor.vertices = second_triangle;
      ++plan.instability.repair_flip_count;
      ++g_recover_stats.swapped;
      swapped = true;
      break;
    }

    if(!swapped) {
      ++g_recover_stats.no_swap_choice_succeeded;
      return false;
    }
  }

  ++g_recover_stats.swap_iter_limit_hit;
  return edge_is_present_in_valid_mesh(plan, target_edge);
}

void recover_missing_constrained_edges(
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  FaceMeshPlan &plan
)
{
  std::vector<LocalEdgeKey> missing_edges;
  collect_missing_constrained_edges(plan, missing_edges);
  if(missing_edges.empty()) {
    return;
  }

  log_constrained_edge_owner_breakdown(
    face_preprocess,
    plan,
    missing_edges,
    "missing constrained-edge owner-breakdown"
  );

  g_recover_stats = {};

  std::size_t recovered_edge_count = 0U;
  std::size_t stalled_edge_count = 0U;
  for(std::size_t edge_index = 0U; edge_index < missing_edges.size(); ++edge_index) {
    const auto edge = missing_edges[edge_index];
    if(try_recover_constrained_edge_by_swap(plan, edge)) {
      ++recovered_edge_count;
    }
    else {
      ++stalled_edge_count;
    }
  }

  SQMESH_LOG_WARN(
    "constrained-edge swap recovery on {}: missing_edges={}, recovered={}, "
    "stalled={}, swap_calls={}, swapped={}, target_already_present={}, "
    "empty_intersected={}, no_swap_choice={}, swap_iter_limit_hit={}, "
    "intersected_total={}, constrained_blocking={}, reject_re_intersect={}, "
    "reject_neighbor_slot={}, reject_same_opposite={}, "
    "reject_new_diag_intersects={}, reject_non_convex={}, reject_orient={}",
    topology_entity_debug_label(face_preprocess.face),
    missing_edges.size(),
    recovered_edge_count,
    stalled_edge_count,
    g_recover_stats.calls,
    g_recover_stats.swapped,
    g_recover_stats.target_already_present,
    g_recover_stats.empty_intersected,
    g_recover_stats.no_swap_choice_succeeded,
    g_recover_stats.swap_iter_limit_hit,
    g_recover_stats.intersected_total,
    g_recover_stats.constrained_blocking_total,
    g_recover_stats.reject_re_intersect,
    g_recover_stats.reject_neighbor_slot,
    g_recover_stats.reject_same_opposite,
    g_recover_stats.reject_new_diag_intersects,
    g_recover_stats.reject_non_convex,
    g_recover_stats.reject_orient);
}

struct DelaunayRepairPassStats final {
  std::size_t executed_passes = 0U;
  std::size_t total_flips = 0U;
  bool hit_pass_limit = false;
};

[[nodiscard]] double triangle_min_angle_3d(
  const FaceMeshPlan &plan,
  const std::array<std::uint32_t, 3> &triangle
) noexcept
{
  if(triangle[0] >= plan.vertices.size() ||
     triangle[1] >= plan.vertices.size() ||
     triangle[2] >= plan.vertices.size()) {
    return 0.0;
  }
  return triangle_min_angle_3d(
    plan.vertices[triangle[0]].position,
    plan.vertices[triangle[1]].position,
    plan.vertices[triangle[2]].position
  );
}

[[nodiscard]] bool replacement_triangles_preserve_area(
  const FaceMeshPlan &plan,
  const LocalTriangle &first_before,
  const LocalTriangle &second_before,
  const std::array<std::uint32_t, 3> &first_after,
  const std::array<std::uint32_t, 3> &second_after
) noexcept
{
  const double before_area =
    std::abs(signed_triangle_area_twice(
      plan.vertices[first_before.vertices[0]].uv,
      plan.vertices[first_before.vertices[1]].uv,
      plan.vertices[first_before.vertices[2]].uv
    )) +
    std::abs(signed_triangle_area_twice(
      plan.vertices[second_before.vertices[0]].uv,
      plan.vertices[second_before.vertices[1]].uv,
      plan.vertices[second_before.vertices[2]].uv
    ));
  const double after_area =
    std::abs(signed_triangle_area_twice(
      plan.vertices[first_after[0]].uv,
      plan.vertices[first_after[1]].uv,
      plan.vertices[first_after[2]].uv
    )) +
    std::abs(signed_triangle_area_twice(
      plan.vertices[second_after[0]].uv,
      plan.vertices[second_after[1]].uv,
      plan.vertices[second_after[2]].uv
    ));
  return before_area > kUvAreaTolerance &&
         std::abs(before_area - after_area) / before_area <= 1.0e-8;
}

[[nodiscard]] bool replacement_flip_is_convex(
  const FaceMeshPlan &plan,
  std::uint32_t shared_a,
  std::uint32_t shared_b,
  std::uint32_t opposite_c,
  std::uint32_t opposite_d
) noexcept
{
  if(shared_a >= plan.vertices.size() ||
     shared_b >= plan.vertices.size() ||
     opposite_c >= plan.vertices.size() ||
     opposite_d >= plan.vertices.size()) {
    return false;
  }

  const auto diagonal = uv_subtract(
    plan.vertices[opposite_d].uv,
    plan.vertices[opposite_c].uv
  );
  const auto to_shared_a = uv_subtract(
    plan.vertices[shared_a].uv,
    plan.vertices[opposite_c].uv
  );
  const auto to_shared_b = uv_subtract(
    plan.vertices[shared_b].uv,
    plan.vertices[opposite_c].uv
  );
  const double cross_shared_a = uv_cross(diagonal, to_shared_a);
  const double cross_shared_b = uv_cross(diagonal, to_shared_b);
  return cross_shared_a * cross_shared_b < 0.0;
}

[[nodiscard]] bool replacement_flip_preserves_quality(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceSizingFieldState &field,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  const FaceMeshPlan &plan,
  const std::array<std::uint32_t, 3> &first_before,
  const std::array<std::uint32_t, 3> &second_before,
  const std::array<std::uint32_t, 3> &first_after,
  const std::array<std::uint32_t, 3> &second_after
)
{
  const auto metric_first_before =
    evaluate_triangle_metric_quality(face_view, field, face_preprocess, plan, first_before);
  const auto metric_second_before =
    evaluate_triangle_metric_quality(face_view, field, face_preprocess, plan, second_before);
  const auto metric_first_after =
    evaluate_triangle_metric_quality(face_view, field, face_preprocess, plan, first_after);
  const auto metric_second_after =
    evaluate_triangle_metric_quality(face_view, field, face_preprocess, plan, second_after);

  if(metric_first_before.valid && metric_second_before.valid &&
     metric_first_after.valid && metric_second_after.valid) {
    const bool before_good =
      triangle_is_good_enough(metric_first_before) &&
      triangle_is_good_enough(metric_second_before);
    const bool after_good =
      triangle_is_good_enough(metric_first_after) &&
      triangle_is_good_enough(metric_second_after);
    const double before_badness =
      triangle_metric_badness(metric_first_before) +
      triangle_metric_badness(metric_second_before);
    const double after_badness =
      triangle_metric_badness(metric_first_after) +
      triangle_metric_badness(metric_second_after);
    if(after_good && !before_good) {
      return true;
    }
    return after_badness + 1.0e-6 < before_badness;
  }

  const double before_min_angle = std::min(
    triangle_min_angle_3d(plan, first_before),
    triangle_min_angle_3d(plan, second_before)
  );
  const double after_min_angle = std::min(
    triangle_min_angle_3d(plan, first_after),
    triangle_min_angle_3d(plan, second_after)
  );
  return after_min_angle > before_min_angle + 1.0e-3;
}

[[nodiscard]] DelaunayRepairPassStats run_delaunay_edge_flip_passes(
    const geo::FaceView& face_view,
    const AutoCfdSurfaceSizingFieldState& field,
    const AutoCfdSurfaceFacePreprocessState& face_preprocess,
    FaceMeshPlan& plan,
    std::size_t max_passes
)
{
  
  DelaunayRepairPassStats stats;
  std::vector<std::array<std::uint32_t, 3>> neighbors;
  std::vector<bool> dirty_triangles;

  for(std::size_t pass = 0U; pass < max_passes; ++pass) {
    build_triangle_neighbors(plan.triangles, neighbors);
    dirty_triangles.assign(plan.triangles.size(), false);

    std::size_t flip_count = 0U;
    for(std::uint32_t triangle_index = 0U;
        triangle_index < plan.triangles.size();
        ++triangle_index) {
      if(triangle_index >= plan.triangles.size() ||
         !plan.triangles[triangle_index].valid ||
         dirty_triangles[triangle_index]) {
        continue;
      }

      for(std::uint8_t edge_slot = 0U; edge_slot < 3U; ++edge_slot) {
        const auto neighbor_index = neighbors[triangle_index][edge_slot];
        if(neighbor_index == invalid_index ||
           neighbor_index <= triangle_index ||
           neighbor_index >= plan.triangles.size() ||
           !plan.triangles[neighbor_index].valid ||
           dirty_triangles[neighbor_index]) {
          continue;
        }

        auto &triangle = plan.triangles[triangle_index];
        auto &neighbor = plan.triangles[neighbor_index];
        const auto shared_a = triangle.vertices[edge_slot];
        const auto shared_b = triangle.vertices[(edge_slot + 1U) % 3U];
        const auto shared_edge = canonical_edge(shared_a, shared_b);
        if(plan.constrained_edges.find(shared_edge) != plan.constrained_edges.end()) {
          continue;
        }

        const auto opposite_c = triangle.vertices[(edge_slot + 2U) % 3U];
        const auto neighbor_edge_slot = find_triangle_edge_slot(neighbor, shared_edge);
        if(neighbor_edge_slot >= 3U) {
          continue;
        }

        const auto opposite_d = neighbor.vertices[(neighbor_edge_slot + 2U) % 3U];
        if(opposite_c == opposite_d ||
           !replacement_flip_is_convex(
             plan, shared_a, shared_b, opposite_c, opposite_d
           )) {
          continue;
        }

        const std::array<double, 2> tri_centroid_uv {
          (plan.vertices[shared_a].uv[0] + plan.vertices[shared_b].uv[0] +
           plan.vertices[opposite_c].uv[0]) / 3.0,
          (plan.vertices[shared_a].uv[1] + plan.vertices[shared_b].uv[1] +
           plan.vertices[opposite_c].uv[1]) / 3.0,
        };
        AutoCfdSurfaceFaceMetricTensor flip_metric;
        const bool have_metric =
          sample_metric_or_reference(
            face_view, field, face_preprocess, &plan,
            tri_centroid_uv[0], tri_centroid_uv[1], flip_metric
          ) == base::StatusCode::ok;

        bool inside = false;
        if(have_metric && flip_metric.usable &&
           flip_metric.first_fundamental_form_defined) {
          inside = point_in_circumcircle_aniso(
            plan.vertices[shared_a].uv,
            plan.vertices[shared_b].uv,
            plan.vertices[opposite_c].uv,
            plan.vertices[opposite_d].uv,
            flip_metric.first_fundamental_form
          );
        }
        else {
          inside = point_in_circumcircle(
            plan.vertices[shared_a].uv,
            plan.vertices[shared_b].uv,
            plan.vertices[opposite_c].uv,
            plan.vertices[opposite_d].uv
          );
        }
        if(!inside) {
          continue;
        }

        const auto first_before = triangle.vertices;
        const auto second_before = neighbor.vertices;
        std::array<std::uint32_t, 3> first_after {
          opposite_c, opposite_d, shared_b,
        };
        std::array<std::uint32_t, 3> second_after {
          opposite_d, opposite_c, shared_a,
        };
        if(!orient_triangle_ccw(plan.vertices, first_after) ||
           !orient_triangle_ccw(plan.vertices, second_after) ||
           !replacement_triangles_preserve_area(
             plan, triangle, neighbor, first_after, second_after
           ) ||
           !replacement_flip_preserves_quality(
             face_view,
             field,
             face_preprocess,
             plan,
             first_before,
             second_before,
             first_after,
             second_after
           )) {
          continue;
        }

        triangle.vertices = first_after;
        neighbor.vertices = second_after;
        dirty_triangles[triangle_index] = true;
        dirty_triangles[neighbor_index] = true;
        ++flip_count;
        ++stats.total_flips;
        ++plan.instability.repair_flip_count;
        break;
      }
    }

    ++stats.executed_passes;
   
    if(flip_count == 0U) {
      return stats;
    }
  }

  stats.hit_pass_limit = true;
  return stats;
}

void run_delaunay_repair(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceSizingFieldState &field,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  FaceMeshPlan &plan,
  bool recover_constrained_edges
)
{
  const auto valid_triangle_count_before =
    std::count_if(
      plan.triangles.begin(),
      plan.triangles.end(),
      [](const LocalTriangle &triangle) noexcept {
        return triangle.valid;
      }
    );
  
  if(recover_constrained_edges) {
    recover_missing_constrained_edges(face_preprocess, plan);
  }

  const auto repair_stats = run_delaunay_edge_flip_passes(
    face_view, field, face_preprocess, plan, kMaximumRepairPasses
  );
}

void accept_all_valid_triangles(FaceMeshPlan &plan) noexcept
{
  for(auto &triangle : plan.triangles) {
    if(triangle.valid) {
      triangle.accepted = true;
    }
  }
}

void clear_all_triangle_acceptance(FaceMeshPlan &plan) noexcept
{
  for(auto &triangle : plan.triangles) {
    triangle.accepted = false;
  }
}

void collect_unrecovered_boundary_retry_requests(
  const FaceMeshPlan &plan,
  const std::vector<LocalEdgeKey> &missing_edges,
  std::vector<AutoCfdSurfaceBoundaryRetryRequest> &requests
)
{
  requests.clear();
  requests.reserve(missing_edges.size());
  for(const auto edge : missing_edges) {
    const auto owner_it = plan.constrained_edge_owners.find(edge);
    if(owner_it == plan.constrained_edge_owners.end() ||
       !geo::is_valid(owner_it->second)) {
      continue;
    }

    if(edge.first >= plan.vertices.size() ||
       edge.second >= plan.vertices.size()) {
      continue;
    }

    const auto first_node_index = plan.vertices[edge.first].boundary_node_index;
    const auto second_node_index = plan.vertices[edge.second].boundary_node_index;
    if(first_node_index == invalid_index ||
       second_node_index == invalid_index ||
       first_node_index == second_node_index) {
      continue;
    }

    requests.push_back(
      {
        owner_it->second,
        std::min(first_node_index, second_node_index),
        std::max(first_node_index, second_node_index),
      }
    );
  }

  std::sort(
    requests.begin(),
    requests.end(),
    [](const AutoCfdSurfaceBoundaryRetryRequest &lhs,
       const AutoCfdSurfaceBoundaryRetryRequest &rhs) noexcept {
      if(lhs.edge.dimension != rhs.edge.dimension) {
        return lhs.edge.dimension < rhs.edge.dimension;
      }
      if(lhs.edge.index != rhs.edge.index) {
        return lhs.edge.index < rhs.edge.index;
      }
      if(lhs.first_node_index != rhs.first_node_index) {
        return lhs.first_node_index < rhs.first_node_index;
      }
      return lhs.second_node_index < rhs.second_node_index;
    }
  );
  requests.erase(
    std::unique(
      requests.begin(),
      requests.end(),
      [](const AutoCfdSurfaceBoundaryRetryRequest &lhs,
         const AutoCfdSurfaceBoundaryRetryRequest &rhs) noexcept {
        return lhs.edge.dimension == rhs.edge.dimension &&
               lhs.edge.index == rhs.edge.index &&
               lhs.first_node_index == rhs.first_node_index &&
               lhs.second_node_index == rhs.second_node_index;
      }
    ),
    requests.end()
  );
}

[[nodiscard]] base::StatusCode run_initial_mesh_boundary_recovery(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceSizingFieldState &field,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  FaceMeshPlan &plan,
  std::vector<AutoCfdSurfaceBoundaryRetryRequest> *unrecovered_boundary_edges = nullptr
)
{

  accept_all_valid_triangles(plan);

  run_delaunay_repair(face_view, field, face_preprocess, plan, true);

  std::vector<LocalEdgeKey> missing_edges;
  collect_missing_constrained_edges(plan, missing_edges);
  if(missing_edges.empty()) {
    return core::detail::clear_error_state();
  }

  log_constrained_edge_owner_breakdown(
    face_preprocess,
    plan,
    missing_edges,
    "unrecovered constrained-edge owner-breakdown"
  );

  
  if(unrecovered_boundary_edges != nullptr) {
    collect_unrecovered_boundary_retry_requests(
      plan, missing_edges, *unrecovered_boundary_edges
    );
    if(!unrecovered_boundary_edges->empty()) {
      return core::detail::publish_error(
        base::StatusCode::unsupported,
        "Auto CFD Surface Mesher: boundary edges unrecovered on face " +
          std::to_string(face_preprocess.face.index) +
          "; requesting shared 1D boundary retry."
      );
    }
  }

  return core::detail::clear_error_state();
}

[[nodiscard]] base::StatusCode filter_accepted_face_triangles(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  FaceMeshPlan &plan
)
{
  std::vector<LocalTriangle> accepted_triangles;
  accepted_triangles.reserve(plan.triangles.size());
  std::size_t processed_triangle_count = 0U;
  for(const auto &triangle : plan.triangles) {
    if(!triangle.valid || !triangle.accepted) {
      continue;
    }
    ++processed_triangle_count;

    LocalTriangle screened_triangle = triangle;
    if(should_flip_triangle_winding(face_view, plan, screened_triangle.vertices)) {
      std::swap(screened_triangle.vertices[1], screened_triangle.vertices[2]);
    }
    if(std::abs(signed_triangle_area_twice(
         plan.vertices[screened_triangle.vertices[0]].uv,
         plan.vertices[screened_triangle.vertices[1]].uv,
         plan.vertices[screened_triangle.vertices[2]].uv
       )) <= kUvAreaTolerance) {
      continue;
    }
    accepted_triangles.push_back(screened_triangle);
  }

  if(accepted_triangles.empty()) {
    return core::detail::publish_error(
      base::StatusCode::unsupported,
      "Auto CFD Surface Mesher: all triangles on face " +
        std::to_string(face_preprocess.face.index) +
        " were degenerate in UV space (zero area). No valid triangles remain."
    );
  }

  plan.triangles = std::move(accepted_triangles);
  log_plan_open_edge_breakdown(face_preprocess, plan, "post_filter");
  return core::detail::clear_error_state();
}


void run_laplacian_smoothing(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceSizingFieldState &field,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  FaceMeshPlan &plan,
  std::size_t niter
)
{
  
  if(niter == 0U || plan.vertices.empty() || plan.triangles.empty()) {
    return;
  }

  // Build vertex-to-neighbor-vertices adjacency from accepted triangles.
  std::vector<std::vector<std::uint32_t>> vertex_neighbors(
    plan.vertices.size()
  );
  // Also build vertex-to-incident-triangles for quality validation.
  std::vector<std::vector<std::uint32_t>> vertex_triangles(
    plan.vertices.size()
  );
  for(std::uint32_t triangle_index = 0U;
      triangle_index < plan.triangles.size();
      ++triangle_index) {
    const auto &tri = plan.triangles[triangle_index];
    if(!tri.valid || !tri.accepted) {
      continue;
    }
    for(int k = 0; k < 3; ++k) {
      const auto v = tri.vertices[k];
      const auto vn = tri.vertices[(k + 1) % 3];
      const auto vp = tri.vertices[(k + 2) % 3];
      vertex_neighbors[v].push_back(vn);
      vertex_neighbors[v].push_back(vp);
      vertex_triangles[v].push_back(triangle_index);
    }
  }

  // Deduplicate neighbor lists.
  for(auto &nbrs : vertex_neighbors) {
    std::sort(nbrs.begin(), nbrs.end());
    nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
  }

  for(std::size_t iter = 0U; iter < niter; ++iter) {
    std::size_t moved_count = 0U;
    for(std::uint32_t vi = 0U; vi < plan.vertices.size(); ++vi) {
      auto &vertex = plan.vertices[vi];
      // Only smooth interior (non-boundary) vertices.
      if(vertex.boundary) {
        continue;
      }
      const auto &nbrs = vertex_neighbors[vi];
      if(nbrs.empty()) {
        continue;
      }

      // Compute centroid of neighbors in UV space.
      double avg_u = 0.0;
      double avg_v = 0.0;
      for(const auto ni : nbrs) {
        avg_u += plan.vertices[ni].uv[0];
        avg_v += plan.vertices[ni].uv[1];
      }
      avg_u /= static_cast<double>(nbrs.size());
      avg_v /= static_cast<double>(nbrs.size());

      // Save old position.
      const auto old_uv = vertex.uv;
      const auto old_position = vertex.position;

      // Compute new position.
      const std::array<double, 2> new_uv {avg_u, avg_v};
      geo::Point3 new_position;
      if(plan.planar_fallback) {
        new_position = planar_world_point(plan, new_uv);
      }
      else {
        geo::FaceSample sample;
        if(geo::sample_face(face_view, new_uv[0], new_uv[1], sample) !=
           base::StatusCode::ok) {
          continue;
        }
        new_position = sample.position;
      }

      // Tentatively apply the new position.
      vertex.uv = new_uv;
      vertex.position = new_position;

      // Validate: check that no incident triangle has zero or negative area.
      bool valid_move = true;
      for(const auto tri_idx : vertex_triangles[vi]) {
        const auto &tri = plan.triangles[tri_idx];
        if(!tri.valid || !tri.accepted) {
          continue;
        }
        const double area = signed_triangle_area_twice(
          plan.vertices[tri.vertices[0]].uv,
          plan.vertices[tri.vertices[1]].uv,
          plan.vertices[tri.vertices[2]].uv
        );
        if(area <= kUvAreaTolerance) {
          valid_move = false;
          break;
        }
      }

      if(!valid_move) {
        // Revert.
        vertex.uv = old_uv;
        vertex.position = old_position;
      }
      else {
        ++moved_count;
      }
    }

    if(moved_count == 0U) {
      break;
    }
  }
}

[[nodiscard]] base::StatusCode generate_face_mesh_plan(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceSizingFieldState &field,
  const AutoCfdSurfaceFacePreprocessState &face_preprocess,
  FaceMeshPlan &plan,
  std::vector<AutoCfdSurfaceBoundaryRetryRequest> *unrecovered_boundary_edges = nullptr
)
{
  auto status = build_face_seed_mesh(face_view, face_preprocess, field, plan);
  if(status != base::StatusCode::ok) {
    return status;
  }
  plan.seed_vertex_count = plan.vertices.size();
  plan.seed_triangle_count = std::count_if(
    plan.triangles.begin(),
    plan.triangles.end(),
    [](const LocalTriangle &triangle) {
      return triangle.valid;
    }
  );
  // ── Sizing field diagnostic output ──
  {
    // Query sizing field at a few boundary vertices
    for(std::size_t vi = 0U; vi < std::min<std::size_t>(plan.vertices.size(), 4U); ++vi) {
      const auto &v = plan.vertices[vi];
      const double owner_size =
        query_auto_cfd_surface_sizing_field(field, face_preprocess.face, v.position);
      const double uv_size =
        query_auto_cfd_surface_sizing_field(field, face_view, v.uv[0], v.uv[1], v.position);
      AutoCfdSurfaceFaceMetricTensor metric;
      MetricSampleTrace metric_trace;
      static_cast<void>(sample_metric_or_reference(
        face_view,
        field,
        face_preprocess,
        &plan,
        v.uv[0],
        v.uv[1],
        metric,
        &metric_trace
      ));
    }
    // Print background octree info
    // Print first fundamental form from reference metric
  }

  const bool allow_frontal_refinement =
    face_preprocess.disposition !=
      AutoCfdSurfaceFacePreprocessDisposition::fallback_only &&
    (!face_preprocess.boundary.has_seams ||
     uses_supported_seam_material_screen(face_preprocess) ||
     face_preprocess.seam_support ==
       AutoCfdSurfaceSeamSupportKind::periodic);

  auto t0 = std::chrono::steady_clock::now();

  status =
    run_initial_mesh_boundary_recovery(
      face_view,
      field,
      face_preprocess,
      plan,
      unrecovered_boundary_edges
    );
  if(status != base::StatusCode::ok) {
    return status;
  }

  auto t1 = std::chrono::steady_clock::now();

  if(allow_frontal_refinement) {
    clear_all_triangle_acceptance(plan);
    run_front_refinement(face_view, field, face_preprocess, plan);
  }
  else {
    accept_all_valid_triangles(plan);
  }

  auto t2 = std::chrono::steady_clock::now();

  log_plan_open_edge_breakdown(face_preprocess, plan, "post_bw");
  run_delaunay_repair(face_view, field, face_preprocess, plan, false);

  auto t3 = std::chrono::steady_clock::now();

  run_laplacian_smoothing(face_view, field, face_preprocess, plan, 3U);

  auto t4 = std::chrono::steady_clock::now();

  status = filter_accepted_face_triangles(face_view, face_preprocess, plan);

  auto t5 = std::chrono::steady_clock::now();

  log_metric_sampling_summary(face_preprocess, plan, "post_filter");

  const auto ms = [](auto d) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(d).count();
  };

  return status;
}

[[nodiscard]] bool should_flip_triangle_winding(
  const geo::FaceView &face_view,
  const FaceMeshPlan &plan,
  const std::array<std::uint32_t, 3> &triangle
) noexcept
{
  if(!plan.reference_normal_defined) {
    return false;
  }
  if(triangle[0] >= plan.vertices.size() ||
     triangle[1] >= plan.vertices.size() ||
     triangle[2] >= plan.vertices.size()) {
    return false;
  }

  const auto &a = plan.vertices[triangle[0]].position;
  const auto &b = plan.vertices[triangle[1]].position;
  const auto &c = plan.vertices[triangle[2]].position;
  const auto ab = point_subtract(b, a);
  const auto ac = point_subtract(c, a);
  auto triangle_normal = cross_product(ab, ac);

  const std::array<double, 2> centroid_uv {
    (plan.vertices[triangle[0]].uv[0] +
     plan.vertices[triangle[1]].uv[0] +
     plan.vertices[triangle[2]].uv[0]) / 3.0,
    (plan.vertices[triangle[0]].uv[1] +
     plan.vertices[triangle[1]].uv[1] +
     plan.vertices[triangle[2]].uv[1]) / 3.0,
  };

  geo::FaceSample sample;
  if(plan.planar_fallback) {
    geo::FaceProjection projection;
    const auto centroid_point = planar_world_point(plan, centroid_uv);
    if(geo::project_point_to_face(face_view, centroid_point, projection) ==
         base::StatusCode::ok &&
       projection.normal_defined) {
      return dot_product(triangle_normal, projection.normal) < 0.0;
    }
  }
  else if(geo::sample_face(face_view, centroid_uv[0], centroid_uv[1], sample) ==
            base::StatusCode::ok &&
          sample.normal_defined) {
    return dot_product(triangle_normal, sample.normal) < 0.0;
  }

  return dot_product(triangle_normal, plan.reference_normal) < 0.0;
}

[[nodiscard]] base::StatusCode build_face_mesh_plans(
  const geo::ModelView &model_view,
  AutoCfdSurfacePipelineState &pipeline_state,
  const MeshingRequest *fallback_request,
  const AutoCfdSurfaceParameters *fallback_parameters,
  std::vector<FaceMeshPlan> &face_meshes
)
{
  constexpr std::size_t kMaximumSharedBoundaryRecoveryRetries = 4U;
  static_cast<void>(fallback_request);
  static_cast<void>(fallback_parameters);

  for(std::size_t retry = 0U;
      retry <= kMaximumSharedBoundaryRecoveryRetries;
      ++retry) {
    face_meshes.clear();
    face_meshes.reserve(model_view.faces.size());

    bool restarted_after_shared_boundary_retry = false;
    for(const auto &face_view : model_view.faces) {
      const auto *face_preprocess =
        find_auto_cfd_surface_face_preprocess_state(pipeline_state, face_view.entity);
      if(face_preprocess == nullptr) {
        const bool in_work_items = std::any_of(
          pipeline_state.face_work_items.begin(),
          pipeline_state.face_work_items.end(),
          [&face_view](const AutoCfdSurfaceFaceWorkItem &work_item) {
            return work_item.face == face_view.entity;
          });
        if(!in_work_items) {
          continue;
        }
        return core::detail::publish_error(
          base::StatusCode::internal_error,
          "Auto CFD Surface Mesher could not resolve face preprocessing state for one of the geometry faces."
        );
      }


      FaceMeshPlan face_mesh;
      std::vector<AutoCfdSurfaceBoundaryRetryRequest> unrecovered_boundary_edges;
      const auto status = generate_face_mesh_plan(
        face_view,
        pipeline_state.sizing_field,
        *face_preprocess,
        face_mesh,
        &unrecovered_boundary_edges
      );
      if(status == base::StatusCode::ok) {
        if(face_mesh.triangles.empty() ||
           std::none_of(
             face_mesh.triangles.begin(),
             face_mesh.triangles.end(),
             [](const LocalTriangle &triangle) {
               return triangle.valid && triangle.accepted;
             })) {
          SQMESH_LOG_WARN(
            "Auto CFD Surface Mesher produced an empty plan for {} "
            "(seed_vertices={}, seed_triangles={}, total_triangles={}, "
            "inserted_vertices={}, front_iterations={}). Face will have "
            "no triangles in the final mesh.",
            topology_entity_debug_label(face_view.entity),
            face_mesh.seed_vertex_count,
            face_mesh.seed_triangle_count,
            face_mesh.triangles.size(),
            face_mesh.inserted_vertex_count,
            face_mesh.front_iteration_count);
        }
        face_meshes.push_back(std::move(face_mesh));
        continue;
      }

      // Unsupported face with no boundary-retry edges: skip with a
      // warning rather than aborting the entire pipeline. This enables
      // partial meshing on models where some faces use configurations
      // the mesher cannot handle yet (e.g. multi-loop periodic faces
      // on missile/torus geometries).
      if(status == base::StatusCode::unsupported &&
         unrecovered_boundary_edges.empty()) {
        const auto reason = sqmesh::base::last_error_message();
        SQMESH_LOG_WARN(
          "Auto CFD Surface Mesher skipping {} (no triangles will be "
          "produced for this face): {}",
          topology_entity_debug_label(face_view.entity),
          reason.empty() ? std::string("unsupported configuration") : reason);
        continue;
      }

      // Non-unsupported failure → hard-abort.
      if(status != base::StatusCode::unsupported) {
        return status;
      }

      // status == unsupported with non-empty retry edges → shared
      // boundary retry path below.

      if(retry >= kMaximumSharedBoundaryRecoveryRetries) {
        const auto reason = sqmesh::base::last_error_message();
        SQMESH_LOG_WARN(
          "Auto CFD Surface Mesher exhausted shared-boundary retries on "
          "{} (no triangles will be produced for this face): {}",
          topology_entity_debug_label(face_view.entity),
          reason.empty() ? std::string("unsupported configuration") : reason);
        continue;
      }

      {
        std::string owner_edges_str;
        geo::TopologyEntityId previous_edge {};
        bool have_previous_edge = false;
        for(std::size_t owner_index = 0U;
            owner_index < unrecovered_boundary_edges.size();
            ++owner_index) {
          const auto edge = unrecovered_boundary_edges[owner_index].edge;
          if(have_previous_edge &&
             edge.dimension == previous_edge.dimension &&
             edge.index == previous_edge.index) {
            continue;
          }
          if(owner_index > 0U) {
            owner_edges_str += ",";
          }
          owner_edges_str += topology_entity_debug_label(edge);
          previous_edge = edge;
          have_previous_edge = true;
        }
      }

      const auto retry_status = split_auto_cfd_surface_boundary_edges_for_retry(
        model_view,
        unrecovered_boundary_edges,
        pipeline_state
      );
      if(retry_status != base::StatusCode::ok) {
        return retry_status;
      }

      restarted_after_shared_boundary_retry = true;
      break;
    }

    if(!restarted_after_shared_boundary_retry) {
      return core::detail::clear_error_state();
    }
  }

  return core::detail::publish_error(
    base::StatusCode::internal_error,
    "Auto CFD Surface Mesher shared/global boundary retry loop terminated unexpectedly."
  );
}

void populate_candidate_stats(
  const std::vector<FaceMeshPlan> &face_meshes,
  testing::AutoCfdSurfaceCandidateStats &stats
)
{
  stats = {};
  stats.face_count = face_meshes.size();
  stats.inserted_vertex_cap_per_face = kMaximumInsertedVerticesPerFace;
  stats.front_iteration_cap_per_face = kMaximumFrontIterationsPerFace;
  stats.faces.reserve(face_meshes.size());

  for(const auto &plan : face_meshes) {
    testing::AutoCfdSurfaceFaceMeshingStats face_stats;
    face_stats.face = plan.face;
    face_stats.seed_vertex_count = plan.seed_vertex_count;
    face_stats.seed_triangle_count = plan.seed_triangle_count;
    face_stats.final_vertex_count = plan.vertices.size();
    face_stats.accepted_triangle_count = plan.triangles.size();
    face_stats.inserted_vertex_count = plan.inserted_vertex_count;
    face_stats.front_iteration_count = plan.front_iteration_count;
    face_stats.inserted_vertex_limit_reached = plan.inserted_vertex_limit_reached;
    face_stats.front_iteration_limit_reached = plan.front_iteration_limit_reached;
    face_stats.instability = plan.instability;

    stats.total_accepted_triangle_count += face_stats.accepted_triangle_count;
    stats.total_inserted_vertex_count += face_stats.inserted_vertex_count;
    stats.total_front_iteration_count += face_stats.front_iteration_count;
    stats.max_inserted_vertex_count =
      std::max(stats.max_inserted_vertex_count, face_stats.inserted_vertex_count);
    stats.max_front_iteration_count =
      std::max(stats.max_front_iteration_count, face_stats.front_iteration_count);
    stats.any_inserted_vertex_limit_reached =
      stats.any_inserted_vertex_limit_reached || face_stats.inserted_vertex_limit_reached;
    stats.any_front_iteration_limit_reached =
      stats.any_front_iteration_limit_reached || face_stats.front_iteration_limit_reached;
    if(face_stats.instability.seed_source_face_failures.has_failures()) {
      ++stats.faces_with_seed_source_face_instability;
    }
    if(face_stats.instability.front_candidate_source_face_failures.has_failures()) {
      ++stats.faces_with_front_source_face_instability;
    }
    if(face_stats.instability.repair_candidate_source_face_failures.has_failures()) {
      ++stats.faces_with_repair_source_face_instability;
    }
    if(face_stats.instability.final_filter_source_face_failures.has_failures()) {
      ++stats.faces_with_final_filter_source_face_instability;
    }
    if(face_stats.instability.has_failures()) {
      ++stats.faces_with_any_meshing_instability;
    }
    stats.instability.seed_source_face_failures.node_reprojection_failure_count +=
      face_stats.instability.seed_source_face_failures.node_reprojection_failure_count;
    stats.instability.seed_source_face_failures.source_face_containment_failure_count +=
      face_stats.instability.seed_source_face_failures.source_face_containment_failure_count;
    stats.instability.seed_source_face_failures.orientation_flip_count +=
      face_stats.instability.seed_source_face_failures.orientation_flip_count;
    stats.instability.front_candidate_source_face_failures.node_reprojection_failure_count +=
      face_stats.instability.front_candidate_source_face_failures.node_reprojection_failure_count;
    stats.instability.front_candidate_source_face_failures.source_face_containment_failure_count +=
      face_stats.instability.front_candidate_source_face_failures.source_face_containment_failure_count;
    stats.instability.front_candidate_source_face_failures.orientation_flip_count +=
      face_stats.instability.front_candidate_source_face_failures.orientation_flip_count;
    stats.instability.repair_candidate_source_face_failures.node_reprojection_failure_count +=
      face_stats.instability.repair_candidate_source_face_failures.node_reprojection_failure_count;
    stats.instability.repair_candidate_source_face_failures.source_face_containment_failure_count +=
      face_stats.instability.repair_candidate_source_face_failures.source_face_containment_failure_count;
    stats.instability.repair_candidate_source_face_failures.orientation_flip_count +=
      face_stats.instability.repair_candidate_source_face_failures.orientation_flip_count;
    stats.instability.final_filter_source_face_failures.node_reprojection_failure_count +=
      face_stats.instability.final_filter_source_face_failures.node_reprojection_failure_count;
    stats.instability.final_filter_source_face_failures.source_face_containment_failure_count +=
      face_stats.instability.final_filter_source_face_failures.source_face_containment_failure_count;
    stats.instability.final_filter_source_face_failures.orientation_flip_count +=
      face_stats.instability.final_filter_source_face_failures.orientation_flip_count;
    stats.instability.repair_flip_count += face_stats.instability.repair_flip_count;
    stats.faces.push_back(face_stats);
  }
}

[[nodiscard]] base::StatusCode build_auto_cfd_surface_candidate_domain(
  const geo::ModelView &model_view,
  AutoCfdSurfacePipelineState &pipeline_state,
  const MeshingRequest *fallback_request,
  const AutoCfdSurfaceParameters *fallback_parameters,
  Domain &candidate_domain
)
{
  std::vector<FaceMeshPlan> face_meshes;
  auto status = build_face_mesh_plans(
    model_view,
    pipeline_state,
    fallback_request,
    fallback_parameters,
    face_meshes
  );
  if(status != base::StatusCode::ok) {
    return status;
  }

  return inject_face_meshes(model_view, pipeline_state, face_meshes, candidate_domain);
}

[[nodiscard]] base::StatusCode inject_face_meshes(
  const geo::ModelView &model_view,
  const AutoCfdSurfacePipelineState &pipeline_state,
  const std::vector<FaceMeshPlan> &face_meshes,
  Domain &output
)
{
  std::size_t estimated_node_upper_bound =
    pipeline_state.boundary.nodes.size();
  std::size_t estimated_triangle_count = 0U;
  for(const auto &face_mesh : face_meshes) {
    estimated_node_upper_bound += face_mesh.vertices.size();
    estimated_triangle_count += face_mesh.triangles.size();
  }

  const auto node_entity_group = output.create_entity_group({EntityOrder::node, "surface_nodes"});
  output.reserve_entity_group_storage(node_entity_group, estimated_node_upper_bound);

  // Create one face EntityGroup per geometric face (PID).
  std::unordered_map<std::uint32_t, EntityGroupIndex> geo_face_entity_groups;
  for(const auto &fm : face_meshes) {
    EntityGroupDefinition def;
    def.order = EntityOrder::face;
    def.name = "face_" + std::to_string(fm.face.index);
    def.zone_id = fm.face.index;
    def.source_entity_tag = fm.face.index;
    def.boundary = true;
    def.default_kind = EntityKind::face_triangle;
    def.semantic = EntityGroupSemantic::boundary;
    const auto tid = output.create_entity_group(def);
    output.reserve_entity_group_storage(tid, fm.triangles.size(), fm.triangles.size() * 3U);
    geo_face_entity_groups[fm.face.index] = tid;
  }

  std::unordered_map<std::uint32_t, EntityRef> boundary_node_refs;
  boundary_node_refs.reserve(pipeline_state.boundary.nodes.size());

  std::unordered_map<MeshEdgeKey, MeshEdgeIncidence, MeshEdgeKeyHash> mesh_edges;
  mesh_edges.reserve(estimated_triangle_count * 2U);

  for(const auto &face_mesh : face_meshes) {
    const auto *face_view = model_view.find_face(face_mesh.face);
    if(face_view == nullptr) {
      return core::detail::publish_error(
        base::StatusCode::internal_error,
        "Auto CFD Surface Mesher could not resolve a generated face back to the geometry view."
      );
    }

    std::size_t emitted_triangle_count = 0U;
    std::size_t collapsed_triangle_count = 0U;

    std::vector<EntityRef> local_node_refs(
      face_mesh.vertices.size(),
      EntityRef {}
    );
    for(std::size_t vertex_index = 0U; vertex_index < face_mesh.vertices.size(); ++vertex_index) {
      const auto &vertex = face_mesh.vertices[vertex_index];
      if(vertex.boundary && vertex.boundary_node_index != invalid_index) {
        auto [it, inserted] = boundary_node_refs.emplace(vertex.boundary_node_index, EntityRef {});
        if(inserted) {
          if(vertex.boundary_node_index >= pipeline_state.boundary.nodes.size()) {
            return core::detail::publish_error(
              base::StatusCode::internal_error,
              "Auto CFD Surface Mesher encountered an out-of-range boundary node while injecting the mesh."
            );
          }
          it->second = output.add_node(
            node_entity_group,
            pipeline_state.boundary.nodes[vertex.boundary_node_index].position
          );
        }
        local_node_refs[vertex_index] = it->second;
        continue;
      }

      local_node_refs[vertex_index] = output.add_node(node_entity_group, vertex.position);
    }

    for(const auto &triangle : face_mesh.triangles) {
      if(triangle.vertices[0] >= local_node_refs.size() ||
         triangle.vertices[1] >= local_node_refs.size() ||
         triangle.vertices[2] >= local_node_refs.size()) {
        return core::detail::publish_error(
          base::StatusCode::internal_error,
          "Auto CFD Surface Mesher encountered an out-of-range local triangle vertex while injecting the mesh."
        );
      }

      std::array<EntityRef, 3> refs {
        local_node_refs[triangle.vertices[0]],
        local_node_refs[triangle.vertices[1]],
        local_node_refs[triangle.vertices[2]],
      };


      if(refs[0] == refs[1] || refs[0] == refs[2] || refs[1] == refs[2]) {
        ++collapsed_triangle_count;
        continue;
      }

      auto oriented_triangle = triangle.vertices;
      if(should_flip_triangle_winding(*face_view, face_mesh, oriented_triangle)) {
        std::swap(oriented_triangle[1], oriented_triangle[2]);
        std::swap(refs[1], refs[2]);
      }

      const auto face_ref = output.add_triangle_face(geo_face_entity_groups[face_mesh.face.index], refs);
      output.set_face_topology_owner(face_ref, face_mesh.face);
      ++emitted_triangle_count;
      const auto edge01_it = face_mesh.constrained_edge_owners.find(
        canonical_edge(oriented_triangle[0], oriented_triangle[1])
      );
      const auto edge12_it = face_mesh.constrained_edge_owners.find(
        canonical_edge(oriented_triangle[1], oriented_triangle[2])
      );
      const auto edge20_it = face_mesh.constrained_edge_owners.find(
        canonical_edge(oriented_triangle[2], oriented_triangle[0])
      );
      record_mesh_edge(
        mesh_edges,
        refs[0],
        refs[1],
        face_ref,
        edge01_it != face_mesh.constrained_edge_owners.end()
          ? edge01_it->second
          : geo::TopologyEntityId {}
      );
      record_mesh_edge(
        mesh_edges,
        refs[1],
        refs[2],
        face_ref,
        edge12_it != face_mesh.constrained_edge_owners.end()
          ? edge12_it->second
          : geo::TopologyEntityId {}
      );
      record_mesh_edge(
        mesh_edges,
        refs[2],
        refs[0],
        face_ref,
        edge20_it != face_mesh.constrained_edge_owners.end()
          ? edge20_it->second
          : geo::TopologyEntityId {}
      );
    }

    if(emitted_triangle_count == 0U) {
      SQMESH_LOG_WARN(
        "inject_face_meshes emitted no triangles for {} (plan_triangles={}, "
        "plan_vertices={}, collapsed_to_shared_node={}). Final mesh will "
        "be missing this face.",
        topology_entity_debug_label(face_mesh.face),
        face_mesh.triangles.size(),
        face_mesh.vertices.size(),
        collapsed_triangle_count);
    }
  }

  // Create one edge EntityGroup per geometric edge; interior edges get a shared entity_group.
  std::unordered_map<std::uint32_t, EntityGroupIndex> geo_edge_entity_groups;
  EntityGroupIndex interior_edge_entity_group = invalid_index;

  for(const auto &entry : mesh_edges) {
    const auto &inc = entry.second;
    const bool has_owner = geo::is_valid(inc.topology_owner) &&
                           !inc.topology_owner_conflicted;
    EntityGroupIndex target;
    if(has_owner) {
      auto it = geo_edge_entity_groups.find(inc.topology_owner.index);
      if(it == geo_edge_entity_groups.end()) {
        EntityGroupDefinition def;
        def.order = EntityOrder::edge;
        def.name = "edge_" + std::to_string(inc.topology_owner.index);
        def.zone_id = inc.topology_owner.index;
        def.source_entity_tag = inc.topology_owner.index;
        def.default_kind = EntityKind::edge_line;
        target = output.create_entity_group(def);
        geo_edge_entity_groups[inc.topology_owner.index] = target;
      }
      else {
        target = it->second;
      }
    }
    else {
      if(interior_edge_entity_group == invalid_index) {
        EntityGroupDefinition def;
        def.order = EntityOrder::edge;
        def.name = "interior_edges";
        def.default_kind = EntityKind::edge_line;
        interior_edge_entity_group = output.create_entity_group(def);
      }
      target = interior_edge_entity_group;
    }

    const auto edge_ref = output.add_edge(target, inc.nodes);
    output.set_edge_faces(edge_ref, inc.left_face, inc.right_face);
    if(has_owner) {
      output.set_edge_topology_owner(edge_ref, inc.topology_owner);
    }
  }

  return core::detail::clear_error_state();
}

[[nodiscard]] double triangle_longest_edge_length(
  const Domain &domain,
  EntityRef face_ref
) noexcept
{
  const auto nodes = domain.face_nodes(face_ref);
  if(nodes.size != 3U) {
    return 0.0;
  }

  const auto &a = domain.node(nodes[0]).coordinates;
  const auto &b = domain.node(nodes[1]).coordinates;
  const auto &c = domain.node(nodes[2]).coordinates;
  const auto ab = point_subtract(b, a);
  const auto bc = point_subtract(c, b);
  const auto ca = point_subtract(a, c);
  return std::max({
    std::sqrt(dot_product(ab, ab)),
    std::sqrt(dot_product(bc, bc)),
    std::sqrt(dot_product(ca, ca)),
  });
}

[[nodiscard]] geo::Point3 triangle_centroid(
  const Domain &domain,
  EntityRef face_ref
) noexcept
{
  const auto nodes = domain.face_nodes(face_ref);
  if(nodes.size != 3U) {
    return {0.0, 0.0, 0.0};
  }

  const auto &a = domain.node(nodes[0]).coordinates;
  const auto &b = domain.node(nodes[1]).coordinates;
  const auto &c = domain.node(nodes[2]).coordinates;
  return {
    (a[0] + b[0] + c[0]) / 3.0,
    (a[1] + b[1] + c[1]) / 3.0,
    (a[2] + b[2] + c[2]) / 3.0,
  };
}

[[nodiscard]] geo::Vector3 triangle_normal(
  const Domain &domain,
  EntityRef face_ref
) noexcept
{
  const auto nodes = domain.face_nodes(face_ref);
  if(nodes.size != 3U) {
    return {0.0, 0.0, 0.0};
  }

  const auto &a = domain.node(nodes[0]).coordinates;
  const auto &b = domain.node(nodes[1]).coordinates;
  const auto &c = domain.node(nodes[2]).coordinates;
  return cross_product(point_subtract(b, a), point_subtract(c, a));
}

[[nodiscard]] bool recover_face_triangle_uvs(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceFacePreprocessState *face_preprocess,
  const std::array<geo::Point3, 3> &triangle_points,
  std::array<std::array<double, 2>, 3> &triangle_uvs
)
{
  for(std::size_t point_index = 0U; point_index < triangle_points.size(); ++point_index) {
    geo::FaceUvMapping mapping;
    if(geo::recover_face_uv(face_view, triangle_points[point_index], mapping) !=
       base::StatusCode::ok) {
      return false;
    }

    triangle_uvs[point_index] = {mapping.u, mapping.v};
    if(point_index > 0U &&
       face_preprocess != nullptr &&
       uses_supported_seam_material_screen(*face_preprocess)) {
      align_uv_to_reference_periodically(
        face_preprocess->uv_bounds,
        triangle_uvs[point_index - 1U],
        triangle_uvs[point_index]
      );
    }
  }

  return true;
}

void record_quality_failure_details(
  testing::AutoCfdSurfaceFinalScreenDiagnostics &diagnostics,
  const ElementQuality &quality,
  const QualityGateFailures &failures
) noexcept
{
  if(failures.minimum_angle) {
    ++diagnostics.minimum_angle_rejection_count;
    if(std::isfinite(quality.min_angle)) {
      diagnostics.rejected_min_angle_minimum =
        std::min(diagnostics.rejected_min_angle_minimum, quality.min_angle);
    }
  }
  if(failures.maximum_angle) {
    ++diagnostics.maximum_angle_rejection_count;
    if(std::isfinite(quality.max_angle)) {
      diagnostics.rejected_max_angle_maximum =
        std::max(diagnostics.rejected_max_angle_maximum, quality.max_angle);
    }
  }
  if(failures.aspect_ratio) {
    ++diagnostics.aspect_ratio_rejection_count;
    if(std::isfinite(quality.aspect_ratio)) {
      diagnostics.rejected_aspect_ratio_maximum =
        std::max(diagnostics.rejected_aspect_ratio_maximum, quality.aspect_ratio);
    }
  }
  if(failures.radius_ratio) {
    ++diagnostics.radius_ratio_rejection_count;
    if(std::isfinite(quality.radius_ratio)) {
      diagnostics.rejected_radius_ratio_minimum =
        std::min(diagnostics.rejected_radius_ratio_minimum, quality.radius_ratio);
    }
  }
  if(failures.skewness) {
    ++diagnostics.skewness_rejection_count;
    if(std::isfinite(quality.skewness)) {
      diagnostics.rejected_skewness_maximum =
        std::max(diagnostics.rejected_skewness_maximum, quality.skewness);
    }
  }
}

void record_final_screen_failure(
  testing::AutoCfdSurfaceFinalScreenDiagnostics &diagnostics,
  testing::AutoCfdSurfaceFinalScreenFailureKind kind
) noexcept
{
  if(diagnostics.first_failure_kind ==
     testing::AutoCfdSurfaceFinalScreenFailureKind::none) {
    diagnostics.first_failure_kind = kind;
  }

  switch(kind) {
  case testing::AutoCfdSurfaceFinalScreenFailureKind::none:
    break;
  case testing::AutoCfdSurfaceFinalScreenFailureKind::non_triangle_face:
    ++diagnostics.non_triangle_face_count;
    break;
  case testing::AutoCfdSurfaceFinalScreenFailureKind::duplicate_triangle:
    ++diagnostics.duplicate_triangle_count;
    break;
  case testing::AutoCfdSurfaceFinalScreenFailureKind::non_manifold_edge:
    ++diagnostics.non_manifold_edge_count;
    break;
  case testing::AutoCfdSurfaceFinalScreenFailureKind::degenerate_triangle:
    ++diagnostics.degenerate_triangle_count;
    break;
  case testing::AutoCfdSurfaceFinalScreenFailureKind::quality_gate_rejection:
    ++diagnostics.quality_gate_rejection_count;
    break;
  case testing::AutoCfdSurfaceFinalScreenFailureKind::missing_face_owner:
    ++diagnostics.missing_face_owner_count;
    break;
  case testing::AutoCfdSurfaceFinalScreenFailureKind::missing_face_view:
    ++diagnostics.missing_face_view_count;
    break;
  case testing::AutoCfdSurfaceFinalScreenFailureKind::node_reprojection_failure:
    ++diagnostics.node_reprojection_failure_count;
    break;
  case testing::AutoCfdSurfaceFinalScreenFailureKind::centroid_left_source_face:
    ++diagnostics.centroid_left_source_face_count;
    break;
  case testing::AutoCfdSurfaceFinalScreenFailureKind::orientation_flip:
    ++diagnostics.orientation_flip_count;
    break;
  case testing::AutoCfdSurfaceFinalScreenFailureKind::no_surface_triangles:
    break;
  case testing::AutoCfdSurfaceFinalScreenFailureKind::internal_quality_guardrail_rejection:
    ++diagnostics.internal_quality_guardrail_rejection_count;
    break;
  case testing::AutoCfdSurfaceFinalScreenFailureKind::boundary_topology_mismatch:
    ++diagnostics.boundary_topology_mismatch_count;
    break;
  }
}

[[nodiscard]] base::StatusCode final_screen_failure_status(
  testing::AutoCfdSurfaceFinalScreenFailureKind kind
) noexcept
{
  switch(kind) {
  case testing::AutoCfdSurfaceFinalScreenFailureKind::missing_face_owner:
  case testing::AutoCfdSurfaceFinalScreenFailureKind::missing_face_view:
    return base::StatusCode::internal_error;
  case testing::AutoCfdSurfaceFinalScreenFailureKind::none:
    return base::StatusCode::ok;
  case testing::AutoCfdSurfaceFinalScreenFailureKind::non_triangle_face:
  case testing::AutoCfdSurfaceFinalScreenFailureKind::duplicate_triangle:
  case testing::AutoCfdSurfaceFinalScreenFailureKind::non_manifold_edge:
  case testing::AutoCfdSurfaceFinalScreenFailureKind::degenerate_triangle:
  case testing::AutoCfdSurfaceFinalScreenFailureKind::quality_gate_rejection:
  case testing::AutoCfdSurfaceFinalScreenFailureKind::node_reprojection_failure:
  case testing::AutoCfdSurfaceFinalScreenFailureKind::centroid_left_source_face:
  case testing::AutoCfdSurfaceFinalScreenFailureKind::orientation_flip:
  case testing::AutoCfdSurfaceFinalScreenFailureKind::no_surface_triangles:
  case testing::AutoCfdSurfaceFinalScreenFailureKind::internal_quality_guardrail_rejection:
  case testing::AutoCfdSurfaceFinalScreenFailureKind::boundary_topology_mismatch:
    return base::StatusCode::unsupported;
  }

  return base::StatusCode::internal_error;
}

[[nodiscard]] std::string final_screen_failure_message(
  const testing::AutoCfdSurfaceFinalScreenDiagnostics &diagnostics
)
{
  using FailureKind = testing::AutoCfdSurfaceFinalScreenFailureKind;

  const auto append_metric_summary =
    [&diagnostics](std::string prefix) -> std::string {
    prefix += " (min-angle=" + std::to_string(diagnostics.minimum_angle_rejection_count);
    if(std::isfinite(diagnostics.rejected_min_angle_minimum)) {
      prefix += ", worst_min_angle_deg=" +
                std::to_string(diagnostics.rejected_min_angle_minimum);
    }
    prefix += ", max-angle=" + std::to_string(diagnostics.maximum_angle_rejection_count);
    if(diagnostics.rejected_max_angle_maximum > 0.0) {
      prefix += ", worst_max_angle_deg=" +
                std::to_string(diagnostics.rejected_max_angle_maximum);
    }
    prefix += ", aspect=" + std::to_string(diagnostics.aspect_ratio_rejection_count);
    if(diagnostics.rejected_aspect_ratio_maximum > 0.0) {
      prefix += ", worst_aspect_ratio=" +
                std::to_string(diagnostics.rejected_aspect_ratio_maximum);
    }
    prefix += ", radius-ratio=" + std::to_string(diagnostics.radius_ratio_rejection_count);
    if(std::isfinite(diagnostics.rejected_radius_ratio_minimum)) {
      prefix += ", worst_radius_ratio=" +
                std::to_string(diagnostics.rejected_radius_ratio_minimum);
    }
    prefix += ", skewness=" + std::to_string(diagnostics.skewness_rejection_count);
    if(diagnostics.rejected_skewness_maximum > 0.0) {
      prefix += ", worst_skewness=" +
                std::to_string(diagnostics.rejected_skewness_maximum);
    }
    prefix += ").";
    return prefix;
  };

  switch(diagnostics.first_failure_kind) {
  case FailureKind::none:
    return {};
  case FailureKind::non_triangle_face:
    return "Auto CFD Surface Mesher final screening requires triangle-only face entity_groups.";
  case FailureKind::duplicate_triangle:
    return "Auto CFD Surface Mesher rejected " +
           std::to_string(diagnostics.duplicate_triangle_count) +
           " duplicate surface triangle candidate(s) during final screening.";
  case FailureKind::non_manifold_edge:
    return "Auto CFD Surface Mesher rejected " +
           std::to_string(diagnostics.non_manifold_edge_count) +
           " non-manifold edge candidate(s) during final screening.";
  case FailureKind::degenerate_triangle:
    return "Auto CFD Surface Mesher rejected " +
           std::to_string(diagnostics.degenerate_triangle_count) +
           " degenerate surface triangle(s) during final screening.";
  case FailureKind::quality_gate_rejection:
    return append_metric_summary(
      "Auto CFD Surface Mesher rejected " +
      std::to_string(diagnostics.quality_gate_rejection_count) +
      " surface triangle(s) at the delivered quality gate during final screening"
    );
  case FailureKind::missing_face_owner:
    return "Auto CFD Surface Mesher final screening requires a face topology owner for every generated triangle.";
  case FailureKind::missing_face_view:
    return "Auto CFD Surface Mesher final screening could not resolve a generated face owner back to the geometry view.";
  case FailureKind::node_reprojection_failure:
    return "Auto CFD Surface Mesher rejected " +
           std::to_string(diagnostics.node_reprojection_failure_count) +
           " triangle(s) whose node could not be reprojected onto the source face.";
  case FailureKind::centroid_left_source_face:
    return "Auto CFD Surface Mesher rejected " +
           std::to_string(diagnostics.centroid_left_source_face_count) +
           " triangle(s) whose centroid left the source face during final screening.";
  case FailureKind::orientation_flip:
    return "Auto CFD Surface Mesher rejected " +
           std::to_string(diagnostics.orientation_flip_count) +
           " triangle(s) whose normal flipped against the source face.";
  case FailureKind::no_surface_triangles:
    return "Auto CFD Surface Mesher produced no surface triangles during final screening.";
  case FailureKind::internal_quality_guardrail_rejection:
    return append_metric_summary(
      "Auto CFD Surface Mesher rejected " +
      std::to_string(diagnostics.internal_quality_guardrail_rejection_count) +
      " surface triangle(s) at the internal development quality guardrail during final screening"
    );
  case FailureKind::boundary_topology_mismatch:
    return "Auto CFD Surface Mesher rejected " +
           std::to_string(diagnostics.boundary_topology_mismatch_count) +
           " face owner(s) whose delivered boundary topology no longer matched the expected loop count or closed-loop boundary structure.";
  }

  return "Auto CFD Surface Mesher final screening failed unexpectedly.";
}

void inspect_auto_cfd_surface_quality_gate_impl(
  const Domain &domain,
  const testing::AutoCfdSurfaceQualityGate &gate,
  const testing::AutoCfdSurfaceQualityGate *internal_guardrail,
  testing::AutoCfdSurfaceFinalScreenDiagnostics &diagnostics
)
{
  std::unordered_map<MeshEdgeKey, std::size_t, MeshEdgeKeyHash> edge_face_counts;
  std::unordered_set<MeshFaceKey, MeshFaceKeyHash> unique_faces;

  for(const auto &entity_group : domain.entity_groups()) {
    if(entity_group.order() != EntityOrder::face ||
       entity_group.role() != EntityGroupRole::computational) {
      continue;
    }

    for(std::uint32_t face_index = 0U; face_index < entity_group.faces().size(); ++face_index) {
      const EntityRef face_ref {entity_group.id(), face_index};
      const auto face_nodes = domain.face_nodes(face_ref);
      if(face_nodes.size != 3U) {
        record_final_screen_failure(
          diagnostics,
          testing::AutoCfdSurfaceFinalScreenFailureKind::non_triangle_face
        );
        continue;
      }

      ++diagnostics.total_triangle_count;
      const auto face_key =
        canonical_face(face_nodes[0], face_nodes[1], face_nodes[2]);
      if(!unique_faces.emplace(face_key).second) {
        record_final_screen_failure(
          diagnostics,
          testing::AutoCfdSurfaceFinalScreenFailureKind::duplicate_triangle
        );
      }

      for(std::size_t node_index = 0U; node_index < face_nodes.size; ++node_index) {
        const auto edge_key = MeshEdgeKey {
          std::min(
            pack_entity_ref(face_nodes[node_index]),
            pack_entity_ref(face_nodes[(node_index + 1U) % face_nodes.size])
          ),
          std::max(
            pack_entity_ref(face_nodes[node_index]),
            pack_entity_ref(face_nodes[(node_index + 1U) % face_nodes.size])
          ),
        };
        const auto count = ++edge_face_counts[edge_key];
        if(count > 2U) {
          record_final_screen_failure(
            diagnostics,
            testing::AutoCfdSurfaceFinalScreenFailureKind::non_manifold_edge
          );
        }
      }

      const auto quality = domain.element_quality(face_ref);
      if(!quality.supported || quality.status != QualityStatus::valid) {
        record_final_screen_failure(
          diagnostics,
          testing::AutoCfdSurfaceFinalScreenFailureKind::degenerate_triangle
        );
        continue;
      }
      if(internal_guardrail != nullptr) {
        const auto guardrail_failures =
          evaluate_quality_gate_failures(quality, *internal_guardrail);
        if(guardrail_failures.any()) {
          record_quality_failure_details(diagnostics, quality, guardrail_failures);
          record_final_screen_failure(
            diagnostics,
            testing::AutoCfdSurfaceFinalScreenFailureKind::internal_quality_guardrail_rejection
          );
          continue;
        }
      }

      const auto delivered_failures = evaluate_quality_gate_failures(quality, gate);
      if(delivered_failures.any()) {
        record_quality_failure_details(diagnostics, quality, delivered_failures);
        record_final_screen_failure(
          diagnostics,
          testing::AutoCfdSurfaceFinalScreenFailureKind::quality_gate_rejection
        );
      }
    }
  }

  if(diagnostics.total_triangle_count == 0U) {
    record_final_screen_failure(
      diagnostics,
      testing::AutoCfdSurfaceFinalScreenFailureKind::no_surface_triangles
    );
  }
}

[[nodiscard]] base::StatusCode validate_auto_cfd_surface_quality_gate_impl(
  const Domain &domain,
  const testing::AutoCfdSurfaceQualityGate &gate
)
{
  testing::AutoCfdSurfaceFinalScreenDiagnostics diagnostics;
  inspect_auto_cfd_surface_quality_gate_impl(domain, gate, nullptr, diagnostics);
  if(!diagnostics.has_failures()) {
    return core::detail::clear_error_state();
  }

  return core::detail::publish_error(
    final_screen_failure_status(diagnostics.first_failure_kind),
    final_screen_failure_message(diagnostics)
  );
}

void inspect_auto_cfd_surface_projection_and_orientation_impl(
  const geo::ModelView &model_view,
  const AutoCfdSurfacePipelineState *pipeline_state,
  const Domain &domain,
  testing::AutoCfdSurfaceFinalScreenDiagnostics &diagnostics
)
{
  for(const auto &entity_group : domain.entity_groups()) {
    if(entity_group.order() != EntityOrder::face ||
       entity_group.role() != EntityGroupRole::computational) {
      continue;
    }

    for(std::uint32_t face_index = 0U; face_index < entity_group.faces().size(); ++face_index) {
      const EntityRef face_ref {entity_group.id(), face_index};
      const auto owner = domain.face_topology_owner(face_ref);
      if(!geo::is_valid(owner) || owner.dimension != geo::TopologyDimension::face) {
        record_final_screen_failure(
          diagnostics,
          testing::AutoCfdSurfaceFinalScreenFailureKind::missing_face_owner
        );
        continue;
      }

      const auto *face_view = model_view.find_face(owner);
      if(face_view == nullptr) {
        record_final_screen_failure(
          diagnostics,
          testing::AutoCfdSurfaceFinalScreenFailureKind::missing_face_view
        );
        continue;
      }
      const auto *face_preprocess =
        pipeline_state == nullptr
          ? nullptr
          : find_auto_cfd_surface_face_preprocess_state(*pipeline_state, owner);

      const auto face_nodes = domain.face_nodes(face_ref);
      const double tolerance =
        projection_tolerance(triangle_longest_edge_length(domain, face_ref));
      std::array<geo::Point3, 3> projected_points {};
      bool node_projection_failed = false;
      for(std::size_t node_index = 0U; node_index < face_nodes.size; ++node_index) {
        geo::FaceProjection projection;
        if(geo::project_point_to_face(
             *face_view,
             domain.node(face_nodes[node_index]).coordinates,
             projection
           ) != base::StatusCode::ok ||
           projection.distance > tolerance) {
          record_final_screen_failure(
            diagnostics,
            testing::AutoCfdSurfaceFinalScreenFailureKind::node_reprojection_failure
          );
          node_projection_failed = true;
          break;
        }
        projected_points[node_index] = projection.projected_point;
      }
      if(node_projection_failed) {
        continue;
      }

      if(face_preprocess != nullptr &&
         uses_uv_material_face_screen(*face_preprocess)) {
        std::array<std::array<double, 2>, 3> triangle_uvs {};
        if(!recover_face_triangle_uvs(
             *face_view,
             face_preprocess,
             projected_points,
             triangle_uvs
           )) {
          record_final_screen_failure(
            diagnostics,
            testing::AutoCfdSurfaceFinalScreenFailureKind::node_reprojection_failure
          );
          continue;
        }

        LocalTriangleSourceFaceFailureKind failure =
          LocalTriangleSourceFaceFailureKind::none;
        if(!triangle_passes_uv_material_face_screen(
             *face_view,
             *face_preprocess,
             projected_points,
             triangle_uvs,
             &failure
           )) {
          switch(failure) {
          case LocalTriangleSourceFaceFailureKind::orientation_flip:
            record_final_screen_failure(
              diagnostics,
              testing::AutoCfdSurfaceFinalScreenFailureKind::orientation_flip
            );
            break;
          case LocalTriangleSourceFaceFailureKind::node_reprojection_failure:
            record_final_screen_failure(
              diagnostics,
              testing::AutoCfdSurfaceFinalScreenFailureKind::node_reprojection_failure
            );
            break;
          case LocalTriangleSourceFaceFailureKind::source_face_containment_failure:
          default:
            record_final_screen_failure(
              diagnostics,
              testing::AutoCfdSurfaceFinalScreenFailureKind::centroid_left_source_face
            );
            break;
          }
        }
        continue;
      }

      geo::FaceProjection centroid_projection;
      if(geo::project_point_to_face(
           *face_view,
           triangle_centroid(domain, face_ref),
           centroid_projection
         ) != base::StatusCode::ok ||
         centroid_projection.distance > tolerance) {
        record_final_screen_failure(
          diagnostics,
          testing::AutoCfdSurfaceFinalScreenFailureKind::centroid_left_source_face
        );
        continue;
      }

      if(centroid_projection.normal_defined) {
        const auto normal = normalized_vector(triangle_normal(domain, face_ref));
        if(dot_product(normal, centroid_projection.normal) <= 0.0) {
          record_final_screen_failure(
            diagnostics,
            testing::AutoCfdSurfaceFinalScreenFailureKind::orientation_flip
          );
        }
      }
    }
  }
}

[[nodiscard]] base::StatusCode validate_auto_cfd_surface_projection_and_orientation(
  const geo::ModelView &model_view,
  const AutoCfdSurfacePipelineState *pipeline_state,
  const Domain &domain
)
{
  testing::AutoCfdSurfaceFinalScreenDiagnostics diagnostics;
  inspect_auto_cfd_surface_projection_and_orientation_impl(
    model_view,
    pipeline_state,
    domain,
    diagnostics
  );
  if(!diagnostics.has_failures()) {
    return core::detail::clear_error_state();
  }

  return core::detail::publish_error(
    final_screen_failure_status(diagnostics.first_failure_kind),
    final_screen_failure_message(diagnostics)
  );
}

[[nodiscard]] std::size_t expected_owner_boundary_component_count(
  const AutoCfdSurfaceFacePreprocessState &face_preprocess
) noexcept
{
  if(uses_supported_seam_material_screen(face_preprocess)) {
    return 2U;
  }
  if(face_preprocess.seam_support ==
       AutoCfdSurfaceSeamSupportKind::periodic) {
    // For a periodic face, the number of OPEN boundary components in
    // the final mesh is the number of loops whose edges are NOT all
    // seam-or-degenerate. A loop whose every edge is either a seam
    // edge or a degenerate edge gets "stitched" back together during
    // the populate_candidate_stats collapse (boundary_node_index
    // dedup + the v[0] != v[1] guard), leaving no open boundary. This correctly
    // returns 0 for a fully closed sphere (1 loop with 2 seam + 2
    // degenerate edges) or torus (1 loop with 2 seam edges in each
    // periodic direction), and returns >0 for a cylinder-without-caps
    // or cone-without-caps where some loops still have real boundary
    // edges after the stitch.
    std::size_t open_loops = 0U;
    for(const auto &loop : face_preprocess.boundary.loops) {
      const auto total = loop.edge_uses.size();
      const auto collapsed =
        loop.seam_edge_use_count + loop.degenerate_edge_use_count;
      if(total > collapsed) {
        ++open_loops;
      }
    }
    return open_loops;
  }
  if(face_preprocess.boundary.outer_loop_count == 1U &&
     face_preprocess.boundary.unknown_loop_count == 0U) {
    return 1U + face_preprocess.boundary.inner_loop_count;
  }

  return 0U;
}

void inspect_auto_cfd_surface_boundary_topology_impl(
  const AutoCfdSurfacePipelineState *pipeline_state,
  const Domain &domain,
  testing::AutoCfdSurfaceFinalScreenDiagnostics &diagnostics
)
{
  if(pipeline_state == nullptr) {
    return;
  }

  std::unordered_map<std::uint64_t, std::unordered_map<std::uint64_t, std::vector<std::uint64_t>>>
    adjacency_by_owner;

  const auto append_owner_boundary_edge =
    [&adjacency_by_owner](
      geo::TopologyEntityId owner,
      EntityRef first,
      EntityRef second
    ) {
      if(!geo::is_valid(owner) || owner.dimension != geo::TopologyDimension::face) {
        return;
      }

      auto &adjacency = adjacency_by_owner[pack_topology_entity_id(owner)];
      const auto first_key = pack_entity_ref(first);
      const auto second_key = pack_entity_ref(second);
      adjacency[first_key].push_back(second_key);
      adjacency[second_key].push_back(first_key);
    };

  for(const auto &entity_group : domain.entity_groups()) {
    if(entity_group.order() != EntityOrder::edge ||
       entity_group.role() != EntityGroupRole::computational) {
      continue;
    }

    for(std::uint32_t edge_index = 0U; edge_index < entity_group.edges().size(); ++edge_index) {
      const EntityRef edge_ref {entity_group.id(), edge_index};
      const auto edge_nodes = domain.edge_nodes(edge_ref);
      if(edge_nodes.size != 2U) {
        continue;
      }

      const auto left_face = domain.adjacent_face(edge_ref, FaceSide::left);
      const auto right_face = domain.adjacent_face(edge_ref, FaceSide::right);
      const auto left_owner =
        is_valid(left_face)
          ? domain.face_topology_owner(left_face)
          : geo::TopologyEntityId {};
      const auto right_owner =
        is_valid(right_face)
          ? domain.face_topology_owner(right_face)
          : geo::TopologyEntityId {};

      if(geo::is_valid(left_owner) &&
         (!geo::is_valid(right_owner) || right_owner != left_owner)) {
        append_owner_boundary_edge(left_owner, edge_nodes[0], edge_nodes[1]);
      }
      if(geo::is_valid(right_owner) &&
         (!geo::is_valid(left_owner) || left_owner != right_owner)) {
        append_owner_boundary_edge(right_owner, edge_nodes[0], edge_nodes[1]);
      }
    }
  }

  for(const auto &face_preprocess : pipeline_state->face_preprocess_states) {
    const auto expected_components =
      expected_owner_boundary_component_count(face_preprocess);
    if(expected_components == 0U) {
      continue;
    }

    const auto owner_key = pack_topology_entity_id(face_preprocess.face);
    const auto owner_it = adjacency_by_owner.find(owner_key);
    if(owner_it == adjacency_by_owner.end()) {
      record_final_screen_failure(
        diagnostics,
        testing::AutoCfdSurfaceFinalScreenFailureKind::boundary_topology_mismatch
      );
      SQMESH_LOG_WARN(
        "Boundary topology mismatch on {}: expected {} component(s) "
        "(outer={}, inner={}, unknown={}, seam={}, degenerate={}), "
        "but no boundary edges were emitted in the final mesh.",
        topology_entity_debug_label(face_preprocess.face),
        expected_components,
        face_preprocess.boundary.outer_loop_count,
        face_preprocess.boundary.inner_loop_count,
        face_preprocess.boundary.unknown_loop_count,
        face_preprocess.boundary.seam_loop_count,
        face_preprocess.boundary.degenerate_loop_count);
      continue;
    }

    auto unvisited = std::unordered_set<std::uint64_t> {};
    unvisited.reserve(owner_it->second.size());
    for(const auto &[node_key, neighbors] : owner_it->second) {
      unvisited.insert(node_key);
    }

    std::size_t component_count = 0U;
    std::vector<std::uint64_t> stack;
    while(!unvisited.empty()) {
      ++component_count;
      stack.push_back(*unvisited.begin());
      unvisited.erase(unvisited.begin());

      while(!stack.empty()) {
        const auto current = stack.back();
        stack.pop_back();
        const auto adjacency_it = owner_it->second.find(current);
        if(adjacency_it == owner_it->second.end()) {
          continue;
        }

        for(const auto neighbor : adjacency_it->second) {
          const auto unvisited_it = unvisited.find(neighbor);
          if(unvisited_it == unvisited.end()) {
            continue;
          }
          stack.push_back(neighbor);
          unvisited.erase(unvisited_it);
        }
      }
    }

    if(component_count != expected_components) {
      record_final_screen_failure(
        diagnostics,
        testing::AutoCfdSurfaceFinalScreenFailureKind::boundary_topology_mismatch
      );
      SQMESH_LOG_WARN(
        "Boundary topology mismatch on {}: expected {} component(s) "
        "(outer={}, inner={}, unknown={}, seam={}, degenerate={}), "
        "but the final mesh delivered {} boundary component(s).",
        topology_entity_debug_label(face_preprocess.face),
        expected_components,
        face_preprocess.boundary.outer_loop_count,
        face_preprocess.boundary.inner_loop_count,
        face_preprocess.boundary.unknown_loop_count,
        face_preprocess.boundary.seam_loop_count,
        face_preprocess.boundary.degenerate_loop_count,
        component_count);
    }
  }
}

[[nodiscard]] base::StatusCode validate_auto_cfd_surface_boundary_topology(
  const AutoCfdSurfacePipelineState *pipeline_state,
  const Domain &domain
)
{
  testing::AutoCfdSurfaceFinalScreenDiagnostics diagnostics;
  inspect_auto_cfd_surface_boundary_topology_impl(
    pipeline_state,
    domain,
    diagnostics
  );
  if(!diagnostics.has_failures()) {
    return core::detail::clear_error_state();
  }

  return core::detail::publish_error(
    final_screen_failure_status(diagnostics.first_failure_kind),
    final_screen_failure_message(diagnostics)
  );
}

void inspect_auto_cfd_surface_final_screen_impl(
  const geo::ModelView &model_view,
  const AutoCfdSurfacePipelineState *pipeline_state,
  const Domain &domain,
  testing::AutoCfdSurfaceFinalScreenDiagnostics &diagnostics
)
{
  diagnostics = {};
  const auto guardrail = internal_quality_guardrail();
  inspect_auto_cfd_surface_quality_gate_impl(
    domain,
    delivered_quality_gate(),
    &guardrail,
    diagnostics
  );
  inspect_auto_cfd_surface_projection_and_orientation_impl(
    model_view,
    pipeline_state,
    domain,
    diagnostics
  );
  inspect_auto_cfd_surface_boundary_topology_impl(
    pipeline_state,
    domain,
    diagnostics
  );
}

[[nodiscard]] base::StatusCode validate_generated_auto_cfd_surface_meshes(
  const geo::ModelView &model_view,
  const AutoCfdSurfacePipelineState &pipeline_state,
  const Domain &domain,
  bool allow_quality_gate_failure,
  bool allow_final_screen_failure
)
{
  if(allow_final_screen_failure) {
    auto status =
      validate_auto_cfd_surface_quality_gate_impl(domain, delivered_quality_gate());
    if(status != base::StatusCode::ok) {
      SQMESH_LOG_WARN("{} Keeping generated mesh for inspection because"
        " 'allow_final_screen_failure' is enabled.",
        sqmesh::base::last_error_message());
    }

    testing::AutoCfdSurfaceFinalScreenDiagnostics boundary_diagnostics;
    inspect_auto_cfd_surface_boundary_topology_impl(
      &pipeline_state,
      domain,
      boundary_diagnostics
    );
    if(boundary_diagnostics.has_failures()) {
      SQMESH_LOG_WARN("{} Keeping generated mesh for inspection because"
        " 'allow_final_screen_failure' is enabled.",
        final_screen_failure_message(boundary_diagnostics));
    }
    log_open_mesh_edge_owner_breakdown(domain);

    // Projection/orientation still stays skipped here because it requires many
    // OCC evaluations per triangle and is materially more expensive than the
    // topology-only debug path above.

    return core::detail::clear_error_state();
  }

  testing::AutoCfdSurfaceFinalScreenDiagnostics quality_diagnostics;
  inspect_auto_cfd_surface_quality_gate_impl(
    domain,
    delivered_quality_gate(),
    nullptr,
    quality_diagnostics
  );
  if(quality_diagnostics.has_failures()) {
    if(!(allow_quality_gate_failure &&
         quality_diagnostics.first_failure_kind ==
           testing::AutoCfdSurfaceFinalScreenFailureKind::quality_gate_rejection)) {
      return core::detail::publish_error(
        final_screen_failure_status(quality_diagnostics.first_failure_kind),
        final_screen_failure_message(quality_diagnostics)
      );
    }

    SQMESH_LOG_WARN("{} Keeping generated mesh for inspection because"
      " 'allow_quality_gate_failure' is enabled.",
      final_screen_failure_message(quality_diagnostics));
  }

  auto status = validate_auto_cfd_surface_projection_and_orientation(
    model_view,
    &pipeline_state,
    domain
  );
  if(status != base::StatusCode::ok) {
    return status;
  }

  return validate_auto_cfd_surface_boundary_topology(&pipeline_state, domain);
}

class AutoCfdSurfaceMesher final : public BottomUpMeshingPipeline
{
public:
  [[nodiscard]] std::string_view name() const noexcept override
  {
    return kAutoCfdSurfaceMesherName;
  }

  [[nodiscard]] std::string_view entity_group_prefix() const noexcept override
  {
    return "surface_";
  }

protected:
  void reset() noexcept override
  {
    BottomUpMeshingPipeline::reset();
    parameters_ = {};
    pipeline_state_ = {};
  }

  [[nodiscard]] base::StatusCode on_configure(
    const ParameterDictionary &parameters
  ) override
  {
    return resolve_auto_cfd_surface_parameters(parameters, parameters_);
  }

  [[nodiscard]] MeshingDimension pipeline_dimension() const noexcept override
  {
    return MeshingDimension::surface;
  }

  [[nodiscard]] double pipeline_global_target_size() const noexcept override
  {
    return parameters_.maximum_length;
  }

  [[nodiscard]] base::StatusCode on_prepare_pipeline() override
  {
    return build_auto_cfd_surface_pipeline_state(
      model_view(),
      size_controls(),
      parameters_,
      pipeline_state_
    );
  }

  [[nodiscard]] base::StatusCode on_stage_faces(Domain &output) override
  {
    auto status =
      build_auto_cfd_surface_candidate_domain(
        model_view(),
        pipeline_state_,
        &request(),
        &parameters_,
        output
      );
    if(status != base::StatusCode::ok) {
      return status;
    }

    status = validate_generated_auto_cfd_surface_meshes(
      model_view(),
      pipeline_state_,
      output,
      parameters_.allow_quality_gate_failure,
      parameters_.allow_final_screen_failure
    );
    if(status != base::StatusCode::ok) {
      if(!entity_group_prefix().empty()) {
  output.remove_entity_groups_with_prefix(entity_group_prefix());
      }
      return status;
    }

    return core::detail::clear_error_state();
  }

private:
  AutoCfdSurfaceParameters parameters_ {};
  AutoCfdSurfacePipelineState pipeline_state_ {};
};

} // namespace

MeshingAlgorithmPtr create_auto_cfd_surface_mesher()
{
  return std::make_unique<AutoCfdSurfaceMesher>();
}

testing::AutoCfdSurfaceQualityGate auto_cfd_surface_internal_quality_guardrail_for_tests() noexcept
{
  return internal_quality_guardrail();
}

testing::AutoCfdSurfaceQualityGate auto_cfd_surface_quality_gate_for_tests() noexcept
{
  return delivered_quality_gate();
}

base::StatusCode inspect_auto_cfd_surface_quality_gate_for_tests(
  const Domain &domain
) noexcept
{
  return validate_auto_cfd_surface_quality_gate_impl(domain, delivered_quality_gate());
}

base::StatusCode build_auto_cfd_surface_candidate_domain_for_tests(
  geo::ModelHandle model_handle,
  const ParameterDictionary &parameters,
  Domain &domain,
  base::ContextHandle context
)
{
  domain = Domain("");

  geo::ModelView model_view;
  auto status = geo::model_view(model_handle, model_view, context);
  if(status != base::StatusCode::ok) {
    return status;
  }

  AutoCfdSurfaceParameters resolved;
  status = resolve_auto_cfd_surface_parameters(parameters, resolved);
  if(status != base::StatusCode::ok) {
    return status;
  }

  ResolvedMeshSizeControls size_controls;
  status = resolve_mesh_size_controls(
    model_view,
    MeshSizeControls {},
    resolved.maximum_length,
    size_controls
  );
  if(status != base::StatusCode::ok) {
    return status;
  }

  AutoCfdSurfacePipelineState pipeline_state;
  status = build_auto_cfd_surface_pipeline_state(
    model_view,
    size_controls,
    resolved,
    pipeline_state
  );
  if(status != base::StatusCode::ok) {
    return status;
  }


  const MeshingRequest fallback_request {
    context,
    base::current_session(context),
    model_handle,
    MeshingDimension::surface,
    MeshSizeControls {},
  };
  return build_auto_cfd_surface_candidate_domain(
    model_view,
    pipeline_state,
    &fallback_request,
    &resolved,
    domain
  );
}

base::StatusCode build_auto_cfd_surface_boundary_domain_for_tests(
  geo::ModelHandle model_handle,
  const ParameterDictionary &parameters,
  Domain &domain,
  base::ContextHandle context
)
{
  domain = Domain("");

  geo::ModelView model_view;
  auto status = geo::model_view(model_handle, model_view, context);
  if(status != base::StatusCode::ok) {
    return status;
  }

  AutoCfdSurfaceParameters resolved;
  status = resolve_auto_cfd_surface_parameters(parameters, resolved);
  if(status != base::StatusCode::ok) {
    return status;
  }

  ResolvedMeshSizeControls size_controls;
  status = resolve_mesh_size_controls(
    model_view,
    MeshSizeControls {},
    resolved.maximum_length,
    size_controls
  );
  if(status != base::StatusCode::ok) {
    return status;
  }

  AutoCfdSurfacePipelineState pipeline_state;
  status = build_auto_cfd_surface_pipeline_state(
    model_view,
    size_controls,
    resolved,
    pipeline_state
  );
  if(status != base::StatusCode::ok) {
    return status;
  }
  domain.set_source_topology_revision(pipeline_state.topology_revision);
  const auto node_entity_group = domain.create_entity_group(
    {EntityOrder::node, "boundary_nodes_1d"}
  );
  const auto edge_entity_group = domain.create_entity_group(
    {EntityOrder::edge, "boundary_edges_1d", invalid_index, false, EntityKind::edge_line}
  );

  std::size_t segment_count = 0U;
  for(const auto &edge_discretization : pipeline_state.boundary.edge_discretizations) {
    if(edge_discretization.points.size() >= 2U) {
      segment_count += edge_discretization.points.size() - 1U;
    }
  }

  domain.reserve_entity_group_storage(
    node_entity_group,
    pipeline_state.boundary.nodes.size()
  );
  domain.reserve_entity_group_storage(
    edge_entity_group,
    segment_count,
    segment_count * 2U
  );

  std::vector<EntityRef> node_refs(
    pipeline_state.boundary.nodes.size(),
    EntityRef {}
  );
  for(std::size_t node_index = 0U;
      node_index < pipeline_state.boundary.nodes.size();
      ++node_index) {
    node_refs[node_index] = domain.add_node(
      node_entity_group,
      pipeline_state.boundary.nodes[node_index].position
    );
  }

  for(const auto &edge_discretization : pipeline_state.boundary.edge_discretizations) {
    if(!geo::is_valid(edge_discretization.edge) ||
       edge_discretization.points.size() < 2U) {
      continue;
    }

    for(std::size_t point_index = 0U;
        point_index + 1U < edge_discretization.points.size();
        ++point_index) {
      const auto first = edge_discretization.points[point_index].node_index;
      const auto second = edge_discretization.points[point_index + 1U].node_index;
      if(first >= node_refs.size() ||
         second >= node_refs.size() ||
         !is_valid(node_refs[first]) ||
         !is_valid(node_refs[second]) ||
         node_refs[first] == node_refs[second]) {
        continue;
      }

      const auto edge_ref = domain.add_edge(
        edge_entity_group,
        {node_refs[first], node_refs[second]}
      );
      domain.set_edge_topology_owner(edge_ref, edge_discretization.edge);
    }
  }

  return core::detail::clear_error_state();
}

base::StatusCode inspect_auto_cfd_surface_candidate_stats_for_tests(
  geo::ModelHandle model_handle,
  const ParameterDictionary &parameters,
  testing::AutoCfdSurfaceCandidateStats &stats,
  base::ContextHandle context
)
{
  stats = {};

  geo::ModelView model_view;
  auto status = geo::model_view(model_handle, model_view, context);
  if(status != base::StatusCode::ok) {
    return status;
  }

  AutoCfdSurfaceParameters resolved;
  status = resolve_auto_cfd_surface_parameters(parameters, resolved);
  if(status != base::StatusCode::ok) {
    return status;
  }

  ResolvedMeshSizeControls size_controls;
  status = resolve_mesh_size_controls(
    model_view,
    MeshSizeControls {},
    resolved.maximum_length,
    size_controls
  );
  if(status != base::StatusCode::ok) {
    return status;
  }

  AutoCfdSurfacePipelineState pipeline_state;
  status = build_auto_cfd_surface_pipeline_state(
    model_view,
    size_controls,
    resolved,
    pipeline_state
  );
  if(status != base::StatusCode::ok) {
    return status;
  }

  std::vector<FaceMeshPlan> face_meshes;
  status = build_face_mesh_plans(
    model_view,
    pipeline_state,
    nullptr,
    nullptr,
    face_meshes
  );
  if(status != base::StatusCode::ok) {
    stats = {};
    return status;
  }

  populate_candidate_stats(face_meshes, stats);
  return core::detail::clear_error_state();
}

base::StatusCode inspect_auto_cfd_surface_final_screen_for_tests(
  geo::ModelHandle model_handle,
  const ParameterDictionary &parameters,
  const Domain &domain,
  testing::AutoCfdSurfaceFinalScreenDiagnostics &diagnostics,
  base::ContextHandle context
)
{
  geo::ModelView model_view;
  auto status = geo::model_view(model_handle, model_view, context);
  if(status != base::StatusCode::ok) {
    diagnostics = {};
    return status;
  }

  AutoCfdSurfacePipelineState pipeline_state;
  if(parameters.size() > 0U) {
    AutoCfdSurfaceParameters resolved;
    status = resolve_auto_cfd_surface_parameters(parameters, resolved);
    if(status != base::StatusCode::ok) {
      diagnostics = {};
      return status;
    }

    ResolvedMeshSizeControls size_controls;
    status = resolve_mesh_size_controls(
      model_view,
      MeshSizeControls {},
      resolved.maximum_length,
      size_controls
    );
    if(status != base::StatusCode::ok) {
      diagnostics = {};
      return status;
    }

    status = build_auto_cfd_surface_pipeline_state(
      model_view,
      size_controls,
      resolved,
      pipeline_state
    );
    if(status != base::StatusCode::ok) {
      diagnostics = {};
      return status;
    }
  }

  inspect_auto_cfd_surface_final_screen_impl(
    model_view,
    parameters.size() == 0U ? nullptr : &pipeline_state,
    domain,
    diagnostics
  );
  if(!diagnostics.has_failures()) {
    return core::detail::clear_error_state();
  }

  return core::detail::publish_error(
    final_screen_failure_status(diagnostics.first_failure_kind),
    final_screen_failure_message(diagnostics)
  );
}

} // namespace sqmesh::mesh::detail

namespace sqmesh::mesh::testing {

AutoCfdSurfaceQualityGate auto_cfd_surface_internal_quality_guardrail() noexcept
{
  return detail::auto_cfd_surface_internal_quality_guardrail_for_tests();
}

AutoCfdSurfaceQualityGate auto_cfd_surface_quality_gate() noexcept
{
  return detail::auto_cfd_surface_quality_gate_for_tests();
}

base::StatusCode inspect_auto_cfd_surface_quality_gate(
  const Domain &domain
) noexcept
{
  return detail::inspect_auto_cfd_surface_quality_gate_for_tests(domain);
}

base::StatusCode build_auto_cfd_surface_candidate_domain(
  geo::ModelHandle model_handle,
  const ParameterDictionary &parameters,
  Domain &domain,
  base::ContextHandle context
)
{
  return detail::build_auto_cfd_surface_candidate_domain_for_tests(
    model_handle,
    parameters,
    domain,
    context
  );
}

base::StatusCode build_auto_cfd_surface_boundary_domain(
  geo::ModelHandle model_handle,
  const ParameterDictionary &parameters,
  Domain &domain,
  base::ContextHandle context
)
{
  return detail::build_auto_cfd_surface_boundary_domain_for_tests(
    model_handle,
    parameters,
    domain,
    context
  );
}

base::StatusCode inspect_auto_cfd_surface_candidate_stats(
  geo::ModelHandle model_handle,
  const ParameterDictionary &parameters,
  AutoCfdSurfaceCandidateStats &stats,
  base::ContextHandle context
)
{
  return detail::inspect_auto_cfd_surface_candidate_stats_for_tests(
    model_handle,
    parameters,
    stats,
    context
  );
}

base::StatusCode inspect_auto_cfd_surface_final_screen(
  geo::ModelHandle model_handle,
  const ParameterDictionary &parameters,
  const Domain &domain,
  AutoCfdSurfaceFinalScreenDiagnostics &diagnostics,
  base::ContextHandle context
)
{
  return detail::inspect_auto_cfd_surface_final_screen_for_tests(
    model_handle,
    parameters,
    domain,
    diagnostics,
    context
  );
}

} // namespace sqmesh::mesh::testing
