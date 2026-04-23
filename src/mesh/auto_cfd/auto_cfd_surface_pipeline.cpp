// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "auto_cfd_surface_pipeline.hpp"
#include "auto_cfd_surface_proximity_index.hpp"
#include "../sizing/proxy_mesh_geometry.hpp"

#include "core/log.hpp"
#include "core/runtime_registry.hpp"
#include "model/geometry_model_storage.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <queue>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace sqmesh::mesh::detail {
namespace {

constexpr std::string_view kAutoCfdSurfaceMesherName = "Auto CFD Surface Mesher";
constexpr double kPi = 3.14159265358979323846;
constexpr double kCurvatureTolerance = 1.0e-12;
constexpr double kTangentContinuityToleranceDegrees = 10.0;
constexpr double kFeatureAngleDegrees = 30.0;
constexpr std::size_t kMinimumEdgeCurvatureSegmentCount = 8U;
constexpr std::size_t kMinimumFaceSampleCount = 2U;
constexpr std::size_t kMaximumFaceSampleCount = 16U;
constexpr std::size_t kMaximumEdgeDiscretizationPoints = 4096U;
constexpr std::size_t kMinimumEdgeIntegrationDepth = 6U;
constexpr std::size_t kMaximumEdgeIntegrationDepth = 20U;
constexpr std::size_t kProximityLayerCount = 3U;
constexpr double kEdgeIntegrationPrecision = 1.0e-3;
constexpr double kEdgeParameterToleranceScale = 1.0e-12;

constexpr double kEdgeSizeSmoothingRatio = 1.8;
constexpr std::size_t kEdgeSizeSmoothingMaxIterations = 2000U;
constexpr double kEdgeFilterCloseRatio = 0.3;
constexpr double kProximityNormalAngleThresholdDegrees = 140.0;
constexpr double kPeriodicUvWrapTolerance = 1.0e-9;
constexpr double kMetricEigenvalueFloor = 1.0e-12;
constexpr double kMetricEigenvalueCeiling = 1.0e12;
constexpr double kMetricDeterminantTolerance = 1.0e-18;

constexpr double kMaximumMetricEigenvalueRatio = 1.0e12;
constexpr double kBackgroundGridPaddingFactor = 1.0;
constexpr double kBackgroundGridValueTolerance = 1.0e-12;
constexpr std::size_t kBackgroundGridMaximumAxisNodes = 96U;
constexpr std::size_t kBackgroundGridMaximumNodeCount = 262144U;
constexpr std::size_t kBackgroundOctreeMaximumLeafSourceCount = 16U;
constexpr std::size_t kBackgroundOctreeMaximumDepth = 12U;
constexpr double kBackgroundOctreeMinimumExtentScale = 0.25;

struct ProximityHitCandidate final {
  bool found = false;
  AutoCfdSurfaceSizingSourceKind kind =
    AutoCfdSurfaceSizingSourceKind::self_proximity;
  geo::TopologyEntityId hit_owner {};
  double hit_distance = std::numeric_limits<double>::infinity();
  double hit_normal_angle_degrees = 0.0;
};

struct ProximityBucketUnionFind final {
  std::vector<std::uint32_t> parents {};
  std::vector<std::uint32_t> ranks {};

  explicit ProximityBucketUnionFind(std::size_t count)
    : parents(count, invalid_auto_cfd_surface_part_index),
      ranks(count, 0U)
  {
    for(std::uint32_t index = 0U; index < parents.size(); ++index) {
      parents[index] = index;
    }
  }

  [[nodiscard]] std::uint32_t find(std::uint32_t index)
  {
    if(index >= parents.size()) {
      return invalid_auto_cfd_surface_part_index;
    }
    if(parents[index] == index) {
      return index;
    }
    parents[index] = find(parents[index]);
    return parents[index];
  }

  void unite(std::uint32_t lhs, std::uint32_t rhs)
  {
    lhs = find(lhs);
    rhs = find(rhs);
    if(lhs == invalid_auto_cfd_surface_part_index ||
       rhs == invalid_auto_cfd_surface_part_index ||
       lhs == rhs) {
      return;
    }

    if(ranks[lhs] < ranks[rhs]) {
      std::swap(lhs, rhs);
    }
    parents[rhs] = lhs;
    if(ranks[lhs] == ranks[rhs]) {
      ++ranks[lhs];
    }
  }
};

struct EdgeDiscretizationSample final {
  double parameter = 0.0;
  geo::Point3 position {0.0, 0.0, 0.0};
  double target_size = 0.0;
  double speed = 0.0;
  double primitive = 0.0;
};

struct MutableBackgroundBounds final {
  geo::Point3 minimum {0.0, 0.0, 0.0};
  geo::Point3 maximum {0.0, 0.0, 0.0};
  bool defined = false;
};

struct BackgroundGridQueueEntry final {
  double value = std::numeric_limits<double>::infinity();
  std::size_t node_index = 0U;
};

struct BackgroundGridQueueGreater final {
  [[nodiscard]] bool operator()(
    const BackgroundGridQueueEntry &lhs,
    const BackgroundGridQueueEntry &rhs
  ) const noexcept
  {
    return lhs.value > rhs.value;
  }
};

[[nodiscard]] base::StatusCode require_numeric_parameter(
  const ParameterDictionary &parameters,
  std::string_view key,
  double &value
)
{
  const auto *entry = parameters.find(key);
  if(entry == nullptr) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      std::string(kAutoCfdSurfaceMesherName) + " requires a '" + std::string(key) +
        "' parameter."
    );
  }
  if(!entry->is_number()) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Parameter '" + std::string(key) + "' must be numeric."
    );
  }

  value = entry->number();
  if(!std::isfinite(value)) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Parameter '" + std::string(key) + "' must be finite."
    );
  }

  return core::detail::clear_error_state();
}

[[nodiscard]] base::StatusCode try_get_optional_numeric_parameter(
  const ParameterDictionary &parameters,
  std::string_view key,
  double &value,
  bool &specified
)
{
  specified = false;
  const auto *entry = parameters.find(key);
  if(entry == nullptr) {
    return core::detail::clear_error_state();
  }
  if(!entry->is_number()) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Parameter '" + std::string(key) + "' must be numeric."
    );
  }

  value = entry->number();
  if(!std::isfinite(value)) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Parameter '" + std::string(key) + "' must be finite."
    );
  }

  specified = true;
  return core::detail::clear_error_state();
}

[[nodiscard]] base::StatusCode try_get_optional_boolean_parameter(
  const ParameterDictionary &parameters,
  std::string_view key,
  bool &value,
  bool &specified
)
{
  specified = false;
  const auto *entry = parameters.find(key);
  if(entry == nullptr) {
    return core::detail::clear_error_state();
  }
  if(!entry->is_boolean()) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Parameter '" + std::string(key) + "' must be boolean."
    );
  }

  value = entry->boolean();
  specified = true;
  return core::detail::clear_error_state();
}

[[nodiscard]] double point_distance(
  const geo::Point3 &lhs,
  const geo::Point3 &rhs
) noexcept
{
  const double dx = lhs[0] - rhs[0];
  const double dy = lhs[1] - rhs[1];
  const double dz = lhs[2] - rhs[2];
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

void extend_bounds(MutableBackgroundBounds &bounds, const geo::Point3 &point) noexcept
{
  if(!bounds.defined) {
    bounds.minimum = point;
    bounds.maximum = point;
    bounds.defined = true;
    return;
  }

  for(std::size_t axis = 0U; axis < 3U; ++axis) {
    bounds.minimum[axis] = std::min(bounds.minimum[axis], point[axis]);
    bounds.maximum[axis] = std::max(bounds.maximum[axis], point[axis]);
  }
}

[[nodiscard]] MutableBackgroundBounds effective_background_bounds(
  const AutoCfdSurfaceSizingFieldState &field
) noexcept
{
  MutableBackgroundBounds bounds;
  if(field.background_bounds.defined) {
    bounds.minimum = field.background_bounds.minimum;
    bounds.maximum = field.background_bounds.maximum;
    bounds.defined = true;
  }

  const auto include_sources =
    [&](const std::vector<AutoCfdSurfaceSizingSource> &sources) {
      for(const auto &source : sources) {
        extend_bounds(bounds, source.position);
      }
    };
  include_sources(field.curvature_sources);
  include_sources(field.proximity_sources);
  return bounds;
}

[[nodiscard]] geo::Vector3 subtract_points(
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

[[nodiscard]] double dot_product(
  const geo::Vector3 &lhs,
  const geo::Vector3 &rhs
) noexcept
{
  return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2];
}

[[nodiscard]] double vector_norm(const geo::Vector3 &vector) noexcept
{
  return std::sqrt(dot_product(vector, vector));
}

[[nodiscard]] geo::Vector3 cross_product(
  const geo::Vector3 &lhs,
  const geo::Vector3 &rhs
) noexcept
{
  return {
    lhs[1] * rhs[2] - lhs[2] * rhs[1],
    lhs[2] * rhs[0] - lhs[0] * rhs[2],
    lhs[0] * rhs[1] - lhs[1] * rhs[0],
  };
}

[[nodiscard]] geo::Vector3 normalized(const geo::Vector3 &vector) noexcept
{
  const double length = vector_norm(vector);
  if(length <= 0.0 || !std::isfinite(length)) {
    return {0.0, 0.0, 0.0};
  }

  return {
    vector[0] / length,
    vector[1] / length,
    vector[2] / length,
  };
}

[[nodiscard]] geo::Vector3 scaled_vector(
  const geo::Vector3 &vector,
  double scale
) noexcept
{
  return {
    vector[0] * scale,
    vector[1] * scale,
    vector[2] * scale,
  };
}

[[nodiscard]] std::size_t background_grid_node_index(
  const AutoCfdSurfaceBackgroundGrid &grid,
  std::size_t i,
  std::size_t j,
  std::size_t k
) noexcept
{
  return (k * grid.node_counts[1] + j) * grid.node_counts[0] + i;
}

[[nodiscard]] std::size_t background_grid_cell_index(
  const AutoCfdSurfaceBackgroundGrid &grid,
  std::size_t i,
  std::size_t j,
  std::size_t k
) noexcept
{
  const std::size_t cell_count_x =
    grid.node_counts[0] > 0U ? grid.node_counts[0] - 1U : 0U;
  const std::size_t cell_count_y =
    grid.node_counts[1] > 0U ? grid.node_counts[1] - 1U : 0U;
  return (k * cell_count_y + j) * cell_count_x + i;
}

[[nodiscard]] std::array<std::size_t, 3> background_grid_node_ijk(
  const AutoCfdSurfaceBackgroundGrid &grid,
  std::size_t node_index
) noexcept
{
  const std::size_t plane_size = grid.node_counts[0] * grid.node_counts[1];
  const std::size_t k = plane_size == 0U ? 0U : node_index / plane_size;
  const std::size_t plane_offset =
    plane_size == 0U ? 0U : node_index - k * plane_size;
  const std::size_t j =
    grid.node_counts[0] == 0U ? 0U : plane_offset / grid.node_counts[0];
  const std::size_t i =
    grid.node_counts[0] == 0U ? 0U : plane_offset % grid.node_counts[0];
  return {i, j, k};
}

[[nodiscard]] std::size_t clamped_background_grid_node_count(
  double extent,
  double spacing
) noexcept
{
  if(!std::isfinite(extent) || extent <= 0.0 || !std::isfinite(spacing) || spacing <= 0.0) {
    return 2U;
  }

  const auto count = static_cast<std::size_t>(std::ceil(extent / spacing)) + 1U;
  return std::max<std::size_t>(2U, count);
}

[[nodiscard]] double distance_point_to_aabb(
  const geo::Point3 &point,
  const geo::Point3 &bounds_min,
  const geo::Point3 &bounds_max
) noexcept
{
  double squared_distance = 0.0;
  for(std::size_t axis = 0U; axis < 3U; ++axis) {
    const double coordinate = point[axis];
    if(coordinate < bounds_min[axis]) {
      const double delta = bounds_min[axis] - coordinate;
      squared_distance += delta * delta;
    }
    else if(coordinate > bounds_max[axis]) {
      const double delta = coordinate - bounds_max[axis];
      squared_distance += delta * delta;
    }
  }
  return std::sqrt(squared_distance);
}

[[nodiscard]] std::uint8_t background_octree_child_index(
  const geo::Point3 &bounds_min,
  const geo::Point3 &bounds_max,
  const geo::Point3 &point
) noexcept
{
  const geo::Point3 center {
    0.5 * (bounds_min[0] + bounds_max[0]),
    0.5 * (bounds_min[1] + bounds_max[1]),
    0.5 * (bounds_min[2] + bounds_max[2]),
  };

  std::uint8_t child = 0U;
  if(point[0] >= center[0]) {
    child |= 1U;
  }
  if(point[1] >= center[1]) {
    child |= 2U;
  }
  if(point[2] >= center[2]) {
    child |= 4U;
  }
  return child;
}

void background_octree_child_bounds(
  const geo::Point3 &bounds_min,
  const geo::Point3 &bounds_max,
  std::uint8_t child_index,
  geo::Point3 &child_min,
  geo::Point3 &child_max
) noexcept
{
  const geo::Point3 center {
    0.5 * (bounds_min[0] + bounds_max[0]),
    0.5 * (bounds_min[1] + bounds_max[1]),
    0.5 * (bounds_min[2] + bounds_max[2]),
  };

  for(std::size_t axis = 0U; axis < 3U; ++axis) {
    const bool upper = ((child_index >> axis) & 1U) != 0U;
    child_min[axis] = upper ? center[axis] : bounds_min[axis];
    child_max[axis] = upper ? bounds_max[axis] : center[axis];
  }
}

[[nodiscard]] std::uint32_t build_background_octree_node(
  AutoCfdSurfaceBackgroundOctree &octree,
  const AutoCfdSurfaceSizingFieldState &field,
  std::uint32_t first,
  std::uint32_t count,
  const geo::Point3 &bounds_min,
  const geo::Point3 &bounds_max,
  std::size_t depth
)
{
  const std::uint32_t node_index = static_cast<std::uint32_t>(octree.nodes.size());
  octree.nodes.push_back({});
  octree.nodes[node_index].bounds_min = bounds_min;
  octree.nodes[node_index].bounds_max = bounds_max;
  octree.nodes[node_index].first_source = first;
  octree.nodes[node_index].source_count = count;
  octree.maximum_depth = std::max(octree.maximum_depth, depth);

  if(count == 0U) {
    ++octree.leaf_count;
    return node_index;
  }

  octree.nodes[node_index].minimum_target_size = std::numeric_limits<double>::infinity();
  for(std::uint32_t offset = 0U; offset < count; ++offset) {
    const auto source_index = octree.source_indices[first + offset];
    if(source_index >= octree.seed_sources.size()) {
      continue;
    }
    octree.nodes[node_index].minimum_target_size = std::min(
      octree.nodes[node_index].minimum_target_size,
      octree.seed_sources[source_index].target_size
    );
  }

  const double extent_x = bounds_max[0] - bounds_min[0];
  const double extent_y = bounds_max[1] - bounds_min[1];
  const double extent_z = bounds_max[2] - bounds_min[2];
  const double maximum_extent = std::max({extent_x, extent_y, extent_z});
  const double minimum_extent =
    std::max(field.minimum_length * kBackgroundOctreeMinimumExtentScale, 1.0e-9);

  if(count <= kBackgroundOctreeMaximumLeafSourceCount ||
     depth >= kBackgroundOctreeMaximumDepth ||
     maximum_extent <= minimum_extent) {
    ++octree.leaf_count;
    return node_index;
  }

  std::array<std::vector<std::uint32_t>, 8U> child_sources;
  for(std::uint32_t offset = 0U; offset < count; ++offset) {
    const auto source_index = octree.source_indices[first + offset];
    if(source_index >= octree.seed_sources.size()) {
      continue;
    }
    const auto child_index = background_octree_child_index(
      bounds_min,
      bounds_max,
      octree.seed_sources[source_index].position
    );
    child_sources[child_index].push_back(source_index);
  }

  std::size_t non_empty_child_count = 0U;
  for(const auto &sources : child_sources) {
    if(!sources.empty()) {
      ++non_empty_child_count;
    }
  }
  if(non_empty_child_count == 0U) {
    ++octree.leaf_count;
    return node_index;
  }

  std::uint32_t write_offset = first;
  for(const auto &sources : child_sources) {
    for(const auto source_index : sources) {
      octree.source_indices[write_offset++] = source_index;
    }
  }

  octree.nodes[node_index].leaf = false;
  std::uint32_t child_first = first;
  for(std::uint8_t child_index = 0U; child_index < child_sources.size(); ++child_index) {
    const auto child_count = static_cast<std::uint32_t>(child_sources[child_index].size());
    if(child_count == 0U) {
      continue;
    }

    geo::Point3 child_min {0.0, 0.0, 0.0};
    geo::Point3 child_max {0.0, 0.0, 0.0};
    background_octree_child_bounds(
      bounds_min,
      bounds_max,
      child_index,
      child_min,
      child_max
    );
    octree.nodes[node_index].children[child_index] = build_background_octree_node(
      octree,
      field,
      child_first,
      child_count,
      child_min,
      child_max,
      depth + 1U
    );
    child_first += child_count;
  }

  return node_index;
}

void configure_background_octree(
  const AutoCfdSurfaceSizingFieldState &field,
  const MutableBackgroundBounds &support_bounds,
  AutoCfdSurfaceBackgroundOctree &octree
)
{
  octree = {};
  if(!support_bounds.defined) {
    return;
  }

  octree.minimum = support_bounds.minimum;
  octree.maximum = support_bounds.maximum;

  const auto append_sources =
    [&](const std::vector<AutoCfdSurfaceSizingSource> &sources) {
      for(const auto &source : sources) {
        if(!std::isfinite(source.target_size) || source.target_size <= 0.0) {
          continue;
        }
        octree.seed_sources.push_back(source);
      }
    };
  append_sources(field.curvature_sources);
  append_sources(field.proximity_sources);
  if(octree.seed_sources.empty()) {
    return;
  }

  octree.source_indices.resize(octree.seed_sources.size());
  for(std::uint32_t index = 0U; index < octree.source_indices.size(); ++index) {
    octree.source_indices[index] = index;
  }

  geo::Point3 bounds_min = support_bounds.minimum;
  geo::Point3 bounds_max = support_bounds.maximum;
  const double padding = std::max(field.minimum_length, 1.0e-6);
  for(std::size_t axis = 0U; axis < 3U; ++axis) {
    bounds_min[axis] -= padding;
    bounds_max[axis] += padding;
  }

  static_cast<void>(build_background_octree_node(
    octree,
    field,
    0U,
    static_cast<std::uint32_t>(octree.source_indices.size()),
    bounds_min,
    bounds_max,
    0U
  ));
  octree.built = !octree.nodes.empty();
}

[[nodiscard]] double trilinear_interpolate(
  double c000,
  double c100,
  double c010,
  double c110,
  double c001,
  double c101,
  double c011,
  double c111,
  double tx,
  double ty,
  double tz
) noexcept
{
  const double c00 = c000 * (1.0 - tx) + c100 * tx;
  const double c10 = c010 * (1.0 - tx) + c110 * tx;
  const double c01 = c001 * (1.0 - tx) + c101 * tx;
  const double c11 = c011 * (1.0 - tx) + c111 * tx;
  const double c0 = c00 * (1.0 - ty) + c10 * ty;
  const double c1 = c01 * (1.0 - ty) + c11 * ty;
  return c0 * (1.0 - tz) + c1 * tz;
}

[[nodiscard]] double curve_length_tolerance(double lhs_length, double rhs_length) noexcept
{
  const double reference = std::max({1.0, std::abs(lhs_length), std::abs(rhs_length)});
  return 1.0e-6 * reference;
}

[[nodiscard]] double distortion_angle_radians(double distortion_angle_degrees) noexcept
{
  return distortion_angle_degrees * (kPi / 180.0);
}

// Centroid factor: msize = radius × centroid_factor × crit_ang,
// where crit_ang = distortion_angle in radians. Produces target sizes
// 1.5× larger than the 2*sin(θ/2)/κ formula, allowing coarser elements
// on gently curved surfaces.
constexpr double kCurvatureCentroidFactor = 1.5;

[[nodiscard]] double curvature_to_target_size(
  double curvature,
  double distortion_angle_degrees
) noexcept
{
  if(!std::isfinite(curvature) || curvature <= kCurvatureTolerance) {
    return std::numeric_limits<double>::infinity();
  }

  const double theta = distortion_angle_radians(distortion_angle_degrees);
  return kCurvatureCentroidFactor * theta / curvature;
}

[[nodiscard]] double face_sample_fraction(
  std::size_t index,
  std::size_t sample_count
) noexcept
{
  return static_cast<double>(index + 1U) / static_cast<double>(sample_count + 1U);
}

[[nodiscard]] double squared_value(double value) noexcept
{
  return value * value;
}

[[nodiscard]] std::array<double, 2> metric_eigenvalues(
  double a,
  double b,
  double c
) noexcept
{
  const double trace = a + c;
  const double determinant = a * c - b * b;
  const double discriminant = std::sqrt(std::max(0.0, trace * trace - 4.0 * determinant));
  const double first = 0.5 * (trace + discriminant);
  const double second = 0.5 * (trace - discriminant);
  return {
    std::min(first, second),
    std::max(first, second),
  };
}

[[nodiscard]] std::array<double, 2> metric_primary_eigenvector(
  double a,
  double b,
  double c,
  double eigenvalue
) noexcept
{
  std::array<double, 2> vector {0.0, 0.0};
  if(std::abs(b) > kMetricDeterminantTolerance) {
    vector = {b, eigenvalue - a};
  }
  else if(a >= c) {
    vector = {1.0, 0.0};
  }
  else {
    vector = {0.0, 1.0};
  }

  const double length =
    std::sqrt(squared_value(vector[0]) + squared_value(vector[1]));
  if(length <= 0.0 || !std::isfinite(length)) {
    return {1.0, 0.0};
  }

  return {
    vector[0] / length,
    vector[1] / length,
  };
}

[[nodiscard]] std::array<double, 3> reconstruct_metric_tensor(
  const std::array<double, 2> &primary_eigenvector,
  double primary_eigenvalue,
  double secondary_eigenvalue
) noexcept
{
  const double qx = primary_eigenvector[0];
  const double qy = primary_eigenvector[1];
  return {
    primary_eigenvalue * qx * qx + secondary_eigenvalue * qy * qy,
    (primary_eigenvalue - secondary_eigenvalue) * qx * qy,
    primary_eigenvalue * qy * qy + secondary_eigenvalue * qx * qx,
  };
}

[[nodiscard]] std::uint64_t pack_entity_id(geo::TopologyEntityId entity) noexcept
{
  return (static_cast<std::uint64_t>(static_cast<std::uint8_t>(entity.dimension)) << 32U) |
         static_cast<std::uint64_t>(entity.index);
}

[[nodiscard]] bool should_skip_face_curvature_status(base::StatusCode status) noexcept
{
  return status == base::StatusCode::unsupported ||
         status == base::StatusCode::invalid_argument;
}

[[nodiscard]] bool supports_auto_cfd_surface_seam_unwrap(
  const geo::FaceBoundaryLoops &boundary
) noexcept
{
  return boundary.loops.size() == 1U &&
         boundary.outer_loop_count == 1U &&
         boundary.inner_loop_count == 0U &&
         boundary.unknown_loop_count == 0U &&
         boundary.closed_loop_count == 1U &&
         boundary.open_loop_count == 0U &&
         boundary.non_continuous_loop_count == 0U &&
         boundary.seam_loop_count == 1U &&
         boundary.seam_edge_use_count == 2U &&
         boundary.degenerate_loop_count == 0U &&
         boundary.degenerate_edge_use_count == 0U &&
         boundary.repeated_edge_loop_count == 1U &&
         boundary.repeated_edge_use_count == 1U;
}

[[nodiscard]] double fallback_metric_scale_from_first_form(
  const std::array<double, 3> &first_form,
  double target_size
) noexcept
{
  const double target_size_squared = std::max(
    squared_value(target_size),
    kMetricEigenvalueFloor
  );
  const double average_scale = 0.5 * (std::max(0.0, first_form[0]) + std::max(0.0, first_form[2]));
  if(std::isfinite(average_scale) && average_scale > kMetricEigenvalueFloor) {
    return std::clamp(
      average_scale / target_size_squared,
      kMetricEigenvalueFloor,
      kMetricEigenvalueCeiling
    );
  }

  return std::clamp(
    1.0 / target_size_squared,
    kMetricEigenvalueFloor,
    kMetricEigenvalueCeiling
  );
}

void build_metric_fallback_from_position(
  geo::TopologyEntityId face,
  const AutoCfdSurfaceSizingFieldState &field,
  const geo::Point3 &position,
  double u,
  double v,
  const std::array<double, 3> *first_form,
  AutoCfdSurfaceMetricFallbackKind fallback_kind,
  AutoCfdSurfaceFaceMetricTensor &metric
)
{
  metric = {};
  metric.face = face;
  metric.u = u;
  metric.v = v;
  metric.target_size = query_auto_cfd_surface_sizing_field(field, face, position);
  if(first_form != nullptr) {
    metric.first_fundamental_form = *first_form;
    metric.first_fundamental_form_defined = true;
  }

  const double eigenvalue =
    fallback_metric_scale_from_first_form(
      metric.first_fundamental_form,
      metric.target_size
    );
  metric.tensor = {eigenvalue, 0.0, eigenvalue};
  metric.raw_eigenvalues = {0.0, 0.0};
  metric.clamped_eigenvalues = {eigenvalue, eigenvalue};
  metric.determinant = eigenvalue * eigenvalue;
  metric.usable = true;
  metric.fallback_kind = fallback_kind;
}

void unwrap_periodic_loop_uvs(
  const geo::FaceUvBounds &bounds,
  std::vector<AutoCfdSurfaceFaceLoopPoint> &points
) noexcept
{
  if(points.size() < 2U) {
    return;
  }

  const std::array<double, 2> spans {
    bounds.u_max - bounds.u_min,
    bounds.v_max - bounds.v_min,
  };

  for(std::size_t dimension = 0U; dimension < spans.size(); ++dimension) {
    const double span = spans[dimension];
    if(!std::isfinite(span) || span <= kPeriodicUvWrapTolerance) {
      continue;
    }

    double previous = points.front().uv[dimension];
    for(std::size_t point_index = 1U; point_index < points.size(); ++point_index) {
      const double wrap_count =
        std::round((points[point_index].uv[dimension] - previous) / span);
      points[point_index].uv[dimension] -= wrap_count * span;
      previous = points[point_index].uv[dimension];
    }
  }
}

[[nodiscard]] base::StatusCode append_recovered_loop_point(
  const AutoCfdSurfaceBoundaryState &boundary_state,
  const AutoCfdSurfaceEdgePoint &edge_point,
  std::vector<AutoCfdSurfaceFaceLoopPoint> &points
)
{
  if(edge_point.node_index >= boundary_state.nodes.size()) {
    return core::detail::publish_error(
      base::StatusCode::internal_error,
      "Auto CFD face preprocessing encountered an out-of-range boundary node index."
    );
  }

  AutoCfdSurfaceFaceLoopPoint point;
  point.node_index = edge_point.node_index;
  point.parameter = edge_point.parameter;
  point.boundary_target_size = edge_point.target_size;
  point.position = boundary_state.nodes[edge_point.node_index].position;
  points.push_back(std::move(point));
  return core::detail::clear_error_state();
}

[[nodiscard]] base::StatusCode append_oriented_edge_discretization_points(
  const AutoCfdSurfaceBoundaryState &boundary_state,
  const AutoCfdSurfaceEdgeDiscretization &edge_discretization,
  const geo::FaceBoundaryEdgeUse &edge_use,
  bool include_terminal_point,
  std::vector<AutoCfdSurfaceFaceLoopPoint> &points
)
{
  if(edge_discretization.points.empty()) {
    return core::detail::clear_error_state();
  }

  if(edge_use.same_orientation_as_edge) {
    const std::size_t limit =
      include_terminal_point
        ? edge_discretization.points.size()
        : edge_discretization.points.size() - 1U;
    for(std::size_t point_index = 0U; point_index < limit; ++point_index) {
      const auto status = append_recovered_loop_point(
        boundary_state,
        edge_discretization.points[point_index],
        points
      );
      if(status != base::StatusCode::ok) {
        return status;
      }
    }
    return core::detail::clear_error_state();
  }

  const std::size_t floor = include_terminal_point ? 0U : 1U;
  for(std::size_t point_index = edge_discretization.points.size();
      point_index-- > floor;) {
    const auto status = append_recovered_loop_point(
      boundary_state,
      edge_discretization.points[point_index],
      points
    );
    if(status != base::StatusCode::ok) {
      return status;
    }
  }

  return core::detail::clear_error_state();
}

[[nodiscard]] base::StatusCode recover_face_loop_from_boundary_dictionary(
  const AutoCfdSurfaceBoundaryState &boundary_state,
  const geo::FaceBoundaryLoop &loop,
  AutoCfdSurfaceFaceLoopState &recovered_loop,
  bool face_has_seams = false
)
{
  recovered_loop = {};
  recovered_loop.kind = loop.kind;
  recovered_loop.closed = loop.closed;
  recovered_loop.continuous = loop.continuous;
  recovered_loop.seam_edge_use_count = loop.seam_edge_use_count;
  recovered_loop.degenerate_edge_use_count = loop.degenerate_edge_use_count;
  recovered_loop.segments.reserve(loop.edge_uses.size());

  for(std::size_t edge_use_index = 0U; edge_use_index < loop.edge_uses.size(); ++edge_use_index) {
    const auto &edge_use = loop.edge_uses[edge_use_index];
    if(edge_use.edge.dimension != geo::TopologyDimension::edge ||
       edge_use.edge.index >= boundary_state.edge_discretizations.size()) {
      return core::detail::publish_error(
        base::StatusCode::internal_error,
        "Auto CFD face preprocessing could not resolve one of the recovered boundary edge discretizations."
      );
    }

    const auto &edge_discretization =
      boundary_state.edge_discretizations[edge_use.edge.index];
    AutoCfdSurfaceFaceLoopSegment segment;
    segment.edge_use = edge_use;
    
    const bool include_terminal_point =
      face_has_seams ||
      (!loop.closed && edge_use_index + 1U == loop.edge_uses.size());

    auto corrected_edge_use = edge_use;
    if(geo::is_valid(edge_use.start_vertex) &&
       edge_discretization.points.size() >= 2U) {
      const auto first_node_idx = edge_discretization.points.front().node_index;
      const auto last_node_idx  = edge_discretization.points.back().node_index;
      if(first_node_idx < boundary_state.nodes.size() &&
         last_node_idx  < boundary_state.nodes.size()) {
        const auto first_topo_vertex =
          boundary_state.nodes[first_node_idx].topology_vertex;
        const auto last_topo_vertex =
          boundary_state.nodes[last_node_idx].topology_vertex;
        if(edge_use.start_vertex == first_topo_vertex) {
          corrected_edge_use.same_orientation_as_edge = true;
        } else if(edge_use.start_vertex == last_topo_vertex) {
          corrected_edge_use.same_orientation_as_edge = false;
        }
        // else: degenerate or seam edge — keep original flag untouched.
      }
    }

    auto status = append_oriented_edge_discretization_points(
      boundary_state,
      edge_discretization,
      corrected_edge_use,
      include_terminal_point,
      segment.points
    );
    if(status != base::StatusCode::ok) {
      return status;
    }

    recovered_loop.points.insert(
      recovered_loop.points.end(),
      segment.points.begin(),
      segment.points.end()
    );
    recovered_loop.segments.push_back(std::move(segment));
  }

  return core::detail::clear_error_state();
}

void populate_recovered_face_uv_extents(
  AutoCfdSurfaceFacePreprocessState &face_preprocess
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

  for(const auto &loop : face_preprocess.loops) {
    for(const auto &point : loop.points) {
      if(!point.uv_defined) {
        continue;
      }
      minima[0] = std::min(minima[0], point.uv[0]);
      minima[1] = std::min(minima[1], point.uv[1]);
      maxima[0] = std::max(maxima[0], point.uv[0]);
      maxima[1] = std::max(maxima[1], point.uv[1]);
    }
  }

  if(std::isfinite(minima[0]) && std::isfinite(minima[1]) &&
     std::isfinite(maxima[0]) && std::isfinite(maxima[1])) {
    face_preprocess.recovered_uv_min = minima;
    face_preprocess.recovered_uv_max = maxima;
  }
}

[[nodiscard]] bool interpolate_missing_segment_uvs(
  std::vector<AutoCfdSurfaceFaceLoopPoint> &points
) noexcept
{
  if(points.empty()) {
    return true;
  }

  std::vector<std::size_t> defined_indices;
  defined_indices.reserve(points.size());
  for(std::size_t point_index = 0U; point_index < points.size(); ++point_index) {
    if(points[point_index].uv_defined) {
      defined_indices.push_back(point_index);
    }
  }

  if(defined_indices.size() == points.size()) {
    return true;
  }
  if(defined_indices.empty()) {
    return false;
  }

  const auto assign_interpolated_uv =
    [&points](
      std::size_t target_index,
      std::size_t first_index,
      std::size_t second_index
    ) noexcept {
      const double first_parameter = points[first_index].parameter;
      const double second_parameter = points[second_index].parameter;
      const double parameter_span = second_parameter - first_parameter;

      if(std::abs(parameter_span) <= 1.0e-12) {
        points[target_index].uv = points[first_index].uv;
        points[target_index].uv_defined = true;
        return;
      }

      const double t =
        (points[target_index].parameter - first_parameter) / parameter_span;
      points[target_index].uv = {
        points[first_index].uv[0] +
          (points[second_index].uv[0] - points[first_index].uv[0]) * t,
        points[first_index].uv[1] +
          (points[second_index].uv[1] - points[first_index].uv[1]) * t,
      };
      points[target_index].uv_defined = true;
    };

  for(std::size_t point_index = 0U; point_index < points.size(); ++point_index) {
    if(points[point_index].uv_defined) {
      continue;
    }

    const auto upper_it =
      std::lower_bound(defined_indices.begin(), defined_indices.end(), point_index);
    const bool has_upper = upper_it != defined_indices.end();
    const bool has_lower = upper_it != defined_indices.begin();

    if(has_lower && has_upper) {
      assign_interpolated_uv(point_index, *(upper_it - 1U), *upper_it);
      continue;
    }

    if(defined_indices.size() == 1U) {
      points[point_index].uv = points[defined_indices.front()].uv;
      points[point_index].uv_defined = true;
      continue;
    }

    if(has_upper) {
      assign_interpolated_uv(point_index, defined_indices[0], defined_indices[1]);
      continue;
    }

    assign_interpolated_uv(
      point_index,
      defined_indices[defined_indices.size() - 2U],
      defined_indices.back()
    );
  }

  return std::all_of(
    points.begin(),
    points.end(),
    [](const AutoCfdSurfaceFaceLoopPoint &point) noexcept {
      return point.uv_defined;
    }
  );
}

[[nodiscard]] double edge_parameter_tolerance(
  const geo::EdgeCurveInfo &curve
) noexcept
{
  const double parameter_range =
    std::max(1.0, std::abs(curve.parameter_max - curve.parameter_min));
  return parameter_range * kEdgeParameterToleranceScale;
}

[[nodiscard]] geo::TopologyEntityId edge_endpoint_vertex(
  const AutoCfdSurfaceCurveWorkItem &curve_work_item,
  bool start
) noexcept
{
  if(curve_work_item.vertex_ids.size() < 2U) {
    return {};
  }

  return start ? curve_work_item.vertex_ids.front() : curve_work_item.vertex_ids.back();
}

[[nodiscard]] const AutoCfdSurfaceAnchorCandidate *find_anchor_candidate(
  const std::vector<AutoCfdSurfaceAnchorCandidate> &anchor_candidates,
  geo::TopologyEntityId vertex
) noexcept
{
  if(vertex.dimension != geo::TopologyDimension::vertex ||
     vertex.index >= anchor_candidates.size()) {
    return nullptr;
  }

  const auto &candidate = anchor_candidates[vertex.index];
  if(candidate.vertex != vertex) {
    return nullptr;
  }

  return &candidate;
}

[[nodiscard]] std::uint32_t add_boundary_node(
  AutoCfdSurfaceBoundaryState &boundary_state,
  const geo::Point3 &position,
  geo::TopologyEntityId topology_vertex,
  const AutoCfdSurfaceAnchorCandidate *anchor_candidate
)
{
  if(topology_vertex.dimension == geo::TopologyDimension::vertex &&
     topology_vertex.index < boundary_state.vertex_node_indices.size()) {
    auto &node_index = boundary_state.vertex_node_indices[topology_vertex.index];
    if(node_index != invalid_index) {
      return node_index;
    }

    node_index = static_cast<std::uint32_t>(boundary_state.nodes.size());
  }

  AutoCfdSurfaceBoundaryNode node;
  node.position = position;
  node.topology_vertex = topology_vertex;
  node.topology_endpoint = geo::is_valid(topology_vertex);
  if(anchor_candidate != nullptr) {
    node.anchor_source = anchor_candidate->source;
    node.frozen_anchor = anchor_candidate->frozen;
  }

  boundary_state.nodes.push_back(std::move(node));

  if(topology_vertex.dimension == geo::TopologyDimension::vertex &&
     topology_vertex.index < boundary_state.vertex_node_indices.size()) {
    return boundary_state.vertex_node_indices[topology_vertex.index];
  }

  return static_cast<std::uint32_t>(boundary_state.nodes.size() - 1U);
}

[[nodiscard]] base::StatusCode sample_edge_discretization_point(
  const geo::EdgeView &edge_view,
  geo::TopologyEntityId owner,
  const AutoCfdSurfaceSizingFieldState &field,
  double parameter,
  EdgeDiscretizationSample &sample
)
{
  geo::EdgeTangentSample tangent_sample;
  const auto status = geo::sample_edge_tangent(edge_view, parameter, tangent_sample);
  if(status != base::StatusCode::ok) {
    sample = {};
    return status;
  }

  sample.parameter = parameter;
  sample.position = tangent_sample.position;
  sample.target_size = query_auto_cfd_surface_sizing_field(
    field,
    owner,
    tangent_sample.position
  );
  sample.speed = std::max(0.0, tangent_sample.speed);
  if(std::getenv("SQMESH_SIZE_PROBE") != nullptr) {
    std::printf("[SIZE] edge=%u t=%.10g xyz=(%.10g %.10g %.10g) size=%.10g\n",
                owner.index,
                parameter,
                tangent_sample.position[0],
                tangent_sample.position[1],
                tangent_sample.position[2],
                sample.target_size);
  }
  return core::detail::clear_error_state();
}

[[nodiscard]] bool should_use_short_edge_fallback(
  const geo::EdgeCurveInfo &curve,
  const AutoCfdSurfaceSizingFieldState &field
) noexcept
{
  const double tolerance = curve_length_tolerance(
    curve.approximate_length,
    field.minimum_length
  );

  return !std::isfinite(curve.approximate_length) ||
         curve.approximate_length <= std::max(field.minimum_length, tolerance) ||
         curve.parameter_max - curve.parameter_min <= edge_parameter_tolerance(curve);
}

[[nodiscard]] double edge_discretization_density(
  const EdgeDiscretizationSample &sample
) noexcept
{
  if(!std::isfinite(sample.speed) || sample.speed <= 0.0 ||
     !std::isfinite(sample.target_size) || sample.target_size <= 0.0) {
    return 0.0;
  }

  return sample.speed / sample.target_size;
}

[[nodiscard]] double edge_discretization_trapezoidal_integral(
  const EdgeDiscretizationSample &start,
  const EdgeDiscretizationSample &end
) noexcept
{
  return 0.5 * (edge_discretization_density(start) + edge_discretization_density(end)) *
         (end.parameter - start.parameter);
}

[[nodiscard]] base::StatusCode integrate_edge_discretization_primitive(
  const geo::EdgeView &edge_view,
  geo::TopologyEntityId owner,
  const geo::EdgeCurveInfo &curve,
  const AutoCfdSurfaceSizingFieldState &field,
  const EdgeDiscretizationSample &start,
  const EdgeDiscretizationSample &end,
  std::size_t depth,
  std::vector<EdgeDiscretizationSample> &samples
)
{
  const double parameter_span = end.parameter - start.parameter;
  if(parameter_span <= edge_parameter_tolerance(curve)) {
    EdgeDiscretizationSample end_sample = end;
    end_sample.primitive =
      samples.empty() ? 0.0 :
                        samples.back().primitive +
                          edge_discretization_trapezoidal_integral(start, end);
    samples.push_back(end_sample);
    return core::detail::clear_error_state();
  }

  const double midpoint_parameter = 0.5 * (start.parameter + end.parameter);
  if(midpoint_parameter <= start.parameter || midpoint_parameter >= end.parameter) {
    EdgeDiscretizationSample end_sample = end;
    end_sample.primitive =
      samples.empty() ? 0.0 :
                        samples.back().primitive +
                          edge_discretization_trapezoidal_integral(start, end);
    samples.push_back(end_sample);
    return core::detail::clear_error_state();
  }

  EdgeDiscretizationSample midpoint;
  auto status = sample_edge_discretization_point(
    edge_view,
    owner,
    field,
    midpoint_parameter,
    midpoint
  );
  if(status != base::StatusCode::ok) {
    return status;
  }

  const double geometry_tolerance =
    curve_length_tolerance(curve.approximate_length, point_distance(start.position, end.position));
  if(point_distance(start.position, midpoint.position) <= geometry_tolerance ||
     point_distance(midpoint.position, end.position) <= geometry_tolerance) {
    EdgeDiscretizationSample end_sample = end;
    end_sample.primitive =
      samples.empty() ? 0.0 :
                        samples.back().primitive +
                          edge_discretization_trapezoidal_integral(start, end);
    samples.push_back(end_sample);
    return core::detail::clear_error_state();
  }

  const double coarse_value = edge_discretization_trapezoidal_integral(start, end);
  const double left_value = edge_discretization_trapezoidal_integral(start, midpoint);
  const double right_value = edge_discretization_trapezoidal_integral(midpoint, end);
  const double integration_error = std::abs(coarse_value - left_value - right_value);

  if(((integration_error <= kEdgeIntegrationPrecision) &&
      (depth >= kMinimumEdgeIntegrationDepth)) ||
     depth >= kMaximumEdgeIntegrationDepth) {
    if(samples.size() + 2U >= kMaximumEdgeDiscretizationPoints) {
      EdgeDiscretizationSample end_sample = end;
      end_sample.primitive =
        samples.empty() ? 0.0 : samples.back().primitive + coarse_value;
      samples.push_back(end_sample);
      return core::detail::clear_error_state();
    }

    EdgeDiscretizationSample midpoint_sample = midpoint;
    midpoint_sample.primitive =
      samples.empty() ? left_value : samples.back().primitive + left_value;
    samples.push_back(midpoint_sample);

    EdgeDiscretizationSample end_sample = end;
    end_sample.primitive = samples.back().primitive + right_value;
    samples.push_back(end_sample);
    return core::detail::clear_error_state();
  }

  status = integrate_edge_discretization_primitive(
    edge_view,
    owner,
    curve,
    field,
    start,
    midpoint,
    depth + 1U,
    samples
  );
  if(status != base::StatusCode::ok) {
    return status;
  }

  return integrate_edge_discretization_primitive(
    edge_view,
    owner,
    curve,
    field,
    midpoint,
    end,
    depth + 1U,
    samples
  );
}


double smooth_edge_discretization_primitive(
  double alpha,
  std::vector<EdgeDiscretizationSample> &samples
) noexcept
{
  if(samples.size() < 2U) {
    return samples.empty() ? 0.0 : samples.back().primitive;
  }

  std::size_t iteration = 0U;
  while(true) {
    std::size_t count1 = 0U;
    std::size_t count2 = 0U;

    for(std::size_t i = 1U; i < samples.size(); ++i) {
      auto &prev = samples[i - 1U];
      auto &curr = samples[i];
      if(prev.target_size <= 0.0 || curr.target_size <= 0.0 ||
         curr.speed <= 0.0) {
        continue;
      }
      const double dt = curr.parameter - prev.parameter;
      if(dt <= 0.0) {
        continue;
      }
      const double dh = curr.target_size - prev.target_size;
      const double dhdt = dh / dt;
      if(dhdt / curr.speed > (alpha - 1.0) * 1.01) {
        const double hnew =
          prev.target_size + dt * (alpha - 1.0) * curr.speed;
        curr.target_size = hnew;
        ++count1;
      }
    }

    for(std::size_t i = samples.size() - 1U; i > 0U; --i) {
      auto &prev = samples[i - 1U];
      auto &curr = samples[i];
      if(prev.target_size <= 0.0 || curr.target_size <= 0.0 ||
         prev.speed <= 0.0) {
        continue;
      }
      const double dt = std::abs(prev.parameter - curr.parameter);
      if(dt <= 0.0) {
        continue;
      }
      const double dh = prev.target_size - curr.target_size;
      const double dhdt = dh / dt;
      if(dhdt / prev.speed > (alpha - 1.0) * 1.01) {
        const double hnew =
          curr.target_size + dt * (alpha - 1.0) * prev.speed;
        prev.target_size = hnew;
        ++count2;
      }
    }

    ++iteration;
    if(iteration > kEdgeSizeSmoothingMaxIterations) {
      break;
    }
    if((count1 + count2) == 0U) {
      break;
    }
  }

  samples.front().primitive = 0.0;
  for(std::size_t i = 1U; i < samples.size(); ++i) {
    auto &prev = samples[i - 1U];
    auto &curr = samples[i];
    curr.primitive =
      prev.primitive +
      (curr.parameter - prev.parameter) * 0.5 *
        (edge_discretization_density(prev) + edge_discretization_density(curr));
  }
  return samples.back().primitive;
}

[[nodiscard]] base::StatusCode distribute_edge_discretization_samples(
  const geo::EdgeView &edge_view,
  geo::TopologyEntityId owner,
  const geo::EdgeCurveInfo &curve,
  const AutoCfdSurfaceSizingFieldState &field,
  const EdgeDiscretizationSample &start,
  const EdgeDiscretizationSample &end,
  std::vector<EdgeDiscretizationSample> &samples
)
{
  std::vector<EdgeDiscretizationSample> primitive_samples;
  primitive_samples.reserve(32U);

  EdgeDiscretizationSample primitive_start = start;
  primitive_start.primitive = 0.0;
  primitive_samples.push_back(primitive_start);

  auto status = integrate_edge_discretization_primitive(
    edge_view,
    owner,
    curve,
    field,
    primitive_start,
    end,
    0U,
    primitive_samples
  );
  if(status != base::StatusCode::ok) {
    return status;
  }
  if(primitive_samples.size() < 2U) {
    samples.push_back(end);
    return core::detail::clear_error_state();
  }

 
  smooth_edge_discretization_primitive(
    std::sqrt(kEdgeSizeSmoothingRatio),
    primitive_samples
  );

  const double total_primitive = primitive_samples.back().primitive;
  if(!std::isfinite(total_primitive) || total_primitive <= 1.0) {
    samples.push_back(end);
    return core::detail::clear_error_state();
  }

  std::size_t segment_count =
    std::max<std::size_t>(1U, static_cast<std::size_t>(total_primitive + 0.99));
  segment_count =
    std::min<std::size_t>(segment_count, kMaximumEdgeDiscretizationPoints - 1U);
  if(segment_count <= 1U) {
    samples.push_back(end);
    return core::detail::clear_error_state();
  }

  const double primitive_step =
    total_primitive / static_cast<double>(segment_count);
  std::size_t primitive_index = 1U;

  for(std::size_t segment_index = 1U; segment_index < segment_count; ++segment_index) {
    const double target_primitive = primitive_step * static_cast<double>(segment_index);
    while(primitive_index < primitive_samples.size() &&
          primitive_samples[primitive_index].primitive < target_primitive) {
      ++primitive_index;
    }
    if(primitive_index >= primitive_samples.size()) {
      break;
    }

    const auto &lhs = primitive_samples[primitive_index - 1U];
    const auto &rhs = primitive_samples[primitive_index];
    const double primitive_span = rhs.primitive - lhs.primitive;
    if(!std::isfinite(primitive_span) || primitive_span <= 0.0) {
      continue;
    }

    const double alpha = std::clamp(
      (target_primitive - lhs.primitive) / primitive_span,
      0.0,
      1.0
    );
    const double parameter =
      lhs.parameter + alpha * (rhs.parameter - lhs.parameter);
    if(parameter <= samples.back().parameter + edge_parameter_tolerance(curve) ||
       parameter >= end.parameter - edge_parameter_tolerance(curve)) {
      continue;
    }

    EdgeDiscretizationSample sample;
    status = sample_edge_discretization_point(
      edge_view,
      owner,
      field,
      parameter,
      sample
    );
    if(status != base::StatusCode::ok) {
      return status;
    }

    if(point_distance(samples.back().position, sample.position) <=
       curve_length_tolerance(curve.approximate_length, field.minimum_length)) {
      continue;
    }

    samples.push_back(sample);
  }

  samples.push_back(end);

  return core::detail::clear_error_state();
}


void filter_edge_discretization_samples(
  const geo::EdgeView &edge_view,
  geo::TopologyEntityId owner,
  const AutoCfdSurfaceSizingFieldState &field,
  std::size_t minimum_interior_count,
  std::vector<EdgeDiscretizationSample> &samples
)
{
  if(samples.size() <= 2U) {
    return;
  }

  struct FilterCandidate final {
    double key = 0.0;
    std::size_t index = 0U;
  };

  std::vector<FilterCandidate> candidates;
  candidates.reserve(samples.size());

  std::size_t anchor_index = 0U;
  for(std::size_t i = 1U; i + 1U < samples.size(); ++i) {
    const auto &anchor = samples[anchor_index];
    const auto &current = samples[i];
    const double d = point_distance(anchor.position, current.position);

    const double mid_parameter = 0.5 * (anchor.parameter + current.parameter);
    EdgeDiscretizationSample mid_sample;
    double lc = 0.0;
    if(sample_edge_discretization_point(
         edge_view, owner, field, mid_parameter, mid_sample) ==
       base::StatusCode::ok) {
      lc = mid_sample.target_size;
    }
    if(lc <= 0.0 || !std::isfinite(lc)) {
      lc = current.target_size;
    }
    if(lc <= 0.0 || !std::isfinite(lc)) {
      anchor_index = i;
      continue;
    }

    if(d < lc * kEdgeFilterCloseRatio) {
      candidates.push_back({lc / d, i});
    }
    else {
      anchor_index = i;
    }
  }

  if(candidates.empty()) {
    return;
  }

  std::sort(
    candidates.begin(),
    candidates.end(),
    [](const FilterCandidate &a, const FilterCandidate &b) noexcept {
      return a.key < b.key;
    }
  );

  const std::size_t interior_count = samples.size() - 2U;
  if(interior_count < candidates.size() ||
     (interior_count - candidates.size()) < minimum_interior_count) {
    return;
  }

  std::vector<std::size_t> removal_indices;
  removal_indices.reserve(candidates.size());
  for(const auto &candidate : candidates) {
    removal_indices.push_back(candidate.index);
  }
  std::sort(
    removal_indices.begin(),
    removal_indices.end(),
    std::greater<std::size_t>()
  );
  for(const std::size_t index : removal_indices) {
    samples.erase(samples.begin() + static_cast<std::ptrdiff_t>(index));
  }
}

[[nodiscard]] geo::TopologyEntityId resolve_edge_endpoint_vertex(
  const AutoCfdSurfaceCurveWorkItem &curve_work_item,
  const AutoCfdSurfaceBoundaryState &boundary_state,
  const geo::Point3 &position
) noexcept
{
  geo::TopologyEntityId best_vertex;
  double best_distance = std::numeric_limits<double>::max();

  for(const auto vertex_id : curve_work_item.vertex_ids) {
    if(vertex_id.dimension != geo::TopologyDimension::vertex ||
       vertex_id.index >= boundary_state.vertex_node_indices.size()) {
      continue;
    }
    const auto node_index = boundary_state.vertex_node_indices[vertex_id.index];
    if(node_index == invalid_index || node_index >= boundary_state.nodes.size()) {
      continue;
    }
    
    const double distance = point_distance(position, boundary_state.nodes[node_index].position);
    if(distance < best_distance) {
      best_distance = distance;
      best_vertex = vertex_id;
    }
  }

  return best_vertex;
}

[[nodiscard]] base::StatusCode discretize_curve_work_item(
  const geo::ModelView &model_view,
  const std::vector<AutoCfdSurfaceAnchorCandidate> &anchor_candidates,
  const AutoCfdSurfaceSizingFieldState &field,
  const AutoCfdSurfaceCurveWorkItem &curve_work_item,
  AutoCfdSurfaceBoundaryState &boundary_state
)
{
  const auto *edge_view = model_view.find_edge(curve_work_item.edge);
  if(edge_view == nullptr) {
    return core::detail::publish_error(
      base::StatusCode::internal_error,
      "Auto CFD edge discretization could not resolve the requested geometry edge."
    );
  }
  if(curve_work_item.edge.index >= boundary_state.edge_discretizations.size()) {
    return core::detail::publish_error(
      base::StatusCode::internal_error,
      "Auto CFD edge discretization encountered an out-of-range geometry edge index."
    );
  }

  EdgeDiscretizationSample start;
  auto status = sample_edge_discretization_point(
    *edge_view,
    curve_work_item.edge,
    field,
    curve_work_item.curve.parameter_min,
    start
  );
  if(status != base::StatusCode::ok) {
    return status;
  }

  EdgeDiscretizationSample end;
  status = sample_edge_discretization_point(
    *edge_view,
    curve_work_item.edge,
    field,
    curve_work_item.curve.parameter_max,
    end
  );
  if(status != base::StatusCode::ok) {
    return status;
  }

  auto &edge_discretization = boundary_state.edge_discretizations[curve_work_item.edge.index];
  edge_discretization = {};
  edge_discretization.edge = curve_work_item.edge;
  edge_discretization.start_vertex = resolve_edge_endpoint_vertex(curve_work_item, boundary_state, start.position);
  edge_discretization.end_vertex = resolve_edge_endpoint_vertex(curve_work_item, boundary_state, end.position);
  edge_discretization.approximate_length =
    std::max(0.0, curve_work_item.curve.approximate_length);

  std::vector<EdgeDiscretizationSample> samples;
  samples.reserve(8U);
  samples.push_back(start);

  if(should_use_short_edge_fallback(curve_work_item.curve, field)) {
    edge_discretization.short_edge_fallback = true;
    samples.push_back(end);
  }
  else {
    status = distribute_edge_discretization_samples(
      *edge_view,
      curve_work_item.edge,
      curve_work_item.curve,
      field,
      start,
      end,
      samples
    );
    if(status != base::StatusCode::ok) {
      return status;
    }
    
    filter_edge_discretization_samples(
      *edge_view,
      curve_work_item.edge,
      field,
      0U,
      samples
    );
  }

  edge_discretization.points.reserve(samples.size());
  for(std::size_t sample_index = 0U; sample_index < samples.size(); ++sample_index) {
    geo::TopologyEntityId topology_vertex;
    if(sample_index == 0U) {
      topology_vertex = edge_discretization.start_vertex;
    }
    else if(sample_index + 1U == samples.size()) {
      topology_vertex = edge_discretization.end_vertex;
    }

    const auto *anchor_candidate =
      find_anchor_candidate(anchor_candidates, topology_vertex);
    const auto node_index = add_boundary_node(
      boundary_state,
      samples[sample_index].position,
      topology_vertex,
      anchor_candidate
    );
    edge_discretization.points.push_back(
      {
        node_index,
        samples[sample_index].parameter,
        samples[sample_index].target_size,
      }
    );
  }

  return core::detail::clear_error_state();
}

[[nodiscard]] base::StatusCode populate_boundary_state(
  const geo::ModelView &model_view,
  const std::vector<AutoCfdSurfaceAnchorCandidate> &anchor_candidates,
  const std::vector<AutoCfdSurfaceCurveWorkItem> &curve_work_items,
  const std::vector<AutoCfdSurfaceFaceWorkItem> &face_work_items,
  const AutoCfdSurfaceSizingFieldState &field,
  AutoCfdSurfaceBoundaryState &boundary_state
)
{
  boundary_state = {};
  boundary_state.vertex_node_indices.assign(model_view.vertices.size(), invalid_index);
  boundary_state.edge_discretizations.assign(model_view.edges.size(), {});
  boundary_state.nodes.reserve(model_view.vertices.size());

  // ── Pre-register all topology vertices ────────────────────────────────────
  //
  // resolve_edge_endpoint_vertex matches sampled endpoint positions against
  // already-registered boundary nodes to determine which topology vertex
  // corresponds to param_min vs param_max.  If a vertex has not been
  // registered yet (vertex_node_indices[v] == invalid_index), the function
  // skips it and may return an invalid topology_vertex.  add_boundary_node
  // then creates an unregistered node, so subsequent edges sharing the same
  // geometric vertex create duplicates — producing open edges at every
  // shared corner in the assembled mesh.
  //
  // Fix: use the face boundary edge-use data (FaceBoundaryEdgeUse), which
  // carries wire-oriented start_vertex / end_vertex determined by the
  // topology backend (OCC TopExp::Vertices on the oriented edge).  For each
  // edge use we evaluate the curve at param_min and param_max, then register
  // start_vertex and end_vertex with the correct positions.
  //
  // Because multiple face loops reference the same edge, a vertex may be
  // encountered more than once; add_boundary_node deduplicates via
  // vertex_node_indices so only the first registration takes effect.
  std::unordered_map<std::uint32_t, const AutoCfdSurfaceCurveWorkItem *> edge_to_curve;
  edge_to_curve.reserve(curve_work_items.size());
  for(const auto &cwi : curve_work_items) {
    edge_to_curve.emplace(cwi.edge.index, &cwi);
  }

  for(const auto &face_work_item : face_work_items) {
    for(const auto &loop : face_work_item.boundary.loops) {
      for(const auto &edge_use : loop.edge_uses) {
        if(edge_use.edge.dimension != geo::TopologyDimension::edge ||
           edge_use.edge.index >= boundary_state.edge_discretizations.size()) {
          continue;
        }

        // Skip if both vertices are already registered.
        const bool sv_needs = geo::is_valid(edge_use.start_vertex) &&
          edge_use.start_vertex.dimension == geo::TopologyDimension::vertex &&
          edge_use.start_vertex.index < boundary_state.vertex_node_indices.size() &&
          boundary_state.vertex_node_indices[edge_use.start_vertex.index] == invalid_index;
        const bool ev_needs = geo::is_valid(edge_use.end_vertex) &&
          edge_use.end_vertex.dimension == geo::TopologyDimension::vertex &&
          edge_use.end_vertex.index < boundary_state.vertex_node_indices.size() &&
          boundary_state.vertex_node_indices[edge_use.end_vertex.index] == invalid_index;
        if(!sv_needs && !ev_needs) {
          continue;
        }

        const auto *edge_view = model_view.find_edge(edge_use.edge);
        if(edge_view == nullptr) {
          continue;
        }

        auto curve_it = edge_to_curve.find(edge_use.edge.index);
        if(curve_it == edge_to_curve.end()) {
          continue;
        }
        const auto *curve_item = curve_it->second;

        // Sample the edge at param_min and param_max.
        geo::EdgeTangentSample start_sample;
        if(geo::sample_edge_tangent(*edge_view, curve_item->curve.parameter_min, start_sample) !=
           base::StatusCode::ok) {
          continue;
        }
        geo::EdgeTangentSample end_sample;
        if(geo::sample_edge_tangent(*edge_view, curve_item->curve.parameter_max, end_sample) !=
           base::StatusCode::ok) {
          continue;
        }

        // edge_use.start_vertex is the vertex at the START of the oriented
        // edge in the wire (from OCC TopExp::Vertices with cumOri=true).
        // same_orientation_as_edge == true means the wire direction matches
        // the edge parameter direction, so start_vertex → param_min and
        // end_vertex → param_max.  When false, they are swapped.
        if(sv_needs) {
          const auto &sv_pos = edge_use.same_orientation_as_edge
            ? start_sample.position
            : end_sample.position;
          const auto *anchor = find_anchor_candidate(anchor_candidates, edge_use.start_vertex);
          static_cast<void>(add_boundary_node(
            boundary_state,
            sv_pos,
            edge_use.start_vertex,
            anchor
          ));
        }

        if(ev_needs) {
          const auto &ev_pos = edge_use.same_orientation_as_edge
            ? end_sample.position
            : start_sample.position;
          const auto *anchor = find_anchor_candidate(anchor_candidates, edge_use.end_vertex);
          static_cast<void>(add_boundary_node(
            boundary_state,
            ev_pos,
            edge_use.end_vertex,
            anchor
          ));
        }
      }
    }
  }

  for(const auto &curve_work_item : curve_work_items) {
    const auto status = discretize_curve_work_item(
      model_view,
      anchor_candidates,
      field,
      curve_work_item,
      boundary_state
    );
    if(status != base::StatusCode::ok) {
      return status;
    }
  }

  // Diagnostic: print per-edge discretization size range and gradation.
  for(const auto &ed : boundary_state.edge_discretizations) {
    if(ed.points.size() < 2U) {
      continue;
    }
    double min_size = ed.points.front().target_size;
    double max_size = ed.points.front().target_size;
    double max_ratio = 1.0;
    for(std::size_t pi = 0U; pi < ed.points.size(); ++pi) {
      const double s = ed.points[pi].target_size;
      min_size = std::min(min_size, s);
      max_size = std::max(max_size, s);
      if(pi > 0U) {
        const double prev = ed.points[pi - 1U].target_size;
        if(prev > 0.0 && s > 0.0) {
          max_ratio = std::max(max_ratio, std::max(s / prev, prev / s));
        }
      }
    }
    if(max_ratio > 1.5 || max_size / std::max(1.0e-30, min_size) > 3.0) {
    }
  }

  return core::detail::clear_error_state();
}

void resolve_face_preprocess_reference_metric(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceSizingFieldState &field,
  AutoCfdSurfaceFacePreprocessState &face_preprocess
)
{
  const auto adopt_first_boundary_metric =
    [&]() {
      for(const auto &loop : face_preprocess.loops) {
        for(const auto &point : loop.points) {
          if(point.metric.usable) {
            face_preprocess.reference_metric = point.metric;
            return true;
          }
        }
      }
      return false;
    };

  if(!face_preprocess.uv_reconstruction_available) {
    if(!adopt_first_boundary_metric() &&
       !face_preprocess.boundary.loops.empty() &&
       !face_preprocess.loops.empty() &&
       !face_preprocess.loops.front().points.empty()) {
      const auto &point = face_preprocess.loops.front().points.front();
      build_metric_fallback_from_position(
        face_view.entity,
        field,
        point.position,
        point.uv[0],
        point.uv[1],
        nullptr,
        AutoCfdSurfaceMetricFallbackKind::isotropic_from_uv_identity,
        face_preprocess.reference_metric
      );
      face_preprocess.metric_fallback_used = true;
    }
    return;
  }

  const auto *outer_loop = geo::primary_outer_boundary_loop(face_preprocess.boundary);
  const AutoCfdSurfaceFaceLoopState *reference_loop = nullptr;
  if(outer_loop != nullptr && face_preprocess.boundary.primary_outer_loop_index <
                               face_preprocess.loops.size()) {
    reference_loop =
      &face_preprocess.loops[face_preprocess.boundary.primary_outer_loop_index];
  }
  else if(!face_preprocess.loops.empty()) {
    reference_loop = &face_preprocess.loops.front();
  }

  if(reference_loop == nullptr || reference_loop->points.empty()) {
    static_cast<void>(adopt_first_boundary_metric());
    return;
  }

  std::array<double, 2> centroid {0.0, 0.0};
  std::size_t point_count = 0U;
  for(const auto &point : reference_loop->points) {
    if(!point.uv_defined) {
      continue;
    }
    centroid[0] += point.uv[0];
    centroid[1] += point.uv[1];
    ++point_count;
  }

  if(point_count == 0U) {
    static_cast<void>(adopt_first_boundary_metric());
    return;
  }

  centroid[0] /= static_cast<double>(point_count);
  centroid[1] /= static_cast<double>(point_count);

  AutoCfdSurfaceFaceMetricTensor metric;
  const auto status = sample_auto_cfd_surface_face_metric(
    face_view,
    field,
    centroid[0],
    centroid[1],
    metric
  );
  if(status == base::StatusCode::ok) {
    face_preprocess.reference_metric = metric;
    face_preprocess.metric_fallback_used |=
      metric.fallback_kind != AutoCfdSurfaceMetricFallbackKind::none;
    return;
  }

  static_cast<void>(adopt_first_boundary_metric());
}

[[nodiscard]] base::StatusCode populate_single_face_preprocess_state(
  const geo::ModelView &model_view,
  const AutoCfdSurfaceSizingFieldState &field,
  const AutoCfdSurfaceBoundaryState &boundary_state,
  const AutoCfdSurfaceFaceWorkItem &face_work_item,
  AutoCfdSurfaceFacePreprocessState &face_preprocess
)
{
  const auto *face_view = model_view.find_face(face_work_item.face);
  if(face_view == nullptr) {
    return core::detail::publish_error(
      base::StatusCode::internal_error,
      "Auto CFD face preprocessing could not resolve one of the face work items from the geometry model view."
    );
  }

  face_preprocess = {};
  face_preprocess.face = face_work_item.face;
  face_preprocess.boundary = face_work_item.boundary;
  face_preprocess.loops.reserve(face_work_item.boundary.loops.size());

  auto bounds_status = geo::face_uv_bounds(*face_view, face_preprocess.uv_bounds);
  if(bounds_status == base::StatusCode::ok) {
    face_preprocess.uv_bounds_defined = true;
  }
  else if(face_work_item.boundary.has_seams &&
          !should_skip_face_curvature_status(bounds_status)) {
    return bounds_status;
  }

  if(face_work_item.boundary.has_seams) {
    
    if(supports_auto_cfd_surface_seam_unwrap(face_work_item.boundary)) {
      face_preprocess.seam_support =
        AutoCfdSurfaceSeamSupportKind::supported_subset;
    }
    else {
      face_preprocess.seam_support =
        AutoCfdSurfaceSeamSupportKind::periodic;
    }
  }

  bool uv_reconstruction_available = true;
  std::size_t recovered_uv_point_count = 0U;
  std::size_t interpolated_uv_point_count = 0U;
  for(auto loop : face_work_item.boundary.loops) {
    AutoCfdSurfaceFaceLoopState recovered_loop;
    auto status = recover_face_loop_from_boundary_dictionary(
      boundary_state,
      loop,
      recovered_loop,
      face_work_item.boundary.has_seams
    );
    if(status != base::StatusCode::ok) {
      return status;
    }

    for(auto &segment : recovered_loop.segments) {
      std::size_t segment_recovered_count = 0U;
      for(auto &point : segment.points) {
        geo::FaceUvMapping mapping;

        
        if(geo::is_valid(segment.edge_use.edge)) {
          status = geo::recover_face_uv_from_edge_use(
            *face_view,
            segment.edge_use,
            point.parameter,
            mapping
          );
        }
        else {
          status = geo::recover_face_uv(*face_view, point.position, mapping);
        }

        if(status != base::StatusCode::ok) {
          // Fallback: if pcurve lookup failed, try 3D projection
          if(geo::is_valid(segment.edge_use.edge)) {
            status = geo::recover_face_uv(*face_view, point.position, mapping);
          }
          if(status != base::StatusCode::ok) {
            if(should_skip_face_curvature_status(status)) {
              continue;
            }
            return status;
          }
        }

        point.uv = {mapping.u, mapping.v};
        point.uv_defined = true;
        ++segment_recovered_count;
      }

      recovered_uv_point_count += segment_recovered_count;
      if(segment_recovered_count != segment.points.size()) {
        const auto defined_before_interpolation = segment_recovered_count;
        if(!interpolate_missing_segment_uvs(segment.points)) {
          uv_reconstruction_available = false;
          break;
        }
        interpolated_uv_point_count +=
          segment.points.size() - defined_before_interpolation;
      }
      if(!uv_reconstruction_available) {
        break;
      }
    }

    if(uv_reconstruction_available) {
      recovered_loop.points.clear();
      for(auto &segment : recovered_loop.segments) {
        if(face_work_item.boundary.has_seams &&
           face_preprocess.seam_support ==
             AutoCfdSurfaceSeamSupportKind::supported_subset &&
           face_preprocess.uv_bounds_defined) {
          unwrap_periodic_loop_uvs(face_preprocess.uv_bounds, segment.points);
        }

        recovered_loop.points.insert(
          recovered_loop.points.end(),
          segment.points.begin(),
          segment.points.end()
        );
      }

      if(face_work_item.boundary.has_seams &&
         face_preprocess.seam_support ==
           AutoCfdSurfaceSeamSupportKind::supported_subset &&
         face_preprocess.uv_bounds_defined) {
        unwrap_periodic_loop_uvs(face_preprocess.uv_bounds, recovered_loop.points);
        face_preprocess.seam_unwrap_applied = true;
      }
    }

    face_preprocess.loops.push_back(std::move(recovered_loop));
  }

  face_preprocess.uv_reconstruction_available = uv_reconstruction_available;
  if(uv_reconstruction_available) {
    populate_recovered_face_uv_extents(face_preprocess);
  }

  if(face_preprocess.uv_reconstruction_available) {
    for(auto &loop : face_preprocess.loops) {
      // Sync UV from segment.points → loop.points. The UV recovery above
      // updates segment.points but loop.points is a stale copy from
      // recover_face_loop_from_boundary_dictionary. Without this sync,
      // the loop.points → segment.points writeback below would overwrite
      // the correct pcurve-based UV with the stale initial values.
      {
        std::size_t sync_offset = 0U;
        for(const auto &segment : loop.segments) {
          for(const auto &pt : segment.points) {
            if(sync_offset < loop.points.size()) {
              loop.points[sync_offset].uv = pt.uv;
              loop.points[sync_offset].uv_defined = pt.uv_defined;
            }
            ++sync_offset;
          }
        }
      }

      for(auto &point : loop.points) {
        if(!point.uv_defined) {
          continue;
        }

        auto status = sample_auto_cfd_surface_face_metric(
          *face_view,
          field,
          point.uv[0],
          point.uv[1],
          point.metric
        );
        if(status != base::StatusCode::ok) {
          if(!should_skip_face_curvature_status(status)) {
            return status;
          }

          build_metric_fallback_from_position(
            face_view->entity,
            field,
            point.position,
            point.uv[0],
            point.uv[1],
            nullptr,
            AutoCfdSurfaceMetricFallbackKind::isotropic_from_uv_identity,
            point.metric
          );
        }

        face_preprocess.metric_fallback_used |=
          point.metric.fallback_kind != AutoCfdSurfaceMetricFallbackKind::none;
      }

      std::size_t point_offset = 0U;
      for(auto &segment : loop.segments) {
        for(auto &point : segment.points) {
          if(point_offset >= loop.points.size()) {
            return core::detail::publish_error(
              base::StatusCode::internal_error,
              "Auto CFD face preprocessing lost segment-to-loop point alignment while rebuilding UV-local boundary data."
            );
          }
          point = loop.points[point_offset++];
        }
      }
    }
  }

  resolve_face_preprocess_reference_metric(*face_view, field, face_preprocess);

  face_preprocess.disposition =
    AutoCfdSurfaceFacePreprocessDisposition::fallback_only;
  if(face_preprocess.uv_reconstruction_available &&
     (!face_work_item.boundary.has_seams ||
      (face_preprocess.seam_support ==
         AutoCfdSurfaceSeamSupportKind::supported_subset &&
       face_preprocess.seam_unwrap_applied) ||
      face_preprocess.seam_support ==
        AutoCfdSurfaceSeamSupportKind::periodic)) {
    face_preprocess.disposition =
      face_preprocess.metric_fallback_used
        ? AutoCfdSurfaceFacePreprocessDisposition::uv_ready_with_metric_fallback
        : AutoCfdSurfaceFacePreprocessDisposition::uv_ready;
  }

  if(face_preprocess.disposition ==
       AutoCfdSurfaceFacePreprocessDisposition::fallback_only &&
     !face_preprocess.metric_fallback_used) {
    bool applied_fallback = false;
    for(auto &loop : face_preprocess.loops) {
      for(auto &point : loop.points) {
        build_metric_fallback_from_position(
          face_view->entity,
          field,
          point.position,
          point.uv[0],
          point.uv[1],
          nullptr,
          AutoCfdSurfaceMetricFallbackKind::isotropic_from_uv_identity,
          point.metric
        );
        face_preprocess.reference_metric = point.metric;
        face_preprocess.metric_fallback_used = true;
        applied_fallback = true;
        break;
      }
      if(applied_fallback) {
        std::size_t point_offset = 0U;
        for(auto &segment : loop.segments) {
          for(auto &point : segment.points) {
            if(point_offset >= loop.points.size()) {
              return core::detail::publish_error(
                base::StatusCode::internal_error,
                "Auto CFD face preprocessing lost segment-to-loop point alignment while applying a bounded face fallback metric."
              );
            }
            point = loop.points[point_offset++];
          }
        }
        break;
      }
    }
  }


  return core::detail::clear_error_state();
}

[[nodiscard]] base::StatusCode populate_face_preprocess_states(
  const geo::ModelView &model_view,
  const AutoCfdSurfaceSizingFieldState &field,
  const AutoCfdSurfaceBoundaryState &boundary_state,
  const std::vector<AutoCfdSurfaceFaceWorkItem> &face_work_items,
  std::vector<AutoCfdSurfaceFacePreprocessState> &face_preprocess_states
)
{
  face_preprocess_states.assign(model_view.faces.size(), {});

  for(const auto &face_work_item : face_work_items) {
    if(face_work_item.face.dimension != geo::TopologyDimension::face ||
       face_work_item.face.index >= face_preprocess_states.size()) {
      return core::detail::publish_error(
        base::StatusCode::internal_error,
        "Auto CFD face preprocessing encountered an out-of-range face index."
      );
    }

    auto &face_preprocess = face_preprocess_states[face_work_item.face.index];
    auto status = populate_single_face_preprocess_state(
      model_view, field, boundary_state, face_work_item, face_preprocess
    );
    if(status != base::StatusCode::ok) {
      return status;
    }
  }

  return core::detail::clear_error_state();
}

// Rebuild only the specified face's preprocess state. Other faces'
// entries in `face_preprocess_states` are left untouched. The index
// must be < face_preprocess_states.size() and must match the matching
// face_work_item's face.index.
[[nodiscard]] base::StatusCode rebuild_face_preprocess_states_for_faces(
  const geo::ModelView &model_view,
  const AutoCfdSurfaceSizingFieldState &field,
  const AutoCfdSurfaceBoundaryState &boundary_state,
  const std::vector<AutoCfdSurfaceFaceWorkItem> &face_work_items,
  const std::unordered_set<std::uint32_t> &affected_face_indices,
  std::vector<AutoCfdSurfaceFacePreprocessState> &face_preprocess_states
)
{
  for(const auto &face_work_item : face_work_items) {
    if(face_work_item.face.dimension != geo::TopologyDimension::face ||
       face_work_item.face.index >= face_preprocess_states.size()) {
      return core::detail::publish_error(
        base::StatusCode::internal_error,
        "Auto CFD face preprocessing encountered an out-of-range face index."
      );
    }
    if(affected_face_indices.find(face_work_item.face.index) ==
         affected_face_indices.end()) {
      continue;
    }

    auto &face_preprocess = face_preprocess_states[face_work_item.face.index];
    auto status = populate_single_face_preprocess_state(
      model_view, field, boundary_state, face_work_item, face_preprocess
    );
    if(status != base::StatusCode::ok) {
      return status;
    }
  }

  return core::detail::clear_error_state();
}

[[nodiscard]] const AutoCfdSurfaceCurveWorkItem *find_curve_work_item(
  const std::vector<AutoCfdSurfaceCurveWorkItem> &curve_work_items,
  geo::TopologyEntityId edge
) noexcept
{
  if(edge.dimension != geo::TopologyDimension::edge) {
    return nullptr;
  }

  const auto it = std::find_if(
    curve_work_items.begin(),
    curve_work_items.end(),
    [edge](const AutoCfdSurfaceCurveWorkItem &work_item) {
      return work_item.edge == edge;
    }
  );
  if(it == curve_work_items.end()) {
    return nullptr;
  }

  return &(*it);
}

[[nodiscard]] bool lookup_shared_endpoint_tangents(
  const geo::EdgeView &lhs,
  const AutoCfdSurfaceCurveWorkItem &lhs_curve,
  const geo::EdgeView &rhs,
  const AutoCfdSurfaceCurveWorkItem &rhs_curve,
  geo::Vector3 &lhs_tangent,
  geo::Vector3 &rhs_tangent
)
{
  lhs_tangent = {0.0, 0.0, 0.0};
  rhs_tangent = {0.0, 0.0, 0.0};

  geo::EdgeTangentSample lhs_start;
  auto status =
    geo::sample_edge_tangent(lhs, lhs_curve.curve.parameter_min, lhs_start);
  if(status != base::StatusCode::ok || !lhs_start.tangent_defined) {
    return false;
  }

  geo::EdgeTangentSample lhs_end;
  status = geo::sample_edge_tangent(lhs, lhs_curve.curve.parameter_max, lhs_end);
  if(status != base::StatusCode::ok || !lhs_end.tangent_defined) {
    return false;
  }

  geo::EdgeTangentSample rhs_start;
  status = geo::sample_edge_tangent(rhs, rhs_curve.curve.parameter_min, rhs_start);
  if(status != base::StatusCode::ok || !rhs_start.tangent_defined) {
    return false;
  }

  geo::EdgeTangentSample rhs_end;
  status = geo::sample_edge_tangent(rhs, rhs_curve.curve.parameter_max, rhs_end);
  if(status != base::StatusCode::ok || !rhs_end.tangent_defined) {
    return false;
  }

  struct EndpointPair final {
    const geo::EdgeTangentSample *lhs_sample = nullptr;
    const geo::EdgeTangentSample *rhs_sample = nullptr;
    double distance = std::numeric_limits<double>::infinity();
  };

  EndpointPair best_pair;
  const EndpointPair candidates[] {
    {&lhs_start, &rhs_start, point_distance(lhs_curve.curve.start_point, rhs_curve.curve.start_point)},
    {&lhs_start, &rhs_end, point_distance(lhs_curve.curve.start_point, rhs_curve.curve.end_point)},
    {&lhs_end, &rhs_start, point_distance(lhs_curve.curve.end_point, rhs_curve.curve.start_point)},
    {&lhs_end, &rhs_end, point_distance(lhs_curve.curve.end_point, rhs_curve.curve.end_point)},
  };

  for(const auto &candidate : candidates) {
    if(candidate.distance < best_pair.distance) {
      best_pair = candidate;
    }
  }

  if(best_pair.lhs_sample == nullptr || best_pair.rhs_sample == nullptr) {
    return false;
  }

  const double tolerance =
    curve_length_tolerance(lhs_curve.curve.approximate_length, rhs_curve.curve.approximate_length);
  if(!std::isfinite(best_pair.distance) || best_pair.distance > tolerance) {
    return false;
  }

  lhs_tangent = normalized(best_pair.lhs_sample->tangent);
  rhs_tangent = normalized(best_pair.rhs_sample->tangent);
  return vector_norm(lhs_tangent) > 0.0 && vector_norm(rhs_tangent) > 0.0;
}

[[nodiscard]] bool feature_edge_vertex_is_tangent_continuous(
  const geo::ModelView &model_view,
  const std::vector<AutoCfdSurfaceCurveWorkItem> &curve_work_items,
  geo::TopologyEntityId first_edge,
  geo::TopologyEntityId second_edge
)
{
  const auto *first_edge_view = model_view.find_edge(first_edge);
  const auto *second_edge_view = model_view.find_edge(second_edge);
  const auto *first_curve = find_curve_work_item(curve_work_items, first_edge);
  const auto *second_curve = find_curve_work_item(curve_work_items, second_edge);
  if(first_edge_view == nullptr || second_edge_view == nullptr ||
     first_curve == nullptr || second_curve == nullptr) {
    return false;
  }

  geo::Vector3 first_tangent;
  geo::Vector3 second_tangent;
  if(!lookup_shared_endpoint_tangents(
       *first_edge_view,
       *first_curve,
       *second_edge_view,
       *second_curve,
       first_tangent,
       second_tangent
     )) {
    return false;
  }

  const double alignment = std::clamp(
    std::abs(dot_product(first_tangent, second_tangent)),
    0.0,
    1.0
  );
  const double continuity_threshold =
    std::cos(kTangentContinuityToleranceDegrees * (kPi / 180.0));
  return alignment >= continuity_threshold;
}

void initialize_sizing_field(
  AutoCfdSurfaceSizingFieldState &field,
  double minimum_length,
  double maximum_length,
  double distortion_angle,
  double growth_rate,
  bool proximity_enabled,
  bool self_proximity_enabled,
  bool pid_proximity_enabled,
  double maximum_normals_angle_degrees,
  double length_to_gap_ratio,
  double proximity_minimum_length
)
{
  field = {};
  field.minimum_length = minimum_length;
  field.maximum_length = maximum_length;
  field.distortion_angle = distortion_angle;
  field.growth_rate = growth_rate;
  field.proximity_enabled = proximity_enabled;
  field.self_proximity_enabled = self_proximity_enabled;
  field.pid_proximity_enabled = pid_proximity_enabled;
  field.proximity_layer_count = kProximityLayerCount;
  field.proximity_normal_angle_threshold_degrees =
    180.0 - maximum_normals_angle_degrees;
  field.proximity_length_to_gap_ratio = length_to_gap_ratio;
  field.proximity_minimum_length =
    std::max(0.0, proximity_minimum_length);
}

void populate_owner_target_sizes(
  const geo::ModelView &model_view,
  const ResolvedMeshSizeControls &size_controls,
  AutoCfdSurfaceSizingFieldState &field
)
{
  field.face_target_sizes.assign(model_view.faces.size(), field.maximum_length);
  for(const auto &face_view : model_view.faces) {
    field.face_target_sizes[face_view.entity.index] = clamp_auto_cfd_surface_size(
      field,
      effective_face_target_size(size_controls, face_view.entity)
    );
  }

  field.edge_target_sizes.assign(model_view.edges.size(), field.maximum_length);
  for(const auto &edge_view : model_view.edges) {
    field.edge_target_sizes[edge_view.entity.index] = clamp_auto_cfd_surface_size(
      field,
      effective_edge_target_size(size_controls, edge_view)
    );
  }
}

void populate_face_proximity_buckets(
  const geo::ModelView &model_view,
  AutoCfdSurfaceSizingFieldState &field
)
{
  field.face_proximity_shell_ids.assign(
    model_view.faces.size(),
    invalid_auto_cfd_surface_part_index
  );
  field.face_proximity_part_ids.assign(
    model_view.faces.size(),
    invalid_auto_cfd_surface_part_index
  );
  field.proximity_stats.shell_count = 0U;
  field.proximity_stats.part_count = 0U;

  std::vector<geo::TopologyEntityId> stack;
  std::uint32_t next_shell_index = 0U;

  for(const auto &face_view : model_view.faces) {
    if(face_view.entity.index >= field.face_proximity_shell_ids.size() ||
       field.face_proximity_shell_ids[face_view.entity.index] !=
         invalid_auto_cfd_surface_part_index) {
      continue;
    }

    stack.clear();
    stack.push_back(face_view.entity);
    field.face_proximity_shell_ids[face_view.entity.index] = next_shell_index;

    while(!stack.empty()) {
      const auto current_face = stack.back();
      stack.pop_back();

      const auto *current_face_view = model_view.find_face(current_face);
      if(current_face_view == nullptr) {
        continue;
      }

      for(const auto edge : current_face_view->edge_ids) {
        const auto *edge_view = model_view.find_edge(edge);
        if(edge_view == nullptr) {
          continue;
        }

        for(const auto adjacent_face : edge_view->face_ids) {
          if(adjacent_face.dimension != geo::TopologyDimension::face ||
             adjacent_face.index >= field.face_proximity_shell_ids.size()) {
            continue;
          }
          if(field.face_proximity_shell_ids[adjacent_face.index] !=
             invalid_auto_cfd_surface_part_index) {
            continue;
          }

          field.face_proximity_shell_ids[adjacent_face.index] = next_shell_index;
          stack.push_back(adjacent_face);
        }
      }
    }

    ++next_shell_index;
  }

  field.proximity_stats.shell_count = next_shell_index;
  if(next_shell_index == 0U) {
    return;
  }

  ProximityBucketUnionFind part_union(next_shell_index);
  for(const auto &region_view : model_view.regions) {
    std::uint32_t region_shell_index = invalid_auto_cfd_surface_part_index;
    for(const auto region_face : region_view.face_ids) {
      if(region_face.dimension != geo::TopologyDimension::face ||
         region_face.index >= field.face_proximity_shell_ids.size()) {
        continue;
      }

      const auto shell_index = field.face_proximity_shell_ids[region_face.index];
      if(shell_index == invalid_auto_cfd_surface_part_index) {
        continue;
      }

      if(region_shell_index == invalid_auto_cfd_surface_part_index) {
        region_shell_index = shell_index;
        continue;
      }
      part_union.unite(region_shell_index, shell_index);
    }
  }

  std::unordered_map<std::uint32_t, std::uint32_t> compact_part_indices;
  compact_part_indices.reserve(next_shell_index);
  std::uint32_t next_part_index = 0U;

  for(const auto &face_view : model_view.faces) {
    if(face_view.entity.index >= field.face_proximity_shell_ids.size() ||
       face_view.entity.index >= field.face_proximity_part_ids.size()) {
      continue;
    }

    const auto shell_index = field.face_proximity_shell_ids[face_view.entity.index];
    if(shell_index == invalid_auto_cfd_surface_part_index) {
      continue;
    }

    const auto root_part_index = part_union.find(shell_index);
    auto [part_it, inserted] =
      compact_part_indices.emplace(root_part_index, next_part_index);
    if(inserted) {
      ++next_part_index;
    }
    field.face_proximity_part_ids[face_view.entity.index] = part_it->second;
  }

  field.proximity_stats.part_count = next_part_index;
}

[[nodiscard]] base::StatusCode collect_proxy_triangle_query_data(
  const mesh::Domain &proxy_domain,
  const AutoCfdSurfaceSizingFieldState &field,
  std::vector<AutoCfdSurfaceProxyTriangle> &triangles
)
{
  triangles.clear();

  for(const auto &entity_group : proxy_domain.entity_groups()) {
    if(entity_group.role() != mesh::EntityGroupRole::geometric_proxy ||
       entity_group.order() != mesh::EntityOrder::face) {
      continue;
    }

    for(std::uint32_t face_index = 0U; face_index < entity_group.faces().size(); ++face_index) {
      const mesh::EntityRef face_ref {entity_group.id(), face_index};
      const auto owner = proxy_domain.face_topology_owner(face_ref);
      if(owner.dimension != geo::TopologyDimension::face ||
         owner.index >= field.face_proximity_shell_ids.size() ||
         owner.index >= field.face_proximity_part_ids.size()) {
        continue;
      }

      const auto nodes = proxy_domain.face_nodes(face_ref);
      if(nodes.size != 3U) {
        continue;
      }

      AutoCfdSurfaceProxyTriangle triangle;
      triangle.owner = owner;
      triangle.shell_index = field.face_proximity_shell_ids[owner.index];
      triangle.part_index = field.face_proximity_part_ids[owner.index];
      triangle.vertices[0] = proxy_domain.node(nodes[0]).coordinates;
      triangle.vertices[1] = proxy_domain.node(nodes[1]).coordinates;
      triangle.vertices[2] = proxy_domain.node(nodes[2]).coordinates;

      triangle.bounds_min = triangle.vertices[0];
      triangle.bounds_max = triangle.vertices[0];
      for(std::size_t vertex_index = 1U; vertex_index < triangle.vertices.size();
          ++vertex_index) {
        for(std::size_t axis = 0U; axis < 3U; ++axis) {
          triangle.bounds_min[axis] =
            std::min(triangle.bounds_min[axis], triangle.vertices[vertex_index][axis]);
          triangle.bounds_max[axis] =
            std::max(triangle.bounds_max[axis], triangle.vertices[vertex_index][axis]);
        }
      }
      triangle.centroid = {
        (triangle.vertices[0][0] + triangle.vertices[1][0] + triangle.vertices[2][0]) / 3.0,
        (triangle.vertices[0][1] + triangle.vertices[1][1] + triangle.vertices[2][1]) / 3.0,
        (triangle.vertices[0][2] + triangle.vertices[1][2] + triangle.vertices[2][2]) / 3.0,
      };

      const auto face_normal = normalized(cross_product(
        subtract_points(triangle.vertices[1], triangle.vertices[0]),
        subtract_points(triangle.vertices[2], triangle.vertices[0])
      ));
      if(vector_norm(face_normal) <= 0.0) {
        continue;
      }
      triangle.normal = face_normal;

      triangles.push_back(std::move(triangle));
    }
  }

  return core::detail::clear_error_state();
}

[[nodiscard]] base::StatusCode populate_curve_work_items(
  const geo::ModelView &model_view,
  const AutoCfdSurfaceSizingFieldState &field,
  std::vector<AutoCfdSurfaceCurveWorkItem> &curve_work_items
)
{
  curve_work_items.clear();
  curve_work_items.reserve(model_view.edges.size());

  for(const auto &edge_view : model_view.edges) {
    AutoCfdSurfaceCurveWorkItem work_item;
    work_item.edge = edge_view.entity;
    work_item.face_ids = edge_view.face_ids;
    work_item.vertex_ids = edge_view.vertex_ids;
    work_item.target_size = effective_auto_cfd_surface_owner_target_size(
      field,
      edge_view.entity
    );

    const auto status = geo::edge_curve_info(edge_view, work_item.curve);
    if(status != base::StatusCode::ok) {
      return status;
    }

    curve_work_items.push_back(std::move(work_item));
  }

  return core::detail::clear_error_state();
}

void populate_face_work_items(
  const geo::ModelView &model_view,
  const AutoCfdSurfaceSizingFieldState &field,
  std::vector<AutoCfdSurfaceFaceWorkItem> &face_work_items
)
{
  face_work_items.clear();
  face_work_items.reserve(model_view.faces.size());

 
  for(const auto &face_view : model_view.faces) {
  
    AutoCfdSurfaceFaceWorkItem work_item;
    work_item.face = face_view.entity;
    work_item.boundary = geo::ordered_boundary_loops(face_view);
    work_item.target_size = effective_auto_cfd_surface_owner_target_size(
      field,
      face_view.entity
    );
    if(face_view.entity.index < field.face_proximity_part_ids.size()) {
      work_item.proximity_part_index = field.face_proximity_part_ids[face_view.entity.index];
    }
    face_work_items.push_back(std::move(work_item));
  }
}

void populate_background_support_bounds(
  const geo::ModelView &model_view,
  const std::vector<AutoCfdSurfaceCurveWorkItem> &curve_work_items,
  AutoCfdSurfaceSizingFieldState &field
)
{
  MutableBackgroundBounds bounds;

  for(const auto &curve_work_item : curve_work_items) {
    extend_bounds(bounds, curve_work_item.curve.start_point);
    extend_bounds(bounds, curve_work_item.curve.end_point);
  }

  for(const auto &face_view : model_view.faces) {
    geo::FaceUvBounds uv_bounds;
    if(geo::face_uv_bounds(face_view, uv_bounds) != base::StatusCode::ok) {
      continue;
    }

    const std::array<std::array<double, 2>, 5> uv_samples = {{
      {uv_bounds.u_min, uv_bounds.v_min},
      {uv_bounds.u_min, uv_bounds.v_max},
      {uv_bounds.u_max, uv_bounds.v_min},
      {uv_bounds.u_max, uv_bounds.v_max},
      {
        0.5 * (uv_bounds.u_min + uv_bounds.u_max),
        0.5 * (uv_bounds.v_min + uv_bounds.v_max),
      },
    }};

    for(const auto &uv : uv_samples) {
      geo::FaceSample face_sample;
      if(geo::sample_face(face_view, uv[0], uv[1], face_sample) != base::StatusCode::ok) {
        continue;
      }
      extend_bounds(bounds, face_sample.position);
    }
  }

  if(bounds.defined) {
    field.background_bounds.minimum = bounds.minimum;
    field.background_bounds.maximum = bounds.maximum;
    field.background_bounds.defined = true;
  }
}

void append_curvature_source(
  AutoCfdSurfaceSizingFieldState &field,
  geo::TopologyEntityId owner,
  AutoCfdSurfaceSizingSourceKind kind,
  const geo::Point3 &position,
  double primary_parameter,
  double secondary_parameter,
  double raw_target_size
);

void append_edge_curvature_sources(
  const geo::ModelView &model_view,
  const AutoCfdSurfaceCurveWorkItem &curve_work_item,
  const AutoCfdSurfaceSizingFieldState &field,
  AutoCfdSurfaceSizingFieldState &mutable_field
);

[[nodiscard]] base::StatusCode populate_anchor_candidates(
  const geo::ModelView &model_view,
  const std::vector<AutoCfdSurfaceCurveWorkItem> &curve_work_items,
  double feature_angle_degrees,
  std::vector<AutoCfdSurfaceAnchorCandidate> &anchor_candidates
)
{
  geo::FeatureEdgeReport feature_report;
  geo::FeatureEdgeOptions feature_options;
  feature_options.feature_angle_degrees = feature_angle_degrees;
  feature_options.include_boundary_edges = true;
  feature_options.include_non_manifold_edges = true;

  auto status = geo::feature_edges(
    model_view.model_handle,
    feature_report,
    feature_options,
    model_view.context_handle
  );
  if(status != base::StatusCode::ok) {
    return status;
  }

  std::unordered_set<std::uint64_t> feature_edge_keys;
  feature_edge_keys.reserve(feature_report.edges.size());
  for(const auto edge : feature_report.edges) {
    feature_edge_keys.insert(pack_entity_id(edge));
  }

  anchor_candidates.clear();
  anchor_candidates.reserve(model_view.vertices.size());

  for(const auto &vertex_view : model_view.vertices) {
    AutoCfdSurfaceAnchorCandidate anchor;
    anchor.vertex = vertex_view.entity;
    anchor.incident_edges = vertex_view.edge_ids;

    for(const auto edge : vertex_view.edge_ids) {
      if(feature_edge_keys.find(pack_entity_id(edge)) != feature_edge_keys.end()) {
        anchor.incident_feature_edges.push_back(edge);
      }
    }

    if(anchor.incident_feature_edges.empty()) {
      anchor_candidates.push_back(std::move(anchor));
      continue;
    }

    if(anchor.incident_feature_edges.size() != 2U) {
      anchor.frozen = true;
      anchor.source = AutoCfdSurfaceAnchorSource::hard_point;
      anchor_candidates.push_back(std::move(anchor));
      continue;
    }

    anchor.tangent_continuous = feature_edge_vertex_is_tangent_continuous(
      model_view,
      curve_work_items,
      anchor.incident_feature_edges[0],
      anchor.incident_feature_edges[1]
    );
    if(!anchor.tangent_continuous) {
      anchor.frozen = true;
      anchor.source =
        AutoCfdSurfaceAnchorSource::non_tangent_feature_intersection;
    }

    anchor_candidates.push_back(std::move(anchor));
  }

  return core::detail::clear_error_state();
}

[[nodiscard]] const geo::FaceBoundaryEdgeUse *find_face_edge_use_for_vertex(
  const std::vector<AutoCfdSurfaceFaceWorkItem> &face_work_items,
  geo::TopologyEntityId edge,
  geo::TopologyEntityId vertex
) noexcept
{
  for(const auto &face_work_item : face_work_items) {
    for(const auto &loop : face_work_item.boundary.loops) {
      for(const auto &edge_use : loop.edge_uses) {
        if(edge_use.edge != edge) {
          continue;
        }
        if(edge_use.start_vertex == vertex || edge_use.end_vertex == vertex) {
          return &edge_use;
        }
      }
    }
  }

  return nullptr;
}

[[nodiscard]] bool lookup_vertex_oriented_edge_tangent(
  const geo::ModelView &model_view,
  const std::vector<AutoCfdSurfaceFaceWorkItem> &face_work_items,
  const AutoCfdSurfaceCurveWorkItem &curve_work_item,
  geo::TopologyEntityId vertex,
  geo::Point3 &position,
  geo::Vector3 &tangent
)
{
  position = {0.0, 0.0, 0.0};
  tangent = {0.0, 0.0, 0.0};

  const auto *edge_view = model_view.find_edge(curve_work_item.edge);
  const auto *edge_use =
    find_face_edge_use_for_vertex(face_work_items, curve_work_item.edge, vertex);
  if(edge_view == nullptr || edge_use == nullptr) {
    return false;
  }

  bool vertex_at_parameter_min = false;
  if(edge_use->start_vertex == vertex) {
    vertex_at_parameter_min = edge_use->same_orientation_as_edge;
  }
  else if(edge_use->end_vertex == vertex) {
    vertex_at_parameter_min = !edge_use->same_orientation_as_edge;
  }
  else {
    return false;
  }

  geo::EdgeTangentSample tangent_sample;
  const double parameter = vertex_at_parameter_min
                             ? curve_work_item.curve.parameter_min
                             : curve_work_item.curve.parameter_max;
  if(geo::sample_edge_tangent(*edge_view, parameter, tangent_sample) !=
       base::StatusCode::ok ||
     !tangent_sample.tangent_defined) {
    return false;
  }

  tangent = normalized(tangent_sample.tangent);
  if(vector_norm(tangent) <= 0.0) {
    return false;
  }
  if(!vertex_at_parameter_min) {
    tangent = {-tangent[0], -tangent[1], -tangent[2]};
  }
  position = tangent_sample.position;
  return true;
}

void populate_edge_curvature_sources(
  const geo::ModelView &model_view,
  const std::vector<AutoCfdSurfaceCurveWorkItem> &curve_work_items,
  AutoCfdSurfaceSizingFieldState &field
)
{
  field.curvature_sources.clear();
  for(const auto &curve_work_item : curve_work_items) {
    append_edge_curvature_sources(model_view, curve_work_item, field, field);
  }
}

void populate_sharp_corner_sources(
  const geo::ModelView &model_view,
  const std::vector<AutoCfdSurfaceFaceWorkItem> &face_work_items,
  const std::vector<AutoCfdSurfaceCurveWorkItem> &curve_work_items,
  const std::vector<AutoCfdSurfaceAnchorCandidate> &anchor_candidates,
  double sharp_angle_limit_degrees,
  double sharp_angle_length,
  AutoCfdSurfaceSizingFieldState &field
)
{
  if(sharp_angle_limit_degrees <= 0.0 || sharp_angle_limit_degrees >= 180.0 ||
     sharp_angle_length <= 0.0) {
    return;
  }

  for(const auto &anchor : anchor_candidates) {
    if(anchor.vertex.dimension != geo::TopologyDimension::vertex ||
       anchor.incident_feature_edges.empty()) {
      continue;
    }

    if(anchor.incident_feature_edges.size() == 2U) {
      const auto *first_curve =
        find_curve_work_item(curve_work_items, anchor.incident_feature_edges[0]);
      const auto *second_curve =
        find_curve_work_item(curve_work_items, anchor.incident_feature_edges[1]);
      if(first_curve == nullptr || second_curve == nullptr) {
        continue;
      }

      geo::Point3 first_position;
      geo::Point3 second_position;
      geo::Vector3 first_tangent;
      geo::Vector3 second_tangent;
      if(!lookup_vertex_oriented_edge_tangent(
           model_view,
           face_work_items,
           *first_curve,
           anchor.vertex,
           first_position,
           first_tangent
         ) ||
         !lookup_vertex_oriented_edge_tangent(
           model_view,
           face_work_items,
           *second_curve,
           anchor.vertex,
           second_position,
           second_tangent
         )) {
        continue;
      }

      const double angle_radians = std::acos(std::clamp(
        dot_product(first_tangent, second_tangent),
        -1.0,
        1.0
      ));
      const double angle_degrees = angle_radians * (180.0 / kPi);
      if(!std::isfinite(angle_degrees) ||
         angle_degrees > sharp_angle_limit_degrees) {
        continue;
      }

      const geo::Point3 source_position {
        0.5 * (first_position[0] + second_position[0]),
        0.5 * (first_position[1] + second_position[1]),
        0.5 * (first_position[2] + second_position[2]),
      };
      append_curvature_source(
        field,
        anchor.vertex,
        AutoCfdSurfaceSizingSourceKind::edge_curvature,
        source_position,
        angle_degrees,
        0.0,
        sharp_angle_length
      );
      continue;
    }

    if(anchor.incident_feature_edges.size() < 3U) {
      continue;
    }

    const auto *curve_work_item =
      find_curve_work_item(curve_work_items, anchor.incident_feature_edges.front());
    if(curve_work_item == nullptr) {
      continue;
    }

    geo::Point3 source_position;
    geo::Vector3 source_tangent;
    if(!lookup_vertex_oriented_edge_tangent(
         model_view,
         face_work_items,
         *curve_work_item,
         anchor.vertex,
         source_position,
         source_tangent
       )) {
      continue;
    }

    append_curvature_source(
      field,
      anchor.vertex,
      AutoCfdSurfaceSizingSourceKind::edge_curvature,
      source_position,
      sharp_angle_limit_degrees,
      0.0,
      sharp_angle_length
    );
  }
}

void append_curvature_source(
  AutoCfdSurfaceSizingFieldState &field,
  geo::TopologyEntityId owner,
  AutoCfdSurfaceSizingSourceKind kind,
  const geo::Point3 &position,
  double primary_parameter,
  double secondary_parameter,
  double raw_target_size
)
{
  if(!std::isfinite(raw_target_size) || raw_target_size <= 0.0) {
    return;
  }

  AutoCfdSurfaceSizingSource source;
  source.owner = owner;
  source.kind = kind;
  source.position = position;
  source.primary_parameter = primary_parameter;
  source.secondary_parameter = secondary_parameter;
  source.raw_target_size = raw_target_size;
  source.target_size = clamp_auto_cfd_surface_size(field, raw_target_size);
  source.relaxed_target_size = source.target_size;
  field.curvature_sources.push_back(std::move(source));
}

void append_proximity_source(
  AutoCfdSurfaceSizingFieldState &field,
  geo::TopologyEntityId owner,
  AutoCfdSurfaceSizingSourceKind kind,
  geo::TopologyEntityId hit_owner,
  const geo::Point3 &position,
  double primary_parameter,
  double secondary_parameter,
  double hit_distance,
  double hit_normal_angle_degrees,
  double curvature_size
)
{
  if(!std::isfinite(hit_distance) || hit_distance <= 0.0) {
    return;
  }

  AutoCfdSurfaceSizingSource source;
  source.owner = owner;
  source.kind = kind;
  source.hit_owner = hit_owner;
  source.position = position;
  source.primary_parameter = primary_parameter;
  source.secondary_parameter = secondary_parameter;
  source.hit_distance = hit_distance;
  source.hit_normal_angle_degrees = hit_normal_angle_degrees;
  source.raw_target_size = std::max(
    hit_distance * field.proximity_length_to_gap_ratio,
    field.proximity_minimum_length
  );
  source.target_size = clamp_auto_cfd_surface_size(field, source.raw_target_size);
  // Relaxation: proximity should not reduce the local size below a fraction
  // of the curvature-driven field. We cap the effective proximity size so
  // it does not dominate the curvature field in sensitive regions.
  constexpr double kProximityRelaxationFraction = 0.2;
  if(std::isfinite(curvature_size) && curvature_size > 0.0) {
    const double relaxation_floor = clamp_auto_cfd_surface_size(
      field, kProximityRelaxationFraction * curvature_size);
    source.relaxed_target_size = std::max(source.target_size, relaxation_floor);
  } else {
    source.relaxed_target_size = source.target_size;
  }
  field.proximity_sources.push_back(std::move(source));
}

void consider_proximity_query_hit(
  const AutoCfdSurfaceProximityRayHit &hit,
  AutoCfdSurfaceSizingSourceKind kind,
  ProximityHitCandidate &candidate
)
{
  if(!hit.found ||
     !std::isfinite(hit.hit_distance) ||
     hit.hit_distance <= 0.0) {
    return;
  }

  if(candidate.found && hit.hit_distance >= candidate.hit_distance) {
    return;
  }

  candidate.found = true;
  candidate.kind = kind;
  candidate.hit_owner = hit.hit_owner;
  candidate.hit_distance = hit.hit_distance;
  candidate.hit_normal_angle_degrees = hit.hit_normal_angle_degrees;
}

void query_proximity_candidates(
  const AutoCfdSurfaceProximityIndex &proximity_index,
  const geo::Point3 &origin,
  const geo::Vector3 &source_normal,
  geo::TopologyEntityId excluded_owner,
  std::uint32_t source_part_index,
  double max_hit_distance,
  double minimum_normal_angle_degrees,
  bool allow_self_proximity,
  bool allow_pid_proximity,
  ProximityHitCandidate &self_hit,
  ProximityHitCandidate &pid_hit
)
{
  if(proximity_index.empty() ||
     !std::isfinite(max_hit_distance) ||
     max_hit_distance <= 0.0) {
    return;
  }

  const std::array<geo::Vector3, 2> ray_directions = {
    source_normal,
    scaled_vector(source_normal, -1.0),
  };

  for(const auto &direction : ray_directions) {
    AutoCfdSurfaceProximityRayRequest request;
    request.origin = origin;
    request.direction = direction;
    request.source_normal = source_normal;
    request.excluded_owner = excluded_owner;
    request.source_part_index = source_part_index;
    request.max_distance = max_hit_distance;
    request.minimum_normal_angle_degrees = minimum_normal_angle_degrees;

    if(allow_self_proximity) {
      request.require_same_part = true;
      request.require_different_part = false;

      AutoCfdSurfaceProximityRayHit self_query_hit;
      if(proximity_index.query_nearest_hit(request, self_query_hit)) {
        consider_proximity_query_hit(
          self_query_hit,
          AutoCfdSurfaceSizingSourceKind::self_proximity,
          self_hit
        );
      }
    }

    if(allow_pid_proximity) {
      request.require_same_part = false;
      request.require_different_part = true;

      AutoCfdSurfaceProximityRayHit pid_query_hit;
      if(proximity_index.query_nearest_hit(request, pid_query_hit)) {
        consider_proximity_query_hit(
          pid_query_hit,
          AutoCfdSurfaceSizingSourceKind::pid_proximity,
          pid_hit
        );
      }
    }
  }
}

void append_proximity_candidate_source(
  AutoCfdSurfaceSizingFieldState &field,
  geo::TopologyEntityId owner,
  const geo::Point3 &position,
  double primary_parameter,
  double secondary_parameter,
  double source_size,
  const ProximityHitCandidate &candidate
)
{
  if(!candidate.found) {
    return;
  }
  // min_gap filter: gaps smaller than min_size * layer_count are too narrow
  // to mesh and should be ignored.
  const double proximity_floor = std::max(
    field.minimum_length,
    field.proximity_minimum_length
  );
  const double min_gap = proximity_floor * static_cast<double>(kProximityLayerCount);
  if(candidate.hit_distance < min_gap) {
    return;
  }
  // Significance filter: proximity target must be below source_size.
  const double candidate_target_size = std::max(
    candidate.hit_distance * field.proximity_length_to_gap_ratio,
    field.proximity_minimum_length
  );
  if(candidate_target_size >= source_size) {
    return;
  }
  append_proximity_source(
    field,
    owner,
    candidate.kind,
    candidate.hit_owner,
    position,
    primary_parameter,
    secondary_parameter,
    candidate.hit_distance,
    candidate.hit_normal_angle_degrees,
    source_size
  );
}

[[nodiscard]] bool sample_support_face_normal(
  const geo::FaceView &face_view,
  const geo::Point3 &position,
  geo::Vector3 &normal
)
{
  geo::FaceProjection projection;
  if(geo::project_point_to_face(face_view, position, projection) != base::StatusCode::ok) {
    return false;
  }

  if(projection.normal_defined) {
    normal = projection.normal;
    return true;
  }

  geo::FaceSample face_sample;
  if(geo::sample_face(face_view, projection.u, projection.v, face_sample) !=
       base::StatusCode::ok ||
     !face_sample.normal_defined) {
    return false;
  }

  normal = face_sample.normal;
  return true;
}

void append_face_proximity_sources(
  const geo::FaceSample &face_sample,
  const AutoCfdSurfaceFaceWorkItem &face_work_item,
  const AutoCfdSurfaceSizingFieldState &base_field,
  const AutoCfdSurfaceProximityIndex &proximity_index,
  AutoCfdSurfaceSizingFieldState &mutable_field
)
{
  if(!face_sample.normal_defined || proximity_index.empty() ||
     mutable_field.proximity_layer_count == 0U) {
    return;
  }
  ++mutable_field.proximity_stats.face_sample_count;

  const auto source_normal = normalized(face_sample.normal);
  if(vector_norm(source_normal) <= 0.0) {
    return;
  }

  const double source_size = query_auto_cfd_surface_sizing_field(
    base_field,
    face_work_item.face,
    face_sample.position
  );
  const double proximity_floor = std::max(
    mutable_field.minimum_length,
    mutable_field.proximity_minimum_length
  );
  if(!std::isfinite(source_size) || source_size <= proximity_floor) {
    return;
  }

  const double max_hit_distance =
    source_size / mutable_field.proximity_length_to_gap_ratio;
  if(!std::isfinite(max_hit_distance) || max_hit_distance <= 0.0) {
    return;
  }

  ProximityHitCandidate self_hit;
  ProximityHitCandidate pid_hit;

  query_proximity_candidates(
    proximity_index,
    face_sample.position,
    source_normal,
    face_work_item.face,
    face_work_item.proximity_part_index,
    max_hit_distance,
    mutable_field.proximity_normal_angle_threshold_degrees,
    mutable_field.self_proximity_enabled,
    mutable_field.pid_proximity_enabled,
    self_hit,
    pid_hit
  );

  append_proximity_candidate_source(
    mutable_field,
    face_work_item.face,
    face_sample.position,
    face_sample.u,
    face_sample.v,
    source_size,
    self_hit
  );
  append_proximity_candidate_source(
    mutable_field,
    face_work_item.face,
    face_sample.position,
    face_sample.u,
    face_sample.v,
    source_size,
    pid_hit
  );
}

void append_edge_proximity_sources(
  const geo::ModelView &model_view,
  const AutoCfdSurfaceCurveWorkItem &curve_work_item,
  const AutoCfdSurfaceSizingFieldState &base_field,
  const AutoCfdSurfaceProximityIndex &proximity_index,
  AutoCfdSurfaceSizingFieldState &mutable_field
)
{
  if(proximity_index.empty() ||
     mutable_field.proximity_layer_count == 0U ||
     curve_work_item.face_ids.empty()) {
    return;
  }

  const auto *edge_view = model_view.find_edge(curve_work_item.edge);
  if(edge_view == nullptr) {
    return;
  }
  if(curve_work_item.curve.approximate_length <= kCurvatureTolerance ||
     curve_work_item.curve.parameter_max <= curve_work_item.curve.parameter_min) {
    return;
  }

  const double target_segment_length = std::max(
    base_field.minimum_length,
    std::min(
      curve_work_item.target_size,
      std::max(curve_work_item.curve.approximate_length / 8.0, base_field.minimum_length)
    )
  );

  geo::EdgeCurveSamples sampled_curve;
  const auto status = geo::sample_edge_curve(
    *edge_view,
    {
      target_segment_length,
      4U,
    },
    sampled_curve
  );
  if(status != base::StatusCode::ok || sampled_curve.samples.empty()) {
    return;
  }

  for(const auto &edge_sample : sampled_curve.samples) {
    const double source_size = query_auto_cfd_surface_sizing_field(
      base_field,
      curve_work_item.edge,
      edge_sample.position
    );
    const double proximity_floor = std::max(
      mutable_field.minimum_length,
      mutable_field.proximity_minimum_length
    );
    if(!std::isfinite(source_size) || source_size <= proximity_floor) {
      continue;
    }

    const double max_hit_distance =
      source_size / mutable_field.proximity_length_to_gap_ratio;
    if(!std::isfinite(max_hit_distance) || max_hit_distance <= 0.0) {
      continue;
    }

    ++mutable_field.proximity_stats.edge_sample_count;
    ProximityHitCandidate self_hit;
    ProximityHitCandidate pid_hit;

    for(const auto support_face_id : curve_work_item.face_ids) {
      if(support_face_id.dimension != geo::TopologyDimension::face ||
         support_face_id.index >= base_field.face_proximity_part_ids.size()) {
        continue;
      }

      const auto *support_face_view = model_view.find_face(support_face_id);
      if(support_face_view == nullptr) {
        continue;
      }

      geo::Vector3 support_normal;
      if(!sample_support_face_normal(
           *support_face_view,
           edge_sample.position,
           support_normal
         )) {
        continue;
      }

      const auto source_normal = normalized(support_normal);
      if(vector_norm(source_normal) <= 0.0) {
        continue;
      }

      query_proximity_candidates(
        proximity_index,
        edge_sample.position,
        source_normal,
        support_face_id,
        base_field.face_proximity_part_ids[support_face_id.index],
        max_hit_distance,
        mutable_field.proximity_normal_angle_threshold_degrees,
        mutable_field.self_proximity_enabled,
        mutable_field.pid_proximity_enabled,
        self_hit,
        pid_hit
      );
    }

    append_proximity_candidate_source(
      mutable_field,
      curve_work_item.edge,
      edge_sample.position,
      edge_sample.parameter,
      0.0,
      source_size,
      self_hit
    );
    append_proximity_candidate_source(
      mutable_field,
      curve_work_item.edge,
      edge_sample.position,
      edge_sample.parameter,
      0.0,
      source_size,
      pid_hit
    );
  }
}

void append_edge_curvature_sources(
  const geo::ModelView &model_view,
  const AutoCfdSurfaceCurveWorkItem &curve_work_item,
  const AutoCfdSurfaceSizingFieldState &field,
  AutoCfdSurfaceSizingFieldState &mutable_field
)
{
  const auto *edge_view = model_view.find_edge(curve_work_item.edge);
  if(edge_view == nullptr) {
    return;
  }

  if(curve_work_item.curve.approximate_length <= kCurvatureTolerance ||
     curve_work_item.curve.parameter_max <= curve_work_item.curve.parameter_min) {
    return;
  }

  const double target_segment_length = std::max(
    field.minimum_length,
    std::min(
      curve_work_item.target_size,
      std::max(curve_work_item.curve.approximate_length / 8.0, field.minimum_length)
    )
  );

  geo::EdgeCurveSamples sampled_curve;
  const auto status = geo::sample_edge_curve(
    *edge_view,
    {
      target_segment_length,
      kMinimumEdgeCurvatureSegmentCount,
    },
    sampled_curve
  );
  if(status != base::StatusCode::ok || sampled_curve.samples.size() < 3U) {
    return;
  }

  // Per-sample target size. Interior samples get curvature from the
  // tangent angle change over the prev+next span; endpoints inherit
  // from their adjacent interior sample so consecutive-pair edge
  // sources span the full curve.
  const auto sample_count = sampled_curve.samples.size();
  std::vector<double> sample_size(sample_count, field.maximum_length);
  for(std::size_t sample_index = 1U;
      sample_index + 1U < sample_count;
      ++sample_index) {
    const auto &previous = sampled_curve.samples[sample_index - 1U];
    const auto &current = sampled_curve.samples[sample_index];
    const auto &next = sampled_curve.samples[sample_index + 1U];
    if(!previous.tangent_defined || !current.tangent_defined || !next.tangent_defined) {
      continue;
    }

    const double previous_distance = point_distance(previous.position, current.position);
    const double next_distance = point_distance(current.position, next.position);
    const double span_length = previous_distance + next_distance;
    if(!std::isfinite(span_length) || span_length <= kCurvatureTolerance) {
      continue;
    }

    const auto previous_tangent = normalized(previous.tangent);
    const auto next_tangent = normalized(next.tangent);
    if(vector_norm(previous_tangent) <= 0.0 || vector_norm(next_tangent) <= 0.0) {
      continue;
    }

    const double angle = std::acos(std::clamp(
      dot_product(previous_tangent, next_tangent),
      -1.0,
      1.0
    ));
    const double curvature = angle / span_length;
    const double raw_target = curvature_to_target_size(curvature, field.distortion_angle);
    if(!std::isfinite(raw_target) || raw_target <= 0.0) {
      continue;  // near-zero curvature, leave sample_size[i] at max
    }
    sample_size[sample_index] = std::clamp(
      raw_target, field.minimum_length, field.maximum_length);
  }

  // Extend interior sample sizes to the curve endpoints so edge sources
  // cover the full curve length.
  if(sample_count >= 3U) {
    sample_size[0] = sample_size[1];
    sample_size[sample_count - 1U] = sample_size[sample_count - 2U];
  }

  const double flat_threshold = 0.9 * field.maximum_length;
  for(std::size_t i = 0U; i + 1U < sample_count; ++i) {
    const double s0 = sample_size[i];
    const double s1 = sample_size[i + 1U];
    if(s0 >= flat_threshold && s1 >= flat_threshold) {
      continue;
    }
    AutoCfdSurfaceCurvatureEdgeSource edge_source;
    edge_source.positions[0] = sampled_curve.samples[i].position;
    edge_source.positions[1] = sampled_curve.samples[i + 1U].position;
    edge_source.sizes[0] = s0;
    edge_source.sizes[1] = s1;
    mutable_field.curvature_edge_sources.push_back(edge_source);
  }
}

[[nodiscard]] std::size_t face_sample_count(
  const AutoCfdSurfaceFaceWorkItem &face_work_item,
  const std::unordered_map<std::uint32_t, double> &curve_lengths
) noexcept
{
  std::unordered_set<std::uint32_t> visited_edges;
  double boundary_length = 0.0;

  for(const auto &loop : face_work_item.boundary.loops) {
    for(const auto &edge_use : loop.edge_uses) {
      if(edge_use.edge.dimension != geo::TopologyDimension::edge) {
        continue;
      }
      if(!visited_edges.insert(edge_use.edge.index).second) {
        continue;
      }

      const auto it = curve_lengths.find(edge_use.edge.index);
      if(it != curve_lengths.end()) {
        boundary_length += it->second;
      }
    }
  }

  if(boundary_length <= 0.0 || face_work_item.target_size <= 0.0) {
    return kMinimumFaceSampleCount;
  }

  const auto sample_count = static_cast<std::size_t>(std::ceil(
    boundary_length / (4.0 * face_work_item.target_size)
  ));
  return std::clamp(
    sample_count,
    kMinimumFaceSampleCount,
    kMaximumFaceSampleCount
  );
}

/// Extract proxy mesh node positions, triangle connectivity, and per-node
/// topology owners from the proxy domain.  Returns false if the domain has
/// no usable proxy triangles.
[[nodiscard]] bool extract_proxy_mesh_arrays(
  const mesh::Domain &proxy_domain,
  std::vector<geo::Point3> &out_nodes,
  std::vector<std::array<std::uint32_t, 3>> &out_triangles,
  std::vector<geo::TopologyEntityId> &out_node_owners
)
{
  out_nodes.clear();
  out_triangles.clear();
  out_node_owners.clear();

  // Map from (entity_group_id, local_index) to flat node index.
  std::unordered_map<std::uint64_t, std::uint32_t> node_map;

  auto get_or_insert_node = [&](mesh::EntityRef ref,
                                const mesh::Domain &dom) -> std::uint32_t {
    const std::uint64_t key =
      (static_cast<std::uint64_t>(ref.entity_group) << 32U) | ref.index;
    auto [it, inserted] = node_map.emplace(key, static_cast<std::uint32_t>(out_nodes.size()));
    if(inserted) {
      out_nodes.push_back(dom.node(ref).coordinates);
      // Owner will be assigned later from incident triangle face owners.
      out_node_owners.push_back({});
    }
    return it->second;
  };

  for(const auto &entity_group : proxy_domain.entity_groups()) {
    if(entity_group.role() != mesh::EntityGroupRole::geometric_proxy ||
       entity_group.order() != mesh::EntityOrder::face) {
      continue;
    }

    for(std::uint32_t fi = 0U; fi < entity_group.faces().size(); ++fi) {
      const mesh::EntityRef face_ref {entity_group.id(), fi};
      const auto face_nodes = proxy_domain.face_nodes(face_ref);
      if(face_nodes.size != 3U) {
        continue;
      }

      const auto n0 = get_or_insert_node(face_nodes[0], proxy_domain);
      const auto n1 = get_or_insert_node(face_nodes[1], proxy_domain);
      const auto n2 = get_or_insert_node(face_nodes[2], proxy_domain);
      out_triangles.push_back({n0, n1, n2});

      // Assign face topology owner to nodes (first-write wins for shared nodes).
      const auto face_owner = proxy_domain.face_topology_owner(face_ref);
      if(out_node_owners[n0].index == geo::invalid_topology_index) {
        out_node_owners[n0] = face_owner;
      }
      if(out_node_owners[n1].index == geo::invalid_topology_index) {
        out_node_owners[n1] = face_owner;
      }
      if(out_node_owners[n2].index == geo::invalid_topology_index) {
        out_node_owners[n2] = face_owner;
      }
    }
  }

  return !out_nodes.empty() && !out_triangles.empty();
}

[[nodiscard]] base::StatusCode populate_curvature_sources(
  const geo::ModelView &model_view,
  const std::vector<AutoCfdSurfaceCurveWorkItem> &curve_work_items,
  const std::vector<AutoCfdSurfaceFaceWorkItem> & /*face_work_items*/,
  AutoCfdSurfaceSizingFieldState &field
)
{
  field.curvature_sources.clear();
  field.curvature_face_sources.clear();
  field.curvature_edge_sources.clear();

  // Curve curvature: each GEdge contributes consecutive-pair line-segment
  // edge sources for continuous coverage along the curve.
  for(const auto &curve_work_item : curve_work_items) {
    append_edge_curvature_sources(model_view, curve_work_item, field, field);
  }

  const auto edge_source_count = field.curvature_edge_sources.size();

  auto storage = core::detail::lookup_model_storage(
    model_view.model_handle, model_view.context_handle);
  if(!storage) {
    return core::detail::clear_error_state();
  }

  model::detail::GeometryCoarseProxyMesh proxy_data;
  if(storage->coarse_proxy_mesh(proxy_data) != base::StatusCode::ok ||
     proxy_data.empty()) {
    return core::detail::clear_error_state();
  }

  ProxyMeshGeometry geom;
  geom.nodes = &proxy_data.nodes;
  geom.triangles = &proxy_data.triangles;
  geom.face_triangle_ranges = &proxy_data.face_triangle_ranges;
  geom.compute_normals_and_areas();
  geom.compute_chord_curvature(field.distortion_angle);

  const auto nv = proxy_data.nodes.size();
  const auto nt = proxy_data.triangles.size();

  const double theta = distortion_angle_radians(field.distortion_angle);
  const double kappa_scale = kCurvatureCentroidFactor * theta;

  // Per-triangle target size. Flat triangles (κ ≤ tolerance or dihedral
  // filtered out upstream) keep maximum_length, which is a no-op in the
  // per-node MIN step below.
  std::vector<double> tri_msize(nt, field.maximum_length);
  for(std::size_t fi = 0; fi < nt; ++fi) {
    const double kmax = geom.triangle_curvatures[fi];
    if(kmax <= kCurvatureTolerance || !std::isfinite(kmax)) {
      continue;
    }
    const double radius = 1.0 / kmax;
    tri_msize[fi] = std::clamp(
      radius * kappa_scale, field.minimum_length, field.maximum_length);
  }

  // Per-node MIN across incident triangles.
  std::vector<double> node_msize(nv, field.maximum_length);
  for(std::size_t fi = 0; fi < nt; ++fi) {
    const double ms = tri_msize[fi];
    if(ms >= field.maximum_length) {
      continue;
    }
    const auto &tri = proxy_data.triangles[fi];
    for(int k = 0; k < 3; ++k) {
      const auto vi = tri[k];
      if(ms < node_msize[vi]) {
        node_msize[vi] = ms;
      }
    }
  }

  // Emit face source with 3 per-vertex sizes; skip triangles where every
  // vertex is at/near max (nothing to refine).
  const double flat_threshold = 0.9 * field.maximum_length;
  field.curvature_face_sources.reserve(nt);
  std::size_t curvature_tri_count = 0U;
  for(std::size_t fi = 0; fi < nt; ++fi) {
    const auto &tri = proxy_data.triangles[fi];
    const double s0 = node_msize[tri[0]];
    const double s1 = node_msize[tri[1]];
    const double s2 = node_msize[tri[2]];
    if(s0 >= flat_threshold && s1 >= flat_threshold && s2 >= flat_threshold) {
      continue;
    }
    AutoCfdSurfaceCurvatureFaceSource face_source;
    face_source.positions[0] = proxy_data.nodes[tri[0]];
    face_source.positions[1] = proxy_data.nodes[tri[1]];
    face_source.positions[2] = proxy_data.nodes[tri[2]];
    face_source.sizes[0] = s0;
    face_source.sizes[1] = s1;
    face_source.sizes[2] = s2;
    field.curvature_face_sources.push_back(face_source);
    ++curvature_tri_count;
  }

  SQMESH_LOG_INFO("Face curvature: {} triangles, {} curvature "
    "triangles, {} face sources, {} edge sources",
    nt, curvature_tri_count,
    field.curvature_face_sources.size(), edge_source_count);

  return core::detail::clear_error_state();
}

/// Growth-rate-aware proximity source deduplication.
///
/// Incremental processing: sources are processed one at a time, updating the
/// background grid after each. Once a small proximity source is processed, its
/// effect propagates outward via growth rate. Subsequent queries at nearby
/// points see the reduced size → smaller search range → fewer detections.
/// This self-limits proximity to the regions where it truly matters (e.g.,
/// trailing edge) and prevents it from spreading across the entire surface.
///
/// We simulate this by sorting sources by target_size ascending, then discarding
/// any source that is already "covered" by an accepted source's growth-rate
/// propagation.  If an accepted source at distance D with target S can reach
/// the candidate position with effective size S + growth_slope * D <= T, the
/// candidate (target T) is redundant.
void deduplicate_proximity_sources(
  std::vector<AutoCfdSurfaceSizingSource> &sources,
  double growth_rate
)
{
  if(sources.size() <= 1U) {
    return;
  }

  // Sort ascending by relaxed_target_size: smallest (most important) first.
  std::sort(sources.begin(), sources.end(),
    [](const AutoCfdSurfaceSizingSource &a, const AutoCfdSurfaceSizingSource &b) {
      return a.relaxed_target_size < b.relaxed_target_size;
    }
  );

  // Spatial deduplication: an accepted source at size S covers a radius of S
  // around itself.  Any candidate within this radius is redundant because the
  // octree size function will already produce a similar-or-smaller size at
  // that location via gradient propagation from the accepted source.
  static_cast<void>(growth_rate);
  std::vector<AutoCfdSurfaceSizingSource> accepted;
  accepted.reserve(sources.size());

  for(auto &source : sources) {
    bool covered = false;
    for(const auto &acc : accepted) {
      // Coverage radius uses relaxed_target_size (the effective size fed into
      // the size function) for accurate spatial coverage estimation.
      const double cover_radius = acc.relaxed_target_size;
      const double cover_radius_sq = cover_radius * cover_radius;

      const double dist_sq =
        (source.position[0] - acc.position[0]) * (source.position[0] - acc.position[0]) +
        (source.position[1] - acc.position[1]) * (source.position[1] - acc.position[1]) +
        (source.position[2] - acc.position[2]) * (source.position[2] - acc.position[2]);

      if(dist_sq <= cover_radius_sq) {
        covered = true;
        break;
      }
    }

    if(!covered) {
      accepted.push_back(std::move(source));
    }
  }

  sources = std::move(accepted);
}

[[nodiscard]] base::StatusCode populate_proximity_sources(
  const geo::ModelView &model_view,
  const std::vector<AutoCfdSurfaceCurveWorkItem> &curve_work_items,
  const std::vector<AutoCfdSurfaceFaceWorkItem> &face_work_items,
  const AutoCfdSurfaceSizingFieldState &base_field,
  AutoCfdSurfaceSizingFieldState &field
)
{
  field.proximity_sources.clear();
  field.proximity_stats.proxy_triangle_count = 0U;
  field.proximity_stats.bvh_node_count = 0U;
  field.proximity_stats.bvh_leaf_count = 0U;
  field.proximity_stats.bvh_max_leaf_triangle_count = 0U;
  field.proximity_stats.face_sample_count = 0U;
  field.proximity_stats.edge_sample_count = 0U;
  if(!field.proximity_enabled ||
     (!field.self_proximity_enabled && !field.pid_proximity_enabled) ||
     !std::isfinite(field.proximity_length_to_gap_ratio) ||
     field.proximity_length_to_gap_ratio <= 0.0) {
    return core::detail::clear_error_state();
  }

  mesh::MeshHandle proxy_mesh = sqmesh::invalid_handle;
  auto status = core::detail::model_proxy_mesh(
    model_view.model_handle,
    proxy_mesh,
    model_view.context_handle
  );
  if(status != base::StatusCode::ok) {
    if(status == base::StatusCode::unsupported) {
      return core::detail::clear_error_state();
    }
    return status;
  }
  if(proxy_mesh == sqmesh::invalid_handle) {
    return core::detail::clear_error_state();
  }

  std::vector<AutoCfdSurfaceProxyTriangle> proxy_triangles;
  status = core::detail::with_mesh_domain(
    proxy_mesh,
    [&](const mesh::Domain &proxy_domain) {
      if(proxy_domain.source_topology_revision() != model_view.snapshot.topology_revision) {
        return core::detail::publish_error(
          base::StatusCode::internal_error,
          "Auto CFD proximity queries require a proxy domain synchronized to the current topology revision."
        );
      }
      return collect_proxy_triangle_query_data(proxy_domain, field, proxy_triangles);
    },
    model_view.context_handle
  );
  if(status != base::StatusCode::ok || proxy_triangles.empty()) {
    return status;
  }

  AutoCfdSurfaceProximityIndex proximity_index;
  proximity_index.build(std::move(proxy_triangles));
  field.proximity_stats.proxy_triangle_count = proximity_index.stats().triangle_count;
  field.proximity_stats.bvh_node_count = proximity_index.stats().node_count;
  field.proximity_stats.bvh_leaf_count = proximity_index.stats().leaf_count;
  field.proximity_stats.bvh_max_leaf_triangle_count =
    proximity_index.stats().max_leaf_triangle_count;
  if(proximity_index.empty()) {
    return core::detail::clear_error_state();
  }

  // Face proximity: iterate proxy triangles, raycast from centroid,
  // inject whole triangle as face source for continuous coverage.
  const double proximity_floor = std::max(
    field.minimum_length,
    field.proximity_minimum_length
  );
  for(const auto &tri : proximity_index.triangles()) {
    const double tri_normal_norm = vector_norm(tri.normal);
    if(tri_normal_norm <= 0.0) {
      continue;
    }

    const double source_size = query_auto_cfd_surface_sizing_field(
      base_field,
      tri.owner,
      tri.centroid
    );
    if(!std::isfinite(source_size) || source_size <= proximity_floor) {
      continue;
    }

    // Global search range: max_msize * cell_per_gap * 1.5
    // (not per-face source_size * 3, which is too small near high-curvature
    // edges like leading edges where curvature shrinks source_size).
    // Cap by source_size * 3 to avoid excessive mid-chord detections.
    const double global_range =
      field.maximum_length / field.proximity_length_to_gap_ratio * 1.5;
    const double local_range =
      source_size / field.proximity_length_to_gap_ratio;
    const double max_hit_distance = std::min(global_range, local_range);
    if(!std::isfinite(max_hit_distance) || max_hit_distance <= 0.0) {
      continue;
    }

    const geo::Vector3 source_normal = {
      tri.normal[0] / tri_normal_norm,
      tri.normal[1] / tri_normal_norm,
      tri.normal[2] / tri_normal_norm,
    };

    ++field.proximity_stats.face_sample_count;

    ProximityHitCandidate self_hit;
    ProximityHitCandidate pid_hit;

    query_proximity_candidates(
      proximity_index,
      tri.centroid,
      source_normal,
      tri.owner,
      tri.part_index,
      max_hit_distance,
      field.proximity_normal_angle_threshold_degrees,
      field.self_proximity_enabled,
      field.pid_proximity_enabled,
      self_hit,
      pid_hit
    );

    // Debug: count ray hits before filtering.
    if(self_hit.found || pid_hit.found) {
      ++field.proximity_stats.ray_hit_count;
    }

    // min_gap filter: ignore gaps smaller than min_size * layer_count.
    // This prevents proximity from generating target sizes below the minimum
    // element size (min_gap = min_msize * layer_factor).
    const double min_gap = proximity_floor * static_cast<double>(kProximityLayerCount);

    for(const auto *candidate : {&self_hit, &pid_hit}) {
      if(!candidate->found) {
        continue;
      }
      if(candidate->hit_distance < min_gap) {
        continue;
      }
      // Significance filter: proximity target must be below source_size.
      const double candidate_target_size = std::max(
        candidate->hit_distance * field.proximity_length_to_gap_ratio,
        field.proximity_minimum_length
      );
      if(candidate_target_size >= source_size) {
        continue;
      }
      append_proximity_source(
        field,
        tri.owner,
        candidate->kind,
        candidate->hit_owner,
        tri.centroid,
        0.0,
        0.0,
        candidate->hit_distance,
        candidate->hit_normal_angle_degrees,
        source_size
      );
    }
  }

  // Edge proximity: keep raycast-based point sources for boundary coverage.
  for(const auto &curve_work_item : curve_work_items) {
    append_edge_proximity_sources(
      model_view,
      curve_work_item,
      base_field,
      proximity_index,
      field
    );
  }

  // Deduplication: approximate incremental source processing by removing
  // spatially redundant sources. This prevents multiple proxy triangles
  // detecting the same gap from creating duplicate octree seeds.
  const std::size_t raw_count = field.proximity_sources.size();
  deduplicate_proximity_sources(field.proximity_sources, field.growth_rate);

  // Compute proximity source statistics for diagnostics.
  double prox_min_raw = std::numeric_limits<double>::infinity();
  double prox_max_raw = 0.0;
  double prox_min_relaxed = std::numeric_limits<double>::infinity();
  double prox_max_relaxed = 0.0;
  double prox_sum_relaxed = 0.0;
  double prox_min_gap = std::numeric_limits<double>::infinity();
  double prox_max_gap = 0.0;
  for(const auto &src : field.proximity_sources) {
    prox_min_raw = std::min(prox_min_raw, src.target_size);
    prox_max_raw = std::max(prox_max_raw, src.target_size);
    prox_min_relaxed = std::min(prox_min_relaxed, src.relaxed_target_size);
    prox_max_relaxed = std::max(prox_max_relaxed, src.relaxed_target_size);
    prox_sum_relaxed += src.relaxed_target_size;
    prox_min_gap = std::min(prox_min_gap, src.hit_distance);
    prox_max_gap = std::max(prox_max_gap, src.hit_distance);
  }
  const double prox_avg_relaxed = field.proximity_sources.empty()
    ? 0.0 : prox_sum_relaxed / static_cast<double>(field.proximity_sources.size());

  SQMESH_LOG_INFO("raw sources: {}, after dedup: {}, "
    "proxy triangles: {}, face samples: {}, edge samples: {}, "
    "ray hits: {}",
    raw_count,
    field.proximity_sources.size(),
    field.proximity_stats.proxy_triangle_count,
    field.proximity_stats.face_sample_count,
    field.proximity_stats.edge_sample_count,
    field.proximity_stats.ray_hit_count
  );
  SQMESH_LOG_INFO("raw_target: min={:.4f} max={:.4f}, "
    "relaxed_target: min={:.4f} avg={:.4f} max={:.4f}, "
    "gap: min={:.4f} max={:.4f}",
    prox_min_raw, prox_max_raw,
    prox_min_relaxed, prox_avg_relaxed, prox_max_relaxed,
    prox_min_gap, prox_max_gap
  );

  // Dump first few proximity sources for position debugging.
  {
    const std::size_t dump_count = std::min<std::size_t>(
      field.proximity_sources.size(), 10U);
    for(std::size_t i = 0U; i < dump_count; ++i) {
      const auto &src = field.proximity_sources[i];
    }
  }

  return core::detail::clear_error_state();
}

void configure_background_grid(
  const AutoCfdSurfaceSizingFieldState &field,
  const MutableBackgroundBounds &support_bounds,
  AutoCfdSurfaceBackgroundGrid &grid
)
{
  grid = {};
  if(!support_bounds.defined) {
    return;
  }

  const double padding =
    std::max(field.minimum_length, field.maximum_length * kBackgroundGridPaddingFactor);
  const double minimum_spacing = std::max(field.minimum_length, 1.0e-9);

  grid.minimum = support_bounds.minimum;
  grid.maximum = support_bounds.maximum;
  for(std::size_t axis = 0U; axis < 3U; ++axis) {
    const double axis_extent = std::max(
      field.minimum_length,
      support_bounds.maximum[axis] - support_bounds.minimum[axis]
    );
    grid.origin[axis] = support_bounds.minimum[axis] - padding;
    const double padded_extent = axis_extent + 2.0 * padding;
    grid.spacing[axis] = minimum_spacing;
    grid.node_counts[axis] =
      clamped_background_grid_node_count(padded_extent, minimum_spacing);
  }

  double scale = 1.0;
  for(std::size_t axis = 0U; axis < 3U; ++axis) {
    if(grid.node_counts[axis] > kBackgroundGridMaximumAxisNodes) {
      scale = std::max(
        scale,
        static_cast<double>(grid.node_counts[axis] - 1U) /
          static_cast<double>(kBackgroundGridMaximumAxisNodes - 1U)
      );
    }
  }

  const double raw_node_count = static_cast<double>(grid.node_counts[0]) *
                                static_cast<double>(grid.node_counts[1]) *
                                static_cast<double>(grid.node_counts[2]);
  if(raw_node_count > static_cast<double>(kBackgroundGridMaximumNodeCount)) {
    scale = std::max(
      scale,
      std::cbrt(raw_node_count / static_cast<double>(kBackgroundGridMaximumNodeCount))
    );
  }

  if(scale > 1.0) {
    for(std::size_t axis = 0U; axis < 3U; ++axis) {
      const double axis_extent = std::max(
        field.minimum_length,
        support_bounds.maximum[axis] - support_bounds.minimum[axis]
      ) + 2.0 * padding;
      grid.spacing[axis] = minimum_spacing * scale;
      grid.node_counts[axis] =
        std::min<std::size_t>(
          kBackgroundGridMaximumAxisNodes,
          clamped_background_grid_node_count(axis_extent, grid.spacing[axis])
        );
    }
  }

  const std::size_t total_node_count =
    grid.node_counts[0] * grid.node_counts[1] * grid.node_counts[2];
  grid.node_values.assign(total_node_count, field.maximum_length);
  const std::size_t total_cell_count =
    (grid.node_counts[0] - 1U) * (grid.node_counts[1] - 1U) * (grid.node_counts[2] - 1U);
  grid.seed_sources.clear();
  grid.cell_source_indices.assign(total_cell_count, {});
  grid.built = !grid.node_values.empty();
}

void seed_background_grid_from_sources(
  const AutoCfdSurfaceSizingFieldState &field,
  AutoCfdSurfaceBackgroundGrid &grid,
  std::priority_queue<
    BackgroundGridQueueEntry,
    std::vector<BackgroundGridQueueEntry>,
    BackgroundGridQueueGreater
  > &queue
)
{
  if(!grid.built) {
    return;
  }

  const auto seed_sources =
    [&](const std::vector<AutoCfdSurfaceSizingSource> &sources) {
      for(const auto &source : sources) {
        if(!std::isfinite(source.target_size) || source.target_size <= 0.0) {
          continue;
        }

        ++grid.seed_count;
        const auto seed_axis =
          [&](std::size_t axis) {
            const double spacing = grid.spacing[axis];
            if(grid.node_counts[axis] <= 1U || !std::isfinite(spacing) || spacing <= 0.0) {
              return std::pair<std::size_t, std::size_t> {0U, 0U};
            }

            const double coordinate = std::clamp(
              (source.position[axis] - grid.origin[axis]) / spacing,
              0.0,
              static_cast<double>(grid.node_counts[axis] - 1U)
            );
            const std::size_t lower = static_cast<std::size_t>(std::floor(coordinate));
            const std::size_t upper = std::min(lower + 1U, grid.node_counts[axis] - 1U);
            return std::pair<std::size_t, std::size_t> {lower, upper};
          };
        const auto cell_axis =
          [&](std::size_t axis) -> std::size_t {
            const double spacing = grid.spacing[axis];
            if(grid.node_counts[axis] <= 1U || !std::isfinite(spacing) || spacing <= 0.0) {
              return 0U;
            }

            const double coordinate = std::clamp(
              (source.position[axis] - grid.origin[axis]) / spacing,
              0.0,
              static_cast<double>(grid.node_counts[axis] - 1U - 1U)
            );
            return static_cast<std::size_t>(std::floor(coordinate));
          };

        const auto x = seed_axis(0U);
        const auto y = seed_axis(1U);
        const auto z = seed_axis(2U);
        const std::size_t cell_index = background_grid_cell_index(
          grid,
          cell_axis(0U),
          cell_axis(1U),
          cell_axis(2U)
        );
        if(cell_index < grid.cell_source_indices.size()) {
          grid.cell_source_indices[cell_index].push_back(
            static_cast<std::uint32_t>(grid.seed_sources.size())
          );
        }
        grid.seed_sources.push_back(source);

        for(const std::size_t i : {x.first, x.second}) {
          for(const std::size_t j : {y.first, y.second}) {
            for(const std::size_t k : {z.first, z.second}) {
              const std::size_t node_index = background_grid_node_index(grid, i, j, k);
              const double candidate =
                clamp_auto_cfd_surface_size(field, source.target_size);
              if(candidate + kBackgroundGridValueTolerance >= grid.node_values[node_index]) {
                continue;
              }
              grid.node_values[node_index] = candidate;
              queue.push({candidate, node_index});
            }
          }
        }
      }
    };

  seed_sources(field.curvature_sources);
  seed_sources(field.proximity_sources);
}

void diffuse_background_grid(
  const AutoCfdSurfaceSizingFieldState &field,
  AutoCfdSurfaceBackgroundGrid &grid
)
{
  if(!grid.built || grid.node_values.empty()) {
    return;
  }

  const double growth_slope = std::max(0.0, field.growth_rate - 1.0);
  std::priority_queue<
    BackgroundGridQueueEntry,
    std::vector<BackgroundGridQueueEntry>,
    BackgroundGridQueueGreater
  > queue;
  seed_background_grid_from_sources(field, grid, queue);
  if(queue.empty()) {
    return;
  }

  while(!queue.empty()) {
    const auto current = queue.top();
    queue.pop();
    if(current.node_index >= grid.node_values.size() ||
       current.value > grid.node_values[current.node_index] + kBackgroundGridValueTolerance) {
      continue;
    }

    const auto ijk = background_grid_node_ijk(grid, current.node_index);
    for(int dk = -1; dk <= 1; ++dk) {
      for(int dj = -1; dj <= 1; ++dj) {
        for(int di = -1; di <= 1; ++di) {
          if(di == 0 && dj == 0 && dk == 0) {
            continue;
          }

          const int neighbor_i = static_cast<int>(ijk[0]) + di;
          const int neighbor_j = static_cast<int>(ijk[1]) + dj;
          const int neighbor_k = static_cast<int>(ijk[2]) + dk;
          if(neighbor_i < 0 || neighbor_j < 0 || neighbor_k < 0 ||
             neighbor_i >= static_cast<int>(grid.node_counts[0]) ||
             neighbor_j >= static_cast<int>(grid.node_counts[1]) ||
             neighbor_k >= static_cast<int>(grid.node_counts[2])) {
            continue;
          }

          const std::size_t neighbor_index = background_grid_node_index(
            grid,
            static_cast<std::size_t>(neighbor_i),
            static_cast<std::size_t>(neighbor_j),
            static_cast<std::size_t>(neighbor_k)
          );
          const double step_distance = std::sqrt(
            squared_value(grid.spacing[0] * static_cast<double>(di)) +
            squared_value(grid.spacing[1] * static_cast<double>(dj)) +
            squared_value(grid.spacing[2] * static_cast<double>(dk))
          );
          const double candidate = clamp_auto_cfd_surface_size(
            field,
            current.value + growth_slope * step_distance
          );
          if(candidate + kBackgroundGridValueTolerance >= grid.node_values[neighbor_index]) {
            continue;
          }
          grid.node_values[neighbor_index] = candidate;
          queue.push({candidate, neighbor_index});
        }
      }
    }
  }
}

void populate_relaxed_source_sizes(AutoCfdSurfaceSizingFieldState &field)
{
  for(auto &source : field.curvature_sources) {
    source.relaxed_target_size = query_auto_cfd_surface_sizing_field(
      field,
      source.owner,
      source.position
    );
  }

  for(auto &source : field.proximity_sources) {
    source.relaxed_target_size = query_auto_cfd_surface_sizing_field(
      field,
      source.owner,
      source.position
    );
  }
}

} // namespace

base::StatusCode resolve_auto_cfd_surface_parameters(
  const ParameterDictionary &parameters,
  AutoCfdSurfaceParameters &resolved
)
{
  resolved = {};

  auto status =
    require_numeric_parameter(parameters, "minimum_length", resolved.minimum_length);
  if(status != base::StatusCode::ok) {
    return status;
  }

  status = require_numeric_parameter(parameters, "maximum_length", resolved.maximum_length);
  if(status != base::StatusCode::ok) {
    return status;
  }

  status =
    require_numeric_parameter(parameters, "distortion_angle", resolved.distortion_angle);
  if(status != base::StatusCode::ok) {
    return status;
  }

  status = require_numeric_parameter(parameters, "growth_rate", resolved.growth_rate);
  if(status != base::StatusCode::ok) {
    return status;
  }

  if(resolved.minimum_length <= 0.0) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Parameter 'minimum_length' must be positive."
    );
  }
  if(resolved.maximum_length <= 0.0) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Parameter 'maximum_length' must be positive."
    );
  }
  if(resolved.maximum_length < resolved.minimum_length) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Parameter 'maximum_length' must be greater than or equal to 'minimum_length'."
    );
  }
  if(resolved.distortion_angle <= 0.0 || resolved.distortion_angle >= 180.0) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Parameter 'distortion_angle' must stay within the open interval (0, 180)."
    );
  }
  if(resolved.growth_rate < 1.0) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Parameter 'growth_rate' must be greater than or equal to 1.0."
    );
  }

  bool proximity_specified = false;
  status = try_get_optional_boolean_parameter(
    parameters,
    "proximity",
    resolved.proximity,
    proximity_specified
  );
  if(status != base::StatusCode::ok) {
    return status;
  }
  resolved.self_proximity = resolved.proximity;
  resolved.pid_proximity = resolved.proximity;
  bool self_proximity_specified = false;
  status = try_get_optional_boolean_parameter(
    parameters,
    "self_proximity",
    resolved.self_proximity,
    self_proximity_specified
  );
  if(status != base::StatusCode::ok) {
    return status;
  }
  bool pid_proximity_specified = false;
  status = try_get_optional_boolean_parameter(
    parameters,
    "pid_proximity",
    resolved.pid_proximity,
    pid_proximity_specified
  );
  if(status != base::StatusCode::ok) {
    return status;
  }
  if(self_proximity_specified || pid_proximity_specified) {
    if(!self_proximity_specified) {
      resolved.self_proximity = resolved.proximity;
    }
    if(!pid_proximity_specified) {
      resolved.pid_proximity = resolved.proximity;
    }
    resolved.proximity = resolved.self_proximity || resolved.pid_proximity;
  }

  resolved.proximity_maximum_normals_angle =
    180.0 - kProximityNormalAngleThresholdDegrees;
  resolved.proximity_length_to_gap_ratio =
    1.0 / static_cast<double>(kProximityLayerCount);
  resolved.proximity_minimum_length = resolved.minimum_length;

  bool proximity_maximum_normals_angle_specified = false;
  status = try_get_optional_numeric_parameter(
    parameters,
    "maximum_normals_angle",
    resolved.proximity_maximum_normals_angle,
    proximity_maximum_normals_angle_specified
  );
  if(status != base::StatusCode::ok) {
    return status;
  }
  bool proximity_length_to_gap_ratio_specified = false;
  status = try_get_optional_numeric_parameter(
    parameters,
    "length_to_gap_ratio",
    resolved.proximity_length_to_gap_ratio,
    proximity_length_to_gap_ratio_specified
  );
  if(status != base::StatusCode::ok) {
    return status;
  }
  bool proximity_minimum_length_specified = false;
  status = try_get_optional_numeric_parameter(
    parameters,
    "proximity_minimum_length",
    resolved.proximity_minimum_length,
    proximity_minimum_length_specified
  );
  if(status != base::StatusCode::ok) {
    return status;
  }
  if(resolved.proximity_maximum_normals_angle <= 0.0 ||
     resolved.proximity_maximum_normals_angle >= 180.0) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Parameter 'maximum_normals_angle' must stay within the open interval (0, 180)."
    );
  }
  if(resolved.proximity_length_to_gap_ratio <= 0.0) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Parameter 'length_to_gap_ratio' must be positive."
    );
  }
  if(resolved.proximity_minimum_length <= 0.0) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Parameter 'proximity_minimum_length' must be positive."
    );
  }

  resolved.spacing_growth_rate = resolved.growth_rate;
  resolved.spacing_minimum_length = resolved.minimum_length;
  resolved.spacing_maximum_length = resolved.maximum_length;
  resolved.spacing_feature_angle = resolved.distortion_angle;
  resolved.spacing_sharp_angle_limit = 0.0;
  resolved.spacing_sharp_angle_length = 0.0;
  resolved.spacing_proximity = false;
  resolved.spacing_self_proximity = false;
  resolved.spacing_pid_proximity = false;
  resolved.spacing_maximum_normals_angle =
    180.0 - kProximityNormalAngleThresholdDegrees;
  resolved.spacing_length_to_gap_ratio =
    1.0 / static_cast<double>(kProximityLayerCount);
  resolved.spacing_proximity_minimum_length = resolved.spacing_minimum_length;

  bool spacing_growth_rate_specified = false;
  status = try_get_optional_numeric_parameter(
    parameters,
    "spacing_growth_rate",
    resolved.spacing_growth_rate,
    spacing_growth_rate_specified
  );
  if(status != base::StatusCode::ok) {
    return status;
  }

  bool spacing_proximity_specified = false;
  status = try_get_optional_boolean_parameter(
    parameters,
    "spacing_proximity",
    resolved.spacing_proximity,
    spacing_proximity_specified
  );
  if(status != base::StatusCode::ok) {
    return status;
  }
  bool spacing_self_proximity_specified = false;
  status = try_get_optional_boolean_parameter(
    parameters,
    "spacing_self_proximity",
    resolved.spacing_self_proximity,
    spacing_self_proximity_specified
  );
  if(status != base::StatusCode::ok) {
    return status;
  }
  bool spacing_pid_proximity_specified = false;
  status = try_get_optional_boolean_parameter(
    parameters,
    "spacing_pid_proximity",
    resolved.spacing_pid_proximity,
    spacing_pid_proximity_specified
  );
  if(status != base::StatusCode::ok) {
    return status;
  }
  if(spacing_self_proximity_specified || spacing_pid_proximity_specified) {
    if(!spacing_self_proximity_specified) {
      resolved.spacing_self_proximity = resolved.spacing_proximity;
    }
    if(!spacing_pid_proximity_specified) {
      resolved.spacing_pid_proximity = resolved.spacing_proximity;
    }
    resolved.spacing_proximity =
      resolved.spacing_self_proximity || resolved.spacing_pid_proximity;
  }

  bool spacing_maximum_normals_angle_specified = false;
  status = try_get_optional_numeric_parameter(
    parameters,
    "spacing_maximum_normals_angle",
    resolved.spacing_maximum_normals_angle,
    spacing_maximum_normals_angle_specified
  );
  if(status != base::StatusCode::ok) {
    return status;
  }
  bool spacing_length_to_gap_ratio_specified = false;
  status = try_get_optional_numeric_parameter(
    parameters,
    "spacing_length_to_gap_ratio",
    resolved.spacing_length_to_gap_ratio,
    spacing_length_to_gap_ratio_specified
  );
  if(status != base::StatusCode::ok) {
    return status;
  }
  bool spacing_proximity_minimum_length_specified = false;
  status = try_get_optional_numeric_parameter(
    parameters,
    "spacing_proximity_minimum_length",
    resolved.spacing_proximity_minimum_length,
    spacing_proximity_minimum_length_specified
  );
  if(status != base::StatusCode::ok) {
    return status;
  }

  bool spacing_feature_angle_specified = false;
  status = try_get_optional_numeric_parameter(
    parameters,
    "spacing_feature_angle",
    resolved.spacing_feature_angle,
    spacing_feature_angle_specified
  );
  if(status != base::StatusCode::ok) {
    return status;
  }
  if(!spacing_feature_angle_specified) {
    status = try_get_optional_numeric_parameter(
      parameters,
      "feature_angle",
      resolved.spacing_feature_angle,
      spacing_feature_angle_specified
    );
    if(status != base::StatusCode::ok) {
      return status;
    }
  }
  resolved.spacing_feature_angle_specified = spacing_feature_angle_specified;

  bool spacing_minimum_length_specified = false;
  status = try_get_optional_numeric_parameter(
    parameters,
    "spacing_minimum_length",
    resolved.spacing_minimum_length,
    spacing_minimum_length_specified
  );
  if(status != base::StatusCode::ok) {
    return status;
  }

  bool spacing_maximum_length_specified = false;
  status = try_get_optional_numeric_parameter(
    parameters,
    "spacing_maximum_length",
    resolved.spacing_maximum_length,
    spacing_maximum_length_specified
  );
  if(status != base::StatusCode::ok) {
    return status;
  }
  if(!spacing_proximity_minimum_length_specified) {
    resolved.spacing_proximity_minimum_length = resolved.spacing_minimum_length;
  }

  bool spacing_sharp_angle_limit_specified = false;
  status = try_get_optional_numeric_parameter(
    parameters,
    "spacing_sharp_angle_limit",
    resolved.spacing_sharp_angle_limit,
    spacing_sharp_angle_limit_specified
  );
  if(status != base::StatusCode::ok) {
    return status;
  }
  if(!spacing_sharp_angle_limit_specified) {
    status = try_get_optional_numeric_parameter(
      parameters,
      "sharp_angle_limit",
      resolved.spacing_sharp_angle_limit,
      spacing_sharp_angle_limit_specified
    );
    if(status != base::StatusCode::ok) {
      return status;
    }
  }

  bool spacing_sharp_angle_length_specified = false;
  status = try_get_optional_numeric_parameter(
    parameters,
    "spacing_sharp_angle_length",
    resolved.spacing_sharp_angle_length,
    spacing_sharp_angle_length_specified
  );
  if(status != base::StatusCode::ok) {
    return status;
  }
  if(!spacing_sharp_angle_length_specified) {
    status = try_get_optional_numeric_parameter(
      parameters,
      "sharp_angle_length",
      resolved.spacing_sharp_angle_length,
      spacing_sharp_angle_length_specified
    );
    if(status != base::StatusCode::ok) {
      return status;
    }
  }

  resolved.spacing_sharp_angle_specified =
    spacing_sharp_angle_limit_specified || spacing_sharp_angle_length_specified;
  if(spacing_sharp_angle_limit_specified != spacing_sharp_angle_length_specified) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Parameters 'spacing_sharp_angle_limit' and 'spacing_sharp_angle_length' must be provided together."
    );
  }

  if(resolved.spacing_minimum_length <= 0.0) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Parameter 'spacing_minimum_length' must be positive."
    );
  }
  if(resolved.spacing_maximum_length <= 0.0) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Parameter 'spacing_maximum_length' must be positive."
    );
  }
  if(resolved.spacing_maximum_length < resolved.spacing_minimum_length) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Parameter 'spacing_maximum_length' must be greater than or equal to 'spacing_minimum_length'."
    );
  }
  if(resolved.spacing_growth_rate < 1.0) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Parameter 'spacing_growth_rate' must be greater than or equal to 1.0."
    );
  }
  if(spacing_feature_angle_specified &&
     (resolved.spacing_feature_angle <= 0.0 ||
      resolved.spacing_feature_angle >= 180.0)) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Parameter 'spacing_feature_angle' must stay within the open interval (0, 180)."
    );
  }
  if(spacing_sharp_angle_limit_specified &&
     (resolved.spacing_sharp_angle_limit <= 0.0 ||
      resolved.spacing_sharp_angle_limit >= 180.0)) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Parameter 'spacing_sharp_angle_limit' must stay within the open interval (0, 180)."
    );
  }
  if(spacing_sharp_angle_length_specified &&
     resolved.spacing_sharp_angle_length <= 0.0) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Parameter 'spacing_sharp_angle_length' must be positive."
    );
  }
  if(resolved.spacing_maximum_normals_angle <= 0.0 ||
     resolved.spacing_maximum_normals_angle >= 180.0) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Parameter 'spacing_maximum_normals_angle' must stay within the open interval (0, 180)."
    );
  }
  if(resolved.spacing_length_to_gap_ratio <= 0.0) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Parameter 'spacing_length_to_gap_ratio' must be positive."
    );
  }
  if(resolved.spacing_proximity_minimum_length <= 0.0) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Parameter 'spacing_proximity_minimum_length' must be positive."
    );
  }

  bool allow_quality_gate_failure_specified = false;
  status = try_get_optional_boolean_parameter(
    parameters,
    "allow_quality_gate_failure",
    resolved.allow_quality_gate_failure,
    allow_quality_gate_failure_specified
  );
  if(status != base::StatusCode::ok) {
    return status;
  }

  bool allow_final_screen_failure_specified = false;
  status = try_get_optional_boolean_parameter(
    parameters,
    "allow_final_screen_failure",
    resolved.allow_final_screen_failure,
    allow_final_screen_failure_specified
  );
  if(status != base::StatusCode::ok) {
    return status;
  }

  std::string_view element_type;
  if(parameters.try_get_text("element_type", element_type) &&
     element_type != "tri" &&
     element_type != "tria" &&
     element_type != "triangle" &&
     element_type != "tri3") {
    return core::detail::publish_error(
      base::StatusCode::unsupported,
      "Auto CFD Surface Mesher only supports triangular surface elements."
    );
  }

  return core::detail::clear_error_state();
}

base::StatusCode build_auto_cfd_surface_pipeline_state(
  const geo::ModelView &model_view,
  const ResolvedMeshSizeControls &size_controls,
  const AutoCfdSurfaceParameters &parameters,
  AutoCfdSurfacePipelineState &state
)
{
  if(model_view.faces.empty()) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Auto CFD Surface Mesher requires a geometry model with at least one face."
    );
  }

  state = {};
  state.topology_revision = model_view.snapshot.topology_revision;
  state.parameters = parameters;

  // ── Phase 1: Build anchor candidates and curve metadata ─────────────
  //
  // We need edge_spacing_field temporarily for anchor candidates and sharp
  // corner detection, but edge discretization will use state.sizing_field.
  AutoCfdSurfaceSizingFieldState edge_spacing_field;
  auto edge_spacing_controls = size_controls;
  edge_spacing_controls.global_target_size = parameters.spacing_maximum_length;
  initialize_sizing_field(
    edge_spacing_field,
    parameters.spacing_minimum_length,
    parameters.spacing_maximum_length,
    parameters.spacing_feature_angle_specified
      ? parameters.spacing_feature_angle
      : parameters.distortion_angle,
    parameters.spacing_growth_rate,
    parameters.spacing_proximity,
    parameters.spacing_self_proximity,
    parameters.spacing_pid_proximity,
    parameters.spacing_maximum_normals_angle,
    parameters.spacing_length_to_gap_ratio,
    parameters.spacing_proximity_minimum_length
  );
  populate_owner_target_sizes(model_view, edge_spacing_controls, edge_spacing_field);
  populate_face_proximity_buckets(model_view, edge_spacing_field);

  std::vector<AutoCfdSurfaceCurveWorkItem> edge_curve_work_items;
  auto status = populate_curve_work_items(
    model_view,
    edge_spacing_field,
    edge_curve_work_items
  );
  if(status != base::StatusCode::ok) {
    return status;
  }

  std::vector<AutoCfdSurfaceFaceWorkItem> edge_face_work_items;
  populate_face_work_items(model_view, edge_spacing_field, edge_face_work_items);

  status = populate_anchor_candidates(
    model_view,
    edge_curve_work_items,
    parameters.spacing_feature_angle_specified
      ? parameters.spacing_feature_angle
      : kFeatureAngleDegrees,
    state.anchor_candidates
  );
  if(status != base::StatusCode::ok) {
    return status;
  }

  // ── Phase 2: Build main sizing field ────────────────────────────────
  initialize_sizing_field(
    state.sizing_field,
    parameters.minimum_length,
    parameters.maximum_length,
    parameters.distortion_angle,
    parameters.growth_rate,
    parameters.proximity,
    parameters.self_proximity,
    parameters.pid_proximity,
    parameters.proximity_maximum_normals_angle,
    parameters.proximity_length_to_gap_ratio,
    parameters.proximity_minimum_length
  );
  populate_owner_target_sizes(model_view, size_controls, state.sizing_field);
  populate_face_proximity_buckets(model_view, state.sizing_field);

  status = populate_curve_work_items(
    model_view,
    state.sizing_field,
    state.curve_work_items
  );
  if(status != base::StatusCode::ok) {
    return status;
  }

  populate_face_work_items(model_view, state.sizing_field, state.face_work_items);
  populate_background_support_bounds(
    model_view,
    state.curve_work_items,
    state.sizing_field
  );

  status = populate_curvature_sources(
    model_view,
    state.curve_work_items,
    state.face_work_items,
    state.sizing_field
  );
  if(status != base::StatusCode::ok) {
    return status;
  }
  if(parameters.spacing_sharp_angle_specified) {
    populate_sharp_corner_sources(
      model_view,
      state.face_work_items,
      state.curve_work_items,
      state.anchor_candidates,
      parameters.spacing_sharp_angle_limit,
      parameters.spacing_sharp_angle_length,
      state.sizing_field
    );
  }

  status = rebuild_auto_cfd_surface_background_field(state.sizing_field);
  if(status != base::StatusCode::ok) {
    return status;
  }

  const auto curvature_only_field = state.sizing_field;
  status = populate_proximity_sources(
    model_view,
    state.curve_work_items,
    state.face_work_items,
    curvature_only_field,
    state.sizing_field
  );
  if(status != base::StatusCode::ok) {
    return status;
  }

  if(!state.sizing_field.proximity_sources.empty()) {
    status = rebuild_auto_cfd_surface_background_field(state.sizing_field);
    if(status != base::StatusCode::ok) {
      return status;
    }
  }

  // ── Phase 3: Edge discretization using main sizing field ────────────
  //
  // Use state.sizing_field (not edge_spacing_field) so boundary node
  // spacing matches the sizing the surface mesher will use in the interior.
  status = populate_boundary_state(
    model_view,
    state.anchor_candidates,
    state.curve_work_items,
    state.face_work_items,
    state.sizing_field,
    state.boundary
  );
  if(status != base::StatusCode::ok) {
    return status;
  }

  populate_relaxed_source_sizes(state.sizing_field);

  status = populate_face_preprocess_states(
    model_view,
    state.sizing_field,
    state.boundary,
    state.face_work_items,
    state.face_preprocess_states
  );
  if(status != base::StatusCode::ok) {
    return status;
  }

  return core::detail::clear_error_state();
}

base::StatusCode rebuild_auto_cfd_surface_face_preprocess_states(
  const geo::ModelView &model_view,
  AutoCfdSurfacePipelineState &state
)
{
  return populate_face_preprocess_states(
    model_view,
    state.sizing_field,
    state.boundary,
    state.face_work_items,
    state.face_preprocess_states
  );
}

base::StatusCode split_auto_cfd_surface_boundary_edges_for_retry(
  const geo::ModelView &model_view,
  const std::vector<AutoCfdSurfaceBoundaryRetryRequest> &requests,
  AutoCfdSurfacePipelineState &state
)
{
  if(requests.empty()) {
    return core::detail::clear_error_state();
  }

  std::vector<AutoCfdSurfaceBoundaryRetryRequest> unique_requests;
  unique_requests.reserve(requests.size());
  for(const auto &request : requests) {
    const auto edge = request.edge;
    if(edge.dimension != geo::TopologyDimension::edge ||
       edge.index >= state.boundary.edge_discretizations.size() ||
       request.first_node_index == invalid_index ||
       request.second_node_index == invalid_index ||
       request.first_node_index == request.second_node_index) {
      return core::detail::publish_error(
        base::StatusCode::internal_error,
        "Auto CFD shared boundary retry encountered an invalid unrecovered boundary sub-edge request."
      );
    }

    unique_requests.push_back(
      {
        edge,
        std::min(request.first_node_index, request.second_node_index),
        std::max(request.first_node_index, request.second_node_index),
      }
    );
  }
  std::sort(
    unique_requests.begin(),
    unique_requests.end(),
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
  unique_requests.erase(
    std::unique(
      unique_requests.begin(),
      unique_requests.end(),
      [](const AutoCfdSurfaceBoundaryRetryRequest &lhs,
         const AutoCfdSurfaceBoundaryRetryRequest &rhs) noexcept {
        return lhs.edge.dimension == rhs.edge.dimension &&
               lhs.edge.index == rhs.edge.index &&
               lhs.first_node_index == rhs.first_node_index &&
               lhs.second_node_index == rhs.second_node_index;
      }
    ),
    unique_requests.end()
  );

  std::vector<std::size_t> split_segment_indices;
  for(const auto &request : unique_requests) {
    const auto edge = request.edge;
    auto *edge_view = model_view.find_edge(edge);
    if(edge_view == nullptr) {
      return core::detail::publish_error(
        base::StatusCode::internal_error,
        "Auto CFD shared boundary retry could not resolve a topology edge back to the geometry view."
      );
    }

    auto &edge_discretization = state.boundary.edge_discretizations[edge.index];
    if(edge_discretization.points.size() < 2U) {
      return core::detail::publish_error(
        base::StatusCode::unsupported,
        "Auto CFD shared boundary retry cannot split a degenerate boundary edge discretization."
      );
    }

    split_segment_indices.clear();
    for(std::size_t point_index = 0U;
        point_index + 1U < edge_discretization.points.size();
        ++point_index) {
      const auto first_node = std::min(
        edge_discretization.points[point_index].node_index,
        edge_discretization.points[point_index + 1U].node_index
      );
      const auto second_node = std::max(
        edge_discretization.points[point_index].node_index,
        edge_discretization.points[point_index + 1U].node_index
      );
      if(first_node == request.first_node_index &&
         second_node == request.second_node_index) {
        split_segment_indices.push_back(point_index);
      }
    }
    if(split_segment_indices.empty()) {
      return core::detail::publish_error(
        base::StatusCode::unsupported,
        "Auto CFD shared boundary retry could not map an unrecovered constrained sub-edge back to the shared 1D discretization."
      );
    }

    std::vector<AutoCfdSurfaceEdgePoint> refined_points;
    refined_points.reserve(edge_discretization.points.size() +
                           split_segment_indices.size());
    std::size_t split_segment_cursor = 0U;
    for(std::size_t point_index = 0U;
        point_index + 1U < edge_discretization.points.size();
        ++point_index) {
      const auto first = edge_discretization.points[point_index];
      const auto second = edge_discretization.points[point_index + 1U];
      refined_points.push_back(first);

      if(split_segment_cursor >= split_segment_indices.size() ||
         split_segment_indices[split_segment_cursor] != point_index) {
        continue;
      }
      ++split_segment_cursor;

      const double midpoint_parameter = 0.5 * (first.parameter + second.parameter);
      if(!std::isfinite(midpoint_parameter) ||
         midpoint_parameter <= first.parameter ||
         midpoint_parameter >= second.parameter) {
        continue;
      }

      geo::EdgeTangentSample tangent_sample;
      auto status = geo::sample_edge_tangent(*edge_view, midpoint_parameter, tangent_sample);
      if(status != base::StatusCode::ok) {
        return status;
      }

      const auto node_index = add_boundary_node(
        state.boundary,
        tangent_sample.position,
        {},
        nullptr
      );
      refined_points.push_back(
        {
          node_index,
          midpoint_parameter,
          query_auto_cfd_surface_sizing_field(
            state.sizing_field,
            edge,
            tangent_sample.position
          ),
        }
      );
    }

    refined_points.push_back(edge_discretization.points.back());
    edge_discretization.points = std::move(refined_points);
    edge_discretization.short_edge_fallback = false;
  }

  // Only the faces whose boundary uses one of the split edges need their
  // preprocess state rebuilt. Rebuilding every face's preprocess (as the
  // full `populate_face_preprocess_states` does) dominates retry wall time
  // on CRM (~70 s per retry iteration for 30+ faces).
  std::unordered_set<std::uint32_t> split_edge_indices;
  split_edge_indices.reserve(unique_requests.size());
  for(const auto &request : unique_requests) {
    if(request.edge.dimension == geo::TopologyDimension::edge) {
      split_edge_indices.insert(request.edge.index);
    }
  }

  std::unordered_set<std::uint32_t> affected_face_indices;
  for(const auto &face_work_item : state.face_work_items) {
    if(face_work_item.face.dimension != geo::TopologyDimension::face) {
      continue;
    }
    for(const auto &loop : face_work_item.boundary.loops) {
      bool found = false;
      for(const auto &edge_use : loop.edge_uses) {
        if(edge_use.edge.dimension == geo::TopologyDimension::edge &&
           split_edge_indices.find(edge_use.edge.index) !=
             split_edge_indices.end()) {
          affected_face_indices.insert(face_work_item.face.index);
          found = true;
          break;
        }
      }
      if(found) {
        break;
      }
    }
  }

  return rebuild_face_preprocess_states_for_faces(
    model_view,
    state.sizing_field,
    state.boundary,
    state.face_work_items,
    affected_face_indices,
    state.face_preprocess_states
  );
}

const AutoCfdSurfaceFacePreprocessState *find_auto_cfd_surface_face_preprocess_state(
  const AutoCfdSurfacePipelineState &state,
  geo::TopologyEntityId face
) noexcept
{
  if(face.dimension != geo::TopologyDimension::face ||
     face.index >= state.face_preprocess_states.size()) {
    return nullptr;
  }

  const auto &face_preprocess = state.face_preprocess_states[face.index];
  if(face_preprocess.face != face) {
    return nullptr;
  }

  return &face_preprocess;
}

double clamp_auto_cfd_surface_size(
  const AutoCfdSurfaceSizingFieldState &field,
  double size
) noexcept
{
  if(!std::isfinite(size)) {
    return field.maximum_length;
  }

  return std::clamp(size, field.minimum_length, field.maximum_length);
}

double effective_auto_cfd_surface_owner_target_size(
  const AutoCfdSurfaceSizingFieldState &field,
  geo::TopologyEntityId owner
) noexcept
{
  double target_size = field.maximum_length;

  switch(owner.dimension) {
  case geo::TopologyDimension::face:
    if(owner.index < field.face_target_sizes.size()) {
      target_size = field.face_target_sizes[owner.index];
    }
    break;
  case geo::TopologyDimension::edge:
    if(owner.index < field.edge_target_sizes.size()) {
      target_size = field.edge_target_sizes[owner.index];
    }
    break;
  case geo::TopologyDimension::vertex:
  case geo::TopologyDimension::region:
  default:
    break;
  }

  return clamp_auto_cfd_surface_size(field, target_size);
}


base::StatusCode rebuild_auto_cfd_surface_background_field(
  AutoCfdSurfaceSizingFieldState &field
)
{
  const auto bounds = effective_background_bounds(field);
  field.background_grid = {};
  configure_background_octree(field, bounds, field.background_octree);

  // Build the adaptive octree size function for smooth gradation.
  field.size_function.clear();
  SizeFunctionConfig sf_config;
  sf_config.minimum_size = field.minimum_length;
  sf_config.maximum_size = field.maximum_length;
  sf_config.growth_rate = field.growth_rate;
  field.size_function.configure(sf_config);

  if(bounds.defined) {
    field.size_function.set_domain(bounds.minimum, bounds.maximum);
  }

  // Feed face sources for continuous coverage across each triangle.
  for(const auto &face_src : field.curvature_face_sources) {
    field.size_function.add_face_source(
      face_src.positions[0], face_src.sizes[0],
      face_src.positions[1], face_src.sizes[1],
      face_src.positions[2], face_src.sizes[2]
    );
  }
  // Feed line-segment edge sources for continuous coverage along each curve.
  for(const auto &edge_src : field.curvature_edge_sources) {
    field.size_function.add_edge_source(
      edge_src.positions[0], edge_src.sizes[0],
      edge_src.positions[1], edge_src.sizes[1]
    );
  }
  // Point sources (per-sample curve curvature, sharp corners).
  for(const auto &source : field.curvature_sources) {
    if(std::isfinite(source.target_size) && source.target_size > 0.0) {
      field.size_function.add_source(source.position, source.target_size);
    }
  }
  // Proximity point sources — one per detected gap location.
  // Unlike curvature (which needs face sources for continuous coverage),
  // proximity only needs point seeds; the octree growth rate handles
  // smooth transitions from the locally-reduced size back to normal.
  // Proximity sources use relaxed_target_size (capped by curvature-driven
  // relaxation) to limit their influence zone by raising sizes near larger
  // sources.
  for(const auto &source : field.proximity_sources) {
    if(std::isfinite(source.relaxed_target_size) && source.relaxed_target_size > 0.0) {
      field.size_function.add_source(source.position, source.relaxed_target_size);
    }
  }

  field.size_function.build();

  return core::detail::clear_error_state();
}

double query_auto_cfd_surface_sizing_field(
  const AutoCfdSurfaceSizingFieldState &field,
  geo::TopologyEntityId owner,
  const geo::Point3 &position
) noexcept
{
  double target_size = effective_auto_cfd_surface_owner_target_size(field, owner);

  // Primary path: adaptive octree size function (smooth gradation).
  if(field.size_function.built()) {
    const double sf_size = field.size_function.query(position);
    if(std::isfinite(sf_size) && sf_size > 0.0) {
      target_size = std::min(target_size, sf_size);
    }
    return clamp_auto_cfd_surface_size(field, target_size);
  }

  // Legacy fallback: old cone-based octree query.
  const auto &octree = field.background_octree;
  if(octree.built && !octree.nodes.empty() && !octree.seed_sources.empty()) {
    struct OctreeQueueEntry final {
      double lower_bound = std::numeric_limits<double>::infinity();
      std::uint32_t node_index = invalid_auto_cfd_surface_octree_node_index;
    };

    struct OctreeQueueGreater final {
      [[nodiscard]] bool operator()(
        const OctreeQueueEntry &lhs,
        const OctreeQueueEntry &rhs
      ) const noexcept
      {
        return lhs.lower_bound > rhs.lower_bound;
      }
    };

    const double growth_slope = std::max(0.0, field.growth_rate - 1.0);
    std::priority_queue<
      OctreeQueueEntry,
      std::vector<OctreeQueueEntry>,
      OctreeQueueGreater
    > queue;

    const auto root_lower_bound =
      octree.nodes.front().minimum_target_size +
      growth_slope * distance_point_to_aabb(
        position,
        octree.nodes.front().bounds_min,
        octree.nodes.front().bounds_max
      );
    queue.push({root_lower_bound, 0U});

    while(!queue.empty()) {
      const auto current = queue.top();
      queue.pop();
      if(current.node_index >= octree.nodes.size() ||
         current.lower_bound >= target_size - kBackgroundGridValueTolerance) {
        continue;
      }

      const auto &node = octree.nodes[current.node_index];
      if(node.leaf) {
        for(std::uint32_t offset = 0U; offset < node.source_count; ++offset) {
          const std::uint32_t source_slot = node.first_source + offset;
          if(source_slot >= octree.source_indices.size()) {
            continue;
          }
          const auto source_index = octree.source_indices[source_slot];
          if(source_index >= octree.seed_sources.size()) {
            continue;
          }
          const auto &source = octree.seed_sources[source_index];
          const double candidate = clamp_auto_cfd_surface_size(
            field,
            source.target_size + growth_slope * point_distance(source.position, position)
          );
          if(std::isfinite(candidate) && candidate > 0.0) {
            target_size = std::min(target_size, candidate);
          }
        }
        continue;
      }

      for(const auto child_index : node.children) {
        if(child_index == invalid_auto_cfd_surface_octree_node_index ||
           child_index >= octree.nodes.size()) {
          continue;
        }
        const auto &child = octree.nodes[child_index];
        const double lower_bound =
          child.minimum_target_size +
          growth_slope * distance_point_to_aabb(
            position,
            child.bounds_min,
            child.bounds_max
          );
        if(lower_bound >= target_size - kBackgroundGridValueTolerance) {
          continue;
        }
        queue.push({lower_bound, child_index});
      }
    }
  }
  return clamp_auto_cfd_surface_size(field, target_size);
}

double query_auto_cfd_surface_sizing_field(
  const AutoCfdSurfaceSizingFieldState &field,
  const geo::FaceView &face_view,
  double u,
  double v,
  const geo::Point3 &position
) noexcept
{
  // Curvature is now fully captured in pre-computed seed sources and
  // gradient-smoothed.  No on-the-fly OCC curvature query needed.
  return query_auto_cfd_surface_sizing_field(
    field, face_view.entity, position
  );
}

base::StatusCode sample_auto_cfd_surface_face_metric(
  const geo::FaceView &face_view,
  const AutoCfdSurfaceSizingFieldState &field,
  double u,
  double v,
  AutoCfdSurfaceFaceMetricTensor &metric
)
{
  metric = {};

  geo::FaceDerivatives derivatives;
  const auto status = geo::sample_face_derivatives(face_view, u, v, derivatives);
  if(status != base::StatusCode::ok) {
    return status;
  }

  metric.face = face_view.entity;
  metric.u = u;
  metric.v = v;
  metric.first_fundamental_form = {
    dot_product(derivatives.du, derivatives.du),
    dot_product(derivatives.du, derivatives.dv),
    dot_product(derivatives.dv, derivatives.dv),
  };
  metric.first_fundamental_form_defined = derivatives.first_derivatives_defined;
  metric.target_size =
    query_auto_cfd_surface_sizing_field(field, face_view, u, v, derivatives.position);

  if(!derivatives.first_derivatives_defined) {
    build_metric_fallback_from_position(
      face_view.entity,
      field,
      derivatives.position,
      u,
      v,
      &metric.first_fundamental_form,
      AutoCfdSurfaceMetricFallbackKind::isotropic_from_uv_identity,
      metric
    );
    return core::detail::clear_error_state();
  }

  const double inverse_target_size_squared =
    1.0 / std::max(squared_value(metric.target_size), kMetricEigenvalueFloor);
  const double a = metric.first_fundamental_form[0] * inverse_target_size_squared;
  const double b = metric.first_fundamental_form[1] * inverse_target_size_squared;
  const double c = metric.first_fundamental_form[2] * inverse_target_size_squared;
  metric.tensor = {a, b, c};
  metric.raw_eigenvalues = metric_eigenvalues(a, b, c);
  metric.determinant = a * c - b * b;

  if(!std::isfinite(metric.raw_eigenvalues[0]) ||
     !std::isfinite(metric.raw_eigenvalues[1]) ||
     !std::isfinite(metric.determinant) ||
     metric.raw_eigenvalues[1] <= kMetricEigenvalueFloor ||
     metric.raw_eigenvalues[0] <= kMetricEigenvalueFloor ||
     metric.determinant <= kMetricDeterminantTolerance) {
    build_metric_fallback_from_position(
      face_view.entity,
      field,
      derivatives.position,
      u,
      v,
      &metric.first_fundamental_form,
      AutoCfdSurfaceMetricFallbackKind::isotropic_from_jacobian_average,
      metric
    );
    return core::detail::clear_error_state();
  }

  
  auto clamped_eigenvalues = metric.raw_eigenvalues;
  clamped_eigenvalues[1] = std::clamp(
    clamped_eigenvalues[1],
    kMetricEigenvalueFloor,
    kMetricEigenvalueCeiling
  );
  clamped_eigenvalues[0] = std::clamp(
    clamped_eigenvalues[0],
    kMetricEigenvalueFloor,
    clamped_eigenvalues[1]
  );

  metric.clamped =
    std::abs(clamped_eigenvalues[0] - metric.raw_eigenvalues[0]) > kMetricEigenvalueFloor ||
    std::abs(clamped_eigenvalues[1] - metric.raw_eigenvalues[1]) > kMetricEigenvalueFloor;
  metric.clamped_eigenvalues = clamped_eigenvalues;

  const auto primary_eigenvector =
    metric_primary_eigenvector(a, b, c, metric.raw_eigenvalues[1]);
  metric.tensor = reconstruct_metric_tensor(
    primary_eigenvector,
    clamped_eigenvalues[1],
    clamped_eigenvalues[0]
  );
  metric.determinant =
    metric.tensor[0] * metric.tensor[2] - metric.tensor[1] * metric.tensor[1];
  metric.usable =
    std::isfinite(metric.determinant) &&
    metric.determinant > kMetricDeterminantTolerance;
  return core::detail::clear_error_state();
}

} // namespace sqmesh::mesh::detail
