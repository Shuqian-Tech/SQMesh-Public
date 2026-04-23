// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "auto_cfd_surface_proximity_index.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace sqmesh::mesh::detail {
namespace {

constexpr std::size_t kMaximumLeafTriangleCount = 8U;
constexpr std::size_t kMaximumTraversalDepth = 128U;
constexpr double kRayTolerance = 1.0e-9;
constexpr double kPi = 3.14159265358979323846;

[[nodiscard]] double dot_product(
  const geo::Vector3 &lhs,
  const geo::Vector3 &rhs
) noexcept
{
  return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2];
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

[[nodiscard]] double vector_norm(const geo::Vector3 &vector) noexcept
{
  return std::sqrt(dot_product(vector, vector));
}

[[nodiscard]] double angle_degrees_between(
  const geo::Vector3 &lhs,
  const geo::Vector3 &rhs
) noexcept
{
  const double lhs_norm = vector_norm(lhs);
  const double rhs_norm = vector_norm(rhs);
  if(lhs_norm <= 0.0 || rhs_norm <= 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double cosine = std::clamp(dot_product(lhs, rhs) / (lhs_norm * rhs_norm), -1.0, 1.0);
  return std::acos(cosine) * (180.0 / kPi);
}

struct RayBoxInterval final {
  bool hit = false;
  double near_distance = 0.0;
  double far_distance = 0.0;
};

[[nodiscard]] RayBoxInterval ray_box_interval(
  const geo::Point3 &origin,
  const geo::Vector3 &direction,
  double max_distance,
  const geo::Point3 &bounds_min,
  const geo::Point3 &bounds_max
) noexcept
{
  double lower = 0.0;
  double upper = max_distance;

  for(std::size_t axis = 0U; axis < 3U; ++axis) {
    if(std::abs(direction[axis]) <= kRayTolerance) {
      if(origin[axis] < bounds_min[axis] - kRayTolerance ||
         origin[axis] > bounds_max[axis] + kRayTolerance) {
        return {};
      }
      continue;
    }

    const double inverse_direction = 1.0 / direction[axis];
    double axis_lower = (bounds_min[axis] - origin[axis]) * inverse_direction;
    double axis_upper = (bounds_max[axis] - origin[axis]) * inverse_direction;
    if(axis_lower > axis_upper) {
      std::swap(axis_lower, axis_upper);
    }

    lower = std::max(lower, axis_lower);
    upper = std::min(upper, axis_upper);
    if(lower > upper) {
      return {};
    }
  }

  if(upper < 0.0 || lower > max_distance) {
    return {};
  }

  return {
    true,
    std::max(0.0, lower),
    upper,
  };
}

[[nodiscard]] bool intersect_ray_with_triangle(
  const geo::Point3 &origin,
  const geo::Vector3 &direction,
  double max_distance,
  const AutoCfdSurfaceProxyTriangle &triangle,
  double &hit_distance
) noexcept
{
  hit_distance = std::numeric_limits<double>::infinity();

  const auto edge01 = subtract_points(triangle.vertices[1], triangle.vertices[0]);
  const auto edge02 = subtract_points(triangle.vertices[2], triangle.vertices[0]);
  const auto p_vector = cross_product(direction, edge02);
  const double determinant = dot_product(edge01, p_vector);
  if(std::abs(determinant) <= kRayTolerance) {
    return false;
  }

  const double inverse_determinant = 1.0 / determinant;
  const auto t_vector = subtract_points(origin, triangle.vertices[0]);
  const double barycentric_u = dot_product(t_vector, p_vector) * inverse_determinant;
  if(barycentric_u < -kRayTolerance || barycentric_u > 1.0 + kRayTolerance) {
    return false;
  }

  const auto q_vector = cross_product(t_vector, edge01);
  const double barycentric_v = dot_product(direction, q_vector) * inverse_determinant;
  if(barycentric_v < -kRayTolerance ||
     barycentric_u + barycentric_v > 1.0 + kRayTolerance) {
    return false;
  }

  const double distance = dot_product(edge02, q_vector) * inverse_determinant;
  if(!std::isfinite(distance) || distance <= kRayTolerance ||
     distance > max_distance + kRayTolerance) {
    return false;
  }

  hit_distance = distance;
  return true;
}

} // namespace

void AutoCfdSurfaceProximityIndex::clear() noexcept
{
  triangles_.clear();
  triangle_indices_.clear();
  nodes_.clear();
  stats_ = {};
}

void AutoCfdSurfaceProximityIndex::build(std::vector<AutoCfdSurfaceProxyTriangle> triangles)
{
  clear();
  if(triangles.empty()) {
    return;
  }

  triangles_ = std::move(triangles);
  triangle_indices_.resize(triangles_.size());
  for(std::uint32_t index = 0U; index < triangle_indices_.size(); ++index) {
    triangle_indices_[index] = index;
  }

  const auto root_index =
    build_node(0U, static_cast<std::uint32_t>(triangle_indices_.size()));
  static_cast<void>(root_index);
  stats_.triangle_count = triangles_.size();
  stats_.node_count = nodes_.size();
}

bool AutoCfdSurfaceProximityIndex::empty() const noexcept
{
  return triangles_.empty() || nodes_.empty();
}

const AutoCfdSurfaceProximityIndexStats &AutoCfdSurfaceProximityIndex::stats() const noexcept
{
  return stats_;
}

const std::vector<AutoCfdSurfaceProxyTriangle> &AutoCfdSurfaceProximityIndex::triangles() const noexcept
{
  return triangles_;
}

bool AutoCfdSurfaceProximityIndex::query_nearest_hit(
  const AutoCfdSurfaceProximityRayRequest &request,
  AutoCfdSurfaceProximityRayHit &hit
) const noexcept
{
  hit = {};
  if(empty() ||
     !std::isfinite(request.max_distance) ||
     request.max_distance <= 0.0 ||
     (request.require_same_part && request.require_different_part) ||
     ((request.require_same_part || request.require_different_part) &&
      request.source_part_index == invalid_auto_cfd_surface_proximity_bucket_index)) {
    return false;
  }

  struct StackEntry final {
    std::uint32_t node_index = 0U;
    double near_distance = 0.0;
  };

  std::array<StackEntry, kMaximumTraversalDepth> stack {};
  std::size_t stack_size = 0U;

  const auto root_interval = ray_box_interval(
    request.origin,
    request.direction,
    request.max_distance,
    nodes_.front().bounds_min,
    nodes_.front().bounds_max
  );
  if(!root_interval.hit) {
    return false;
  }

  stack[stack_size++] = {0U, root_interval.near_distance};
  double best_distance = request.max_distance;

  while(stack_size > 0U) {
    const auto current = stack[--stack_size];
    if(current.node_index >= nodes_.size() ||
       current.near_distance > best_distance + kRayTolerance) {
      continue;
    }

    const auto &node = nodes_[current.node_index];
    if(node.leaf) {
      for(std::uint32_t offset = 0U; offset < node.count; ++offset) {
        const std::uint32_t triangle_slot = node.first + offset;
        if(triangle_slot >= triangle_indices_.size()) {
          continue;
        }

        const auto &triangle = triangles_[triangle_indices_[triangle_slot]];
        if(triangle.owner == request.excluded_owner ||
           triangle.part_index == invalid_auto_cfd_surface_proximity_bucket_index) {
          continue;
        }
        if(request.require_same_part &&
           triangle.part_index != request.source_part_index) {
          continue;
        }
        if(request.require_different_part &&
           triangle.part_index == request.source_part_index) {
          continue;
        }

        // Two proximity patterns need detection:
        //   1. Opposing faces (trailing edge): source and hit normals are
        //      anti-parallel — check angle(source_normal, hit_normal).
        //   2. Fold/crease (leading edge): source and hit normals are roughly
        //      aligned (both point outward) — check angle(direction, hit_normal)
        //      where direction = -source_normal for the backward ray.
        // Accept if EITHER check passes (OR logic).
        const double angle_vs_source =
          angle_degrees_between(request.source_normal, triangle.normal);
        const double angle_vs_direction =
          angle_degrees_between(request.direction, triangle.normal);
        const double normal_angle_degrees =
          std::max(angle_vs_source, angle_vs_direction);
        if(!std::isfinite(normal_angle_degrees) ||
           normal_angle_degrees < request.minimum_normal_angle_degrees) {
          continue;
        }

        double candidate_distance = std::numeric_limits<double>::infinity();
        if(!intersect_ray_with_triangle(
             request.origin,
             request.direction,
             best_distance,
             triangle,
             candidate_distance
           )) {
          continue;
        }

        best_distance = candidate_distance;
        hit.found = true;
        hit.hit_owner = triangle.owner;
        hit.hit_distance = candidate_distance;
        hit.hit_normal_angle_degrees = normal_angle_degrees;
      }
      continue;
    }

    const auto push_child =
      [&](std::uint32_t child_index, double near_distance) {
        if(child_index >= nodes_.size() || stack_size >= stack.size()) {
          return;
        }
        stack[stack_size++] = {child_index, near_distance};
      };

    const auto left_interval = ray_box_interval(
      request.origin,
      request.direction,
      best_distance,
      nodes_[node.left_child].bounds_min,
      nodes_[node.left_child].bounds_max
    );
    const auto right_interval = ray_box_interval(
      request.origin,
      request.direction,
      best_distance,
      nodes_[node.right_child].bounds_min,
      nodes_[node.right_child].bounds_max
    );

    if(left_interval.hit && right_interval.hit) {
      if(left_interval.near_distance <= right_interval.near_distance) {
        push_child(node.right_child, right_interval.near_distance);
        push_child(node.left_child, left_interval.near_distance);
      }
      else {
        push_child(node.left_child, left_interval.near_distance);
        push_child(node.right_child, right_interval.near_distance);
      }
    }
    else if(left_interval.hit) {
      push_child(node.left_child, left_interval.near_distance);
    }
    else if(right_interval.hit) {
      push_child(node.right_child, right_interval.near_distance);
    }
  }

  return hit.found;
}

std::uint32_t AutoCfdSurfaceProximityIndex::build_node(
  std::uint32_t first,
  std::uint32_t count
)
{
  const std::uint32_t node_index = static_cast<std::uint32_t>(nodes_.size());
  nodes_.push_back({});
  nodes_[node_index].first = first;
  nodes_[node_index].count = count;

  if(count == 0U) {
    nodes_[node_index].leaf = true;
    ++stats_.leaf_count;
    return node_index;
  }

  const auto &first_triangle = triangles_[triangle_indices_[first]];
  nodes_[node_index].bounds_min = first_triangle.bounds_min;
  nodes_[node_index].bounds_max = first_triangle.bounds_max;
  geo::Point3 centroid_min = first_triangle.centroid;
  geo::Point3 centroid_max = first_triangle.centroid;

  for(std::uint32_t offset = 1U; offset < count; ++offset) {
    const auto &triangle = triangles_[triangle_indices_[first + offset]];
    for(std::size_t axis = 0U; axis < 3U; ++axis) {
      nodes_[node_index].bounds_min[axis] = std::min(
        nodes_[node_index].bounds_min[axis],
        triangle.bounds_min[axis]
      );
      nodes_[node_index].bounds_max[axis] = std::max(
        nodes_[node_index].bounds_max[axis],
        triangle.bounds_max[axis]
      );
      centroid_min[axis] = std::min(centroid_min[axis], triangle.centroid[axis]);
      centroid_max[axis] = std::max(centroid_max[axis], triangle.centroid[axis]);
    }
  }

  std::size_t split_axis = 0U;
  double split_extent = centroid_max[0] - centroid_min[0];
  for(std::size_t axis = 1U; axis < 3U; ++axis) {
    const double extent = centroid_max[axis] - centroid_min[axis];
    if(extent > split_extent) {
      split_axis = axis;
      split_extent = extent;
    }
  }

  if(count <= kMaximumLeafTriangleCount || split_extent <= kRayTolerance) {
    nodes_[node_index].leaf = true;
    ++stats_.leaf_count;
    stats_.max_leaf_triangle_count = std::max<std::size_t>(
      stats_.max_leaf_triangle_count,
      count
    );
    return node_index;
  }

  const std::uint32_t midpoint = first + count / 2U;
  std::nth_element(
    triangle_indices_.begin() + first,
    triangle_indices_.begin() + midpoint,
    triangle_indices_.begin() + first + count,
    [&](std::uint32_t lhs_index, std::uint32_t rhs_index) {
      return triangles_[lhs_index].centroid[split_axis] <
             triangles_[rhs_index].centroid[split_axis];
    }
  );

  nodes_[node_index].leaf = false;
  nodes_[node_index].count = 0U;
  nodes_[node_index].left_child = build_node(first, midpoint - first);
  nodes_[node_index].right_child = build_node(midpoint, first + count - midpoint);
  return node_index;
}

} // namespace sqmesh::mesh::detail
