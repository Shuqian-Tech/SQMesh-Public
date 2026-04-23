// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#pragma once

#include "sqmesh/geo/api.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace sqmesh::mesh::detail {

inline constexpr std::uint32_t invalid_auto_cfd_surface_proximity_bucket_index =
  std::numeric_limits<std::uint32_t>::max();

struct AutoCfdSurfaceProxyTriangle final {
  geo::TopologyEntityId owner {};
  std::uint32_t shell_index = invalid_auto_cfd_surface_proximity_bucket_index;
  std::uint32_t part_index = invalid_auto_cfd_surface_proximity_bucket_index;
  std::array<geo::Point3, 3> vertices {};
  geo::Point3 bounds_min {0.0, 0.0, 0.0};
  geo::Point3 bounds_max {0.0, 0.0, 0.0};
  geo::Point3 centroid {0.0, 0.0, 0.0};
  geo::Vector3 normal {0.0, 0.0, 0.0};
};

struct AutoCfdSurfaceProximityIndexStats final {
  std::size_t triangle_count = 0U;
  std::size_t node_count = 0U;
  std::size_t leaf_count = 0U;
  std::size_t max_leaf_triangle_count = 0U;
};

struct AutoCfdSurfaceProximityRayRequest final {
  geo::Point3 origin {0.0, 0.0, 0.0};
  geo::Vector3 direction {0.0, 0.0, 0.0};
  geo::Vector3 source_normal {0.0, 0.0, 0.0};
  geo::TopologyEntityId excluded_owner {};
  std::uint32_t source_part_index = invalid_auto_cfd_surface_proximity_bucket_index;
  double max_distance = 0.0;
  double minimum_normal_angle_degrees = 0.0;
  bool require_same_part = false;
  bool require_different_part = false;
};

struct AutoCfdSurfaceProximityRayHit final {
  bool found = false;
  geo::TopologyEntityId hit_owner {};
  double hit_distance = std::numeric_limits<double>::infinity();
  double hit_normal_angle_degrees = 0.0;
};

class AutoCfdSurfaceProximityIndex final {
public:
  void clear() noexcept;

  void build(std::vector<AutoCfdSurfaceProxyTriangle> triangles);

  [[nodiscard]] bool empty() const noexcept;

  [[nodiscard]] const AutoCfdSurfaceProximityIndexStats &stats() const noexcept;

  [[nodiscard]] const std::vector<AutoCfdSurfaceProxyTriangle> &triangles() const noexcept;

  [[nodiscard]] bool query_nearest_hit(
    const AutoCfdSurfaceProximityRayRequest &request,
    AutoCfdSurfaceProximityRayHit &hit
  ) const noexcept;

private:
  struct Node final {
    geo::Point3 bounds_min {0.0, 0.0, 0.0};
    geo::Point3 bounds_max {0.0, 0.0, 0.0};
    std::uint32_t first = 0U;
    std::uint32_t count = 0U;
    std::uint32_t left_child = invalid_auto_cfd_surface_proximity_bucket_index;
    std::uint32_t right_child = invalid_auto_cfd_surface_proximity_bucket_index;
    bool leaf = false;
  };

  [[nodiscard]] std::uint32_t build_node(std::uint32_t first, std::uint32_t count);

  std::vector<AutoCfdSurfaceProxyTriangle> triangles_ {};
  std::vector<std::uint32_t> triangle_indices_ {};
  std::vector<Node> nodes_ {};
  AutoCfdSurfaceProximityIndexStats stats_ {};
};

} // namespace sqmesh::mesh::detail
