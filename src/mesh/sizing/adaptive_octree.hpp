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
#include <unordered_map>
#include <vector>

namespace sqmesh::mesh::detail {

inline constexpr std::uint32_t invalid_octree_node_index =
  std::numeric_limits<std::uint32_t>::max();

inline constexpr std::uint32_t invalid_octree_vertex_index =
  std::numeric_limits<std::uint32_t>::max();

/// Shared vertex at an octree corner.  Multiple leaf nodes reference the same
/// vertex via integer-coordinate hashing, guaranteeing C0 continuity across
/// cell boundaries during trilinear interpolation.
struct AdaptiveOctreeVertex final {
  double mesh_size = std::numeric_limits<double>::infinity();
};

struct AdaptiveOctreeNode final {
  geo::Point3 center {0.0, 0.0, 0.0};
  double half_extent = 0.0;
  double size_value = std::numeric_limits<double>::infinity();

  std::uint32_t parent = invalid_octree_node_index;
  std::array<std::uint32_t, 8> children {
    invalid_octree_node_index, invalid_octree_node_index,
    invalid_octree_node_index, invalid_octree_node_index,
    invalid_octree_node_index, invalid_octree_node_index,
    invalid_octree_node_index, invalid_octree_node_index,
  };

  /// Indices into the shared vertex array.  One per corner (Morton order).
  std::array<std::uint32_t, 8> vertex_indices {
    invalid_octree_vertex_index, invalid_octree_vertex_index,
    invalid_octree_vertex_index, invalid_octree_vertex_index,
    invalid_octree_vertex_index, invalid_octree_vertex_index,
    invalid_octree_vertex_index, invalid_octree_vertex_index,
  };

  std::uint8_t depth = 0;
  bool leaf = true;
};

struct AdaptiveOctreeStats final {
  std::size_t node_count = 0U;
  std::size_t leaf_count = 0U;
  std::size_t max_depth_reached = 0U;
  std::size_t seed_count = 0U;
  std::size_t refinement_splits = 0U;
  std::size_t balance_splits = 0U;
  std::size_t smoothing_updates = 0U;
  std::size_t vertex_count = 0U;
  std::size_t diffusion_passes = 0U;
};

struct AdaptiveOctreeConfig final {
  std::uint8_t max_depth = 14;
  double minimum_node_extent = 1.0e-9;
};

class AdaptiveOctree final {
public:
  void clear() noexcept;

  /// Initialize the root node to a cubic domain enclosing [domain_min, domain_max].
  void initialize(
    const geo::Point3 &domain_min,
    const geo::Point3 &domain_max,
    const AdaptiveOctreeConfig &config
  );

  /// Insert a sizing seed. The containing leaf's size_value is reduced to
  /// min(current, size_value), and the leaf is subdivided if its extent
  /// is much larger than the seed size.
  void insert_size_seed(const geo::Point3 &position, double size_value);

  /// Adaptive refinement: split leaf nodes where adjacent leaves have
  /// size values that differ by more than growth_slope * separation.
  void refine_for_gradation(double growth_rate);

  /// Enforce 2:1 balance: no two adjacent leaves may differ by more than
  /// one level of depth.  This creates the intermediate cells needed for
  /// smooth gradient propagation.
  void balance_2_to_1();

  /// Dijkstra-based gradient smoothing on the octree adjacency graph.
  /// Enforces |h(A) - h(B)| <= (growth_rate - 1) * dist(A, B) for all
  /// adjacent leaf pairs.
  void smooth_gradients(
    double growth_rate,
    double minimum_size,
    double maximum_size
  );

  /// Create shared corner vertices for all leaf nodes and propagate
  /// size values via 28-pair vertex diffusion.  Uses vertex hashing
  /// so adjacent cells share the same vertex objects, guaranteeing C0
  /// continuity for trilinear queries.  Must be called after smooth_gradients().
  void compute_leaf_gradients();

  /// Query the interpolated size value at an arbitrary point.
  /// Uses trilinear interpolation from the 8 shared corner vertices.
  [[nodiscard]] double query(const geo::Point3 &position) const noexcept;

  [[nodiscard]] bool empty() const noexcept;
  [[nodiscard]] bool built() const noexcept;
  [[nodiscard]] const AdaptiveOctreeStats &stats() const noexcept;

  /// Mutable access to vertices for geometric source size updates.
  [[nodiscard]] std::vector<AdaptiveOctreeVertex> &mutable_vertices() noexcept
  {
    return vertices_;
  }

  /// Get 3D positions of all vertices (computed from integer coordinates).
  void get_vertex_positions(std::vector<geo::Point3> &positions) const;

private:
  /// Find the leaf node containing a point. Returns node index.
  [[nodiscard]] std::uint32_t find_leaf(const geo::Point3 &position) const noexcept;

  /// Subdivide a leaf into 8 children. Children inherit the parent's size_value.
  void subdivide(std::uint32_t node_index);

  /// Determine which child octant (0..7) a point falls in relative to a node center.
  [[nodiscard]] static std::uint8_t child_octant(
    const geo::Point3 &center,
    const geo::Point3 &position
  ) noexcept;

  /// Compute child center from parent center, half_extent, and octant index.
  [[nodiscard]] static geo::Point3 child_center(
    const geo::Point3 &parent_center,
    double parent_half_extent,
    std::uint8_t octant
  ) noexcept;

  /// Find the face neighbor of a node across the given face direction (0..5).
  /// Faces: 0=-X, 1=+X, 2=-Y, 3=+Y, 4=-Z, 5=+Z
  /// Returns invalid_octree_node_index if at domain boundary.
  [[nodiscard]] std::uint32_t find_face_neighbor(
    std::uint32_t node_index,
    std::uint8_t face
  ) const noexcept;

  /// Collect all leaf descendants of a node that touch a given incoming face.
  void collect_leaf_neighbors_across_face(
    std::uint32_t node_index,
    std::uint8_t incoming_face,
    std::vector<std::uint32_t> &result
  ) const;

  /// Compute the min-corner integer coordinate for a node by tracing from root.
  void compute_node_min_icoord(
    std::uint32_t node_index,
    std::array<int, 3> &coord
  ) const;

  /// Pack integer coordinates into a hash key.
  [[nodiscard]] static std::uint64_t pack_vertex_key(int ix, int iy, int iz) noexcept;

  /// Look up or create a shared vertex at the given integer coordinates.
  [[nodiscard]] std::uint32_t get_or_create_vertex(int ix, int iy, int iz);

  std::vector<AdaptiveOctreeNode> nodes_ {};
  std::vector<AdaptiveOctreeVertex> vertices_ {};
  std::unordered_map<std::uint64_t, std::uint32_t> vertex_map_ {};
  AdaptiveOctreeConfig config_ {};
  AdaptiveOctreeStats stats_ {};
  double growth_rate_ = 1.2;
  bool built_ = false;
};

} // namespace sqmesh::mesh::detail
