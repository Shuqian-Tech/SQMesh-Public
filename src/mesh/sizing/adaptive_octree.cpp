// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "adaptive_octree.hpp"

#include "core/log.hpp"

#include <algorithm>
#include <cmath>
#include <queue>
#include <utility>

namespace sqmesh::mesh::detail {
namespace {

constexpr double kOctreePadding = 1.0e-6;

[[nodiscard]] double point_distance(
  const geo::Point3 &a,
  const geo::Point3 &b
) noexcept
{
  const double dx = a[0] - b[0];
  const double dy = a[1] - b[1];
  const double dz = a[2] - b[2];
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

/// The axis (0=X,1=Y,2=Z) and sign direction for each face index.
/// Face 0=-X, 1=+X, 2=-Y, 3=+Y, 4=-Z, 5=+Z.
[[nodiscard]] constexpr std::uint8_t face_axis(std::uint8_t face) noexcept
{
  return face / 2U;
}

[[nodiscard]] constexpr bool face_positive(std::uint8_t face) noexcept
{
  return (face & 1U) != 0U;
}

/// The opposite face direction.
[[nodiscard]] constexpr std::uint8_t opposite_face(std::uint8_t face) noexcept
{
  return face ^ 1U;
}

/// Check if octant is on the positive side of the given axis.
/// Morton order: bit 0 = X, bit 1 = Y, bit 2 = Z.
[[nodiscard]] constexpr bool octant_on_positive_side(
  std::uint8_t octant,
  std::uint8_t axis
) noexcept
{
  return (octant & (1U << axis)) != 0U;
}

/// Flip the axis bit in an octant index (mirror across axis).
[[nodiscard]] constexpr std::uint8_t octant_flip_axis(
  std::uint8_t octant,
  std::uint8_t axis
) noexcept
{
  return octant ^ static_cast<std::uint8_t>(1U << axis);
}

/// Which child octant index does node_index occupy in its parent?
[[nodiscard]] std::uint8_t which_child(
  const std::vector<AdaptiveOctreeNode> &nodes,
  std::uint32_t node_index
) noexcept
{
  const auto parent_index = nodes[node_index].parent;
  if(parent_index == invalid_octree_node_index) {
    return 0U;
  }
  const auto &parent = nodes[parent_index];
  for(std::uint8_t i = 0U; i < 8U; ++i) {
    if(parent.children[i] == node_index) {
      return i;
    }
  }
  return 0U;
}

// ---------------------------------------------------------------------------
// 28 vertex pair table for octant diffusion.
// All C(8,2) = 28 combinations of 8 corner vertices.
// ---------------------------------------------------------------------------

constexpr int kVertexPairCount = 28;

constexpr int kVertexPairs[28][2] = {
  // 12 edges (distance = 1.0 * side_length)
  {0, 1}, {2, 3}, {4, 5}, {6, 7},   // X-aligned
  {0, 2}, {1, 3}, {4, 6}, {5, 7},   // Y-aligned
  {0, 4}, {1, 5}, {2, 6}, {3, 7},   // Z-aligned
  // 12 face diagonals (distance = sqrt(2) * side_length)
  {0, 6}, {2, 4}, {1, 7}, {3, 5},
  {0, 5}, {1, 4}, {2, 7}, {3, 6},
  {0, 3}, {1, 2}, {4, 7}, {5, 6},
  // 4 body diagonals (distance = sqrt(3) * side_length)
  {0, 7}, {1, 6}, {3, 4}, {2, 5},
};

constexpr double kVertexPairDist[28] = {
  1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
  1.4142135623730951, 1.4142135623730951, 1.4142135623730951, 1.4142135623730951,
  1.4142135623730951, 1.4142135623730951, 1.4142135623730951, 1.4142135623730951,
  1.4142135623730951, 1.4142135623730951, 1.4142135623730951, 1.4142135623730951,
  1.7320508075688772, 1.7320508075688772, 1.7320508075688772, 1.7320508075688772,
};

} // namespace

// ---------------------------------------------------------------------------
// clear / initialize
// ---------------------------------------------------------------------------

void AdaptiveOctree::clear() noexcept
{
  nodes_.clear();
  vertices_.clear();
  vertex_map_.clear();
  stats_ = {};
  built_ = false;
}

void AdaptiveOctree::initialize(
  const geo::Point3 &domain_min,
  const geo::Point3 &domain_max,
  const AdaptiveOctreeConfig &config
)
{
  clear();
  config_ = config;

  // Build a cubic root that encloses the domain.
  geo::Point3 center {0.0, 0.0, 0.0};
  double max_half = 0.0;
  for(std::size_t axis = 0U; axis < 3U; ++axis) {
    center[axis] = 0.5 * (domain_min[axis] + domain_max[axis]);
    const double half = 0.5 * (domain_max[axis] - domain_min[axis]) + kOctreePadding;
    max_half = std::max(max_half, half);
  }
  max_half = std::max(max_half, kOctreePadding);

  nodes_.push_back({});
  nodes_[0].center = center;
  nodes_[0].half_extent = max_half;
  nodes_[0].depth = 0;
  nodes_[0].leaf = true;
  nodes_[0].parent = invalid_octree_node_index;

  stats_.node_count = 1U;
  stats_.leaf_count = 1U;
  built_ = true;
}

// ---------------------------------------------------------------------------
// child_octant / child_center
// ---------------------------------------------------------------------------

std::uint8_t AdaptiveOctree::child_octant(
  const geo::Point3 &center,
  const geo::Point3 &position
) noexcept
{
  std::uint8_t octant = 0U;
  if(position[0] >= center[0]) { octant |= 1U; }
  if(position[1] >= center[1]) { octant |= 2U; }
  if(position[2] >= center[2]) { octant |= 4U; }
  return octant;
}

geo::Point3 AdaptiveOctree::child_center(
  const geo::Point3 &parent_center,
  double parent_half_extent,
  std::uint8_t octant
) noexcept
{
  const double quarter = parent_half_extent * 0.5;
  return {
    parent_center[0] + ((octant & 1U) ? quarter : -quarter),
    parent_center[1] + ((octant & 2U) ? quarter : -quarter),
    parent_center[2] + ((octant & 4U) ? quarter : -quarter),
  };
}

// ---------------------------------------------------------------------------
// subdivide
// ---------------------------------------------------------------------------

void AdaptiveOctree::subdivide(std::uint32_t node_index)
{
  if(!nodes_[node_index].leaf) {
    return;
  }

  // Cache parent properties before push_back can invalidate references.
  const geo::Point3 parent_center = nodes_[node_index].center;
  const double parent_half = nodes_[node_index].half_extent;
  const auto child_depth = static_cast<std::uint8_t>(nodes_[node_index].depth + 1U);
  const double child_half = parent_half * 0.5;

  nodes_[node_index].leaf = false;
  --stats_.leaf_count;

  // Reserve to prevent reallocation during the loop.
  nodes_.reserve(nodes_.size() + 8U);

  for(std::uint8_t octant = 0U; octant < 8U; ++octant) {
    const auto child_idx = static_cast<std::uint32_t>(nodes_.size());

    AdaptiveOctreeNode child_node;
    child_node.center = child_center(parent_center, parent_half, octant);
    child_node.half_extent = child_half;
    child_node.size_value = std::numeric_limits<double>::infinity();
    child_node.parent = node_index;
    child_node.depth = child_depth;
    child_node.leaf = true;

    nodes_.push_back(child_node);
    nodes_[node_index].children[octant] = child_idx;

    ++stats_.leaf_count;
    ++stats_.node_count;
    if(child_depth > stats_.max_depth_reached) {
      stats_.max_depth_reached = child_depth;
    }
  }
}

// ---------------------------------------------------------------------------
// find_leaf
// ---------------------------------------------------------------------------

std::uint32_t AdaptiveOctree::find_leaf(const geo::Point3 &position) const noexcept
{
  if(nodes_.empty()) {
    return invalid_octree_node_index;
  }

  std::uint32_t current = 0U;
  while(!nodes_[current].leaf) {
    const auto octant = child_octant(nodes_[current].center, position);
    const auto child = nodes_[current].children[octant];
    if(child == invalid_octree_node_index || child >= nodes_.size()) {
      return current;
    }
    current = child;
  }
  return current;
}

// ---------------------------------------------------------------------------
// insert_size_seed
// ---------------------------------------------------------------------------

void AdaptiveOctree::insert_size_seed(const geo::Point3 &position, double size_value)
{
  if(nodes_.empty() || !std::isfinite(size_value) || size_value <= 0.0) {
    return;
  }

  ++stats_.seed_count;

  std::uint32_t leaf = find_leaf(position);
  if(leaf == invalid_octree_node_index) {
    return;
  }

  // Subdivide while the leaf is much coarser than the seed size.
  while(nodes_[leaf].half_extent > 2.0 * size_value &&
        nodes_[leaf].depth < config_.max_depth &&
        nodes_[leaf].half_extent > config_.minimum_node_extent) {
    subdivide(leaf);
    leaf = find_leaf(position);
    if(leaf == invalid_octree_node_index) {
      return;
    }
  }

  // Only set size_value on the FINAL leaf containing the seed point.
  nodes_[leaf].size_value = std::min(nodes_[leaf].size_value, size_value);
}

// ---------------------------------------------------------------------------
// find_face_neighbor (Morton-code-based tree traversal)
// ---------------------------------------------------------------------------

std::uint32_t AdaptiveOctree::find_face_neighbor(
  std::uint32_t node_index,
  std::uint8_t face
) const noexcept
{
  if(node_index >= nodes_.size()) {
    return invalid_octree_node_index;
  }

  const std::uint8_t axis = face_axis(face);
  const bool positive = face_positive(face);

  const std::uint32_t parent = nodes_[node_index].parent;
  if(parent == invalid_octree_node_index) {
    return invalid_octree_node_index;
  }

  const std::uint8_t octant = which_child(nodes_, node_index);
  const bool on_boundary = (octant_on_positive_side(octant, axis) == positive);

  if(!on_boundary) {
    // Sibling across this axis within the same parent.
    const std::uint8_t sibling_octant = octant_flip_axis(octant, axis);
    const auto sibling = nodes_[parent].children[sibling_octant];
    if(sibling == invalid_octree_node_index) {
      return invalid_octree_node_index;
    }
    if(nodes_[sibling].leaf) {
      return sibling;
    }
    const auto child = nodes_[sibling].children[octant];
    return (child != invalid_octree_node_index) ? child : sibling;
  }

  // On boundary: recurse up to find the parent's neighbor.
  const auto parent_neighbor = find_face_neighbor(parent, face);
  if(parent_neighbor == invalid_octree_node_index) {
    return invalid_octree_node_index;
  }

  if(nodes_[parent_neighbor].leaf) {
    return parent_neighbor;
  }

  // Descend into the parent neighbor's child that mirrors our octant.
  const std::uint8_t mirror_octant = octant_flip_axis(octant, axis);
  const auto child = nodes_[parent_neighbor].children[mirror_octant];
  return (child != invalid_octree_node_index) ? child : parent_neighbor;
}

// ---------------------------------------------------------------------------
// collect_leaf_neighbors_across_face
// ---------------------------------------------------------------------------

void AdaptiveOctree::collect_leaf_neighbors_across_face(
  std::uint32_t node_index,
  std::uint8_t incoming_face,
  std::vector<std::uint32_t> &result
) const
{
  if(node_index == invalid_octree_node_index || node_index >= nodes_.size()) {
    return;
  }

  if(nodes_[node_index].leaf) {
    result.push_back(node_index);
    return;
  }

  const std::uint8_t axis = face_axis(incoming_face);
  const bool want_positive = face_positive(incoming_face);

  for(std::uint8_t octant = 0U; octant < 8U; ++octant) {
    if(octant_on_positive_side(octant, axis) != want_positive) {
      continue;
    }
    const auto child = nodes_[node_index].children[octant];
    if(child == invalid_octree_node_index) {
      continue;
    }
    collect_leaf_neighbors_across_face(child, incoming_face, result);
  }
}

// ---------------------------------------------------------------------------
// refine_for_gradation
// ---------------------------------------------------------------------------

constexpr std::size_t kMaxOctreeNodes = 500000U;
constexpr std::size_t kMaxRefinementPasses = 8U;

void AdaptiveOctree::refine_for_gradation(double growth_rate)
{
  if(nodes_.empty() || growth_rate <= 1.0) {
    return;
  }

  const double alpha = growth_rate - 1.0;

  for(std::size_t pass = 0U; pass < kMaxRefinementPasses; ++pass) {
    if(nodes_.size() >= kMaxOctreeNodes) {
      break;
    }

    bool changed = false;
    const auto count = static_cast<std::uint32_t>(nodes_.size());
    for(std::uint32_t i = 0U; i < count; ++i) {
      if(!nodes_[i].leaf || !std::isfinite(nodes_[i].size_value)) {
        continue;
      }

      for(std::uint8_t face = 0U; face < 6U; ++face) {
        const auto nb = find_face_neighbor(i, face);
        if(nb == invalid_octree_node_index || nb >= nodes_.size()) {
          continue;
        }
        if(!nodes_[nb].leaf || !std::isfinite(nodes_[nb].size_value)) {
          continue;
        }

        const double dist = point_distance(nodes_[i].center, nodes_[nb].center);
        const double allowed_diff = alpha * dist;
        const double actual_diff = std::abs(
          nodes_[i].size_value - nodes_[nb].size_value
        );

        if(actual_diff > allowed_diff * 2.0) {
          const auto to_split =
            (nodes_[i].half_extent >= nodes_[nb].half_extent) ? i : nb;
          if(nodes_[to_split].leaf &&
             nodes_[to_split].depth < config_.max_depth &&
             nodes_[to_split].half_extent > config_.minimum_node_extent &&
             nodes_.size() < kMaxOctreeNodes) {
            subdivide(to_split);
            ++stats_.refinement_splits;
            changed = true;
          }
        }
      }
    }

    if(!changed) {
      break;
    }
  }
}

// ---------------------------------------------------------------------------
// balance_2_to_1
// ---------------------------------------------------------------------------

void AdaptiveOctree::balance_2_to_1()
{
  if(nodes_.empty()) {
    return;
  }

  std::vector<std::uint32_t> work_queue;
  work_queue.reserve(nodes_.size());

  for(std::uint32_t i = 0U; i < nodes_.size(); ++i) {
    if(nodes_[i].leaf) {
      work_queue.push_back(i);
    }
  }

  constexpr std::size_t kMaxBalanceNodes = 500000U;
  std::size_t head = 0U;

  while(head < work_queue.size() && nodes_.size() < kMaxBalanceNodes) {
    const auto idx = work_queue[head++];
    if(idx >= nodes_.size() || !nodes_[idx].leaf) {
      continue;
    }

    const auto my_depth = nodes_[idx].depth;

    for(std::uint8_t face = 0U; face < 6U; ++face) {
      const auto nb = find_face_neighbor(idx, face);
      if(nb == invalid_octree_node_index || nb >= nodes_.size()) {
        continue;
      }
      if(!nodes_[nb].leaf) {
        continue;
      }

      if(nodes_[nb].depth + 1 < my_depth &&
         nodes_[nb].depth < config_.max_depth &&
         nodes_[nb].half_extent > config_.minimum_node_extent &&
         nodes_.size() < kMaxBalanceNodes) {
        subdivide(nb);
        ++stats_.balance_splits;

        for(std::uint8_t c = 0U; c < 8U; ++c) {
          const auto child = nodes_[nb].children[c];
          if(child != invalid_octree_node_index) {
            work_queue.push_back(child);
          }
        }
      }
    }
  }
}

// ---------------------------------------------------------------------------
// smooth_gradients (Dijkstra)
// ---------------------------------------------------------------------------

void AdaptiveOctree::smooth_gradients(
  double growth_rate,
  double minimum_size,
  double maximum_size
)
{
  growth_rate_ = growth_rate;

  if(nodes_.empty() || growth_rate <= 1.0) {
    return;
  }

  const double alpha = growth_rate - 1.0;

  struct QueueEntry final {
    double size = std::numeric_limits<double>::infinity();
    std::uint32_t index = 0U;
    [[nodiscard]] bool operator>(const QueueEntry &rhs) const noexcept
    {
      return size > rhs.size;
    }
  };

  std::priority_queue<QueueEntry, std::vector<QueueEntry>, std::greater<QueueEntry>> queue;

  for(std::uint32_t i = 0U; i < nodes_.size(); ++i) {
    if(!nodes_[i].leaf) {
      continue;
    }
    if(!std::isfinite(nodes_[i].size_value) || nodes_[i].size_value <= 0.0) {
      nodes_[i].size_value = maximum_size;
    }
    queue.push({nodes_[i].size_value, i});
  }

  std::vector<std::uint32_t> neighbor_buffer;

  while(!queue.empty()) {
    const auto current = queue.top();
    queue.pop();

    if(current.index >= nodes_.size() || !nodes_[current.index].leaf) {
      continue;
    }
    if(current.size > nodes_[current.index].size_value + 1.0e-12) {
      continue;
    }

    for(std::uint8_t face = 0U; face < 6U; ++face) {
      const auto neighbor_index = find_face_neighbor(current.index, face);
      if(neighbor_index == invalid_octree_node_index) {
        continue;
      }

      neighbor_buffer.clear();
      if(nodes_[neighbor_index].leaf) {
        neighbor_buffer.push_back(neighbor_index);
      } else {
        collect_leaf_neighbors_across_face(
          neighbor_index, opposite_face(face), neighbor_buffer
        );
      }

      for(const auto nb : neighbor_buffer) {
        if(nb >= nodes_.size() || !nodes_[nb].leaf) {
          continue;
        }

        const double dist = point_distance(
          nodes_[current.index].center, nodes_[nb].center
        );
        const double limited = nodes_[current.index].size_value + alpha * dist;
        const double clamped = std::max(minimum_size, std::min(maximum_size, limited));

        if(clamped < nodes_[nb].size_value - 1.0e-12) {
          nodes_[nb].size_value = clamped;
          queue.push({clamped, nb});
          ++stats_.smoothing_updates;
        }
      }
    }
  }
}

// ---------------------------------------------------------------------------
// Vertex management (integer-coordinate hashing)
// ---------------------------------------------------------------------------

std::uint64_t AdaptiveOctree::pack_vertex_key(int ix, int iy, int iz) noexcept
{
  // max_depth <= 20, so coordinates fit in 21 bits each (0..2^20 = 1048576).
  const auto ux = static_cast<std::uint64_t>(ix) & 0x1FFFFFULL;
  const auto uy = static_cast<std::uint64_t>(iy) & 0x1FFFFFULL;
  const auto uz = static_cast<std::uint64_t>(iz) & 0x1FFFFFULL;
  return ux | (uy << 21ULL) | (uz << 42ULL);
}

void AdaptiveOctree::compute_node_min_icoord(
  std::uint32_t node_index,
  std::array<int, 3> &coord
) const
{
  // Build path from this node to root.
  // max_depth <= 20 so the path is short.
  std::uint8_t path[24];
  int path_len = 0;

  std::uint32_t cur = node_index;
  while(cur != 0U && nodes_[cur].parent != invalid_octree_node_index) {
    path[path_len++] = which_child(nodes_, cur);
    cur = nodes_[cur].parent;
    if(path_len >= 24) {
      break;
    }
  }

  // Replay from root to node.
  coord = {0, 0, 0};
  int span = 1 << config_.max_depth;

  for(int i = path_len - 1; i >= 0; --i) {
    span >>= 1;
    const std::uint8_t octant = path[i];
    if(octant & 1U) { coord[0] += span; }
    if(octant & 2U) { coord[1] += span; }
    if(octant & 4U) { coord[2] += span; }
  }
}

std::uint32_t AdaptiveOctree::get_or_create_vertex(int ix, int iy, int iz)
{
  const auto key = pack_vertex_key(ix, iy, iz);
  const auto it = vertex_map_.find(key);
  if(it != vertex_map_.end()) {
    return it->second;
  }

  const auto idx = static_cast<std::uint32_t>(vertices_.size());
  vertices_.push_back({});
  vertex_map_[key] = idx;
  return idx;
}

// ---------------------------------------------------------------------------
// compute_leaf_gradients  (vertex creation + 28-pair diffusion)
// ---------------------------------------------------------------------------

void AdaptiveOctree::compute_leaf_gradients()
{
  if(nodes_.empty() || growth_rate_ <= 1.0) {
    return;
  }

  const double alpha = growth_rate_ - 1.0;

  // ── Phase 1: Create shared vertices for all leaf corners ──────────────

  vertices_.clear();
  vertex_map_.clear();
  vertices_.reserve(stats_.leaf_count * 2U); // rough estimate

  for(std::uint32_t i = 0U; i < nodes_.size(); ++i) {
    if(!nodes_[i].leaf) {
      continue;
    }

    std::array<int, 3> min_coord;
    compute_node_min_icoord(i, min_coord);
    const int span = 1 << (config_.max_depth - nodes_[i].depth);

    for(std::uint8_t corner = 0U; corner < 8U; ++corner) {
      const int ix = min_coord[0] + ((corner & 1U) ? span : 0);
      const int iy = min_coord[1] + ((corner & 2U) ? span : 0);
      const int iz = min_coord[2] + ((corner & 4U) ? span : 0);
      nodes_[i].vertex_indices[corner] = get_or_create_vertex(ix, iy, iz);
    }
  }

  stats_.vertex_count = vertices_.size();

  // ── Phase 2: Initialize vertex sizes from Dijkstra-smoothed cell centers ─

  // Each vertex gets the minimum size_value of all leaf cells that share it.
  // This is conservative and ensures small-size cells dominate at shared corners.
  for(std::uint32_t i = 0U; i < nodes_.size(); ++i) {
    if(!nodes_[i].leaf) {
      continue;
    }
    const double cell_size = nodes_[i].size_value;
    for(std::uint8_t corner = 0U; corner < 8U; ++corner) {
      const auto vi = nodes_[i].vertex_indices[corner];
      if(vi < vertices_.size()) {
        vertices_[vi].mesh_size = std::min(vertices_[vi].mesh_size, cell_size);
      }
    }
  }

  // ── Phase 3: 28-pair diffusion ────────────────────────────────────────
  //
  // For each leaf octant, enforce growth rate between all 28 vertex pairs.
  // Repeat until convergence (max 10 passes).

  constexpr int kMaxDiffusionPasses = 10;

  for(int pass = 0; pass < kMaxDiffusionPasses; ++pass) {
    bool changed = false;

    for(std::uint32_t i = 0U; i < nodes_.size(); ++i) {
      if(!nodes_[i].leaf) {
        continue;
      }

      const double side_length = nodes_[i].half_extent * 2.0;

      // Read current vertex sizes for this octant.
      double msize[8];
      for(int n = 0; n < 8; ++n) {
        const auto vi = nodes_[i].vertex_indices[n];
        msize[n] = (vi < vertices_.size()) ? vertices_[vi].mesh_size
                                           : nodes_[i].size_value;
      }

      // Enforce growth rate on all 28 pairs.
      for(int p = 0; p < kVertexPairCount; ++p) {
        const int a = kVertexPairs[p][0];
        const int b = kVertexPairs[p][1];
        const double dist = side_length * kVertexPairDist[p];

        if(msize[b] <= msize[a]) {
          const double limit = msize[b] + alpha * dist;
          if(limit < msize[a]) {
            msize[a] = limit;
          }
        } else {
          const double limit = msize[a] + alpha * dist;
          if(limit < msize[b]) {
            msize[b] = limit;
          }
        }
      }

      // Write back.  Only update if the value decreased.
      for(int n = 0; n < 8; ++n) {
        const auto vi = nodes_[i].vertex_indices[n];
        if(vi < vertices_.size() && msize[n] < vertices_[vi].mesh_size - 1.0e-12) {
          vertices_[vi].mesh_size = msize[n];
          changed = true;
        }
      }
    }

    ++stats_.diffusion_passes;

    if(!changed) {
      break;
    }
  }

}

// ---------------------------------------------------------------------------
// query (trilinear interpolation from shared corner vertices)
// ---------------------------------------------------------------------------

double AdaptiveOctree::query(const geo::Point3 &position) const noexcept
{
  if(nodes_.empty()) {
    return std::numeric_limits<double>::infinity();
  }

  const auto leaf = find_leaf(position);
  if(leaf == invalid_octree_node_index) {
    return std::numeric_limits<double>::infinity();
  }

  const auto &node = nodes_[leaf];

  // Check that vertices are assigned (compute_leaf_gradients was called).
  if(node.vertex_indices[0] == invalid_octree_vertex_index) {
    return node.size_value;
  }

  // Normalized coordinates r, s, t in [-1, 1] within the cell.
  const double h = node.half_extent;
  if(h <= 0.0) {
    return node.size_value;
  }

  double r = (position[0] - node.center[0]) / h;
  double s = (position[1] - node.center[1]) / h;
  double t = (position[2] - node.center[2]) / h;

  // Clamp to cell bounds.
  r = std::max(-1.0, std::min(1.0, r));
  s = std::max(-1.0, std::min(1.0, s));
  t = std::max(-1.0, std::min(1.0, t));

  // Trilinear interpolation weights.
  // Corner ordering: bit 0 = X, bit 1 = Y, bit 2 = Z (Morton).
  const double rp = (r + 1.0) * 0.5;
  const double rn = 1.0 - rp;
  const double sp = (s + 1.0) * 0.5;
  const double sn = 1.0 - sp;
  const double tp = (t + 1.0) * 0.5;
  const double tn = 1.0 - tp;

  const double psi[8] = {
    rn * sn * tn,   // corner 0: -X -Y -Z
    rp * sn * tn,   // corner 1: +X -Y -Z
    rn * sp * tn,   // corner 2: -X +Y -Z
    rp * sp * tn,   // corner 3: +X +Y -Z
    rn * sn * tp,   // corner 4: -X -Y +Z
    rp * sn * tp,   // corner 5: +X -Y +Z
    rn * sp * tp,   // corner 6: -X +Y +Z
    rp * sp * tp,   // corner 7: +X +Y +Z
  };

  double result = 0.0;
  for(int n = 0; n < 8; ++n) {
    const auto vi = node.vertex_indices[n];
    const double vtx_size = (vi < vertices_.size())
                              ? vertices_[vi].mesh_size
                              : node.size_value;
    result += psi[n] * vtx_size;
  }

  if(!std::isfinite(result) || result <= 0.0) {
    return node.size_value;
  }
  return result;
}

// ---------------------------------------------------------------------------
// accessors
// ---------------------------------------------------------------------------

bool AdaptiveOctree::empty() const noexcept
{
  return nodes_.empty();
}

bool AdaptiveOctree::built() const noexcept
{
  return built_ && !nodes_.empty();
}

const AdaptiveOctreeStats &AdaptiveOctree::stats() const noexcept
{
  return stats_;
}

void AdaptiveOctree::get_vertex_positions(
  std::vector<geo::Point3> &positions
) const
{
  positions.assign(vertices_.size(), {0.0, 0.0, 0.0});
  std::vector<bool> assigned(vertices_.size(), false);

  for(std::uint32_t ni = 0U; ni < nodes_.size(); ++ni) {
    if(!nodes_[ni].leaf) continue;
    const auto &node = nodes_[ni];
    const double h = node.half_extent;

    for(std::uint8_t corner = 0U; corner < 8U; ++corner) {
      const auto vi = node.vertex_indices[corner];
      if(vi >= vertices_.size() || assigned[vi]) continue;

      positions[vi] = {
        node.center[0] + ((corner & 1U) ? h : -h),
        node.center[1] + ((corner & 2U) ? h : -h),
        node.center[2] + ((corner & 4U) ? h : -h),
      };
      assigned[vi] = true;
    }
  }
}

} // namespace sqmesh::mesh::detail
