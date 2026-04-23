// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "size_function.hpp"

#include "core/log.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace sqmesh::mesh::detail {
namespace {

constexpr double kEpsilon = 1.0e-30;

[[nodiscard]] double point_distance(
  const geo::Point3 &a, const geo::Point3 &b
) noexcept
{
  const double dx = b[0] - a[0];
  const double dy = b[1] - a[1];
  const double dz = b[2] - a[2];
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

[[nodiscard]] double dot3(
  const geo::Vector3 &a, const geo::Vector3 &b
) noexcept
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// Compute closest point on edge (p0, p1) to query point, return distance
// and interpolated size.
void closest_point_on_edge(
  const geo::Point3 &p0, const geo::Point3 &p1,
  double size0, double size1,
  const geo::Point3 &query,
  double &dist_out, double &size_out
) noexcept
{
  const double dx = p1[0] - p0[0];
  const double dy = p1[1] - p0[1];
  const double dz = p1[2] - p0[2];
  const double len_sq = dx * dx + dy * dy + dz * dz;

  if(len_sq <= kEpsilon) {
    dist_out = point_distance(p0, query);
    size_out = size0;
    return;
  }

  double t = ((query[0] - p0[0]) * dx +
              (query[1] - p0[1]) * dy +
              (query[2] - p0[2]) * dz) / len_sq;
  t = std::max(0.0, std::min(1.0, t));

  const geo::Point3 closest = {
    p0[0] + t * dx,
    p0[1] + t * dy,
    p0[2] + t * dz,
  };
  dist_out = point_distance(closest, query);
  size_out = size0 + t * (size1 - size0);
}

// Compute closest point on triangle (p0, p1, p2) to query point,
// return distance and barycentric-interpolated size.
void closest_point_on_triangle(
  const geo::Point3 &p0, const geo::Point3 &p1, const geo::Point3 &p2,
  double size0, double size1, double size2,
  const geo::Vector3 &normal, double area,
  const geo::Point3 &query,
  double &dist_out, double &size_out
) noexcept
{
  if(area <= kEpsilon) {
    // Degenerate triangle — fall back to closest vertex
    const double d0 = point_distance(p0, query);
    const double d1 = point_distance(p1, query);
    const double d2 = point_distance(p2, query);
    if(d0 <= d1 && d0 <= d2) { dist_out = d0; size_out = size0; }
    else if(d1 <= d2) { dist_out = d1; size_out = size1; }
    else { dist_out = d2; size_out = size2; }
    return;
  }

  // Project query point onto triangle plane
  const geo::Vector3 vec = {
    query[0] - p0[0], query[1] - p0[1], query[2] - p0[2]
  };
  const double dotp = dot3(vec, normal);
  const geo::Point3 proj = {
    query[0] - normal[0] * dotp,
    query[1] - normal[1] * dotp,
    query[2] - normal[2] * dotp,
  };

  // Barycentric coordinates via sub-triangle areas.
  const geo::Point3 *verts[3] = {&p0, &p1, &p2};
  double bary[3];
  double bary_sum = 0.0;

  // area_con: opposite vertex pairs for each barycentric coordinate
  constexpr int area_con[3][2] = {{1, 2}, {2, 0}, {0, 1}};
  for(int n = 0; n < 3; ++n) {
    const int i0 = area_con[n][0];
    const int i1 = area_con[n][1];
    const geo::Vector3 d0 = {
      (*verts[i0])[0] - proj[0],
      (*verts[i0])[1] - proj[1],
      (*verts[i0])[2] - proj[2],
    };
    const geo::Vector3 d1 = {
      (*verts[i1])[0] - proj[0],
      (*verts[i1])[1] - proj[1],
      (*verts[i1])[2] - proj[2],
    };
    const geo::Vector3 cross = {
      d0[1] * d1[2] - d1[1] * d0[2],
      d0[2] * d1[0] - d1[2] * d0[0],
      d0[0] * d1[1] - d1[0] * d0[1],
    };
    bary[n] = dot3(cross, normal) / area;
    if(bary[n] < 0.0) bary[n] = 0.0;
    bary_sum += bary[n];
  }

  if(bary_sum > 1.0) {
    bary[0] /= bary_sum;
    bary[1] /= bary_sum;
    bary[2] /= bary_sum;
  }

  // Closest point on triangle (clamped barycentric)
  const geo::Point3 closest = {
    p0[0] * bary[0] + p1[0] * bary[1] + p2[0] * bary[2],
    p0[1] * bary[0] + p1[1] * bary[1] + p2[1] * bary[2],
    p0[2] * bary[0] + p1[2] * bary[1] + p2[2] * bary[2],
  };

  dist_out = point_distance(closest, query);
  size_out = size0 * bary[0] + size1 * bary[1] + size2 * bary[2];
}

} // anonymous namespace

void SizeFunction::clear() noexcept
{
  sources_.clear();
  edge_sources_.clear();
  face_sources_.clear();
  octree_.clear();
  stats_ = {};
  built_ = false;
  domain_set_ = false;
}

void SizeFunction::configure(const SizeFunctionConfig &config)
{
  config_ = config;
}

void SizeFunction::set_domain(
  const geo::Point3 &domain_min,
  const geo::Point3 &domain_max
)
{
  domain_min_ = domain_min;
  domain_max_ = domain_max;
  domain_set_ = true;
}

void SizeFunction::add_source(const geo::Point3 &position, double size_value)
{
  if(!std::isfinite(size_value) || size_value <= 0.0) {
    return;
  }
  sources_.push_back({position, size_value});
}

void SizeFunction::add_sources(
  const geo::Point3 *positions,
  const double *size_values,
  std::size_t count
)
{
  for(std::size_t i = 0U; i < count; ++i) {
    add_source(positions[i], size_values[i]);
  }
}

void SizeFunction::add_edge_source(
  const geo::Point3 &p0, double size0,
  const geo::Point3 &p1, double size1
)
{
  if(!std::isfinite(size0) || size0 <= 0.0 ||
     !std::isfinite(size1) || size1 <= 0.0) {
    return;
  }
  edge_sources_.push_back({p0, p1, size0, size1});
}

void SizeFunction::add_face_source(
  const geo::Point3 &p0, double size0,
  const geo::Point3 &p1, double size1,
  const geo::Point3 &p2, double size2
)
{
  if(!std::isfinite(size0) || size0 <= 0.0 ||
     !std::isfinite(size1) || size1 <= 0.0 ||
     !std::isfinite(size2) || size2 <= 0.0) {
    return;
  }

  // Compute normal and area for distance computation.
  const geo::Vector3 e01 = {p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]};
  const geo::Vector3 e02 = {p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]};
  geo::Vector3 normal = {
    e01[1]*e02[2] - e01[2]*e02[1],
    e01[2]*e02[0] - e01[0]*e02[2],
    e01[0]*e02[1] - e01[1]*e02[0],
  };
  const double area = std::sqrt(dot3(normal, normal));
  if(area > kEpsilon) {
    normal[0] /= area;
    normal[1] /= area;
    normal[2] /= area;
  }

  face_sources_.push_back({p0, p1, p2, size0, size1, size2, normal, area});
  
}

void SizeFunction::apply_geometric_source_sizes()
{
  if(face_sources_.empty() && edge_sources_.empty()) {
    return;
  }

  const double alpha = std::max(0.0, config_.growth_rate - 1.0);


  auto &vertices = octree_.mutable_vertices();
  const auto vertex_count = vertices.size();
  if(vertex_count == 0U) {
    return;
  }

  // Get vertex positions from the octree.
  std::vector<geo::Point3> vertex_positions(vertex_count);
  octree_.get_vertex_positions(vertex_positions);

  std::size_t updates = 0U;

  // Process face sources with bounding box pruning.
  // Each source only affects vertices within its influence radius:
  // max_influence = max_msize / (growth_rate - 1) + source_size
  const double max_influence_base =
    (alpha > kEpsilon) ? config_.maximum_size / alpha : config_.maximum_size * 10.0;

  for(const auto &fs : face_sources_) {
    const double source_max_size = std::max({fs.size0, fs.size1, fs.size2});
    const double influence = max_influence_base + source_max_size;

    // Source bounding box expanded by influence radius
    const double xmin = std::min({fs.p0[0], fs.p1[0], fs.p2[0]}) - influence;
    const double xmax = std::max({fs.p0[0], fs.p1[0], fs.p2[0]}) + influence;
    const double ymin = std::min({fs.p0[1], fs.p1[1], fs.p2[1]}) - influence;
    const double ymax = std::max({fs.p0[1], fs.p1[1], fs.p2[1]}) + influence;
    const double zmin = std::min({fs.p0[2], fs.p1[2], fs.p2[2]}) - influence;
    const double zmax = std::max({fs.p0[2], fs.p1[2], fs.p2[2]}) + influence;

    for(std::size_t vi = 0U; vi < vertex_count; ++vi) {
      const auto &vp = vertex_positions[vi];
      if(vp[0] < xmin || vp[0] > xmax ||
         vp[1] < ymin || vp[1] > ymax ||
         vp[2] < zmin || vp[2] > zmax) {
        continue;
      }

      double dist = 0.0;
      double source_size = 0.0;
      closest_point_on_triangle(
        fs.p0, fs.p1, fs.p2,
        fs.size0, fs.size1, fs.size2,
        fs.normal, fs.area,
        vp, dist, source_size
      );

      dist = std::max(0.0, dist - source_size);
      double diffused = source_size + alpha * dist;
      diffused = std::max(config_.minimum_size,
                          std::min(config_.maximum_size, diffused));

      if(diffused < vertices[vi].mesh_size) {
        vertices[vi].mesh_size = diffused;
        ++updates;
      }
    }
  }

  // Process edge sources with bounding box pruning.
  for(const auto &es : edge_sources_) {
    const double source_max_size = std::max(es.size0, es.size1);
    const double influence = max_influence_base + source_max_size;

    const double xmin = std::min(es.p0[0], es.p1[0]) - influence;
    const double xmax = std::max(es.p0[0], es.p1[0]) + influence;
    const double ymin = std::min(es.p0[1], es.p1[1]) - influence;
    const double ymax = std::max(es.p0[1], es.p1[1]) + influence;
    const double zmin = std::min(es.p0[2], es.p1[2]) - influence;
    const double zmax = std::max(es.p0[2], es.p1[2]) + influence;

    for(std::size_t vi = 0U; vi < vertex_count; ++vi) {
      const auto &vp = vertex_positions[vi];
      if(vp[0] < xmin || vp[0] > xmax ||
         vp[1] < ymin || vp[1] > ymax ||
         vp[2] < zmin || vp[2] > zmax) {
        continue;
      }

      double dist = 0.0;
      double source_size = 0.0;
      closest_point_on_edge(
        es.p0, es.p1, es.size0, es.size1,
        vp, dist, source_size
      );

      dist = std::max(0.0, dist - source_size);
      double diffused = source_size + alpha * dist;
      diffused = std::max(config_.minimum_size,
                          std::min(config_.maximum_size, diffused));

      if(diffused < vertices[vi].mesh_size) {
        vertices[vi].mesh_size = diffused;
        ++updates;
      }
    }
  }

}

void SizeFunction::build()
{
  built_ = false;
  octree_.clear();

  if(sources_.empty() && face_sources_.empty() && edge_sources_.empty()) {
    return;
  }

  // Compute domain bounds from all sources (point + face + edge) and domain.
  const bool has_point = !sources_.empty();
  geo::Point3 bounds_min = has_point
    ? sources_.front().position
    : (!face_sources_.empty()
       ? face_sources_.front().p0
       : edge_sources_.front().p0);
  geo::Point3 bounds_max = bounds_min;

  const auto extend_bounds = [&](const geo::Point3 &p) {
    for(std::size_t axis = 0U; axis < 3U; ++axis) {
      bounds_min[axis] = std::min(bounds_min[axis], p[axis]);
      bounds_max[axis] = std::max(bounds_max[axis], p[axis]);
    }
  };
  for(const auto &s : sources_)       { extend_bounds(s.position); }
  for(const auto &fs : face_sources_) {
    extend_bounds(fs.p0); extend_bounds(fs.p1); extend_bounds(fs.p2);
  }
  for(const auto &es : edge_sources_) {
    extend_bounds(es.p0); extend_bounds(es.p1);
  }

  if(domain_set_) {
    for(std::size_t axis = 0U; axis < 3U; ++axis) {
      bounds_min[axis] = std::min(bounds_min[axis], domain_min_[axis]);
      bounds_max[axis] = std::max(bounds_max[axis], domain_max_[axis]);
    }
  }

  // Add padding based on maximum_size to allow field to propagate beyond sources.
  const double padding = std::max(config_.maximum_size, 1.0e-6);
  for(std::size_t axis = 0U; axis < 3U; ++axis) {
    bounds_min[axis] -= padding;
    bounds_max[axis] += padding;
  }

  // Initialize octree.
  AdaptiveOctreeConfig octree_config;
  octree_config.max_depth = config_.max_octree_depth;
  octree_config.minimum_node_extent = std::max(config_.minimum_size * 0.1, 1.0e-12);
  octree_.initialize(bounds_min, bounds_max, octree_config);

  // Phase 1: insert sizing seeds from every source type. The octree
  // structure is driven by these seeds, so we need coverage across the
  // surface (face-source vertices + centroid) and along edges
  // (edge-source endpoints) — not just point sources.
  for(const auto &source : sources_) {
    octree_.insert_size_seed(source.position, source.size_value);
  }
  for(const auto &fs : face_sources_) {
    octree_.insert_size_seed(fs.p0, fs.size0);
    octree_.insert_size_seed(fs.p1, fs.size1);
    octree_.insert_size_seed(fs.p2, fs.size2);
    const geo::Point3 centroid = {
      (fs.p0[0] + fs.p1[0] + fs.p2[0]) / 3.0,
      (fs.p0[1] + fs.p1[1] + fs.p2[1]) / 3.0,
      (fs.p0[2] + fs.p1[2] + fs.p2[2]) / 3.0,
    };
    const double size_centroid = (fs.size0 + fs.size1 + fs.size2) / 3.0;
    octree_.insert_size_seed(centroid, size_centroid);
  }
  for(const auto &es : edge_sources_) {
    octree_.insert_size_seed(es.p0, es.size0);
    octree_.insert_size_seed(es.p1, es.size1);
  }

  auto seed_stats = octree_.stats();

  // Phase 2: 2:1 balance the octree so adjacent leaves differ by at most 1 level.
  octree_.balance_2_to_1();

  auto balance_stats = octree_.stats();

  // Phase 3: Dijkstra gradient smoothing.
  octree_.smooth_gradients(config_.growth_rate, config_.minimum_size, config_.maximum_size);

  // Phase 4: Compute per-leaf gradient vectors for interpolation.
  octree_.compute_leaf_gradients();

  // Phase 5: distance-based size updates from face/edge sources.
  // Distance is measured in global 3D space, so a source triangle can
  // bleed size into an unrelated nearby face (e.g. a wing source can
  // reach the symmetry plane through 3D proximity).
  apply_geometric_source_sizes();

  stats_.source_count = sources_.size();
  stats_.edge_source_count = edge_sources_.size();
  stats_.face_source_count = face_sources_.size();
  stats_.octree_stats = octree_.stats();
  built_ = true;

  SQMESH_LOG_INFO("SizeFunction built: {} point sources, {} face sources "
    "({} octree seeds), {} edge sources",
    sources_.size(),
    face_sources_.size(),
    sources_.size() + 4U * face_sources_.size() + 2U * edge_sources_.size(),
    edge_sources_.size());
}

double SizeFunction::query(const geo::Point3 &position) const noexcept
{
  if(!built_) {
    return config_.maximum_size;
  }

  const double raw = octree_.query(position);
  return clamp(raw);
}

double SizeFunction::clamp(double size) const noexcept
{
  if(!std::isfinite(size) || size <= 0.0) {
    return config_.maximum_size;
  }
  return std::max(config_.minimum_size, std::min(config_.maximum_size, size));
}

bool SizeFunction::built() const noexcept
{
  return built_;
}

const SizeFunctionStats &SizeFunction::stats() const noexcept
{
  return stats_;
}

} // namespace sqmesh::mesh::detail
