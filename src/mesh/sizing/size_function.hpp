// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#pragma once

#include "adaptive_octree.hpp"
#include "sqmesh/geo/api.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace sqmesh::mesh::detail {

struct SizeFunctionConfig final {
  double minimum_size = 0.0;
  double maximum_size = 0.0;
  double growth_rate = 1.2;
  std::uint8_t max_octree_depth = 14;
};

struct SizeFunctionStats final {
  std::size_t source_count = 0U;
  std::size_t edge_source_count = 0U;
  std::size_t face_source_count = 0U;
  AdaptiveOctreeStats octree_stats {};
};

class SizeFunction final {
public:
  void clear() noexcept;

  /// Configure size bounds and growth rate. Must be called before build().
  void configure(const SizeFunctionConfig &config);

  /// Set the spatial domain. Must be called before build().
  void set_domain(
    const geo::Point3 &domain_min,
    const geo::Point3 &domain_max
  );

  /// Add a single sizing source.
  void add_source(const geo::Point3 &position, double size_value);

  /// Add multiple sizing sources at once.
  void add_sources(
    const geo::Point3 *positions,
    const double *size_values,
    std::size_t count
  );

  /// Add an edge source (two endpoints with sizes).
  /// Stored as a geometric entity for distance-based size computation.
  void add_edge_source(
    const geo::Point3 &p0, double size0,
    const geo::Point3 &p1, double size1
  );

  /// Add a face source (triangle with per-vertex sizes).
  /// Stored as a geometric entity for distance-based size computation.
  void add_face_source(
    const geo::Point3 &p0, double size0,
    const geo::Point3 &p1, double size1,
    const geo::Point3 &p2, double size2
  );

  /// Build the adaptive octree and apply gradient smoothing.
  /// After this call, query() is valid.
  void build();

  /// Query the size function at an arbitrary point.
  /// Returns the interpolated mesh size, clamped to [min_size, max_size].
  [[nodiscard]] double query(const geo::Point3 &position) const noexcept;

  /// Clamp a raw size value to the configured bounds.
  [[nodiscard]] double clamp(double size) const noexcept;

  [[nodiscard]] bool built() const noexcept;
  [[nodiscard]] const SizeFunctionStats &stats() const noexcept;

private:
  struct Source final {
    geo::Point3 position {0.0, 0.0, 0.0};
    double size_value = 0.0;
  };

  struct EdgeSource final {
    geo::Point3 p0 {0.0, 0.0, 0.0};
    geo::Point3 p1 {0.0, 0.0, 0.0};
    double size0 = 0.0;
    double size1 = 0.0;
  };

  struct FaceSource final {
    geo::Point3 p0 {0.0, 0.0, 0.0};
    geo::Point3 p1 {0.0, 0.0, 0.0};
    geo::Point3 p2 {0.0, 0.0, 0.0};
    double size0 = 0.0;
    double size1 = 0.0;
    double size2 = 0.0;
    geo::Vector3 normal {0.0, 0.0, 0.0};
    double area = 0.0;
  };

  /// Apply geometric entity distance-based size updates to octree vertices.
  void apply_geometric_source_sizes();

  SizeFunctionConfig config_ {};
  geo::Point3 domain_min_ {0.0, 0.0, 0.0};
  geo::Point3 domain_max_ {0.0, 0.0, 0.0};
  bool domain_set_ = false;
  std::vector<Source> sources_ {};
  std::vector<EdgeSource> edge_sources_ {};
  std::vector<FaceSource> face_sources_ {};
  AdaptiveOctree octree_ {};
  SizeFunctionStats stats_ {};
  bool built_ = false;
};

} // namespace sqmesh::mesh::detail
