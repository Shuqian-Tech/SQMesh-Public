// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "sqmesh/mesh/api.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace sqmesh::mesh {

namespace {

[[nodiscard]] EntityKind default_kind_for_order(EntityOrder order)
{
  switch(order) {
  case EntityOrder::node:
    return EntityKind::node_point;
  case EntityOrder::edge:
    return EntityKind::edge_line;
  case EntityOrder::face:
    return EntityKind::face_triangle;
  case EntityOrder::cell:
    return EntityKind::cell_tetra;
  default:
    throw std::invalid_argument("Unsupported mesh entity order.");
  }
}

[[nodiscard]] EntityKind normalize_kind(EntityOrder order, EntityKind kind)
{
  const auto normalized = kind == EntityKind::invalid ? default_kind_for_order(order) : kind;
  if(entity_order(normalized) != order) {
    throw std::invalid_argument("EntityGroup/entity kind does not match the requested entity_group order.");
  }
  return normalized;
}

[[nodiscard]] EntityGroupSemantic default_semantic_for_entity_group(
  EntityOrder order,
  bool boundary
) noexcept
{
  switch(order) {
  case EntityOrder::node:
    return EntityGroupSemantic::node;
  case EntityOrder::cell:
    return EntityGroupSemantic::region;
  case EntityOrder::edge:
  case EntityOrder::face:
    return boundary ? EntityGroupSemantic::boundary : EntityGroupSemantic::interior;
  default:
    return EntityGroupSemantic::unspecified;
  }
}

[[nodiscard]] EntityGroupSemantic normalize_entity_group_semantic(
  EntityOrder order,
  bool boundary,
  EntityGroupSemantic semantic
)
{
  const auto normalized =
    semantic == EntityGroupSemantic::unspecified
      ? default_semantic_for_entity_group(order, boundary)
      : semantic;

  switch(normalized) {
  case EntityGroupSemantic::node:
    if(order != EntityOrder::node) {
      throw std::invalid_argument("Only node entity_groups can use node entity_group semantics.");
    }
    return normalized;
  case EntityGroupSemantic::region:
    if(order != EntityOrder::cell) {
      throw std::invalid_argument("Only cell entity_groups can use region entity_group semantics.");
    }
    return normalized;
  case EntityGroupSemantic::boundary:
    if((order != EntityOrder::edge && order != EntityOrder::face) || !boundary) {
      throw std::invalid_argument(
        "Boundary entity_group semantics require boundary edge or face entity_groups."
      );
    }
    return normalized;
  case EntityGroupSemantic::interior:
    if((order != EntityOrder::edge && order != EntityOrder::face) || boundary) {
      throw std::invalid_argument(
        "Interior entity_group semantics require non-boundary edge or face entity_groups."
      );
    }
    return normalized;
  case EntityGroupSemantic::interface:
    if((order != EntityOrder::edge && order != EntityOrder::face) || boundary) {
      throw std::invalid_argument(
        "Interface entity_group semantics require non-boundary edge or face entity_groups."
      );
    }
    return normalized;
  default:
    throw std::invalid_argument("Unsupported mesh entity_group semantics.");
  }
}

void normalize_entity_group_region_zone_ids(
  EntityGroupSemantic semantic,
  std::uint32_t zone_id,
  std::uint32_t &primary_region_zone_id,
  std::uint32_t &secondary_region_zone_id
)
{
  switch(semantic) {
  case EntityGroupSemantic::node:
    if(primary_region_zone_id != invalid_index ||
       secondary_region_zone_id != invalid_index) {
      throw std::invalid_argument("Node entity_groups cannot carry region-zone bindings.");
    }
    return;
  case EntityGroupSemantic::interior:
    if(secondary_region_zone_id != invalid_index) {
      throw std::invalid_argument(
        "Interior entity_groups can reference at most one owning region-zone id."
      );
    }
    return;
  case EntityGroupSemantic::boundary:
    if(secondary_region_zone_id != invalid_index) {
      throw std::invalid_argument(
        "Boundary entity_groups can reference at most one adjacent region-zone id."
      );
    }
    return;
  case EntityGroupSemantic::interface:
    if(primary_region_zone_id == invalid_index ||
       secondary_region_zone_id == invalid_index) {
      throw std::invalid_argument(
        "Interface entity_groups require two adjacent region-zone ids."
      );
    }
    if(primary_region_zone_id == secondary_region_zone_id) {
      throw std::invalid_argument(
        "Interface entity_groups require distinct adjacent region-zone ids."
      );
    }
    if(secondary_region_zone_id < primary_region_zone_id) {
      std::swap(primary_region_zone_id, secondary_region_zone_id);
    }
    return;
  case EntityGroupSemantic::region:
    if(primary_region_zone_id == invalid_index) {
      primary_region_zone_id = zone_id;
    }
    if(secondary_region_zone_id != invalid_index) {
      throw std::invalid_argument("Region entity_groups can bind only one region-zone id.");
    }
    return;
  default:
    throw std::invalid_argument("Unsupported mesh entity_group semantics.");
  }
}

[[nodiscard]] std::uint32_t make_entity_flags(
  EntityOrder order,
  bool boundary,
  std::uint32_t extra_flags
) noexcept
{
  const auto sanitized_flags =
    extra_flags & ~(entity_order_mask | entity_boundary_mask);
  return sanitized_flags | entity_live_mask | entity_meshed_mask |
         static_cast<std::uint32_t>(order) |
         (boundary ? entity_boundary_mask : 0U);
}

[[nodiscard]] bool has_entity_group_prefix(
  std::string_view name,
  std::string_view prefix
) noexcept
{
  return name.size() >= prefix.size() && name.substr(0, prefix.size()) == prefix;
}

[[nodiscard]] EntityRef remap_entity_ref(
  EntityRef ref,
  const std::vector<EntityGroupIndex> &entity_group_remap
) noexcept
{
  if(!is_valid(ref) || ref.entity_group >= entity_group_remap.size()) {
    return {};
  }

  const auto remapped_entity_group = entity_group_remap[ref.entity_group];
  if(remapped_entity_group == invalid_index) {
    return {};
  }

  ref.entity_group = remapped_entity_group;
  return ref;
}

[[nodiscard]] const EntityGroup &require_entity_group(
  const Domain &domain,
  EntityGroupIndex entity_group_index,
  EntityOrder order
)
{
  const auto &entity_group = domain.entity_group(entity_group_index);
  if(entity_group.order() != order) {
    throw std::invalid_argument("EntityGroup order mismatch for requested mesh entity access.");
  }
  return entity_group;
}

void validate_entity_group_kind(const EntityGroup &entity_group, EntityKind expected_kind)
{
  if(entity_group.default_kind() != expected_kind) {
    throw std::invalid_argument("EntityGroup kind does not match the requested entity type.");
  }
}

void validate_node_ref(const Domain &domain, EntityRef node_ref)
{
  if(!is_valid(node_ref)) {
    throw std::invalid_argument("Node connectivity requires valid node references.");
  }
  static_cast<void>(domain.node(node_ref));
}

void validate_face_ref(const Domain &domain, EntityRef face_ref)
{
  if(!is_valid(face_ref)) {
    throw std::invalid_argument("Cell connectivity requires valid face references.");
  }
  static_cast<void>(domain.face(face_ref));
}

void validate_cell_ref(const Domain &domain, EntityRef cell_ref)
{
  if(!is_valid(cell_ref)) {
    throw std::invalid_argument("Face adjacency requires valid cell references.");
  }
  static_cast<void>(domain.cell(cell_ref));
}

[[nodiscard]] EntityRefRange make_range(
  const std::vector<EntityRef> &channel,
  ConnectivitySpan span
) noexcept
{
  return {
    span.count == 0U ? nullptr : channel.data() + span.offset,
    span.count,
  };
}

using Point = std::array<double, 3>;

constexpr double kPi = 3.14159265358979323846264338327950288;
constexpr double kTriangleIdealAngle = kPi / 3.0;
constexpr double kTriangleIdealSin = 0.86602540378443864676372317075293618;
constexpr double kTriangleAspectNormalization = 0.86602540378443864676372317075293618;
constexpr double kTetraIdealDihedral = 1.230959417340774682134929450118454;
constexpr double kTetraAspectNormalization = 0.8164965809277260327324280249019638;
constexpr double kTetraJacobianNormalization = 1.4142135623730950488016887242096981;
constexpr double kRadiansToDegrees = 180.0 / kPi;
constexpr double kQualityTolerance = 1.0e-12;

struct QualityMetricAccumulator final {
  std::size_t count = 0;
  double minimum = 0.0;
  double maximum = 0.0;
  double sum = 0.0;

  void observe(double value) noexcept
  {
    if(!std::isfinite(value)) {
      return;
    }

    if(count == 0U) {
      minimum = value;
      maximum = value;
    }
    else {
      minimum = std::min(minimum, value);
      maximum = std::max(maximum, value);
    }

    sum += value;
    ++count;
  }

  [[nodiscard]] QualityMetricSummary summary() const noexcept
  {
    QualityMetricSummary metric;
    metric.count = count;
    if(count == 0U) {
      return metric;
    }

    metric.minimum = minimum;
    metric.maximum = maximum;
    metric.average = sum / static_cast<double>(count);
    return metric;
  }
};

struct KindQualityAccumulator final {
  QualityMetricAccumulator jacobian {};
  QualityMetricAccumulator skewness {};
  QualityMetricAccumulator aspect_ratio {};
  QualityMetricAccumulator radius_ratio {};
  QualityMetricAccumulator min_angle {};
  QualityMetricAccumulator max_angle {};
};

[[nodiscard]] double clamp_unit(double value) noexcept
{
  return std::max(-1.0, std::min(1.0, value));
}

[[nodiscard]] double clamp_unit_interval(double value) noexcept
{
  return std::max(0.0, std::min(1.0, value));
}

[[nodiscard]] Point subtract_points(
  const Point &lhs,
  const Point &rhs
) noexcept
{
  return {
    lhs[0] - rhs[0],
    lhs[1] - rhs[1],
    lhs[2] - rhs[2],
  };
}

[[nodiscard]] double dot_product(
  const Point &lhs,
  const Point &rhs
) noexcept
{
  return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2];
}

[[nodiscard]] Point cross_product(
  const Point &lhs,
  const Point &rhs
) noexcept
{
  return {
    lhs[1] * rhs[2] - lhs[2] * rhs[1],
    lhs[2] * rhs[0] - lhs[0] * rhs[2],
    lhs[0] * rhs[1] - lhs[1] * rhs[0],
  };
}

[[nodiscard]] Point scale_point(
  const Point &value,
  double scale
) noexcept
{
  return {
    value[0] * scale,
    value[1] * scale,
    value[2] * scale,
  };
}

[[nodiscard]] Point add_points(
  const Point &lhs,
  const Point &rhs
) noexcept
{
  return {
    lhs[0] + rhs[0],
    lhs[1] + rhs[1],
    lhs[2] + rhs[2],
  };
}

[[nodiscard]] double vector_norm(const Point &value) noexcept
{
  return std::sqrt(dot_product(value, value));
}

[[nodiscard]] double triple_product(
  const Point &a,
  const Point &b,
  const Point &c
) noexcept
{
  return dot_product(a, cross_product(b, c));
}

[[nodiscard]] double geometric_tolerance(
  double length_scale,
  int dimension
) noexcept
{
  double tolerance = kQualityTolerance;
  const double scale = std::max(length_scale, 1.0);
  for(int index = 0; index < dimension; ++index) {
    tolerance *= scale;
  }
  return tolerance;
}

[[nodiscard]] double corner_angle(
  const Point &origin,
  const Point &first,
  const Point &second
) noexcept
{
  const auto a = subtract_points(first, origin);
  const auto b = subtract_points(second, origin);
  const auto a_norm = vector_norm(a);
  const auto b_norm = vector_norm(b);
  if(a_norm <= 0.0 || b_norm <= 0.0) {
    return 0.0;
  }

  return std::acos(clamp_unit(dot_product(a, b) / (a_norm * b_norm)));
}

[[nodiscard]] double corner_jacobian_sine(
  const Point &origin,
  const Point &first,
  const Point &second
) noexcept
{
  const auto a = subtract_points(first, origin);
  const auto b = subtract_points(second, origin);
  const auto denominator = vector_norm(a) * vector_norm(b);
  if(denominator <= 0.0) {
    return 0.0;
  }

  return vector_norm(cross_product(a, b)) / denominator;
}

[[nodiscard]] Point outward_triangle_normal(
  const Point &a,
  const Point &b,
  const Point &c,
  const Point &opposite_vertex
) noexcept
{
  auto normal = cross_product(subtract_points(b, a), subtract_points(c, a));
  if(dot_product(normal, subtract_points(opposite_vertex, a)) > 0.0) {
    normal[0] = -normal[0];
    normal[1] = -normal[1];
    normal[2] = -normal[2];
  }

  const auto magnitude = vector_norm(normal);
  if(magnitude <= 0.0) {
    return {0.0, 0.0, 0.0};
  }

  return {
    normal[0] / magnitude,
    normal[1] / magnitude,
    normal[2] / magnitude,
  };
}

[[nodiscard]] double tetra_circumradius(
  const Point &a,
  const Point &b,
  const Point &c,
  double signed_det,
  double volume_tolerance
) noexcept
{
  if(std::abs(signed_det) <= volume_tolerance) {
    return 0.0;
  }

  const auto numerator = add_points(
    add_points(
      scale_point(cross_product(b, c), dot_product(a, a)),
      scale_point(cross_product(c, a), dot_product(b, b))
    ),
    scale_point(cross_product(a, b), dot_product(c, c))
  );
  return vector_norm(scale_point(numerator, 1.0 / (2.0 * signed_det)));
}

[[nodiscard]] std::size_t find_or_append_kind_summary(
  std::vector<KindQualitySummary> &summaries,
  std::vector<KindQualityAccumulator> &accumulators,
  EntityKind kind
)
{
  for(std::size_t index = 0; index < summaries.size(); ++index) {
    if(summaries[index].kind == kind) {
      return index;
    }
  }

  summaries.push_back(KindQualitySummary {});
  summaries.back().kind = kind;
  accumulators.push_back(KindQualityAccumulator {});
  return summaries.size() - 1U;
}

[[nodiscard]] ElementQuality make_triangle_quality(
  const Domain &domain,
  EntityRef face_ref
)
{
  ElementQuality report;
  report.entity = face_ref;
  report.kind = EntityKind::face_triangle;
  report.supported = true;

  const auto nodes = domain.face_nodes(face_ref);
  const auto &p0 = domain.node(nodes[0]).coordinates;
  const auto &p1 = domain.node(nodes[1]).coordinates;
  const auto &p2 = domain.node(nodes[2]).coordinates;

  const auto e01 = subtract_points(p1, p0);
  const auto e12 = subtract_points(p2, p1);
  const auto e20 = subtract_points(p0, p2);
  const auto l01 = vector_norm(e01);
  const auto l12 = vector_norm(e12);
  const auto l20 = vector_norm(e20);
  const auto longest_edge = std::max({l01, l12, l20});
  const auto area = 0.5 * vector_norm(cross_product(e01, subtract_points(p2, p0)));
  const auto area_tolerance = geometric_tolerance(longest_edge, 2);

  const auto angle0 = corner_angle(p0, p1, p2);
  const auto angle1 = corner_angle(p1, p2, p0);
  const auto angle2 = corner_angle(p2, p0, p1);
  const auto min_angle = std::min({angle0, angle1, angle2});
  const auto max_angle = std::max({angle0, angle1, angle2});

  const auto corner0 = corner_jacobian_sine(p0, p1, p2);
  const auto corner1 = corner_jacobian_sine(p1, p2, p0);
  const auto corner2 = corner_jacobian_sine(p2, p0, p1);
  report.jacobian =
    clamp_unit_interval(std::min({corner0, corner1, corner2}) / kTriangleIdealSin);
  report.jacobian_is_signed = false;

  report.skewness = clamp_unit_interval(
    std::max(
      (max_angle - kTriangleIdealAngle) / (kPi - kTriangleIdealAngle),
      (kTriangleIdealAngle - min_angle) / kTriangleIdealAngle
    )
  );

  if(area <= area_tolerance || longest_edge <= area_tolerance) {
    report.degenerate = true;
    report.status = QualityStatus::degenerate;
    report.jacobian = 0.0;
    report.skewness = 1.0;
    report.aspect_ratio = std::numeric_limits<double>::infinity();
    report.radius_ratio = 0.0;
    return report;
  }

  const auto shortest_altitude = (2.0 * area) / longest_edge;
  if(shortest_altitude <= area_tolerance) {
    report.degenerate = true;
    report.status = QualityStatus::degenerate;
    report.jacobian = 0.0;
    report.skewness = 1.0;
    report.aspect_ratio = std::numeric_limits<double>::infinity();
    report.radius_ratio = 0.0;
    return report;
  }

  const auto perimeter = l01 + l12 + l20;
  const auto radius_ratio_denominator = perimeter * l01 * l12 * l20;
  report.aspect_ratio =
    kTriangleAspectNormalization * longest_edge / shortest_altitude;
  report.radius_ratio =
    radius_ratio_denominator <= geometric_tolerance(longest_edge, 4)
      ? 0.0
      : clamp_unit_interval((16.0 * area * area) / radius_ratio_denominator);
  report.min_angle = min_angle * kRadiansToDegrees;
  report.max_angle = max_angle * kRadiansToDegrees;
  report.status = QualityStatus::valid;
  return report;
}

[[nodiscard]] ElementQuality make_tetra_quality(
  const Domain &domain,
  EntityRef cell_ref
)
{
  ElementQuality report;
  report.entity = cell_ref;
  report.kind = EntityKind::cell_tetra;
  report.supported = true;
  report.jacobian_is_signed = true;

  const auto nodes = domain.cell_nodes(cell_ref);
  const auto &p0 = domain.node(nodes[0]).coordinates;
  const auto &p1 = domain.node(nodes[1]).coordinates;
  const auto &p2 = domain.node(nodes[2]).coordinates;
  const auto &p3 = domain.node(nodes[3]).coordinates;

  const auto e01 = subtract_points(p1, p0);
  const auto e02 = subtract_points(p2, p0);
  const auto e03 = subtract_points(p3, p0);
  const auto e12 = subtract_points(p2, p1);
  const auto e13 = subtract_points(p3, p1);
  const auto e23 = subtract_points(p3, p2);

  const auto edge_lengths = std::array<double, 6> {
    vector_norm(e01),
    vector_norm(e02),
    vector_norm(e03),
    vector_norm(e12),
    vector_norm(e13),
    vector_norm(e23),
  };
  const auto longest_edge = std::max({
    edge_lengths[0],
    edge_lengths[1],
    edge_lengths[2],
    edge_lengths[3],
    edge_lengths[4],
    edge_lengths[5],
  });
  const auto volume_tolerance = geometric_tolerance(longest_edge, 3);
  const auto signed_det = triple_product(e01, e02, e03);
  const auto absolute_det = std::abs(signed_det);
  const auto orientation_sign = signed_det < 0.0 ? -1.0 : 1.0;

  const auto face_areas = std::array<double, 4> {
    0.5 * vector_norm(cross_product(subtract_points(p2, p1), subtract_points(p3, p1))),
    0.5 * vector_norm(cross_product(subtract_points(p3, p0), subtract_points(p2, p0))),
    0.5 * vector_norm(cross_product(subtract_points(p1, p0), subtract_points(p3, p0))),
    0.5 * vector_norm(cross_product(subtract_points(p2, p0), subtract_points(p1, p0))),
  };

  double min_scaled_corner_jacobian = std::numeric_limits<double>::infinity();
  const auto update_scaled_corner = [&](const Point &a, const Point &b, const Point &c) {
    const auto denominator = vector_norm(a) * vector_norm(b) * vector_norm(c);
    if(denominator <= volume_tolerance) {
      min_scaled_corner_jacobian = 0.0;
      return;
    }

    const auto scaled =
      kTetraJacobianNormalization * std::abs(triple_product(a, b, c)) / denominator;
    min_scaled_corner_jacobian = std::min(min_scaled_corner_jacobian, scaled);
  };

  update_scaled_corner(subtract_points(p1, p0), subtract_points(p2, p0), subtract_points(p3, p0));
  update_scaled_corner(subtract_points(p0, p1), subtract_points(p3, p1), subtract_points(p2, p1));
  update_scaled_corner(subtract_points(p0, p2), subtract_points(p1, p2), subtract_points(p3, p2));
  update_scaled_corner(subtract_points(p0, p3), subtract_points(p2, p3), subtract_points(p1, p3));
  report.jacobian =
    std::isfinite(min_scaled_corner_jacobian) ? orientation_sign * min_scaled_corner_jacobian : 0.0;

  const auto normals = std::array<Point, 4> {
    outward_triangle_normal(p1, p2, p3, p0),
    outward_triangle_normal(p0, p3, p2, p1),
    outward_triangle_normal(p0, p1, p3, p2),
    outward_triangle_normal(p0, p2, p1, p3),
  };
  const auto dihedral_pairs = std::array<std::array<int, 2>, 6> {{
    {0, 1},
    {0, 2},
    {0, 3},
    {1, 2},
    {1, 3},
    {2, 3},
  }};

  double min_dihedral = std::numeric_limits<double>::infinity();
  double max_dihedral = 0.0;
  for(const auto &pair : dihedral_pairs) {
    const auto cosine = clamp_unit(dot_product(normals[pair[0]], normals[pair[1]]));
    const auto dihedral = kPi - std::acos(cosine);
    min_dihedral = std::min(min_dihedral, dihedral);
    max_dihedral = std::max(max_dihedral, dihedral);
  }

  report.skewness = clamp_unit_interval(
    std::max(
      (max_dihedral - kTetraIdealDihedral) / (kPi - kTetraIdealDihedral),
      (kTetraIdealDihedral - min_dihedral) / kTetraIdealDihedral
    )
  );

  double shortest_altitude = std::numeric_limits<double>::infinity();
  const double surface_area =
    face_areas[0] + face_areas[1] + face_areas[2] + face_areas[3];
  for(const auto face_area : face_areas) {
    if(face_area <= volume_tolerance) {
      shortest_altitude = 0.0;
      break;
    }
    shortest_altitude = std::min(shortest_altitude, absolute_det / (2.0 * face_area));
  }

  const bool degenerate =
    absolute_det <= volume_tolerance ||
    longest_edge <= volume_tolerance ||
    shortest_altitude <= volume_tolerance;
  if(degenerate) {
    report.degenerate = true;
    report.status = QualityStatus::degenerate;
    report.jacobian = 0.0;
    report.skewness = 1.0;
    report.aspect_ratio = std::numeric_limits<double>::infinity();
    report.radius_ratio = 0.0;
    return report;
  }

  const auto circumradius =
    tetra_circumradius(e01, e02, e03, signed_det, volume_tolerance);
  const auto inradius =
    surface_area <= volume_tolerance ? 0.0 : absolute_det / (2.0 * surface_area);
  report.aspect_ratio =
    kTetraAspectNormalization * longest_edge / shortest_altitude;
  report.radius_ratio =
    circumradius <= volume_tolerance
      ? 0.0
      : clamp_unit_interval((3.0 * inradius) / circumradius);
  report.min_angle = min_dihedral * kRadiansToDegrees;
  report.max_angle = max_dihedral * kRadiansToDegrees;

  if(signed_det < -volume_tolerance) {
    report.inverted = true;
    report.status = QualityStatus::inverted;
    return report;
  }

  report.status = QualityStatus::valid;
  return report;
}

} // namespace

const EntityGroupInfo &EntityGroup::info() const noexcept
{
  return info_;
}

const EntityGroupImportInfo &EntityGroup::import_info() const noexcept
{
  return info_.import_info;
}

EntityGroupIndex EntityGroup::id() const noexcept
{
  return info_.id;
}

std::uint32_t EntityGroup::zone_id() const noexcept
{
  return info_.zone_id;
}

std::uint32_t EntityGroup::source_entity_tag() const noexcept
{
  return info_.source_entity_tag;
}

EntityOrder EntityGroup::order() const noexcept
{
  return info_.order;
}

bool EntityGroup::is_boundary() const noexcept
{
  return info_.boundary;
}

EntityGroupSemantic EntityGroup::semantic() const noexcept
{
  return info_.semantic;
}

EntityGroupRole EntityGroup::role() const noexcept
{
  return info_.role;
}

bool EntityGroup::is_interface() const noexcept
{
  return info_.semantic == EntityGroupSemantic::interface;
}

std::uint32_t EntityGroup::primary_region_zone_id() const noexcept
{
  return info_.primary_region_zone_id;
}

std::uint32_t EntityGroup::secondary_region_zone_id() const noexcept
{
  return info_.secondary_region_zone_id;
}

EntityKind EntityGroup::default_kind() const noexcept
{
  return info_.default_kind;
}

std::string_view EntityGroup::name() const noexcept
{
  return info_.name;
}

std::size_t EntityGroup::entity_count() const noexcept
{
  switch(info_.order) {
  case EntityOrder::node:
    return nodes_.size();
  case EntityOrder::edge:
    return edges_.size();
  case EntityOrder::face:
    return faces_.size();
  case EntityOrder::cell:
    return cells_.size();
  default:
    return 0U;
  }
}

const std::vector<Node> &EntityGroup::nodes() const noexcept
{
  return nodes_;
}

const std::vector<Edge> &EntityGroup::edges() const noexcept
{
  return edges_;
}

const std::vector<Face> &EntityGroup::faces() const noexcept
{
  return faces_;
}

const std::vector<Cell> &EntityGroup::cells() const noexcept
{
  return cells_;
}

const std::vector<EntityRef> &EntityGroup::node_channel() const noexcept
{
  return node_channel_;
}

const std::vector<EntityRef> &EntityGroup::face_channel() const noexcept
{
  return face_channel_;
}

EntityGroup::EntityGroup(EntityGroupInfo info) : info_(std::move(info))
{
}

void EntityGroup::reserve(
  std::size_t entity_count,
  std::size_t node_channel_count,
  std::size_t face_channel_count
)
{
  switch(info_.order) {
  case EntityOrder::node:
    nodes_.reserve(entity_count);
    break;
  case EntityOrder::edge:
    edges_.reserve(entity_count);
    node_channel_.reserve(node_channel_count);
    edge_topology_owner_channel_.reserve(entity_count);
    break;
  case EntityOrder::face:
    faces_.reserve(entity_count);
    node_channel_.reserve(node_channel_count);
    face_topology_owner_channel_.reserve(entity_count);
    face_source_entity_tag_channel_.reserve(entity_count);
    break;
  case EntityOrder::cell:
    cells_.reserve(entity_count);
    node_channel_.reserve(node_channel_count);
    face_channel_.reserve(face_channel_count);
    break;
  default:
    throw std::invalid_argument("Unsupported mesh entity_group storage order.");
  }
}

Domain::Domain(std::string name) : name_(std::move(name))
{
}

std::string_view Domain::name() const noexcept
{
  return name_;
}

const std::vector<EntityGroup> &Domain::entity_groups() const noexcept
{
  return entity_groups_;
}

const EntityGroup &Domain::entity_group(EntityGroupIndex entity_group_index) const
{
  if(entity_group_index >= entity_groups_.size()) {
    throw std::out_of_range("EntityGroup index is outside the domain entity_group table.");
  }
  return entity_groups_[entity_group_index];
}

std::size_t Domain::entity_group_count() const noexcept
{
  return entity_groups_.size();
}

std::size_t Domain::entity_group_count(EntityOrder order) const noexcept
{
  std::size_t count = 0;
  for(const auto &entity_group_entry : entity_groups_) {
    if(entity_group_entry.order() == order) {
      ++count;
    }
  }
  return count;
}

std::size_t Domain::entity_group_count(EntityGroupRole role) const noexcept
{
  std::size_t count = 0;
  for(const auto &entity_group_entry : entity_groups_) {
    if(entity_group_entry.role() == role) {
      ++count;
    }
  }
  return count;
}

void Domain::remove_entity_groups_with_prefix(std::string_view prefix)
{
  if(prefix.empty()) {
    return;
  }

  std::vector<EntityGroup> original_entity_groups = std::move(entity_groups_);
  std::vector<EntityGroupIndex> entity_group_remap(original_entity_groups.size(), invalid_index);

  entity_groups_.clear();
  entity_groups_.reserve(original_entity_groups.size());

  for(std::size_t old_index = 0; old_index < original_entity_groups.size(); ++old_index) {
    if(has_entity_group_prefix(original_entity_groups[old_index].name(), prefix)) {
      continue;
    }

    const auto new_index = static_cast<EntityGroupIndex>(entity_groups_.size());
    entity_group_remap[old_index] = new_index;
    entity_groups_.push_back(std::move(original_entity_groups[old_index]));
    entity_groups_.back().info_.id = new_index;
  }

  for(auto &entity_group_entry : entity_groups_) {
    const auto entity_group_id = entity_group_entry.info_.id;

    for(auto &node : entity_group_entry.nodes_) {
      node.header.entity_group = entity_group_id;
    }
    for(auto &edge : entity_group_entry.edges_) {
      edge.header.entity_group = entity_group_id;
      edge.left_face = remap_entity_ref(edge.left_face, entity_group_remap);
      edge.right_face = remap_entity_ref(edge.right_face, entity_group_remap);
    }
    for(auto &face : entity_group_entry.faces_) {
      face.header.entity_group = entity_group_id;
      face.left_cell = remap_entity_ref(face.left_cell, entity_group_remap);
      face.right_cell = remap_entity_ref(face.right_cell, entity_group_remap);
    }
    for(auto &cell : entity_group_entry.cells_) {
      cell.header.entity_group = entity_group_id;
    }
    for(auto &node_ref : entity_group_entry.node_channel_) {
      node_ref = remap_entity_ref(node_ref, entity_group_remap);
    }
    for(auto &face_ref : entity_group_entry.face_channel_) {
      face_ref = remap_entity_ref(face_ref, entity_group_remap);
    }
  }
}

EntityGroupIndex Domain::create_entity_group(EntityGroupDefinition definition)
{
  EntityGroupInfo info;
  info.id = static_cast<EntityGroupIndex>(entity_groups_.size());
  info.zone_id = definition.zone_id == invalid_index ? info.id : definition.zone_id;
  info.source_entity_tag = definition.source_entity_tag;
  info.order = definition.order;
  info.boundary = definition.boundary;
  info.default_kind = normalize_kind(definition.order, definition.default_kind);
  info.semantic = normalize_entity_group_semantic(
    info.order,
    info.boundary,
    definition.semantic
  );
  info.role = definition.role;
  info.primary_region_zone_id = definition.primary_region_zone_id;
  info.secondary_region_zone_id = definition.secondary_region_zone_id;
  normalize_entity_group_region_zone_ids(
    info.semantic,
    info.zone_id,
    info.primary_region_zone_id,
    info.secondary_region_zone_id
  );
  info.name = std::move(definition.name);
  info.import_info = std::move(definition.import_info);

  entity_groups_.emplace_back(std::move(info));
  return entity_groups_.back().id();
}

void Domain::reserve_entity_group_storage(
  EntityGroupIndex entity_group_index,
  std::size_t entity_count,
  std::size_t node_channel_count,
  std::size_t face_channel_count
)
{
  entity_group_mutable(entity_group_index).reserve(entity_count, node_channel_count, face_channel_count);
}

EntityRef Domain::add_node(
  EntityGroupIndex entity_group_index,
  const std::array<double, 3> &coordinates,
  std::uint32_t extra_flags
)
{
  auto &entity_group_entry = entity_group_mutable(entity_group_index);
  if(entity_group_entry.order() != EntityOrder::node) {
    throw std::invalid_argument("EntityGroup order mismatch for requested mesh entity operation.");
  }
  validate_entity_group_kind(entity_group_entry, EntityKind::node_point);

  EntityRef ref {entity_group_index, static_cast<std::uint32_t>(entity_group_entry.nodes_.size())};
  Node node;
  node.header.id = next_entity_id_++;
  node.header.flags = make_entity_flags(EntityOrder::node, entity_group_entry.is_boundary(), extra_flags);
  node.header.index = ref.index;
  node.header.kind = EntityKind::node_point;
  node.header.entity_group = entity_group_index;
  node.coordinates = coordinates;

  entity_group_entry.nodes_.push_back(node);
  return ref;
}

EntityRef Domain::add_edge(
  EntityGroupIndex entity_group_index,
  const std::array<EntityRef, 2> &nodes,
  std::uint32_t extra_flags
)
{
  auto &entity_group_entry = entity_group_mutable(entity_group_index);
  if(entity_group_entry.order() != EntityOrder::edge) {
    throw std::invalid_argument("EntityGroup order mismatch for requested mesh entity operation.");
  }
  validate_entity_group_kind(entity_group_entry, EntityKind::edge_line);
  for(const auto node_ref : nodes) {
    validate_node_ref(*this, node_ref);
  }

  const ConnectivitySpan span {
    static_cast<std::uint32_t>(entity_group_entry.node_channel_.size()),
    static_cast<std::uint32_t>(nodes.size())
  };
  entity_group_entry.node_channel_.insert(entity_group_entry.node_channel_.end(), nodes.begin(), nodes.end());

  EntityRef ref {entity_group_index, static_cast<std::uint32_t>(entity_group_entry.edges_.size())};
  Edge edge;
  edge.header.id = next_entity_id_++;
  const auto semantic_flags =
    entity_group_entry.is_interface() ? entity_interface_mask : 0U;
  edge.header.flags = make_entity_flags(
    EntityOrder::edge,
    entity_group_entry.is_boundary(),
    extra_flags | semantic_flags
  );
  edge.header.index = ref.index;
  edge.header.kind = EntityKind::edge_line;
  edge.header.entity_group = entity_group_index;
  edge.node_span = span;

  entity_group_entry.edges_.push_back(edge);
  entity_group_entry.edge_topology_owner_channel_.push_back(geo::invalid_topology_index);
  return ref;
}

EntityRef Domain::add_triangle_face(
  EntityGroupIndex entity_group_index,
  const std::array<EntityRef, 3> &nodes,
  std::uint32_t extra_flags
)
{
  auto &entity_group_entry = entity_group_mutable(entity_group_index);
  if(entity_group_entry.order() != EntityOrder::face) {
    throw std::invalid_argument("EntityGroup order mismatch for requested mesh entity operation.");
  }
  validate_entity_group_kind(entity_group_entry, EntityKind::face_triangle);
  for(const auto node_ref : nodes) {
    validate_node_ref(*this, node_ref);
  }

  const ConnectivitySpan span {
    static_cast<std::uint32_t>(entity_group_entry.node_channel_.size()),
    static_cast<std::uint32_t>(nodes.size())
  };
  entity_group_entry.node_channel_.insert(entity_group_entry.node_channel_.end(), nodes.begin(), nodes.end());

  EntityRef ref {entity_group_index, static_cast<std::uint32_t>(entity_group_entry.faces_.size())};
  Face face;
  face.header.id = next_entity_id_++;
  const auto semantic_flags =
    entity_group_entry.is_interface() ? entity_interface_mask : 0U;
  face.header.flags = make_entity_flags(
    EntityOrder::face,
    entity_group_entry.is_boundary(),
    extra_flags | semantic_flags
  );
  face.header.index = ref.index;
  face.header.kind = EntityKind::face_triangle;
  face.header.entity_group = entity_group_index;
  face.node_span = span;

  entity_group_entry.faces_.push_back(face);
  entity_group_entry.face_topology_owner_channel_.push_back(geo::invalid_topology_index);
  entity_group_entry.face_source_entity_tag_channel_.push_back(invalid_index);
  return ref;
}

EntityRef Domain::add_quad_face(
  EntityGroupIndex entity_group_index,
  const std::array<EntityRef, 4> &nodes,
  std::uint32_t extra_flags
)
{
  auto &entity_group_entry = entity_group_mutable(entity_group_index);
  if(entity_group_entry.order() != EntityOrder::face) {
    throw std::invalid_argument("EntityGroup order mismatch for requested mesh entity operation.");
  }
  validate_entity_group_kind(entity_group_entry, EntityKind::face_quad);
  for(const auto node_ref : nodes) {
    validate_node_ref(*this, node_ref);
  }

  const ConnectivitySpan span {
    static_cast<std::uint32_t>(entity_group_entry.node_channel_.size()),
    static_cast<std::uint32_t>(nodes.size())
  };
  entity_group_entry.node_channel_.insert(entity_group_entry.node_channel_.end(), nodes.begin(), nodes.end());

  EntityRef ref {entity_group_index, static_cast<std::uint32_t>(entity_group_entry.faces_.size())};
  Face face;
  face.header.id = next_entity_id_++;
  const auto semantic_flags =
    entity_group_entry.is_interface() ? entity_interface_mask : 0U;
  face.header.flags = make_entity_flags(
    EntityOrder::face,
    entity_group_entry.is_boundary(),
    extra_flags | semantic_flags
  );
  face.header.index = ref.index;
  face.header.kind = EntityKind::face_quad;
  face.header.entity_group = entity_group_index;
  face.node_span = span;

  entity_group_entry.faces_.push_back(face);
  entity_group_entry.face_topology_owner_channel_.push_back(geo::invalid_topology_index);
  entity_group_entry.face_source_entity_tag_channel_.push_back(invalid_index);
  return ref;
}

EntityRef Domain::add_tetra_cell(
  EntityGroupIndex entity_group_index,
  const std::array<EntityRef, 4> &nodes,
  const std::array<EntityRef, 4> &faces,
  std::uint32_t extra_flags
)
{
  auto &entity_group_entry = entity_group_mutable(entity_group_index);
  if(entity_group_entry.order() != EntityOrder::cell) {
    throw std::invalid_argument("EntityGroup order mismatch for requested mesh entity operation.");
  }
  validate_entity_group_kind(entity_group_entry, EntityKind::cell_tetra);
  for(const auto node_ref : nodes) {
    validate_node_ref(*this, node_ref);
  }
  for(const auto face_ref : faces) {
    validate_face_ref(*this, face_ref);
  }

  const ConnectivitySpan node_span {
    static_cast<std::uint32_t>(entity_group_entry.node_channel_.size()),
    static_cast<std::uint32_t>(nodes.size())
  };
  entity_group_entry.node_channel_.insert(entity_group_entry.node_channel_.end(), nodes.begin(), nodes.end());

  const ConnectivitySpan face_span {
    static_cast<std::uint32_t>(entity_group_entry.face_channel_.size()),
    static_cast<std::uint32_t>(faces.size())
  };
  entity_group_entry.face_channel_.insert(entity_group_entry.face_channel_.end(), faces.begin(), faces.end());

  EntityRef ref {entity_group_index, static_cast<std::uint32_t>(entity_group_entry.cells_.size())};
  Cell cell;
  cell.header.id = next_entity_id_++;
  cell.header.flags = make_entity_flags(EntityOrder::cell, entity_group_entry.is_boundary(), extra_flags);
  cell.header.index = ref.index;
  cell.header.kind = EntityKind::cell_tetra;
  cell.header.entity_group = entity_group_index;
  cell.node_span = node_span;
  cell.face_span = face_span;

  entity_group_entry.cells_.push_back(cell);
  return ref;
}

EntityRef Domain::add_pyramid_cell(
  EntityGroupIndex entity_group_index,
  const std::array<EntityRef, 5> &nodes,
  const std::array<EntityRef, 5> &faces,
  std::uint32_t extra_flags
)
{
  auto &entity_group_entry = entity_group_mutable(entity_group_index);
  if(entity_group_entry.order() != EntityOrder::cell) {
    throw std::invalid_argument("EntityGroup order mismatch for requested mesh entity operation.");
  }
  validate_entity_group_kind(entity_group_entry, EntityKind::cell_pyramid);
  for(const auto node_ref : nodes) {
    validate_node_ref(*this, node_ref);
  }
  for(const auto face_ref : faces) {
    validate_face_ref(*this, face_ref);
  }

  const ConnectivitySpan node_span {
    static_cast<std::uint32_t>(entity_group_entry.node_channel_.size()),
    static_cast<std::uint32_t>(nodes.size())
  };
  entity_group_entry.node_channel_.insert(entity_group_entry.node_channel_.end(), nodes.begin(), nodes.end());

  const ConnectivitySpan face_span {
    static_cast<std::uint32_t>(entity_group_entry.face_channel_.size()),
    static_cast<std::uint32_t>(faces.size())
  };
  entity_group_entry.face_channel_.insert(entity_group_entry.face_channel_.end(), faces.begin(), faces.end());

  EntityRef ref {entity_group_index, static_cast<std::uint32_t>(entity_group_entry.cells_.size())};
  Cell cell;
  cell.header.id = next_entity_id_++;
  cell.header.flags = make_entity_flags(EntityOrder::cell, entity_group_entry.is_boundary(), extra_flags);
  cell.header.index = ref.index;
  cell.header.kind = EntityKind::cell_pyramid;
  cell.header.entity_group = entity_group_index;
  cell.node_span = node_span;
  cell.face_span = face_span;

  entity_group_entry.cells_.push_back(cell);
  return ref;
}

EntityRef Domain::add_prism_cell(
  EntityGroupIndex entity_group_index,
  const std::array<EntityRef, 6> &nodes,
  const std::array<EntityRef, 5> &faces,
  std::uint32_t extra_flags
)
{
  auto &entity_group_entry = entity_group_mutable(entity_group_index);
  if(entity_group_entry.order() != EntityOrder::cell) {
    throw std::invalid_argument("EntityGroup order mismatch for requested mesh entity operation.");
  }
  validate_entity_group_kind(entity_group_entry, EntityKind::cell_prism);
  for(const auto node_ref : nodes) {
    validate_node_ref(*this, node_ref);
  }
  for(const auto face_ref : faces) {
    validate_face_ref(*this, face_ref);
  }

  const ConnectivitySpan node_span {
    static_cast<std::uint32_t>(entity_group_entry.node_channel_.size()),
    static_cast<std::uint32_t>(nodes.size())
  };
  entity_group_entry.node_channel_.insert(entity_group_entry.node_channel_.end(), nodes.begin(), nodes.end());

  const ConnectivitySpan face_span {
    static_cast<std::uint32_t>(entity_group_entry.face_channel_.size()),
    static_cast<std::uint32_t>(faces.size())
  };
  entity_group_entry.face_channel_.insert(entity_group_entry.face_channel_.end(), faces.begin(), faces.end());

  EntityRef ref {entity_group_index, static_cast<std::uint32_t>(entity_group_entry.cells_.size())};
  Cell cell;
  cell.header.id = next_entity_id_++;
  cell.header.flags = make_entity_flags(EntityOrder::cell, entity_group_entry.is_boundary(), extra_flags);
  cell.header.index = ref.index;
  cell.header.kind = EntityKind::cell_prism;
  cell.header.entity_group = entity_group_index;
  cell.node_span = node_span;
  cell.face_span = face_span;

  entity_group_entry.cells_.push_back(cell);
  return ref;
}

EntityRef Domain::add_hexa_cell(
  EntityGroupIndex entity_group_index,
  const std::array<EntityRef, 8> &nodes,
  const std::array<EntityRef, 6> &faces,
  std::uint32_t extra_flags
)
{
  auto &entity_group_entry = entity_group_mutable(entity_group_index);
  if(entity_group_entry.order() != EntityOrder::cell) {
    throw std::invalid_argument("EntityGroup order mismatch for requested mesh entity operation.");
  }
  validate_entity_group_kind(entity_group_entry, EntityKind::cell_hexa);
  for(const auto node_ref : nodes) {
    validate_node_ref(*this, node_ref);
  }
  for(const auto face_ref : faces) {
    validate_face_ref(*this, face_ref);
  }

  const ConnectivitySpan node_span {
    static_cast<std::uint32_t>(entity_group_entry.node_channel_.size()),
    static_cast<std::uint32_t>(nodes.size())
  };
  entity_group_entry.node_channel_.insert(entity_group_entry.node_channel_.end(), nodes.begin(), nodes.end());

  const ConnectivitySpan face_span {
    static_cast<std::uint32_t>(entity_group_entry.face_channel_.size()),
    static_cast<std::uint32_t>(faces.size())
  };
  entity_group_entry.face_channel_.insert(entity_group_entry.face_channel_.end(), faces.begin(), faces.end());

  EntityRef ref {entity_group_index, static_cast<std::uint32_t>(entity_group_entry.cells_.size())};
  Cell cell;
  cell.header.id = next_entity_id_++;
  cell.header.flags = make_entity_flags(EntityOrder::cell, entity_group_entry.is_boundary(), extra_flags);
  cell.header.index = ref.index;
  cell.header.kind = EntityKind::cell_hexa;
  cell.header.entity_group = entity_group_index;
  cell.node_span = node_span;
  cell.face_span = face_span;

  entity_group_entry.cells_.push_back(cell);
  return ref;
}

void Domain::set_face_cells(
  EntityRef face_ref,
  EntityRef left_cell,
  EntityRef right_cell
)
{
  auto &entity_group_entry = entity_group_mutable(face_ref.entity_group);
  if(entity_group_entry.order() != EntityOrder::face) {
    throw std::invalid_argument("EntityGroup order mismatch for requested mesh entity operation.");
  }
  if(face_ref.index >= entity_group_entry.faces_.size()) {
    throw std::out_of_range("Face reference index is outside the owning face entity_group.");
  }
  validate_cell_ref(*this, left_cell);
  if(is_valid(right_cell)) {
    validate_cell_ref(*this, right_cell);
  }

  auto &face_entry = entity_group_entry.faces_[face_ref.index];
  face_entry.left_cell = left_cell;
  face_entry.right_cell = right_cell;
}

void Domain::set_edge_topology_owner(
  EntityRef edge_ref,
  geo::TopologyEntityId owner
)
{
  auto &entity_group_entry = entity_group_mutable(edge_ref.entity_group);
  if(entity_group_entry.order() != EntityOrder::edge) {
    throw std::invalid_argument("EntityGroup order mismatch for requested mesh entity operation.");
  }
  if(edge_ref.index >= entity_group_entry.edges_.size()) {
    throw std::out_of_range("Edge reference index is outside the owning edge entity_group.");
  }
  if(geo::is_valid(owner) && owner.dimension != geo::TopologyDimension::edge) {
    throw std::invalid_argument("Edge topology ownership currently supports only geometry-edge ids.");
  }

  entity_group_entry.edge_topology_owner_channel_[edge_ref.index] = owner.index;
}

geo::TopologyEntityId Domain::edge_topology_owner(EntityRef edge_ref) const
{
  const auto &entity_group_entry = require_entity_group(*this, edge_ref.entity_group, EntityOrder::edge);
  static_cast<void>(edge(edge_ref));

  if(edge_ref.index >= entity_group_entry.edge_topology_owner_channel_.size()) {
    return {geo::TopologyDimension::edge, geo::invalid_topology_index};
  }

  return {
    geo::TopologyDimension::edge,
    entity_group_entry.edge_topology_owner_channel_[edge_ref.index],
  };
}

void Domain::set_edge_faces(
  EntityRef edge_ref,
  EntityRef left_face,
  EntityRef right_face
)
{
  auto &entity_group_entry = entity_group_mutable(edge_ref.entity_group);
  if(entity_group_entry.order() != EntityOrder::edge) {
    throw std::invalid_argument("EntityGroup order mismatch for requested mesh entity operation.");
  }
  if(edge_ref.index >= entity_group_entry.edges_.size()) {
    throw std::out_of_range("Edge reference index is outside the owning edge entity_group.");
  }
  validate_face_ref(*this, left_face);
  if(is_valid(right_face)) {
    validate_face_ref(*this, right_face);
  }

  auto &edge_entry = entity_group_entry.edges_[edge_ref.index];
  edge_entry.left_face = left_face;
  edge_entry.right_face = right_face;
}

const Node &Domain::node(EntityRef node_ref) const
{
  const auto &entity_group_entry = require_entity_group(*this, node_ref.entity_group, EntityOrder::node);
  if(node_ref.index >= entity_group_entry.nodes().size()) {
    throw std::out_of_range("Node reference index is outside the owning node entity_group.");
  }
  return entity_group_entry.nodes()[node_ref.index];
}

const Edge &Domain::edge(EntityRef edge_ref) const
{
  const auto &entity_group_entry = require_entity_group(*this, edge_ref.entity_group, EntityOrder::edge);
  if(edge_ref.index >= entity_group_entry.edges().size()) {
    throw std::out_of_range("Edge reference index is outside the owning edge entity_group.");
  }
  return entity_group_entry.edges()[edge_ref.index];
}

const Face &Domain::face(EntityRef face_ref) const
{
  const auto &entity_group_entry = require_entity_group(*this, face_ref.entity_group, EntityOrder::face);
  if(face_ref.index >= entity_group_entry.faces().size()) {
    throw std::out_of_range("Face reference index is outside the owning face entity_group.");
  }
  return entity_group_entry.faces()[face_ref.index];
}

const Cell &Domain::cell(EntityRef cell_ref) const
{
  const auto &entity_group_entry = require_entity_group(*this, cell_ref.entity_group, EntityOrder::cell);
  if(cell_ref.index >= entity_group_entry.cells().size()) {
    throw std::out_of_range("Cell reference index is outside the owning cell entity_group.");
  }
  return entity_group_entry.cells()[cell_ref.index];
}

EntityRefRange Domain::edge_nodes(EntityRef edge_ref) const
{
  const auto &entity_group_entry = require_entity_group(*this, edge_ref.entity_group, EntityOrder::edge);
  const auto &edge_entry = edge(edge_ref);
  return make_range(entity_group_entry.node_channel(), edge_entry.node_span);
}

EntityRefRange Domain::face_nodes(EntityRef face_ref) const
{
  const auto &entity_group_entry = require_entity_group(*this, face_ref.entity_group, EntityOrder::face);
  const auto &face_entry = face(face_ref);
  return make_range(entity_group_entry.node_channel(), face_entry.node_span);
}

EntityRefRange Domain::cell_nodes(EntityRef cell_ref) const
{
  const auto &entity_group_entry = require_entity_group(*this, cell_ref.entity_group, EntityOrder::cell);
  const auto &cell_entry = cell(cell_ref);
  return make_range(entity_group_entry.node_channel(), cell_entry.node_span);
}

EntityRefRange Domain::cell_faces(EntityRef cell_ref) const
{
  const auto &entity_group_entry = require_entity_group(*this, cell_ref.entity_group, EntityOrder::cell);
  const auto &cell_entry = cell(cell_ref);
  return make_range(entity_group_entry.face_channel(), cell_entry.face_span);
}

EntityRef Domain::adjacent_face(EntityRef edge_ref, FaceSide side) const
{
  const auto &edge_entry = edge(edge_ref);
  return side == FaceSide::left ? edge_entry.left_face : edge_entry.right_face;
}

EntityRef Domain::adjacent_cell(EntityRef face_ref, FaceSide side) const
{
  const auto &face_entry = face(face_ref);
  return side == FaceSide::left ? face_entry.left_cell : face_entry.right_cell;
}

void Domain::set_face_topology_owner(
  EntityRef face_ref,
  geo::TopologyEntityId owner
)
{
  auto &entity_group_entry = entity_group_mutable(face_ref.entity_group);
  if(entity_group_entry.order() != EntityOrder::face) {
    throw std::invalid_argument("EntityGroup order mismatch for requested mesh entity operation.");
  }
  if(face_ref.index >= entity_group_entry.faces_.size()) {
    throw std::out_of_range("Face reference index is outside the owning face entity_group.");
  }
  if(geo::is_valid(owner) && owner.dimension != geo::TopologyDimension::face) {
    throw std::invalid_argument("Face topology ownership currently supports only geometry-face ids.");
  }

  entity_group_entry.face_topology_owner_channel_[face_ref.index] = owner.index;
}

geo::TopologyEntityId Domain::face_topology_owner(EntityRef face_ref) const
{
  const auto &entity_group_entry = require_entity_group(*this, face_ref.entity_group, EntityOrder::face);
  static_cast<void>(face(face_ref));

  if(face_ref.index >= entity_group_entry.face_topology_owner_channel_.size()) {
    return {};
  }

  return {
    geo::TopologyDimension::face,
    entity_group_entry.face_topology_owner_channel_[face_ref.index],
  };
}

void Domain::set_face_source_entity_tag(
  EntityRef face_ref,
  std::uint32_t source_entity_tag
)
{
  auto &entity_group_entry = entity_group_mutable(face_ref.entity_group);
  if(entity_group_entry.order() != EntityOrder::face) {
    throw std::invalid_argument("EntityGroup order mismatch for requested mesh entity operation.");
  }
  if(face_ref.index >= entity_group_entry.faces_.size()) {
    throw std::out_of_range("Face reference index is outside the owning face entity_group.");
  }

  entity_group_entry.face_source_entity_tag_channel_[face_ref.index] = source_entity_tag;
}

std::uint32_t Domain::face_source_entity_tag(EntityRef face_ref) const
{
  const auto &entity_group_entry = require_entity_group(*this, face_ref.entity_group, EntityOrder::face);
  static_cast<void>(face(face_ref));

  if(face_ref.index >= entity_group_entry.face_source_entity_tag_channel_.size()) {
    return invalid_index;
  }

  return entity_group_entry.face_source_entity_tag_channel_[face_ref.index];
}

std::size_t Domain::node_count() const noexcept
{
  std::size_t count = 0;
  for(const auto &entity_group_entry : entity_groups_) {
    if(entity_group_entry.role() == EntityGroupRole::computational) {
      count += entity_group_entry.nodes().size();
    }
  }
  return count;
}

std::size_t Domain::edge_count() const noexcept
{
  std::size_t count = 0;
  for(const auto &entity_group_entry : entity_groups_) {
    if(entity_group_entry.role() == EntityGroupRole::computational) {
      count += entity_group_entry.edges().size();
    }
  }
  return count;
}

std::size_t Domain::face_count() const noexcept
{
  std::size_t count = 0;
  for(const auto &entity_group_entry : entity_groups_) {
    if(entity_group_entry.role() == EntityGroupRole::computational) {
      count += entity_group_entry.faces().size();
    }
  }
  return count;
}

std::size_t Domain::cell_count() const noexcept
{
  std::size_t count = 0;
  for(const auto &entity_group_entry : entity_groups_) {
    if(entity_group_entry.role() == EntityGroupRole::computational) {
      count += entity_group_entry.cells().size();
    }
  }
  return count;
}

MeshSummary Domain::summary() const noexcept
{
  return {
    node_count(),
    edge_count(),
    face_count(),
    cell_count(),
    source_topology_revision_,
  };
}

DomainStatistics Domain::statistics() const noexcept
{
  DomainStatistics stats;
  stats.summary = summary();
  stats.entity_group_count = 0;

  for(const auto &entity_group_entry : entity_groups_) {
    if(entity_group_entry.role() == EntityGroupRole::computational) {
      ++stats.entity_group_count;
    }

    // Skip entity_group metric accumulation for non-computational entity_groups like geometric proxies
    if(entity_group_entry.role() != EntityGroupRole::computational) continue;

    for(const auto &edge_entry : entity_group_entry.edges()) {
      if(edge_entry.header.kind == EntityKind::edge_line) {
        ++stats.line_edge_count;
      }
      if(is_boundary_entity(edge_entry.header.flags) || !is_valid(edge_entry.right_face)) {
        ++stats.boundary_edge_count;
      }
      else {
        ++stats.interior_edge_count;
      }
    }
    for(const auto &face_entry : entity_group_entry.faces()) {
      if(face_entry.header.kind == EntityKind::face_triangle) {
        ++stats.triangle_face_count;
      }
      if(is_boundary_entity(face_entry.header.flags) || !is_valid(face_entry.right_cell)) {
        ++stats.boundary_face_count;
      }
      else {
        ++stats.interior_face_count;
      }
    }
    for(const auto &cell_entry : entity_group_entry.cells()) {
      if(cell_entry.header.kind == EntityKind::cell_tetra) {
        ++stats.tetra_cell_count;
      }
    }
  }

  return stats;
}

ElementQuality Domain::element_quality(EntityRef entity_ref) const
{
  ElementQuality report;
  report.entity = entity_ref;

  const auto &owner_entity_group = entity_group(entity_ref.entity_group);
  switch(owner_entity_group.order()) {
  case EntityOrder::face: {
    const auto &face_entry = face(entity_ref);
    report.kind = face_entry.header.kind;
    if(face_entry.header.kind == EntityKind::face_triangle) {
      return make_triangle_quality(*this, entity_ref);
    }
    return report;
  }
  case EntityOrder::cell: {
    const auto &cell_entry = cell(entity_ref);
    report.kind = cell_entry.header.kind;
    if(cell_entry.header.kind == EntityKind::cell_tetra) {
      return make_tetra_quality(*this, entity_ref);
    }
    return report;
  }
  default:
    report.kind = owner_entity_group.default_kind();
    return report;
  }
}

MeshQualityReport Domain::quality_report() const
{
  MeshQualityReport report;
  report.summary = summary();
  report.elements.reserve(face_count() + cell_count());

  std::vector<KindQualityAccumulator> accumulators;
  accumulators.reserve(2U);

  const auto record_quality = [&](const ElementQuality &element_quality) {
    report.elements.push_back(element_quality);
    if(!element_quality.supported) {
      ++report.unsupported_element_count;
      return;
    }

    ++report.supported_element_count;
    switch(element_quality.status) {
    case QualityStatus::valid:
      ++report.valid_element_count;
      break;
    case QualityStatus::degenerate:
      ++report.degenerate_element_count;
      break;
    case QualityStatus::inverted:
      ++report.inverted_element_count;
      break;
    case QualityStatus::unsupported:
    default:
      ++report.unsupported_element_count;
      return;
    }

    const auto summary_index =
      find_or_append_kind_summary(report.kinds, accumulators, element_quality.kind);
    auto &kind_summary = report.kinds[summary_index];
    auto &kind_accumulator = accumulators[summary_index];

    ++kind_summary.supported_element_count;
    switch(element_quality.status) {
    case QualityStatus::valid:
      ++kind_summary.valid_element_count;
      break;
    case QualityStatus::degenerate:
      ++kind_summary.degenerate_element_count;
      break;
    case QualityStatus::inverted:
      ++kind_summary.inverted_element_count;
      break;
    case QualityStatus::unsupported:
    default:
      break;
    }

    kind_accumulator.jacobian.observe(element_quality.jacobian);
    kind_accumulator.skewness.observe(element_quality.skewness);
    kind_accumulator.aspect_ratio.observe(element_quality.aspect_ratio);
    kind_accumulator.radius_ratio.observe(element_quality.radius_ratio);
    kind_accumulator.min_angle.observe(element_quality.min_angle);
    kind_accumulator.max_angle.observe(element_quality.max_angle);
  };

  for(const auto &entity_group_entry : entity_groups_) {
    if(entity_group_entry.role() != EntityGroupRole::computational) continue;

    if(entity_group_entry.order() == EntityOrder::face) {
      for(std::uint32_t index = 0; index < entity_group_entry.faces().size(); ++index) {
        record_quality(element_quality({entity_group_entry.id(), index}));
      }
    }
    else if(entity_group_entry.order() == EntityOrder::cell) {
      for(std::uint32_t index = 0; index < entity_group_entry.cells().size(); ++index) {
        record_quality(element_quality({entity_group_entry.id(), index}));
      }
    }
  }

  for(std::size_t index = 0; index < report.kinds.size(); ++index) {
    report.kinds[index].jacobian = accumulators[index].jacobian.summary();
    report.kinds[index].skewness = accumulators[index].skewness.summary();
    report.kinds[index].aspect_ratio = accumulators[index].aspect_ratio.summary();
    report.kinds[index].radius_ratio = accumulators[index].radius_ratio.summary();
    report.kinds[index].min_angle = accumulators[index].min_angle.summary();
    report.kinds[index].max_angle = accumulators[index].max_angle.summary();
  }

  return report;
}

EntityGroup &Domain::entity_group_mutable(EntityGroupIndex entity_group_index)
{
  if(entity_group_index >= entity_groups_.size()) {
    throw std::out_of_range("EntityGroup index is outside the domain entity_group table.");
  }
  return entity_groups_[entity_group_index];
}

std::string_view module_name() noexcept
{
  return "sqmesh::mesh";
}

Domain make_dummy_tetra_domain()
{
  Domain domain("dummy_tetra_domain");

  const auto node_entity_group = domain.create_entity_group({EntityOrder::node, "nodes"});
  const auto face_entity_group = domain.create_entity_group(
    {EntityOrder::face, "boundary_faces", invalid_index, true, EntityKind::face_triangle}
  );
  const auto cell_entity_group = domain.create_entity_group(
    {EntityOrder::cell, "volume_cells", invalid_index, false, EntityKind::cell_tetra}
  );

  domain.reserve_entity_group_storage(node_entity_group, 4U);
  domain.reserve_entity_group_storage(face_entity_group, 4U, 12U);
  domain.reserve_entity_group_storage(cell_entity_group, 1U, 4U, 4U);

  const auto n0 = domain.add_node(node_entity_group, {0.0, 0.0, 0.0});
  const auto n1 = domain.add_node(node_entity_group, {1.0, 0.0, 0.0});
  const auto n2 = domain.add_node(node_entity_group, {0.0, 1.0, 0.0});
  const auto n3 = domain.add_node(node_entity_group, {0.0, 0.0, 1.0});

  const auto f0 = domain.add_triangle_face(face_entity_group, {n0, n2, n1});
  const auto f1 = domain.add_triangle_face(face_entity_group, {n0, n1, n3});
  const auto f2 = domain.add_triangle_face(face_entity_group, {n1, n2, n3});
  const auto f3 = domain.add_triangle_face(face_entity_group, {n2, n0, n3});

  const auto c0 = domain.add_tetra_cell(cell_entity_group, {n0, n1, n2, n3}, {f0, f1, f2, f3});

  domain.set_face_cells(f0, c0);
  domain.set_face_cells(f1, c0);
  domain.set_face_cells(f2, c0);
  domain.set_face_cells(f3, c0);

  return domain;
}

} // namespace sqmesh::mesh
