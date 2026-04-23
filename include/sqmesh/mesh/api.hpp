// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#pragma once

#include "sqmesh/base/api.hpp"
#include "sqmesh/geo/api.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>
#include <variant>

namespace sqmesh::mesh {

using MeshHandle = sqmesh::Handle;
using EntityGroupIndex = std::uint32_t;
inline constexpr EntityGroupIndex invalid_index = std::numeric_limits<EntityGroupIndex>::max();

enum class EntityOrder : std::uint8_t {
  node = 0,
  edge = 1,
  face = 2,
  cell = 3,
};

enum class FaceSide : std::uint8_t {
  left = 0,
  right = 1,
};

inline constexpr std::uint32_t entity_order_mask = 0x3U;
inline constexpr std::uint32_t entity_boundary_mask = 0x4U;
inline constexpr std::uint32_t entity_interface_mask = 0x20U;
inline constexpr std::uint32_t entity_meshed_mask = 0x00100000U;
inline constexpr std::uint32_t entity_live_mask = 0x40000000U;

[[nodiscard]] constexpr std::uint16_t encode_entity_kind(
  EntityOrder order,
  std::uint8_t arity
) noexcept
{
  return static_cast<std::uint16_t>(
    (static_cast<std::uint16_t>(order) << 8U) | static_cast<std::uint16_t>(arity)
  );
}

enum class EntityKind : std::uint16_t {
  invalid = 0,
  node_point = encode_entity_kind(EntityOrder::node, 1U),
  edge_line = encode_entity_kind(EntityOrder::edge, 2U),
  face_triangle = encode_entity_kind(EntityOrder::face, 3U),
  face_quad = encode_entity_kind(EntityOrder::face, 4U),
  cell_tetra = encode_entity_kind(EntityOrder::cell, 4U),
  cell_pyramid = encode_entity_kind(EntityOrder::cell, 5U),
  cell_prism = encode_entity_kind(EntityOrder::cell, 6U),
  cell_hexa = encode_entity_kind(EntityOrder::cell, 8U),
};

[[nodiscard]] constexpr EntityOrder entity_order(EntityKind kind) noexcept
{
  return static_cast<EntityOrder>(static_cast<std::uint16_t>(kind) >> 8U);
}

[[nodiscard]] constexpr std::uint8_t entity_arity(EntityKind kind) noexcept
{
  return static_cast<std::uint8_t>(static_cast<std::uint16_t>(kind) & 0xFFU);
}

[[nodiscard]] constexpr bool is_boundary_entity(std::uint32_t flags) noexcept
{
  return (flags & entity_boundary_mask) != 0U;
}

[[nodiscard]] constexpr bool is_interface_entity(std::uint32_t flags) noexcept
{
  return (flags & entity_interface_mask) != 0U;
}

[[nodiscard]] constexpr EntityOrder entity_order(std::uint32_t flags) noexcept
{
  return static_cast<EntityOrder>(flags & entity_order_mask);
}

struct EntityRef final {
  EntityGroupIndex entity_group = invalid_index;
  std::uint32_t index = invalid_index;
};

[[nodiscard]] constexpr bool operator==(EntityRef lhs, EntityRef rhs) noexcept
{
  return lhs.entity_group == rhs.entity_group && lhs.index == rhs.index;
}

[[nodiscard]] constexpr bool operator!=(EntityRef lhs, EntityRef rhs) noexcept
{
  return !(lhs == rhs);
}

[[nodiscard]] constexpr bool is_valid(EntityRef ref) noexcept
{
  return ref.entity_group != invalid_index && ref.index != invalid_index;
}

struct ConnectivitySpan final {
  std::uint32_t offset = 0;
  std::uint32_t count = 0;
};

struct EntityHeader final {
  std::uint64_t id = 0;
  std::uint32_t flags = 0;
  std::uint32_t index = invalid_index;
  EntityKind kind = EntityKind::invalid;
  std::uint16_t reserved = 0;
  EntityGroupIndex entity_group = invalid_index;
};

struct Node final {
  EntityHeader header {};
  std::array<double, 3> coordinates {0.0, 0.0, 0.0};
};

struct Edge final {
  EntityHeader header {};
  ConnectivitySpan node_span {};
  EntityRef left_face {};
  EntityRef right_face {};
};

struct Face final {
  EntityHeader header {};
  ConnectivitySpan node_span {};
  EntityRef left_cell {};
  EntityRef right_cell {};
};

struct Cell final {
  EntityHeader header {};
  ConnectivitySpan node_span {};
  ConnectivitySpan face_span {};
};

/*
 * Reviewer-facing layout snapshot for the mesh core.
 *
 * These values intentionally stay explicit instead of being reconstructed ad
 * hoc in every smoke test or benchmark. They make entity-size drift and the
 * adjacency-relevant member offsets reviewable in one place.
 */
struct MeshCoreLayout final {
  std::size_t entity_ref_size = 0;
  std::size_t entity_ref_alignment = 0;
  std::size_t connectivity_span_size = 0;
  std::size_t connectivity_span_alignment = 0;
  std::size_t entity_header_size = 0;
  std::size_t entity_header_alignment = 0;
  std::size_t node_size = 0;
  std::size_t node_alignment = 0;
  std::size_t edge_size = 0;
  std::size_t edge_alignment = 0;
  std::size_t face_size = 0;
  std::size_t face_alignment = 0;
  std::size_t cell_size = 0;
  std::size_t cell_alignment = 0;
  std::size_t node_header_offset = 0;
  std::size_t edge_header_offset = 0;
  std::size_t edge_node_span_offset = 0;
  std::size_t edge_left_face_offset = 0;
  std::size_t edge_right_face_offset = 0;
  std::size_t face_header_offset = 0;
  std::size_t face_node_span_offset = 0;
  std::size_t face_left_cell_offset = 0;
  std::size_t face_right_cell_offset = 0;
  std::size_t cell_header_offset = 0;
  std::size_t cell_node_span_offset = 0;
  std::size_t cell_face_span_offset = 0;
};

[[nodiscard]] constexpr MeshCoreLayout mesh_core_layout() noexcept
{
  return {
    sizeof(EntityRef),
    alignof(EntityRef),
    sizeof(ConnectivitySpan),
    alignof(ConnectivitySpan),
    sizeof(EntityHeader),
    alignof(EntityHeader),
    sizeof(Node),
    alignof(Node),
    sizeof(Edge),
    alignof(Edge),
    sizeof(Face),
    alignof(Face),
    sizeof(Cell),
    alignof(Cell),
    offsetof(Node, header),
    offsetof(Edge, header),
    offsetof(Edge, node_span),
    offsetof(Edge, left_face),
    offsetof(Edge, right_face),
    offsetof(Face, header),
    offsetof(Face, node_span),
    offsetof(Face, left_cell),
    offsetof(Face, right_cell),
    offsetof(Cell, header),
    offsetof(Cell, node_span),
    offsetof(Cell, face_span),
  };
}

[[nodiscard]] constexpr MeshCoreLayout expected_mesh_core_layout() noexcept
{
  return {
    8U,
    4U,
    8U,
    4U,
    24U,
    8U,
    48U,
    8U,
    48U,
    8U,
    48U,
    8U,
    40U,
    8U,
    0U,
    0U,
    24U,
    32U,
    40U,
    0U,
    24U,
    32U,
    40U,
    0U,
    24U,
    32U,
  };
}

[[nodiscard]] constexpr bool same_mesh_core_layout(
  const MeshCoreLayout &lhs,
  const MeshCoreLayout &rhs
) noexcept
{
  return lhs.entity_ref_size == rhs.entity_ref_size &&
         lhs.entity_ref_alignment == rhs.entity_ref_alignment &&
         lhs.connectivity_span_size == rhs.connectivity_span_size &&
         lhs.connectivity_span_alignment == rhs.connectivity_span_alignment &&
         lhs.entity_header_size == rhs.entity_header_size &&
         lhs.entity_header_alignment == rhs.entity_header_alignment &&
         lhs.node_size == rhs.node_size &&
         lhs.node_alignment == rhs.node_alignment &&
         lhs.edge_size == rhs.edge_size &&
         lhs.edge_alignment == rhs.edge_alignment &&
         lhs.face_size == rhs.face_size &&
         lhs.face_alignment == rhs.face_alignment &&
         lhs.cell_size == rhs.cell_size &&
         lhs.cell_alignment == rhs.cell_alignment &&
         lhs.node_header_offset == rhs.node_header_offset &&
         lhs.edge_header_offset == rhs.edge_header_offset &&
         lhs.edge_node_span_offset == rhs.edge_node_span_offset &&
         lhs.edge_left_face_offset == rhs.edge_left_face_offset &&
         lhs.edge_right_face_offset == rhs.edge_right_face_offset &&
         lhs.face_header_offset == rhs.face_header_offset &&
         lhs.face_node_span_offset == rhs.face_node_span_offset &&
         lhs.face_left_cell_offset == rhs.face_left_cell_offset &&
         lhs.face_right_cell_offset == rhs.face_right_cell_offset &&
         lhs.cell_header_offset == rhs.cell_header_offset &&
         lhs.cell_node_span_offset == rhs.cell_node_span_offset &&
         lhs.cell_face_span_offset == rhs.cell_face_span_offset;
}

[[nodiscard]] constexpr bool matches_mesh_core_layout_baseline(
  const MeshCoreLayout &layout
) noexcept
{
  return same_mesh_core_layout(layout, expected_mesh_core_layout());
}

/*
 * `mesh_layout_benchmark` currently builds one disconnected tetrahedron per
 * benchmark cell. Each generated cell therefore owns:
 * - 4 `Node` records
 * - 4 boundary `Face` records
 * - 8 face-scoped `uint32_t` side-channel entries
 * - 1 `Cell` record
 * - 12 face-to-node side-channel `EntityRef`s
 * - 4 cell-to-face side-channel `EntityRef`s
 *
 * Keeping the derivation here makes the memory regression gate explicit and
 * easy to update if the synthetic topology ever changes.
 */
[[nodiscard]] constexpr std::size_t layout_benchmark_bytes_per_cell(
  const MeshCoreLayout &layout
) noexcept
{
  return 4U * layout.node_size + 4U * layout.face_size + layout.cell_size +
         16U * layout.entity_ref_size + 8U * sizeof(std::uint32_t);
}

[[nodiscard]] constexpr std::size_t expected_layout_benchmark_bytes_per_cell() noexcept
{
  return layout_benchmark_bytes_per_cell(expected_mesh_core_layout());
}

struct MeshSummary final {
  std::size_t node_count = 0;
  std::size_t edge_count = 0;
  std::size_t face_count = 0;
  std::size_t cell_count = 0;
  std::uint64_t source_topology_revision = 0U;
};

enum class ParameterType : std::uint8_t {
  empty = 0,
  integer = 1,
  number = 2,
  boolean = 3,
  text = 4,
};

class ParameterValue final
{
public:
  ParameterValue() = default;
  ParameterValue(std::int64_t value);
  ParameterValue(double value);
  ParameterValue(bool value);
  ParameterValue(std::string value);
  ParameterValue(const char *value);

  [[nodiscard]] ParameterType type() const noexcept;
  [[nodiscard]] bool empty() const noexcept;
  [[nodiscard]] bool is_integer() const noexcept;
  [[nodiscard]] bool is_number() const noexcept;
  [[nodiscard]] bool is_boolean() const noexcept;
  [[nodiscard]] bool is_text() const noexcept;
  [[nodiscard]] std::int64_t integer(std::int64_t fallback = 0) const noexcept;
  [[nodiscard]] double number(double fallback = 0.0) const noexcept;
  [[nodiscard]] bool boolean(bool fallback = false) const noexcept;
  [[nodiscard]] std::string_view text() const noexcept;

private:
  using Storage = std::variant<std::monostate, std::int64_t, double, bool, std::string>;

  Storage storage_ {};
};

class ParameterDictionary final
{
public:
  void set(std::string key, ParameterValue value);
  void set_integer(std::string key, std::int64_t value);
  void set_number(std::string key, double value);
  void set_boolean(std::string key, bool value);
  void set_text(std::string key, std::string value);

  [[nodiscard]] bool contains(std::string_view key) const;
  [[nodiscard]] const ParameterValue *find(std::string_view key) const noexcept;
  [[nodiscard]] bool try_get_integer(
    std::string_view key,
    std::int64_t &value
  ) const noexcept;
  [[nodiscard]] bool try_get_number(
    std::string_view key,
    double &value
  ) const noexcept;
  [[nodiscard]] bool try_get_boolean(
    std::string_view key,
    bool &value
  ) const noexcept;
  [[nodiscard]] bool try_get_text(
    std::string_view key,
    std::string_view &value
  ) const noexcept;
  [[nodiscard]] std::size_t size() const noexcept;

private:
  std::map<std::string, ParameterValue, std::less<>> entries_ {};
};

struct MeshSizeControl final {
  geo::TopologyEntityId entity {};
  double target_size = 0.0;
};

struct MeshSizeControls final {
  // When non-zero, meshing requires the local-size entity ids to come from the
  // same geometry topology revision.
  std::uint64_t topology_revision = 0;
  // Edge and face ids are currently supported. If multiple controls apply to
  // the same entity or edge neighborhood, the smallest target size wins.
  std::vector<MeshSizeControl> local_sizes {};

  void add_local_size(geo::TopologyEntityId entity, double target_size)
  {
    local_sizes.push_back({entity, target_size});
  }

  [[nodiscard]] bool empty() const noexcept
  {
    return local_sizes.empty();
  }
};

struct MeshingOptions final {
  ParameterDictionary parameters {};
  MeshSizeControls size_controls {};
};

struct AutoCfdSpacingOptions final {
  double growth_rate = 1.2;
  double feature_angle = 15.0;
  double minimum_length = 1.0;
  double maximum_length = 1000.0;
  double sharp_angle_limit = 15.0;
  double sharp_angle_length = 1.0;
  bool self_proximity = false;
  bool pid_proximity = false;
  double maximum_normals_angle = 40.0;
  double length_to_gap_ratio = 0.3;
  double proximity_minimum_length = 1.0;
};

struct AutoCfdSurfaceDefaults final {
  std::string element_type = "tria";
  double growth_rate = 1.2;
  double distortion_angle = 15.0;
  double minimum_length = 1.0;
  double maximum_length = 1000.0;
  bool self_proximity = false;
  bool pid_proximity = false;
  double maximum_normals_angle = 40.0;
  double length_to_gap_ratio = 0.3;
  double proximity_minimum_length = 1.0;
};

struct DomainStatistics final {
  MeshSummary summary {};
  std::size_t entity_group_count = 0;
  std::size_t boundary_edge_count = 0;
  std::size_t interior_edge_count = 0;
  std::size_t boundary_face_count = 0;
  std::size_t interior_face_count = 0;
  std::size_t line_edge_count = 0;
  std::size_t triangle_face_count = 0;
  std::size_t tetra_cell_count = 0;
};

enum class QualityStatus : std::uint8_t {
  unsupported = 0,
  valid = 1,
  degenerate = 2,
  inverted = 3,
};

struct QualityMetricSummary final {
  std::size_t count = 0;
  double minimum = 0.0;
  double maximum = 0.0;
  double average = 0.0;
};

/*
 * Phase 1 quality conventions:
 * - `face_triangle`
 *   - `jacobian` is an unsigned corner-Jacobian quality:
 *     `min_i(sin(theta_i)) / sin(60 deg)`.
 *     This stays in `[0, 1]` for non-degenerate triangles. SQMesh does not
 *     report standalone triangle inversion because the mesh core does not store
 *     a reference surface normal.
 *   - `skewness` is equiangular skewness:
 *     `max((theta_max - 60 deg) / 120 deg, (60 deg - theta_min) / 60 deg)`.
 *     An equilateral triangle is `0`, and a degenerate triangle tends to `1`.
 *   - `aspect_ratio` is the normalized longest-edge to shortest-altitude ratio:
 *     `(sqrt(3) / 2) * longest_edge / shortest_altitude`.
 *     An equilateral triangle is `1`.
 *   - `radius_ratio` is the normalized inscribed-to-circumscribed circle ratio:
 *     `2 * r / R = 16 * area^2 / (perimeter * l01 * l12 * l20)`.
 *     An equilateral triangle is `1`, and a degenerate triangle is `0`.
 *   - `min_angle` / `max_angle` are the minimum and maximum interior corner
 *     angles over the 3 triangle corners, reported in degrees.
 *     Degenerate triangles leave both angle fields as `NaN`.
 * - `cell_tetra`
 *   - `jacobian` is the minimum signed scaled corner Jacobian over the 4
 *     vertices:
 *     `sqrt(2) * det(a_i, b_i, c_i) / (|a_i| |b_i| |c_i|)`.
 *     A regular tetrahedron is `1`, a degenerate tetrahedron is `0`, and an
 *     inverted tetrahedron is negative.
 *   - `skewness` is equihedral skewness from the 6 interior dihedral angles:
 *     `max((phi_max - phi_ideal) / (180 deg - phi_ideal),
 *          (phi_ideal - phi_min) / phi_ideal)`,
 *     where `phi_ideal = acos(1 / 3)`.
 *   - `aspect_ratio` is the normalized longest-edge to shortest-altitude ratio:
 *     `sqrt(2 / 3) * longest_edge / shortest_altitude`.
 *     A regular tetrahedron is `1`.
 *   - `radius_ratio` is the normalized inscribed-to-circumscribed sphere ratio:
 *     `3 * r / R`, where `r = 3 * volume / surface_area` and `R` is computed
 *     from the circumsphere through the 4 tetrahedron vertices.
 *     A regular tetrahedron is `1`. This measure is geometric and remains
 *     unsigned; tetrahedron inversion is still reported separately through the
 *     signed `jacobian` and `status` fields.
 *   - `min_angle` / `max_angle` are the minimum and maximum interior dihedral
 *     angles over the 6 tetrahedron edges, reported in degrees.
 *     Degenerate tetrahedra leave both angle fields as `NaN`.
 */
struct ElementQuality final {
  EntityRef entity {};
  EntityKind kind = EntityKind::invalid;
  QualityStatus status = QualityStatus::unsupported;
  bool supported = false;
  bool degenerate = false;
  bool inverted = false;
  bool jacobian_is_signed = false;
  double jacobian = std::numeric_limits<double>::quiet_NaN();
  double skewness = std::numeric_limits<double>::quiet_NaN();
  double aspect_ratio = std::numeric_limits<double>::quiet_NaN();
  double radius_ratio = std::numeric_limits<double>::quiet_NaN();
  double min_angle = std::numeric_limits<double>::quiet_NaN();
  double max_angle = std::numeric_limits<double>::quiet_NaN();
};

struct KindQualitySummary final {
  EntityKind kind = EntityKind::invalid;
  std::size_t supported_element_count = 0;
  std::size_t valid_element_count = 0;
  std::size_t degenerate_element_count = 0;
  std::size_t inverted_element_count = 0;
  QualityMetricSummary jacobian {};
  QualityMetricSummary skewness {};
  QualityMetricSummary aspect_ratio {};
  QualityMetricSummary radius_ratio {};
  QualityMetricSummary min_angle {};
  QualityMetricSummary max_angle {};
};

struct MeshQualityReport final {
  MeshSummary summary {};
  std::size_t supported_element_count = 0;
  std::size_t unsupported_element_count = 0;
  std::size_t valid_element_count = 0;
  std::size_t degenerate_element_count = 0;
  std::size_t inverted_element_count = 0;
  std::vector<KindQualitySummary> kinds {};
  std::vector<ElementQuality> elements {};
};

/*
 * Reviewer-facing mesh-core entity_group semantics:
 * - `node`: point storage only.
 * - `interior`: non-boundary edge/face topology inside one region.
 * - `boundary`: an edge/face entity_group that bounds a single adjacent region.
 * - `interface`: an edge/face entity_group shared by two distinct adjacent regions.
 * - `region`: a cell entity_group representing one solver-facing volume region.
 *
 * `primary_region_zone_id` / `secondary_region_zone_id` expose the adjacent
 * region-zone ids when they are known from import or construction time.
 */
enum class EntityGroupSemantic : std::uint8_t {
  unspecified = 0,
  node = 1,
  interior = 2,
  boundary = 3,
  interface = 4,
  region = 5,
};

enum class EntityGroupRole : std::uint8_t {
  computational = 0,
  geometric_proxy = 1,
  annotation = 2,
};

enum class EntityGroupImportFormat : std::uint8_t {
  none = 0,
  cgns = 1,
  nastran = 2,
};

enum class NastranEntityGroupSourceCard : std::uint8_t {
  unspecified = 0,
  cbar = 1,
  cbeam = 2,
  crod = 3,
  ctria3 = 4,
  cquad4 = 5,
  ctetra = 6,
  mixed = 7,
};

enum class NastranPropertyCard : std::uint8_t {
  unspecified = 0,
  prod = 1,
  pshell = 2,
  psolid = 3,
};

enum class NastranMaterialCard : std::uint8_t {
  unspecified = 0,
  mat1 = 1,
};

struct CgnsEntityGroupImportInfo final {
  std::uint32_t base_index = invalid_index;
  std::string base_name {};
  std::uint32_t zone_index = invalid_index;
  std::string zone_name {};
  std::string local_name {};
  std::int32_t bc_type_value = -1;
};

struct NastranEntityGroupImportInfo final {
  NastranEntityGroupSourceCard source_card = NastranEntityGroupSourceCard::unspecified;
  NastranPropertyCard property_card = NastranPropertyCard::unspecified;
  NastranMaterialCard material_card = NastranMaterialCard::unspecified;
  std::uint32_t material_id = invalid_index;
  double rod_area = std::numeric_limits<double>::quiet_NaN();
  double shell_thickness = std::numeric_limits<double>::quiet_NaN();
  double youngs_modulus = std::numeric_limits<double>::quiet_NaN();
  double shear_modulus = std::numeric_limits<double>::quiet_NaN();
  double poisson_ratio = std::numeric_limits<double>::quiet_NaN();
};

struct EntityGroupImportInfo final {
  EntityGroupImportFormat format = EntityGroupImportFormat::none;
  CgnsEntityGroupImportInfo cgns {};
  NastranEntityGroupImportInfo nastran {};
};

struct EntityGroupDefinition final {
  EntityOrder order = EntityOrder::node;
  std::string name {};
  std::uint32_t zone_id = invalid_index;
  bool boundary = false;
  EntityKind default_kind = EntityKind::invalid;
  EntityGroupSemantic semantic = EntityGroupSemantic::unspecified;
  EntityGroupRole role = EntityGroupRole::computational;
  std::uint32_t primary_region_zone_id = invalid_index;
  std::uint32_t secondary_region_zone_id = invalid_index;
  // Preserves an imported source entity/block tag when the format exposes one.
  std::uint32_t source_entity_tag = invalid_index;
  EntityGroupImportInfo import_info {};
};

struct EntityGroupInfo final {
  EntityGroupIndex id = invalid_index;
  std::uint32_t zone_id = invalid_index;
  std::uint32_t source_entity_tag = invalid_index;
  EntityOrder order = EntityOrder::node;
  bool boundary = false;
  EntityKind default_kind = EntityKind::invalid;
  EntityGroupSemantic semantic = EntityGroupSemantic::unspecified;
  EntityGroupRole role = EntityGroupRole::computational;
  std::uint32_t primary_region_zone_id = invalid_index;
  std::uint32_t secondary_region_zone_id = invalid_index;
  std::string name {};
  EntityGroupImportInfo import_info {};
};

/*
 * Phase 1.1 metadata split:
 * - `EntityGroupImportInfo` carries narrow imported entity_group provenance when the
 *   source data is honestly entity_group-scoped.
 * - per-edge source topology ownership lives in a typed side channel on the
 *   owning edge entity_group and is accessed through `Domain::edge_topology_owner()`.
 * - per-face source topology ownership lives in a typed side channel on the
 *   owning face entity_group and is accessed through `Domain::face_topology_owner()`.
 * - per-face imported source-entity tags live in a second typed side channel
 *   and are accessed through `Domain::face_source_entity_tag()`.
 */
struct EntityRefRange final {
  const EntityRef *data = nullptr;
  std::size_t size = 0;

  [[nodiscard]] const EntityRef *begin() const noexcept { return data; }
  [[nodiscard]] const EntityRef *end() const noexcept { return data + size; }
  [[nodiscard]] bool empty() const noexcept { return size == 0U; }
  [[nodiscard]] const EntityRef &operator[](std::size_t offset) const noexcept
  {
    return data[offset];
  }
};

/*
 * Mesh-core layout:
 * - `Domain -> EntityGroup` ownership with compact shared entity prefixes, pooled
 *   contiguous storage, and explicit `Face -> Left/Right Cell` adjacency. This
 *   favors batched mesh construction and hot adjacency queries over
 *   pointer-linked OO objects with virtual interfaces.
 * - One `Domain` owns homogeneous `EntityGroup` zones; each entity_group owns contiguous
 *   entity pools plus entity_group-local side channels. Entities are indexed by a
 *   stable `EntityGroupIndex` so they stay trivially movable inside vectors and
 *   future pool growth does not invalidate entity headers.
 */
class EntityGroup final
{
public:
  EntityGroup() = default;
  explicit EntityGroup(EntityGroupInfo info);

  [[nodiscard]] const EntityGroupInfo &info() const noexcept;
  [[nodiscard]] const EntityGroupImportInfo &import_info() const noexcept;
  [[nodiscard]] EntityGroupIndex id() const noexcept;
  [[nodiscard]] std::uint32_t zone_id() const noexcept;
  [[nodiscard]] std::uint32_t source_entity_tag() const noexcept;
  [[nodiscard]] EntityOrder order() const noexcept;
  [[nodiscard]] bool is_boundary() const noexcept;
  [[nodiscard]] EntityGroupSemantic semantic() const noexcept;
  [[nodiscard]] EntityGroupRole role() const noexcept;
  [[nodiscard]] bool is_interface() const noexcept;
  [[nodiscard]] std::uint32_t primary_region_zone_id() const noexcept;
  [[nodiscard]] std::uint32_t secondary_region_zone_id() const noexcept;
  [[nodiscard]] EntityKind default_kind() const noexcept;
  [[nodiscard]] std::string_view name() const noexcept;
  [[nodiscard]] std::size_t entity_count() const noexcept;

  [[nodiscard]] const std::vector<Node> &nodes() const noexcept;
  [[nodiscard]] const std::vector<Edge> &edges() const noexcept;
  [[nodiscard]] const std::vector<Face> &faces() const noexcept;
  [[nodiscard]] const std::vector<Cell> &cells() const noexcept;
  [[nodiscard]] const std::vector<EntityRef> &node_channel() const noexcept;
  [[nodiscard]] const std::vector<EntityRef> &face_channel() const noexcept;

private:
  friend class Domain;
  void reserve(
    std::size_t entity_count,
    std::size_t node_channel_count,
    std::size_t face_channel_count
  );

  EntityGroupInfo info_ {};
  std::vector<Node> nodes_ {};
  std::vector<Edge> edges_ {};
  std::vector<Face> faces_ {};
  std::vector<Cell> cells_ {};
  std::vector<EntityRef> node_channel_ {};
  std::vector<EntityRef> face_channel_ {};
  std::vector<std::uint32_t> edge_topology_owner_channel_ {};
  std::vector<std::uint32_t> face_topology_owner_channel_ {};
  std::vector<std::uint32_t> face_source_entity_tag_channel_ {};
};

class Domain final
{
public:
  explicit Domain(std::string name = {});

  [[nodiscard]] std::string_view name() const noexcept;
  [[nodiscard]] const std::vector<EntityGroup> &entity_groups() const noexcept;
  [[nodiscard]] const EntityGroup &entity_group(EntityGroupIndex entity_group_index) const;
  [[nodiscard]] std::size_t entity_group_count() const noexcept;
  [[nodiscard]] std::size_t entity_group_count(EntityOrder order) const noexcept;
  [[nodiscard]] std::size_t entity_group_count(EntityGroupRole role) const noexcept;

  void remove_entity_groups_with_prefix(std::string_view prefix);

  [[nodiscard]] EntityGroupIndex create_entity_group(EntityGroupDefinition definition);
  void reserve_entity_group_storage(
    EntityGroupIndex entity_group_index,
    std::size_t entity_count,
    std::size_t node_channel_count = 0,
    std::size_t face_channel_count = 0
  );

  [[nodiscard]] EntityRef add_node(
    EntityGroupIndex entity_group_index,
    const std::array<double, 3> &coordinates,
    std::uint32_t extra_flags = 0
  );
  [[nodiscard]] EntityRef add_edge(
    EntityGroupIndex entity_group_index,
    const std::array<EntityRef, 2> &nodes,
    std::uint32_t extra_flags = 0
  );
  [[nodiscard]] EntityRef add_triangle_face(
    EntityGroupIndex entity_group_index,
    const std::array<EntityRef, 3> &nodes,
    std::uint32_t extra_flags = 0
  );
  [[nodiscard]] EntityRef add_quad_face(
    EntityGroupIndex entity_group_index,
    const std::array<EntityRef, 4> &nodes,
    std::uint32_t extra_flags = 0
  );
  [[nodiscard]] EntityRef add_tetra_cell(
    EntityGroupIndex entity_group_index,
    const std::array<EntityRef, 4> &nodes,
    const std::array<EntityRef, 4> &faces,
    std::uint32_t extra_flags = 0
  );
  [[nodiscard]] EntityRef add_pyramid_cell(
    EntityGroupIndex entity_group_index,
    const std::array<EntityRef, 5> &nodes,
    const std::array<EntityRef, 5> &faces,
    std::uint32_t extra_flags = 0
  );
  [[nodiscard]] EntityRef add_prism_cell(
    EntityGroupIndex entity_group_index,
    const std::array<EntityRef, 6> &nodes,
    const std::array<EntityRef, 5> &faces,
    std::uint32_t extra_flags = 0
  );
  [[nodiscard]] EntityRef add_hexa_cell(
    EntityGroupIndex entity_group_index,
    const std::array<EntityRef, 8> &nodes,
    const std::array<EntityRef, 6> &faces,
    std::uint32_t extra_flags = 0
  );

  void set_face_cells(
    EntityRef face_ref,
    EntityRef left_cell,
    EntityRef right_cell = {}
  );
  void set_edge_faces(
    EntityRef edge_ref,
    EntityRef left_face,
    EntityRef right_face = {}
  );

  [[nodiscard]] const Node &node(EntityRef node_ref) const;
  [[nodiscard]] const Edge &edge(EntityRef edge_ref) const;
  [[nodiscard]] const Face &face(EntityRef face_ref) const;
  [[nodiscard]] const Cell &cell(EntityRef cell_ref) const;

  [[nodiscard]] EntityRefRange edge_nodes(EntityRef edge_ref) const;
  [[nodiscard]] EntityRefRange face_nodes(EntityRef face_ref) const;
  [[nodiscard]] EntityRefRange cell_nodes(EntityRef cell_ref) const;
  [[nodiscard]] EntityRefRange cell_faces(EntityRef cell_ref) const;
  [[nodiscard]] EntityRef adjacent_face(EntityRef edge_ref, FaceSide side) const;
  [[nodiscard]] EntityRef adjacent_cell(EntityRef face_ref, FaceSide side) const;
  void set_source_topology_revision(std::uint64_t topology_revision) noexcept
  {
    source_topology_revision_ = topology_revision;
  }
  [[nodiscard]] std::uint64_t source_topology_revision() const noexcept
  {
    return source_topology_revision_;
  }
  // Preserves the first narrow per-edge ownership channel for source geometry
  // edge ids without overloading entity_group metadata. Compare the returned ids
  // against the mesh summary's `source_topology_revision` before reusing them
  // after geometry-topology updates.
  void set_edge_topology_owner(
    EntityRef edge_ref,
    geo::TopologyEntityId owner
  );
  [[nodiscard]] geo::TopologyEntityId edge_topology_owner(EntityRef edge_ref) const;
  // Preserves the first narrow per-face ownership channel for source geometry
  // face ids without overloading entity_group metadata. Compare the returned ids
  // against the mesh summary's `source_topology_revision` before reusing them
  // after geometry-topology updates.
  void set_face_topology_owner(
    EntityRef face_ref,
    geo::TopologyEntityId owner
  );
  [[nodiscard]] geo::TopologyEntityId face_topology_owner(EntityRef face_ref) const;
  void set_face_source_entity_tag(
    EntityRef face_ref,
    std::uint32_t source_entity_tag
  );
  [[nodiscard]] std::uint32_t face_source_entity_tag(EntityRef face_ref) const;

  // Mesh summaries and aggregate statistics intentionally count only
  // `EntityGroupRole::computational` entities. Proxy-only or annotation-only entity_groups
  // remain visible through `domain_snapshot()` and `entity_groups()`.
  [[nodiscard]] std::size_t node_count() const noexcept;
  [[nodiscard]] std::size_t edge_count() const noexcept;
  [[nodiscard]] std::size_t face_count() const noexcept;
  [[nodiscard]] std::size_t cell_count() const noexcept;
  [[nodiscard]] MeshSummary summary() const noexcept;
  [[nodiscard]] DomainStatistics statistics() const noexcept;
  [[nodiscard]] ElementQuality element_quality(EntityRef entity_ref) const;
  [[nodiscard]] MeshQualityReport quality_report() const;

private:
  [[nodiscard]] EntityGroup &entity_group_mutable(EntityGroupIndex entity_group_index);

  std::string name_ {};
  std::vector<EntityGroup> entity_groups_ {};
  std::uint64_t next_entity_id_ = 1;
  std::uint64_t source_topology_revision_ = 0U;
};

/*
 * Phase 1 mesh-IO subsets:
 * - MSH import supports Gmsh 2.2 ASCII/binary and Gmsh 4.1 ASCII/binary files
 *   with finite node coordinates plus line (1), triangle (2), and tetrahedron
 *   (4) elements.
 *   `PhysicalNames` are mapped onto entity_group `zone_id` and `name` when present;
 *   on Gmsh 4.1 files, `$Entities` and `$PartitionedEntities` physical tags
 *   are preserved when available and otherwise fall back to per-entity entity_group
 *   labels. Distinct Gmsh geometrical entity tags remain visible on
 *   `EntityGroupInfo::source_entity_tag` and keep separate entity blocks distinct
 *   even when they share one physical group. Gmsh 4.1 import also accepts
 *   `$Periodic`, `$GhostElements`, and `$Parametrizations` sections around
 *   that supported topology subset; periodic/ghost/partition ownership and
 *   parametrization payloads are syntax-checked then ignored, while 1D/2D/3D
 *   parametric node coordinates in `$Nodes` blocks are accepted and discarded.
 *   Binary Gmsh 4.1 import accepts either 4-byte or 8-byte `MeshFormat`
 *   data-size fields for the supported tag/count records. When tetrahedral cell
 *   groups on opposite sides of a
 *   shared triangle resolve to distinct region entity_groups, SQMesh stores that face
 *   entity_group explicitly as `EntityGroupSemantic::interface` instead of collapsing it
 *   into a generic interior/boundary label; synthesized boundary/interface face
 *   entity_group names are qualified by adjacent region zone ids when the source does
 *   not provide them. MSH export writes only the core `$PhysicalNames`,
 *   `$Entities`, `$Nodes`, and `$Elements` sections for the same supported
 *   line/triangle/tetra element families, selected explicitly with
 *   `MshExportOptions::format_version`; when a entity_group carries a preserved
 *   `source_entity_tag`, export reuses it as the Gmsh geometrical entity tag
 *   instead of collapsing back onto the physical tag. Export does not
 *   currently re-emit partitioned, periodic, ghost, or parametrization
 *   metadata.
 * - OBJ import accepts surface meshes with `v`, polyline `l`, and simple
 *   planar polygonal `f` records whose vertex references are pairwise distinct
 *   and contain no three consecutive collinear vertices, plus `g`/`o`
 *   grouping labels. Supported polygons are triangulated on import through a
 *   bounded ear-clipping path, including simple concave polygons inside that
 *   non-collinear subset. Repeated-vertex, self-intersecting, non-planar, and
 *   otherwise degenerate polygon faces fail cleanly with `unsupported`.
 *   Positive and negative face-token references are supported on the bounded
 *   `v`, `v/vt`, `v//vn`, and `v/vt/vn` forms. `vt`, `vn`, `vp`, `mtllib`,
 *   and `usemtl` payloads are accepted and kept on the internal
 *   parser/import-description path, but current `Domain`/`EntityGroup` storage does
 *   not retain that per-corner/material fidelity.
 *   OBJ export stays on the bounded surface subset: point-node `v`, line `l`,
 *   and triangle `f` records, plus optional `g`/`o` labels. Export does not
 *   re-emit parser-side `vt`/`vn`/`vp` or material-library/material-use
 *   metadata. Group labels are preserved as entity_group names for the supported
 *   mesh subset.
 * - CGNS import/export is optional behind `SQMESH_ENABLE_CGNS`.
 *   CGNS import accepts the current bounded unstructured subset only:
 *   one or more bases that all share one cell dimension and one physical
 *   coordinate dimension; one or more unstructured zones per base; canonical
 *   `CoordinateX/Y/Z` arrays in each zone; homogeneous `BAR_2`, `TRI_3`,
 *   `QUAD_4`, and `TETRA_4` sections; and `MIXED` sections carrying only that
 *   same element-family subset. `QUAD_4` surface elements are fan-triangulated
 *   on import to match the current triangle-based face storage. Section names
 *   map onto entity_group names when available; when more than one zone or more than
 *   one base is imported, the flattened `Domain` qualifies section/boundary
 *   entity_group names as `z<zone-index>_<zone-name>__<name>` on the single-base
 *   multi-zone path and
 *   `b<base-index>_<base-name>__z<zone-index>_<zone-name>__<name>` on the
 *   multi-base path to preserve provenance and avoid cross-base entity_group merges,
 *   while the unqualified local name plus base/zone provenance remain
 *   available on `EntityGroupImportInfo::cgns`. `BC_t` boundary metadata is
 *   represented for the supported edge/face boundary subset
 *   (`PointList` / `PointRange` with `EdgeCenter` plus `FaceCenter` or
 *   `CellCenter` as appropriate for the zone dimension, and the same bounded
 *   edge/face subset through legacy `ElementList` / `ElementRange`);
 *   boundary-condition names take precedence over section-name fallback for the
 *   addressed boundary elements. SQMesh-written `SQMeshZoneId` descriptors map
 *   back onto entity_group `zone_id`, and import otherwise falls back to the section
 *   index. When a supported edge/face entity_group is named from a `BC_t` node,
 *   `EntityGroupImportInfo::cgns.bc_type_value` preserves the raw `BCType_t`
 *   enumerator as an integer for lossless passthrough.
 *   CGNS export still defaults to one unstructured base with one unstructured
 *   zone containing `BAR_2`, `TRI_3`, and `TETRA_4` sections plus `BC_t`
 *   records for explicit edge entity_groups and boundary face entity_groups. That default
 *   path rewrites any imported multi-zone or multi-base domain as one written
 *   zone: it keeps the current entity_group names as section/BC labels subject to
 *   CGNS's 32-character label limit and uniqueness suffixing, but it does not
 *   preserve CGNS zone/base hierarchy or `GridConnectivity_t`.
 *   When `CgnsExportOptions::write_multi_zone` is enabled, export can instead
 *   preserve multiple unstructured zones inside one CGNS base on the same
 *   bounded surface/tetra subset if the mesh entity_groups retain CGNS import
 *   metadata or the importer-style `z<zone-index>_<zone-name>__<name>`
 *   qualification. Within that one-base preserved path, export can also
 *   reconstruct two-way `GridConnectivity_t` records for shared-node interface
 *   face entity_groups on volume zones and shared-node interface edge entity_groups on
 *   surface zones. Imported multi-base provenance is preserved for review and
 *   clean rejection, but current export does not claim to reconstruct more
 *   than one CGNS base.
 *   Structured zones, parent-data sections, broader element families,
 *   non-canonical coordinate arrays on the supported 2D/3D subset, preserved
 *   multi-base export, and wider interface shapes including non-shared-node
 *   interfaces still fail cleanly with `unsupported`. Export re-emits the
 *   preserved raw `BCType_t` value when available and otherwise falls back to
 *   `BCGeneral`.
 * - NASTRAN BDF/DAT import supports a bounded mainstream bulk-data subset:
 *   single-line free-field, single-line 8-character small-field, and canonical
 *   `*`-continued 16-character large-field `GRID`, `CBAR`, `CBEAM`, `CROD`,
 *   `CTRIA3`, `CQUAD4`, and first-order `CTETRA` cards, plus the same syntax
 *   envelope for `MAT1`, `PROD`, `PSHELL`, and `PSOLID` metadata cards.
 *   SQMesh preserves a bounded typed subset on `EntityGroupImportInfo::nastran`:
 *   `MAT1` keeps `MID` plus bounded `E/G/NU`, `PROD` keeps `PID/MID/A`,
 *   `PSHELL` keeps `PID/MID1/T`, and `PSOLID` keeps `PID/MID`. `GRID`
 *   requires blank/zero coordinate-system and superelement fields;
 *   `CBAR`/`CBEAM` orientation payload and `CTRIA3`/`CQUAD4` shell extras are
 *   ignored on import.
 *   `CQUAD4` faces are fan-triangulated to match the current triangle-face
 *   storage. Imported `PID` values map onto entity_group `zone_id`; supported
 *   line/shell/volume card identity is also preserved in `EntityGroupImportInfo`
 *   where that card family remains honestly entity_group-scoped. When deterministic
 *   `line_pid_<pid>` entity_grouping legitimately combines imported `CBAR`/`CBEAM`/
 *   `CROD` records on one line entity_group, `EntityGroupImportInfo::nastran.source_card`
 *   is recorded as `mixed` instead of pretending one family owns the whole
 *   entity_group. `CQUAD4` shells keep a distinct `cquad4_pid_<pid>` face-entity_group name
 *   for reviewer-visible determinism, while the other supported families keep
 *   deterministic `line_pid_<pid>`, `surface_pid_<pid>`, and
 *   `volume_pid_<pid>` naming.
 *   Arbitrary continuation identifiers, unsupported topology cards,
 *   non-default `GRID` coordinate systems, and higher-order volume elements
 *   fail cleanly with `unsupported`.
 * - NASTRAN export writes the supported SQMesh line/triangle/tetra topology as
 *   free-field, 8-character small-field, or canonical `*`-continued
 *   16-character large-field bulk data using `GRID`, `CROD`, `CTRIA3`,
 *   `CQUAD4`, and first-order `CTETRA` cards. When NASTRAN-imported entity_group
 *   metadata remains self-consistent, the exporter also re-emits the bounded
 *   preserved metadata subset via `MAT1`, `PROD`, `PSHELL`, and `PSOLID`.
 *   Conflicting metadata for one `MID` or `PID` fails truthfully instead of
 *   guessing. EntityGroups imported from `CQUAD4` are re-emitted as `CQUAD4` only
 *   while they still contain the importer-style paired triangles
 *   `(G1,G2,G3)` plus `(G1,G3,G4)`; the exporter prefers
 *   `EntityGroupImportInfo::nastran` for that distinction and otherwise falls back
 *   to the deterministic entity_group-name convention. All other face entity_groups
 *   continue to export as `CTRIA3`.
 */
struct MshImportOptions final {
  bool read_physical_names = true;
};

enum class MshFormatVersion : std::uint8_t {
  gmsh22_ascii = 0,
  gmsh22_binary = 1,
  gmsh41_ascii = 2,
  gmsh41_binary = 3,
};

struct MshExportOptions final {
  bool write_physical_names = true;
  MshFormatVersion format_version = MshFormatVersion::gmsh22_ascii;
};

struct ObjImportOptions final {
  // Preserves bounded `o`/`g` labels as the imported domain and entity_group names.
  bool read_groups = true;
};

struct ObjExportOptions final {
  // Writes only `o`/`g` labels plus `v`, `l`, and triangle `f` records on the
  // supported surface subset; parser-side OBJ attribute/material metadata are
  // not re-emitted.
  bool write_groups = true;
  bool write_object_name = true;
  // Controls which entity_group role is exported. Defaults to `computational` so
  // that proxy/annotation entity_groups in a shared Domain are excluded.
  EntityGroupRole export_role = EntityGroupRole::computational;
};

struct CgnsImportOptions final {
  bool read_section_names = true;
};

struct CgnsExportOptions final {
  std::string base_name {};
  std::string zone_name {};
  // Keeps the existing one-base/one-zone export by default, even for imported
  // multi-zone or multi-base domains. When enabled, export can preserve
  // multiple unstructured CGNS zones only inside one CGNS base on the current
  // bounded subset: every non-node entity_group must carry the importer-style
  // `z<zone-index>_<zone-name>__<name>` qualification or equivalent typed CGNS
  // import metadata, all zones must stay in one surface-or-volume dimension,
  // imported multi-base provenance still rejects cleanly, and interface
  // connectivity is reconstructed only for shared-node interface face entity_groups
  // on volume zones or shared-node interface edge entity_groups on surface zones.
  bool write_multi_zone = false;
};

struct NastranImportOptions final {
  bool allow_free_fields = true;
  bool allow_small_fixed_fields = true;
  bool allow_large_fixed_fields = true;
};

enum class NastranFieldFormat : std::uint8_t {
  free_field = 0,
  small_fixed = 1,
  large_fixed = 2,
};

struct NastranExportOptions final {
  NastranFieldFormat field_format = NastranFieldFormat::free_field;
  bool write_begin_bulk = true;
  bool write_enddata = true;
};

[[nodiscard]] std::string_view module_name() noexcept;
[[nodiscard]] Domain make_dummy_tetra_domain();
// Built-in meshing algorithms keep explicit subset boundaries:
// - `Dummy Mesher`: compatibility path that always returns one tetrahedron.
// - `Auto CFD Surface Mesher`: bounded in-repo CFD surface path that validates
//   `minimum_length`, `maximum_length`, `distortion_angle`, `growth_rate`, and
//   optional boolean `proximity`, optional split `self_proximity` /
//   `pid_proximity`, optional proximity controls (`maximum_normals_angle`,
//   `length_to_gap_ratio`, `proximity_minimum_length`), while also accepting
//   optional edge-spacing overrides (`spacing_growth_rate`, `spacing_feature_angle`,
//   `spacing_minimum_length`, `spacing_maximum_length`,
//   `spacing_sharp_angle_limit`, `spacing_sharp_angle_length`,
//   `spacing_proximity`, `spacing_self_proximity`, `spacing_pid_proximity`,
//   `spacing_maximum_normals_angle`, `spacing_length_to_gap_ratio`,
//   `spacing_proximity_minimum_length` plus `feature_angle` /
//   `sharp_angle_*` aliases). The mesher freezes supported hard/non-tangent
//   0D anchors, builds native curvature/proximity sizing sources, injects them
//   into a bounded native global background grid with growth-limited diffusion
//   while preserving the current query contract. The delivered proximity path
//   now uses a spatially indexed proxy-triangle BVH/AABB-style query, keeps
//   explicit shell/part buckets for Self vs. PID semantics, and adds
//   edge-owned proximity sources when face-interior samples would miss a
//   boundary-local narrow gap. The mesher then carries a durable
//   size-field-driven 1D boundary discretization state for later loop
//   recovery/shared-edge reuse, runs a bounded UV-space
//   frontal-Delaunay-style triangle path with local constrained repair, and
//   applies final duplicate/non-manifold/reprojection plus
//   source-face-containment/orientation and near-degenerate screening before
//   runtime-owned `Domain` injection. The current public meshing-success
//   subset is simple planar faces, one-hole trimmed planar faces, the
//   supported unwrap subset of seam-bearing single-loop faces, and the
//   bounded smoke-verified low-distortion subset of non-seam singly trimmed
//   cylindrical and conical patch faces. The current
//   public delivered triangle gate is measured rather than CFD-red-line
//   complete: `min_angle >= 2.6 deg`, `max_angle <= 170 deg`,
//   `aspect_ratio <= 20.5`, `radius_ratio >= 0.018`, `skewness <= 0.96`.
//   Internal review-only diagnostics also keep a broader development
//   guardrail (`1.5 / 175 / 40 / 0.005 / 0.98`) so tests can distinguish
//   severe meshing-quality breakage from delivered-quality screening.
//   Post-Part-2 re-measurement only justified tightening min-angle and
//   aspect-ratio; the supported seam-bearing plus bounded assembly evidence
//   still bottoms out or tops out near `168.75 / 0.0191683 / 0.955466` for
//   max-angle / radius-ratio / skewness, so those public limits remain
//   conservative.
//   Bounded assembly-style smoke/benchmark evidence now covers both the
//   original six-face OCC compound composed from the delivered
//   simple/trimmed/seam-bearing subset and a repeated 12-face / 24-face
//   evidence path built from the same supported cell. The current per-face
//   hard stops remain explicit:
//   `kMaximumInsertedVerticesPerFace = 100000`,
//   `kMaximumFrontIterationsPerFace = 500000`. On the measured repeated path,
//   per-face usage still stayed far below those caps (`0 / 100000` inserted
//   vertices and `64 / 500000` front iterations), but the larger 12-face / 24-face
//   candidate domains now drift to roughly `min_angle = 2.48506` and
//   `aspect_ratio = 20.9289`, which the delivered final screen rejects as
//   quality-gate failures. The current larger-scale ceiling therefore appears
//   before per-face cap exhaustion and remains evidence for bounded mixed-part
//   correctness only, not a claim of assembly-scale readiness.
//   The widened non-planar trimmed claim remains intentionally narrow: it
//   does not imply broad arbitrary trimmed-surface success outside those
//   smoke-verified cylindrical/conical patch subsets.
//   Degenerate seam trims and discrete STL public meshing remain
//   `unsupported`.
// - `Simple Surface Mesher`: geometry-driven triangle mesher for simple faces
//   plus explicit broader subsets: multi-hole trimmed faces and seam-bearing
//   single-loop faces without degenerate trim edges.
// - `Tetrahedral Volume Mesher`: geometry-driven Tet4 mesher for single-region
//   OCC solids bounded by straight edges and planar faces on the currently
//   tested subset: orthogonal boxes, triangular prisms, mixed-face pyramids,
//   and the bounded star-shaped L-prism case covered by the smoke suite. When
//   `optimization_pass = "tetra_relocate"` is requested through the parameter
//   dictionary, this mesher also applies a bounded Tet4 quality-improvement
//   pass on that same supported subset by relocating the generated interior
//   kernel node and planar face-center nodes only; boundary edge/vertex samples
//   remain fixed on the source geometry.
[[nodiscard]] bool is_algorithm_registered(std::string_view algorithm_name) noexcept;
[[nodiscard]] bool cgns_io_available() noexcept;
[[nodiscard]] base::StatusCode import_msh(
  std::string_view path,
  MeshHandle &mesh_handle,
  const MshImportOptions &options = {},
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode export_msh(
  MeshHandle mesh_handle,
  std::string_view path,
  const MshExportOptions &options = {},
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
// Imports the bounded OBJ subset documented above, including simple planar
// concave polygons with pairwise-distinct, non-collinear consecutive face
// vertices plus parser-level preservation of supported face-token and
// material-library/material-use records.
[[nodiscard]] base::StatusCode import_obj(
  std::string_view path,
  MeshHandle &mesh_handle,
  const ObjImportOptions &options = {},
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
// Exports only the bounded surface subset: point-node `v`, line `l`, and
// triangle `f` records with optional `o`/`g` labels. Parser-side OBJ
// attribute/material metadata are not re-emitted.
[[nodiscard]] base::StatusCode export_obj(
  MeshHandle mesh_handle,
  std::string_view path,
  const ObjExportOptions &options = {},
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
// Imports the bounded unstructured CGNS subset documented above, flattening
// supported multi-zone and same-dimension multi-base files into one `Domain`.
[[nodiscard]] base::StatusCode import_cgns(
  std::string_view path,
  MeshHandle &mesh_handle,
  const CgnsImportOptions &options = {},
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
// Exports the bounded CGNS subset documented above: one written zone by
// default, or one-base preserved multi-zone output when `write_multi_zone` is
// enabled on the supported subset.
[[nodiscard]] base::StatusCode export_cgns(
  MeshHandle mesh_handle,
  std::string_view path,
  const CgnsExportOptions &options = {},
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode import_nastran(
  std::string_view path,
  MeshHandle &mesh_handle,
  const NastranImportOptions &options = {},
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode export_nastran(
  MeshHandle mesh_handle,
  std::string_view path,
  const NastranExportOptions &options = {},
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
void apply_auto_cfd_spacing(
  MeshingOptions &options,
  const AutoCfdSpacingOptions &spacing
) noexcept;
void apply_auto_cfd_surface_defaults(
  MeshingOptions &options,
  const AutoCfdSurfaceDefaults &defaults
) noexcept;
[[nodiscard]] base::StatusCode create_surface_mesh(
  geo::ModelHandle model_handle,
  std::string_view algorithm_name,
  const MeshingOptions &options,
  MeshHandle &mesh_handle,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode create_surface_mesh(
  geo::ModelHandle model_handle,
  std::string_view algorithm_name,
  const ParameterDictionary &parameters,
  MeshHandle &mesh_handle,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode create_volume_mesh(
  geo::ModelHandle model_handle,
  std::string_view algorithm_name,
  const MeshingOptions &options,
  MeshHandle &mesh_handle,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode create_volume_mesh(
  geo::ModelHandle model_handle,
  std::string_view algorithm_name,
  const ParameterDictionary &parameters,
  MeshHandle &mesh_handle,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode mesh_summary(
  MeshHandle mesh_handle,
  MeshSummary &summary,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode mesh_quality_report(
  MeshHandle mesh_handle,
  MeshQualityReport &report,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode domain_snapshot(
  MeshHandle mesh_handle,
  Domain &domain,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode nodes_count(
  MeshHandle mesh_handle,
  std::size_t &count,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;
[[nodiscard]] base::StatusCode cells_count(
  MeshHandle mesh_handle,
  std::size_t &count,
  base::ContextHandle context_handle = sqmesh::invalid_handle
) noexcept;

static_assert(sizeof(EntityRef) == 8U, "EntityRef must stay compact for channel storage.");
static_assert(
  sizeof(ConnectivitySpan) == 8U,
  "ConnectivitySpan must stay compact for entity_group-local side arrays."
);
static_assert(sizeof(EntityHeader) == 24U, "EntityHeader should match the shared compact prefix.");
static_assert(offsetof(Node, header) == 0U, "Node must start with EntityHeader.");
static_assert(offsetof(Edge, header) == 0U, "Edge must start with EntityHeader.");
static_assert(offsetof(Face, header) == 0U, "Face must start with EntityHeader.");
static_assert(offsetof(Cell, header) == 0U, "Cell must start with EntityHeader.");
static_assert(sizeof(Node) == 48U, "Node size drift breaks compact batch assumptions.");
static_assert(sizeof(Edge) == 48U, "Edge size drift breaks adjacency cache assumptions.");
static_assert(sizeof(Face) == 48U, "Face size drift breaks adjacency cache assumptions.");
static_assert(sizeof(Cell) == 40U, "Cell size drift breaks compact connectivity assumptions.");
static_assert(std::is_standard_layout<EntityHeader>::value, "EntityHeader must stay standard-layout.");
static_assert(std::is_standard_layout<Node>::value, "Node must stay standard-layout.");
static_assert(std::is_standard_layout<Edge>::value, "Edge must stay standard-layout.");
static_assert(std::is_standard_layout<Face>::value, "Face must stay standard-layout.");
static_assert(std::is_standard_layout<Cell>::value, "Cell must stay standard-layout.");

} // namespace sqmesh::mesh
