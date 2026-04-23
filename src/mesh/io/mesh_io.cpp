// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "sqmesh/mesh/api.hpp"

#include "mesh_io_test_hook.hpp"

#include "core/runtime_registry.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#if defined(sqmesh_HAVE_CGNS)
#include <cgnslib.h>
#endif

namespace sqmesh::mesh {
namespace {

class MeshIoFailure final : public std::runtime_error
{
public:
  MeshIoFailure(base::StatusCode code, std::string message)
    : std::runtime_error(std::move(message)), code_(code)
  {
  }

  [[nodiscard]] base::StatusCode code() const noexcept
  {
    return code_;
  }

private:
  base::StatusCode code_ = base::StatusCode::internal_error;
};

[[noreturn]] void fail(base::StatusCode code, std::string message)
{
  throw MeshIoFailure(code, std::move(message));
}

void validate_non_empty_path(
  std::string_view path,
  std::string_view message
);

[[nodiscard]] std::map<std::uint64_t, std::size_t> build_node_indices(
  const Domain &domain,
  EntityGroupRole role_filter = EntityGroupRole::computational
);

[[nodiscard]] std::uint64_t pack_entity_ref(EntityRef ref) noexcept;

void validate_domain_for_msh_export(const Domain &domain);

[[nodiscard]] std::string format_nastran_small_field_integer(
  std::uint64_t value,
  std::string_view context
);

[[nodiscard]] std::string format_nastran_small_field_real(
  double value,
  std::string_view context
);

[[nodiscard]] std::string format_nastran_large_field_integer(
  std::uint64_t value,
  std::string_view context
);

[[nodiscard]] std::string format_nastran_large_field_real(
  double value,
  std::string_view context
);

void write_nastran_small_field_record(
  std::ostream &output,
  std::string_view card,
  const std::vector<std::string> &fields
);

void write_nastran_large_field_record(
  std::ostream &output,
  std::string_view card,
  const std::vector<std::string> &fields
);

[[nodiscard]] bool starts_with(
  std::string_view value,
  std::string_view prefix
) noexcept
{
  return value.size() >= prefix.size() &&
         value.compare(0U, prefix.size(), prefix) == 0;
}

[[nodiscard]] std::string trim_copy(std::string_view value)
{
  std::size_t first = 0U;
  while(first < value.size() &&
        std::isspace(static_cast<unsigned char>(value[first])) != 0) {
    ++first;
  }

  std::size_t last = value.size();
  while(last > first &&
        std::isspace(static_cast<unsigned char>(value[last - 1U])) != 0) {
    --last;
  }

  return std::string(value.substr(first, last - first));
}

[[nodiscard]] std::string rtrim_copy(std::string_view value)
{
  std::size_t last = value.size();
  while(last > 0U &&
        (value[last - 1U] == '\r' ||
         std::isspace(static_cast<unsigned char>(value[last - 1U])) != 0)) {
    --last;
  }
  return std::string(value.substr(0U, last));
}

[[nodiscard]] std::string uppercase_copy(std::string_view value)
{
  std::string upper(value);
  for(auto &ch : upper) {
    ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
  }
  return upper;
}

[[nodiscard]] std::string trim_upper_copy(std::string_view value)
{
  return uppercase_copy(trim_copy(value));
}

template <typename T>
[[nodiscard]] std::string format_value(const T &value)
{
  std::ostringstream stream;
  stream << value;
  return stream.str();
}

[[nodiscard]] bool is_finite_point(const std::array<double, 3> &point) noexcept
{
  return std::isfinite(point[0]) &&
         std::isfinite(point[1]) &&
         std::isfinite(point[2]);
}

template <typename NodeId>
void validate_msh_node_coordinates(
  const std::array<double, 3> &point,
  const NodeId &node_id
)
{
  if(is_finite_point(point)) {
    return;
  }

  fail(
    base::StatusCode::io_error,
    "MSH import found a non-finite coordinate on node " + format_value(node_id) + '.'
  );
}

[[nodiscard]] std::string msh_element_family_label(long long element_type)
{
  switch(element_type) {
  case 1:
    return "line";
  case 2:
    return "triangle";
  case 3:
    return "quadrangle";
  case 4:
    return "tetrahedron";
  case 5:
    return "hexahedron";
  case 6:
    return "prism";
  case 7:
    return "pyramid";
  case 8:
    return "second-order line";
  case 9:
    return "second-order triangle";
  case 10:
    return "second-order quadrangle";
  case 11:
    return "second-order tetrahedron";
  case 15:
    return "point";
  default:
    return {};
  }
}

[[noreturn]] void fail_unsupported_msh_element_type(long long element_type)
{
  const auto label = msh_element_family_label(element_type);
  const auto type_text =
    label.empty()
      ? "element type (" + std::to_string(element_type) + ')'
      : label + " (" + std::to_string(element_type) + ')';

  fail(
    base::StatusCode::unsupported,
    "MSH import encountered unsupported " + type_text +
      "; the current bounded subset accepts only line (1), triangle (2), and tetrahedron (4) elements."
  );
}

template <typename T>
void swap_msh_binary_value(T &value) noexcept
{
  auto *bytes = reinterpret_cast<unsigned char *>(&value);
  std::reverse(bytes, bytes + sizeof(T));
}

template <typename T>
void read_msh_binary_value(
  std::istream &input,
  T &value,
  bool swap_bytes,
  std::string_view context
)
{
  input.read(reinterpret_cast<char *>(&value), static_cast<std::streamsize>(sizeof(T)));
  if(!input) {
    fail(
      base::StatusCode::io_error,
      "MSH import found an incomplete " + std::string(context) + '.'
    );
  }

  if(swap_bytes) {
    swap_msh_binary_value(value);
  }
}

template <typename T>
void read_msh_binary_values(
  std::istream &input,
  T *values,
  std::size_t count,
  bool swap_bytes,
  std::string_view context
)
{
  if(count == 0U) {
    return;
  }

  input.read(
    reinterpret_cast<char *>(values),
    static_cast<std::streamsize>(sizeof(T) * count)
  );
  if(!input) {
    fail(
      base::StatusCode::io_error,
      "MSH import found an incomplete " + std::string(context) + '.'
    );
  }

  if(swap_bytes) {
    for(std::size_t index = 0U; index < count; ++index) {
      swap_msh_binary_value(values[index]);
    }
  }
}

template <typename T>
void write_msh_binary_value(
  std::ostream &output,
  const T &value
)
{
  output.write(reinterpret_cast<const char *>(&value), static_cast<std::streamsize>(sizeof(T)));
  if(!output) {
    fail(base::StatusCode::io_error, "MSH export could not write the output file.");
  }
}

template <typename T>
void write_msh_binary_values(
  std::ostream &output,
  const T *values,
  std::size_t count
)
{
  if(count == 0U) {
    return;
  }

  output.write(
    reinterpret_cast<const char *>(values),
    static_cast<std::streamsize>(sizeof(T) * count)
  );
  if(!output) {
    fail(base::StatusCode::io_error, "MSH export could not write the output file.");
  }
}

[[nodiscard]] int checked_msh_binary_int(
  std::uint64_t value,
  std::string_view context
)
{
  if(value > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
    fail(base::StatusCode::unsupported, std::string(context));
  }
  return static_cast<int>(value);
}

[[nodiscard]] std::size_t checked_msh_binary_size(
  std::uint64_t value,
  std::string_view context
)
{
  if(value > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
    fail(base::StatusCode::unsupported, std::string(context));
  }
  return static_cast<std::size_t>(value);
}

[[nodiscard]] long long checked_msh_binary_tag(
  std::uint64_t value,
  std::string_view context
)
{
  if(value > static_cast<std::uint64_t>(std::numeric_limits<long long>::max())) {
    fail(base::StatusCode::unsupported, std::string(context));
  }
  return static_cast<long long>(value);
}

[[nodiscard]] bool msh41_binary_data_size_supported(int data_size) noexcept
{
  return data_size == 4 || data_size == 8;
}

[[nodiscard]] std::size_t msh41_parametric_coordinate_count(
  int entity_dim,
  int parametric
) noexcept
{
  if(parametric == 0 || entity_dim <= 0) {
    return 0U;
  }

  if(entity_dim >= 3) {
    return 3U;
  }

  return static_cast<std::size_t>(entity_dim);
}

void read_msh41_binary_size_value(
  std::istream &input,
  std::uint64_t &value,
  bool swap_bytes,
  int data_size,
  std::string_view context
)
{
  if(data_size == 4) {
    std::uint32_t raw_value = 0U;
    read_msh_binary_value(input, raw_value, swap_bytes, context);
    value = raw_value;
    return;
  }

  if(data_size == 8) {
    std::uint64_t raw_value = 0U;
    read_msh_binary_value(input, raw_value, swap_bytes, context);
    value = raw_value;
    return;
  }

  fail(
    base::StatusCode::unsupported,
    "MSH import currently supports binary Gmsh 4.1 only with 4-byte or 8-byte MeshFormat data sizes."
  );
}

void read_msh41_binary_size_values(
  std::istream &input,
  std::uint64_t *values,
  std::size_t count,
  bool swap_bytes,
  int data_size,
  std::string_view context
)
{
  if(count == 0U) {
    return;
  }

  if(data_size == 4) {
    std::vector<std::uint32_t> raw_values(count, 0U);
    read_msh_binary_values(input, raw_values.data(), raw_values.size(), swap_bytes, context);
    for(std::size_t index = 0U; index < count; ++index) {
      values[index] = raw_values[index];
    }
    return;
  }

  if(data_size == 8) {
    std::vector<std::uint64_t> raw_values(count, 0U);
    read_msh_binary_values(input, raw_values.data(), raw_values.size(), swap_bytes, context);
    for(std::size_t index = 0U; index < count; ++index) {
      values[index] = raw_values[index];
    }
    return;
  }

  fail(
    base::StatusCode::unsupported,
    "MSH import currently supports binary Gmsh 4.1 only with 4-byte or 8-byte MeshFormat data sizes."
  );
}

[[nodiscard]] std::string default_domain_name(std::string_view path)
{
  const auto stem = std::filesystem::path(std::string(path)).stem().string();
  return stem.empty() ? std::string("imported_mesh") : stem;
}

[[nodiscard]] std::string default_entity_group_name(
  EntityOrder order,
  bool boundary,
  EntityGroupSemantic semantic = EntityGroupSemantic::unspecified,
  std::uint32_t primary_region_zone_id = invalid_index,
  std::uint32_t secondary_region_zone_id = invalid_index
) noexcept
{
  const auto normalized_semantic =
    semantic == EntityGroupSemantic::unspecified
      ? (order == EntityOrder::node
           ? EntityGroupSemantic::node
           : (order == EntityOrder::cell
                ? EntityGroupSemantic::region
                : (boundary ? EntityGroupSemantic::boundary : EntityGroupSemantic::interior)))
      : semantic;

  if(normalized_semantic == EntityGroupSemantic::interface &&
     primary_region_zone_id != invalid_index &&
     secondary_region_zone_id != invalid_index) {
    auto first_zone_id = primary_region_zone_id;
    auto second_zone_id = secondary_region_zone_id;
    if(second_zone_id < first_zone_id) {
      std::swap(first_zone_id, second_zone_id);
    }
    switch(order) {
    case EntityOrder::edge:
      return "z" + std::to_string(first_zone_id) + "__z" +
             std::to_string(second_zone_id) + "__interface_edges";
    case EntityOrder::face:
      return "z" + std::to_string(first_zone_id) + "__z" +
             std::to_string(second_zone_id) + "__interface_faces";
    default:
      break;
    }
  }

  if(primary_region_zone_id != invalid_index) {
    switch(normalized_semantic) {
    case EntityGroupSemantic::boundary:
      switch(order) {
      case EntityOrder::edge:
        return "z" + std::to_string(primary_region_zone_id) + "__boundary_edges";
      case EntityOrder::face:
        return "z" + std::to_string(primary_region_zone_id) + "__boundary_faces";
      default:
        break;
      }
      break;
    case EntityGroupSemantic::region:
      if(order == EntityOrder::cell) {
        return "z" + std::to_string(primary_region_zone_id) + "__volume_cells";
      }
      break;
    default:
      break;
    }
  }

  switch(order) {
  case EntityOrder::node:
    return "nodes";
  case EntityOrder::edge:
    return boundary ? "boundary_edges" : "edges";
  case EntityOrder::face:
    return boundary ? "boundary_faces" : "interior_faces";
  case EntityOrder::cell:
    return "volume_cells";
  default:
    return "mesh_entities";
  }
}

[[nodiscard]] EntityKind default_kind(EntityOrder order)
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
    fail(
      base::StatusCode::unsupported,
      "Unsupported mesh entity order in mesh IO."
    );
  }
}

[[nodiscard]] NastranEntityGroupSourceCard nastran_entity_group_source_card(
  std::string_view card
) noexcept
{
  if(card == "CBAR") {
    return NastranEntityGroupSourceCard::cbar;
  }
  if(card == "CBEAM") {
    return NastranEntityGroupSourceCard::cbeam;
  }
  if(card == "CROD") {
    return NastranEntityGroupSourceCard::crod;
  }
  if(card == "CTRIA3") {
    return NastranEntityGroupSourceCard::ctria3;
  }
  if(card == "CQUAD4") {
    return NastranEntityGroupSourceCard::cquad4;
  }
  if(card == "CTETRA") {
    return NastranEntityGroupSourceCard::ctetra;
  }
  return NastranEntityGroupSourceCard::unspecified;
}

[[nodiscard]] NastranEntityGroupSourceCard combine_nastran_entity_group_source_card(
  NastranEntityGroupSourceCard current,
  NastranEntityGroupSourceCard incoming
) noexcept
{
  if(current == NastranEntityGroupSourceCard::unspecified) {
    return incoming;
  }
  if(incoming == NastranEntityGroupSourceCard::unspecified || current == incoming) {
    return current;
  }
  return NastranEntityGroupSourceCard::mixed;
}

struct NastranMaterialInfo final {
  NastranMaterialCard card = NastranMaterialCard::unspecified;
  std::uint32_t material_id = invalid_index;
  double youngs_modulus = std::numeric_limits<double>::quiet_NaN();
  double shear_modulus = std::numeric_limits<double>::quiet_NaN();
  double poisson_ratio = std::numeric_limits<double>::quiet_NaN();
};

struct NastranPropertyInfo final {
  NastranPropertyCard card = NastranPropertyCard::unspecified;
  std::uint32_t material_id = invalid_index;
  double rod_area = std::numeric_limits<double>::quiet_NaN();
  double shell_thickness = std::numeric_limits<double>::quiet_NaN();
};

[[nodiscard]] EntityGroupImportInfo make_nastran_entity_group_import_info(
  std::string_view card,
  const NastranPropertyInfo *property = nullptr,
  const NastranMaterialInfo *material = nullptr
)
{
  EntityGroupImportInfo import_info;
  import_info.format = EntityGroupImportFormat::nastran;
  import_info.nastran.source_card = nastran_entity_group_source_card(card);
  if(property != nullptr) {
    import_info.nastran.property_card = property->card;
    import_info.nastran.material_id = property->material_id;
    import_info.nastran.rod_area = property->rod_area;
    import_info.nastran.shell_thickness = property->shell_thickness;
  }
  if(material != nullptr) {
    import_info.nastran.material_card = material->card;
    import_info.nastran.material_id = material->material_id;
    import_info.nastran.youngs_modulus = material->youngs_modulus;
    import_info.nastran.shear_modulus = material->shear_modulus;
    import_info.nastran.poisson_ratio = material->poisson_ratio;
  }
  return import_info;
}

[[nodiscard]] bool same_nastran_optional_real(double lhs, double rhs) noexcept
{
  if(!std::isfinite(lhs) && !std::isfinite(rhs)) {
    return true;
  }
  return lhs == rhs;
}

struct GroupInfo final {
  std::uint32_t zone_id = invalid_index;
  std::uint32_t source_entity_tag = invalid_index;
  std::string name {};
  EntityGroupImportInfo import_info {};
};

struct ObjOptionalIndex final {
  bool present = false;
  std::size_t value = 0U;
};

struct ObjCoordinate final {
  std::array<double, 3> values {0.0, 0.0, 0.0};
  std::uint8_t component_count = 0U;
};

struct ObjMaterialLibrary final {
  std::vector<std::string> entries {};
};

struct ObjFaceVertex final {
  std::size_t vertex_index = 0U;
  ObjOptionalIndex texture_index {};
  ObjOptionalIndex normal_index {};
};

struct ObjImportData final {
  std::vector<ObjCoordinate> texture_coordinates {};
  std::vector<std::array<double, 3>> normals {};
  std::vector<ObjCoordinate> parameter_vertices {};
  std::vector<ObjMaterialLibrary> material_libraries {};
};

struct ImportedEdge final {
  std::array<std::size_t, 2> nodes {};
  GroupInfo group {};
  std::string material_name {};
};

struct ImportedFace final {
  std::array<std::size_t, 3> nodes {};
  std::array<ObjFaceVertex, 3> obj_face_vertices {};
  GroupInfo group {};
  std::string material_name {};
  std::uint32_t source_entity_tag = invalid_index;
};

struct ImportedCell final {
  std::array<std::size_t, 4> nodes {};
  GroupInfo group {};
};

struct ImportedMeshDescription final {
  std::string name {};
  std::vector<std::array<double, 3>> nodes {};
  std::vector<ImportedEdge> edges {};
  std::vector<ImportedFace> faces {};
  std::vector<ImportedCell> cells {};
  ObjImportData obj {};
};

void append_imported_face(
  std::vector<ImportedFace> &faces,
  const std::array<std::size_t, 3> &nodes,
  GroupInfo group = {},
  std::uint32_t source_entity_tag = invalid_index
)
{
  ImportedFace face;
  face.nodes = nodes;
  face.group = std::move(group);
  face.source_entity_tag = source_entity_tag;
  faces.push_back(std::move(face));
}

struct NastranImportedGroupKey final {
  EntityOrder order = EntityOrder::node;
  std::uint32_t zone_id = invalid_index;
  std::string name {};

  [[nodiscard]] bool operator<(const NastranImportedGroupKey &other) const noexcept
  {
    if(order != other.order) {
      return order < other.order;
    }
    if(zone_id != other.zone_id) {
      return zone_id < other.zone_id;
    }
    return name < other.name;
  }
};

void finalize_nastran_entity_group_import_info(ImportedMeshDescription &description)
{
  std::map<NastranImportedGroupKey, NastranEntityGroupSourceCard> source_cards;

  const auto accumulate_group = [&](EntityOrder order, const GroupInfo &group) {
    if(group.import_info.format != EntityGroupImportFormat::nastran) {
      return;
    }

    NastranImportedGroupKey key;
    key.order = order;
    key.zone_id = group.zone_id;
    key.name = group.name;

    auto &combined = source_cards[key];
    combined = combine_nastran_entity_group_source_card(
      combined,
      group.import_info.nastran.source_card
    );
  };

  for(const auto &edge : description.edges) {
    accumulate_group(EntityOrder::edge, edge.group);
  }
  for(const auto &face : description.faces) {
    accumulate_group(EntityOrder::face, face.group);
  }
  for(const auto &cell : description.cells) {
    accumulate_group(EntityOrder::cell, cell.group);
  }

  const auto apply_group = [&](EntityOrder order, GroupInfo &group) {
    if(group.import_info.format != EntityGroupImportFormat::nastran) {
      return;
    }

    NastranImportedGroupKey key;
    key.order = order;
    key.zone_id = group.zone_id;
    key.name = group.name;
    const auto card_it = source_cards.find(key);
    if(card_it != source_cards.end()) {
      group.import_info.nastran.source_card = card_it->second;
    }
  };

  for(auto &edge : description.edges) {
    apply_group(EntityOrder::edge, edge.group);
  }
  for(auto &face : description.faces) {
    apply_group(EntityOrder::face, face.group);
  }
  for(auto &cell : description.cells) {
    apply_group(EntityOrder::cell, cell.group);
  }
}

struct ObjProjectedPoint final {
  double u = 0.0;
  double v = 0.0;
};

[[nodiscard]] std::array<double, 3> subtract_points(
  const std::array<double, 3> &lhs,
  const std::array<double, 3> &rhs
) noexcept
{
  return {
    lhs[0] - rhs[0],
    lhs[1] - rhs[1],
    lhs[2] - rhs[2],
  };
}

[[nodiscard]] double dot_product(
  const std::array<double, 3> &lhs,
  const std::array<double, 3> &rhs
) noexcept
{
  return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2];
}

[[nodiscard]] double vector_length(const std::array<double, 3> &value) noexcept
{
  return std::sqrt(dot_product(value, value));
}

[[nodiscard]] ObjProjectedPoint project_obj_point(
  const std::array<double, 3> &point,
  std::size_t dropped_axis
) noexcept
{
  switch(dropped_axis) {
  case 0U:
    return {point[1], point[2]};
  case 1U:
    return {point[0], point[2]};
  default:
    return {point[0], point[1]};
  }
}

[[nodiscard]] double obj_cross_2d(
  const ObjProjectedPoint &a,
  const ObjProjectedPoint &b,
  const ObjProjectedPoint &c
) noexcept
{
  const auto ab_u = b.u - a.u;
  const auto ab_v = b.v - a.v;
  const auto bc_u = c.u - b.u;
  const auto bc_v = c.v - b.v;
  return ab_u * bc_v - ab_v * bc_u;
}

[[nodiscard]] bool obj_point_in_segment_bounds_2d(
  const ObjProjectedPoint &point,
  const ObjProjectedPoint &a,
  const ObjProjectedPoint &b,
  double coordinate_tolerance
) noexcept
{
  return point.u >= std::min(a.u, b.u) - coordinate_tolerance &&
         point.u <= std::max(a.u, b.u) + coordinate_tolerance &&
         point.v >= std::min(a.v, b.v) - coordinate_tolerance &&
         point.v <= std::max(a.v, b.v) + coordinate_tolerance;
}

[[nodiscard]] bool obj_point_on_segment_2d(
  const ObjProjectedPoint &point,
  const ObjProjectedPoint &a,
  const ObjProjectedPoint &b,
  double area_tolerance,
  double coordinate_tolerance
) noexcept
{
  const auto cross = obj_cross_2d(a, b, point);
  return std::isfinite(cross) &&
         std::abs(cross) <= area_tolerance &&
         obj_point_in_segment_bounds_2d(point, a, b, coordinate_tolerance);
}

[[nodiscard]] int obj_orient_2d(
  const ObjProjectedPoint &a,
  const ObjProjectedPoint &b,
  const ObjProjectedPoint &c,
  double area_tolerance
) noexcept
{
  const auto cross = obj_cross_2d(a, b, c);
  if(!std::isfinite(cross) || std::abs(cross) <= area_tolerance) {
    return 0;
  }
  return cross > 0.0 ? 1 : -1;
}

[[nodiscard]] bool obj_segments_intersect_2d(
  const ObjProjectedPoint &a0,
  const ObjProjectedPoint &a1,
  const ObjProjectedPoint &b0,
  const ObjProjectedPoint &b1,
  double area_tolerance,
  double coordinate_tolerance
) noexcept
{
  const auto a_orientation_0 = obj_orient_2d(a0, a1, b0, area_tolerance);
  const auto a_orientation_1 = obj_orient_2d(a0, a1, b1, area_tolerance);
  const auto b_orientation_0 = obj_orient_2d(b0, b1, a0, area_tolerance);
  const auto b_orientation_1 = obj_orient_2d(b0, b1, a1, area_tolerance);

  if(a_orientation_0 == 0 &&
     obj_point_on_segment_2d(b0, a0, a1, area_tolerance, coordinate_tolerance)) {
    return true;
  }
  if(a_orientation_1 == 0 &&
     obj_point_on_segment_2d(b1, a0, a1, area_tolerance, coordinate_tolerance)) {
    return true;
  }
  if(b_orientation_0 == 0 &&
     obj_point_on_segment_2d(a0, b0, b1, area_tolerance, coordinate_tolerance)) {
    return true;
  }
  if(b_orientation_1 == 0 &&
     obj_point_on_segment_2d(a1, b0, b1, area_tolerance, coordinate_tolerance)) {
    return true;
  }

  return a_orientation_0 != a_orientation_1 &&
         b_orientation_0 != b_orientation_1;
}

[[nodiscard]] bool obj_point_in_triangle_2d(
  const ObjProjectedPoint &point,
  const ObjProjectedPoint &a,
  const ObjProjectedPoint &b,
  const ObjProjectedPoint &c,
  int orientation_sign,
  double area_tolerance
) noexcept
{
  const auto ab = orientation_sign * obj_cross_2d(a, b, point);
  const auto bc = orientation_sign * obj_cross_2d(b, c, point);
  const auto ca = orientation_sign * obj_cross_2d(c, a, point);

  return std::isfinite(ab) && std::isfinite(bc) && std::isfinite(ca) &&
         ab >= -area_tolerance &&
         bc >= -area_tolerance &&
         ca >= -area_tolerance;
}

struct ObjPreparedPolygon final {
  std::vector<std::size_t> refs {};
  std::vector<ObjProjectedPoint> projected {};
  double area_tolerance = 0.0;
  double coordinate_tolerance = 0.0;
  int orientation_sign = 1;
};

[[nodiscard]] ObjPreparedPolygon prepare_obj_face_polygon(
  std::vector<std::size_t> refs,
  const std::vector<std::array<double, 3>> &nodes,
  std::size_t line_number
)
{
  if(refs.size() >= 4U && refs.front() == refs.back()) {
    refs.pop_back();
  }

  if(refs.size() < 3U) {
    fail(
      base::StatusCode::unsupported,
      "OBJ import requires at least three distinct vertices in an `f` record."
    );
  }

  std::set<std::size_t> unique_refs(refs.begin(), refs.end());
  if(unique_refs.size() != refs.size()) {
    fail(
      base::StatusCode::unsupported,
      "OBJ import currently supports only simple polygon faces without repeated vertex references; line " +
        std::to_string(line_number) + " is not supported."
    );
  }

  ObjPreparedPolygon polygon;
  polygon.refs = std::move(refs);

  std::array<double, 3> minimum = nodes[polygon.refs.front()];
  std::array<double, 3> maximum = minimum;
  std::array<double, 3> normal = {0.0, 0.0, 0.0};
  for(std::size_t index = 0U; index < polygon.refs.size(); ++index) {
    const auto &current = nodes[polygon.refs[index]];
    const auto &next = nodes[polygon.refs[(index + 1U) % polygon.refs.size()]];
    normal[0] += (current[1] - next[1]) * (current[2] + next[2]);
    normal[1] += (current[2] - next[2]) * (current[0] + next[0]);
    normal[2] += (current[0] - next[0]) * (current[1] + next[1]);

    for(std::size_t axis = 0U; axis < 3U; ++axis) {
      minimum[axis] = std::min(minimum[axis], current[axis]);
      maximum[axis] = std::max(maximum[axis], current[axis]);
    }
  }

  const auto max_extent = std::max({
    maximum[0] - minimum[0],
    maximum[1] - minimum[1],
    maximum[2] - minimum[2],
  });
  const auto plane_tolerance = 1.0e-9 * std::max(1.0, max_extent);
  const auto normal_length = vector_length(normal);
  if(!std::isfinite(normal_length) || normal_length <= plane_tolerance) {
    fail(
      base::StatusCode::unsupported,
      "OBJ import currently supports only non-degenerate polygon faces; line " +
        std::to_string(line_number) + " is degenerate."
    );
  }

  const std::array<double, 3> unit_normal = {
    normal[0] / normal_length,
    normal[1] / normal_length,
    normal[2] / normal_length,
  };
  const auto &origin = nodes[polygon.refs.front()];
  for(std::size_t index = 1U; index < polygon.refs.size(); ++index) {
    const auto offset = subtract_points(nodes[polygon.refs[index]], origin);
    if(std::abs(dot_product(offset, unit_normal)) > plane_tolerance) {
      fail(
        base::StatusCode::unsupported,
        "OBJ import currently supports only planar polygon faces; line " +
          std::to_string(line_number) + " is not coplanar."
      );
    }
  }

  std::size_t dropped_axis = 2U;
  if(std::abs(unit_normal[0]) >= std::abs(unit_normal[1]) &&
     std::abs(unit_normal[0]) >= std::abs(unit_normal[2])) {
    dropped_axis = 0U;
  }
  else if(std::abs(unit_normal[1]) >= std::abs(unit_normal[2])) {
    dropped_axis = 1U;
  }

  polygon.projected.reserve(polygon.refs.size());
  for(const auto ref : polygon.refs) {
    polygon.projected.push_back(project_obj_point(nodes[ref], dropped_axis));
  }

  polygon.area_tolerance = 1.0e-12 * std::max(1.0, max_extent * max_extent);
  polygon.coordinate_tolerance = 1.0e-9 * std::max(1.0, max_extent);
  for(std::size_t index = 0U; index < polygon.projected.size(); ++index) {
    const auto cross =
      obj_cross_2d(
        polygon.projected[index],
        polygon.projected[(index + 1U) % polygon.projected.size()],
        polygon.projected[(index + 2U) % polygon.projected.size()]
      );
    if(!std::isfinite(cross) || std::abs(cross) <= polygon.area_tolerance) {
      fail(
        base::StatusCode::unsupported,
        "OBJ import currently supports only simple planar polygon faces with non-collinear vertices; line " +
          std::to_string(line_number) + " is degenerate."
      );
    }
  }

  for(std::size_t first = 0U; first < polygon.projected.size(); ++first) {
    const auto first_next = (first + 1U) % polygon.projected.size();
    for(std::size_t second = first + 1U; second < polygon.projected.size(); ++second) {
      const auto second_next = (second + 1U) % polygon.projected.size();
      if(first == second || first_next == second || second_next == first) {
        continue;
      }

      if(obj_segments_intersect_2d(
           polygon.projected[first],
           polygon.projected[first_next],
           polygon.projected[second],
           polygon.projected[second_next],
           polygon.area_tolerance,
           polygon.coordinate_tolerance
         )) {
        fail(
          base::StatusCode::unsupported,
          "OBJ import currently supports only simple planar polygon faces; line " +
            std::to_string(line_number) + " is self-intersecting."
        );
      }
    }
  }

  double signed_area = 0.0;
  for(std::size_t index = 0U; index < polygon.projected.size(); ++index) {
    const auto &current = polygon.projected[index];
    const auto &next = polygon.projected[(index + 1U) % polygon.projected.size()];
    signed_area += current.u * next.v - next.u * current.v;
  }
  if(!std::isfinite(signed_area) || std::abs(signed_area) <= polygon.area_tolerance) {
    fail(
      base::StatusCode::unsupported,
      "OBJ import currently supports only non-degenerate polygon faces; line " +
        std::to_string(line_number) + " collapses in projection."
    );
  }

  polygon.orientation_sign = signed_area > 0.0 ? 1 : -1;

  return polygon;
}

[[nodiscard]] std::vector<std::array<std::size_t, 3>> triangulate_obj_face_refs(
  std::vector<std::size_t> refs,
  const std::vector<std::array<double, 3>> &nodes,
  std::size_t line_number
)
{
  auto polygon = prepare_obj_face_polygon(std::move(refs), nodes, line_number);

  std::vector<std::size_t> active_indices(polygon.refs.size(), 0U);
  std::iota(active_indices.begin(), active_indices.end(), 0U);

  std::vector<std::array<std::size_t, 3>> triangles;
  triangles.reserve(polygon.refs.size() - 2U);

  while(active_indices.size() > 3U) {
    bool found_ear = false;
    for(std::size_t position = 0U; position < active_indices.size(); ++position) {
      const auto previous_index =
        active_indices[(position + active_indices.size() - 1U) % active_indices.size()];
      const auto current_index = active_indices[position];
      const auto next_index = active_indices[(position + 1U) % active_indices.size()];

      const auto ear_cross =
        polygon.orientation_sign * obj_cross_2d(
          polygon.projected[previous_index],
          polygon.projected[current_index],
          polygon.projected[next_index]
        );
      if(!std::isfinite(ear_cross) || ear_cross <= polygon.area_tolerance) {
        continue;
      }

      bool contains_other_vertex = false;
      for(const auto candidate_index : active_indices) {
        if(candidate_index == previous_index ||
           candidate_index == current_index ||
           candidate_index == next_index) {
          continue;
        }

        if(obj_point_in_triangle_2d(
             polygon.projected[candidate_index],
             polygon.projected[previous_index],
             polygon.projected[current_index],
             polygon.projected[next_index],
             polygon.orientation_sign,
             polygon.area_tolerance
           )) {
          contains_other_vertex = true;
          break;
        }
      }

      if(contains_other_vertex) {
        continue;
      }

      triangles.push_back({
        polygon.refs[previous_index],
        polygon.refs[current_index],
        polygon.refs[next_index],
      });
      active_indices.erase(
        active_indices.begin() + static_cast<std::vector<std::size_t>::difference_type>(position)
      );
      found_ear = true;
      break;
    }

    if(!found_ear) {
      fail(
        base::StatusCode::unsupported,
        "OBJ import currently supports only simple planar polygon faces with non-collinear vertices; line " +
          std::to_string(line_number) + " could not be ear-clipped cleanly."
      );
    }
  }

  const auto final_cross =
    polygon.orientation_sign * obj_cross_2d(
      polygon.projected[active_indices[0]],
      polygon.projected[active_indices[1]],
      polygon.projected[active_indices[2]]
    );
  if(!std::isfinite(final_cross) || final_cross <= polygon.area_tolerance) {
    fail(
      base::StatusCode::unsupported,
      "OBJ import currently supports only non-degenerate polygon faces; line " +
        std::to_string(line_number) + " is degenerate."
    );
  }

  triangles.push_back({
    polygon.refs[active_indices[0]],
    polygon.refs[active_indices[1]],
    polygon.refs[active_indices[2]],
  });
  return triangles;
}

struct MshEntityInfo final {
  std::vector<int> physical_tags {};
};

void merge_msh_entity_info(
  std::map<std::pair<int, int>, MshEntityInfo> &entity_info,
  int dimension,
  int tag,
  MshEntityInfo info
)
{
  auto [entity_it, inserted] =
    entity_info.emplace(std::make_pair(dimension, tag), MshEntityInfo {});
  if(inserted) {
    entity_it->second = std::move(info);
    return;
  }

  for(const auto physical_tag : info.physical_tags) {
    if(std::find(
         entity_it->second.physical_tags.begin(),
         entity_it->second.physical_tags.end(),
         physical_tag
       ) == entity_it->second.physical_tags.end()) {
      entity_it->second.physical_tags.push_back(physical_tag);
    }
  }
}

struct EntityGroupKey final {
  EntityOrder order = EntityOrder::node;
  bool boundary = false;
  EntityKind kind = EntityKind::invalid;
  EntityGroupSemantic semantic = EntityGroupSemantic::unspecified;
  std::uint32_t zone_id = invalid_index;
  std::uint32_t source_entity_tag = invalid_index;
  std::uint32_t primary_region_zone_id = invalid_index;
  std::uint32_t secondary_region_zone_id = invalid_index;
  std::string name {};
  EntityGroupImportInfo import_info {};

  [[nodiscard]] bool operator<(const EntityGroupKey &other) const noexcept
  {
    if(order != other.order) {
      return order < other.order;
    }
    if(boundary != other.boundary) {
      return boundary < other.boundary;
    }
    if(kind != other.kind) {
      return kind < other.kind;
    }
    if(semantic != other.semantic) {
      return semantic < other.semantic;
    }
    if(zone_id != other.zone_id) {
      return zone_id < other.zone_id;
    }
    if(source_entity_tag != other.source_entity_tag) {
      return source_entity_tag < other.source_entity_tag;
    }
    if(primary_region_zone_id != other.primary_region_zone_id) {
      return primary_region_zone_id < other.primary_region_zone_id;
    }
    if(secondary_region_zone_id != other.secondary_region_zone_id) {
      return secondary_region_zone_id < other.secondary_region_zone_id;
    }
    if(import_info.format != other.import_info.format) {
      return import_info.format < other.import_info.format;
    }
    if(import_info.cgns.base_index != other.import_info.cgns.base_index) {
      return import_info.cgns.base_index < other.import_info.cgns.base_index;
    }
    if(import_info.cgns.base_name != other.import_info.cgns.base_name) {
      return import_info.cgns.base_name < other.import_info.cgns.base_name;
    }
    if(import_info.cgns.zone_index != other.import_info.cgns.zone_index) {
      return import_info.cgns.zone_index < other.import_info.cgns.zone_index;
    }
    if(import_info.cgns.zone_name != other.import_info.cgns.zone_name) {
      return import_info.cgns.zone_name < other.import_info.cgns.zone_name;
    }
    if(import_info.cgns.local_name != other.import_info.cgns.local_name) {
      return import_info.cgns.local_name < other.import_info.cgns.local_name;
    }
    if(import_info.cgns.bc_type_value != other.import_info.cgns.bc_type_value) {
      return import_info.cgns.bc_type_value < other.import_info.cgns.bc_type_value;
    }
    if(import_info.nastran.source_card != other.import_info.nastran.source_card) {
      return import_info.nastran.source_card < other.import_info.nastran.source_card;
    }
    if(import_info.nastran.property_card != other.import_info.nastran.property_card) {
      return import_info.nastran.property_card < other.import_info.nastran.property_card;
    }
    if(import_info.nastran.material_card != other.import_info.nastran.material_card) {
      return import_info.nastran.material_card < other.import_info.nastran.material_card;
    }
    if(import_info.nastran.material_id != other.import_info.nastran.material_id) {
      return import_info.nastran.material_id < other.import_info.nastran.material_id;
    }
    return name < other.name;
  }
};

[[nodiscard]] GroupInfo normalize_group(
  GroupInfo group,
  EntityOrder order,
  bool boundary,
  EntityGroupSemantic semantic,
  std::uint32_t primary_region_zone_id,
  std::uint32_t secondary_region_zone_id
)
{
  if(group.name.empty()) {
    group.name = default_entity_group_name(
      order,
      boundary,
      semantic,
      primary_region_zone_id,
      secondary_region_zone_id
    );
  }
  return group;
}

[[nodiscard]] EntityGroupIndex ensure_entity_group(
  Domain &domain,
  std::map<EntityGroupKey, EntityGroupIndex> &entity_group_map,
  EntityOrder order,
  bool boundary,
  GroupInfo group,
  EntityGroupSemantic semantic = EntityGroupSemantic::unspecified,
  std::uint32_t primary_region_zone_id = invalid_index,
  std::uint32_t secondary_region_zone_id = invalid_index
)
{
  const auto normalized_semantic =
    semantic == EntityGroupSemantic::unspecified
      ? (order == EntityOrder::node
           ? EntityGroupSemantic::node
           : (order == EntityOrder::cell
                ? EntityGroupSemantic::region
                : (boundary ? EntityGroupSemantic::boundary : EntityGroupSemantic::interior)))
      : semantic;
  group = normalize_group(
    std::move(group),
    order,
    boundary,
    normalized_semantic,
    primary_region_zone_id,
    secondary_region_zone_id
  );

  EntityGroupKey key;
  key.order = order;
  key.boundary = boundary;
  key.kind = default_kind(order);
  key.semantic = normalized_semantic;
  key.zone_id = group.zone_id;
  key.source_entity_tag = group.source_entity_tag;
  key.primary_region_zone_id = primary_region_zone_id;
  key.secondary_region_zone_id = secondary_region_zone_id;
  key.import_info = group.import_info;
  key.name = group.name;

  const auto entity_group_it = entity_group_map.find(key);
  if(entity_group_it != entity_group_map.end()) {
    return entity_group_it->second;
  }

  EntityGroupDefinition definition;
  definition.order = order;
  definition.name = group.name;
  definition.zone_id = group.zone_id;
  definition.source_entity_tag = group.source_entity_tag;
  definition.boundary = boundary;
  definition.default_kind = key.kind;
  definition.semantic = normalized_semantic;
  definition.primary_region_zone_id = primary_region_zone_id;
  definition.secondary_region_zone_id = secondary_region_zone_id;
  definition.import_info = group.import_info;

  const auto entity_group = domain.create_entity_group(std::move(definition));
  entity_group_map.emplace(std::move(key), entity_group);
  return entity_group;
}

using EdgeKey = std::array<std::size_t, 2>;
using FaceKey = std::array<std::size_t, 3>;

[[nodiscard]] EdgeKey make_edge_key(std::array<std::size_t, 2> nodes) noexcept
{
  if(nodes[1] < nodes[0]) {
    std::swap(nodes[0], nodes[1]);
  }
  return nodes;
}

[[nodiscard]] FaceKey make_face_key(std::array<std::size_t, 3> nodes) noexcept
{
  std::sort(nodes.begin(), nodes.end());
  return nodes;
}

struct FaceBuildRecord final {
  std::array<std::size_t, 3> nodes {};
  GroupInfo group {};
  std::uint32_t source_entity_tag = invalid_index;
  bool explicit_record = false;
  std::vector<std::size_t> cell_indices {};
  EntityRef ref {};
};

[[nodiscard]] EntityGroupSemantic classify_face_entity_group_semantic(
  const FaceBuildRecord &record,
  const std::vector<EntityGroupIndex> &cell_entity_groups,
  const Domain &domain,
  std::uint32_t &primary_region_zone_id,
  std::uint32_t &secondary_region_zone_id
)
{
  primary_region_zone_id = invalid_index;
  secondary_region_zone_id = invalid_index;

  if(record.cell_indices.empty()) {
    return EntityGroupSemantic::boundary;
  }

  primary_region_zone_id =
    domain.entity_group(cell_entity_groups[record.cell_indices.front()]).primary_region_zone_id();
  if(record.cell_indices.size() == 1U) {
    return EntityGroupSemantic::boundary;
  }

  secondary_region_zone_id =
    domain.entity_group(cell_entity_groups[record.cell_indices.back()]).primary_region_zone_id();
  if(primary_region_zone_id != invalid_index &&
     secondary_region_zone_id != invalid_index &&
     primary_region_zone_id != secondary_region_zone_id) {
    return EntityGroupSemantic::interface;
  }

  secondary_region_zone_id = invalid_index;
  return EntityGroupSemantic::interior;
}

[[nodiscard]] std::array<std::array<std::size_t, 3>, 4> tetra_faces(
  const std::array<std::size_t, 4> &nodes
) noexcept
{
  return {{
    {nodes[0], nodes[2], nodes[1]},
    {nodes[0], nodes[1], nodes[3]},
    {nodes[1], nodes[2], nodes[3]},
    {nodes[2], nodes[0], nodes[3]},
  }};
}

[[nodiscard]] Domain build_domain_from_description(
  const ImportedMeshDescription &description
)
{
  Domain domain(description.name);
  std::map<EntityGroupKey, EntityGroupIndex> entity_group_map;

  const auto node_entity_group = ensure_entity_group(
    domain,
    entity_group_map,
    EntityOrder::node,
    false,
    {}
  );

  std::vector<EntityRef> node_refs;
  node_refs.reserve(description.nodes.size());
  for(const auto &node : description.nodes) {
    node_refs.push_back(domain.add_node(node_entity_group, node));
  }

  std::map<FaceKey, FaceBuildRecord> face_records;
  for(const auto &face : description.faces) {
    const auto key = make_face_key(face.nodes);
    auto [record_it, inserted] = face_records.emplace(key, FaceBuildRecord {});
    if(!inserted && record_it->second.explicit_record) {
      fail(
        base::StatusCode::unsupported,
        "Mesh IO does not support duplicate triangle records with the same connectivity."
      );
    }
    record_it->second.nodes = face.nodes;
    record_it->second.group = face.group;
    record_it->second.source_entity_tag =
      face.source_entity_tag != invalid_index
        ? face.source_entity_tag
        : face.group.source_entity_tag;
    record_it->second.explicit_record = true;
  }

  std::vector<std::array<FaceKey, 4>> cell_faces;
  cell_faces.reserve(description.cells.size());
  for(std::size_t cell_index = 0; cell_index < description.cells.size(); ++cell_index) {
    const auto local_faces = tetra_faces(description.cells[cell_index].nodes);
    std::array<FaceKey, 4> face_keys {};

    for(std::size_t face_offset = 0; face_offset < local_faces.size(); ++face_offset) {
      const auto key = make_face_key(local_faces[face_offset]);
      face_keys[face_offset] = key;
      auto [record_it, inserted] = face_records.emplace(key, FaceBuildRecord {});
      if(inserted) {
        record_it->second.nodes = local_faces[face_offset];
      }
      record_it->second.cell_indices.push_back(cell_index);
    }

    cell_faces.push_back(face_keys);
  }

  std::vector<EntityGroupIndex> cell_entity_groups(description.cells.size(), invalid_index);
  for(std::size_t cell_index = 0; cell_index < description.cells.size(); ++cell_index) {
    cell_entity_groups[cell_index] = ensure_entity_group(
      domain,
      entity_group_map,
      EntityOrder::cell,
      false,
      description.cells[cell_index].group,
      EntityGroupSemantic::region,
      description.cells[cell_index].group.zone_id
    );
  }

  for(auto &entry : face_records) {
    auto &record = entry.second;
    if(record.cell_indices.size() > 2U) {
      fail(
        base::StatusCode::unsupported,
        "Mesh IO does not support non-manifold tetrahedral face adjacency."
      );
    }

    const bool boundary = record.cell_indices.size() <= 1U;
    std::uint32_t primary_region_zone_id = invalid_index;
    std::uint32_t secondary_region_zone_id = invalid_index;
    const auto semantic = classify_face_entity_group_semantic(
      record,
      cell_entity_groups,
      domain,
      primary_region_zone_id,
      secondary_region_zone_id
    );
    GroupInfo group = record.explicit_record ? record.group : GroupInfo {};
    const auto face_entity_group = ensure_entity_group(
      domain,
      entity_group_map,
      EntityOrder::face,
      boundary,
      std::move(group),
      semantic,
      primary_region_zone_id,
      secondary_region_zone_id
    );

    record.ref = domain.add_triangle_face(
      face_entity_group,
      {
        node_refs[record.nodes[0]],
        node_refs[record.nodes[1]],
        node_refs[record.nodes[2]],
      }
    );
    if(record.source_entity_tag != invalid_index) {
      domain.set_face_source_entity_tag(record.ref, record.source_entity_tag);
    }
  }

  std::vector<EntityRef> cell_refs;
  cell_refs.reserve(description.cells.size());
  for(std::size_t cell_index = 0; cell_index < description.cells.size(); ++cell_index) {
    const auto &cell = description.cells[cell_index];
    const auto cell_entity_group = cell_entity_groups[cell_index];

    std::array<EntityRef, 4> face_refs {};
    for(std::size_t face_offset = 0; face_offset < face_refs.size(); ++face_offset) {
      face_refs[face_offset] = face_records.at(cell_faces[cell_index][face_offset]).ref;
    }

    cell_refs.push_back(
      domain.add_tetra_cell(
        cell_entity_group,
        {
          node_refs[cell.nodes[0]],
          node_refs[cell.nodes[1]],
          node_refs[cell.nodes[2]],
          node_refs[cell.nodes[3]],
        },
        face_refs
      )
    );
  }

  for(const auto &entry : face_records) {
    const auto &record = entry.second;
    if(record.cell_indices.empty()) {
      continue;
    }

    const auto left_cell = cell_refs[record.cell_indices[0]];
    const auto right_cell =
      record.cell_indices.size() == 2U ? cell_refs[record.cell_indices[1]] : EntityRef {};
    domain.set_face_cells(record.ref, left_cell, right_cell);
  }

  std::map<EdgeKey, std::vector<EntityRef>> edge_face_adjacency;
  for(const auto &entry : face_records) {
    const auto &record = entry.second;
    edge_face_adjacency[make_edge_key({record.nodes[0], record.nodes[1]})].push_back(record.ref);
    edge_face_adjacency[make_edge_key({record.nodes[1], record.nodes[2]})].push_back(record.ref);
    edge_face_adjacency[make_edge_key({record.nodes[2], record.nodes[0]})].push_back(record.ref);
  }

  std::map<EdgeKey, bool> seen_edges;
  for(const auto &edge : description.edges) {
    const auto key = make_edge_key(edge.nodes);
    if(!seen_edges.emplace(key, true).second) {
      fail(
        base::StatusCode::unsupported,
        "Mesh IO does not support duplicate line records with the same connectivity."
      );
    }

    const auto adjacency_it = edge_face_adjacency.find(key);
    const auto adjacent_face_count =
      adjacency_it == edge_face_adjacency.end() ? 0U : adjacency_it->second.size();
    if(adjacent_face_count > 2U) {
      fail(
        base::StatusCode::unsupported,
        "Mesh IO does not support non-manifold line-to-face adjacency."
      );
    }

    const bool boundary = adjacent_face_count <= 1U;
    std::uint32_t primary_region_zone_id = invalid_index;
    if(boundary && adjacency_it != edge_face_adjacency.end() && !adjacency_it->second.empty()) {
      primary_region_zone_id =
        domain.entity_group(adjacency_it->second.front().entity_group).primary_region_zone_id();
    }
    const auto edge_entity_group = ensure_entity_group(
      domain,
      entity_group_map,
      EntityOrder::edge,
      boundary,
      edge.group,
      boundary ? EntityGroupSemantic::boundary : EntityGroupSemantic::interior,
      primary_region_zone_id
    );

    const auto edge_ref = domain.add_edge(
      edge_entity_group,
      {
        node_refs[edge.nodes[0]],
        node_refs[edge.nodes[1]],
      }
    );

    if(adjacent_face_count > 0U) {
      const auto right_face =
        adjacent_face_count == 2U ? adjacency_it->second[1] : EntityRef {};
      domain.set_edge_faces(edge_ref, adjacency_it->second[0], right_face);
    }
  }

  return domain;
}

[[nodiscard]] base::ContextHandle require_context(base::ContextHandle context_handle)
{
  if(context_handle != sqmesh::invalid_handle) {
    return context_handle;
  }

  const auto current = base::current_context();
  if(current == sqmesh::invalid_handle) {
    fail(
      base::last_error_code(),
      std::string(base::last_error_message())
    );
  }
  return current;
}

[[nodiscard]] base::StatusCode store_imported_domain(
  Domain domain,
  std::string_view algorithm_name,
  MeshHandle &mesh_handle,
  base::ContextHandle context_handle
)
{
  mesh_handle = sqmesh::invalid_handle;

  geo::ModelHandle model_handle = sqmesh::invalid_handle;
  const auto status = geo::create_placeholder_model(model_handle, context_handle);
  if(status != base::StatusCode::ok) {
    return status;
  }

  auto domain_storage = std::make_shared<Domain>(std::move(domain));
  const auto summary = domain_storage->summary();
  return core::detail::store_mesh(
    model_handle,
    algorithm_name,
    summary,
    std::move(domain_storage),
    mesh_handle,
    context_handle
  );
}

void validate_non_empty_path(
  std::string_view path,
  std::string_view message
)
{
  if(path.empty()) {
    fail(base::StatusCode::invalid_argument, std::string(message));
  }
}

[[nodiscard]] long long parse_integer(
  std::string_view token,
  std::string_view context
)
{
  try {
    std::size_t consumed = 0U;
    const auto value = std::stoll(std::string(token), &consumed);
    if(consumed != token.size()) {
      fail(
        base::StatusCode::io_error,
        std::string(context) + " contains a non-integer token."
      );
    }
    return value;
  }
  catch(const std::exception &) {
    fail(
      base::StatusCode::io_error,
      std::string(context) + " contains a malformed integer token."
    );
  }
}

enum class NastranRecordFormat : std::uint8_t {
  free_field = 0,
  small_fixed = 1,
  large_fixed = 2,
};

struct NastranPhysicalLine final {
  std::string raw {};
  std::size_t line_number = 0U;
};

struct NastranRecord final {
  std::string card {};
  std::vector<std::string> fields {};
  std::size_t line_number = 0U;
  NastranRecordFormat format = NastranRecordFormat::free_field;
};

[[nodiscard]] std::vector<std::string> split_nastran_free_fields(std::string_view value)
{
  std::vector<std::string> fields;
  std::size_t begin = 0U;
  while(begin <= value.size()) {
    const auto comma = value.find(',', begin);
    if(comma == std::string_view::npos) {
      fields.push_back(trim_copy(value.substr(begin)));
      break;
    }
    fields.push_back(trim_copy(value.substr(begin, comma - begin)));
    begin = comma + 1U;
  }
  return fields;
}

[[nodiscard]] std::string normalize_nastran_card_name(std::string_view token)
{
  auto normalized = trim_upper_copy(token);
  while(!normalized.empty() && normalized.back() == '*') {
    normalized.pop_back();
  }
  return normalized;
}

[[nodiscard]] bool is_nastran_continuation_line(std::string_view trimmed) noexcept
{
  return !trimmed.empty() && (trimmed.front() == '+' || trimmed.front() == '*');
}

void trim_trailing_empty_fields(
  std::vector<std::string> &fields
)
{
  while(!fields.empty() && fields.back().empty()) {
    fields.pop_back();
  }
}

void append_nastran_small_fixed_fields(
  std::string_view line,
  std::vector<std::string> &fields
)
{
  for(std::size_t begin = 8U; begin < line.size(); begin += 8U) {
    fields.push_back(
      trim_copy(line.substr(begin, std::min<std::size_t>(8U, line.size() - begin)))
    );
  }
}

void append_nastran_large_fixed_fields(
  std::string_view line,
  std::vector<std::string> &fields
)
{
  for(const std::size_t begin : {8U, 24U, 40U, 56U}) {
    if(begin >= line.size()) {
      break;
    }
    fields.push_back(
      trim_copy(line.substr(begin, std::min<std::size_t>(16U, line.size() - begin)))
    );
  }
}

[[nodiscard]] bool has_non_empty_fields(
  const std::vector<std::string> &fields,
  std::size_t first
) noexcept
{
  for(std::size_t index = first; index < fields.size(); ++index) {
    if(!fields[index].empty()) {
      return true;
    }
  }
  return false;
}

[[nodiscard]] bool is_zero_or_blank(std::string_view token)
{
  if(trim_copy(token).empty()) {
    return true;
  }
  return parse_integer(trim_copy(token), "NASTRAN integer field") == 0;
}

[[nodiscard]] long long parse_positive_integer(
  std::string_view token,
  std::string_view context
)
{
  const auto value = parse_integer(trim_copy(token), context);
  if(value <= 0) {
    fail(
      base::StatusCode::unsupported,
      std::string(context) + " requires a positive integer value."
    );
  }
  return value;
}

[[nodiscard]] std::uint32_t parse_nastran_property_id(
  std::string_view token,
  std::string_view context
)
{
  const auto value = parse_positive_integer(token, context);
  if(value >= static_cast<long long>(invalid_index)) {
    fail(
      base::StatusCode::unsupported,
      std::string(context) + " exceeds the supported SQMesh entity_group zone-id range."
    );
  }
  return static_cast<std::uint32_t>(value);
}

[[nodiscard]] std::uint32_t parse_nastran_material_id(
  std::string_view token,
  std::string_view context
)
{
  const auto value = parse_positive_integer(token, context);
  if(value >= static_cast<long long>(invalid_index)) {
    fail(
      base::StatusCode::unsupported,
      std::string(context) + " exceeds the supported SQMesh material-id range."
    );
  }
  return static_cast<std::uint32_t>(value);
}

[[nodiscard]] std::uint32_t parse_nastran_source_entity_tag(
  long long value,
  std::string_view context
)
{
  if(value >= static_cast<long long>(invalid_index)) {
    fail(
      base::StatusCode::unsupported,
      std::string(context) + " exceeds the supported SQMesh source-entity-tag range."
    );
  }
  return static_cast<std::uint32_t>(value);
}

[[nodiscard]] double parse_nastran_real(
  std::string_view token,
  std::string_view context
)
{
  const auto trimmed = trim_copy(token);
  if(trimmed.empty()) {
    fail(
      base::StatusCode::io_error,
      std::string(context) + " is missing a real-valued field."
    );
  }

  std::string normalized(trimmed);
  bool has_exponent = false;
  for(auto &ch : normalized) {
    if(ch == 'D' || ch == 'd') {
      ch = 'E';
      has_exponent = true;
    }
    else if(ch == 'E' || ch == 'e') {
      has_exponent = true;
    }
  }

  if(!has_exponent) {
    std::string compact;
    compact.reserve(normalized.size() + 2U);
    bool seen_non_space = false;
    for(std::size_t index = 0; index < normalized.size(); ++index) {
      const auto ch = normalized[index];
      if(!seen_non_space && ch != ' ') {
        seen_non_space = true;
      }
      if(index != 0U && seen_non_space && (ch == '+' || ch == '-')) {
        compact.push_back('E');
      }
      compact.push_back(ch);
    }
    normalized = std::move(compact);
  }

  try {
    std::size_t consumed = 0U;
    const auto value = std::stod(normalized, &consumed);
    if(consumed != normalized.size()) {
      fail(
        base::StatusCode::io_error,
        std::string(context) + " contains a malformed real token."
      );
    }
    return value;
  }
  catch(const std::exception &) {
    fail(
      base::StatusCode::io_error,
      std::string(context) + " contains a malformed real token."
    );
  }
}

[[nodiscard]] double parse_optional_nastran_real(
  std::string_view token,
  std::string_view context
)
{
  if(trim_copy(token).empty()) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return parse_nastran_real(token, context);
}

[[nodiscard]] bool is_supported_nastran_card(std::string_view card) noexcept
{
  return card == "GRID" || card == "CBAR" || card == "CBEAM" ||
         card == "CROD" || card == "CTRIA3" || card == "CQUAD4" ||
         card == "CTETRA" || card == "MAT1" || card == "PROD" ||
         card == "PSHELL" || card == "PSOLID";
}

[[nodiscard]] bool is_supported_nastran_mesh_card(std::string_view card) noexcept
{
  return card == "GRID" || card == "CBAR" || card == "CBEAM" ||
         card == "CROD" || card == "CTRIA3" || card == "CQUAD4" ||
         card == "CTETRA";
}

[[nodiscard]] bool are_blank_or_zero_fields(
  const std::vector<std::string> &fields,
  std::size_t first
)
{
  for(std::size_t index = first; index < fields.size(); ++index) {
    if(!is_zero_or_blank(fields[index])) {
      return false;
    }
  }
  return true;
}

[[nodiscard]] std::string detect_nastran_card(std::string_view line)
{
  const auto comma = line.find(',');
  if(comma != std::string_view::npos) {
    return normalize_nastran_card_name(line.substr(0U, comma));
  }
  return normalize_nastran_card_name(line.substr(0U, std::min<std::size_t>(8U, line.size())));
}

[[nodiscard]] NastranRecord parse_nastran_record(
  const std::vector<NastranPhysicalLine> &lines,
  const NastranImportOptions &options
)
{
  if(!options.allow_free_fields && !options.allow_small_fixed_fields &&
     !options.allow_large_fixed_fields) {
    fail(
      base::StatusCode::invalid_argument,
      "NASTRAN import requires at least one enabled field format."
    );
  }

  if(lines.empty()) {
    fail(base::StatusCode::io_error, "NASTRAN import found an empty logical record.");
  }

  const auto &line = lines.front().raw;
  const auto line_number = lines.front().line_number;
  const auto context = "NASTRAN line " + std::to_string(line_number);
  const auto trimmed = trim_copy(line);
  if(trimmed.empty()) {
    fail(base::StatusCode::io_error, context + " is empty.");
  }
  if(is_nastran_continuation_line(trimmed)) {
    fail(
      base::StatusCode::unsupported,
      context + " starts with a continuation marker without a preceding supported card."
    );
  }

  NastranRecord record;
  record.line_number = line_number;

  if(line.find(',') != std::string_view::npos) {
    if(!options.allow_free_fields) {
      fail(
        base::StatusCode::unsupported,
        context + " uses free-field syntax, which this import call disabled."
      );
    }

    auto tokens = split_nastran_free_fields(line);
    if(tokens.empty() || tokens[0].empty()) {
      fail(base::StatusCode::io_error, context + " is missing a card name.");
    }

    record.card = normalize_nastran_card_name(tokens[0]);
    if(!is_supported_nastran_card(record.card)) {
      fail(
        base::StatusCode::unsupported,
        context + " contains unsupported card `" + record.card + "`."
      );
    }

    record.format = NastranRecordFormat::free_field;
    record.fields.assign(tokens.begin() + 1, tokens.end());
    for(std::size_t index = 1U; index < lines.size(); ++index) {
      const auto continuation = trim_copy(lines[index].raw);
      if(!is_nastran_continuation_line(continuation)) {
        fail(
          base::StatusCode::io_error,
          context + " has a malformed free-field continuation record."
        );
      }

      auto continuation_tokens = split_nastran_free_fields(lines[index].raw);
      if(continuation_tokens.empty() || continuation_tokens[0].empty() ||
         !is_nastran_continuation_line(trim_copy(continuation_tokens[0]))) {
        fail(
          base::StatusCode::unsupported,
          "NASTRAN line " + std::to_string(lines[index].line_number) +
            " uses an unsupported free-field continuation identifier."
        );
      }

      record.fields.insert(
        record.fields.end(),
        continuation_tokens.begin() + 1,
        continuation_tokens.end()
      );
    }
    trim_trailing_empty_fields(record.fields);
    return record;
  }

  const auto card_field =
    line.substr(0U, std::min<std::size_t>(8U, line.size()));
  const auto large_fixed = card_field.find('*') != std::string_view::npos;

  if(large_fixed) {
    if(!options.allow_large_fixed_fields) {
      fail(
        base::StatusCode::unsupported,
        context + " uses large-field syntax, which this import call disabled."
      );
    }

    record.card = normalize_nastran_card_name(card_field);
    if(record.card.empty()) {
      fail(base::StatusCode::io_error, context + " is missing a card name.");
    }
    if(!is_supported_nastran_card(record.card)) {
      fail(
        base::StatusCode::unsupported,
        context + " contains unsupported card `" + record.card + "`."
      );
    }

    record.format = NastranRecordFormat::large_fixed;
    append_nastran_large_fixed_fields(line, record.fields);
    for(std::size_t index = 1U; index < lines.size(); ++index) {
      const auto continuation = trim_copy(lines[index].raw);
      if(continuation.empty() || continuation.front() != '*') {
        fail(
          base::StatusCode::unsupported,
          "NASTRAN line " + std::to_string(lines[index].line_number) +
            " uses an unsupported large-field continuation identifier."
        );
      }
      append_nastran_large_fixed_fields(lines[index].raw, record.fields);
    }
    trim_trailing_empty_fields(record.fields);
    return record;
  }

  if(!options.allow_small_fixed_fields) {
    fail(
      base::StatusCode::unsupported,
      context + " uses small-field syntax, which this import call disabled."
    );
  }

  record.card = normalize_nastran_card_name(card_field);
  if(record.card.empty()) {
    fail(base::StatusCode::io_error, context + " is missing a card name.");
  }
  if(!is_supported_nastran_card(record.card)) {
    fail(
      base::StatusCode::unsupported,
      context + " contains unsupported card `" + record.card + "`."
    );
  }
  if(lines.size() > 1U) {
    fail(
      base::StatusCode::unsupported,
      context + " uses a continued 8-character small-field record, which SQMesh does not support."
    );
  }

  record.format = NastranRecordFormat::small_fixed;
  append_nastran_small_fixed_fields(line, record.fields);
  trim_trailing_empty_fields(record.fields);
  return record;
}

[[nodiscard]] GroupInfo nastran_group(
  EntityOrder order,
  std::uint32_t property_id,
  std::string_view card = {},
  const NastranPropertyInfo *property = nullptr,
  const NastranMaterialInfo *material = nullptr
)
{
  GroupInfo group;
  group.zone_id = property_id;
  switch(order) {
  case EntityOrder::edge:
    group.name = "line_pid_" + std::to_string(property_id);
    group.import_info = make_nastran_entity_group_import_info(card, property, material);
    break;
  case EntityOrder::face:
    group.name =
      card == "CQUAD4"
        ? "cquad4_pid_" + std::to_string(property_id)
        : "surface_pid_" + std::to_string(property_id);
    group.import_info = make_nastran_entity_group_import_info(card, property, material);
    break;
  case EntityOrder::cell:
    group.name = "volume_pid_" + std::to_string(property_id);
    group.import_info = make_nastran_entity_group_import_info(card, property, material);
    break;
  default:
    break;
  }
  return group;
}

[[nodiscard]] ImportedMeshDescription parse_nastran_description(
  std::string_view path,
  const NastranImportOptions &options
)
{
  validate_non_empty_path(path, "NASTRAN import requires a non-empty path.");

  const std::string file_path(path);
  std::ifstream input(file_path);
  if(!input) {
    fail(base::StatusCode::io_error, "NASTRAN import could not open the input file.");
  }

  ImportedMeshDescription description;
  description.name = default_domain_name(path);

  std::vector<NastranRecord> records;
  std::vector<NastranPhysicalLine> pending_record;
  std::string line;
  bool in_bulk = false;
  std::size_t line_number = 0U;
  const auto flush_pending_record = [&]() {
    if(pending_record.empty()) {
      return;
    }
    records.push_back(parse_nastran_record(pending_record, options));
    pending_record.clear();
  };

  while(std::getline(input, line)) {
    ++line_number;
    const auto raw = rtrim_copy(line);
    const auto trimmed = trim_copy(raw);
    if(trimmed.empty() || trimmed[0] == '$') {
      continue;
    }

    const auto upper = trim_upper_copy(trimmed);
    if(upper == "ENDDATA") {
      flush_pending_record();
      break;
    }

    if(!in_bulk) {
      if(upper == "BEGIN BULK") {
        in_bulk = true;
        continue;
      }
      if(is_nastran_continuation_line(trimmed)) {
        continue;
      }

      const auto card = detect_nastran_card(raw);
      if(!is_supported_nastran_card(card)) {
        continue;
      }
      in_bulk = true;
    }

    if(upper == "BEGIN BULK") {
      flush_pending_record();
      continue;
    }

    if(is_nastran_continuation_line(trimmed)) {
      if(pending_record.empty()) {
        fail(
          base::StatusCode::unsupported,
          "NASTRAN line " + std::to_string(line_number) +
            " starts with a continuation marker without a preceding supported card."
        );
      }
      pending_record.push_back({raw, line_number});
      continue;
    }

    flush_pending_record();
    pending_record.push_back({raw, line_number});
  }
  flush_pending_record();

  std::map<long long, std::size_t> node_map;
  std::map<std::uint32_t, NastranMaterialInfo> material_info_by_id;
  std::map<std::uint32_t, NastranPropertyInfo> property_info_by_id;
  for(const auto &record : records) {
    if(record.card == "MAT1") {
      if(record.fields.empty()) {
        fail(base::StatusCode::io_error, "NASTRAN MAT1 record is missing the material id.");
      }
      const auto material_id =
        parse_nastran_material_id(record.fields[0], "NASTRAN MAT1 material id");
      NastranMaterialInfo material_info;
      material_info.card = NastranMaterialCard::mat1;
      material_info.material_id = material_id;
      if(record.fields.size() > 1U) {
        material_info.youngs_modulus =
          parse_optional_nastran_real(record.fields[1], "NASTRAN MAT1 E");
      }
      if(record.fields.size() > 2U) {
        material_info.shear_modulus =
          parse_optional_nastran_real(record.fields[2], "NASTRAN MAT1 G");
      }
      if(record.fields.size() > 3U) {
        material_info.poisson_ratio =
          parse_optional_nastran_real(record.fields[3], "NASTRAN MAT1 NU");
      }
      if(!material_info_by_id.emplace(material_id, material_info).second) {
        fail(base::StatusCode::io_error, "NASTRAN import found a duplicate MAT1 material id.");
      }
      continue;
    }

    if(record.card == "PROD") {
      if(record.fields.size() < 2U) {
        fail(base::StatusCode::io_error, "NASTRAN PROD record is missing required fields.");
      }
      const auto property_id =
        parse_nastran_property_id(record.fields[0], "NASTRAN PROD property id");
      NastranPropertyInfo property_info;
      property_info.card = NastranPropertyCard::prod;
      if(!trim_copy(record.fields[1]).empty()) {
        property_info.material_id =
          parse_nastran_material_id(record.fields[1], "NASTRAN PROD material id");
      }
      if(record.fields.size() > 2U) {
        property_info.rod_area =
          parse_optional_nastran_real(record.fields[2], "NASTRAN PROD A");
      }
      if(!property_info_by_id.emplace(property_id, property_info).second) {
        fail(base::StatusCode::io_error, "NASTRAN import found a duplicate PROD property id.");
      }
      continue;
    }

    if(record.card == "PSHELL") {
      if(record.fields.size() < 2U) {
        fail(base::StatusCode::io_error, "NASTRAN PSHELL record is missing required fields.");
      }
      const auto property_id =
        parse_nastran_property_id(record.fields[0], "NASTRAN PSHELL property id");
      NastranPropertyInfo property_info;
      property_info.card = NastranPropertyCard::pshell;
      if(!trim_copy(record.fields[1]).empty()) {
        property_info.material_id =
          parse_nastran_material_id(record.fields[1], "NASTRAN PSHELL MID1");
      }
      if(record.fields.size() > 2U) {
        property_info.shell_thickness =
          parse_optional_nastran_real(record.fields[2], "NASTRAN PSHELL T");
      }
      if(!property_info_by_id.emplace(property_id, property_info).second) {
        fail(base::StatusCode::io_error, "NASTRAN import found a duplicate PSHELL property id.");
      }
      continue;
    }

    if(record.card == "PSOLID") {
      if(record.fields.size() < 2U) {
        fail(base::StatusCode::io_error, "NASTRAN PSOLID record is missing required fields.");
      }
      const auto property_id =
        parse_nastran_property_id(record.fields[0], "NASTRAN PSOLID property id");
      NastranPropertyInfo property_info;
      property_info.card = NastranPropertyCard::psolid;
      if(!trim_copy(record.fields[1]).empty()) {
        property_info.material_id =
          parse_nastran_material_id(record.fields[1], "NASTRAN PSOLID material id");
      }
      if(!property_info_by_id.emplace(property_id, property_info).second) {
        fail(base::StatusCode::io_error, "NASTRAN import found a duplicate PSOLID property id.");
      }
      continue;
    }

    if(record.card != "GRID") {
      continue;
    }

    if(record.fields.size() < 5U) {
      fail(
        base::StatusCode::io_error,
        "NASTRAN GRID line " + std::to_string(record.line_number) +
          " is missing required fields."
      );
    }
    if(!is_zero_or_blank(record.fields[1])) {
      fail(
        base::StatusCode::unsupported,
        "NASTRAN import currently supports only GRID cards with blank or zero CP."
      );
    }
    if(record.fields.size() > 5U && !is_zero_or_blank(record.fields[5])) {
      fail(
        base::StatusCode::unsupported,
        "NASTRAN import currently supports only GRID cards with blank or zero CD."
      );
    }
    if(!are_blank_or_zero_fields(record.fields, 6U)) {
      fail(
        base::StatusCode::unsupported,
        "NASTRAN import currently supports only GRID cards with blank or zero PS/SEID fields."
      );
    }

    const auto node_id = parse_positive_integer(
      record.fields[0],
      "NASTRAN GRID node id"
    );
    if(node_map.find(node_id) != node_map.end()) {
      fail(base::StatusCode::io_error, "NASTRAN import found a duplicate GRID id.");
    }

    node_map.emplace(node_id, description.nodes.size());
    description.nodes.push_back({
      parse_nastran_real(record.fields[2], "NASTRAN GRID X"),
      parse_nastran_real(record.fields[3], "NASTRAN GRID Y"),
      parse_nastran_real(record.fields[4], "NASTRAN GRID Z"),
    });
  }

  if(description.nodes.empty()) {
    fail(base::StatusCode::io_error, "NASTRAN import found no GRID nodes.");
  }

  std::map<long long, bool> seen_element_ids;
  for(const auto &record : records) {
    if(!is_supported_nastran_mesh_card(record.card) || record.card == "GRID") {
      continue;
    }

    const auto element_id = parse_positive_integer(
      record.fields.empty() ? std::string_view {} : std::string_view(record.fields[0]),
      "NASTRAN element id"
    );
    if(!seen_element_ids.emplace(element_id, true).second) {
      fail(base::StatusCode::io_error, "NASTRAN import found a duplicate element id.");
    }

    const auto property_id = parse_nastran_property_id(
      record.fields.size() > 1U ? std::string_view(record.fields[1]) : std::string_view {},
      "NASTRAN property id"
    );
    const auto property_it = property_info_by_id.find(property_id);
    const NastranPropertyInfo *property_info =
      property_it != property_info_by_id.end() ? &property_it->second : nullptr;
    const NastranMaterialInfo *material_info = nullptr;
    if(property_info != nullptr && property_info->material_id != invalid_index) {
      const auto material_it = material_info_by_id.find(property_info->material_id);
      if(material_it != material_info_by_id.end()) {
        material_info = &material_it->second;
      }
    }

    auto map_node = [&](std::size_t field_index, std::string_view context) {
      if(field_index >= record.fields.size()) {
        fail(
          base::StatusCode::io_error,
          std::string(context) + " is missing a node reference."
        );
      }
      const auto node_id = parse_positive_integer(record.fields[field_index], context);
      const auto node_it = node_map.find(node_id);
      if(node_it == node_map.end()) {
        fail(
          base::StatusCode::io_error,
          std::string(context) + " references an undefined GRID id."
        );
      }
      return node_it->second;
    };

    if(record.card == "CBAR" || record.card == "CBEAM" || record.card == "CROD") {
      if(record.fields.size() < 4U) {
        fail(
          base::StatusCode::io_error,
          "NASTRAN " + record.card + " record is missing required fields."
        );
      }
      description.edges.push_back({
        {
          map_node(2U, "NASTRAN " + record.card + " G1"),
          map_node(3U, "NASTRAN " + record.card + " G2"),
        },
        nastran_group(EntityOrder::edge, property_id, record.card, property_info, material_info),
      });
      continue;
    }

    if(record.card == "CTRIA3") {
      if(record.fields.size() < 5U) {
        fail(base::StatusCode::io_error, "NASTRAN CTRIA3 record is missing required fields.");
      }
      append_imported_face(
        description.faces,
        {
          map_node(2U, "NASTRAN CTRIA3 G1"),
          map_node(3U, "NASTRAN CTRIA3 G2"),
          map_node(4U, "NASTRAN CTRIA3 G3"),
        },
        nastran_group(EntityOrder::face, property_id, record.card, property_info, material_info),
        parse_nastran_source_entity_tag(element_id, "NASTRAN CTRIA3 element id")
      );
      continue;
    }

    if(record.card == "CQUAD4") {
      if(record.fields.size() < 6U) {
        fail(base::StatusCode::io_error, "NASTRAN CQUAD4 record is missing required fields.");
      }

      const auto group =
        nastran_group(EntityOrder::face, property_id, record.card, property_info, material_info);
      const auto g1 = map_node(2U, "NASTRAN CQUAD4 G1");
      const auto g2 = map_node(3U, "NASTRAN CQUAD4 G2");
      const auto g3 = map_node(4U, "NASTRAN CQUAD4 G3");
      const auto g4 = map_node(5U, "NASTRAN CQUAD4 G4");
      append_imported_face(
        description.faces,
        {g1, g2, g3},
        group,
        parse_nastran_source_entity_tag(element_id, "NASTRAN CQUAD4 element id")
      );
      append_imported_face(
        description.faces,
        {g1, g3, g4},
        group,
        parse_nastran_source_entity_tag(element_id, "NASTRAN CQUAD4 element id")
      );
      continue;
    }

    if(record.card == "CTETRA") {
      if(record.fields.size() < 6U) {
        fail(base::StatusCode::io_error, "NASTRAN CTETRA record is missing required fields.");
      }
      if(has_non_empty_fields(record.fields, 6U)) {
        fail(
          base::StatusCode::unsupported,
          "NASTRAN import currently supports only first-order 4-node CTETRA cards."
        );
      }
      description.cells.push_back({
        {
          map_node(2U, "NASTRAN CTETRA G1"),
          map_node(3U, "NASTRAN CTETRA G2"),
          map_node(4U, "NASTRAN CTETRA G3"),
          map_node(5U, "NASTRAN CTETRA G4"),
        },
        nastran_group(EntityOrder::cell, property_id, record.card, property_info, material_info),
      });
      continue;
    }
  }

  finalize_nastran_entity_group_import_info(description);
  return description;
}

[[nodiscard]] std::vector<std::string> split_tokens(std::string_view value)
{
  const std::string text(value);
  std::istringstream stream(text);
  std::vector<std::string> tokens;
  std::string token;
  while(stream >> token) {
    tokens.push_back(token);
  }
  return tokens;
}

[[nodiscard]] double parse_real_token(
  std::string_view token,
  std::string_view context
)
{
  const auto trimmed = trim_copy(token);
  if(trimmed.empty()) {
    fail(
      base::StatusCode::io_error,
      std::string(context) + " is missing a real-valued token."
    );
  }

  try {
    std::size_t consumed = 0U;
    const auto value = std::stod(trimmed, &consumed);
    if(consumed != trimmed.size()) {
      fail(
        base::StatusCode::io_error,
        std::string(context) + " contains a malformed real token."
      );
    }
    return value;
  }
  catch(const std::exception &) {
    fail(
      base::StatusCode::io_error,
      std::string(context) + " contains a malformed real token."
    );
  }
}

[[nodiscard]] std::size_t resolve_obj_index(
  std::string_view token,
  std::size_t entry_count,
  std::string_view context,
  std::string_view label
)
{
  const auto trimmed = trim_copy(token);
  if(trimmed.empty()) {
    fail(
      base::StatusCode::io_error,
      std::string(context) + " contains an empty OBJ " + std::string(label) + " reference."
    );
  }

  const auto value = parse_integer(
    trimmed,
    std::string(context) + " " + std::string(label) + " reference"
  );
  if(value == 0) {
    fail(
      base::StatusCode::unsupported,
      std::string(context) + " does not support OBJ " + std::string(label) + " index 0."
    );
  }

  const auto resolved =
    value > 0
      ? value - 1
      : static_cast<long long>(entry_count) + value;
  if(resolved < 0 || static_cast<std::size_t>(resolved) >= entry_count) {
    fail(
      base::StatusCode::io_error,
      std::string(context) + " references a " + std::string(label) +
        " that has not been defined."
    );
  }

  return static_cast<std::size_t>(resolved);
}

[[nodiscard]] std::size_t parse_obj_vertex_ref(
  std::string_view token,
  std::size_t vertex_count,
  std::string_view context
)
{
  const auto slash = token.find('/');
  return resolve_obj_index(token.substr(0U, slash), vertex_count, context, "vertex");
}

[[nodiscard]] ObjCoordinate parse_obj_texture_coordinate(
  std::string_view payload
)
{
  const auto tokens = split_tokens(payload);
  if(tokens.empty()) {
    fail(base::StatusCode::io_error, "OBJ texture-coordinate record is malformed.");
  }
  if(tokens.size() > 3U) {
    fail(
      base::StatusCode::unsupported,
      "OBJ import currently supports only `vt u`, `vt u v`, or `vt u v w` texture-coordinate records."
    );
  }

  ObjCoordinate coordinate;
  coordinate.component_count = static_cast<std::uint8_t>(tokens.size());
  for(std::size_t index = 0U; index < tokens.size(); ++index) {
    coordinate.values[index] =
      parse_real_token(tokens[index], "OBJ texture-coordinate record");
  }
  return coordinate;
}

[[nodiscard]] std::array<double, 3> parse_obj_normal(
  std::string_view payload
)
{
  const auto tokens = split_tokens(payload);
  if(tokens.size() != 3U) {
    fail(base::StatusCode::io_error, "OBJ normal record is malformed.");
  }

  return {
    parse_real_token(tokens[0], "OBJ normal record"),
    parse_real_token(tokens[1], "OBJ normal record"),
    parse_real_token(tokens[2], "OBJ normal record"),
  };
}

[[nodiscard]] ObjCoordinate parse_obj_parameter_vertex(
  std::string_view payload
)
{
  const auto tokens = split_tokens(payload);
  if(tokens.empty()) {
    fail(base::StatusCode::io_error, "OBJ parameter-vertex record is malformed.");
  }
  if(tokens.size() > 3U) {
    fail(
      base::StatusCode::unsupported,
      "OBJ import currently supports only `vp u`, `vp u v`, or `vp u v w` parameter-vertex records."
    );
  }

  ObjCoordinate coordinate;
  coordinate.component_count = static_cast<std::uint8_t>(tokens.size());
  for(std::size_t index = 0U; index < tokens.size(); ++index) {
    coordinate.values[index] =
      parse_real_token(tokens[index], "OBJ parameter-vertex record");
  }
  return coordinate;
}

[[noreturn]] void fail_malformed_obj_face_token(
  std::string_view token
)
{
  fail(
    base::StatusCode::io_error,
    "OBJ face record contains malformed token `" + std::string(token) +
      "`; supported forms are `v`, `v/vt`, `v//vn`, and `v/vt/vn`."
  );
}

[[nodiscard]] ObjFaceVertex parse_obj_face_vertex(
  std::string_view token,
  std::size_t vertex_count,
  std::size_t texture_coordinate_count,
  std::size_t normal_count
)
{
  const auto first_slash = token.find('/');
  if(first_slash == std::string_view::npos) {
    ObjFaceVertex vertex;
    vertex.vertex_index = resolve_obj_index(token, vertex_count, "OBJ face record", "vertex");
    return vertex;
  }

  const auto second_slash = token.find('/', first_slash + 1U);
  if(second_slash == std::string_view::npos) {
    const auto texture_token = token.substr(first_slash + 1U);
    if(texture_token.empty()) {
      fail_malformed_obj_face_token(token);
    }

    ObjFaceVertex vertex;
    vertex.vertex_index = resolve_obj_index(
      token.substr(0U, first_slash),
      vertex_count,
      "OBJ face record",
      "vertex"
    );
    vertex.texture_index = {
      true,
      resolve_obj_index(
        texture_token,
        texture_coordinate_count,
        "OBJ face record",
        "texture coordinate"
      ),
    };
    return vertex;
  }

  if(token.find('/', second_slash + 1U) != std::string_view::npos) {
    fail_malformed_obj_face_token(token);
  }

  const auto texture_token =
    token.substr(first_slash + 1U, second_slash - first_slash - 1U);
  const auto normal_token = token.substr(second_slash + 1U);
  if(normal_token.empty()) {
    fail_malformed_obj_face_token(token);
  }

  ObjFaceVertex vertex;
  vertex.vertex_index = resolve_obj_index(
    token.substr(0U, first_slash),
    vertex_count,
    "OBJ face record",
    "vertex"
  );
  if(!texture_token.empty()) {
    vertex.texture_index = {
      true,
      resolve_obj_index(
        texture_token,
        texture_coordinate_count,
        "OBJ face record",
        "texture coordinate"
      ),
    };
  }
  vertex.normal_index = {
    true,
    resolve_obj_index(normal_token, normal_count, "OBJ face record", "normal"),
  };
  return vertex;
}

[[nodiscard]] std::string strip_obj_inline_comment(std::string_view value)
{
  const auto comment = value.find('#');
  if(comment == std::string_view::npos) {
    return std::string(value);
  }
  return std::string(value.substr(0U, comment));
}

[[nodiscard]] const ObjFaceVertex &find_obj_face_vertex(
  const std::vector<ObjFaceVertex> &vertices,
  std::size_t vertex_index
)
{
  // `prepare_obj_face_polygon()` rejects repeated vertex references inside one
  // face, so triangulated corners can recover parser-side metadata by vertex
  // id without ambiguously matching multiple original corners.
  const auto vertex_it =
    std::find_if(
      vertices.begin(),
      vertices.end(),
      [&](const ObjFaceVertex &vertex) {
        return vertex.vertex_index == vertex_index;
      }
    );
  if(vertex_it == vertices.end()) {
    fail(
      base::StatusCode::internal_error,
      "OBJ import lost per-corner face metadata during triangulation."
    );
  }
  return *vertex_it;
}

[[nodiscard]] ImportedMeshDescription parse_obj_description(
  std::string_view path,
  const ObjImportOptions &options
)
{
  validate_non_empty_path(path, "OBJ import requires a non-empty path.");

  const std::string file_path(path);
  std::ifstream input(file_path);
  if(!input) {
    fail(base::StatusCode::io_error, "OBJ import could not open the input file.");
  }

  ImportedMeshDescription description;
  description.name = default_domain_name(path);

  std::string current_object = description.name;
  std::string current_group;
  std::string current_material;
  std::string line;
  std::size_t line_number = 0U;

  while(std::getline(input, line)) {
    ++line_number;
    const auto trimmed = trim_copy(strip_obj_inline_comment(line));
    if(trimmed.empty() || trimmed[0] == '#') {
      continue;
    }

    const auto first_whitespace = trimmed.find_first_of(" \t");
    const auto keyword =
      std::string_view(trimmed).substr(
        0U,
        first_whitespace == std::string::npos ? trimmed.size() : first_whitespace
      );
    const auto payload =
      first_whitespace == std::string::npos
        ? std::string_view {}
        : std::string_view(trimmed).substr(first_whitespace + 1U);
    const auto payload_trimmed = trim_copy(payload);

    if(keyword == "v") {
      std::istringstream stream(payload_trimmed);
      std::array<double, 3> node {};
      if(!(stream >> node[0] >> node[1] >> node[2])) {
        fail(
          base::StatusCode::io_error,
          "OBJ vertex record is malformed."
        );
      }
      std::string extra;
      if(stream >> extra) {
        fail(
          base::StatusCode::unsupported,
          "OBJ import currently supports only `v x y z` vertex records."
        );
      }
      description.nodes.push_back(node);
      continue;
    }

    if(keyword == "vt") {
      description.obj.texture_coordinates.push_back(
        parse_obj_texture_coordinate(payload_trimmed)
      );
      continue;
    }

    if(keyword == "vn") {
      description.obj.normals.push_back(parse_obj_normal(payload_trimmed));
      continue;
    }

    if(keyword == "vp") {
      description.obj.parameter_vertices.push_back(
        parse_obj_parameter_vertex(payload_trimmed)
      );
      continue;
    }

    if(keyword == "mtllib") {
      const auto tokens = split_tokens(payload_trimmed);
      if(tokens.empty()) {
        fail(base::StatusCode::io_error, "OBJ mtllib record is missing a library name.");
      }

      ObjMaterialLibrary library;
      library.entries = tokens;
      description.obj.material_libraries.push_back(std::move(library));
      continue;
    }

    if(keyword == "usemtl") {
      if(payload_trimmed.empty()) {
        fail(base::StatusCode::io_error, "OBJ usemtl record is missing a material name.");
      }
      current_material = payload_trimmed;
      continue;
    }

    if(keyword == "s") {
      continue;
    }

    if(keyword == "o") {
      if(options.read_groups) {
        current_object = payload_trimmed;
        if(description.name.empty() || description.name == default_domain_name(path)) {
          if(!current_object.empty()) {
            description.name = current_object;
          }
        }
      }
      continue;
    }

    if(keyword == "g") {
      if(options.read_groups) {
        current_group = payload_trimmed;
      }
      continue;
    }

    if(keyword == "l") {
      const auto tokens = split_tokens(payload_trimmed);
      if(tokens.size() < 2U) {
        fail(
          base::StatusCode::unsupported,
          "OBJ import requires at least two vertices in an `l` record."
        );
      }

      std::vector<std::size_t> refs;
      refs.reserve(tokens.size());
      for(const auto &token_value : tokens) {
        refs.push_back(
          parse_obj_vertex_ref(token_value, description.nodes.size(), "OBJ line record")
        );
      }

      for(std::size_t index = 1U; index < refs.size(); ++index) {
        ImportedEdge edge;
        edge.nodes = {refs[index - 1U], refs[index]};
        edge.material_name = current_material;
        if(options.read_groups) {
          edge.group.name = current_group.empty() ? current_object : current_group;
        }
        description.edges.push_back(std::move(edge));
      }
      continue;
    }

    if(keyword == "f") {
      const auto tokens = split_tokens(payload_trimmed);
      if(tokens.size() < 3U) {
        fail(
          base::StatusCode::unsupported,
          "OBJ import requires at least three vertices in an `f` record."
        );
      }

      std::vector<ObjFaceVertex> face_vertices;
      face_vertices.reserve(tokens.size());
      std::vector<std::size_t> refs;
      refs.reserve(tokens.size());
      for(const auto &token_value : tokens) {
        const auto face_vertex =
          parse_obj_face_vertex(
            token_value,
            description.nodes.size(),
            description.obj.texture_coordinates.size(),
            description.obj.normals.size()
          );
        refs.push_back(face_vertex.vertex_index);
        face_vertices.push_back(face_vertex);
      }

      const auto triangles =
        triangulate_obj_face_refs(std::move(refs), description.nodes, line_number);

      for(const auto &triangle : triangles) {
        ImportedFace face;
        face.nodes = triangle;
        face.obj_face_vertices = {
          find_obj_face_vertex(face_vertices, triangle[0]),
          find_obj_face_vertex(face_vertices, triangle[1]),
          find_obj_face_vertex(face_vertices, triangle[2]),
        };
        face.material_name = current_material;
        if(options.read_groups) {
          face.group.name = current_group.empty() ? current_object : current_group;
        }
        description.faces.push_back(std::move(face));
      }
      continue;
    }

    fail(
      base::StatusCode::unsupported,
      "OBJ import encountered unsupported record `" + std::string(keyword) +
      "` on line " + std::to_string(line_number) + "."
    );
  }

  if(description.nodes.empty()) {
    fail(base::StatusCode::io_error, "OBJ import found no vertices.");
  }

  return description;
}

[[nodiscard]] testing::ObjReviewOptionalIndex make_obj_review_optional_index(
  const ObjOptionalIndex &index
) noexcept
{
  return {index.present, index.value};
}

[[nodiscard]] testing::ObjReviewCoordinate make_obj_review_coordinate(
  const ObjCoordinate &coordinate
) noexcept
{
  return {coordinate.values, coordinate.component_count};
}

[[nodiscard]] testing::ObjReviewFaceCorner make_obj_review_face_corner(
  const ObjFaceVertex &vertex
) noexcept
{
  return {
    vertex.vertex_index,
    make_obj_review_optional_index(vertex.texture_index),
    make_obj_review_optional_index(vertex.normal_index),
  };
}

[[nodiscard]] testing::ObjImportReview make_obj_import_review(
  const ImportedMeshDescription &description
)
{
  testing::ObjImportReview review;
  review.name = description.name;

  review.texture_coordinates.reserve(description.obj.texture_coordinates.size());
  for(const auto &coordinate : description.obj.texture_coordinates) {
    review.texture_coordinates.push_back(make_obj_review_coordinate(coordinate));
  }

  review.normals = description.obj.normals;

  review.parameter_vertices.reserve(description.obj.parameter_vertices.size());
  for(const auto &coordinate : description.obj.parameter_vertices) {
    review.parameter_vertices.push_back(make_obj_review_coordinate(coordinate));
  }

  review.material_libraries.reserve(description.obj.material_libraries.size());
  for(const auto &library : description.obj.material_libraries) {
    testing::ObjReviewMaterialLibrary review_library;
    review_library.entries = library.entries;
    review.material_libraries.push_back(std::move(review_library));
  }

  review.edges.reserve(description.edges.size());
  for(const auto &edge : description.edges) {
    testing::ObjReviewEdge review_edge;
    review_edge.nodes = edge.nodes;
    review_edge.material_name = edge.material_name;
    review.edges.push_back(std::move(review_edge));
  }

  review.faces.reserve(description.faces.size());
  for(const auto &face : description.faces) {
    testing::ObjReviewFace review_face;
    review_face.material_name = face.material_name;
    for(std::size_t index = 0U; index < face.obj_face_vertices.size(); ++index) {
      review_face.corners[index] =
        make_obj_review_face_corner(face.obj_face_vertices[index]);
    }
    review.faces.push_back(std::move(review_face));
  }

  return review;
}

void expect_msh_end(
  std::istream &input,
  const char *expected_end
)
{
  std::string line;
  if(!std::getline(input >> std::ws, line)) {
    fail(
      base::StatusCode::io_error,
      std::string("MSH import reached EOF before ") + expected_end + '.'
    );
  }
  if(trim_copy(line) != expected_end) {
    fail(
      base::StatusCode::io_error,
      std::string("MSH import expected ") + expected_end + '.'
    );
  }
}

void skip_msh_section(
  std::istream &input,
  const std::string &section
)
{
  const auto end_token = "$End" + section.substr(1U);
  std::string line;
  while(std::getline(input, line)) {
    if(trim_copy(line) == end_token) {
      return;
    }
  }

  fail(
    base::StatusCode::io_error,
    "MSH import reached EOF while skipping section " + section + '.'
  );
}

[[nodiscard]] GroupInfo msh_group_from_entity(
  int dimension,
  int entity_tag,
  const std::vector<int> &physical_tags,
  const std::map<std::pair<int, int>, std::string> &physical_names
)
{
  if(!physical_tags.empty() && physical_tags.front() > 0) {
    GroupInfo group;
    group.zone_id = static_cast<std::uint32_t>(physical_tags.front());
    group.source_entity_tag =
      entity_tag > 0 ? static_cast<std::uint32_t>(entity_tag) : invalid_index;
    const auto name_it = physical_names.find({dimension, physical_tags.front()});
    group.name =
      name_it == physical_names.end()
        ? "physical_" + std::to_string(physical_tags.front())
        : name_it->second;
    return group;
  }

  if(entity_tag > 0) {
    GroupInfo group;
    group.zone_id = static_cast<std::uint32_t>(entity_tag);
    group.source_entity_tag = static_cast<std::uint32_t>(entity_tag);
    group.name = "entity_" + std::to_string(entity_tag);
    return group;
  }

  return {};
}

[[nodiscard]] GroupInfo msh_group(
  int dimension,
  const std::vector<long long> &tags,
  const std::map<std::pair<int, int>, std::string> &physical_names
)
{
  std::vector<int> physical_tags;
  if(!tags.empty() && tags[0] > 0) {
    physical_tags.push_back(static_cast<int>(tags[0]));
  }

  const auto entity_tag =
    tags.size() >= 2U && tags[1] > 0 ? static_cast<int>(tags[1]) : 0;
  return msh_group_from_entity(dimension, entity_tag, physical_tags, physical_names);
}

template <std::size_t N>
[[nodiscard]] std::array<std::size_t, N> map_msh_nodes(
  const std::array<long long, N> &node_ids,
  const std::map<long long, std::size_t> &node_map,
  std::string_view record_name
)
{
  std::array<std::size_t, N> nodes {};
  for(std::size_t offset = 0; offset < N; ++offset) {
    const auto node_it = node_map.find(node_ids[offset]);
    if(node_it == node_map.end()) {
      fail(
        base::StatusCode::io_error,
        std::string(record_name) + " references an undefined node id."
      );
    }
    nodes[offset] = node_it->second;
  }
  return nodes;
}

void parse_msh41_entities(
  std::istream &input,
  std::map<std::pair<int, int>, MshEntityInfo> &entity_info
)
{
  std::size_t point_count = 0U;
  std::size_t curve_count = 0U;
  std::size_t surface_count = 0U;
  std::size_t volume_count = 0U;
  if(!(input >> point_count >> curve_count >> surface_count >> volume_count)) {
    fail(base::StatusCode::io_error, "MSH import found a malformed $Entities section.");
  }

  const std::array<std::size_t, 4> counts = {
    point_count,
    curve_count,
    surface_count,
    volume_count,
  };

  for(int dimension = 0; dimension <= 3; ++dimension) {
    for(std::size_t index = 0U; index < counts[static_cast<std::size_t>(dimension)]; ++index) {
      int tag = 0;
      if(!(input >> tag)) {
        fail(base::StatusCode::io_error, "MSH import found an incomplete entity record.");
      }

      double min_x = 0.0;
      double min_y = 0.0;
      double min_z = 0.0;
      double max_x = 0.0;
      double max_y = 0.0;
      double max_z = 0.0;
      if(dimension == 0) {
        if(!(input >> min_x >> min_y >> min_z)) {
          fail(base::StatusCode::io_error, "MSH import found a malformed point entity record.");
        }
      }
      else if(!(input >> min_x >> min_y >> min_z >> max_x >> max_y >> max_z)) {
        fail(base::StatusCode::io_error, "MSH import found a malformed entity bounding box.");
      }
      static_cast<void>(min_x);
      static_cast<void>(min_y);
      static_cast<void>(min_z);
      static_cast<void>(max_x);
      static_cast<void>(max_y);
      static_cast<void>(max_z);

      std::size_t physical_count = 0U;
      if(!(input >> physical_count)) {
        fail(base::StatusCode::io_error, "MSH import found a malformed entity physical-tag list.");
      }

      MshEntityInfo info;
      info.physical_tags.reserve(physical_count);
      for(std::size_t physical_index = 0U; physical_index < physical_count; ++physical_index) {
        int physical_tag = 0;
        if(!(input >> physical_tag)) {
          fail(base::StatusCode::io_error, "MSH import found an incomplete entity physical-tag list.");
        }
        info.physical_tags.push_back(physical_tag);
      }

      if(dimension > 0) {
        std::size_t boundary_count = 0U;
        if(!(input >> boundary_count)) {
          fail(base::StatusCode::io_error, "MSH import found a malformed entity boundary list.");
        }
        for(std::size_t boundary_index = 0U; boundary_index < boundary_count; ++boundary_index) {
          int ignored_tag = 0;
          if(!(input >> ignored_tag)) {
            fail(base::StatusCode::io_error, "MSH import found an incomplete entity boundary list.");
          }
        }
      }

      if(!entity_info.emplace(std::make_pair(dimension, tag), std::move(info)).second) {
        fail(base::StatusCode::io_error, "MSH import found a duplicate entity tag.");
      }
    }
  }
}

void parse_msh41_entities_binary(
  std::istream &input,
  bool swap_bytes,
  int tag_data_size,
  std::map<std::pair<int, int>, MshEntityInfo> &entity_info
)
{
  std::array<std::uint64_t, 4> count_values = {0U, 0U, 0U, 0U};
  read_msh41_binary_size_values(
    input,
    count_values.data(),
    count_values.size(),
    swap_bytes,
    tag_data_size,
    "$Entities header"
  );
  const std::array<std::size_t, 4> counts = {
    checked_msh_binary_size(
      count_values[0],
      "MSH import found an entity count outside the supported range."
    ),
    checked_msh_binary_size(
      count_values[1],
      "MSH import found an entity count outside the supported range."
    ),
    checked_msh_binary_size(
      count_values[2],
      "MSH import found an entity count outside the supported range."
    ),
    checked_msh_binary_size(
      count_values[3],
      "MSH import found an entity count outside the supported range."
    ),
  };

  for(int dimension = 0; dimension <= 3; ++dimension) {
    for(std::size_t index = 0U; index < counts[static_cast<std::size_t>(dimension)]; ++index) {
      int tag = 0;
      read_msh_binary_value(input, tag, swap_bytes, "entity tag");

      std::array<double, 6> bounds = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      const auto bound_count = dimension == 0 ? 3U : 6U;
      read_msh_binary_values(
        input,
        bounds.data(),
        bound_count,
        swap_bytes,
        "entity bounding box"
      );

      std::uint64_t physical_count_value = 0U;
      read_msh41_binary_size_value(
        input,
        physical_count_value,
        swap_bytes,
        tag_data_size,
        "entity physical-tag count"
      );
      const auto physical_count = checked_msh_binary_size(
        physical_count_value,
        "MSH import found an entity physical-tag count outside the supported range."
      );

      MshEntityInfo info;
      info.physical_tags.resize(physical_count, 0);
      read_msh_binary_values(
        input,
        info.physical_tags.data(),
        info.physical_tags.size(),
        swap_bytes,
        "entity physical-tag list"
      );

      if(dimension > 0) {
        std::uint64_t boundary_count_value = 0U;
        read_msh41_binary_size_value(
          input,
          boundary_count_value,
          swap_bytes,
          tag_data_size,
          "entity boundary count"
        );
        const auto boundary_count = checked_msh_binary_size(
          boundary_count_value,
          "MSH import found an entity boundary count outside the supported range."
        );

        std::vector<int> ignored_tags(boundary_count, 0);
        read_msh_binary_values(
          input,
          ignored_tags.data(),
          ignored_tags.size(),
          swap_bytes,
          "entity boundary list"
        );
      }

      if(!entity_info.emplace(std::make_pair(dimension, tag), std::move(info)).second) {
        fail(base::StatusCode::io_error, "MSH import found a duplicate entity tag.");
      }
    }
  }
}

void parse_msh41_partitioned_entities(
  std::istream &input,
  std::map<std::pair<int, int>, MshEntityInfo> &entity_info
)
{
  std::size_t partition_count = 0U;
  std::size_t ghost_count = 0U;
  if(!(input >> partition_count >> ghost_count)) {
    fail(
      base::StatusCode::io_error,
      "MSH import found a malformed $PartitionedEntities header."
    );
  }
  static_cast<void>(partition_count);

  for(std::size_t ghost_index = 0U; ghost_index < ghost_count; ++ghost_index) {
    int ghost_entity_tag = 0;
    int partition_tag = 0;
    if(!(input >> ghost_entity_tag >> partition_tag)) {
      fail(
        base::StatusCode::io_error,
        "MSH import found an incomplete partitioned-ghost entity record."
      );
    }
  }

  std::size_t point_count = 0U;
  std::size_t curve_count = 0U;
  std::size_t surface_count = 0U;
  std::size_t volume_count = 0U;
  if(!(input >> point_count >> curve_count >> surface_count >> volume_count)) {
    fail(
      base::StatusCode::io_error,
      "MSH import found a malformed $PartitionedEntities entity-count header."
    );
  }

  const std::array<std::size_t, 4> counts = {
    point_count,
    curve_count,
    surface_count,
    volume_count,
  };

  for(int dimension = 0; dimension <= 3; ++dimension) {
    for(std::size_t index = 0U; index < counts[static_cast<std::size_t>(dimension)]; ++index) {
      int tag = 0;
      int parent_dimension = 0;
      int parent_tag = 0;
      std::size_t partition_tag_count = 0U;
      if(!(input >> tag >> parent_dimension >> parent_tag >> partition_tag_count)) {
        fail(
          base::StatusCode::io_error,
          "MSH import found a malformed partitioned entity record."
        );
      }
      static_cast<void>(parent_dimension);
      static_cast<void>(parent_tag);

      for(std::size_t partition_index = 0U;
          partition_index < partition_tag_count;
          ++partition_index) {
        int ignored_partition_tag = 0;
        if(!(input >> ignored_partition_tag)) {
          fail(
            base::StatusCode::io_error,
            "MSH import found an incomplete partitioned entity partition list."
          );
        }
      }

      double min_x = 0.0;
      double min_y = 0.0;
      double min_z = 0.0;
      double max_x = 0.0;
      double max_y = 0.0;
      double max_z = 0.0;
      if(dimension == 0) {
        if(!(input >> min_x >> min_y >> min_z)) {
          fail(
            base::StatusCode::io_error,
            "MSH import found a malformed partitioned point entity record."
          );
        }
      }
      else if(!(input >> min_x >> min_y >> min_z >> max_x >> max_y >> max_z)) {
        fail(
          base::StatusCode::io_error,
          "MSH import found a malformed partitioned entity bounding box."
        );
      }
      static_cast<void>(min_x);
      static_cast<void>(min_y);
      static_cast<void>(min_z);
      static_cast<void>(max_x);
      static_cast<void>(max_y);
      static_cast<void>(max_z);

      std::size_t physical_count = 0U;
      if(!(input >> physical_count)) {
        fail(
          base::StatusCode::io_error,
          "MSH import found a malformed partitioned entity physical-tag list."
        );
      }

      MshEntityInfo info;
      info.physical_tags.reserve(physical_count);
      for(std::size_t physical_index = 0U; physical_index < physical_count; ++physical_index) {
        int physical_tag = 0;
        if(!(input >> physical_tag)) {
          fail(
            base::StatusCode::io_error,
            "MSH import found an incomplete partitioned entity physical-tag list."
          );
        }
        info.physical_tags.push_back(physical_tag);
      }

      if(dimension > 0) {
        std::size_t boundary_count = 0U;
        if(!(input >> boundary_count)) {
          fail(
            base::StatusCode::io_error,
            "MSH import found a malformed partitioned entity boundary list."
          );
        }
        for(std::size_t boundary_index = 0U; boundary_index < boundary_count; ++boundary_index) {
          int ignored_tag = 0;
          if(!(input >> ignored_tag)) {
            fail(
              base::StatusCode::io_error,
              "MSH import found an incomplete partitioned entity boundary list."
            );
          }
        }
      }

      merge_msh_entity_info(entity_info, dimension, tag, std::move(info));
    }
  }
}

void parse_msh41_partitioned_entities_binary(
  std::istream &input,
  bool swap_bytes,
  int tag_data_size,
  std::map<std::pair<int, int>, MshEntityInfo> &entity_info
)
{
  std::uint64_t partition_count_value = 0U;
  std::uint64_t ghost_count_value = 0U;
  read_msh41_binary_size_value(
    input,
    partition_count_value,
    swap_bytes,
    tag_data_size,
    "$PartitionedEntities partition count"
  );
  read_msh41_binary_size_value(
    input,
    ghost_count_value,
    swap_bytes,
    tag_data_size,
    "$PartitionedEntities ghost-entity count"
  );
  static_cast<void>(partition_count_value);

  const auto ghost_count = checked_msh_binary_size(
    ghost_count_value,
    "MSH import found a partitioned-ghost entity count outside the supported range."
  );
  for(std::size_t ghost_index = 0U; ghost_index < ghost_count; ++ghost_index) {
    std::array<int, 2> ghost_record = {0, 0};
    read_msh_binary_values(
      input,
      ghost_record.data(),
      ghost_record.size(),
      swap_bytes,
      "partitioned-ghost entity record"
    );
  }

  std::array<std::uint64_t, 4> count_values = {0U, 0U, 0U, 0U};
  read_msh41_binary_size_values(
    input,
    count_values.data(),
    count_values.size(),
    swap_bytes,
    tag_data_size,
    "$PartitionedEntities entity-count header"
  );
  const std::array<std::size_t, 4> counts = {
    checked_msh_binary_size(
      count_values[0],
      "MSH import found a partitioned point-entity count outside the supported range."
    ),
    checked_msh_binary_size(
      count_values[1],
      "MSH import found a partitioned curve-entity count outside the supported range."
    ),
    checked_msh_binary_size(
      count_values[2],
      "MSH import found a partitioned surface-entity count outside the supported range."
    ),
    checked_msh_binary_size(
      count_values[3],
      "MSH import found a partitioned volume-entity count outside the supported range."
    ),
  };

  for(int dimension = 0; dimension <= 3; ++dimension) {
    for(std::size_t index = 0U; index < counts[static_cast<std::size_t>(dimension)]; ++index) {
      std::array<int, 3> record_header = {0, 0, 0};
      read_msh_binary_values(
        input,
        record_header.data(),
        record_header.size(),
        swap_bytes,
        "partitioned entity header"
      );
      const auto tag = record_header[0];

      std::uint64_t partition_tag_count_value = 0U;
      read_msh41_binary_size_value(
        input,
        partition_tag_count_value,
        swap_bytes,
        tag_data_size,
        "partitioned entity partition-tag count"
      );
      const auto partition_tag_count = checked_msh_binary_size(
        partition_tag_count_value,
        "MSH import found a partitioned entity partition count outside the supported range."
      );
      std::vector<int> ignored_partition_tags(partition_tag_count, 0);
      read_msh_binary_values(
        input,
        ignored_partition_tags.data(),
        ignored_partition_tags.size(),
        swap_bytes,
        "partitioned entity partition-tag list"
      );

      std::array<double, 6> bounds = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      const auto bound_count = dimension == 0 ? 3U : 6U;
      read_msh_binary_values(
        input,
        bounds.data(),
        bound_count,
        swap_bytes,
        "partitioned entity bounding box"
      );

      std::uint64_t physical_count_value = 0U;
      read_msh41_binary_size_value(
        input,
        physical_count_value,
        swap_bytes,
        tag_data_size,
        "partitioned entity physical-tag count"
      );
      const auto physical_count = checked_msh_binary_size(
        physical_count_value,
        "MSH import found a partitioned entity physical-tag count outside the supported range."
      );

      MshEntityInfo info;
      info.physical_tags.resize(physical_count, 0);
      read_msh_binary_values(
        input,
        info.physical_tags.data(),
        info.physical_tags.size(),
        swap_bytes,
        "partitioned entity physical-tag list"
      );

      if(dimension > 0) {
        std::uint64_t boundary_count_value = 0U;
        read_msh41_binary_size_value(
          input,
          boundary_count_value,
          swap_bytes,
          tag_data_size,
          "partitioned entity boundary count"
        );
        const auto boundary_count = checked_msh_binary_size(
          boundary_count_value,
          "MSH import found a partitioned entity boundary count outside the supported range."
        );

        std::vector<int> ignored_boundary_tags(boundary_count, 0);
        read_msh_binary_values(
          input,
          ignored_boundary_tags.data(),
          ignored_boundary_tags.size(),
          swap_bytes,
          "partitioned entity boundary list"
        );
      }

      merge_msh_entity_info(entity_info, dimension, tag, std::move(info));
    }
  }
}

void parse_msh41_periodic(
  std::istream &input
)
{
  std::size_t periodic_link_count = 0U;
  if(!(input >> periodic_link_count)) {
    fail(base::StatusCode::io_error, "MSH import found a malformed $Periodic section.");
  }

  for(std::size_t link_index = 0U; link_index < periodic_link_count; ++link_index) {
    int entity_dimension = 0;
    int entity_tag = 0;
    int master_entity_tag = 0;
    if(!(input >> entity_dimension >> entity_tag >> master_entity_tag)) {
      fail(base::StatusCode::io_error, "MSH import found a malformed periodic-link header.");
    }
    static_cast<void>(entity_dimension);
    static_cast<void>(entity_tag);
    static_cast<void>(master_entity_tag);

    std::size_t affine_value_count = 0U;
    if(!(input >> affine_value_count)) {
      fail(base::StatusCode::io_error, "MSH import found a malformed periodic affine-transform header.");
    }
    for(std::size_t affine_index = 0U; affine_index < affine_value_count; ++affine_index) {
      double ignored_value = 0.0;
      if(!(input >> ignored_value)) {
        fail(base::StatusCode::io_error, "MSH import found an incomplete periodic affine transform.");
      }
    }

    std::size_t corresponding_node_count = 0U;
    if(!(input >> corresponding_node_count)) {
      fail(base::StatusCode::io_error, "MSH import found a malformed periodic node-correspondence header.");
    }
    for(std::size_t node_index = 0U; node_index < corresponding_node_count; ++node_index) {
      std::size_t node_tag = 0U;
      std::size_t master_node_tag = 0U;
      if(!(input >> node_tag >> master_node_tag)) {
        fail(base::StatusCode::io_error, "MSH import found an incomplete periodic node correspondence.");
      }
    }
  }
}

void parse_msh41_periodic_binary(
  std::istream &input,
  bool swap_bytes,
  int tag_data_size
)
{
  std::uint64_t periodic_link_count_value = 0U;
  read_msh41_binary_size_value(
    input,
    periodic_link_count_value,
    swap_bytes,
    tag_data_size,
    "$Periodic link count"
  );
  const auto periodic_link_count = checked_msh_binary_size(
    periodic_link_count_value,
    "MSH import found a periodic-link count outside the supported range."
  );

  for(std::size_t link_index = 0U; link_index < periodic_link_count; ++link_index) {
    std::array<int, 3> link_header = {0, 0, 0};
    read_msh_binary_values(
      input,
      link_header.data(),
      link_header.size(),
      swap_bytes,
      "periodic-link header"
    );

    std::uint64_t affine_value_count_value = 0U;
    read_msh41_binary_size_value(
      input,
      affine_value_count_value,
      swap_bytes,
      tag_data_size,
      "periodic affine-transform count"
    );
    const auto affine_value_count = checked_msh_binary_size(
      affine_value_count_value,
      "MSH import found a periodic affine-transform count outside the supported range."
    );
    std::vector<double> ignored_affine_values(affine_value_count, 0.0);
    read_msh_binary_values(
      input,
      ignored_affine_values.data(),
      ignored_affine_values.size(),
      swap_bytes,
      "periodic affine transform"
    );

    std::uint64_t corresponding_node_count_value = 0U;
    read_msh41_binary_size_value(
      input,
      corresponding_node_count_value,
      swap_bytes,
      tag_data_size,
      "periodic node-correspondence count"
    );
    const auto corresponding_node_count = checked_msh_binary_size(
      corresponding_node_count_value,
      "MSH import found a periodic node-correspondence count outside the supported range."
    );
    std::vector<std::uint64_t> ignored_node_pairs(corresponding_node_count * 2U, 0U);
    read_msh41_binary_size_values(
      input,
      ignored_node_pairs.data(),
      ignored_node_pairs.size(),
      swap_bytes,
      tag_data_size,
      "periodic node-correspondence list"
    );
  }
}

void parse_msh41_ghost_elements(
  std::istream &input
)
{
  std::size_t ghost_element_count = 0U;
  if(!(input >> ghost_element_count)) {
    fail(base::StatusCode::io_error, "MSH import found a malformed $GhostElements section.");
  }

  for(std::size_t ghost_index = 0U; ghost_index < ghost_element_count; ++ghost_index) {
    std::size_t element_tag = 0U;
    int partition_tag = 0;
    std::size_t ghost_partition_count = 0U;
    if(!(input >> element_tag >> partition_tag >> ghost_partition_count)) {
      fail(base::StatusCode::io_error, "MSH import found a malformed ghost-element record.");
    }
    static_cast<void>(element_tag);
    static_cast<void>(partition_tag);

    for(std::size_t partition_index = 0U;
        partition_index < ghost_partition_count;
        ++partition_index) {
      int ignored_ghost_partition = 0;
      if(!(input >> ignored_ghost_partition)) {
        fail(base::StatusCode::io_error, "MSH import found an incomplete ghost-partition list.");
      }
    }
  }
}

void parse_msh41_ghost_elements_binary(
  std::istream &input,
  bool swap_bytes,
  int tag_data_size
)
{
  std::uint64_t ghost_element_count_value = 0U;
  read_msh41_binary_size_value(
    input,
    ghost_element_count_value,
    swap_bytes,
    tag_data_size,
    "$GhostElements count"
  );
  const auto ghost_element_count = checked_msh_binary_size(
    ghost_element_count_value,
    "MSH import found a ghost-element count outside the supported range."
  );

  for(std::size_t ghost_index = 0U; ghost_index < ghost_element_count; ++ghost_index) {
    std::uint64_t element_tag_value = 0U;
    read_msh41_binary_size_value(
      input,
      element_tag_value,
      swap_bytes,
      tag_data_size,
      "ghost element tag"
    );
    int partition_tag = 0;
    read_msh_binary_value(
      input,
      partition_tag,
      swap_bytes,
      "ghost element partition tag"
    );

    std::uint64_t ghost_partition_count_value = 0U;
    read_msh41_binary_size_value(
      input,
      ghost_partition_count_value,
      swap_bytes,
      tag_data_size,
      "ghost-partition count"
    );
    const auto ghost_partition_count = checked_msh_binary_size(
      ghost_partition_count_value,
      "MSH import found a ghost-partition count outside the supported range."
    );
    std::vector<int> ignored_ghost_partitions(ghost_partition_count, 0);
    read_msh_binary_values(
      input,
      ignored_ghost_partitions.data(),
      ignored_ghost_partitions.size(),
      swap_bytes,
      "ghost-partition list"
    );
  }
}

void parse_msh41_parametrizations(
  std::istream &input
)
{
  std::size_t curve_param_count = 0U;
  std::size_t surface_param_count = 0U;
  if(!(input >> curve_param_count >> surface_param_count)) {
    fail(base::StatusCode::io_error, "MSH import found a malformed $Parametrizations header.");
  }

  for(std::size_t curve_index = 0U; curve_index < curve_param_count; ++curve_index) {
    int curve_tag = 0;
    std::size_t node_count = 0U;
    if(!(input >> curve_tag >> node_count)) {
      fail(base::StatusCode::io_error, "MSH import found a malformed curve parametrization header.");
    }
    static_cast<void>(curve_tag);

    for(std::size_t node_index = 0U; node_index < node_count; ++node_index) {
      double x = 0.0;
      double y = 0.0;
      double z = 0.0;
      double u = 0.0;
      if(!(input >> x >> y >> z >> u)) {
        fail(base::StatusCode::io_error, "MSH import found an incomplete curve parametrization node.");
      }
    }
  }

  for(std::size_t surface_index = 0U; surface_index < surface_param_count; ++surface_index) {
    int surface_tag = 0;
    std::size_t node_count = 0U;
    std::size_t triangle_count = 0U;
    if(!(input >> surface_tag >> node_count >> triangle_count)) {
      fail(base::StatusCode::io_error, "MSH import found a malformed surface parametrization header.");
    }
    static_cast<void>(surface_tag);

    for(std::size_t node_index = 0U; node_index < node_count; ++node_index) {
      std::array<double, 11> ignored_node_record = {
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      };
      for(auto &value : ignored_node_record) {
        if(!(input >> value)) {
          fail(base::StatusCode::io_error, "MSH import found an incomplete surface parametrization node.");
        }
      }
    }

    for(std::size_t triangle_index = 0U; triangle_index < triangle_count; ++triangle_index) {
      std::array<int, 3> triangle = {0, 0, 0};
      if(!(input >> triangle[0] >> triangle[1] >> triangle[2])) {
        fail(base::StatusCode::io_error, "MSH import found an incomplete surface parametrization triangle.");
      }
    }
  }
}

void parse_msh41_parametrizations_binary(
  std::istream &input,
  bool swap_bytes,
  int tag_data_size
)
{
  std::array<std::uint64_t, 2> header = {0U, 0U};
  read_msh41_binary_size_values(
    input,
    header.data(),
    header.size(),
    swap_bytes,
    tag_data_size,
    "$Parametrizations header"
  );
  const auto curve_param_count = checked_msh_binary_size(
    header[0],
    "MSH import found a curve-parametrization count outside the supported range."
  );
  const auto surface_param_count = checked_msh_binary_size(
    header[1],
    "MSH import found a surface-parametrization count outside the supported range."
  );

  for(std::size_t curve_index = 0U; curve_index < curve_param_count; ++curve_index) {
    int curve_tag = 0;
    read_msh_binary_value(input, curve_tag, swap_bytes, "curve parametrization tag");

    std::uint64_t node_count_value = 0U;
    read_msh41_binary_size_value(
      input,
      node_count_value,
      swap_bytes,
      tag_data_size,
      "curve parametrization node count"
    );
    const auto node_count = checked_msh_binary_size(
      node_count_value,
      "MSH import found a curve-parametrization node count outside the supported range."
    );

    std::vector<double> ignored_node_records(node_count * 4U, 0.0);
    read_msh_binary_values(
      input,
      ignored_node_records.data(),
      ignored_node_records.size(),
      swap_bytes,
      "curve parametrization nodes"
    );
  }

  for(std::size_t surface_index = 0U; surface_index < surface_param_count; ++surface_index) {
    int surface_tag = 0;
    read_msh_binary_value(input, surface_tag, swap_bytes, "surface parametrization tag");

    std::array<std::uint64_t, 2> counts = {0U, 0U};
    read_msh41_binary_size_values(
      input,
      counts.data(),
      counts.size(),
      swap_bytes,
      tag_data_size,
      "surface parametrization counts"
    );
    const auto node_count = checked_msh_binary_size(
      counts[0],
      "MSH import found a surface-parametrization node count outside the supported range."
    );
    const auto triangle_count = checked_msh_binary_size(
      counts[1],
      "MSH import found a surface-parametrization triangle count outside the supported range."
    );

    std::vector<double> ignored_node_records(node_count * 11U, 0.0);
    read_msh_binary_values(
      input,
      ignored_node_records.data(),
      ignored_node_records.size(),
      swap_bytes,
      "surface parametrization nodes"
    );

    std::vector<int> ignored_triangles(triangle_count * 3U, 0);
    read_msh_binary_values(
      input,
      ignored_triangles.data(),
      ignored_triangles.size(),
      swap_bytes,
      "surface parametrization triangles"
    );
  }
}

void parse_msh22_nodes(
  std::istream &input,
  ImportedMeshDescription &description,
  std::map<long long, std::size_t> &node_map
)
{
  std::string count_line;
  if(!std::getline(input, count_line)) {
    fail(base::StatusCode::io_error, "MSH import found an incomplete $Nodes section.");
  }
  const auto count = parse_integer(trim_copy(count_line), "MSH Nodes");
  description.nodes.reserve(description.nodes.size() + static_cast<std::size_t>(count));
  for(long long index = 0; index < count; ++index) {
    std::string node_line;
    if(!std::getline(input, node_line)) {
      fail(base::StatusCode::io_error, "MSH import found an incomplete node record.");
    }

    std::istringstream stream(node_line);
    long long node_id = 0;
    std::array<double, 3> node {};
    if(!(stream >> node_id >> node[0] >> node[1] >> node[2])) {
      fail(base::StatusCode::io_error, "MSH import found a malformed node record.");
    }

    if(node_map.find(node_id) != node_map.end()) {
      fail(base::StatusCode::io_error, "MSH import found a duplicate node id.");
    }

    validate_msh_node_coordinates(node, node_id);

    node_map.emplace(node_id, description.nodes.size());
    description.nodes.push_back(node);
  }
}

void parse_msh22_nodes_binary(
  std::istream &input,
  ImportedMeshDescription &description,
  std::map<long long, std::size_t> &node_map,
  bool swap_bytes
)
{
  std::string count_line;
  if(!std::getline(input, count_line)) {
    fail(base::StatusCode::io_error, "MSH import found an incomplete $Nodes section.");
  }

  const auto count = parse_integer(trim_copy(count_line), "MSH Nodes");
  description.nodes.reserve(description.nodes.size() + static_cast<std::size_t>(count));
  for(long long index = 0; index < count; ++index) {
    int node_id = 0;
    std::array<double, 3> node = {0.0, 0.0, 0.0};
    read_msh_binary_value(input, node_id, swap_bytes, "node id");
    read_msh_binary_values(input, node.data(), node.size(), swap_bytes, "node coordinates");

    if(node_map.find(node_id) != node_map.end()) {
      fail(base::StatusCode::io_error, "MSH import found a duplicate node id.");
    }

    validate_msh_node_coordinates(node, node_id);

    node_map.emplace(node_id, description.nodes.size());
    description.nodes.push_back(node);
  }
}

void parse_msh41_nodes(
  std::istream &input,
  ImportedMeshDescription &description,
  std::map<long long, std::size_t> &node_map
)
{
  std::size_t block_count = 0U;
  std::size_t total_count = 0U;
  std::size_t min_tag = 0U;
  std::size_t max_tag = 0U;
  if(!(input >> block_count >> total_count >> min_tag >> max_tag)) {
    fail(base::StatusCode::io_error, "MSH import found a malformed $Nodes header.");
  }
  static_cast<void>(min_tag);
  static_cast<void>(max_tag);

  description.nodes.reserve(description.nodes.size() + total_count);
  for(std::size_t block_index = 0U; block_index < block_count; ++block_index) {
    int entity_dim = 0;
    int entity_tag = 0;
    int parametric = 0;
    std::size_t block_size = 0U;
    if(!(input >> entity_dim >> entity_tag >> parametric >> block_size)) {
      fail(base::StatusCode::io_error, "MSH import found a malformed node block header.");
    }
    static_cast<void>(entity_tag);

    std::vector<long long> node_ids(block_size, 0);
    for(std::size_t node_index = 0U; node_index < block_size; ++node_index) {
      std::size_t node_id = 0U;
      if(!(input >> node_id)) {
        fail(base::StatusCode::io_error, "MSH import found an incomplete node-tag block.");
      }
      node_ids[node_index] = static_cast<long long>(node_id);
    }

    const auto parametric_count =
      msh41_parametric_coordinate_count(entity_dim, parametric);

    for(std::size_t node_index = 0U; node_index < block_size; ++node_index) {
      std::array<double, 3> node {};
      if(!(input >> node[0] >> node[1] >> node[2])) {
        fail(base::StatusCode::io_error, "MSH import found a malformed node coordinate record.");
      }
      for(std::size_t param_index = 0U; param_index < parametric_count; ++param_index) {
        double ignored = 0.0;
        if(!(input >> ignored)) {
          fail(base::StatusCode::io_error, "MSH import found an incomplete parametric node record.");
        }
      }

      const auto node_id = node_ids[node_index];
      if(node_map.find(node_id) != node_map.end()) {
        fail(base::StatusCode::io_error, "MSH import found a duplicate node id.");
      }

      validate_msh_node_coordinates(node, node_id);

      node_map.emplace(node_id, description.nodes.size());
      description.nodes.push_back(node);
    }
  }
}

void parse_msh41_nodes_binary(
  std::istream &input,
  ImportedMeshDescription &description,
  std::map<long long, std::size_t> &node_map,
  bool swap_bytes,
  int tag_data_size
)
{
  std::array<std::uint64_t, 4> header = {0U, 0U, 0U, 0U};
  read_msh41_binary_size_values(
    input,
    header.data(),
    header.size(),
    swap_bytes,
    tag_data_size,
    "$Nodes header"
  );

  const auto block_count = checked_msh_binary_size(
    header[0],
    "MSH import found a node-block count outside the supported range."
  );
  const auto total_count = checked_msh_binary_size(
    header[1],
    "MSH import found a node count outside the supported range."
  );

  description.nodes.reserve(description.nodes.size() + total_count);
  for(std::size_t block_index = 0U; block_index < block_count; ++block_index) {
    std::array<int, 3> block_header = {0, 0, 0};
    read_msh_binary_values(
      input,
      block_header.data(),
      block_header.size(),
      swap_bytes,
      "node block header"
    );

    const auto entity_dim = block_header[0];
    const auto parametric = block_header[2];

    std::uint64_t block_size_value = 0U;
    read_msh41_binary_size_value(
      input,
      block_size_value,
      swap_bytes,
      tag_data_size,
      "node block size"
    );
    const auto block_size = checked_msh_binary_size(
      block_size_value,
      "MSH import found a node block size outside the supported range."
    );

    std::vector<std::uint64_t> node_ids(block_size, 0U);
    read_msh41_binary_size_values(
      input,
      node_ids.data(),
      node_ids.size(),
      swap_bytes,
      tag_data_size,
      "node-tag block"
    );

    const auto parametric_count =
      msh41_parametric_coordinate_count(entity_dim, parametric);

    for(std::size_t node_index = 0U; node_index < block_size; ++node_index) {
      std::array<double, 3> node = {0.0, 0.0, 0.0};
      read_msh_binary_values(
        input,
        node.data(),
        node.size(),
        swap_bytes,
        "node coordinate record"
      );

      std::vector<double> ignored_parametric(parametric_count, 0.0);
      read_msh_binary_values(
        input,
        ignored_parametric.data(),
        ignored_parametric.size(),
        swap_bytes,
        "parametric node record"
      );

      const auto node_id = checked_msh_binary_tag(
        node_ids[node_index],
        "MSH import found a node tag outside the supported integer range."
      );
      if(node_map.find(node_id) != node_map.end()) {
        fail(base::StatusCode::io_error, "MSH import found a duplicate node id.");
      }

      validate_msh_node_coordinates(node, node_id);

      node_map.emplace(node_id, description.nodes.size());
      description.nodes.push_back(node);
    }
  }
}

void parse_msh22_elements(
  std::istream &input,
  ImportedMeshDescription &description,
  const std::map<long long, std::size_t> &node_map,
  const std::map<std::pair<int, int>, std::string> &physical_names
)
{
  std::string count_line;
  if(!std::getline(input, count_line)) {
    fail(base::StatusCode::io_error, "MSH import found an incomplete $Elements section.");
  }
  const auto count = parse_integer(trim_copy(count_line), "MSH Elements");
  for(long long index = 0; index < count; ++index) {
    std::string element_line;
    if(!std::getline(input, element_line)) {
      fail(base::StatusCode::io_error, "MSH import found an incomplete element record.");
    }

    std::istringstream stream(element_line);
    long long element_id = 0;
    long long element_type = 0;
    long long tag_count = 0;
    if(!(stream >> element_id >> element_type >> tag_count)) {
      fail(base::StatusCode::io_error, "MSH import found a malformed element record.");
    }
    static_cast<void>(element_id);

    std::vector<long long> tags(static_cast<std::size_t>(tag_count), 0);
    for(long long tag_index = 0; tag_index < tag_count; ++tag_index) {
      if(!(stream >> tags[static_cast<std::size_t>(tag_index)])) {
        fail(base::StatusCode::io_error, "MSH import found a malformed element tag list.");
      }
    }

    switch(element_type) {
    case 1: {
      std::array<long long, 2> node_ids {};
      if(!(stream >> node_ids[0] >> node_ids[1])) {
        fail(base::StatusCode::io_error, "MSH import found a malformed line element.");
      }
      description.edges.push_back({
        map_msh_nodes(node_ids, node_map, "MSH line element"),
        msh_group(1, tags, physical_names),
      });
    } break;
    case 2: {
      std::array<long long, 3> node_ids {};
      if(!(stream >> node_ids[0] >> node_ids[1] >> node_ids[2])) {
        fail(base::StatusCode::io_error, "MSH import found a malformed triangle element.");
      }
      append_imported_face(
        description.faces,
        map_msh_nodes(node_ids, node_map, "MSH triangle element"),
        msh_group(2, tags, physical_names)
      );
    } break;
    case 4: {
      std::array<long long, 4> node_ids {};
      if(!(stream >> node_ids[0] >> node_ids[1] >> node_ids[2] >> node_ids[3])) {
        fail(base::StatusCode::io_error, "MSH import found a malformed tetrahedron element.");
      }
      description.cells.push_back({
        map_msh_nodes(node_ids, node_map, "MSH tetrahedron element"),
        msh_group(3, tags, physical_names),
      });
    } break;
    default:
      fail_unsupported_msh_element_type(element_type);
    }

    std::string extra;
    if(stream >> extra) {
      fail(base::StatusCode::io_error, "MSH import found extra tokens in an element record.");
    }
  }
}

void parse_msh22_elements_binary(
  std::istream &input,
  ImportedMeshDescription &description,
  const std::map<long long, std::size_t> &node_map,
  const std::map<std::pair<int, int>, std::string> &physical_names,
  bool swap_bytes
)
{
  std::string count_line;
  if(!std::getline(input, count_line)) {
    fail(base::StatusCode::io_error, "MSH import found an incomplete $Elements section.");
  }

  const auto count = parse_integer(trim_copy(count_line), "MSH Elements");
  long long parsed_count = 0;
  while(parsed_count < count) {
    std::array<int, 3> header = {0, 0, 0};
    read_msh_binary_values(
      input,
      header.data(),
      header.size(),
      swap_bytes,
      "binary element-block header"
    );

    const auto element_type = header[0];
    const auto block_size = header[1];
    const auto tag_count = header[2];
    if(block_size < 0 || tag_count < 0) {
      fail(base::StatusCode::io_error, "MSH import found a malformed binary element-block header.");
    }

    parsed_count += block_size;
    for(int element_index = 0; element_index < block_size; ++element_index) {
      switch(element_type) {
      case 1: {
        std::vector<int> data(static_cast<std::size_t>(1 + tag_count + 2), 0);
        read_msh_binary_values(input, data.data(), data.size(), swap_bytes, "line element record");
        std::vector<long long> tags(static_cast<std::size_t>(tag_count), 0);
        for(int tag_index = 0; tag_index < tag_count; ++tag_index) {
          tags[static_cast<std::size_t>(tag_index)] =
            data[1 + static_cast<std::size_t>(tag_index)];
        }
        const auto node_offset = static_cast<std::size_t>(1 + tag_count);
        description.edges.push_back({
          map_msh_nodes(
            std::array<long long, 2> {
              data[node_offset],
              data[node_offset + 1U],
            },
            node_map,
            "MSH line element"
          ),
          msh_group(1, tags, physical_names),
        });
      } break;
      case 2: {
        std::vector<int> data(static_cast<std::size_t>(1 + tag_count + 3), 0);
        read_msh_binary_values(
          input,
          data.data(),
          data.size(),
          swap_bytes,
          "triangle element record"
        );
        std::vector<long long> tags(static_cast<std::size_t>(tag_count), 0);
        for(int tag_index = 0; tag_index < tag_count; ++tag_index) {
          tags[static_cast<std::size_t>(tag_index)] =
            data[1 + static_cast<std::size_t>(tag_index)];
        }
        const auto node_offset = static_cast<std::size_t>(1 + tag_count);
        append_imported_face(
          description.faces,
          map_msh_nodes(
            std::array<long long, 3> {
              data[node_offset],
              data[node_offset + 1U],
              data[node_offset + 2U],
            },
            node_map,
            "MSH triangle element"
          ),
          msh_group(2, tags, physical_names)
        );
      } break;
      case 4: {
        std::vector<int> data(static_cast<std::size_t>(1 + tag_count + 4), 0);
        read_msh_binary_values(
          input,
          data.data(),
          data.size(),
          swap_bytes,
          "tetrahedron element record"
        );
        std::vector<long long> tags(static_cast<std::size_t>(tag_count), 0);
        for(int tag_index = 0; tag_index < tag_count; ++tag_index) {
          tags[static_cast<std::size_t>(tag_index)] =
            data[1 + static_cast<std::size_t>(tag_index)];
        }
        const auto node_offset = static_cast<std::size_t>(1 + tag_count);
        description.cells.push_back({
          map_msh_nodes(
            std::array<long long, 4> {
              data[node_offset],
              data[node_offset + 1U],
              data[node_offset + 2U],
              data[node_offset + 3U],
            },
            node_map,
            "MSH tetrahedron element"
          ),
          msh_group(3, tags, physical_names),
        });
      } break;
      default:
        fail_unsupported_msh_element_type(element_type);
      }
    }
  }

  if(parsed_count != count) {
    fail(base::StatusCode::io_error, "MSH import found a malformed binary $Elements section.");
  }
}

void parse_msh41_elements(
  std::istream &input,
  ImportedMeshDescription &description,
  const std::map<long long, std::size_t> &node_map,
  const std::map<std::pair<int, int>, MshEntityInfo> &entity_info,
  const std::map<std::pair<int, int>, std::string> &physical_names
)
{
  std::size_t block_count = 0U;
  std::size_t total_count = 0U;
  std::size_t min_tag = 0U;
  std::size_t max_tag = 0U;
  if(!(input >> block_count >> total_count >> min_tag >> max_tag)) {
    fail(base::StatusCode::io_error, "MSH import found a malformed $Elements header.");
  }
  static_cast<void>(total_count);
  static_cast<void>(min_tag);
  static_cast<void>(max_tag);

  for(std::size_t block_index = 0U; block_index < block_count; ++block_index) {
    int entity_dim = 0;
    int entity_tag = 0;
    int element_type = 0;
    std::size_t block_size = 0U;
    if(!(input >> entity_dim >> entity_tag >> element_type >> block_size)) {
      fail(base::StatusCode::io_error, "MSH import found a malformed element block header.");
    }

    const auto entity_it = entity_info.find({entity_dim, entity_tag});
    const auto group =
      entity_it == entity_info.end()
        ? msh_group_from_entity(entity_dim, entity_tag, {}, physical_names)
        : msh_group_from_entity(
            entity_dim,
            entity_tag,
            entity_it->second.physical_tags,
            physical_names
          );

    for(std::size_t element_index = 0U; element_index < block_size; ++element_index) {
      std::size_t element_id = 0U;
      if(!(input >> element_id)) {
        fail(base::StatusCode::io_error, "MSH import found an incomplete element block.");
      }
      static_cast<void>(element_id);

      switch(element_type) {
      case 1: {
        std::array<long long, 2> node_ids {};
        if(!(input >> node_ids[0] >> node_ids[1])) {
          fail(base::StatusCode::io_error, "MSH import found a malformed line element.");
        }
        description.edges.push_back({
          map_msh_nodes(node_ids, node_map, "MSH line element"),
          group,
        });
      } break;
      case 2: {
        std::array<long long, 3> node_ids {};
        if(!(input >> node_ids[0] >> node_ids[1] >> node_ids[2])) {
          fail(base::StatusCode::io_error, "MSH import found a malformed triangle element.");
        }
        append_imported_face(
          description.faces,
          map_msh_nodes(node_ids, node_map, "MSH triangle element"),
          group
        );
      } break;
      case 4: {
        std::array<long long, 4> node_ids {};
        if(!(input >> node_ids[0] >> node_ids[1] >> node_ids[2] >> node_ids[3])) {
          fail(base::StatusCode::io_error, "MSH import found a malformed tetrahedron element.");
        }
        description.cells.push_back({
          map_msh_nodes(node_ids, node_map, "MSH tetrahedron element"),
          group,
        });
      } break;
      default:
        fail_unsupported_msh_element_type(element_type);
      }
    }
  }
}

void parse_msh41_elements_binary(
  std::istream &input,
  ImportedMeshDescription &description,
  const std::map<long long, std::size_t> &node_map,
  const std::map<std::pair<int, int>, MshEntityInfo> &entity_info,
  const std::map<std::pair<int, int>, std::string> &physical_names,
  bool swap_bytes,
  int tag_data_size
)
{
  std::array<std::uint64_t, 4> header = {0U, 0U, 0U, 0U};
  read_msh41_binary_size_values(
    input,
    header.data(),
    header.size(),
    swap_bytes,
    tag_data_size,
    "$Elements header"
  );

  const auto block_count = checked_msh_binary_size(
    header[0],
    "MSH import found an element-block count outside the supported range."
  );
  const auto total_count = checked_msh_binary_size(
    header[1],
    "MSH import found an element count outside the supported range."
  );
  std::size_t parsed_count = 0U;

  for(std::size_t block_index = 0U; block_index < block_count; ++block_index) {
    std::array<int, 3> block_header = {0, 0, 0};
    read_msh_binary_values(
      input,
      block_header.data(),
      block_header.size(),
      swap_bytes,
      "element block header"
    );

    const auto entity_dim = block_header[0];
    const auto entity_tag = block_header[1];
    const auto element_type = block_header[2];

    std::uint64_t block_size_value = 0U;
    read_msh41_binary_size_value(
      input,
      block_size_value,
      swap_bytes,
      tag_data_size,
      "element block size"
    );
    const auto block_size = checked_msh_binary_size(
      block_size_value,
      "MSH import found an element block size outside the supported range."
    );
    parsed_count += block_size;

    const auto entity_it = entity_info.find({entity_dim, entity_tag});
    const auto group =
      entity_it == entity_info.end()
        ? msh_group_from_entity(entity_dim, entity_tag, {}, physical_names)
        : msh_group_from_entity(
            entity_dim,
            entity_tag,
            entity_it->second.physical_tags,
            physical_names
          );

    for(std::size_t element_index = 0U; element_index < block_size; ++element_index) {
      std::uint64_t element_id = 0U;
      read_msh41_binary_size_value(
        input,
        element_id,
        swap_bytes,
        tag_data_size,
        "element tag"
      );
      static_cast<void>(element_id);

      switch(element_type) {
      case 1: {
        std::array<std::uint64_t, 2> node_ids = {0U, 0U};
        read_msh41_binary_size_values(
          input,
          node_ids.data(),
          node_ids.size(),
          swap_bytes,
          tag_data_size,
          "line element record"
        );
        description.edges.push_back({
          map_msh_nodes(
            std::array<long long, 2> {
              checked_msh_binary_tag(
                node_ids[0],
                "MSH import found a node tag outside the supported integer range."
              ),
              checked_msh_binary_tag(
                node_ids[1],
                "MSH import found a node tag outside the supported integer range."
              ),
            },
            node_map,
            "MSH line element"
          ),
          group,
        });
      } break;
      case 2: {
        std::array<std::uint64_t, 3> node_ids = {0U, 0U, 0U};
        read_msh41_binary_size_values(
          input,
          node_ids.data(),
          node_ids.size(),
          swap_bytes,
          tag_data_size,
          "triangle element record"
        );
        append_imported_face(
          description.faces,
          map_msh_nodes(
            std::array<long long, 3> {
              checked_msh_binary_tag(
                node_ids[0],
                "MSH import found a node tag outside the supported integer range."
              ),
              checked_msh_binary_tag(
                node_ids[1],
                "MSH import found a node tag outside the supported integer range."
              ),
              checked_msh_binary_tag(
                node_ids[2],
                "MSH import found a node tag outside the supported integer range."
              ),
            },
            node_map,
            "MSH triangle element"
          ),
          group
        );
      } break;
      case 4: {
        std::array<std::uint64_t, 4> node_ids = {0U, 0U, 0U, 0U};
        read_msh41_binary_size_values(
          input,
          node_ids.data(),
          node_ids.size(),
          swap_bytes,
          tag_data_size,
          "tetrahedron element record"
        );
        description.cells.push_back({
          map_msh_nodes(
            std::array<long long, 4> {
              checked_msh_binary_tag(
                node_ids[0],
                "MSH import found a node tag outside the supported integer range."
              ),
              checked_msh_binary_tag(
                node_ids[1],
                "MSH import found a node tag outside the supported integer range."
              ),
              checked_msh_binary_tag(
                node_ids[2],
                "MSH import found a node tag outside the supported integer range."
              ),
              checked_msh_binary_tag(
                node_ids[3],
                "MSH import found a node tag outside the supported integer range."
              ),
            },
            node_map,
            "MSH tetrahedron element"
          ),
          group,
        });
      } break;
      default:
        fail_unsupported_msh_element_type(element_type);
      }
    }
  }

  if(parsed_count != total_count) {
    fail(base::StatusCode::io_error, "MSH import found a malformed binary $Elements section.");
  }
}

[[nodiscard]] ImportedMeshDescription parse_msh_description(
  std::string_view path,
  const MshImportOptions &options
)
{
  validate_non_empty_path(path, "MSH import requires a non-empty path.");

  const std::string file_path(path);
  std::ifstream input(file_path, std::ios::binary);
  if(!input) {
    fail(base::StatusCode::io_error, "MSH import could not open the input file.");
  }

  ImportedMeshDescription description;
  description.name = default_domain_name(path);

  bool saw_mesh_format = false;
  bool saw_nodes = false;
  double mesh_version = 0.0;
  bool binary_format = false;
  bool swap_binary_bytes = false;
  int msh41_binary_data_size = 0;
  std::map<std::pair<int, int>, std::string> physical_names;
  std::map<std::pair<int, int>, MshEntityInfo> entity_info;
  std::map<long long, std::size_t> node_map;
  std::string line;

  while(std::getline(input, line)) {
    const auto token = trim_copy(line);
    if(token.empty()) {
      continue;
    }

    if(token == "$MeshFormat") {
      std::string format_line;
      if(!std::getline(input, format_line)) {
        fail(base::StatusCode::io_error, "MSH import found an incomplete $MeshFormat section.");
      }

      double version = 0.0;
      int format = 1;
      int data_size = 0;
      std::istringstream stream(format_line);
      if(!(stream >> version >> format >> data_size)) {
        fail(base::StatusCode::io_error, "MSH import could not parse the MeshFormat record.");
      }
      if(version != 2.2 && version != 4.1) {
        fail(
          base::StatusCode::unsupported,
          "MSH import currently supports only Gmsh 2.2 and Gmsh 4.1 files."
        );
      }
      if(format != 0 && format != 1) {
        fail(
          base::StatusCode::unsupported,
          "MSH import currently supports only ASCII (0) and binary (1) MeshFormat modes."
        );
      }
      if(format == 1) {
        if(version == 2.2 && data_size != static_cast<int>(sizeof(double))) {
          fail(
            base::StatusCode::unsupported,
            "MSH import currently supports binary Gmsh 2.2 only with MeshFormat data size matching sizeof(double)."
          );
        }
        if(version == 4.1 && !msh41_binary_data_size_supported(data_size)) {
          fail(
            base::StatusCode::unsupported,
            "MSH import currently supports binary Gmsh 4.1 only with 4-byte or 8-byte MeshFormat data sizes."
          );
        }
        msh41_binary_data_size = version == 4.1 ? data_size : 0;

        int one = 0;
        input.read(reinterpret_cast<char *>(&one), static_cast<std::streamsize>(sizeof(one)));
        if(!input) {
          fail(base::StatusCode::io_error, "MSH import found an incomplete binary MeshFormat marker.");
        }

        swap_binary_bytes = false;
        if(one != 1) {
          auto swapped_one = one;
          swap_msh_binary_value(swapped_one);
          if(swapped_one != 1) {
            fail(base::StatusCode::io_error, "MSH import found an invalid binary endianness marker.");
          }
          swap_binary_bytes = true;
        }
      }

      expect_msh_end(input, "$EndMeshFormat");
      saw_mesh_format = true;
      mesh_version = version;
      binary_format = format == 1;
      continue;
    }

    if(token == "$PhysicalNames") {
      std::string count_line;
      if(!std::getline(input, count_line)) {
        fail(base::StatusCode::io_error, "MSH import found an incomplete $PhysicalNames section.");
      }
      const auto count = parse_integer(trim_copy(count_line), "MSH PhysicalNames");
      for(long long index = 0; index < count; ++index) {
        std::string record_line;
        if(!std::getline(input, record_line)) {
          fail(base::StatusCode::io_error, "MSH import found an incomplete PhysicalNames record.");
        }

        std::istringstream stream(record_line);
        int dimension = 0;
        int tag = 0;
        if(!(stream >> dimension >> tag)) {
          fail(base::StatusCode::io_error, "MSH import found a malformed PhysicalNames record.");
        }

        if(options.read_physical_names) {
          const auto quote_begin = record_line.find('"');
          const auto quote_end = record_line.rfind('"');
          if(quote_begin == std::string::npos || quote_end == quote_begin) {
            fail(base::StatusCode::io_error, "MSH import found a malformed PhysicalNames string.");
          }
          physical_names[{dimension, tag}] =
            record_line.substr(quote_begin + 1U, quote_end - quote_begin - 1U);
        }
      }
      expect_msh_end(input, "$EndPhysicalNames");
      continue;
    }

    if(token == "$Entities") {
      if(!saw_mesh_format) {
        fail(base::StatusCode::io_error, "MSH import requires a $MeshFormat section before $Entities.");
      }
      if(mesh_version == 4.1) {
        if(binary_format) {
          parse_msh41_entities_binary(
            input,
            swap_binary_bytes,
            msh41_binary_data_size,
            entity_info
          );
        }
        else {
          parse_msh41_entities(input, entity_info);
        }
        expect_msh_end(input, "$EndEntities");
      }
      else {
        skip_msh_section(input, token);
      }
      continue;
    }

    if(token == "$PartitionedEntities") {
      if(!saw_mesh_format) {
        fail(
          base::StatusCode::io_error,
          "MSH import requires a $MeshFormat section before $PartitionedEntities."
        );
      }
      if(mesh_version != 4.1) {
        skip_msh_section(input, token);
        continue;
      }

      if(binary_format) {
        parse_msh41_partitioned_entities_binary(
          input,
          swap_binary_bytes,
          msh41_binary_data_size,
          entity_info
        );
      }
      else {
        parse_msh41_partitioned_entities(input, entity_info);
      }
      expect_msh_end(input, "$EndPartitionedEntities");
      continue;
    }

    if(token == "$Nodes") {
      if(!saw_mesh_format) {
        fail(base::StatusCode::io_error, "MSH import requires a $MeshFormat section before $Nodes.");
      }
      if(mesh_version == 4.1) {
        if(binary_format) {
          parse_msh41_nodes_binary(
            input,
            description,
            node_map,
            swap_binary_bytes,
            msh41_binary_data_size
          );
        }
        else {
          parse_msh41_nodes(input, description, node_map);
        }
      }
      else {
        if(binary_format) {
          parse_msh22_nodes_binary(input, description, node_map, swap_binary_bytes);
        }
        else {
          parse_msh22_nodes(input, description, node_map);
        }
      }
      expect_msh_end(input, "$EndNodes");
      saw_nodes = true;
      continue;
    }

    if(token == "$Elements") {
      if(!saw_mesh_format) {
        fail(base::StatusCode::io_error, "MSH import requires a $MeshFormat section before $Elements.");
      }
      if(mesh_version == 4.1) {
        if(binary_format) {
          parse_msh41_elements_binary(
            input,
            description,
            node_map,
            entity_info,
            physical_names,
            swap_binary_bytes,
            msh41_binary_data_size
          );
        }
        else {
          parse_msh41_elements(input, description, node_map, entity_info, physical_names);
        }
      }
      else {
        if(binary_format) {
          parse_msh22_elements_binary(
            input,
            description,
            node_map,
            physical_names,
            swap_binary_bytes
          );
        }
        else {
          parse_msh22_elements(input, description, node_map, physical_names);
        }
      }
      expect_msh_end(input, "$EndElements");
      continue;
    }

    if(token == "$Periodic") {
      if(!saw_mesh_format) {
        fail(base::StatusCode::io_error, "MSH import requires a $MeshFormat section before $Periodic.");
      }
      if(mesh_version == 4.1) {
        if(binary_format) {
          parse_msh41_periodic_binary(
            input,
            swap_binary_bytes,
            msh41_binary_data_size
          );
        }
        else {
          parse_msh41_periodic(input);
        }
        expect_msh_end(input, "$EndPeriodic");
      }
      else {
        skip_msh_section(input, token);
      }
      continue;
    }

    if(token == "$GhostElements") {
      if(!saw_mesh_format) {
        fail(base::StatusCode::io_error, "MSH import requires a $MeshFormat section before $GhostElements.");
      }
      if(mesh_version == 4.1) {
        if(binary_format) {
          parse_msh41_ghost_elements_binary(
            input,
            swap_binary_bytes,
            msh41_binary_data_size
          );
        }
        else {
          parse_msh41_ghost_elements(input);
        }
        expect_msh_end(input, "$EndGhostElements");
      }
      else {
        skip_msh_section(input, token);
      }
      continue;
    }

    if(token == "$Parametrizations") {
      if(!saw_mesh_format) {
        fail(base::StatusCode::io_error, "MSH import requires a $MeshFormat section before $Parametrizations.");
      }
      if(mesh_version == 4.1) {
        if(binary_format) {
          parse_msh41_parametrizations_binary(
            input,
            swap_binary_bytes,
            msh41_binary_data_size
          );
        }
        else {
          parse_msh41_parametrizations(input);
        }
        expect_msh_end(input, "$EndParametrizations");
      }
      else {
        skip_msh_section(input, token);
      }
      continue;
    }

    if(starts_with(token, "$End")) {
      continue;
    }

    if(starts_with(token, "$")) {
      if(binary_format) {
        fail(
          base::StatusCode::unsupported,
          "MSH import currently supports binary Gmsh files only with "
          "$PhysicalNames, $Entities, $PartitionedEntities, $Nodes, "
          "$Elements, $Periodic, $GhostElements, and $Parametrizations "
          "sections."
        );
      }
      skip_msh_section(input, token);
      continue;
    }

    fail(base::StatusCode::io_error, "MSH import found data outside of a section.");
  }

  if(!saw_mesh_format) {
    fail(base::StatusCode::io_error, "MSH import requires a $MeshFormat section.");
  }
  if(!saw_nodes) {
    fail(base::StatusCode::io_error, "MSH import requires a $Nodes section.");
  }

  return description;
}

[[nodiscard]] std::map<std::uint64_t, std::size_t> build_node_indices(
  const Domain &domain,
  EntityGroupRole role_filter
)
{
  std::map<std::uint64_t, std::size_t> node_indices;
  std::size_t next_index = 1U;

  for(const auto &entity_group : domain.entity_groups()) {
    if(entity_group.order() != EntityOrder::node) {
      continue;
    }
    if(entity_group.role() != role_filter) {
      continue;
    }

    for(const auto &node : entity_group.nodes()) {
      const auto key =
        (static_cast<std::uint64_t>(node.header.entity_group) << 32U) |
        static_cast<std::uint64_t>(node.header.index);
      node_indices.emplace(key, next_index++);
    }
  }

  return node_indices;
}

[[nodiscard]] std::uint64_t pack_entity_ref(EntityRef ref) noexcept
{
  return (static_cast<std::uint64_t>(ref.entity_group) << 32U) |
         static_cast<std::uint64_t>(ref.index);
}

[[nodiscard]] std::string sanitize_msh_name(std::string_view value)
{
  std::string sanitized(value);
  for(auto &character : sanitized) {
    if(character == '"') {
      character = '\'';
    }
  }
  return sanitized;
}

[[nodiscard]] std::string sanitize_obj_label(std::string_view value)
{
  std::string sanitized = trim_copy(value);
  if(sanitized.empty()) {
    return "default";
  }

  for(auto &character : sanitized) {
    if(character == '\n' || character == '\r') {
      character = ' ';
    }
  }
  return sanitized;
}

#if defined(sqmesh_HAVE_CGNS)

[[nodiscard]] EntityOrder highest_entity_order(const Domain &domain) noexcept
{
  if(domain.cell_count() > 0U) {
    return EntityOrder::cell;
  }
  if(domain.face_count() > 0U) {
    return EntityOrder::face;
  }
  if(domain.edge_count() > 0U) {
    return EntityOrder::edge;
  }
  return EntityOrder::node;
}

inline constexpr std::size_t cgns_name_limit = 32U;
inline constexpr char cgns_zone_id_descriptor_name[] = "SQMeshZoneId";

struct CgnsBoundaryRecord final {
  EntityOrder order = EntityOrder::node;
  std::string name {};
  std::int32_t bc_type_value = -1;
  cgsize_t start = 0;
  cgsize_t end = 0;
};

class CgnsFile final
{
public:
  CgnsFile() = default;
  explicit CgnsFile(int index) noexcept : index_(index) {}

  CgnsFile(const CgnsFile &) = delete;
  CgnsFile &operator=(const CgnsFile &) = delete;

  CgnsFile(CgnsFile &&other) noexcept : index_(other.index_)
  {
    other.index_ = -1;
  }

  CgnsFile &operator=(CgnsFile &&other) noexcept
  {
    if(this != &other) {
      close();
      index_ = other.index_;
      other.index_ = -1;
    }
    return *this;
  }

  ~CgnsFile()
  {
    close();
  }

  [[nodiscard]] int get() const noexcept
  {
    return index_;
  }

  void reset(int index = -1) noexcept
  {
    close();
    index_ = index;
  }

  void close() noexcept
  {
    if(index_ >= 0) {
      static_cast<void>(cg_close(index_));
      index_ = -1;
    }
  }

private:
  int index_ = -1;
};

[[noreturn]] void fail_cgns(
  base::StatusCode code,
  std::string_view operation
)
{
  fail(code, std::string(operation) + ": " + cg_get_error());
}

void check_cgns(
  int status,
  base::StatusCode code,
  std::string_view operation
)
{
  if(status != CG_OK) {
    fail_cgns(code, operation);
  }
}

[[nodiscard]] cgsize_t to_cgsize(
  std::size_t value,
  std::string_view context
)
{
  if(value > static_cast<std::size_t>(std::numeric_limits<cgsize_t>::max())) {
    fail(
      base::StatusCode::unsupported,
      std::string(context) + " exceeds the CGNS index range."
    );
  }
  return static_cast<cgsize_t>(value);
}

[[nodiscard]] std::size_t from_cgsize(
  cgsize_t value,
  std::string_view context
)
{
  if(value < 0) {
    fail(
      base::StatusCode::io_error,
      std::string(context) + " is negative."
    );
  }
  return static_cast<std::size_t>(value);
}

[[nodiscard]] std::string sanitize_cgns_label(
  std::string_view value,
  std::string_view fallback
)
{
  std::string sanitized = trim_copy(value);
  if(sanitized.empty()) {
    sanitized = std::string(fallback);
  }

  for(auto &character : sanitized) {
    if(!std::isalnum(static_cast<unsigned char>(character)) &&
       character != '_' && character != '-') {
      character = '_';
    }
  }

  if(sanitized.size() > cgns_name_limit) {
    sanitized.resize(cgns_name_limit);
  }
  if(sanitized.empty()) {
    sanitized = std::string(fallback);
  }
  return sanitized;
}

[[nodiscard]] std::string make_unique_cgns_label(
  std::string_view candidate,
  std::string_view fallback,
  std::map<std::string, std::size_t> &seen_names
)
{
  std::string base_name = sanitize_cgns_label(candidate, fallback);
  std::string unique_name = base_name;

  auto [it, inserted] = seen_names.emplace(unique_name, 0U);
  if(inserted) {
    return unique_name;
  }

  for(std::size_t suffix = ++it->second; ; ++suffix) {
    const auto suffix_text = "_" + std::to_string(suffix);
    const auto prefix_limit =
      cgns_name_limit > suffix_text.size() ? cgns_name_limit - suffix_text.size() : 0U;
    std::string candidate_name = base_name.substr(0U, prefix_limit) + suffix_text;
    if(seen_names.emplace(candidate_name, 0U).second) {
      return candidate_name;
    }
  }
}

struct ParsedCgnsQualifiedName final {
  bool qualified = false;
  std::size_t zone_index = 0U;
  std::string zone_name {};
  std::string local_name {};
};

struct CgnsBaseDescription final {
  int index = 0;
  std::string name {};
  int cell_dim = 0;
  int physical_dim = 0;
  int zone_count = 0;
};

struct CgnsZoneEntityGroupView final {
  const EntityGroup *entity_group = nullptr;
  std::string local_name {};
};

struct CgnsMultiZonePlan final {
  std::size_t zone_index = 0U;
  std::string zone_name {};
  std::vector<CgnsZoneEntityGroupView> entity_groups {};
};

struct CgnsInterfaceEntityGroupView final {
  const EntityGroup *entity_group = nullptr;
  std::string name {};
  std::size_t primary_zone_index = 0U;
  std::size_t secondary_zone_index = 0U;
};

struct CgnsMultiZoneExportPlan final {
  std::vector<CgnsMultiZonePlan> zones {};
  std::vector<CgnsInterfaceEntityGroupView> interfaces {};
};

struct CgnsWrittenZone final {
  int zone_index = 0;
  std::string name {};
  std::map<std::uint64_t, std::size_t> node_indices {};
};

[[nodiscard]] ParsedCgnsQualifiedName parse_cgns_qualified_entity_group_name(
  std::string_view name
)
{
  ParsedCgnsQualifiedName parsed;
  if(name.size() < 4U || name.front() != 'z') {
    return parsed;
  }

  std::size_t digit_end = 1U;
  while(digit_end < name.size() &&
        std::isdigit(static_cast<unsigned char>(name[digit_end])) != 0) {
    ++digit_end;
  }
  if(digit_end == 1U || digit_end >= name.size() || name[digit_end] != '_') {
    return parsed;
  }

  const auto separator = name.find("__", digit_end + 1U);
  if(separator == std::string_view::npos || separator <= digit_end + 1U) {
    return parsed;
  }

  try {
    parsed.zone_index = static_cast<std::size_t>(
      std::stoul(std::string(name.substr(1U, digit_end - 1U)))
    );
  }
  catch(const std::exception &) {
    return parsed;
  }

  parsed.zone_name = std::string(name.substr(digit_end + 1U, separator - digit_end - 1U));
  parsed.local_name = std::string(name.substr(separator + 2U));
  parsed.qualified = !parsed.zone_name.empty() && !parsed.local_name.empty();
  return parsed;
}

[[nodiscard]] std::size_t cgns_section_zone_id(
  int file_index,
  int base_index,
  int zone_index,
  int section_index
)
{
  check_cgns(
    cg_goto(
      file_index,
      base_index,
      "Zone_t",
      zone_index,
      "Elements_t",
      section_index,
      "end"
    ),
    base::StatusCode::io_error,
    "CGNS import could not access section metadata"
  );

  int descriptor_count = 0;
  check_cgns(
    cg_ndescriptors(&descriptor_count),
    base::StatusCode::io_error,
    "CGNS import could not query section descriptors"
  );

  for(int descriptor_index = 1; descriptor_index <= descriptor_count; ++descriptor_index) {
    char descriptor_name[cgns_name_limit + 1] = {};
    char *descriptor_text = nullptr;
    check_cgns(
      cg_descriptor_read(descriptor_index, descriptor_name, &descriptor_text),
      base::StatusCode::io_error,
      "CGNS import could not read section metadata descriptor"
    );

    std::string text = descriptor_text == nullptr ? std::string() : std::string(descriptor_text);
    if(descriptor_text != nullptr) {
      static_cast<void>(cg_free(descriptor_text));
    }

    if(std::string_view(descriptor_name) != cgns_zone_id_descriptor_name) {
      continue;
    }

    try {
      const auto zone_id = std::stoul(text);
      return static_cast<std::size_t>(zone_id);
    }
    catch(const std::exception &) {
      fail(
        base::StatusCode::io_error,
        "CGNS import found a malformed SQMeshZoneId descriptor."
      );
    }
  }

  return static_cast<std::size_t>(section_index);
}

template <std::size_t N>
[[nodiscard]] std::array<std::size_t, N> map_cgns_nodes(
  const cgsize_t *connectivity,
  std::size_t node_count,
  std::string_view record_name
)
{
  std::array<std::size_t, N> nodes {};
  for(std::size_t offset = 0; offset < N; ++offset) {
    const auto raw = from_cgsize(connectivity[offset], record_name);
    if(raw == 0U || raw > node_count) {
      fail(
        base::StatusCode::io_error,
        std::string(record_name) + " references an out-of-range node id."
      );
    }
    nodes[offset] = raw - 1U;
  }
  return nodes;
}

template <std::size_t N>
[[nodiscard]] std::array<std::size_t, N> map_cgns_nodes_with_offset(
  const cgsize_t *connectivity,
  std::size_t node_count,
  std::size_t node_offset,
  std::string_view record_name
)
{
  auto nodes = map_cgns_nodes<N>(connectivity, node_count, record_name);
  for(auto &node : nodes) {
    node += node_offset;
  }
  return nodes;
}

[[nodiscard]] std::string cgns_zone_import_prefix(
  std::string_view zone_name,
  int zone_index,
  std::string_view base_prefix = {}
)
{
  const auto label = [&]() {
    auto trimmed = trim_copy(zone_name);
    if(trimmed.empty()) {
      trimmed = "zone";
    }
    return trimmed;
  }();
  const auto zone_prefix = "z" + std::to_string(zone_index) + "_" + label;
  if(base_prefix.empty()) {
    return zone_prefix;
  }
  return std::string(base_prefix) + "__" + zone_prefix;
}

[[nodiscard]] std::string cgns_base_label(std::string_view base_name)
{
  auto label = trim_copy(base_name);
  if(label.empty()) {
    label = "base";
  }
  return label;
}

[[nodiscard]] std::string cgns_base_import_prefix(
  std::string_view base_name,
  int base_index
)
{
  return "b" + std::to_string(base_index) + "_" + cgns_base_label(base_name);
}

[[nodiscard]] std::string cgns_zone_label(std::string_view zone_name)
{
  auto label = trim_copy(zone_name);
  if(label.empty()) {
    label = "zone";
  }
  return label;
}

[[nodiscard]] EntityGroupImportInfo cgns_entity_group_import_info(
  int base_index,
  std::string_view base_name,
  int zone_index,
  std::string_view zone_name,
  std::string_view local_name,
  std::int32_t bc_type_value = -1
)
{
  EntityGroupImportInfo import_info;
  import_info.format = EntityGroupImportFormat::cgns;
  import_info.cgns.base_index = static_cast<std::uint32_t>(base_index);
  import_info.cgns.base_name = cgns_base_label(base_name);
  import_info.cgns.zone_index = static_cast<std::uint32_t>(zone_index);
  import_info.cgns.zone_name = cgns_zone_label(zone_name);
  import_info.cgns.local_name = std::string(local_name);
  import_info.cgns.bc_type_value = bc_type_value;
  return import_info;
}

[[nodiscard]] EntityGroupImportInfo cgns_entity_group_import_info_for_name(
  int fallback_base_index,
  std::string_view fallback_base_name,
  int fallback_zone_index,
  std::string_view fallback_zone_name,
  std::string_view entity_group_name,
  std::string_view fallback_local_name,
  std::int32_t bc_type_value = -1
)
{
  const auto parsed = parse_cgns_qualified_entity_group_name(entity_group_name);
  if(parsed.qualified &&
     parsed.zone_index <= static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
    EntityGroupImportInfo import_info;
    import_info.format = EntityGroupImportFormat::cgns;
    import_info.cgns.base_index = static_cast<std::uint32_t>(fallback_base_index);
    import_info.cgns.base_name = cgns_base_label(fallback_base_name);
    import_info.cgns.zone_index = static_cast<std::uint32_t>(parsed.zone_index);
    import_info.cgns.zone_name = std::move(parsed.zone_name);
    import_info.cgns.local_name = std::move(parsed.local_name);
    import_info.cgns.bc_type_value = bc_type_value;
    return import_info;
  }

  return cgns_entity_group_import_info(
    fallback_base_index,
    fallback_base_name,
    fallback_zone_index,
    fallback_zone_name,
    fallback_local_name,
    bc_type_value
  );
}

[[nodiscard]] std::string qualify_cgns_import_name(
  std::string_view zone_prefix,
  std::string_view local_name
)
{
  if(zone_prefix.empty()) {
    return std::string(local_name);
  }
  if(local_name.empty()) {
    return std::string(zone_prefix);
  }
  return std::string(zone_prefix) + "__" + std::string(local_name);
}

[[nodiscard]] std::string cgns_section_fallback_label(
  CGNS_ENUMT(ElementType_t) element_type
)
{
  switch(element_type) {
  case CGNS_ENUMV(BAR_2):
    return "edge_section";
  case CGNS_ENUMV(TRI_3):
  case CGNS_ENUMV(QUAD_4):
    return "surface_section";
  case CGNS_ENUMV(TETRA_4):
    return "volume_section";
  case CGNS_ENUMV(MIXED):
    return "mixed_section";
  default:
    return "section";
  }
}

void append_cgns_zone_to_description(
  int file_index,
  int base_index,
  std::string_view base_name,
  int zone_index,
  int cell_dim,
  int physical_dim,
  bool qualify_base_names,
  bool qualify_group_names,
  const CgnsImportOptions &options,
  ImportedMeshDescription &description
)
{
  CGNS_ENUMT(ZoneType_t) zone_type = CGNS_ENUMV(ZoneTypeNull);
  check_cgns(
    cg_zone_type(file_index, base_index, zone_index, &zone_type),
    base::StatusCode::io_error,
    "CGNS import could not read the zone type"
  );
  if(zone_type != CGNS_ENUMV(Unstructured)) {
    fail(
      base::StatusCode::unsupported,
      "CGNS import currently supports only unstructured zones."
    );
  }

  char zone_name[cgns_name_limit + 1] = {};
  cgsize_t zone_size[3] = {0, 0, 0};
  check_cgns(
    cg_zone_read(file_index, base_index, zone_index, zone_name, zone_size),
    base::StatusCode::io_error,
    "CGNS import could not read the zone"
  );

  const auto zone_label = cgns_zone_label(zone_name);
  const auto base_label = cgns_base_label(base_name);
  const auto base_prefix =
    qualify_base_names ? cgns_base_import_prefix(base_name, base_index) : std::string {};
  const auto zone_prefix =
    qualify_group_names
      ? cgns_zone_import_prefix(zone_name, zone_index, base_prefix)
      : std::string {};
  const auto node_offset = description.nodes.size();
  const auto node_count = from_cgsize(zone_size[0], "CGNS zone node count");
  if(node_count == 0U) {
    fail(base::StatusCode::io_error, "CGNS import found no nodes.");
  }
  description.nodes.resize(node_offset + node_count);

  int coordinate_count = 0;
  check_cgns(
    cg_ncoords(file_index, base_index, zone_index, &coordinate_count),
    base::StatusCode::io_error,
    "CGNS import could not query coordinate arrays"
  );
  if(coordinate_count < 2) {
    fail(
      base::StatusCode::unsupported,
      "CGNS import requires CoordinateX and CoordinateY arrays."
    );
  }

  std::map<std::string, int> coordinate_names;
  for(int coordinate_index = 1; coordinate_index <= coordinate_count; ++coordinate_index) {
    CGNS_ENUMT(DataType_t) coordinate_type = CGNS_ENUMV(DataTypeNull);
    char coordinate_name[cgns_name_limit + 1] = {};
    check_cgns(
      cg_coord_info(
        file_index,
        base_index,
        zone_index,
        coordinate_index,
        &coordinate_type,
        coordinate_name
      ),
      base::StatusCode::io_error,
      "CGNS import could not inspect coordinate metadata"
    );
    coordinate_names.emplace(std::string(coordinate_name), coordinate_index);
  }

  if(coordinate_names.find("CoordinateX") == coordinate_names.end() ||
     coordinate_names.find("CoordinateY") == coordinate_names.end()) {
    fail(
      base::StatusCode::unsupported,
      "CGNS import currently requires canonical CoordinateX and CoordinateY arrays."
    );
  }

  std::vector<double> x(node_count, 0.0);
  std::vector<double> y(node_count, 0.0);
  std::vector<double> z(node_count, 0.0);
  const cgsize_t range_min[1] = {1};
  const cgsize_t range_max[1] = {to_cgsize(node_count, "CGNS coordinate count")};
  check_cgns(
    cg_coord_read(
      file_index,
      base_index,
      zone_index,
      "CoordinateX",
      CGNS_ENUMV(RealDouble),
      range_min,
      range_max,
      x.data()
    ),
    base::StatusCode::io_error,
    "CGNS import could not read CoordinateX"
  );
  check_cgns(
    cg_coord_read(
      file_index,
      base_index,
      zone_index,
      "CoordinateY",
      CGNS_ENUMV(RealDouble),
      range_min,
      range_max,
      y.data()
    ),
    base::StatusCode::io_error,
    "CGNS import could not read CoordinateY"
  );
  if(coordinate_names.find("CoordinateZ") != coordinate_names.end()) {
    check_cgns(
      cg_coord_read(
        file_index,
        base_index,
        zone_index,
        "CoordinateZ",
        CGNS_ENUMV(RealDouble),
        range_min,
        range_max,
        z.data()
      ),
      base::StatusCode::io_error,
      "CGNS import could not read CoordinateZ"
    );
  }
  else if(physical_dim == 3) {
    fail(
      base::StatusCode::unsupported,
      "CGNS import currently requires a canonical CoordinateZ array for 3D coordinates."
    );
  }

  for(std::size_t index = 0; index < node_count; ++index) {
    description.nodes[node_offset + index] = {x[index], y[index], z[index]};
  }

  int section_count = 0;
  check_cgns(
    cg_nsections(file_index, base_index, zone_index, &section_count),
    base::StatusCode::io_error,
    "CGNS import could not query element sections"
  );
  if(section_count == 0) {
    fail(base::StatusCode::io_error, "CGNS import found no element sections.");
  }

  std::map<cgsize_t, std::vector<std::size_t>> edge_element_map;
  std::map<cgsize_t, std::vector<std::size_t>> face_element_map;
  std::map<cgsize_t, std::vector<std::size_t>> cell_element_map;

  const auto append_element =
    [&](CGNS_ENUMT(ElementType_t) element_type,
        const cgsize_t *element_connectivity,
        cgsize_t element_id,
        const GroupInfo &group) {
      switch(element_type) {
      case CGNS_ENUMV(BAR_2): {
        const auto edge_index = description.edges.size();
        description.edges.push_back({
          map_cgns_nodes_with_offset<2>(
            element_connectivity,
            node_count,
            node_offset,
            "CGNS BAR_2 element"
          ),
          group,
        });
        edge_element_map[element_id].push_back(edge_index);
        break;
      }
      case CGNS_ENUMV(TRI_3): {
        const auto face_index = description.faces.size();
        append_imported_face(
          description.faces,
          map_cgns_nodes_with_offset<3>(
            element_connectivity,
            node_count,
            node_offset,
            "CGNS TRI_3 element"
          ),
          group
        );
        face_element_map[element_id].push_back(face_index);
        break;
      }
      case CGNS_ENUMV(QUAD_4): {
        const auto nodes = map_cgns_nodes_with_offset<4>(
          element_connectivity,
          node_count,
          node_offset,
          "CGNS QUAD_4 element"
        );
        const auto first_face_index = description.faces.size();
        append_imported_face(description.faces, {nodes[0], nodes[1], nodes[2]}, group);
        append_imported_face(description.faces, {nodes[0], nodes[2], nodes[3]}, group);
        face_element_map[element_id].push_back(first_face_index);
        face_element_map[element_id].push_back(first_face_index + 1U);
        break;
      }
      case CGNS_ENUMV(TETRA_4): {
        const auto cell_index = description.cells.size();
        description.cells.push_back({
          map_cgns_nodes_with_offset<4>(
            element_connectivity,
            node_count,
            node_offset,
            "CGNS TETRA_4 element"
          ),
          group,
        });
        cell_element_map[element_id].push_back(cell_index);
        break;
      }
      default:
        fail(
          base::StatusCode::unsupported,
          "CGNS import currently supports only BAR_2, TRI_3, QUAD_4, TETRA_4, "
          "and MIXED sections containing those families. Unsupported element "
          "type: " + std::string(cg_ElementTypeName(element_type)) + '.'
        );
      }
    };

  for(int section_index = 1; section_index <= section_count; ++section_index) {
    char section_name[cgns_name_limit + 1] = {};
    CGNS_ENUMT(ElementType_t) element_type = CGNS_ENUMV(ElementTypeNull);
    cgsize_t start = 0;
    cgsize_t end = 0;
    int boundary_element_count = 0;
    int parent_flag = 0;
    check_cgns(
      cg_section_read(
        file_index,
        base_index,
        zone_index,
        section_index,
        section_name,
        &element_type,
        &start,
        &end,
        &boundary_element_count,
        &parent_flag
      ),
      base::StatusCode::io_error,
      "CGNS import could not read an element section"
    );
    static_cast<void>(boundary_element_count);

    if(parent_flag != 0) {
      fail(
        base::StatusCode::unsupported,
        "CGNS import does not yet support sections with parent data."
      );
    }

    const auto section_size = from_cgsize(end - start + 1, "CGNS section size");
    if(section_size == 0U) {
      continue;
    }

    GroupInfo group;
    group.zone_id = static_cast<std::uint32_t>(
      cgns_section_zone_id(file_index, base_index, zone_index, section_index)
    );
    const auto trimmed_section_name = trim_copy(section_name);
    const auto imported_local_name =
      trimmed_section_name.empty()
        ? cgns_section_fallback_label(element_type)
        : trimmed_section_name;
    if(qualify_group_names) {
      group.name = qualify_cgns_import_name(
        zone_prefix,
        options.read_section_names && !trimmed_section_name.empty()
          ? trimmed_section_name
          : cgns_section_fallback_label(element_type)
      );
    }
    else if(options.read_section_names && !trimmed_section_name.empty()) {
      group.name = trimmed_section_name;
    }
    group.import_info = cgns_entity_group_import_info_for_name(
      base_index,
      base_label,
      zone_index,
      zone_label,
      group.name,
      imported_local_name
    );

    if(element_type == CGNS_ENUMV(MIXED)) {
      cgsize_t element_data_size = 0;
      check_cgns(
        cg_ElementDataSize(file_index, base_index, zone_index, section_index, &element_data_size),
        base::StatusCode::io_error,
        "CGNS import could not query MIXED section connectivity size"
      );

      std::vector<cgsize_t> connectivity(
        from_cgsize(element_data_size, "CGNS MIXED section connectivity size"),
        0
      );
      std::vector<cgsize_t> connectivity_offsets(section_size + 1U, 0);
      check_cgns(
        cg_poly_elements_read(
          file_index,
          base_index,
          zone_index,
          section_index,
          connectivity.data(),
          connectivity_offsets.data(),
          nullptr
        ),
        base::StatusCode::io_error,
        "CGNS import could not read MIXED section connectivity"
      );

      for(std::size_t element_index = 0; element_index < section_size; ++element_index) {
        const auto begin_offset = from_cgsize(
          connectivity_offsets[element_index],
          "CGNS MIXED section connectivity offset"
        );
        const auto end_offset = from_cgsize(
          connectivity_offsets[element_index + 1U],
          "CGNS MIXED section connectivity offset"
        );
        if(begin_offset >= end_offset || end_offset > connectivity.size()) {
          fail(
            base::StatusCode::io_error,
            "CGNS import found a malformed MIXED section connectivity table."
          );
        }

        const auto mixed_type = static_cast<CGNS_ENUMT(ElementType_t)>(
          connectivity[begin_offset]
        );
        int nodes_per_element = 0;
        check_cgns(
          cg_npe(mixed_type, &nodes_per_element),
          base::StatusCode::io_error,
          "CGNS import could not query MIXED element arity"
        );
        const auto record_size =
          static_cast<std::size_t>(nodes_per_element) + 1U;
        if(end_offset - begin_offset != record_size) {
          fail(
            base::StatusCode::io_error,
            "CGNS import found a malformed MIXED element record."
          );
        }

        append_element(
          mixed_type,
          connectivity.data() + begin_offset + 1U,
          start + static_cast<cgsize_t>(element_index),
          group
        );
      }
      continue;
    }

    int nodes_per_element = 0;
    check_cgns(
      cg_npe(element_type, &nodes_per_element),
      base::StatusCode::io_error,
      "CGNS import could not query section arity"
    );
    std::vector<cgsize_t> connectivity(
      section_size * static_cast<std::size_t>(nodes_per_element),
      0
    );
    check_cgns(
      cg_elements_read(
        file_index,
        base_index,
        zone_index,
        section_index,
        connectivity.data(),
        nullptr
      ),
      base::StatusCode::io_error,
      "CGNS import could not read section connectivity"
    );

    for(std::size_t element_index = 0; element_index < section_size; ++element_index) {
      append_element(
        element_type,
        connectivity.data() +
          element_index * static_cast<std::size_t>(nodes_per_element),
        start + static_cast<cgsize_t>(element_index),
        group
      );
    }
  }

  int boundary_condition_count = 0;
  check_cgns(
    cg_nbocos(file_index, base_index, zone_index, &boundary_condition_count),
    base::StatusCode::io_error,
    "CGNS import could not query boundary-condition nodes"
  );

  if(boundary_condition_count <= 0) {
    return;
  }

  std::vector<std::string> edge_boundary_names(description.edges.size());
  std::vector<std::string> face_boundary_names(description.faces.size());
  std::vector<std::string> edge_boundary_local_names(description.edges.size());
  std::vector<std::string> face_boundary_local_names(description.faces.size());
  std::vector<std::int32_t> edge_boundary_type_values(description.edges.size(), -1);
  std::vector<std::int32_t> face_boundary_type_values(description.faces.size(), -1);

  const auto assign_boundary_name =
    [&](std::vector<std::string> &assigned_names,
        std::vector<std::string> &assigned_local_names,
        std::vector<std::int32_t> &assigned_bc_types,
        const std::vector<std::size_t> &indices,
        std::string_view boundary_name,
        std::string_view boundary_local_name,
        std::int32_t boundary_type_value) {
      for(const auto index : indices) {
        if(index >= assigned_names.size()) {
          fail(
            base::StatusCode::io_error,
            "CGNS import found a boundary-condition element outside the "
            "supported imported subset."
          );
        }
        if(!assigned_names[index].empty() &&
           assigned_names[index] != boundary_name) {
          fail(
            base::StatusCode::unsupported,
            "CGNS import does not support multiple boundary-condition names "
            "attached to the same supported element."
          );
        }
        if(assigned_bc_types[index] != -1 &&
           assigned_bc_types[index] != boundary_type_value) {
          fail(
            base::StatusCode::unsupported,
            "CGNS import does not support multiple BCType_t values attached to the same supported element."
          );
        }
        assigned_names[index] = std::string(boundary_name);
        assigned_local_names[index] = std::string(boundary_local_name);
        assigned_bc_types[index] = boundary_type_value;
      }
    };

  const auto read_boundary_point_set =
    [&](int boundary_index,
        CGNS_ENUMT(PointSetType_t) point_set_type,
        cgsize_t point_count) {
      std::size_t storage_size = 0U;
      switch(point_set_type) {
      case CGNS_ENUMV(PointList):
      case CGNS_ENUMV(ElementList):
        storage_size = from_cgsize(point_count, "CGNS boundary-condition point count");
        if(storage_size == 0U) {
          fail(
            base::StatusCode::io_error,
            "CGNS import found an empty boundary-condition point set."
          );
        }
        break;
      case CGNS_ENUMV(PointRange):
      case CGNS_ENUMV(ElementRange):
        storage_size = 2U;
        break;
      default:
        fail(
          base::StatusCode::unsupported,
          "CGNS import currently supports only PointList, PointRange, "
          "ElementList, and ElementRange boundary-condition point sets."
        );
      }

      std::vector<cgsize_t> raw_points(storage_size, 0);
      check_cgns(
        cg_boco_read(
          file_index,
          base_index,
          zone_index,
          boundary_index,
          raw_points.data(),
          nullptr
        ),
        base::StatusCode::io_error,
        "CGNS import could not read a boundary-condition point set"
      );

      std::vector<cgsize_t> expanded_points;
      switch(point_set_type) {
      case CGNS_ENUMV(PointList):
      case CGNS_ENUMV(ElementList):
        expanded_points = raw_points;
        break;
      case CGNS_ENUMV(PointRange):
      case CGNS_ENUMV(ElementRange): {
        if(raw_points.size() != 2U ||
           raw_points[0] <= 0 ||
           raw_points[1] <= 0) {
          fail(
            base::StatusCode::io_error,
            "CGNS import found a malformed boundary-condition point range."
          );
        }

        const auto first = raw_points[0];
        const auto last = raw_points[1];
        const auto range_begin = first < last ? first : last;
        const auto range_end = first < last ? last : first;
        expanded_points.reserve(
          from_cgsize(
            range_end - range_begin + 1,
            "CGNS boundary-condition expanded point count"
          )
        );
        for(cgsize_t point = range_begin; point <= range_end; ++point) {
          expanded_points.push_back(point);
        }
        break;
      }
      default:
        break;
      }

      return expanded_points;
    };

  for(int boundary_index = 1; boundary_index <= boundary_condition_count; ++boundary_index) {
    char boundary_name[cgns_name_limit + 1] = {};
    CGNS_ENUMT(BCType_t) boundary_type = CGNS_ENUMV(BCTypeNull);
    CGNS_ENUMT(PointSetType_t) point_set_type = CGNS_ENUMV(PointSetTypeNull);
    cgsize_t point_count = 0;
    int normal_index[3] = {0, 0, 0};
    cgsize_t normal_list_size = 0;
    CGNS_ENUMT(DataType_t) normal_data_type = CGNS_ENUMV(DataTypeNull);
    int dataset_count = 0;
    check_cgns(
      cg_boco_info(
        file_index,
        base_index,
        zone_index,
        boundary_index,
        boundary_name,
        &boundary_type,
        &point_set_type,
        &point_count,
        normal_index,
        &normal_list_size,
        &normal_data_type,
        &dataset_count
      ),
      base::StatusCode::io_error,
      "CGNS import could not read boundary-condition metadata"
    );
    static_cast<void>(normal_index);
    static_cast<void>(normal_list_size);
    static_cast<void>(normal_data_type);
    static_cast<void>(dataset_count);

    auto boundary_local_label = trim_copy(boundary_name);
    if(boundary_local_label.empty()) {
      boundary_local_label = "boundary_condition";
    }
    const auto boundary_type_value = static_cast<std::int32_t>(boundary_type);
    auto boundary_label = boundary_local_label;
    if(qualify_group_names) {
      boundary_label = qualify_cgns_import_name(zone_prefix, boundary_label);
    }
    const auto referenced_points =
      read_boundary_point_set(boundary_index, point_set_type, point_count);

    if(point_set_type == CGNS_ENUMV(PointList) ||
       point_set_type == CGNS_ENUMV(PointRange)) {
      CGNS_ENUMT(GridLocation_t) location = CGNS_ENUMV(Vertex);
      const auto location_status =
        cg_boco_gridlocation_read(file_index, base_index, zone_index, boundary_index, &location);
      if(location_status != CG_OK &&
         location_status != CG_NODE_NOT_FOUND) {
        fail_cgns(
          base::StatusCode::io_error,
          "CGNS import could not read a boundary-condition grid location"
        );
      }

      if(location == CGNS_ENUMV(EdgeCenter)) {
        for(const auto point : referenced_points) {
          const auto edge_it = edge_element_map.find(point);
          if(edge_it == edge_element_map.end()) {
            fail(
              base::StatusCode::io_error,
              "CGNS import found an EdgeCenter boundary-condition reference "
              "to an unknown BAR_2 element id."
            );
          }
          assign_boundary_name(
            edge_boundary_names,
            edge_boundary_local_names,
            edge_boundary_type_values,
            edge_it->second,
            boundary_label,
            boundary_local_label,
            boundary_type_value
          );
        }
        continue;
      }

      if(location == CGNS_ENUMV(FaceCenter) ||
         (cell_dim == 2 && location == CGNS_ENUMV(CellCenter))) {
        for(const auto point : referenced_points) {
          const auto face_it = face_element_map.find(point);
          if(face_it == face_element_map.end()) {
            fail(
              base::StatusCode::io_error,
              "CGNS import found a FaceCenter boundary-condition reference "
              "to an unknown supported face element id."
            );
          }
          assign_boundary_name(
            face_boundary_names,
            face_boundary_local_names,
            face_boundary_type_values,
            face_it->second,
            boundary_label,
            boundary_local_label,
            boundary_type_value
          );
        }
        continue;
      }

      fail(
        base::StatusCode::unsupported,
        "CGNS import currently represents PointList/PointRange boundary "
        "conditions only for EdgeCenter plus FaceCenter/CellCenter on the "
        "supported unstructured subset."
      );
    }

    EntityOrder inferred_order = EntityOrder::node;
    bool order_resolved = false;
    for(const auto point : referenced_points) {
      const auto edge_it = edge_element_map.find(point);
      if(edge_it != edge_element_map.end()) {
        if(order_resolved && inferred_order != EntityOrder::edge) {
          fail(
            base::StatusCode::unsupported,
            "CGNS import does not support ElementList/ElementRange boundary "
            "conditions that mix edge and face references."
          );
        }
        inferred_order = EntityOrder::edge;
        order_resolved = true;
        assign_boundary_name(
          edge_boundary_names,
          edge_boundary_local_names,
          edge_boundary_type_values,
          edge_it->second,
          boundary_label,
          boundary_local_label,
          boundary_type_value
        );
        continue;
      }

      const auto face_it = face_element_map.find(point);
      if(face_it != face_element_map.end()) {
        if(order_resolved && inferred_order != EntityOrder::face) {
          fail(
            base::StatusCode::unsupported,
            "CGNS import does not support ElementList/ElementRange boundary "
            "conditions that mix edge and face references."
          );
        }
        inferred_order = EntityOrder::face;
        order_resolved = true;
        assign_boundary_name(
          face_boundary_names,
          face_boundary_local_names,
          face_boundary_type_values,
          face_it->second,
          boundary_label,
          boundary_local_label,
          boundary_type_value
        );
        continue;
      }

      if(cell_element_map.find(point) != cell_element_map.end()) {
        fail(
          base::StatusCode::unsupported,
          "CGNS import currently represents boundary conditions only on "
          "supported boundary BAR_2 and face-element records."
        );
      }

      fail(
        base::StatusCode::io_error,
        "CGNS import found an ElementList/ElementRange boundary-condition "
        "reference to an unknown supported element id."
      );
    }
  }

  for(std::size_t edge_index = 0; edge_index < edge_boundary_names.size(); ++edge_index) {
    if(!edge_boundary_names[edge_index].empty()) {
      description.edges[edge_index].group.name = edge_boundary_names[edge_index];
      description.edges[edge_index].group.import_info =
        cgns_entity_group_import_info_for_name(
          base_index,
          base_label,
          zone_index,
          zone_label,
          description.edges[edge_index].group.name,
          edge_boundary_local_names[edge_index],
          edge_boundary_type_values[edge_index]
        );
    }
  }
  for(std::size_t face_index = 0; face_index < face_boundary_names.size(); ++face_index) {
    if(!face_boundary_names[face_index].empty()) {
      description.faces[face_index].group.name = face_boundary_names[face_index];
      description.faces[face_index].group.import_info =
        cgns_entity_group_import_info_for_name(
          base_index,
          base_label,
          zone_index,
          zone_label,
          description.faces[face_index].group.name,
          face_boundary_local_names[face_index],
          face_boundary_type_values[face_index]
      );
    }
  }
}

[[nodiscard]] CgnsBaseDescription read_cgns_base_description(
  int file_index,
  int base_index
)
{
  char base_name[cgns_name_limit + 1] = {};
  int cell_dim = 0;
  int physical_dim = 0;
  check_cgns(
    cg_base_read(file_index, base_index, base_name, &cell_dim, &physical_dim),
    base::StatusCode::io_error,
    "CGNS import could not read the base node"
  );
  if(cell_dim < 2 || cell_dim > 3) {
    fail(
      base::StatusCode::unsupported,
      "CGNS import currently supports only 2D or 3D unstructured bases."
    );
  }
  if(physical_dim < 2 || physical_dim > 3) {
    fail(
      base::StatusCode::unsupported,
      "CGNS import currently supports only 2D or 3D physical coordinates."
    );
  }

  int zone_count = 0;
  check_cgns(
    cg_nzones(file_index, base_index, &zone_count),
    base::StatusCode::io_error,
    "CGNS import could not query zone count"
  );
  if(zone_count <= 0) {
    fail(base::StatusCode::io_error, "CGNS import found no zones.");
  }

  CgnsBaseDescription description;
  description.index = base_index;
  description.name = std::string(base_name);
  description.cell_dim = cell_dim;
  description.physical_dim = physical_dim;
  description.zone_count = zone_count;
  return description;
}

[[nodiscard]] ImportedMeshDescription parse_cgns_description(
  std::string_view path,
  const CgnsImportOptions &options
)
{
  validate_non_empty_path(path, "CGNS import requires a non-empty path.");

  int file_index = -1;
  check_cgns(
    cg_open(std::string(path).c_str(), CG_MODE_READ, &file_index),
    base::StatusCode::io_error,
    "CGNS import could not open the input file"
  );
  CgnsFile file(file_index);

  int base_count = 0;
  check_cgns(
    cg_nbases(file.get(), &base_count),
    base::StatusCode::io_error,
    "CGNS import could not query base count"
  );
  if(base_count <= 0) {
    fail(base::StatusCode::io_error, "CGNS import found no bases.");
  }

  std::vector<CgnsBaseDescription> bases;
  bases.reserve(static_cast<std::size_t>(base_count));
  std::size_t total_zone_count = 0U;
  for(int base_index = 1; base_index <= base_count; ++base_index) {
    auto base = read_cgns_base_description(file.get(), base_index);
    if(!bases.empty() &&
       (base.cell_dim != bases.front().cell_dim ||
        base.physical_dim != bases.front().physical_dim)) {
      fail(
        base::StatusCode::unsupported,
        "CGNS import currently requires all supported bases in one file to "
        "share one cell dimension and one physical coordinate dimension."
      );
    }
    total_zone_count += static_cast<std::size_t>(base.zone_count);
    bases.push_back(std::move(base));
  }

  ImportedMeshDescription description;
  if(base_count == 1 && total_zone_count == 1U) {
    char zone_name[cgns_name_limit + 1] = {};
    cgsize_t zone_size[3] = {0, 0, 0};
    check_cgns(
      cg_zone_read(file.get(), 1, 1, zone_name, zone_size),
      base::StatusCode::io_error,
      "CGNS import could not read the zone"
    );
    description.name =
      zone_name[0] == '\0' ? bases.front().name : std::string(zone_name);
  }
  else if(base_count == 1) {
    description.name = bases.front().name;
  }
  else {
    description.name = default_domain_name(path);
  }
  if(description.name.empty()) {
    description.name = default_domain_name(path);
  }

  const bool qualify_base_names = base_count > 1;
  const bool qualify_group_names = qualify_base_names || total_zone_count > 1U;
  for(const auto &base : bases) {
    for(int zone_index = 1; zone_index <= base.zone_count; ++zone_index) {
      append_cgns_zone_to_description(
        file.get(),
        base.index,
        base.name,
        zone_index,
        base.cell_dim,
        base.physical_dim,
        qualify_base_names,
        qualify_group_names,
        options,
        description
      );
    }
  }

  if(description.faces.empty() && description.cells.empty()) {
    fail(
      base::StatusCode::unsupported,
      "CGNS import requires at least one supported TRI_3, QUAD_4, or "
      "TETRA_4 element."
    );
  }

  return description;
}

void validate_domain_for_cgns_export(const Domain &domain)
{
  validate_domain_for_msh_export(domain);

  const auto highest_order = highest_entity_order(domain);
  if(highest_order == EntityOrder::node || highest_order == EntityOrder::edge) {
    fail(
      base::StatusCode::unsupported,
      "CGNS export currently supports only surface or volume meshes."
    );
  }
}

[[nodiscard]] int cgns_cell_dimension(const Domain &domain) noexcept
{
  return highest_entity_order(domain) == EntityOrder::cell ? 3 : 2;
}

[[nodiscard]] std::size_t cgns_topology_element_count(const Domain &domain) noexcept
{
  return cgns_cell_dimension(domain) == 3 ? domain.cell_count() : domain.face_count();
}

[[nodiscard]] std::vector<double> collect_coordinate_component(
  const Domain &domain,
  const std::map<std::uint64_t, std::size_t> &node_indices,
  std::size_t component
)
{
  std::vector<double> values(node_indices.size(), 0.0);
  for(const auto &entity_group : domain.entity_groups()) {
    if(entity_group.order() != EntityOrder::node) {
      continue;
    }
    for(const auto &node : entity_group.nodes()) {
      const auto node_it =
        node_indices.find(pack_entity_ref({node.header.entity_group, node.header.index}));
      if(node_it == node_indices.end()) {
        continue;
      }
      const auto index = node_it->second - 1U;
      values[index] = node.coordinates[component];
    }
  }
  return values;
}

[[nodiscard]] std::vector<cgsize_t> collect_entity_group_connectivity(
  const Domain &domain,
  const EntityGroup &entity_group,
  const std::map<std::uint64_t, std::size_t> &node_indices
)
{
  std::vector<cgsize_t> connectivity;

  switch(entity_group.order()) {
  case EntityOrder::edge:
    connectivity.reserve(entity_group.edges().size() * 2U);
    for(const auto &edge : entity_group.edges()) {
      for(const auto node_ref : domain.edge_nodes({edge.header.entity_group, edge.header.index})) {
        connectivity.push_back(
          to_cgsize(node_indices.at(pack_entity_ref(node_ref)), "CGNS edge connectivity")
        );
      }
    }
    break;
  case EntityOrder::face:
    connectivity.reserve(entity_group.faces().size() * 3U);
    for(const auto &face : entity_group.faces()) {
      for(const auto node_ref : domain.face_nodes({face.header.entity_group, face.header.index})) {
        connectivity.push_back(
          to_cgsize(node_indices.at(pack_entity_ref(node_ref)), "CGNS face connectivity")
        );
      }
    }
    break;
  case EntityOrder::cell:
    connectivity.reserve(entity_group.cells().size() * 4U);
    for(const auto &cell : entity_group.cells()) {
      for(const auto node_ref : domain.cell_nodes({cell.header.entity_group, cell.header.index})) {
        connectivity.push_back(
          to_cgsize(node_indices.at(pack_entity_ref(node_ref)), "CGNS cell connectivity")
        );
      }
    }
    break;
  default:
    break;
  }

  return connectivity;
}

[[nodiscard]] int cgns_cell_dimension(
  const CgnsMultiZonePlan &plan
) noexcept
{
  for(const auto &entity_group_view : plan.entity_groups) {
    if(entity_group_view.entity_group->order() == EntityOrder::cell &&
       entity_group_view.entity_group->entity_count() > 0U) {
      return 3;
    }
  }
  return 2;
}

[[nodiscard]] std::size_t cgns_topology_element_count(
  const CgnsMultiZonePlan &plan
) noexcept
{
  std::size_t count = 0U;
  const auto cell_dimension = cgns_cell_dimension(plan);
  for(const auto &entity_group_view : plan.entity_groups) {
    if(cell_dimension == 3) {
      if(entity_group_view.entity_group->order() == EntityOrder::cell) {
        count += entity_group_view.entity_group->entity_count();
      }
    }
    else if(entity_group_view.entity_group->order() == EntityOrder::face) {
      count += entity_group_view.entity_group->entity_count();
    }
  }
  return count;
}

void append_zone_entity_group_node_indices(
  const Domain &domain,
  const EntityGroup &entity_group,
  std::map<std::uint64_t, std::size_t> &node_indices,
  std::size_t &next_index
)
{
  const auto append_node_ref =
    [&](EntityRef node_ref) {
      if(node_indices.emplace(pack_entity_ref(node_ref), next_index).second) {
        ++next_index;
      }
    };

  switch(entity_group.order()) {
  case EntityOrder::edge:
    for(const auto &edge : entity_group.edges()) {
      for(const auto node_ref : domain.edge_nodes({edge.header.entity_group, edge.header.index})) {
        append_node_ref(node_ref);
      }
    }
    break;
  case EntityOrder::face:
    for(const auto &face : entity_group.faces()) {
      for(const auto node_ref : domain.face_nodes({face.header.entity_group, face.header.index})) {
        append_node_ref(node_ref);
      }
    }
    break;
  case EntityOrder::cell:
    for(const auto &cell : entity_group.cells()) {
      for(const auto node_ref : domain.cell_nodes({cell.header.entity_group, cell.header.index})) {
        append_node_ref(node_ref);
      }
    }
    break;
  default:
    break;
  }
}

[[nodiscard]] std::map<std::uint64_t, std::size_t> build_zone_node_indices(
  const Domain &domain,
  const CgnsMultiZonePlan &plan
)
{
  std::map<std::uint64_t, std::size_t> node_indices;
  std::size_t next_index = 1U;
  for(const auto &entity_group_view : plan.entity_groups) {
    append_zone_entity_group_node_indices(
      domain,
      *entity_group_view.entity_group,
      node_indices,
      next_index
    );
  }
  return node_indices;
}

[[nodiscard]] ParsedCgnsQualifiedName cgns_entity_group_export_identity(const EntityGroup &entity_group)
{
  const auto &import_info = entity_group.import_info();
  if(import_info.format == EntityGroupImportFormat::cgns &&
     import_info.cgns.zone_index != invalid_index &&
     !import_info.cgns.zone_name.empty() &&
     !import_info.cgns.local_name.empty()) {
    ParsedCgnsQualifiedName parsed;
    parsed.zone_index = import_info.cgns.zone_index;
    parsed.zone_name = import_info.cgns.zone_name;
    parsed.local_name = import_info.cgns.local_name;
    parsed.qualified = true;
    return parsed;
  }

  return parse_cgns_qualified_entity_group_name(entity_group.name());
}

[[nodiscard]] ParsedCgnsQualifiedName require_cgns_entity_group_export_identity(
  const EntityGroup &entity_group,
  std::string_view context
)
{
  const auto parsed = cgns_entity_group_export_identity(entity_group);
  if(parsed.qualified) {
    return parsed;
  }

  fail(
    base::StatusCode::unsupported,
    "CGNS multi-zone export requires " + std::string(context) +
      " to carry CGNS import metadata or an importer-style zone-qualified name."
  );
}

[[nodiscard]] std::int32_t cgns_entity_group_boundary_type_value(const EntityGroup &entity_group) noexcept
{
  const auto &import_info = entity_group.import_info();
  if(import_info.format != EntityGroupImportFormat::cgns) {
    return -1;
  }
  return import_info.cgns.bc_type_value;
}

[[nodiscard]] std::string cgns_interface_entity_group_name(const EntityGroup &entity_group)
{
  const auto &import_info = entity_group.import_info();
  if(import_info.format == EntityGroupImportFormat::cgns &&
     !import_info.cgns.local_name.empty()) {
    return import_info.cgns.local_name;
  }

  const auto name = trim_copy(entity_group.name());
  if(!name.empty()) {
    return name;
  }

  return entity_group.order() == EntityOrder::edge
           ? std::string("edge_interface")
           : std::string("face_interface");
}

[[nodiscard]] std::pair<std::size_t, std::size_t> cgns_interface_entity_group_zone_pair(
  const Domain &domain,
  const EntityGroup &entity_group
)
{
  bool saw_pair = false;
  std::pair<std::size_t, std::size_t> zone_pair {0U, 0U};

  const auto merge_zone_pair =
    [&](std::size_t first_zone_index, std::size_t second_zone_index) {
      if(first_zone_index == second_zone_index) {
        fail(
          base::StatusCode::unsupported,
          "CGNS multi-zone export found an interface entity_group whose adjacent "
          "zone identities collapse to one zone."
        );
      }

      if(second_zone_index < first_zone_index) {
        std::swap(first_zone_index, second_zone_index);
      }

      if(!saw_pair) {
        zone_pair = {first_zone_index, second_zone_index};
        saw_pair = true;
        return;
      }

      if(zone_pair.first != first_zone_index ||
         zone_pair.second != second_zone_index) {
        fail(
          base::StatusCode::unsupported,
          "CGNS multi-zone export currently requires each interface entity_group "
          "to map to exactly one adjacent zone pair."
        );
      }
    };

  switch(entity_group.order()) {
  case EntityOrder::edge:
    for(const auto &edge : entity_group.edges()) {
      const auto edge_ref = EntityRef {edge.header.entity_group, edge.header.index};
      const auto left_face = domain.adjacent_face(edge_ref, FaceSide::left);
      const auto right_face = domain.adjacent_face(edge_ref, FaceSide::right);
      if(!is_valid(left_face) || !is_valid(right_face)) {
        fail(
          base::StatusCode::unsupported,
          "CGNS multi-zone export currently requires each interface edge to "
          "keep left/right face adjacency."
        );
      }

      const auto left_zone =
        require_cgns_entity_group_export_identity(domain.entity_group(left_face.entity_group), "an interface-adjacent face entity_group");
      const auto right_zone =
        require_cgns_entity_group_export_identity(domain.entity_group(right_face.entity_group), "an interface-adjacent face entity_group");
      merge_zone_pair(left_zone.zone_index, right_zone.zone_index);
    }
    break;
  case EntityOrder::face:
    for(const auto &face : entity_group.faces()) {
      const auto face_ref = EntityRef {face.header.entity_group, face.header.index};
      const auto left_cell = domain.adjacent_cell(face_ref, FaceSide::left);
      const auto right_cell = domain.adjacent_cell(face_ref, FaceSide::right);
      if(!is_valid(left_cell) || !is_valid(right_cell)) {
        fail(
          base::StatusCode::unsupported,
          "CGNS multi-zone export currently requires each interface face to "
          "keep left/right cell adjacency."
        );
      }

      const auto left_zone =
        require_cgns_entity_group_export_identity(domain.entity_group(left_cell.entity_group), "an interface-adjacent region entity_group");
      const auto right_zone =
        require_cgns_entity_group_export_identity(domain.entity_group(right_cell.entity_group), "an interface-adjacent region entity_group");
      merge_zone_pair(left_zone.zone_index, right_zone.zone_index);
    }
    break;
  default:
    fail(
      base::StatusCode::unsupported,
      "CGNS multi-zone export currently reconstructs interface connectivity "
      "only on edge or face entity_groups."
    );
  }

  if(!saw_pair) {
    fail(
      base::StatusCode::unsupported,
      "CGNS multi-zone export found an empty interface entity_group."
    );
  }

  const auto primary_zone_index = entity_group.primary_region_zone_id();
  const auto secondary_zone_index = entity_group.secondary_region_zone_id();
  if(primary_zone_index != invalid_index &&
     secondary_zone_index != invalid_index) {
    std::size_t first_zone_index = primary_zone_index;
    std::size_t second_zone_index = secondary_zone_index;
    if(second_zone_index < first_zone_index) {
      std::swap(first_zone_index, second_zone_index);
    }
    if(first_zone_index != zone_pair.first ||
       second_zone_index != zone_pair.second) {
      fail(
        base::StatusCode::unsupported,
        "CGNS multi-zone export found an interface entity_group whose stored "
        "adjacent zone ids disagree with its actual adjacency."
      );
    }
  }

  return zone_pair;
}

[[nodiscard]] std::vector<EntityRef> collect_cgns_interface_entity_group_nodes(
  const Domain &domain,
  const EntityGroup &entity_group
)
{
  std::vector<EntityRef> node_refs;
  std::set<std::uint64_t> seen_nodes;

  const auto append_node_ref =
    [&](EntityRef node_ref) {
      const auto packed_ref = pack_entity_ref(node_ref);
      if(seen_nodes.insert(packed_ref).second) {
        node_refs.push_back(node_ref);
      }
    };

  switch(entity_group.order()) {
  case EntityOrder::edge:
    for(const auto &edge : entity_group.edges()) {
      for(const auto node_ref :
          domain.edge_nodes({edge.header.entity_group, edge.header.index})) {
        append_node_ref(node_ref);
      }
    }
    break;
  case EntityOrder::face:
    for(const auto &face : entity_group.faces()) {
      for(const auto node_ref :
          domain.face_nodes({face.header.entity_group, face.header.index})) {
        append_node_ref(node_ref);
      }
    }
    break;
  default:
    break;
  }

  return node_refs;
}

[[nodiscard]] std::vector<cgsize_t> map_cgns_interface_nodes_to_zone_points(
  const std::map<std::uint64_t, std::size_t> &node_indices,
  const std::vector<EntityRef> &node_refs,
  std::string_view interface_name,
  std::string_view zone_name
)
{
  std::vector<cgsize_t> points;
  points.reserve(node_refs.size());

  for(const auto node_ref : node_refs) {
    const auto node_it = node_indices.find(pack_entity_ref(node_ref));
    if(node_it == node_indices.end()) {
      fail(
        base::StatusCode::unsupported,
        "CGNS multi-zone export currently reconstructs interface "
        "connectivity only for shared-node interfaces; entity_group '" +
          std::string(interface_name) + "' references nodes that are not part "
          "of zone '" + std::string(zone_name) + "'."
      );
    }
    points.push_back(
      to_cgsize(node_it->second, "CGNS interface point list")
    );
  }

  return points;
}

[[nodiscard]] CgnsMultiZoneExportPlan build_cgns_multi_zone_export_plan(
  const Domain &domain
)
{
  std::map<std::size_t, CgnsMultiZonePlan> zones;
  std::vector<const EntityGroup *> interface_entity_groups;
  std::set<std::pair<std::uint32_t, std::string>> imported_base_ids;
  bool saw_qualified_entity_group = false;
  bool saw_unqualified_non_node_entity_group = false;

  for(const auto &entity_group : domain.entity_groups()) {
    if(entity_group.order() == EntityOrder::node || entity_group.entity_count() == 0U) {
      continue;
    }

    const auto &import_info = entity_group.import_info();
    if(import_info.format == EntityGroupImportFormat::cgns &&
       import_info.cgns.base_index != invalid_index &&
       !import_info.cgns.base_name.empty()) {
      imported_base_ids.emplace(
        import_info.cgns.base_index,
        import_info.cgns.base_name
      );
    }

    if(entity_group.semantic() == EntityGroupSemantic::interface) {
      interface_entity_groups.push_back(&entity_group);
      continue;
    }

    const auto parsed = cgns_entity_group_export_identity(entity_group);
    if(!parsed.qualified) {
      saw_unqualified_non_node_entity_group = true;
      continue;
    }

    saw_qualified_entity_group = true;
    auto [zone_it, inserted] = zones.emplace(parsed.zone_index, CgnsMultiZonePlan {});
    auto &zone = zone_it->second;
    if(inserted) {
      zone.zone_index = parsed.zone_index;
      zone.zone_name = parsed.zone_name;
    }
    else if(zone.zone_name != parsed.zone_name) {
      fail(
        base::StatusCode::unsupported,
        "CGNS multi-zone export found inconsistent zone-qualified entity_group names."
      );
    }
    zone.entity_groups.push_back({&entity_group, parsed.local_name});
  }

  if(imported_base_ids.size() > 1U) {
    fail(
      base::StatusCode::unsupported,
      "CGNS multi-zone export does not yet preserve multiple imported bases."
    );
  }

  if(!saw_qualified_entity_group) {
    if(!interface_entity_groups.empty()) {
      fail(
        base::StatusCode::unsupported,
        "CGNS multi-zone export requires zone-qualified non-interface entity_groups "
        "before it can reconstruct interface connectivity."
      );
    }
    return {};
  }
  if(saw_unqualified_non_node_entity_group) {
    fail(
      base::StatusCode::unsupported,
      "CGNS multi-zone export currently requires CGNS import metadata or "
      "importer-style zone-qualified names on every non-node entity_group."
    );
  }
  if(zones.size() <= 1U) {
    if(!interface_entity_groups.empty()) {
      fail(
        base::StatusCode::unsupported,
        "CGNS multi-zone export requires at least two reconstructed zones "
        "when interface entity_groups are present."
      );
    }
    return {};
  }

  CgnsMultiZoneExportPlan export_plan;
  export_plan.zones.reserve(zones.size());
  std::set<std::size_t> zone_indices;
  int expected_cell_dimension = 0;
  for(auto &[zone_index, zone] : zones) {
    static_cast<void>(zone_index);
    if(cgns_topology_element_count(zone) == 0U) {
      fail(
        base::StatusCode::unsupported,
        "CGNS multi-zone export requires each zone to contain supported face or cell elements."
      );
    }

    const auto cell_dimension = cgns_cell_dimension(zone);
    if(expected_cell_dimension == 0) {
      expected_cell_dimension = cell_dimension;
    }
    else if(cell_dimension != expected_cell_dimension) {
      fail(
        base::StatusCode::unsupported,
        "CGNS multi-zone export currently requires all zones to share one mesh dimension."
      );
    }

    zone_indices.insert(zone.zone_index);
    export_plan.zones.push_back(std::move(zone));
  }

  for(const auto *entity_group : interface_entity_groups) {
    if(expected_cell_dimension == 3 &&
       entity_group->order() != EntityOrder::face) {
      fail(
        base::StatusCode::unsupported,
        "CGNS multi-zone export currently reconstructs interface "
        "connectivity only for face entity_groups on volume zones."
      );
    }
    if(expected_cell_dimension == 2 &&
       entity_group->order() != EntityOrder::edge) {
      fail(
        base::StatusCode::unsupported,
        "CGNS multi-zone export currently reconstructs interface "
        "connectivity only for edge entity_groups on surface zones."
      );
    }

    const auto zone_pair = cgns_interface_entity_group_zone_pair(domain, *entity_group);
    if(zone_indices.find(zone_pair.first) == zone_indices.end() ||
       zone_indices.find(zone_pair.second) == zone_indices.end()) {
      fail(
        base::StatusCode::unsupported,
        "CGNS multi-zone export found an interface entity_group that references a "
        "zone outside the current export set."
      );
    }

    export_plan.interfaces.push_back({
      entity_group,
      cgns_interface_entity_group_name(*entity_group),
      zone_pair.first,
      zone_pair.second,
    });
  }

  return export_plan;
}

void write_cgns_section_zone_id_descriptor(
  int file_index,
  int base_index,
  int zone_index,
  int section_index,
  std::uint32_t zone_id
)
{
  check_cgns(
    cg_goto(
      file_index,
      base_index,
      "Zone_t",
      zone_index,
      "Elements_t",
      section_index,
      "end"
    ),
    base::StatusCode::io_error,
    "CGNS export could not access section metadata"
  );
  const auto zone_id_text = std::to_string(zone_id);
  check_cgns(
    cg_descriptor_write(
      cgns_zone_id_descriptor_name,
      zone_id_text.c_str()
    ),
    base::StatusCode::io_error,
    "CGNS export could not write the SQMeshZoneId descriptor"
  );
}

void write_cgns_boundary_condition(
  int file_index,
  int base_index,
  int zone_index,
  const CgnsBoundaryRecord &record,
  CGNS_ENUMT(GridLocation_t) location,
  std::map<std::string, std::size_t> &used_boundary_names
)
{
  const auto boundary_name = make_unique_cgns_label(
    record.name,
    record.order == EntityOrder::edge ? "edge_boundary" : "face_boundary",
    used_boundary_names
  );
  const cgsize_t point_range[2] = {
    record.start,
    record.end,
  };

  int boundary_index = 0;
  check_cgns(
    cg_boco_write(
      file_index,
      base_index,
      zone_index,
      boundary_name.c_str(),
      static_cast<CGNS_ENUMT(BCType_t)>(
        record.bc_type_value >= 0
          ? record.bc_type_value
          : static_cast<std::int32_t>(CGNS_ENUMV(BCGeneral))
      ),
      CGNS_ENUMV(PointRange),
      2,
      point_range,
      &boundary_index
    ),
    base::StatusCode::io_error,
    "CGNS export could not write a boundary-condition node"
  );
  check_cgns(
    cg_boco_gridlocation_write(
      file_index,
      base_index,
      zone_index,
      boundary_index,
      location
    ),
    base::StatusCode::io_error,
    "CGNS export could not write a boundary-condition grid location"
  );
}

void write_cgns_interface_connectivity(
  int file_index,
  int base_index,
  int zone_index,
  std::string_view interface_name,
  std::string_view fallback_name,
  const std::vector<cgsize_t> &points,
  std::string_view donor_zone_name,
  const std::vector<cgsize_t> &donor_points,
  std::map<std::string, std::size_t> &used_connectivity_names
)
{
  if(points.empty() || donor_points.empty() ||
     points.size() != donor_points.size()) {
    fail(
      base::StatusCode::unsupported,
      "CGNS multi-zone export found an invalid interface point list."
    );
  }

  const auto connectivity_name = make_unique_cgns_label(
    interface_name,
    fallback_name,
    used_connectivity_names
  );

  int connectivity_index = 0;
  check_cgns(
    cg_conn_write(
      file_index,
      base_index,
      zone_index,
      connectivity_name.c_str(),
      CGNS_ENUMV(Vertex),
      CGNS_ENUMV(Abutting1to1),
      CGNS_ENUMV(PointList),
      to_cgsize(points.size(), "CGNS interface point count"),
      points.data(),
      std::string(donor_zone_name).c_str(),
      CGNS_ENUMV(Unstructured),
      CGNS_ENUMV(PointListDonor),
      CGNS_ENUMV(DataTypeNull),
      to_cgsize(donor_points.size(), "CGNS donor interface point count"),
      donor_points.data(),
      &connectivity_index
    ),
    base::StatusCode::io_error,
    "CGNS export could not write a zone interface connectivity node"
  );
}

void write_cgns_zone_contents(
  int file_index,
  int base_index,
  int zone_index,
  int cell_dimension,
  const Domain &domain,
  const std::map<std::uint64_t, std::size_t> &node_indices,
  const std::vector<CgnsZoneEntityGroupView> &entity_group_views
)
{
  std::map<std::string, std::size_t> used_section_names;
  std::vector<CgnsBoundaryRecord> boundary_records;
  cgsize_t next_element = 1;

  const auto write_sections =
    [&](EntityOrder order, CGNS_ENUMT(ElementType_t) element_type) {
      for(const auto &entity_group_view : entity_group_views) {
        const auto &entity_group = *entity_group_view.entity_group;
        if(entity_group.order() != order || entity_group.entity_count() == 0U) {
          continue;
        }

        const auto connectivity = collect_entity_group_connectivity(domain, entity_group, node_indices);
        const auto section_size = to_cgsize(entity_group.entity_count(), "CGNS section size");
        const auto section_name = make_unique_cgns_label(
          entity_group_view.local_name.empty() ? std::string_view(entity_group.name()) : std::string_view(entity_group_view.local_name),
          order == EntityOrder::edge
            ? "bar_section"
            : (order == EntityOrder::face ? "tri_section" : "tet_section"),
          used_section_names
        );

        int section_index = 0;
        check_cgns(
          cg_section_write(
            file_index,
            base_index,
            zone_index,
            section_name.c_str(),
            element_type,
            next_element,
            next_element + section_size - 1,
            0,
            connectivity.data(),
            &section_index
          ),
          base::StatusCode::io_error,
          "CGNS export could not write an element section"
        );
        write_cgns_section_zone_id_descriptor(
          file_index,
          base_index,
          zone_index,
          section_index,
          entity_group.zone_id()
        );
        if(order == EntityOrder::edge ||
           (order == EntityOrder::face && entity_group.is_boundary())) {
          boundary_records.push_back({
            order,
            entity_group_view.local_name.empty() ? std::string(entity_group.name()) : entity_group_view.local_name,
            cgns_entity_group_boundary_type_value(entity_group),
            next_element,
            next_element + section_size - 1,
          });
        }
        next_element += section_size;
      }
    };

  write_sections(EntityOrder::edge, CGNS_ENUMV(BAR_2));
  write_sections(EntityOrder::face, CGNS_ENUMV(TRI_3));
  write_sections(EntityOrder::cell, CGNS_ENUMV(TETRA_4));

  std::map<std::string, std::size_t> used_boundary_names;
  for(const auto &record : boundary_records) {
    write_cgns_boundary_condition(
      file_index,
      base_index,
      zone_index,
      record,
      record.order == EntityOrder::edge
        ? CGNS_ENUMV(EdgeCenter)
        : (cell_dimension == 2
             ? CGNS_ENUMV(CellCenter)
             : CGNS_ENUMV(FaceCenter)),
      used_boundary_names
    );
  }
}

void write_cgns_unstructured(
  const Domain &domain,
  std::string_view path,
  const CgnsExportOptions &options
);

void write_cgns_multi_zone_unstructured(
  const Domain &domain,
  std::string_view path,
  const CgnsExportOptions &options
)
{
  validate_domain_for_cgns_export(domain);

  const auto export_plan = build_cgns_multi_zone_export_plan(domain);
  if(export_plan.zones.empty()) {
    write_cgns_unstructured(domain, path, options);
    return;
  }

  int file_index = -1;
  check_cgns(
    cg_open(std::string(path).c_str(), CG_MODE_WRITE, &file_index),
    base::StatusCode::io_error,
    "CGNS export could not open the output file"
  );
  CgnsFile file(file_index);

  const std::string base_name = sanitize_cgns_label(
    options.base_name.empty() ? default_domain_name(path) : options.base_name,
    "sqmesh_base"
  );
  const int cell_dimension = cgns_cell_dimension(export_plan.zones.front());

  int base_index = 0;
  check_cgns(
    cg_base_write(
      file.get(),
      base_name.c_str(),
      cell_dimension,
      3,
      &base_index
    ),
    base::StatusCode::io_error,
    "CGNS export could not write the base node"
  );

  std::map<std::string, std::size_t> used_zone_names;
  std::map<std::size_t, CgnsWrittenZone> written_zones;
  for(const auto &plan : export_plan.zones) {
    const auto node_indices = build_zone_node_indices(domain, plan);
    const cgsize_t zone_size[3] = {
      to_cgsize(node_indices.size(), "CGNS zone node count"),
      to_cgsize(cgns_topology_element_count(plan), "CGNS topological element count"),
      0,
    };
    const auto zone_name = make_unique_cgns_label(
      plan.zone_name,
      "sqmesh_zone",
      used_zone_names
    );

    int zone_index = 0;
    check_cgns(
      cg_zone_write(
        file.get(),
        base_index,
        zone_name.c_str(),
        zone_size,
        CGNS_ENUMV(Unstructured),
        &zone_index
      ),
      base::StatusCode::io_error,
      "CGNS export could not write the zone"
    );

    const auto x = collect_coordinate_component(domain, node_indices, 0U);
    const auto y = collect_coordinate_component(domain, node_indices, 1U);
    const auto z = collect_coordinate_component(domain, node_indices, 2U);

    int coordinate_index = 0;
    check_cgns(
      cg_coord_write(
        file.get(),
        base_index,
        zone_index,
        CGNS_ENUMV(RealDouble),
        "CoordinateX",
        x.data(),
        &coordinate_index
      ),
      base::StatusCode::io_error,
      "CGNS export could not write CoordinateX"
    );
    check_cgns(
      cg_coord_write(
        file.get(),
        base_index,
        zone_index,
        CGNS_ENUMV(RealDouble),
        "CoordinateY",
        y.data(),
        &coordinate_index
      ),
      base::StatusCode::io_error,
      "CGNS export could not write CoordinateY"
    );
    check_cgns(
      cg_coord_write(
        file.get(),
        base_index,
        zone_index,
        CGNS_ENUMV(RealDouble),
        "CoordinateZ",
        z.data(),
        &coordinate_index
      ),
      base::StatusCode::io_error,
      "CGNS export could not write CoordinateZ"
    );

    write_cgns_zone_contents(
      file.get(),
      base_index,
      zone_index,
      cell_dimension,
      domain,
      node_indices,
      plan.entity_groups
    );

    written_zones.emplace(
      plan.zone_index,
      CgnsWrittenZone {
        zone_index,
        zone_name,
        node_indices,
      }
    );
  }

  std::map<std::size_t, std::map<std::string, std::size_t>> used_connectivity_names;
  for(const auto &interface_view : export_plan.interfaces) {
    const auto primary_zone_it = written_zones.find(interface_view.primary_zone_index);
    const auto secondary_zone_it = written_zones.find(interface_view.secondary_zone_index);
    if(primary_zone_it == written_zones.end() ||
       secondary_zone_it == written_zones.end()) {
      fail(
        base::StatusCode::unsupported,
        "CGNS multi-zone export found an interface entity_group whose adjacent "
        "zones were not written."
      );
    }

    const auto node_refs =
      collect_cgns_interface_entity_group_nodes(domain, *interface_view.entity_group);
    const auto primary_points = map_cgns_interface_nodes_to_zone_points(
      primary_zone_it->second.node_indices,
      node_refs,
      interface_view.name,
      primary_zone_it->second.name
    );
    const auto secondary_points = map_cgns_interface_nodes_to_zone_points(
      secondary_zone_it->second.node_indices,
      node_refs,
      interface_view.name,
      secondary_zone_it->second.name
    );

    write_cgns_interface_connectivity(
      file.get(),
      base_index,
      primary_zone_it->second.zone_index,
      interface_view.name,
      interface_view.entity_group->order() == EntityOrder::edge
        ? "edge_interface"
        : "face_interface",
      primary_points,
      secondary_zone_it->second.name,
      secondary_points,
      used_connectivity_names[interface_view.primary_zone_index]
    );
    write_cgns_interface_connectivity(
      file.get(),
      base_index,
      secondary_zone_it->second.zone_index,
      interface_view.name,
      interface_view.entity_group->order() == EntityOrder::edge
        ? "edge_interface"
        : "face_interface",
      secondary_points,
      primary_zone_it->second.name,
      primary_points,
      used_connectivity_names[interface_view.secondary_zone_index]
    );
  }
}

void write_cgns_unstructured(
  const Domain &domain,
  std::string_view path,
  const CgnsExportOptions &options
)
{
  validate_domain_for_cgns_export(domain);

  int file_index = -1;
  check_cgns(
    cg_open(std::string(path).c_str(), CG_MODE_WRITE, &file_index),
    base::StatusCode::io_error,
    "CGNS export could not open the output file"
  );
  CgnsFile file(file_index);

  const auto node_indices = build_node_indices(domain);
  const std::string base_name = sanitize_cgns_label(
    options.base_name.empty() ? default_domain_name(path) : options.base_name,
    "sqmesh_base"
  );
  const std::string zone_name = sanitize_cgns_label(
    options.zone_name.empty()
      ? (domain.name().empty() ? default_domain_name(path) : std::string(domain.name()))
      : options.zone_name,
    "sqmesh_zone"
  );

  int base_index = 0;
  check_cgns(
    cg_base_write(
      file.get(),
      base_name.c_str(),
      cgns_cell_dimension(domain),
      3,
      &base_index
    ),
    base::StatusCode::io_error,
    "CGNS export could not write the base node"
  );

  const cgsize_t zone_size[3] = {
    to_cgsize(node_indices.size(), "CGNS node count"),
    to_cgsize(cgns_topology_element_count(domain), "CGNS topological element count"),
    0,
  };
  int zone_index = 0;
  check_cgns(
    cg_zone_write(
      file.get(),
      base_index,
      zone_name.c_str(),
      zone_size,
      CGNS_ENUMV(Unstructured),
      &zone_index
    ),
    base::StatusCode::io_error,
    "CGNS export could not write the zone"
  );

  const auto x = collect_coordinate_component(domain, node_indices, 0U);
  const auto y = collect_coordinate_component(domain, node_indices, 1U);
  const auto z = collect_coordinate_component(domain, node_indices, 2U);

  int coordinate_index = 0;
  check_cgns(
    cg_coord_write(
      file.get(),
      base_index,
      zone_index,
      CGNS_ENUMV(RealDouble),
      "CoordinateX",
      x.data(),
      &coordinate_index
    ),
    base::StatusCode::io_error,
    "CGNS export could not write CoordinateX"
  );
  check_cgns(
    cg_coord_write(
      file.get(),
      base_index,
      zone_index,
      CGNS_ENUMV(RealDouble),
      "CoordinateY",
      y.data(),
      &coordinate_index
    ),
    base::StatusCode::io_error,
    "CGNS export could not write CoordinateY"
  );
  check_cgns(
    cg_coord_write(
      file.get(),
      base_index,
      zone_index,
      CGNS_ENUMV(RealDouble),
      "CoordinateZ",
      z.data(),
      &coordinate_index
    ),
    base::StatusCode::io_error,
    "CGNS export could not write CoordinateZ"
  );

  std::map<std::string, std::size_t> used_section_names;
  std::vector<CgnsBoundaryRecord> boundary_records;
  cgsize_t next_element = 1;
  auto write_sections = [&](EntityOrder order, CGNS_ENUMT(ElementType_t) element_type) {
    for(const auto &entity_group : domain.entity_groups()) {
      if(entity_group.order() != order || entity_group.entity_count() == 0U) {
        continue;
      }

      const auto connectivity = collect_entity_group_connectivity(domain, entity_group, node_indices);
      const auto section_size = to_cgsize(entity_group.entity_count(), "CGNS section size");
      const auto section_name = make_unique_cgns_label(
        entity_group.name(),
        order == EntityOrder::edge
          ? "bar_section"
          : (order == EntityOrder::face ? "tri_section" : "tet_section"),
        used_section_names
      );

      int section_index = 0;
      check_cgns(
        cg_section_write(
          file.get(),
          base_index,
          zone_index,
          section_name.c_str(),
          element_type,
          next_element,
          next_element + section_size - 1,
          0,
          connectivity.data(),
          &section_index
        ),
        base::StatusCode::io_error,
        "CGNS export could not write an element section"
      );
      write_cgns_section_zone_id_descriptor(
        file.get(),
        base_index,
        zone_index,
        section_index,
        entity_group.zone_id()
      );
      if(order == EntityOrder::edge ||
         (order == EntityOrder::face && entity_group.is_boundary())) {
        boundary_records.push_back({
          order,
          std::string(entity_group.name()),
          cgns_entity_group_boundary_type_value(entity_group),
          next_element,
          next_element + section_size - 1,
        });
      }
      next_element += section_size;
    }
  };

  write_sections(EntityOrder::edge, CGNS_ENUMV(BAR_2));
  write_sections(EntityOrder::face, CGNS_ENUMV(TRI_3));
  write_sections(EntityOrder::cell, CGNS_ENUMV(TETRA_4));

  std::map<std::string, std::size_t> used_boundary_names;
  for(const auto &record : boundary_records) {
    write_cgns_boundary_condition(
      file.get(),
      base_index,
      zone_index,
      record,
      record.order == EntityOrder::edge
        ? CGNS_ENUMV(EdgeCenter)
        : (cgns_cell_dimension(domain) == 2
             ? CGNS_ENUMV(CellCenter)
             : CGNS_ENUMV(FaceCenter)),
      used_boundary_names
    );
  }
}

#endif

void validate_domain_for_msh_export(const Domain &domain)
{
  for(const auto &entity_group : domain.entity_groups()) {
    switch(entity_group.order()) {
    case EntityOrder::node:
      for(const auto &node : entity_group.nodes()) {
        if(node.header.kind != EntityKind::node_point) {
          fail(
            base::StatusCode::unsupported,
            "MSH export currently supports only point nodes."
          );
        }
      }
      break;
    case EntityOrder::edge:
      for(const auto &edge : entity_group.edges()) {
        if(edge.header.kind != EntityKind::edge_line) {
          fail(
            base::StatusCode::unsupported,
            "MSH export currently supports only line edges."
          );
        }
      }
      break;
    case EntityOrder::face:
      for(const auto &face : entity_group.faces()) {
        if(face.header.kind != EntityKind::face_triangle &&
           face.header.kind != EntityKind::face_quad) {
          fail(
            base::StatusCode::unsupported,
            "MSH export currently supports only triangle and quadrilateral faces."
          );
        }
      }
      break;
    case EntityOrder::cell:
      for(const auto &cell : entity_group.cells()) {
        if(cell.header.kind != EntityKind::cell_tetra &&
           cell.header.kind != EntityKind::cell_prism &&
           cell.header.kind != EntityKind::cell_hexa &&
           cell.header.kind != EntityKind::cell_pyramid) {
          fail(
            base::StatusCode::unsupported,
            "MSH export currently supports only tetrahedron, prism, hexahedron, and pyramid cells."
          );
        }
      }
      break;
    default:
      fail(
        base::StatusCode::unsupported,
        "MSH export found an unsupported entity order."
      );
    }
  }
}

[[nodiscard]] bool entity_group_has_msh_export_elements(const EntityGroup &entity_group) noexcept
{
  switch(entity_group.order()) {
  case EntityOrder::edge:
    return !entity_group.edges().empty();
  case EntityOrder::face:
    return !entity_group.faces().empty();
  case EntityOrder::cell:
    return !entity_group.cells().empty();
  default:
    return false;
  }
}

[[nodiscard]] int msh_dimension(EntityOrder order)
{
  switch(order) {
  case EntityOrder::edge:
    return 1;
  case EntityOrder::face:
    return 2;
  case EntityOrder::cell:
    return 3;
  default:
    fail(
      base::StatusCode::unsupported,
      "MSH export found an unsupported entity order."
    );
  }
}

[[nodiscard]] int msh_element_type(EntityKind kind)
{
  switch(kind) {
  case EntityKind::edge_line:
    return 1;    // Gmsh 2-node line
  case EntityKind::face_triangle:
    return 2;    // Gmsh 3-node triangle
  case EntityKind::face_quad:
    return 3;    // Gmsh 4-node quadrangle
  case EntityKind::cell_tetra:
    return 4;    // Gmsh 4-node tetrahedron
  case EntityKind::cell_hexa:
    return 5;    // Gmsh 8-node hexahedron
  case EntityKind::cell_prism:
    return 6;    // Gmsh 6-node prism
  case EntityKind::cell_pyramid:
    return 7;    // Gmsh 5-node pyramid
  default:
    fail(
      base::StatusCode::unsupported,
      "MSH export found an unsupported entity kind."
    );
  }
}

[[nodiscard]] int msh_element_type_for_order(EntityOrder order)
{
  switch(order) {
  case EntityOrder::edge:
    return 1;
  case EntityOrder::face:
    return 2;
  case EntityOrder::cell:
    return 4;
  default:
    fail(
      base::StatusCode::unsupported,
      "MSH export found an unsupported entity order."
    );
  }
}

struct MshEntityGroupExportRecord final {
  const EntityGroup *entity_group = nullptr;
  int dimension = 0;
  int element_type = 0;
  std::uint32_t entity_tag = 0U;
  bool has_physical = false;
  std::uint32_t physical_tag = 0U;
  std::string physical_name {};
};

[[nodiscard]] std::uint32_t allocate_msh_entity_tag(
  int dimension,
  std::uint32_t preferred_tag,
  std::map<int, std::set<std::uint32_t>> &used_tags,
  std::array<std::uint32_t, 4> &next_tags
)
{
  auto &used = used_tags[dimension];
  if(preferred_tag != 0U && preferred_tag != invalid_index &&
     used.emplace(preferred_tag).second) {
    next_tags[static_cast<std::size_t>(dimension)] =
      std::max(next_tags[static_cast<std::size_t>(dimension)], preferred_tag + 1U);
    return preferred_tag;
  }

  auto candidate = std::max<std::uint32_t>(
    1U,
    next_tags[static_cast<std::size_t>(dimension)]
  );
  while(!used.emplace(candidate).second) {
    ++candidate;
  }
  next_tags[static_cast<std::size_t>(dimension)] = candidate + 1U;
  return candidate;
}

[[nodiscard]] std::vector<MshEntityGroupExportRecord> build_msh_entity_group_export_records(
  const Domain &domain
)
{
  std::vector<MshEntityGroupExportRecord> records;
  std::map<int, std::set<std::uint32_t>> used_tags;
  std::array<std::uint32_t, 4> next_tags = {1U, 1U, 1U, 1U};

  for(const auto &entity_group : domain.entity_groups()) {
    if(entity_group.order() == EntityOrder::node || !entity_group_has_msh_export_elements(entity_group)) {
      continue;
    }

    MshEntityGroupExportRecord record;
    record.entity_group = &entity_group;
    record.dimension = msh_dimension(entity_group.order());
    record.element_type = (entity_group.default_kind() != EntityKind::invalid)
      ? msh_element_type(entity_group.default_kind())
      : msh_element_type_for_order(entity_group.order());
    record.has_physical = entity_group.zone_id() != 0U && entity_group.zone_id() != invalid_index;
    record.physical_tag = record.has_physical ? entity_group.zone_id() : 0U;
    const auto preferred_entity_tag =
      entity_group.source_entity_tag() != invalid_index ? entity_group.source_entity_tag() : record.physical_tag;
    record.entity_tag = allocate_msh_entity_tag(
      record.dimension,
      preferred_entity_tag,
      used_tags,
      next_tags
    );
    record.physical_name = sanitize_msh_name(entity_group.name());
    records.push_back(std::move(record));
  }

  return records;
}

[[nodiscard]] std::map<std::pair<int, std::uint32_t>, std::string>
build_msh_physical_names(const std::vector<MshEntityGroupExportRecord> &records)
{
  std::map<std::pair<int, std::uint32_t>, std::string> physical_names;
  for(const auto &record : records) {
    if(!record.has_physical) {
      continue;
    }
    const auto key = std::make_pair(record.dimension, record.physical_tag);
    const auto [entry_it, inserted] =
      physical_names.emplace(key, record.physical_name);
    if(!inserted && entry_it->second != record.physical_name) {
      fail(
        base::StatusCode::unsupported,
        "MSH export cannot assign two different physical names to the same dimension/tag pair (" +
          std::to_string(record.dimension) + ", " +
          std::to_string(record.physical_tag) + ")."
      );
    }
  }
  return physical_names;
}

struct MshBounds final {
  bool initialized = false;
  std::array<double, 3> min = {0.0, 0.0, 0.0};
  std::array<double, 3> max = {0.0, 0.0, 0.0};
};

void expand_msh_bounds(
  MshBounds &bounds,
  const Node &node
) noexcept
{
  if(!bounds.initialized) {
    bounds.initialized = true;
    bounds.min = node.coordinates;
    bounds.max = node.coordinates;
    return;
  }

  for(std::size_t axis = 0U; axis < 3U; ++axis) {
    bounds.min[axis] = std::min(bounds.min[axis], node.coordinates[axis]);
    bounds.max[axis] = std::max(bounds.max[axis], node.coordinates[axis]);
  }
}

[[nodiscard]] MshBounds compute_msh_entity_group_bounds(
  const Domain &domain,
  const EntityGroup &entity_group
)
{
  MshBounds bounds;

  switch(entity_group.order()) {
  case EntityOrder::edge:
    for(const auto &edge : entity_group.edges()) {
      for(const auto node_ref : domain.edge_nodes({edge.header.entity_group, edge.header.index})) {
        expand_msh_bounds(bounds, domain.node(node_ref));
      }
    }
    break;
  case EntityOrder::face:
    for(const auto &face : entity_group.faces()) {
      for(const auto node_ref : domain.face_nodes({face.header.entity_group, face.header.index})) {
        expand_msh_bounds(bounds, domain.node(node_ref));
      }
    }
    break;
  case EntityOrder::cell:
    for(const auto &cell : entity_group.cells()) {
      for(const auto node_ref : domain.cell_nodes({cell.header.entity_group, cell.header.index})) {
        expand_msh_bounds(bounds, domain.node(node_ref));
      }
    }
    break;
  default:
    break;
  }

  return bounds;
}

void write_msh22_ascii(
  const Domain &domain,
  std::string_view path,
  const MshExportOptions &options
)
{
  validate_domain_for_msh_export(domain);

  const auto node_indices = build_node_indices(domain);
  const auto entity_group_records = build_msh_entity_group_export_records(domain);
  const auto physical_names = build_msh_physical_names(entity_group_records);
  const std::string file_path(path);
  std::ofstream output(file_path, std::ios::trunc);
  if(!output) {
    fail(base::StatusCode::io_error, "MSH export could not open the output file.");
  }

  output << "$MeshFormat\n";
  output << "2.2 0 8\n";
  output << "$EndMeshFormat\n";

  if(options.write_physical_names && !physical_names.empty()) {
    output << "$PhysicalNames\n";
    output << physical_names.size() << '\n';
    for(const auto &entry : physical_names) {
      output << entry.first.first << ' ' << entry.first.second << " \""
             << entry.second << "\"\n";
    }
    output << "$EndPhysicalNames\n";
  }

  output << std::setprecision(17);
  output << "$Nodes\n";
  output << node_indices.size() << '\n';
  for(const auto &entity_group : domain.entity_groups()) {
    if(entity_group.order() != EntityOrder::node) {
      continue;
    }

    for(const auto &node : entity_group.nodes()) {
      const auto index_it = node_indices.find(
        pack_entity_ref({node.header.entity_group, node.header.index})
      );
      output << index_it->second << ' ' << node.coordinates[0] << ' '
             << node.coordinates[1] << ' ' << node.coordinates[2] << '\n';
    }
  }
  output << "$EndNodes\n";

  const auto summary = domain.summary();
  const auto element_count = summary.edge_count + summary.face_count + summary.cell_count;
  output << "$Elements\n";
  output << element_count << '\n';

  std::size_t next_element = 1U;
  for(int dimension = 1; dimension <= 3; ++dimension) {
    for(const auto &record : entity_group_records) {
      if(record.dimension != dimension) {
        continue;
      }

      const auto &entity_group = *record.entity_group;
      switch(entity_group.order()) {
      case EntityOrder::edge:
        for(const auto &edge : entity_group.edges()) {
          const auto nodes = domain.edge_nodes({edge.header.entity_group, edge.header.index});
          output << next_element++ << ' ' << record.element_type << " 2 "
                 << record.physical_tag << ' ' << record.entity_tag;
          for(const auto node_ref : nodes) {
            output << ' ' << node_indices.at(pack_entity_ref(node_ref));
          }
          output << '\n';
        }
        break;
      case EntityOrder::face:
        for(const auto &face : entity_group.faces()) {
          const auto nodes = domain.face_nodes({face.header.entity_group, face.header.index});
          output << next_element++ << ' ' << record.element_type << " 2 "
                 << record.physical_tag << ' ' << record.entity_tag;
          for(const auto node_ref : nodes) {
            output << ' ' << node_indices.at(pack_entity_ref(node_ref));
          }
          output << '\n';
        }
        break;
      case EntityOrder::cell:
        for(const auto &cell : entity_group.cells()) {
          const auto nodes = domain.cell_nodes({cell.header.entity_group, cell.header.index});
          output << next_element++ << ' ' << record.element_type << " 2 "
                 << record.physical_tag << ' ' << record.entity_tag;
          for(const auto node_ref : nodes) {
            output << ' ' << node_indices.at(pack_entity_ref(node_ref));
          }
          output << '\n';
        }
        break;
      default:
        break;
      }
    }
  }
  output << "$EndElements\n";
}

void write_msh22_binary(
  const Domain &domain,
  std::string_view path,
  const MshExportOptions &options
)
{
  validate_domain_for_msh_export(domain);

  const auto node_indices = build_node_indices(domain);
  const auto entity_group_records = build_msh_entity_group_export_records(domain);
  const auto physical_names = build_msh_physical_names(entity_group_records);
  const std::string file_path(path);
  std::ofstream output(file_path, std::ios::trunc | std::ios::binary);
  if(!output) {
    fail(base::StatusCode::io_error, "MSH export could not open the output file.");
  }

  output << "$MeshFormat\n";
  output << "2.2 1 " << sizeof(double) << '\n';
  const int one = 1;
  write_msh_binary_value(output, one);
  output.put('\n');
  output << "$EndMeshFormat\n";

  if(options.write_physical_names && !physical_names.empty()) {
    output << "$PhysicalNames\n";
    output << physical_names.size() << '\n';
    for(const auto &entry : physical_names) {
      output << entry.first.first << ' ' << entry.first.second << " \""
             << entry.second << "\"\n";
    }
    output << "$EndPhysicalNames\n";
  }

  output << "$Nodes\n";
  output << node_indices.size() << '\n';
  for(const auto &entity_group : domain.entity_groups()) {
    if(entity_group.order() != EntityOrder::node) {
      continue;
    }

    for(const auto &node : entity_group.nodes()) {
      const auto index_it = node_indices.find(
        pack_entity_ref({node.header.entity_group, node.header.index})
      );
      const auto node_tag = checked_msh_binary_int(
        index_it->second,
        "MSH export found a node tag outside the supported Gmsh binary integer range."
      );
      write_msh_binary_value(output, node_tag);
      write_msh_binary_values(output, node.coordinates.data(), node.coordinates.size());
    }
  }
  output.put('\n');
  output << "$EndNodes\n";

  const auto summary = domain.summary();
  const auto element_count = summary.edge_count + summary.face_count + summary.cell_count;
  output << "$Elements\n";
  output << element_count << '\n';

  std::size_t next_element = 1U;
  for(int dimension = 1; dimension <= 3; ++dimension) {
    for(const auto &record : entity_group_records) {
      if(record.dimension != dimension) {
        continue;
      }

      const auto &entity_group = *record.entity_group;
      const auto block_size =
        entity_group.order() == EntityOrder::edge
          ? entity_group.edges().size()
          : entity_group.order() == EntityOrder::face ? entity_group.faces().size() : entity_group.cells().size();
      const std::array<int, 3> header = {
        record.element_type,
        checked_msh_binary_int(
          block_size,
          "MSH export found an element block outside the supported Gmsh binary integer range."
        ),
        2,
      };
      write_msh_binary_values(output, header.data(), header.size());

      const auto physical_tag = checked_msh_binary_int(
        record.physical_tag,
        "MSH export found a physical tag outside the supported Gmsh binary integer range."
      );
      const auto entity_tag = checked_msh_binary_int(
        record.entity_tag,
        "MSH export found an entity tag outside the supported Gmsh binary integer range."
      );

      switch(entity_group.order()) {
      case EntityOrder::edge:
        for(const auto &edge : entity_group.edges()) {
          const auto nodes = domain.edge_nodes({edge.header.entity_group, edge.header.index});
          const std::array<int, 5> data = {
            checked_msh_binary_int(
              next_element++,
              "MSH export found an element tag outside the supported Gmsh binary integer range."
            ),
            physical_tag,
            entity_tag,
            checked_msh_binary_int(
              node_indices.at(pack_entity_ref(nodes[0])),
              "MSH export found a node tag outside the supported Gmsh binary integer range."
            ),
            checked_msh_binary_int(
              node_indices.at(pack_entity_ref(nodes[1])),
              "MSH export found a node tag outside the supported Gmsh binary integer range."
            ),
          };
          write_msh_binary_values(output, data.data(), data.size());
        }
        break;
      case EntityOrder::face:
        for(const auto &face : entity_group.faces()) {
          const auto nodes = domain.face_nodes({face.header.entity_group, face.header.index});
          const std::array<int, 6> data = {
            checked_msh_binary_int(
              next_element++,
              "MSH export found an element tag outside the supported Gmsh binary integer range."
            ),
            physical_tag,
            entity_tag,
            checked_msh_binary_int(
              node_indices.at(pack_entity_ref(nodes[0])),
              "MSH export found a node tag outside the supported Gmsh binary integer range."
            ),
            checked_msh_binary_int(
              node_indices.at(pack_entity_ref(nodes[1])),
              "MSH export found a node tag outside the supported Gmsh binary integer range."
            ),
            checked_msh_binary_int(
              node_indices.at(pack_entity_ref(nodes[2])),
              "MSH export found a node tag outside the supported Gmsh binary integer range."
            ),
          };
          write_msh_binary_values(output, data.data(), data.size());
        }
        break;
      case EntityOrder::cell:
        for(const auto &cell : entity_group.cells()) {
          const auto nodes = domain.cell_nodes({cell.header.entity_group, cell.header.index});
          const std::array<int, 7> data = {
            checked_msh_binary_int(
              next_element++,
              "MSH export found an element tag outside the supported Gmsh binary integer range."
            ),
            physical_tag,
            entity_tag,
            checked_msh_binary_int(
              node_indices.at(pack_entity_ref(nodes[0])),
              "MSH export found a node tag outside the supported Gmsh binary integer range."
            ),
            checked_msh_binary_int(
              node_indices.at(pack_entity_ref(nodes[1])),
              "MSH export found a node tag outside the supported Gmsh binary integer range."
            ),
            checked_msh_binary_int(
              node_indices.at(pack_entity_ref(nodes[2])),
              "MSH export found a node tag outside the supported Gmsh binary integer range."
            ),
            checked_msh_binary_int(
              node_indices.at(pack_entity_ref(nodes[3])),
              "MSH export found a node tag outside the supported Gmsh binary integer range."
            ),
          };
          write_msh_binary_values(output, data.data(), data.size());
        }
        break;
      default:
        break;
      }
    }
  }
  output.put('\n');
  output << "$EndElements\n";
}

void write_msh41_ascii(
  const Domain &domain,
  std::string_view path,
  const MshExportOptions &options
)
{
  validate_domain_for_msh_export(domain);

  struct MshNodeBlockKey final {
    int dimension = 0;
    std::uint32_t entity_tag = 0U;

    [[nodiscard]] bool operator<(const MshNodeBlockKey &other) const noexcept
    {
      if(dimension != other.dimension) {
        return dimension < other.dimension;
      }
      return entity_tag < other.entity_tag;
    }
  };

  struct MshPointEntity final {
    std::uint32_t entity_tag = 0U;
    EntityRef node_ref {};
  };

  const auto node_indices = build_node_indices(domain);
  const auto entity_group_records = build_msh_entity_group_export_records(domain);
  const auto physical_names = build_msh_physical_names(entity_group_records);

  std::map<std::uint64_t, MshNodeBlockKey> node_owners;
  const auto claim_nodes = [&](const MshEntityGroupExportRecord &record) {
    const auto &entity_group = *record.entity_group;
    switch(entity_group.order()) {
    case EntityOrder::edge:
      for(const auto &edge : entity_group.edges()) {
        for(const auto node_ref : domain.edge_nodes({edge.header.entity_group, edge.header.index})) {
          node_owners.emplace(
            pack_entity_ref(node_ref),
            MshNodeBlockKey {record.dimension, record.entity_tag}
          );
        }
      }
      break;
    case EntityOrder::face:
      for(const auto &face : entity_group.faces()) {
        for(const auto node_ref : domain.face_nodes({face.header.entity_group, face.header.index})) {
          node_owners.emplace(
            pack_entity_ref(node_ref),
            MshNodeBlockKey {record.dimension, record.entity_tag}
          );
        }
      }
      break;
    case EntityOrder::cell:
      for(const auto &cell : entity_group.cells()) {
        for(const auto node_ref : domain.cell_nodes({cell.header.entity_group, cell.header.index})) {
          node_owners.emplace(
            pack_entity_ref(node_ref),
            MshNodeBlockKey {record.dimension, record.entity_tag}
          );
        }
      }
      break;
    default:
      break;
    }
  };

  for(const auto order : {EntityOrder::cell, EntityOrder::face, EntityOrder::edge}) {
    for(const auto &record : entity_group_records) {
      if(record.entity_group->order() == order) {
        claim_nodes(record);
      }
    }
  }

  std::vector<MshPointEntity> point_entities;
  std::uint32_t next_point_entity_tag = 1U;
  for(const auto &entity_group : domain.entity_groups()) {
    if(entity_group.order() != EntityOrder::node) {
      continue;
    }
    for(const auto &node : entity_group.nodes()) {
      const auto node_ref = EntityRef {node.header.entity_group, node.header.index};
      const auto node_key = pack_entity_ref(node_ref);
      if(node_owners.find(node_key) != node_owners.end()) {
        continue;
      }

      const auto point_tag = next_point_entity_tag++;
      point_entities.push_back({point_tag, node_ref});
      node_owners.emplace(node_key, MshNodeBlockKey {0, point_tag});
    }
  }

  std::map<MshNodeBlockKey, std::vector<EntityRef>> node_blocks;
  for(const auto &entity_group : domain.entity_groups()) {
    if(entity_group.order() != EntityOrder::node) {
      continue;
    }
    for(const auto &node : entity_group.nodes()) {
      const auto node_ref = EntityRef {node.header.entity_group, node.header.index};
      node_blocks[node_owners.at(pack_entity_ref(node_ref))].push_back(node_ref);
    }
  }

  const std::string file_path(path);
  std::ofstream output(file_path, std::ios::trunc);
  if(!output) {
    fail(base::StatusCode::io_error, "MSH export could not open the output file.");
  }

  output << "$MeshFormat\n";
  output << "4.1 0 8\n";
  output << "$EndMeshFormat\n";

  if(options.write_physical_names && !physical_names.empty()) {
    output << "$PhysicalNames\n";
    output << physical_names.size() << '\n';
    for(const auto &entry : physical_names) {
      output << entry.first.first << ' ' << entry.first.second << " \""
             << entry.second << "\"\n";
    }
    output << "$EndPhysicalNames\n";
  }

  output << std::setprecision(17);
  output << "$Entities\n";
  output << point_entities.size() << ' '
         << std::count_if(
              entity_group_records.begin(),
              entity_group_records.end(),
              [](const MshEntityGroupExportRecord &record) { return record.dimension == 1; }
            )
         << ' '
         << std::count_if(
              entity_group_records.begin(),
              entity_group_records.end(),
              [](const MshEntityGroupExportRecord &record) { return record.dimension == 2; }
            )
         << ' '
         << std::count_if(
              entity_group_records.begin(),
              entity_group_records.end(),
              [](const MshEntityGroupExportRecord &record) { return record.dimension == 3; }
            )
         << '\n';

  for(const auto &point : point_entities) {
    const auto &node = domain.node(point.node_ref);
    output << point.entity_tag << ' ' << node.coordinates[0] << ' '
           << node.coordinates[1] << ' ' << node.coordinates[2] << " 0\n";
  }

  for(int dimension = 1; dimension <= 3; ++dimension) {
    for(const auto &record : entity_group_records) {
      if(record.dimension != dimension) {
        continue;
      }

      const auto bounds = compute_msh_entity_group_bounds(domain, *record.entity_group);
      output << record.entity_tag << ' ' << bounds.min[0] << ' ' << bounds.min[1]
             << ' ' << bounds.min[2] << ' ' << bounds.max[0] << ' '
             << bounds.max[1] << ' ' << bounds.max[2] << ' '
             << (record.has_physical ? 1U : 0U);
      if(record.has_physical) {
        output << ' ' << record.physical_tag;
      }
      output << " 0\n";
    }
  }
  output << "$EndEntities\n";

  const auto max_node_tag = node_indices.empty() ? 0U : node_indices.size();
  output << "$Nodes\n";
  output << node_blocks.size() << ' ' << node_indices.size() << ' '
         << (node_indices.empty() ? 0U : 1U) << ' ' << max_node_tag << '\n';
  for(const auto &block : node_blocks) {
    output << block.first.dimension << ' ' << block.first.entity_tag << " 0 "
           << block.second.size() << '\n';
    for(const auto node_ref : block.second) {
      output << node_indices.at(pack_entity_ref(node_ref)) << '\n';
    }
    for(const auto node_ref : block.second) {
      const auto &node = domain.node(node_ref);
      output << node.coordinates[0] << ' ' << node.coordinates[1] << ' '
             << node.coordinates[2] << '\n';
    }
  }
  output << "$EndNodes\n";

  const auto summary = domain.summary();
  const auto element_count = summary.edge_count + summary.face_count + summary.cell_count;
  output << "$Elements\n";
  output << entity_group_records.size() << ' ' << element_count << ' '
         << (element_count == 0U ? 0U : 1U) << ' ' << element_count << '\n';

  std::size_t next_element = 1U;
  for(int dimension = 1; dimension <= 3; ++dimension) {
    for(const auto &record : entity_group_records) {
      if(record.dimension != dimension) {
        continue;
      }

      const auto &entity_group = *record.entity_group;
      const auto block_size =
        entity_group.order() == EntityOrder::edge
          ? entity_group.edges().size()
          : entity_group.order() == EntityOrder::face ? entity_group.faces().size() : entity_group.cells().size();
      output << record.dimension << ' ' << record.entity_tag << ' '
             << record.element_type << ' ' << block_size << '\n';

      switch(entity_group.order()) {
      case EntityOrder::edge:
        for(const auto &edge : entity_group.edges()) {
          const auto nodes = domain.edge_nodes({edge.header.entity_group, edge.header.index});
          output << next_element++;
          for(const auto node_ref : nodes) {
            output << ' ' << node_indices.at(pack_entity_ref(node_ref));
          }
          output << '\n';
        }
        break;
      case EntityOrder::face:
        for(const auto &face : entity_group.faces()) {
          const auto nodes = domain.face_nodes({face.header.entity_group, face.header.index});
          output << next_element++;
          for(const auto node_ref : nodes) {
            output << ' ' << node_indices.at(pack_entity_ref(node_ref));
          }
          output << '\n';
        }
        break;
      case EntityOrder::cell:
        for(const auto &cell : entity_group.cells()) {
          const auto nodes = domain.cell_nodes({cell.header.entity_group, cell.header.index});
          output << next_element++;
          for(const auto node_ref : nodes) {
            output << ' ' << node_indices.at(pack_entity_ref(node_ref));
          }
          output << '\n';
        }
        break;
      default:
        break;
      }
    }
  }
  output << "$EndElements\n";
}

void write_msh41_binary(
  const Domain &domain,
  std::string_view path,
  const MshExportOptions &options
)
{
  validate_domain_for_msh_export(domain);

  struct MshNodeBlockKey final {
    int dimension = 0;
    std::uint32_t entity_tag = 0U;

    [[nodiscard]] bool operator<(const MshNodeBlockKey &other) const noexcept
    {
      if(dimension != other.dimension) {
        return dimension < other.dimension;
      }
      return entity_tag < other.entity_tag;
    }
  };

  struct MshPointEntity final {
    std::uint32_t entity_tag = 0U;
    EntityRef node_ref {};
  };

  const auto node_indices = build_node_indices(domain);
  const auto entity_group_records = build_msh_entity_group_export_records(domain);
  const auto physical_names = build_msh_physical_names(entity_group_records);

  std::map<std::uint64_t, MshNodeBlockKey> node_owners;
  const auto claim_nodes = [&](const MshEntityGroupExportRecord &record) {
    const auto &entity_group = *record.entity_group;
    switch(entity_group.order()) {
    case EntityOrder::edge:
      for(const auto &edge : entity_group.edges()) {
        for(const auto node_ref : domain.edge_nodes({edge.header.entity_group, edge.header.index})) {
          node_owners.emplace(
            pack_entity_ref(node_ref),
            MshNodeBlockKey {record.dimension, record.entity_tag}
          );
        }
      }
      break;
    case EntityOrder::face:
      for(const auto &face : entity_group.faces()) {
        for(const auto node_ref : domain.face_nodes({face.header.entity_group, face.header.index})) {
          node_owners.emplace(
            pack_entity_ref(node_ref),
            MshNodeBlockKey {record.dimension, record.entity_tag}
          );
        }
      }
      break;
    case EntityOrder::cell:
      for(const auto &cell : entity_group.cells()) {
        for(const auto node_ref : domain.cell_nodes({cell.header.entity_group, cell.header.index})) {
          node_owners.emplace(
            pack_entity_ref(node_ref),
            MshNodeBlockKey {record.dimension, record.entity_tag}
          );
        }
      }
      break;
    default:
      break;
    }
  };

  for(const auto order : {EntityOrder::cell, EntityOrder::face, EntityOrder::edge}) {
    for(const auto &record : entity_group_records) {
      if(record.entity_group->order() == order) {
        claim_nodes(record);
      }
    }
  }

  std::vector<MshPointEntity> point_entities;
  std::uint32_t next_point_entity_tag = 1U;
  for(const auto &entity_group : domain.entity_groups()) {
    if(entity_group.order() != EntityOrder::node) {
      continue;
    }
    for(const auto &node : entity_group.nodes()) {
      const auto node_ref = EntityRef {node.header.entity_group, node.header.index};
      const auto node_key = pack_entity_ref(node_ref);
      if(node_owners.find(node_key) != node_owners.end()) {
        continue;
      }

      const auto point_tag = next_point_entity_tag++;
      point_entities.push_back({point_tag, node_ref});
      node_owners.emplace(node_key, MshNodeBlockKey {0, point_tag});
    }
  }

  std::map<MshNodeBlockKey, std::vector<EntityRef>> node_blocks;
  for(const auto &entity_group : domain.entity_groups()) {
    if(entity_group.order() != EntityOrder::node) {
      continue;
    }
    for(const auto &node : entity_group.nodes()) {
      const auto node_ref = EntityRef {node.header.entity_group, node.header.index};
      node_blocks[node_owners.at(pack_entity_ref(node_ref))].push_back(node_ref);
    }
  }

  const std::string file_path(path);
  std::ofstream output(file_path, std::ios::trunc | std::ios::binary);
  if(!output) {
    fail(base::StatusCode::io_error, "MSH export could not open the output file.");
  }

  output << "$MeshFormat\n";
  output << "4.1 1 " << sizeof(std::size_t) << '\n';
  const int one = 1;
  write_msh_binary_value(output, one);
  output.put('\n');
  output << "$EndMeshFormat\n";

  if(options.write_physical_names && !physical_names.empty()) {
    output << "$PhysicalNames\n";
    output << physical_names.size() << '\n';
    for(const auto &entry : physical_names) {
      output << entry.first.first << ' ' << entry.first.second << " \""
             << entry.second << "\"\n";
    }
    output << "$EndPhysicalNames\n";
  }

  output << "$Entities\n";
  const std::array<std::size_t, 4> entity_counts = {
    point_entities.size(),
    static_cast<std::size_t>(std::count_if(
      entity_group_records.begin(),
      entity_group_records.end(),
      [](const MshEntityGroupExportRecord &record) { return record.dimension == 1; }
    )),
    static_cast<std::size_t>(std::count_if(
      entity_group_records.begin(),
      entity_group_records.end(),
      [](const MshEntityGroupExportRecord &record) { return record.dimension == 2; }
    )),
    static_cast<std::size_t>(std::count_if(
      entity_group_records.begin(),
      entity_group_records.end(),
      [](const MshEntityGroupExportRecord &record) { return record.dimension == 3; }
    )),
  };
  write_msh_binary_values(output, entity_counts.data(), entity_counts.size());

  for(const auto &point : point_entities) {
    const auto entity_tag = checked_msh_binary_int(
      point.entity_tag,
      "MSH export found an entity tag outside the supported Gmsh binary integer range."
    );
    write_msh_binary_value(output, entity_tag);

    const auto &node = domain.node(point.node_ref);
    write_msh_binary_values(output, node.coordinates.data(), node.coordinates.size());

    const std::size_t physical_count = 0U;
    write_msh_binary_value(output, physical_count);
  }

  for(int dimension = 1; dimension <= 3; ++dimension) {
    for(const auto &record : entity_group_records) {
      if(record.dimension != dimension) {
        continue;
      }

      const auto entity_tag = checked_msh_binary_int(
        record.entity_tag,
        "MSH export found an entity tag outside the supported Gmsh binary integer range."
      );
      write_msh_binary_value(output, entity_tag);

      const auto bounds = compute_msh_entity_group_bounds(domain, *record.entity_group);
      const std::array<double, 6> box = {
        bounds.min[0],
        bounds.min[1],
        bounds.min[2],
        bounds.max[0],
        bounds.max[1],
        bounds.max[2],
      };
      write_msh_binary_values(output, box.data(), box.size());

      const std::size_t physical_count = record.has_physical ? 1U : 0U;
      write_msh_binary_value(output, physical_count);
      if(record.has_physical) {
        const auto physical_tag = checked_msh_binary_int(
          record.physical_tag,
          "MSH export found a physical tag outside the supported Gmsh binary integer range."
        );
        write_msh_binary_value(output, physical_tag);
      }

      const std::size_t boundary_count = 0U;
      write_msh_binary_value(output, boundary_count);
    }
  }
  output.put('\n');
  output << "$EndEntities\n";

  const auto max_node_tag = node_indices.empty() ? 0U : node_indices.size();
  output << "$Nodes\n";
  const std::array<std::size_t, 4> node_header = {
    node_blocks.size(),
    node_indices.size(),
    node_indices.empty() ? 0U : 1U,
    max_node_tag,
  };
  write_msh_binary_values(output, node_header.data(), node_header.size());
  for(const auto &block : node_blocks) {
    const std::array<int, 3> block_header = {
      block.first.dimension,
      checked_msh_binary_int(
        block.first.entity_tag,
        "MSH export found an entity tag outside the supported Gmsh binary integer range."
      ),
      0,
    };
    write_msh_binary_values(output, block_header.data(), block_header.size());

    const auto block_size = block.second.size();
    write_msh_binary_value(output, block_size);
    for(const auto node_ref : block.second) {
      const auto node_tag = node_indices.at(pack_entity_ref(node_ref));
      write_msh_binary_value(output, node_tag);
    }
    for(const auto node_ref : block.second) {
      const auto &node = domain.node(node_ref);
      write_msh_binary_values(output, node.coordinates.data(), node.coordinates.size());
    }
  }
  output.put('\n');
  output << "$EndNodes\n";

  const auto summary = domain.summary();
  const auto element_count = summary.edge_count + summary.face_count + summary.cell_count;
  output << "$Elements\n";
  const std::array<std::size_t, 4> element_header = {
    entity_group_records.size(),
    element_count,
    element_count == 0U ? 0U : 1U,
    element_count,
  };
  write_msh_binary_values(output, element_header.data(), element_header.size());

  std::size_t next_element = 1U;
  for(int dimension = 1; dimension <= 3; ++dimension) {
    for(const auto &record : entity_group_records) {
      if(record.dimension != dimension) {
        continue;
      }

      const auto &entity_group = *record.entity_group;
      const auto block_size =
        entity_group.order() == EntityOrder::edge
          ? entity_group.edges().size()
          : entity_group.order() == EntityOrder::face ? entity_group.faces().size() : entity_group.cells().size();
      const std::array<int, 3> block_header = {
        record.dimension,
        checked_msh_binary_int(
          record.entity_tag,
          "MSH export found an entity tag outside the supported Gmsh binary integer range."
        ),
        record.element_type,
      };
      write_msh_binary_values(output, block_header.data(), block_header.size());
      write_msh_binary_value(output, block_size);

      switch(entity_group.order()) {
      case EntityOrder::edge:
        for(const auto &edge : entity_group.edges()) {
          const auto nodes = domain.edge_nodes({edge.header.entity_group, edge.header.index});
          const std::array<std::size_t, 3> data = {
            next_element++,
            node_indices.at(pack_entity_ref(nodes[0])),
            node_indices.at(pack_entity_ref(nodes[1])),
          };
          write_msh_binary_values(output, data.data(), data.size());
        }
        break;
      case EntityOrder::face:
        for(const auto &face : entity_group.faces()) {
          const auto nodes = domain.face_nodes({face.header.entity_group, face.header.index});
          const std::array<std::size_t, 4> data = {
            next_element++,
            node_indices.at(pack_entity_ref(nodes[0])),
            node_indices.at(pack_entity_ref(nodes[1])),
            node_indices.at(pack_entity_ref(nodes[2])),
          };
          write_msh_binary_values(output, data.data(), data.size());
        }
        break;
      case EntityOrder::cell:
        for(const auto &cell : entity_group.cells()) {
          const auto nodes = domain.cell_nodes({cell.header.entity_group, cell.header.index});
          const std::array<std::size_t, 5> data = {
            next_element++,
            node_indices.at(pack_entity_ref(nodes[0])),
            node_indices.at(pack_entity_ref(nodes[1])),
            node_indices.at(pack_entity_ref(nodes[2])),
            node_indices.at(pack_entity_ref(nodes[3])),
          };
          write_msh_binary_values(output, data.data(), data.size());
        }
        break;
      default:
        break;
      }
    }
  }
  output.put('\n');
  output << "$EndElements\n";
}

void validate_domain_for_nastran_export(const Domain &domain)
{
  for(const auto &entity_group : domain.entity_groups()) {
    switch(entity_group.order()) {
    case EntityOrder::node:
      for(const auto &node : entity_group.nodes()) {
        if(node.header.kind != EntityKind::node_point) {
          fail(
            base::StatusCode::unsupported,
            "NASTRAN export currently supports only point nodes."
          );
        }
      }
      break;
    case EntityOrder::edge:
      for(const auto &edge : entity_group.edges()) {
        if(edge.header.kind != EntityKind::edge_line) {
          fail(
            base::StatusCode::unsupported,
            "NASTRAN export currently supports only CROD-compatible line edges."
          );
        }
      }
      break;
    case EntityOrder::face:
      for(const auto &face : entity_group.faces()) {
        if(face.header.kind != EntityKind::face_triangle) {
          fail(
            base::StatusCode::unsupported,
            "NASTRAN export currently supports only CTRIA3-compatible triangle faces."
          );
        }
      }
      break;
    case EntityOrder::cell:
      for(const auto &cell : entity_group.cells()) {
        if(cell.header.kind != EntityKind::cell_tetra) {
          fail(
            base::StatusCode::unsupported,
            "NASTRAN export currently supports only first-order CTETRA cells."
          );
        }
      }
      break;
    default:
      fail(
        base::StatusCode::unsupported,
        "NASTRAN export found an unsupported entity order."
      );
    }
  }
}

[[nodiscard]] bool has_nastran_entity_group_prefix(
  std::string_view name,
  std::string_view prefix
) noexcept
{
  return name.size() >= prefix.size() && name.compare(0U, prefix.size(), prefix) == 0;
}

[[nodiscard]] bool is_nastran_cquad4_entity_group(const EntityGroup &entity_group) noexcept
{
  return entity_group.order() == EntityOrder::face &&
         ((entity_group.import_info().format == EntityGroupImportFormat::nastran &&
           entity_group.import_info().nastran.source_card ==
             NastranEntityGroupSourceCard::cquad4) ||
          has_nastran_entity_group_prefix(entity_group.name(), "cquad4_pid_"));
}

[[nodiscard]] bool try_build_nastran_cquad4(
  const Domain &domain,
  EntityRef first_face,
  EntityRef second_face,
  std::array<EntityRef, 4> &quad_nodes
) noexcept
{
  const auto first_source_tag = domain.face_source_entity_tag(first_face);
  if(first_source_tag == invalid_index ||
     first_source_tag != domain.face_source_entity_tag(second_face)) {
    return false;
  }

  const auto first_nodes = domain.face_nodes(first_face);
  const auto second_nodes = domain.face_nodes(second_face);
  if(first_nodes.size != 3U || second_nodes.size != 3U) {
    return false;
  }

  const auto g1 = first_nodes[0];
  const auto g2 = first_nodes[1];
  const auto g3 = first_nodes[2];
  if(g1 == g2 || g2 == g3 || g1 == g3) {
    return false;
  }

  bool saw_g1 = false;
  bool saw_g3 = false;
  EntityRef g4 {};
  bool saw_g4 = false;
  for(const auto candidate : second_nodes) {
    if(candidate == g1) {
      saw_g1 = true;
      continue;
    }
    if(candidate == g3) {
      saw_g3 = true;
      continue;
    }
    if(candidate == g2 || saw_g4) {
      return false;
    }
    g4 = candidate;
    saw_g4 = true;
  }

  if(!saw_g1 || !saw_g3 || !saw_g4 || g4 == g1 || g4 == g2 || g4 == g3) {
    return false;
  }

  quad_nodes = {g1, g2, g3, g4};
  return true;
}

void write_nastran_bulk_data(
  const Domain &domain,
  std::string_view path,
  const NastranExportOptions &options
)
{
  validate_domain_for_nastran_export(domain);

  const auto node_indices = build_node_indices(domain);
  const std::string file_path(path);
  std::ofstream output(file_path, std::ios::trunc);
  if(!output) {
    fail(base::StatusCode::io_error, "NASTRAN export could not open the output file.");
  }

  output << "$ SQMesh NASTRAN bulk-data export\n";
  if(options.write_begin_bulk) {
    output << "BEGIN BULK\n";
  }

  const auto write_grid = [&](std::uint64_t node_id, const Node &node) {
    if(options.field_format == NastranFieldFormat::free_field) {
      output << std::setprecision(17);
      output << "GRID," << node_id << ",," << node.coordinates[0] << ','
             << node.coordinates[1] << ',' << node.coordinates[2] << '\n';
      return;
    }

    const std::vector<std::string> fields = {
      options.field_format == NastranFieldFormat::large_fixed
        ? format_nastran_large_field_integer(node_id, "NASTRAN GRID id")
        : format_nastran_small_field_integer(node_id, "NASTRAN GRID id"),
      "",
      options.field_format == NastranFieldFormat::large_fixed
        ? format_nastran_large_field_real(node.coordinates[0], "NASTRAN GRID X")
        : format_nastran_small_field_real(node.coordinates[0], "NASTRAN GRID X"),
      options.field_format == NastranFieldFormat::large_fixed
        ? format_nastran_large_field_real(node.coordinates[1], "NASTRAN GRID Y")
        : format_nastran_small_field_real(node.coordinates[1], "NASTRAN GRID Y"),
      options.field_format == NastranFieldFormat::large_fixed
        ? format_nastran_large_field_real(node.coordinates[2], "NASTRAN GRID Z")
        : format_nastran_small_field_real(node.coordinates[2], "NASTRAN GRID Z"),
    };
    if(options.field_format == NastranFieldFormat::large_fixed) {
      write_nastran_large_field_record(output, "GRID", fields);
      return;
    }

    write_nastran_small_field_record(output, "GRID", fields);
  };

  const auto write_edge = [&](std::uint64_t element_id, std::uint32_t property_id, EntityRef edge_ref) {
    const auto nodes = domain.edge_nodes(edge_ref);
    if(options.field_format == NastranFieldFormat::free_field) {
      output << "CROD," << element_id << ',' << property_id;
      for(const auto node_ref : nodes) {
        output << ',' << node_indices.at(pack_entity_ref(node_ref));
      }
      output << '\n';
      return;
    }

    const std::vector<std::string> fields = {
      options.field_format == NastranFieldFormat::large_fixed
        ? format_nastran_large_field_integer(element_id, "NASTRAN CROD id")
        : format_nastran_small_field_integer(element_id, "NASTRAN CROD id"),
      options.field_format == NastranFieldFormat::large_fixed
        ? format_nastran_large_field_integer(property_id, "NASTRAN CROD PID")
        : format_nastran_small_field_integer(property_id, "NASTRAN CROD PID"),
      options.field_format == NastranFieldFormat::large_fixed
        ? format_nastran_large_field_integer(
            node_indices.at(pack_entity_ref(nodes[0])),
            "NASTRAN CROD G1"
          )
        : format_nastran_small_field_integer(
            node_indices.at(pack_entity_ref(nodes[0])),
            "NASTRAN CROD G1"
          ),
      options.field_format == NastranFieldFormat::large_fixed
        ? format_nastran_large_field_integer(
            node_indices.at(pack_entity_ref(nodes[1])),
            "NASTRAN CROD G2"
          )
        : format_nastran_small_field_integer(
            node_indices.at(pack_entity_ref(nodes[1])),
            "NASTRAN CROD G2"
          ),
    };
    if(options.field_format == NastranFieldFormat::large_fixed) {
      write_nastran_large_field_record(output, "CROD", fields);
      return;
    }

    write_nastran_small_field_record(output, "CROD", fields);
  };

  const auto write_face = [&](std::uint64_t element_id, std::uint32_t property_id, EntityRef face_ref) {
    const auto nodes = domain.face_nodes(face_ref);
    if(options.field_format == NastranFieldFormat::free_field) {
      output << "CTRIA3," << element_id << ',' << property_id;
      for(const auto node_ref : nodes) {
        output << ',' << node_indices.at(pack_entity_ref(node_ref));
      }
      output << '\n';
      return;
    }

    const std::vector<std::string> fields = {
      options.field_format == NastranFieldFormat::large_fixed
        ? format_nastran_large_field_integer(element_id, "NASTRAN CTRIA3 id")
        : format_nastran_small_field_integer(element_id, "NASTRAN CTRIA3 id"),
      options.field_format == NastranFieldFormat::large_fixed
        ? format_nastran_large_field_integer(property_id, "NASTRAN CTRIA3 PID")
        : format_nastran_small_field_integer(property_id, "NASTRAN CTRIA3 PID"),
      options.field_format == NastranFieldFormat::large_fixed
        ? format_nastran_large_field_integer(
            node_indices.at(pack_entity_ref(nodes[0])),
            "NASTRAN CTRIA3 G1"
          )
        : format_nastran_small_field_integer(
            node_indices.at(pack_entity_ref(nodes[0])),
            "NASTRAN CTRIA3 G1"
          ),
      options.field_format == NastranFieldFormat::large_fixed
        ? format_nastran_large_field_integer(
            node_indices.at(pack_entity_ref(nodes[1])),
            "NASTRAN CTRIA3 G2"
          )
        : format_nastran_small_field_integer(
            node_indices.at(pack_entity_ref(nodes[1])),
            "NASTRAN CTRIA3 G2"
          ),
      options.field_format == NastranFieldFormat::large_fixed
        ? format_nastran_large_field_integer(
            node_indices.at(pack_entity_ref(nodes[2])),
            "NASTRAN CTRIA3 G3"
          )
        : format_nastran_small_field_integer(
            node_indices.at(pack_entity_ref(nodes[2])),
            "NASTRAN CTRIA3 G3"
          ),
    };
    if(options.field_format == NastranFieldFormat::large_fixed) {
      write_nastran_large_field_record(output, "CTRIA3", fields);
      return;
    }

    write_nastran_small_field_record(output, "CTRIA3", fields);
  };

  const auto write_quad =
    [&](std::uint64_t element_id, std::uint32_t property_id, const std::array<EntityRef, 4> &nodes) {
      if(options.field_format == NastranFieldFormat::free_field) {
        output << "CQUAD4," << element_id << ',' << property_id;
        for(const auto node_ref : nodes) {
          output << ',' << node_indices.at(pack_entity_ref(node_ref));
        }
        output << '\n';
        return;
      }

      const std::vector<std::string> fields = {
        options.field_format == NastranFieldFormat::large_fixed
          ? format_nastran_large_field_integer(element_id, "NASTRAN CQUAD4 id")
          : format_nastran_small_field_integer(element_id, "NASTRAN CQUAD4 id"),
        options.field_format == NastranFieldFormat::large_fixed
          ? format_nastran_large_field_integer(property_id, "NASTRAN CQUAD4 PID")
          : format_nastran_small_field_integer(property_id, "NASTRAN CQUAD4 PID"),
        options.field_format == NastranFieldFormat::large_fixed
          ? format_nastran_large_field_integer(
              node_indices.at(pack_entity_ref(nodes[0])),
              "NASTRAN CQUAD4 G1"
            )
          : format_nastran_small_field_integer(
              node_indices.at(pack_entity_ref(nodes[0])),
              "NASTRAN CQUAD4 G1"
            ),
        options.field_format == NastranFieldFormat::large_fixed
          ? format_nastran_large_field_integer(
              node_indices.at(pack_entity_ref(nodes[1])),
              "NASTRAN CQUAD4 G2"
            )
          : format_nastran_small_field_integer(
              node_indices.at(pack_entity_ref(nodes[1])),
              "NASTRAN CQUAD4 G2"
            ),
        options.field_format == NastranFieldFormat::large_fixed
          ? format_nastran_large_field_integer(
              node_indices.at(pack_entity_ref(nodes[2])),
              "NASTRAN CQUAD4 G3"
            )
          : format_nastran_small_field_integer(
              node_indices.at(pack_entity_ref(nodes[2])),
              "NASTRAN CQUAD4 G3"
            ),
        options.field_format == NastranFieldFormat::large_fixed
          ? format_nastran_large_field_integer(
              node_indices.at(pack_entity_ref(nodes[3])),
              "NASTRAN CQUAD4 G4"
            )
          : format_nastran_small_field_integer(
              node_indices.at(pack_entity_ref(nodes[3])),
              "NASTRAN CQUAD4 G4"
            ),
      };
      if(options.field_format == NastranFieldFormat::large_fixed) {
        write_nastran_large_field_record(output, "CQUAD4", fields);
        return;
      }

      write_nastran_small_field_record(output, "CQUAD4", fields);
    };

  const auto write_cell = [&](std::uint64_t element_id, std::uint32_t property_id, EntityRef cell_ref) {
    const auto nodes = domain.cell_nodes(cell_ref);
    if(options.field_format == NastranFieldFormat::free_field) {
      output << "CTETRA," << element_id << ',' << property_id;
      for(const auto node_ref : nodes) {
        output << ',' << node_indices.at(pack_entity_ref(node_ref));
      }
      output << '\n';
      return;
    }

    const std::vector<std::string> fields = {
      options.field_format == NastranFieldFormat::large_fixed
        ? format_nastran_large_field_integer(element_id, "NASTRAN CTETRA id")
        : format_nastran_small_field_integer(element_id, "NASTRAN CTETRA id"),
      options.field_format == NastranFieldFormat::large_fixed
        ? format_nastran_large_field_integer(property_id, "NASTRAN CTETRA PID")
        : format_nastran_small_field_integer(property_id, "NASTRAN CTETRA PID"),
      options.field_format == NastranFieldFormat::large_fixed
        ? format_nastran_large_field_integer(
            node_indices.at(pack_entity_ref(nodes[0])),
            "NASTRAN CTETRA G1"
          )
        : format_nastran_small_field_integer(
            node_indices.at(pack_entity_ref(nodes[0])),
            "NASTRAN CTETRA G1"
          ),
      options.field_format == NastranFieldFormat::large_fixed
        ? format_nastran_large_field_integer(
            node_indices.at(pack_entity_ref(nodes[1])),
            "NASTRAN CTETRA G2"
          )
        : format_nastran_small_field_integer(
            node_indices.at(pack_entity_ref(nodes[1])),
            "NASTRAN CTETRA G2"
          ),
      options.field_format == NastranFieldFormat::large_fixed
        ? format_nastran_large_field_integer(
            node_indices.at(pack_entity_ref(nodes[2])),
            "NASTRAN CTETRA G3"
          )
        : format_nastran_small_field_integer(
            node_indices.at(pack_entity_ref(nodes[2])),
            "NASTRAN CTETRA G3"
          ),
      options.field_format == NastranFieldFormat::large_fixed
        ? format_nastran_large_field_integer(
            node_indices.at(pack_entity_ref(nodes[3])),
            "NASTRAN CTETRA G4"
          )
        : format_nastran_small_field_integer(
            node_indices.at(pack_entity_ref(nodes[3])),
            "NASTRAN CTETRA G4"
          ),
    };
    if(options.field_format == NastranFieldFormat::large_fixed) {
      write_nastran_large_field_record(output, "CTETRA", fields);
      return;
    }

    write_nastran_small_field_record(output, "CTETRA", fields);
  };

  struct NastranMaterialExportRecord final {
    NastranMaterialCard card = NastranMaterialCard::unspecified;
    std::uint32_t material_id = invalid_index;
    double youngs_modulus = std::numeric_limits<double>::quiet_NaN();
    double shear_modulus = std::numeric_limits<double>::quiet_NaN();
    double poisson_ratio = std::numeric_limits<double>::quiet_NaN();
  };

  struct NastranPropertyExportRecord final {
    NastranPropertyCard card = NastranPropertyCard::unspecified;
    std::uint32_t property_id = invalid_index;
    std::uint32_t material_id = invalid_index;
    double rod_area = std::numeric_limits<double>::quiet_NaN();
    double shell_thickness = std::numeric_limits<double>::quiet_NaN();
  };

  std::map<std::uint32_t, NastranMaterialExportRecord> materials_by_id;
  std::map<std::uint32_t, NastranPropertyExportRecord> properties_by_id;

  const auto merge_material_metadata = [&](const EntityGroup &entity_group) {
    if(entity_group.import_info().format != EntityGroupImportFormat::nastran) {
      return;
    }

    const auto &import_info = entity_group.import_info().nastran;
    if(import_info.material_card != NastranMaterialCard::mat1) {
      return;
    }
    if(import_info.material_id == invalid_index) {
      fail(
        base::StatusCode::unsupported,
        "NASTRAN export requires a valid MID when MAT1 metadata is attached to a entity_group."
      );
    }

    NastranMaterialExportRecord candidate;
    candidate.card = import_info.material_card;
    candidate.material_id = import_info.material_id;
    candidate.youngs_modulus = import_info.youngs_modulus;
    candidate.shear_modulus = import_info.shear_modulus;
    candidate.poisson_ratio = import_info.poisson_ratio;

    const auto [it, inserted] = materials_by_id.emplace(candidate.material_id, candidate);
    if(inserted) {
      return;
    }

    const auto &current = it->second;
    if(current.card != candidate.card ||
       !same_nastran_optional_real(current.youngs_modulus, candidate.youngs_modulus) ||
       !same_nastran_optional_real(current.shear_modulus, candidate.shear_modulus) ||
       !same_nastran_optional_real(current.poisson_ratio, candidate.poisson_ratio)) {
      fail(
        base::StatusCode::unsupported,
        "NASTRAN export found conflicting MAT1 metadata for the same MID."
      );
    }
  };

  const auto merge_property_metadata = [&](const EntityGroup &entity_group) {
    if(entity_group.import_info().format != EntityGroupImportFormat::nastran) {
      return;
    }

    const auto &import_info = entity_group.import_info().nastran;
    if(import_info.property_card == NastranPropertyCard::unspecified) {
      return;
    }
    if(entity_group.zone_id() == invalid_index) {
      fail(
        base::StatusCode::unsupported,
        "NASTRAN export requires a valid PID when imported property metadata is attached to a entity_group."
      );
    }

    NastranPropertyExportRecord candidate;
    candidate.card = import_info.property_card;
    candidate.property_id = entity_group.zone_id();
    candidate.material_id = import_info.material_id;
    candidate.rod_area = import_info.rod_area;
    candidate.shell_thickness = import_info.shell_thickness;

    const auto [it, inserted] = properties_by_id.emplace(candidate.property_id, candidate);
    if(inserted) {
      return;
    }

    const auto &current = it->second;
    if(current.card != candidate.card || current.material_id != candidate.material_id ||
       !same_nastran_optional_real(current.rod_area, candidate.rod_area) ||
       !same_nastran_optional_real(current.shell_thickness, candidate.shell_thickness)) {
      fail(
        base::StatusCode::unsupported,
        "NASTRAN export found conflicting property metadata for the same PID."
      );
    }
  };

  const auto write_mat1 = [&](const NastranMaterialExportRecord &material) {
    if(options.field_format == NastranFieldFormat::free_field) {
      output << "MAT1," << material.material_id;
      if(std::isfinite(material.youngs_modulus) || std::isfinite(material.shear_modulus) ||
         std::isfinite(material.poisson_ratio)) {
        output << ',';
        if(std::isfinite(material.youngs_modulus)) {
          output << std::setprecision(17) << material.youngs_modulus;
        }
        output << ',';
        if(std::isfinite(material.shear_modulus)) {
          output << std::setprecision(17) << material.shear_modulus;
        }
        output << ',';
        if(std::isfinite(material.poisson_ratio)) {
          output << std::setprecision(17) << material.poisson_ratio;
        }
      }
      output << '\n';
      return;
    }

    const std::vector<std::string> fields = {
      options.field_format == NastranFieldFormat::large_fixed
        ? format_nastran_large_field_integer(material.material_id, "NASTRAN MAT1 MID")
        : format_nastran_small_field_integer(material.material_id, "NASTRAN MAT1 MID"),
      std::isfinite(material.youngs_modulus)
        ? (options.field_format == NastranFieldFormat::large_fixed
             ? format_nastran_large_field_real(material.youngs_modulus, "NASTRAN MAT1 E")
             : format_nastran_small_field_real(material.youngs_modulus, "NASTRAN MAT1 E"))
        : std::string {},
      std::isfinite(material.shear_modulus)
        ? (options.field_format == NastranFieldFormat::large_fixed
             ? format_nastran_large_field_real(material.shear_modulus, "NASTRAN MAT1 G")
             : format_nastran_small_field_real(material.shear_modulus, "NASTRAN MAT1 G"))
        : std::string {},
      std::isfinite(material.poisson_ratio)
        ? (options.field_format == NastranFieldFormat::large_fixed
             ? format_nastran_large_field_real(material.poisson_ratio, "NASTRAN MAT1 NU")
             : format_nastran_small_field_real(material.poisson_ratio, "NASTRAN MAT1 NU"))
        : std::string {},
    };
    if(options.field_format == NastranFieldFormat::large_fixed) {
      write_nastran_large_field_record(output, "MAT1", fields);
      return;
    }
    write_nastran_small_field_record(output, "MAT1", fields);
  };

  const auto write_property = [&](const NastranPropertyExportRecord &property) {
    const auto material_field =
      property.material_id == invalid_index
        ? std::string {}
        : (options.field_format == NastranFieldFormat::large_fixed
             ? format_nastran_large_field_integer(property.material_id, "NASTRAN property MID")
             : format_nastran_small_field_integer(property.material_id, "NASTRAN property MID"));

    if(property.card == NastranPropertyCard::prod) {
      if(options.field_format == NastranFieldFormat::free_field) {
        output << "PROD," << property.property_id << ',';
        if(property.material_id != invalid_index) {
          output << property.material_id;
        }
        if(std::isfinite(property.rod_area)) {
          output << ',' << std::setprecision(17) << property.rod_area;
        }
        output << '\n';
        return;
      }

      const std::vector<std::string> fields = {
        options.field_format == NastranFieldFormat::large_fixed
          ? format_nastran_large_field_integer(property.property_id, "NASTRAN PROD PID")
          : format_nastran_small_field_integer(property.property_id, "NASTRAN PROD PID"),
        material_field,
        std::isfinite(property.rod_area)
          ? (options.field_format == NastranFieldFormat::large_fixed
               ? format_nastran_large_field_real(property.rod_area, "NASTRAN PROD A")
               : format_nastran_small_field_real(property.rod_area, "NASTRAN PROD A"))
          : std::string {},
      };
      if(options.field_format == NastranFieldFormat::large_fixed) {
        write_nastran_large_field_record(output, "PROD", fields);
        return;
      }
      write_nastran_small_field_record(output, "PROD", fields);
      return;
    }

    if(property.card == NastranPropertyCard::pshell) {
      if(options.field_format == NastranFieldFormat::free_field) {
        output << "PSHELL," << property.property_id << ',';
        if(property.material_id != invalid_index) {
          output << property.material_id;
        }
        if(std::isfinite(property.shell_thickness)) {
          output << ',' << std::setprecision(17) << property.shell_thickness;
        }
        output << '\n';
        return;
      }

      const std::vector<std::string> fields = {
        options.field_format == NastranFieldFormat::large_fixed
          ? format_nastran_large_field_integer(property.property_id, "NASTRAN PSHELL PID")
          : format_nastran_small_field_integer(property.property_id, "NASTRAN PSHELL PID"),
        material_field,
        std::isfinite(property.shell_thickness)
          ? (options.field_format == NastranFieldFormat::large_fixed
               ? format_nastran_large_field_real(property.shell_thickness, "NASTRAN PSHELL T")
               : format_nastran_small_field_real(property.shell_thickness, "NASTRAN PSHELL T"))
          : std::string {},
      };
      if(options.field_format == NastranFieldFormat::large_fixed) {
        write_nastran_large_field_record(output, "PSHELL", fields);
        return;
      }
      write_nastran_small_field_record(output, "PSHELL", fields);
      return;
    }

    if(property.card == NastranPropertyCard::psolid) {
      if(options.field_format == NastranFieldFormat::free_field) {
        output << "PSOLID," << property.property_id << ',';
        if(property.material_id != invalid_index) {
          output << property.material_id;
        }
        output << '\n';
        return;
      }

      const std::vector<std::string> fields = {
        options.field_format == NastranFieldFormat::large_fixed
          ? format_nastran_large_field_integer(property.property_id, "NASTRAN PSOLID PID")
          : format_nastran_small_field_integer(property.property_id, "NASTRAN PSOLID PID"),
        material_field,
      };
      if(options.field_format == NastranFieldFormat::large_fixed) {
        write_nastran_large_field_record(output, "PSOLID", fields);
        return;
      }
      write_nastran_small_field_record(output, "PSOLID", fields);
    }
  };

  for(const auto &entity_group : domain.entity_groups()) {
    merge_material_metadata(entity_group);
    merge_property_metadata(entity_group);
  }

  for(const auto &[material_id, material] : materials_by_id) {
    static_cast<void>(material_id);
    write_mat1(material);
  }
  for(const auto &[property_id, property] : properties_by_id) {
    static_cast<void>(property_id);
    write_property(property);
  }

  for(const auto &entity_group : domain.entity_groups()) {
    if(entity_group.order() != EntityOrder::node) {
      continue;
    }

    for(const auto &node : entity_group.nodes()) {
      write_grid(
        node_indices.at(pack_entity_ref({node.header.entity_group, node.header.index})),
        node
      );
    }
  }

  std::uint64_t next_element = 1U;
  for(const auto &entity_group : domain.entity_groups()) {
    if(entity_group.order() != EntityOrder::edge) {
      continue;
    }
    for(const auto &edge : entity_group.edges()) {
      write_edge(next_element++, entity_group.zone_id(), {edge.header.entity_group, edge.header.index});
    }
  }

  for(const auto &entity_group : domain.entity_groups()) {
    if(entity_group.order() != EntityOrder::face) {
      continue;
    }
    for(std::size_t face_index = 0U; face_index < entity_group.faces().size();) {
      const auto face_ref =
        EntityRef {entity_group.id(), static_cast<std::uint32_t>(face_index)};
      if(is_nastran_cquad4_entity_group(entity_group) && face_index + 1U < entity_group.faces().size()) {
        std::array<EntityRef, 4> quad_nodes {};
        const auto second_face_ref =
          EntityRef {entity_group.id(), static_cast<std::uint32_t>(face_index + 1U)};
        if(try_build_nastran_cquad4(domain, face_ref, second_face_ref, quad_nodes)) {
          write_quad(next_element++, entity_group.zone_id(), quad_nodes);
          face_index += 2U;
          continue;
        }
      }

      write_face(next_element++, entity_group.zone_id(), face_ref);
      ++face_index;
    }
  }

  for(const auto &entity_group : domain.entity_groups()) {
    if(entity_group.order() != EntityOrder::cell) {
      continue;
    }
    for(const auto &cell : entity_group.cells()) {
      write_cell(next_element++, entity_group.zone_id(), {cell.header.entity_group, cell.header.index});
    }
  }

  if(options.write_enddata) {
    output << "ENDDATA\n";
  }
}

[[nodiscard]] std::string format_nastran_small_field_integer(
  std::uint64_t value,
  std::string_view context
)
{
  auto token = std::to_string(value);
  if(token.size() > 8U) {
    fail(
      base::StatusCode::unsupported,
      std::string(context) + " exceeds the 8-character small-field width."
    );
  }
  return token;
}

[[nodiscard]] std::string format_nastran_large_field_integer(
  std::uint64_t value,
  std::string_view context
)
{
  auto token = std::to_string(value);
  if(token.size() > 16U) {
    fail(
      base::StatusCode::unsupported,
      std::string(context) + " exceeds the 16-character large-field width."
    );
  }
  return token;
}

[[nodiscard]] std::string normalize_nastran_scientific(std::string token)
{
  auto exponent = token.find_first_of("eE");
  if(exponent == std::string::npos) {
    return token;
  }

  token[exponent] = 'E';
  while(exponent > 0U && token[exponent - 1U] == '0' &&
        token.find('.', 0U) < exponent) {
    token.erase(exponent - 1U, 1U);
    --exponent;
  }
  if(exponent > 0U && token[exponent - 1U] == '.') {
    token.erase(exponent - 1U, 1U);
    --exponent;
  }

  auto digit = exponent + 2U;
  while(digit + 1U < token.size() && token[digit] == '0') {
    token.erase(digit, 1U);
  }
  return token;
}

[[nodiscard]] std::string format_nastran_small_field_real(
  double value,
  std::string_view context
)
{
  if(value == 0.0) {
    return "0.";
  }

  auto try_format = [&](bool scientific, int precision) {
    std::ostringstream stream;
    if(scientific) {
      stream << std::uppercase << std::scientific;
    }
    stream << std::setprecision(precision) << value;
    auto token = scientific ? normalize_nastran_scientific(stream.str()) : stream.str();
    return token;
  };

  for(int precision = 7; precision >= 1; --precision) {
    const auto token = try_format(false, precision);
    if(token.size() <= 8U) {
      return token;
    }
  }

  for(int precision = 6; precision >= 0; --precision) {
    const auto token = try_format(true, precision);
    if(token.size() <= 8U) {
      return token;
    }
  }

  fail(
    base::StatusCode::unsupported,
    std::string(context) + " cannot be represented in a single 8-character small field."
  );
}

[[nodiscard]] std::string format_nastran_large_field_real(
  double value,
  std::string_view context
)
{
  if(value == 0.0) {
    return "0.";
  }

  auto try_format = [&](bool scientific, int precision) {
    std::ostringstream stream;
    if(scientific) {
      stream << std::uppercase << std::scientific;
    }
    stream << std::setprecision(precision) << value;
    return scientific ? normalize_nastran_scientific(stream.str()) : stream.str();
  };

  for(int precision = 15; precision >= 1; --precision) {
    const auto token = try_format(false, precision);
    if(token.size() <= 16U) {
      return token;
    }
  }

  for(int precision = 14; precision >= 0; --precision) {
    const auto token = try_format(true, precision);
    if(token.size() <= 16U) {
      return token;
    }
  }

  fail(
    base::StatusCode::unsupported,
    std::string(context) + " cannot be represented in a single 16-character large field."
  );
}

void write_nastran_small_field_record(
  std::ostream &output,
  std::string_view card,
  const std::vector<std::string> &fields
)
{
  output << std::left << std::setw(8) << card << std::right;
  for(const auto &field : fields) {
    if(field.size() > 8U) {
      fail(
        base::StatusCode::unsupported,
        "NASTRAN small-field export produced a field wider than 8 characters."
      );
    }
    output << std::setw(8) << field;
  }
  output << '\n';
}

void write_nastran_large_field_record(
  std::ostream &output,
  std::string_view card,
  const std::vector<std::string> &fields
)
{
  const auto first_card = std::string(card) + '*';
  if(first_card.size() > 8U) {
    fail(
      base::StatusCode::unsupported,
      "NASTRAN large-field export produced a card name wider than 8 characters."
    );
  }

  const auto write_line = [&](std::string_view line_card, std::size_t field_offset) {
    output << std::left << std::setw(8) << line_card;
    for(std::size_t index = 0U; index < 4U && field_offset + index < fields.size(); ++index) {
      const auto &field = fields[field_offset + index];
      if(field.size() > 16U) {
        fail(
          base::StatusCode::unsupported,
          "NASTRAN large-field export produced a field wider than 16 characters."
        );
      }
      output << std::setw(16) << field;
    }
    output << '\n';
  };

  write_line(first_card, 0U);
  for(std::size_t offset = 4U; offset < fields.size(); offset += 4U) {
    write_line("*", offset);
  }
}

void validate_domain_for_obj_export(const Domain &domain)
{
  if(domain.cell_count() != 0U) {
    fail(
      base::StatusCode::unsupported,
      "OBJ export currently supports only surface meshes without volume cells."
    );
  }

  for(const auto &entity_group : domain.entity_groups()) {
    switch(entity_group.order()) {
    case EntityOrder::node:
      for(const auto &node : entity_group.nodes()) {
        if(node.header.kind != EntityKind::node_point) {
          fail(
            base::StatusCode::unsupported,
            "OBJ export currently supports only point nodes."
          );
        }
      }
      break;
    case EntityOrder::edge:
      for(const auto &edge : entity_group.edges()) {
        if(edge.header.kind != EntityKind::edge_line) {
          fail(
            base::StatusCode::unsupported,
            "OBJ export currently supports only line edges."
          );
        }
      }
      break;
    case EntityOrder::face:
      for(const auto &face : entity_group.faces()) {
        if(face.header.kind != EntityKind::face_triangle) {
          fail(
            base::StatusCode::unsupported,
            "OBJ export currently supports only triangle faces."
          );
        }
      }
      break;
    case EntityOrder::cell:
      break;
    default:
      fail(
        base::StatusCode::unsupported,
        "OBJ export found an unsupported entity order."
      );
    }
  }
}

void write_obj_surface(
  const Domain &domain,
  std::string_view path,
  const ObjExportOptions &options
)
{
  validate_domain_for_obj_export(domain);

  const auto node_indices = build_node_indices(domain, options.export_role);
  const std::string file_path(path);
  std::ofstream output(file_path, std::ios::trunc);
  if(!output) {
    fail(base::StatusCode::io_error, "OBJ export could not open the output file.");
  }

  output << "# SQMesh OBJ export\n";
  if(options.write_object_name && !domain.name().empty()) {
    output << "o " << sanitize_obj_label(domain.name()) << '\n';
  }

  output << std::setprecision(17);
  for(const auto &entity_group : domain.entity_groups()) {
    if(entity_group.order() != EntityOrder::node) {
      continue;
    }
    if(entity_group.role() != options.export_role) {
      continue;
    }

    for(const auto &node : entity_group.nodes()) {
      output << "v " << node.coordinates[0] << ' ' << node.coordinates[1]
             << ' ' << node.coordinates[2] << '\n';
    }
  }

  for(const auto &entity_group : domain.entity_groups()) {
    if(entity_group.order() != EntityOrder::edge) {
      continue;
    }
    if(entity_group.role() != options.export_role) {
      continue;
    }
    if(options.write_groups) {
      output << "g " << sanitize_obj_label(entity_group.name()) << '\n';
    }
    for(const auto &edge : entity_group.edges()) {
      const auto nodes = domain.edge_nodes({edge.header.entity_group, edge.header.index});
      output << "l";
      for(const auto node_ref : nodes) {
        output << ' ' << node_indices.at(pack_entity_ref(node_ref));
      }
      output << '\n';
    }
  }

  for(const auto &entity_group : domain.entity_groups()) {
    if(entity_group.order() != EntityOrder::face) {
      continue;
    }
    if(entity_group.role() != options.export_role) {
      continue;
    }
    if(options.write_groups) {
      output << "g " << sanitize_obj_label(entity_group.name()) << '\n';
    }
    for(const auto &face : entity_group.faces()) {
      const auto nodes = domain.face_nodes({face.header.entity_group, face.header.index});
      output << "f";
      for(const auto node_ref : nodes) {
        output << ' ' << node_indices.at(pack_entity_ref(node_ref));
      }
      output << '\n';
    }
  }
}

} // namespace

base::StatusCode import_msh(
  std::string_view path,
  MeshHandle &mesh_handle,
  const MshImportOptions &options,
  base::ContextHandle context_handle
) noexcept
{
  mesh_handle = sqmesh::invalid_handle;

  try {
    const auto description = parse_msh_description(path, options);
    auto domain = build_domain_from_description(description);
    const auto resolved_context = require_context(context_handle);
    return store_imported_domain(std::move(domain), "MSH import", mesh_handle, resolved_context);
  }
  catch(const MeshIoFailure &failure) {
    return core::detail::publish_error(failure.code(), failure.what());
  }
  catch(const std::exception &exception) {
    return core::detail::publish_error(
      base::StatusCode::internal_error,
      exception.what()
    );
  }
  catch(...) {
    return core::detail::publish_error(
      base::StatusCode::internal_error,
      "MSH import failed with an unknown error."
    );
  }
}

base::StatusCode export_msh(
  MeshHandle mesh_handle,
  std::string_view path,
  const MshExportOptions &options,
  base::ContextHandle context_handle
) noexcept
{
  try {
    validate_non_empty_path(path, "MSH export requires a non-empty path.");
    return core::detail::with_mesh_domain(
      mesh_handle,
      [&](const Domain &domain) {
        try {
          switch(options.format_version) {
          case MshFormatVersion::gmsh22_ascii:
            write_msh22_ascii(domain, path, options);
            break;
          case MshFormatVersion::gmsh22_binary:
            write_msh22_binary(domain, path, options);
            break;
          case MshFormatVersion::gmsh41_ascii:
            write_msh41_ascii(domain, path, options);
            break;
          case MshFormatVersion::gmsh41_binary:
            write_msh41_binary(domain, path, options);
            break;
          default:
            fail(
              base::StatusCode::invalid_argument,
              "MSH export received an unsupported format-version selector."
            );
          }

          return base::StatusCode::ok;
        }
        catch(const MeshIoFailure &failure) {
          return core::detail::publish_error(failure.code(), failure.what());
        }
        catch(const std::exception &exception) {
          return core::detail::publish_error(
            base::StatusCode::internal_error,
            exception.what()
          );
        }
        catch(...) {
          return core::detail::publish_error(
            base::StatusCode::internal_error,
            "MSH export failed with an unknown error."
          );
        }
      },
      context_handle
    );
  }
  catch(const MeshIoFailure &failure) {
    return core::detail::publish_error(failure.code(), failure.what());
  }
  catch(const std::exception &exception) {
    return core::detail::publish_error(
      base::StatusCode::internal_error,
      exception.what()
    );
  }
  catch(...) {
    return core::detail::publish_error(
      base::StatusCode::internal_error,
      "MSH export failed with an unknown error."
    );
  }
}

base::StatusCode import_obj(
  std::string_view path,
  MeshHandle &mesh_handle,
  const ObjImportOptions &options,
  base::ContextHandle context_handle
) noexcept
{
  mesh_handle = sqmesh::invalid_handle;

  try {
    const auto description = parse_obj_description(path, options);
    auto domain = build_domain_from_description(description);
    const auto resolved_context = require_context(context_handle);
    return store_imported_domain(std::move(domain), "OBJ import", mesh_handle, resolved_context);
  }
  catch(const MeshIoFailure &failure) {
    return core::detail::publish_error(failure.code(), failure.what());
  }
  catch(const std::exception &exception) {
    return core::detail::publish_error(
      base::StatusCode::internal_error,
      exception.what()
    );
  }
  catch(...) {
    return core::detail::publish_error(
      base::StatusCode::internal_error,
      "OBJ import failed with an unknown error."
    );
  }
}

namespace testing {

base::StatusCode inspect_obj_import_review(
  std::string_view path,
  ObjImportReview &review,
  const ObjImportOptions &options
) noexcept
{
  review = {};

  try {
    review = make_obj_import_review(parse_obj_description(path, options));
    return base::StatusCode::ok;
  }
  catch(const MeshIoFailure &failure) {
    return core::detail::publish_error(failure.code(), failure.what());
  }
  catch(const std::exception &exception) {
    return core::detail::publish_error(
      base::StatusCode::internal_error,
      exception.what()
    );
  }
  catch(...) {
    return core::detail::publish_error(
      base::StatusCode::internal_error,
      "OBJ import review failed with an unknown error."
    );
  }
}

}  // namespace testing

base::StatusCode export_obj(
  MeshHandle mesh_handle,
  std::string_view path,
  const ObjExportOptions &options,
  base::ContextHandle context_handle
) noexcept
{
  try {
    validate_non_empty_path(path, "OBJ export requires a non-empty path.");
    return core::detail::with_mesh_domain(
      mesh_handle,
      [&](const Domain &domain) {
        try {
          write_obj_surface(domain, path, options);
          return base::StatusCode::ok;
        }
        catch(const MeshIoFailure &failure) {
          return core::detail::publish_error(failure.code(), failure.what());
        }
        catch(const std::exception &exception) {
          return core::detail::publish_error(
            base::StatusCode::internal_error,
            exception.what()
          );
        }
        catch(...) {
          return core::detail::publish_error(
            base::StatusCode::internal_error,
            "OBJ export failed with an unknown error."
          );
        }
      },
      context_handle
    );
  }
  catch(const MeshIoFailure &failure) {
    return core::detail::publish_error(failure.code(), failure.what());
  }
  catch(const std::exception &exception) {
    return core::detail::publish_error(
      base::StatusCode::internal_error,
      exception.what()
    );
  }
  catch(...) {
    return core::detail::publish_error(
      base::StatusCode::internal_error,
      "OBJ export failed with an unknown error."
    );
  }
}

base::StatusCode import_nastran(
  std::string_view path,
  MeshHandle &mesh_handle,
  const NastranImportOptions &options,
  base::ContextHandle context_handle
) noexcept
{
  mesh_handle = sqmesh::invalid_handle;

  try {
    const auto description = parse_nastran_description(path, options);
    auto domain = build_domain_from_description(description);
    const auto resolved_context = require_context(context_handle);
    return store_imported_domain(
      std::move(domain),
      "NASTRAN import",
      mesh_handle,
      resolved_context
    );
  }
  catch(const MeshIoFailure &failure) {
    return core::detail::publish_error(failure.code(), failure.what());
  }
  catch(const std::exception &exception) {
    return core::detail::publish_error(
      base::StatusCode::internal_error,
      exception.what()
    );
  }
  catch(...) {
    return core::detail::publish_error(
      base::StatusCode::internal_error,
      "NASTRAN import failed with an unknown error."
    );
  }
}

base::StatusCode export_nastran(
  MeshHandle mesh_handle,
  std::string_view path,
  const NastranExportOptions &options,
  base::ContextHandle context_handle
) noexcept
{
  try {
    validate_non_empty_path(path, "NASTRAN export requires a non-empty path.");
    return core::detail::with_mesh_domain(
      mesh_handle,
      [&](const Domain &domain) {
        try {
          write_nastran_bulk_data(domain, path, options);
          return base::StatusCode::ok;
        }
        catch(const MeshIoFailure &failure) {
          return core::detail::publish_error(failure.code(), failure.what());
        }
        catch(const std::exception &exception) {
          return core::detail::publish_error(
            base::StatusCode::internal_error,
            exception.what()
          );
        }
        catch(...) {
          return core::detail::publish_error(
            base::StatusCode::internal_error,
            "NASTRAN export failed with an unknown error."
          );
        }
      },
      context_handle
    );
  }
  catch(const MeshIoFailure &failure) {
    return core::detail::publish_error(failure.code(), failure.what());
  }
  catch(const std::exception &exception) {
    return core::detail::publish_error(
      base::StatusCode::internal_error,
      exception.what()
    );
  }
  catch(...) {
    return core::detail::publish_error(
      base::StatusCode::internal_error,
      "NASTRAN export failed with an unknown error."
    );
  }
}

bool cgns_io_available() noexcept
{
#if defined(sqmesh_HAVE_CGNS)
  return true;
#else
  return false;
#endif
}

base::StatusCode import_cgns(
  std::string_view path,
  MeshHandle &mesh_handle,
  const CgnsImportOptions &options,
  base::ContextHandle context_handle
) noexcept
{
  mesh_handle = sqmesh::invalid_handle;

#if defined(sqmesh_HAVE_CGNS)
  try {
    const auto description = parse_cgns_description(path, options);
    auto domain = build_domain_from_description(description);
    const auto resolved_context = require_context(context_handle);
    return store_imported_domain(std::move(domain), "CGNS import", mesh_handle, resolved_context);
  }
  catch(const MeshIoFailure &failure) {
    return core::detail::publish_error(failure.code(), failure.what());
  }
  catch(const std::exception &exception) {
    return core::detail::publish_error(
      base::StatusCode::internal_error,
      exception.what()
    );
  }
  catch(...) {
    return core::detail::publish_error(
      base::StatusCode::internal_error,
      "CGNS import failed with an unknown error."
    );
  }
#else
  static_cast<void>(path);
  static_cast<void>(options);
  static_cast<void>(context_handle);
  return core::detail::publish_error(
    base::StatusCode::unsupported,
    "CGNS import requires sqmesh_ENABLE_CGNS and a discoverable CGNS library."
  );
#endif
}

base::StatusCode export_cgns(
  MeshHandle mesh_handle,
  std::string_view path,
  const CgnsExportOptions &options,
  base::ContextHandle context_handle
) noexcept
{
#if defined(sqmesh_HAVE_CGNS)
  try {
    validate_non_empty_path(path, "CGNS export requires a non-empty path.");
    return core::detail::with_mesh_domain(
      mesh_handle,
      [&](const Domain &domain) {
        try {
          if(options.write_multi_zone) {
            write_cgns_multi_zone_unstructured(domain, path, options);
          }
          else {
            write_cgns_unstructured(domain, path, options);
          }

          return base::StatusCode::ok;
        }
        catch(const MeshIoFailure &failure) {
          return core::detail::publish_error(failure.code(), failure.what());
        }
        catch(const std::exception &exception) {
          return core::detail::publish_error(
            base::StatusCode::internal_error,
            exception.what()
          );
        }
        catch(...) {
          return core::detail::publish_error(
            base::StatusCode::internal_error,
            "CGNS export failed with an unknown error."
          );
        }
      },
      context_handle
    );
  }
  catch(const MeshIoFailure &failure) {
    return core::detail::publish_error(failure.code(), failure.what());
  }
  catch(const std::exception &exception) {
    return core::detail::publish_error(
      base::StatusCode::internal_error,
      exception.what()
    );
  }
  catch(...) {
    return core::detail::publish_error(
      base::StatusCode::internal_error,
      "CGNS export failed with an unknown error."
    );
  }
#else
  static_cast<void>(mesh_handle);
  static_cast<void>(path);
  static_cast<void>(options);
  static_cast<void>(context_handle);
  return core::detail::publish_error(
    base::StatusCode::unsupported,
    "CGNS export requires sqmesh_ENABLE_CGNS and a discoverable CGNS library."
  );
#endif
}

} // namespace sqmesh::mesh
