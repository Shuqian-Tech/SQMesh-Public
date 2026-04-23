// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "stl_import.hpp"

#include "../core/runtime_registry.hpp"

#include "geometry_model_storage.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <memory>
#include <queue>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace sqmesh::model::detail {
namespace {

constexpr double kMinimumMergeTolerance = 1.0e-12;
constexpr double kPi = 3.14159265358979323846264338327950288;

struct RawFacet final {
  geo::Vector3 file_normal {0.0, 0.0, 0.0};
  std::array<geo::Point3, 3> vertices {};
};

struct BoundingBox final {
  geo::Point3 min {
    std::numeric_limits<double>::max(),
    std::numeric_limits<double>::max(),
    std::numeric_limits<double>::max(),
  };
  geo::Point3 max {
    std::numeric_limits<double>::lowest(),
    std::numeric_limits<double>::lowest(),
    std::numeric_limits<double>::lowest(),
  };
  bool initialized = false;
};

struct TriangleData final {
  std::array<std::uint32_t, 3> vertices {};
  std::array<std::uint32_t, 3> edges {};
  geo::Vector3 normal {0.0, 0.0, 0.0};
};

struct VertexData final {
  geo::Point3 position {0.0, 0.0, 0.0};
  std::vector<std::uint32_t> parent_edges {};
};

struct EdgeData final {
  std::array<std::uint32_t, 2> vertices {};
  std::vector<std::uint32_t> triangles {};
  std::vector<std::uint32_t> parent_faces {};
};

struct FaceData final {
  std::vector<std::uint32_t> triangles {};
  std::vector<std::uint32_t> child_edges {};
  std::vector<geo::FaceBoundaryLoop> boundary_loops {};
};

struct DiscreteModel final {
  std::vector<VertexData> vertices {};
  std::vector<EdgeData> edges {};
  std::vector<TriangleData> triangles {};
  std::vector<FaceData> faces {};
  geo::ModelSummary summary {};
};

struct CellKey final {
  long long x = 0;
  long long y = 0;
  long long z = 0;
};

struct CellKeyHash final {
  [[nodiscard]] std::size_t operator()(const CellKey &key) const noexcept
  {
    const auto hx = std::hash<long long> {}(key.x);
    const auto hy = std::hash<long long> {}(key.y);
    const auto hz = std::hash<long long> {}(key.z);
    return hx ^ (hy << 1U) ^ (hz << 2U);
  }
};

[[nodiscard]] bool operator==(const CellKey &lhs, const CellKey &rhs) noexcept
{
  return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

struct EdgeKey final {
  std::uint32_t v0 = 0U;
  std::uint32_t v1 = 0U;
};

struct EdgeKeyHash final {
  [[nodiscard]] std::size_t operator()(const EdgeKey &key) const noexcept
  {
    const auto combined =
      (static_cast<std::uint64_t>(key.v0) << 32U) |
      static_cast<std::uint64_t>(key.v1);
    return std::hash<std::uint64_t> {}(combined);
  }
};

[[nodiscard]] bool operator==(const EdgeKey &lhs, const EdgeKey &rhs) noexcept
{
  return lhs.v0 == rhs.v0 && lhs.v1 == rhs.v1;
}

struct TriangleKey final {
  std::array<std::uint32_t, 3> vertices {};
};

struct TriangleKeyHash final {
  [[nodiscard]] std::size_t operator()(const TriangleKey &key) const noexcept
  {
    const auto h0 = std::hash<std::uint32_t> {}(key.vertices[0]);
    const auto h1 = std::hash<std::uint32_t> {}(key.vertices[1]);
    const auto h2 = std::hash<std::uint32_t> {}(key.vertices[2]);
    return h0 ^ (h1 << 1U) ^ (h2 << 2U);
  }
};

[[nodiscard]] bool operator==(const TriangleKey &lhs, const TriangleKey &rhs) noexcept
{
  return lhs.vertices == rhs.vertices;
}

[[nodiscard]] geo::TopologyEntityId make_entity_id(
  geo::TopologyDimension dimension,
  std::size_t index
) noexcept
{
  return geo::TopologyEntityId {
    dimension,
    static_cast<std::uint32_t>(index),
  };
}

[[nodiscard]] geo::Point3 subtract(
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

[[nodiscard]] geo::Point3 interpolate(
  const geo::Point3 &lhs,
  const geo::Point3 &rhs,
  double t
) noexcept
{
  return {
    lhs[0] + (rhs[0] - lhs[0]) * t,
    lhs[1] + (rhs[1] - lhs[1]) * t,
    lhs[2] + (rhs[2] - lhs[2]) * t,
  };
}

[[nodiscard]] geo::Vector3 cross_product(
  const geo::Point3 &lhs,
  const geo::Point3 &rhs
) noexcept
{
  return {
    lhs[1] * rhs[2] - lhs[2] * rhs[1],
    lhs[2] * rhs[0] - lhs[0] * rhs[2],
    lhs[0] * rhs[1] - lhs[1] * rhs[0],
  };
}

[[nodiscard]] double dot_product(
  const geo::Vector3 &lhs,
  const geo::Vector3 &rhs
) noexcept
{
  return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2];
}

[[nodiscard]] double squared_norm(const geo::Vector3 &vector) noexcept
{
  return dot_product(vector, vector);
}

[[nodiscard]] double norm(const geo::Vector3 &vector) noexcept
{
  return std::sqrt(squared_norm(vector));
}

[[nodiscard]] geo::Vector3 normalized(const geo::Vector3 &vector) noexcept
{
  const double length = norm(vector);
  if(length <= 0.0) {
    return {0.0, 0.0, 0.0};
  }

  return {
    vector[0] / length,
    vector[1] / length,
    vector[2] / length,
  };
}

[[nodiscard]] double point_distance(
  const geo::Point3 &lhs,
  const geo::Point3 &rhs
) noexcept
{
  return norm(subtract(lhs, rhs));
}

[[nodiscard]] double point_distance_squared(
  const geo::Point3 &lhs,
  const geo::Point3 &rhs
) noexcept
{
  return squared_norm(subtract(lhs, rhs));
}

[[nodiscard]] std::string trim(std::string_view input)
{
  std::size_t begin = 0U;
  while(begin < input.size() && std::isspace(static_cast<unsigned char>(input[begin])) != 0) {
    ++begin;
  }

  std::size_t end = input.size();
  while(end > begin &&
        std::isspace(static_cast<unsigned char>(input[end - 1U])) != 0) {
    --end;
  }

  return std::string(input.substr(begin, end - begin));
}

[[nodiscard]] std::string lower_ascii(std::string_view input)
{
  std::string lowered;
  lowered.reserve(input.size());
  for(const char value : input) {
    lowered.push_back(static_cast<char>(
      std::tolower(static_cast<unsigned char>(value))
    ));
  }
  return lowered;
}

[[nodiscard]] std::vector<std::string> tokenize(std::string_view line)
{
  std::vector<std::string> tokens;
  std::istringstream stream {std::string(line)};
  std::string token;
  while(stream >> token) {
    tokens.push_back(std::move(token));
  }
  return tokens;
}

[[nodiscard]] bool read_non_empty_line(
  std::istream &stream,
  std::string &line,
  std::size_t &line_number
)
{
  while(std::getline(stream, line)) {
    ++line_number;
    if(!trim(line).empty()) {
      return true;
    }
  }
  return false;
}

[[nodiscard]] bool parse_double_token(
  const std::string &token,
  double &value
)
{
  try {
    std::size_t parsed = 0U;
    const double parsed_value = std::stod(token, &parsed);
    if(parsed != token.size()) {
      return false;
    }
    value = parsed_value;
    return true;
  }
  catch(...) {
    return false;
  }
}

[[nodiscard]] bool parse_ascii_facet(
  std::istream &stream,
  const std::vector<std::string> &facet_tokens,
  std::size_t facet_line_number,
  std::size_t &line_number,
  RawFacet &facet,
  std::string &failure_reason
)
{
  failure_reason.clear();
  const auto fail_parse = [&](std::string message) {
    failure_reason = std::move(message);
    return false;
  };

  if(facet_tokens.size() < 5U ||
     lower_ascii(facet_tokens[0]) != "facet" ||
     lower_ascii(facet_tokens[1]) != "normal" ||
     !parse_double_token(facet_tokens[2], facet.file_normal[0]) ||
     !parse_double_token(facet_tokens[3], facet.file_normal[1]) ||
     !parse_double_token(facet_tokens[4], facet.file_normal[2])) {
    return fail_parse(
      "ASCII STL facet header on line " + std::to_string(facet_line_number) +
      " is malformed; expected `facet normal nx ny nz`."
    );
  }

  std::string line;
  if(!read_non_empty_line(stream, line, line_number)) {
    return fail_parse(
      "ASCII STL facet on line " + std::to_string(facet_line_number) +
      " ended before its `outer loop` record."
    );
  }
  const auto outer_tokens = tokenize(line);
  if(outer_tokens.size() < 2U ||
     lower_ascii(outer_tokens[0]) != "outer" ||
     lower_ascii(outer_tokens[1]) != "loop") {
    return fail_parse(
      "ASCII STL facet on line " + std::to_string(facet_line_number) +
      " is missing the `outer loop` record."
    );
  }

  for(std::size_t vertex_index = 0U; vertex_index < facet.vertices.size(); ++vertex_index) {
    if(!read_non_empty_line(stream, line, line_number)) {
      return fail_parse(
        "ASCII STL facet on line " + std::to_string(facet_line_number) +
        " ended before vertex record " + std::to_string(vertex_index + 1U) + '.'
      );
    }

    const auto vertex_tokens = tokenize(line);
    if(vertex_tokens.size() < 4U ||
       lower_ascii(vertex_tokens[0]) != "vertex" ||
       !parse_double_token(vertex_tokens[1], facet.vertices[vertex_index][0]) ||
       !parse_double_token(vertex_tokens[2], facet.vertices[vertex_index][1]) ||
       !parse_double_token(vertex_tokens[3], facet.vertices[vertex_index][2])) {
      return fail_parse(
        "ASCII STL facet vertex record on line " + std::to_string(line_number) +
        " is malformed; expected `vertex x y z`."
      );
    }
  }

  if(!read_non_empty_line(stream, line, line_number)) {
    return fail_parse(
      "ASCII STL facet on line " + std::to_string(facet_line_number) +
      " ended before its `endloop` record."
    );
  }
  if(lower_ascii(trim(line)) != "endloop") {
    return fail_parse(
      "ASCII STL facet on line " + std::to_string(facet_line_number) +
      " is missing the `endloop` record."
    );
  }

  if(!read_non_empty_line(stream, line, line_number)) {
    return fail_parse(
      "ASCII STL facet on line " + std::to_string(facet_line_number) +
      " ended before its `endfacet` record."
    );
  }
  if(lower_ascii(trim(line)) != "endfacet") {
    return fail_parse(
      "ASCII STL facet on line " + std::to_string(facet_line_number) +
      " is missing the `endfacet` record."
    );
  }

  return true;
}

[[nodiscard]] bool looks_like_text_stl_payload(const std::string &bytes) noexcept
{
  bool saw_newline = false;
  const std::size_t probe_size = std::min<std::size_t>(bytes.size(), 256U);
  for(std::size_t index = 0U; index < probe_size; ++index) {
    const unsigned char value = static_cast<unsigned char>(bytes[index]);
    if(value == '\n' || value == '\r') {
      saw_newline = true;
      continue;
    }
    if(value == '\t' || value == ' ') {
      continue;
    }
    if(value < 32U || value > 126U) {
      return false;
    }
  }

  return saw_newline;
}

[[nodiscard]] bool parse_ascii_stl(
  const std::string &bytes,
  std::vector<RawFacet> &facets,
  std::string &failure_reason
)
{
  facets.clear();
  failure_reason.clear();

  std::istringstream stream(bytes);
  std::string line;
  std::size_t line_number = 0U;
  if(!read_non_empty_line(stream, line, line_number)) {
    failure_reason = "ASCII STL import requires a non-empty `solid` block.";
    return false;
  }

  const auto header_tokens = tokenize(line);
  if(header_tokens.empty() || lower_ascii(header_tokens[0]) != "solid") {
    failure_reason = "ASCII STL must begin with a `solid` header.";
    return false;
  }

  bool saw_facet = false;
  while(read_non_empty_line(stream, line, line_number)) {
    const auto tokens = tokenize(line);
    if(tokens.empty()) {
      continue;
    }

    const auto keyword = lower_ascii(tokens[0]);
    if(keyword == "solid" || keyword == "endsolid") {
      continue;
    }
    if(keyword != "facet") {
      facets.clear();
      failure_reason =
        "ASCII STL encountered unexpected record `" + keyword + "` on line " +
        std::to_string(line_number) + "; expected `facet` or `endsolid`.";
      return false;
    }

    RawFacet facet;
    if(!parse_ascii_facet(
         stream,
         tokens,
         line_number,
         line_number,
         facet,
         failure_reason
       )) {
      facets.clear();
      return false;
    }

    facets.push_back(std::move(facet));
    saw_facet = true;
  }

  if(!saw_facet) {
    facets.clear();
    failure_reason = "ASCII STL did not contain any `facet` records.";
  }
  return saw_facet;
}

[[nodiscard]] std::uint32_t read_le_u32(const char *bytes) noexcept
{
  return static_cast<std::uint32_t>(static_cast<unsigned char>(bytes[0])) |
         (static_cast<std::uint32_t>(static_cast<unsigned char>(bytes[1])) << 8U) |
         (static_cast<std::uint32_t>(static_cast<unsigned char>(bytes[2])) << 16U) |
         (static_cast<std::uint32_t>(static_cast<unsigned char>(bytes[3])) << 24U);
}

[[nodiscard]] double read_le_f32(const char *bytes) noexcept
{
  const std::uint32_t packed = read_le_u32(bytes);
  float value = 0.0F;
  std::memcpy(&value, &packed, sizeof(value));
  return static_cast<double>(value);
}

[[nodiscard]] bool parse_binary_stl(
  const std::string &bytes,
  std::vector<RawFacet> &facets,
  std::string &failure_reason
)
{
  facets.clear();
  failure_reason.clear();

  if(bytes.size() < 84U) {
    failure_reason = "Binary STL is smaller than the required 84-byte header.";
    return false;
  }

  const std::uint32_t facet_count = read_le_u32(bytes.data() + 80U);
  if(facet_count == 0U) {
    failure_reason = "Binary STL declares zero facets.";
    return false;
  }
  const std::uint64_t expected_size =
    84ULL + static_cast<std::uint64_t>(facet_count) * 50ULL;
  if(static_cast<std::uint64_t>(bytes.size()) < expected_size) {
    failure_reason =
      "Binary STL declares " + std::to_string(facet_count) +
      " facet(s) but the file is " + std::to_string(bytes.size()) +
      " bytes; expected at least " + std::to_string(expected_size) + " bytes.";
    return false;
  }

  facets.reserve(facet_count);
  std::size_t offset = 84U;
  for(std::uint32_t facet_index = 0U; facet_index < facet_count; ++facet_index) {
    RawFacet facet;
    facet.file_normal = {
      read_le_f32(bytes.data() + offset),
      read_le_f32(bytes.data() + offset + 4U),
      read_le_f32(bytes.data() + offset + 8U),
    };
    offset += 12U;

    for(std::size_t vertex_index = 0U; vertex_index < facet.vertices.size(); ++vertex_index) {
      facet.vertices[vertex_index] = {
        read_le_f32(bytes.data() + offset),
        read_le_f32(bytes.data() + offset + 4U),
        read_le_f32(bytes.data() + offset + 8U),
      };
      offset += 12U;
    }

    offset += 2U;
    facets.push_back(std::move(facet));
  }

  return !facets.empty();
}

[[nodiscard]] bool looks_like_ascii_stl(const std::string &bytes) noexcept
{
  if(bytes.size() < 5U ||
     lower_ascii(std::string_view(bytes.data(), 5U)) != "solid") {
    return false;
  }

  return looks_like_text_stl_payload(bytes);
}

[[nodiscard]] base::StatusCode read_stl_facets(
  std::string_view path,
  std::vector<RawFacet> &facets
) noexcept
{
  facets.clear();

  std::ifstream stream(std::string(path), std::ios::binary);
  if(!stream) {
    return core::detail::publish_error(
      base::StatusCode::io_error,
      "STL import could not open the requested file."
    );
  }

  stream.seekg(0, std::ios::end);
  const auto size = stream.tellg();
  if(size <= 0) {
    return core::detail::publish_error(
      base::StatusCode::io_error,
      "STL import requires a non-empty file."
    );
  }
  stream.seekg(0, std::ios::beg);

  std::string bytes(static_cast<std::size_t>(size), '\0');
  stream.read(bytes.data(), static_cast<std::streamsize>(bytes.size()));
  if(!stream) {
    return core::detail::publish_error(
      base::StatusCode::io_error,
      "STL import failed while reading the input file."
    );
  }

  const bool ascii_hint = looks_like_ascii_stl(bytes);
  const bool text_payload_hint = looks_like_text_stl_payload(bytes);
  std::string ascii_failure;
  if(ascii_hint && parse_ascii_stl(bytes, facets, ascii_failure)) {
    return base::StatusCode::ok;
  }
  std::string binary_failure;
  if(parse_binary_stl(bytes, facets, binary_failure)) {
    return base::StatusCode::ok;
  }
  if(parse_ascii_stl(bytes, facets, ascii_failure)) {
    return base::StatusCode::ok;
  }

  if(!ascii_failure.empty() &&
     (ascii_hint || text_payload_hint || binary_failure.empty())) {
    return core::detail::publish_error(
      base::StatusCode::io_error,
      ascii_failure
    );
  }
  if(!binary_failure.empty()) {
    return core::detail::publish_error(
      base::StatusCode::io_error,
      binary_failure
    );
  }

  return core::detail::publish_error(
    base::StatusCode::io_error,
    "STL import could not parse the file as ASCII or binary STL."
  );
}

void update_bbox(BoundingBox &bbox, const geo::Point3 &point) noexcept
{
  if(!bbox.initialized) {
    bbox.min = point;
    bbox.max = point;
    bbox.initialized = true;
    return;
  }

  for(std::size_t axis = 0U; axis < 3U; ++axis) {
    bbox.min[axis] = std::min(bbox.min[axis], point[axis]);
    bbox.max[axis] = std::max(bbox.max[axis], point[axis]);
  }
}

[[nodiscard]] CellKey cell_key_for_point(
  const geo::Point3 &point,
  double tolerance
) noexcept
{
  return CellKey {
    static_cast<long long>(std::floor(point[0] / tolerance)),
    static_cast<long long>(std::floor(point[1] / tolerance)),
    static_cast<long long>(std::floor(point[2] / tolerance)),
  };
}

[[nodiscard]] std::uint32_t weld_vertex(
  const geo::Point3 &point,
  double tolerance,
  double tolerance_squared,
  std::unordered_map<CellKey, std::vector<std::uint32_t>, CellKeyHash> &buckets,
  std::vector<VertexData> &vertices
)
{
  const CellKey base_key = cell_key_for_point(point, tolerance);
  for(long long dx = -1; dx <= 1; ++dx) {
    for(long long dy = -1; dy <= 1; ++dy) {
      for(long long dz = -1; dz <= 1; ++dz) {
        const CellKey candidate_key {
          base_key.x + dx,
          base_key.y + dy,
          base_key.z + dz,
        };
        const auto bucket_it = buckets.find(candidate_key);
        if(bucket_it == buckets.end()) {
          continue;
        }

        for(const std::uint32_t index : bucket_it->second) {
          if(point_distance_squared(vertices[index].position, point) <= tolerance_squared) {
            return index;
          }
        }
      }
    }
  }

  const auto index = static_cast<std::uint32_t>(vertices.size());
  vertices.push_back(VertexData {point, {}});
  buckets[base_key].push_back(index);
  return index;
}

[[nodiscard]] std::uint32_t other_vertex(
  const EdgeData &edge,
  std::uint32_t vertex
) noexcept
{
  return edge.vertices[0] == vertex ? edge.vertices[1] : edge.vertices[0];
}

void build_boundary_loops(DiscreteModel &model) noexcept
{
  for(std::size_t face_index = 0U; face_index < model.faces.size(); ++face_index) {
    auto &face = model.faces[face_index];

    std::unordered_map<std::uint32_t, std::vector<std::uint32_t>> boundary_adjacency;
    std::vector<std::uint32_t> boundary_edges;
    for(const std::uint32_t edge_index : face.child_edges) {
      const auto &edge = model.edges[edge_index];
      if(edge.triangles.size() != 1U) {
        continue;
      }
      boundary_edges.push_back(edge_index);
      boundary_adjacency[edge.vertices[0]].push_back(edge_index);
      boundary_adjacency[edge.vertices[1]].push_back(edge_index);
    }

    if(boundary_edges.empty()) {
      continue;
    }

    std::unordered_set<std::uint32_t> used_edges;
    while(used_edges.size() < boundary_edges.size()) {
      std::uint32_t start_edge = boundary_edges.front();
      for(const std::uint32_t candidate_edge : boundary_edges) {
        if(used_edges.find(candidate_edge) == used_edges.end()) {
          start_edge = candidate_edge;
          break;
        }
      }

      const auto &edge = model.edges[start_edge];
      std::uint32_t start_vertex = edge.vertices[0];
      const auto left_degree = boundary_adjacency[edge.vertices[0]].size();
      const auto right_degree = boundary_adjacency[edge.vertices[1]].size();
      if(left_degree == 2U && right_degree != 2U) {
        start_vertex = edge.vertices[1];
      } else if(left_degree != 2U) {
        start_vertex = edge.vertices[0];
      }

      geo::FaceBoundaryLoop loop;
      loop.kind = geo::FaceBoundaryLoopKind::unknown;

      std::uint32_t current_vertex = start_vertex;
      std::uint32_t current_edge = start_edge;
      while(true) {
        used_edges.insert(current_edge);

        const auto &current_edge_data = model.edges[current_edge];
        const std::uint32_t next_vertex = other_vertex(current_edge_data, current_vertex);
        loop.edge_uses.push_back(geo::FaceBoundaryEdgeUse {
          make_entity_id(geo::TopologyDimension::edge, current_edge),
          make_entity_id(geo::TopologyDimension::vertex, current_vertex),
          make_entity_id(geo::TopologyDimension::vertex, next_vertex),
          current_vertex == current_edge_data.vertices[0],
        });

        current_vertex = next_vertex;
        std::uint32_t next_edge = geo::invalid_topology_index;
        const auto adjacency_it = boundary_adjacency.find(current_vertex);
        if(adjacency_it != boundary_adjacency.end()) {
          for(const std::uint32_t candidate_edge : adjacency_it->second) {
            if(used_edges.find(candidate_edge) == used_edges.end()) {
              next_edge = candidate_edge;
              break;
            }
          }
        }

        if(next_edge == geo::invalid_topology_index) {
          loop.closed =
            !loop.edge_uses.empty() &&
            loop.edge_uses.front().start_vertex == loop.edge_uses.back().end_vertex;
          break;
        }

        current_edge = next_edge;
      }

      if(!loop.edge_uses.empty()) {
        face.boundary_loops.push_back(std::move(loop));
      }
    }
  }
}

[[nodiscard]] base::StatusCode build_discrete_model(
  const std::vector<RawFacet> &facets,
  const geo::StlImportOptions &options,
  DiscreteModel &model
) noexcept
{
  model = {};

  BoundingBox bbox;
  for(const auto &facet : facets) {
    for(const auto &vertex : facet.vertices) {
      update_bbox(bbox, vertex);
    }
  }
  if(!bbox.initialized) {
    return core::detail::publish_error(
      base::StatusCode::io_error,
      "STL import did not find any facet vertices."
    );
  }

  const double diagonal = point_distance(bbox.max, bbox.min);
  const double weld_tolerance = std::max(
    kMinimumMergeTolerance,
    diagonal * options.relative_merge_tolerance
  );
  const double weld_tolerance_squared = weld_tolerance * weld_tolerance;

  std::unordered_map<CellKey, std::vector<std::uint32_t>, CellKeyHash> vertex_buckets;
  std::unordered_set<TriangleKey, TriangleKeyHash> unique_triangles;

  for(const auto &facet : facets) {
    std::array<std::uint32_t, 3> triangle_vertices {};
    for(std::size_t vertex_index = 0U; vertex_index < triangle_vertices.size(); ++vertex_index) {
      triangle_vertices[vertex_index] = weld_vertex(
        facet.vertices[vertex_index],
        weld_tolerance,
        weld_tolerance_squared,
        vertex_buckets,
        model.vertices
      );
    }

    if(triangle_vertices[0] == triangle_vertices[1] ||
       triangle_vertices[1] == triangle_vertices[2] ||
       triangle_vertices[0] == triangle_vertices[2]) {
      continue;
    }

    TriangleKey key {triangle_vertices};
    std::sort(key.vertices.begin(), key.vertices.end());
    if(!unique_triangles.insert(key).second) {
      continue;
    }

    const auto &p0 = model.vertices[triangle_vertices[0]].position;
    const auto &p1 = model.vertices[triangle_vertices[1]].position;
    const auto &p2 = model.vertices[triangle_vertices[2]].position;
    geo::Vector3 normal = cross_product(subtract(p1, p0), subtract(p2, p0));
    const double normal_length = norm(normal);
    if(normal_length <= weld_tolerance) {
      normal = facet.file_normal;
    }
    normal = normalized(normal);
    if(norm(normal) <= 0.0) {
      continue;
    }

    model.triangles.push_back(TriangleData {
      triangle_vertices,
      {},
      normal,
    });
  }

  if(model.triangles.empty()) {
    return core::detail::publish_error(
      base::StatusCode::io_error,
      "STL import did not find any non-degenerate triangles."
    );
  }

  std::unordered_map<EdgeKey, std::uint32_t, EdgeKeyHash> edge_lookup;
  for(std::size_t triangle_index = 0U; triangle_index < model.triangles.size(); ++triangle_index) {
    auto &triangle = model.triangles[triangle_index];
    for(std::size_t edge_local_index = 0U; edge_local_index < triangle.edges.size(); ++edge_local_index) {
      const std::uint32_t first = triangle.vertices[edge_local_index];
      const std::uint32_t second =
        triangle.vertices[(edge_local_index + 1U) % triangle.vertices.size()];
      const EdgeKey key {
        std::min(first, second),
        std::max(first, second),
      };

      const auto edge_it = edge_lookup.find(key);
      std::uint32_t edge_index = 0U;
      if(edge_it == edge_lookup.end()) {
        edge_index = static_cast<std::uint32_t>(model.edges.size());
        model.edges.push_back(EdgeData {{key.v0, key.v1}, {}, {}});
        edge_lookup.emplace(key, edge_index);
      } else {
        edge_index = edge_it->second;
      }

      triangle.edges[edge_local_index] = edge_index;
      model.edges[edge_index].triangles.push_back(
        static_cast<std::uint32_t>(triangle_index)
      );
    }
  }

  std::vector<bool> visited_triangles(model.triangles.size(), false);
  for(std::size_t triangle_index = 0U; triangle_index < model.triangles.size(); ++triangle_index) {
    if(visited_triangles[triangle_index]) {
      continue;
    }

    FaceData face;
    std::unordered_set<std::uint32_t> face_edges;
    std::queue<std::uint32_t> pending;
    pending.push(static_cast<std::uint32_t>(triangle_index));
    visited_triangles[triangle_index] = true;

    while(!pending.empty()) {
      const std::uint32_t current_triangle = pending.front();
      pending.pop();
      face.triangles.push_back(current_triangle);

      for(const std::uint32_t edge_index : model.triangles[current_triangle].edges) {
        if(face_edges.insert(edge_index).second) {
          face.child_edges.push_back(edge_index);
        }

        for(const std::uint32_t neighbor_triangle : model.edges[edge_index].triangles) {
          if(!visited_triangles[neighbor_triangle]) {
            visited_triangles[neighbor_triangle] = true;
            pending.push(neighbor_triangle);
          }
        }
      }
    }

    const auto face_index = static_cast<std::uint32_t>(model.faces.size());
    for(const std::uint32_t edge_index : face.child_edges) {
      model.edges[edge_index].parent_faces.push_back(face_index);
    }
    model.faces.push_back(std::move(face));
  }

  for(std::size_t edge_index = 0U; edge_index < model.edges.size(); ++edge_index) {
    const auto &edge = model.edges[edge_index];
    model.vertices[edge.vertices[0]].parent_edges.push_back(
      static_cast<std::uint32_t>(edge_index)
    );
    model.vertices[edge.vertices[1]].parent_edges.push_back(
      static_cast<std::uint32_t>(edge_index)
    );
  }

  build_boundary_loops(model);

  model.summary.entity_count =
    model.faces.size() + model.edges.size() + model.vertices.size();
  model.summary.face_count = model.faces.size();
  model.summary.edge_count = model.edges.size();
  model.summary.vertex_count = model.vertices.size();
  for(const auto &face : model.faces) {
    bool closed_component = true;
    for(const std::uint32_t edge_index : face.child_edges) {
      if(model.edges[edge_index].triangles.size() != 2U) {
        closed_component = false;
        break;
      }
    }
    if(closed_component) {
      ++model.summary.shell_count;
    }
  }

  return base::StatusCode::ok;
}

class DiscreteModelStorage final : public GeometryModelStorage {
public:
  explicit DiscreteModelStorage(DiscreteModel model) : model_(std::move(model)) {}

  [[nodiscard]] GeometryKernel kernel() const noexcept override
  {
    return GeometryKernel::discrete;
  }

  [[nodiscard]] base::StatusCode coarse_proxy_mesh(
    GeometryCoarseProxyMesh &proxy_mesh
  ) const noexcept override
  {
    proxy_mesh = {};
    return core::detail::publish_error(
      base::StatusCode::unsupported,
      "Discrete STL-backed models do not expose an internal OCC proxy triangulation cache."
    );
  }

  [[nodiscard]] base::StatusCode topology_snapshot(
    geo::TopologySnapshot &snapshot
  ) const noexcept override
  {
    snapshot = {};
    snapshot.topology_revision = 1U;
    snapshot.faces.reserve(model_.faces.size());
    snapshot.edges.reserve(model_.edges.size());
    snapshot.vertices.reserve(model_.vertices.size());

    for(std::size_t face_index = 0U; face_index < model_.faces.size(); ++face_index) {
      snapshot.faces.push_back(geo::TopologyEntityInfo {
        make_entity_id(geo::TopologyDimension::face, face_index),
        0U,
        model_.faces[face_index].child_edges.size(),
      });
    }
    for(std::size_t edge_index = 0U; edge_index < model_.edges.size(); ++edge_index) {
      snapshot.edges.push_back(geo::TopologyEntityInfo {
        make_entity_id(geo::TopologyDimension::edge, edge_index),
        model_.edges[edge_index].parent_faces.size(),
        2U,
      });
    }
    for(std::size_t vertex_index = 0U; vertex_index < model_.vertices.size(); ++vertex_index) {
      snapshot.vertices.push_back(geo::TopologyEntityInfo {
        make_entity_id(geo::TopologyDimension::vertex, vertex_index),
        model_.vertices[vertex_index].parent_edges.size(),
        0U,
      });
    }

    return base::StatusCode::ok;
  }

  [[nodiscard]] base::StatusCode topology_children(
    geo::TopologyEntityId entity,
    std::vector<geo::TopologyEntityId> &children
  ) const noexcept override
  {
    children.clear();

    switch(entity.dimension) {
    case geo::TopologyDimension::face: {
      const auto *face = find_face(
        entity,
        "Discrete STL face-child lookup requires a face topology entity id.",
        "Discrete STL face-child lookup received an out-of-range face entity id."
      );
      if(face == nullptr) {
        return base::last_error_code();
      }
      children.reserve(face->child_edges.size());
      for(const std::uint32_t edge_index : face->child_edges) {
        children.push_back(make_entity_id(geo::TopologyDimension::edge, edge_index));
      }
      return base::StatusCode::ok;
    }
    case geo::TopologyDimension::edge: {
      const auto *edge = find_edge(
        entity,
        "Discrete STL edge-child lookup requires an edge topology entity id.",
        "Discrete STL edge-child lookup received an out-of-range edge entity id."
      );
      if(edge == nullptr) {
        return base::last_error_code();
      }
      children.reserve(2U);
      children.push_back(make_entity_id(geo::TopologyDimension::vertex, edge->vertices[0]));
      children.push_back(make_entity_id(geo::TopologyDimension::vertex, edge->vertices[1]));
      return base::StatusCode::ok;
    }
    case geo::TopologyDimension::vertex: {
      const auto *vertex = find_vertex(
        entity,
        "Discrete STL vertex-child lookup requires a vertex topology entity id.",
        "Discrete STL vertex-child lookup received an out-of-range vertex entity id."
      );
      if(vertex == nullptr) {
        return base::last_error_code();
      }
      static_cast<void>(vertex);
      return base::StatusCode::ok;
    }
    case geo::TopologyDimension::region:
    default:
      return core::detail::publish_error(
        base::StatusCode::invalid_argument,
        "Discrete STL models do not expose region entities."
      );
    }
  }

  [[nodiscard]] base::StatusCode topology_parents(
    geo::TopologyEntityId entity,
    std::vector<geo::TopologyEntityId> &parents
  ) const noexcept override
  {
    parents.clear();

    switch(entity.dimension) {
    case geo::TopologyDimension::face: {
      const auto *face = find_face(
        entity,
        "Discrete STL face-parent lookup requires a face topology entity id.",
        "Discrete STL face-parent lookup received an out-of-range face entity id."
      );
      if(face == nullptr) {
        return base::last_error_code();
      }
      static_cast<void>(face);
      return base::StatusCode::ok;
    }
    case geo::TopologyDimension::edge: {
      const auto *edge = find_edge(
        entity,
        "Discrete STL edge-parent lookup requires an edge topology entity id.",
        "Discrete STL edge-parent lookup received an out-of-range edge entity id."
      );
      if(edge == nullptr) {
        return base::last_error_code();
      }
      parents.reserve(edge->parent_faces.size());
      for(const std::uint32_t face_index : edge->parent_faces) {
        parents.push_back(make_entity_id(geo::TopologyDimension::face, face_index));
      }
      return base::StatusCode::ok;
    }
    case geo::TopologyDimension::vertex: {
      const auto *vertex = find_vertex(
        entity,
        "Discrete STL vertex-parent lookup requires a vertex topology entity id.",
        "Discrete STL vertex-parent lookup received an out-of-range vertex entity id."
      );
      if(vertex == nullptr) {
        return base::last_error_code();
      }
      parents.reserve(vertex->parent_edges.size());
      for(const std::uint32_t edge_index : vertex->parent_edges) {
        parents.push_back(make_entity_id(geo::TopologyDimension::edge, edge_index));
      }
      return base::StatusCode::ok;
    }
    case geo::TopologyDimension::region:
    default:
      return core::detail::publish_error(
        base::StatusCode::invalid_argument,
        "Discrete STL models do not expose region entities."
      );
    }
  }

  [[nodiscard]] base::StatusCode face_uv_bounds(
    geo::TopologyEntityId face_entity,
    geo::FaceUvBounds &bounds
  ) const noexcept override
  {
    bounds = {};
    if(find_face(
         face_entity,
         "Discrete STL face UV bounds require a face topology entity id.",
         "Discrete STL face UV bounds received an out-of-range face entity id."
       ) == nullptr) {
      return base::last_error_code();
    }
    return unsupported_face_query(
      "Discrete STL-backed models do not expose parametric face UV bounds."
    );
  }

  [[nodiscard]] base::StatusCode sample_face(
    geo::TopologyEntityId face_entity,
    double,
    double,
    geo::FaceSample &sample
  ) const noexcept override
  {
    sample = {};
    if(find_face(
         face_entity,
         "Discrete STL face sampling requires a face topology entity id.",
         "Discrete STL face sampling received an out-of-range face entity id."
       ) == nullptr) {
      return base::last_error_code();
    }
    return unsupported_face_query(
      "Discrete STL-backed models do not support parametric face sampling."
    );
  }

  [[nodiscard]] base::StatusCode sample_face_curvature(
    geo::TopologyEntityId face_entity,
    double,
    double,
    geo::FaceCurvatureSample &sample
  ) const noexcept override
  {
    sample = {};
    if(find_face(
         face_entity,
         "Discrete STL face curvature sampling requires a face topology entity id.",
         "Discrete STL face curvature sampling received an out-of-range face entity id."
       ) == nullptr) {
      return base::last_error_code();
    }
    return unsupported_face_query(
      "Discrete STL-backed models do not expose continuous face curvature."
    );
  }

  [[nodiscard]] base::StatusCode sample_face_derivatives(
    geo::TopologyEntityId face_entity,
    double,
    double,
    geo::FaceDerivatives &sample
  ) const noexcept override
  {
    sample = {};
    if(find_face(
         face_entity,
         "Discrete STL face derivative sampling requires a face topology entity id.",
         "Discrete STL face derivative sampling received an out-of-range face entity id."
       ) == nullptr) {
      return base::last_error_code();
    }
    return unsupported_face_query(
      "Discrete STL-backed models do not expose parametric face derivatives."
    );
  }

  [[nodiscard]] base::StatusCode project_point_to_face(
    geo::TopologyEntityId face_entity,
    const geo::Point3 &,
    geo::FaceProjection &projection
  ) const noexcept override
  {
    projection = {};
    if(find_face(
         face_entity,
         "Discrete STL face projection requires a face topology entity id.",
         "Discrete STL face projection received an out-of-range face entity id."
       ) == nullptr) {
      return base::last_error_code();
    }
    return unsupported_face_query(
      "Discrete STL-backed models do not support point-to-face projection."
    );
  }

  [[nodiscard]] base::StatusCode recover_face_uv(
    geo::TopologyEntityId face_entity,
    const geo::Point3 &,
    geo::FaceUvMapping &mapping
  ) const noexcept override
  {
    mapping = {};
    if(find_face(
         face_entity,
         "Discrete STL face inverse mapping requires a face topology entity id.",
         "Discrete STL face inverse mapping received an out-of-range face entity id."
       ) == nullptr) {
      return base::last_error_code();
    }
    return unsupported_face_query(
      "Discrete STL-backed models do not expose face UV coordinates."
    );
  }

  [[nodiscard]] base::StatusCode recover_face_uv_from_edge(
    geo::TopologyEntityId /*face_entity*/,
    geo::TopologyEntityId /*edge_entity*/,
    double /*edge_parameter*/,
    geo::FaceUvMapping &mapping
  ) const noexcept override
  {
    mapping = {};
    return unsupported_face_query(
      "Discrete STL-backed models do not expose edge pcurve UV coordinates."
    );
  }

  [[nodiscard]] base::StatusCode recover_face_uv_from_edge_use(
    geo::TopologyEntityId /*face_entity*/,
    const geo::FaceBoundaryEdgeUse & /*edge_use*/,
    double /*edge_parameter*/,
    geo::FaceUvMapping &mapping
  ) const noexcept override
  {
    mapping = {};
    return unsupported_face_query(
      "Discrete STL-backed models do not expose face-local edge-use pcurve UV coordinates."
    );
  }

  [[nodiscard]] base::StatusCode edge_curve_info(
    geo::TopologyEntityId edge_entity,
    geo::EdgeCurveInfo &info
  ) const noexcept override
  {
    info = {};

    const auto *edge = find_edge(
      edge_entity,
      "Discrete STL edge geometry lookup requires an edge topology entity id.",
      "Discrete STL edge geometry lookup received an out-of-range edge entity id."
    );
    if(edge == nullptr) {
      return base::last_error_code();
    }

    const auto &start = model_.vertices[edge->vertices[0]].position;
    const auto &end = model_.vertices[edge->vertices[1]].position;
    info.edge = edge_entity;
    info.parameter_min = 0.0;
    info.parameter_max = 1.0;
    info.start_point = start;
    info.end_point = end;
    info.approximate_length = point_distance(start, end);
    return base::StatusCode::ok;
  }

  [[nodiscard]] base::StatusCode sample_edge_tangent(
    geo::TopologyEntityId edge_entity,
    double parameter,
    geo::EdgeTangentSample &sample
  ) const noexcept override
  {
    sample = {};

    const auto *edge = find_edge(
      edge_entity,
      "Discrete STL edge tangent sampling requires an edge topology entity id.",
      "Discrete STL edge tangent sampling received an out-of-range edge entity id."
    );
    if(edge == nullptr) {
      return base::last_error_code();
    }
    if(parameter < 0.0 || parameter > 1.0) {
      return core::detail::publish_error(
        base::StatusCode::invalid_argument,
        "Discrete STL edge tangent sampling requires parameter in [0, 1]."
      );
    }

    const auto &start = model_.vertices[edge->vertices[0]].position;
    const auto &end = model_.vertices[edge->vertices[1]].position;
    const geo::Vector3 derivative = subtract(end, start);
    const double speed = norm(derivative);

    sample.edge = edge_entity;
    sample.parameter = parameter;
    sample.position = interpolate(start, end, parameter);
    sample.derivative = derivative;
    sample.speed = speed;
    if(speed > 0.0) {
      sample.tangent = {
        derivative[0] / speed,
        derivative[1] / speed,
        derivative[2] / speed,
      };
      sample.tangent_defined = true;
    }

    return base::StatusCode::ok;
  }

  [[nodiscard]] base::StatusCode face_boundary_loops(
    geo::TopologyEntityId face_entity,
    geo::FaceBoundaryLoops &boundary
  ) const noexcept override
  {
    boundary = {};

    const auto *face = find_face(
      face_entity,
      "Discrete STL face boundary loop lookup requires a face topology entity id.",
      "Discrete STL face boundary loop lookup received an out-of-range face entity id."
    );
    if(face == nullptr) {
      return base::last_error_code();
    }

    boundary.face = face_entity;
    boundary.loops = face->boundary_loops;
    return base::StatusCode::ok;
  }

  [[nodiscard]] base::StatusCode feature_edges(
    geo::FeatureEdgeReport &report,
    const geo::FeatureEdgeOptions &options
  ) const noexcept override
  {
    report = {};

    for(std::size_t edge_index = 0U; edge_index < model_.edges.size(); ++edge_index) {
      const auto &edge = model_.edges[edge_index];
      const auto edge_id = make_entity_id(geo::TopologyDimension::edge, edge_index);

      if(edge.triangles.size() == 1U) {
        if(options.include_boundary_edges) {
          report.edges.push_back(edge_id);
          ++report.boundary_edge_count;
        }
        continue;
      }

      if(edge.triangles.size() > 2U) {
        if(options.include_non_manifold_edges) {
          report.edges.push_back(edge_id);
          ++report.non_manifold_edge_count;
        }
        continue;
      }

      if(edge.triangles.size() != 2U) {
        continue;
      }

      const auto &triangle0 = model_.triangles[edge.triangles[0]];
      const auto &triangle1 = model_.triangles[edge.triangles[1]];
      const double alignment = std::clamp(
        std::abs(dot_product(triangle0.normal, triangle1.normal)),
        0.0,
        1.0
      );
      const double angle_degrees =
        std::acos(alignment) * (180.0 / kPi);
      if(angle_degrees > options.feature_angle_degrees) {
        report.edges.push_back(edge_id);
        ++report.sharp_edge_count;
      }
    }

    return base::StatusCode::ok;
  }

private:
  [[nodiscard]] const FaceData *find_face(
    geo::TopologyEntityId entity,
    std::string_view wrong_dimension_message,
    std::string_view out_of_range_message
  ) const noexcept
  {
    if(entity.dimension != geo::TopologyDimension::face) {
      static_cast<void>(core::detail::publish_error(
        base::StatusCode::invalid_argument,
        wrong_dimension_message
      ));
      return nullptr;
    }
    if(entity.index >= model_.faces.size()) {
      static_cast<void>(core::detail::publish_error(
        base::StatusCode::invalid_argument,
        out_of_range_message
      ));
      return nullptr;
    }
    return &model_.faces[entity.index];
  }

  [[nodiscard]] const EdgeData *find_edge(
    geo::TopologyEntityId entity,
    std::string_view wrong_dimension_message,
    std::string_view out_of_range_message
  ) const noexcept
  {
    if(entity.dimension != geo::TopologyDimension::edge) {
      static_cast<void>(core::detail::publish_error(
        base::StatusCode::invalid_argument,
        wrong_dimension_message
      ));
      return nullptr;
    }
    if(entity.index >= model_.edges.size()) {
      static_cast<void>(core::detail::publish_error(
        base::StatusCode::invalid_argument,
        out_of_range_message
      ));
      return nullptr;
    }
    return &model_.edges[entity.index];
  }

  [[nodiscard]] const VertexData *find_vertex(
    geo::TopologyEntityId entity,
    std::string_view wrong_dimension_message,
    std::string_view out_of_range_message
  ) const noexcept
  {
    if(entity.dimension != geo::TopologyDimension::vertex) {
      static_cast<void>(core::detail::publish_error(
        base::StatusCode::invalid_argument,
        wrong_dimension_message
      ));
      return nullptr;
    }
    if(entity.index >= model_.vertices.size()) {
      static_cast<void>(core::detail::publish_error(
        base::StatusCode::invalid_argument,
        out_of_range_message
      ));
      return nullptr;
    }
    return &model_.vertices[entity.index];
  }

  [[nodiscard]] static base::StatusCode unsupported_face_query(
    std::string_view message
  ) noexcept
  {
    return core::detail::publish_error(base::StatusCode::unsupported, message);
  }

  DiscreteModel model_ {};
};

} // namespace

base::StatusCode import_stl(
  std::string_view path,
  geo::ModelHandle &model_handle,
  const geo::StlImportOptions &options,
  base::ContextHandle context_handle
) noexcept
{
  model_handle = sqmesh::invalid_handle;

  if(path.empty()) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "STL import requires a non-empty path."
    );
  }
  if(options.relative_merge_tolerance < 0.0) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "STL import requires relative_merge_tolerance to be non-negative."
    );
  }

  std::vector<RawFacet> facets;
  const auto read_status = read_stl_facets(path, facets);
  if(read_status != base::StatusCode::ok) {
    return read_status;
  }

  DiscreteModel model;
  const auto build_status = build_discrete_model(facets, options, model);
  if(build_status != base::StatusCode::ok) {
    return build_status;
  }

  const auto summary = model.summary;
  auto storage = std::make_shared<DiscreteModelStorage>(std::move(model));
  return core::detail::store_model(
    std::move(storage),
    summary,
    model_handle,
    context_handle
  );
}

} // namespace sqmesh::model::detail
