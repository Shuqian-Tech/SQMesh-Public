// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "sqmesh/sqmesh.hpp"

#include "../src/mesh/io/mesh_io_test_hook.hpp"

#include "../src/core/runtime_registry.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#ifdef sqmesh_TEST_CGNS_ENABLED
#include <cgnslib.h>
#endif

namespace {

bool expect(bool condition, const char *message)
{
  if(condition) {
    return true;
  }

  const auto &detail = sqmesh::base::last_error_message();
  if(detail.empty()) {
    std::fprintf(stderr, "mesh_io_smoke: %s\n", message);
  }
  else {
    const std::string detail_text(detail);
    std::fprintf(stderr, "mesh_io_smoke: %s (%s)\n", message, detail_text.c_str());
  }
  return false;
}

bool write_text_file(
  const std::filesystem::path &path,
  const char *contents
)
{
  std::ofstream output(path, std::ios::trunc);
  if(!output) {
    std::fprintf(stderr, "mesh_io_smoke: failed to open %s for writing\n", path.string().c_str());
    return false;
  }

  output << contents;
  return static_cast<bool>(output);
}

bool write_text_file(
  const std::filesystem::path &path,
  const std::string &contents
)
{
  return write_text_file(path, contents.c_str());
}

bool read_text_file(
  const std::filesystem::path &path,
  std::string &contents
)
{
  std::ifstream input(path);
  if(!input) {
    std::fprintf(stderr, "mesh_io_smoke: failed to open %s for reading\n", path.string().c_str());
    return false;
  }

  std::ostringstream stream;
  stream << input.rdbuf();
  if(!input.good() && !input.eof()) {
    std::fprintf(stderr, "mesh_io_smoke: failed to read %s\n", path.string().c_str());
    return false;
  }

  contents = stream.str();
  return true;
}

bool store_manual_mesh(
  sqmesh::mesh::Domain domain,
  const char *algorithm_name,
  sqmesh::mesh::MeshHandle &mesh_handle,
  sqmesh::base::ContextHandle context
)
{
  mesh_handle = sqmesh::invalid_handle;

  sqmesh::geo::ModelHandle model_handle = sqmesh::invalid_handle;
  if(sqmesh::geo::create_placeholder_model(model_handle, context) !=
     sqmesh::base::StatusCode::ok) {
    return false;
  }

  auto domain_storage =
    std::make_shared<sqmesh::mesh::Domain>(std::move(domain));
  const auto summary = domain_storage->summary();
  return sqmesh::core::detail::store_mesh(
           model_handle,
           algorithm_name,
           summary,
           std::move(domain_storage),
           mesh_handle,
           context
         ) == sqmesh::base::StatusCode::ok;
}

bool read_binary_file(
  const std::filesystem::path &path,
  std::string &contents
)
{
  std::ifstream input(path, std::ios::binary);
  if(!input) {
    std::fprintf(stderr, "mesh_io_smoke: failed to open %s for reading\n", path.string().c_str());
    return false;
  }

  std::ostringstream stream;
  stream << input.rdbuf();
  if(!input.good() && !input.eof()) {
    std::fprintf(stderr, "mesh_io_smoke: failed to read %s\n", path.string().c_str());
    return false;
  }

  contents = stream.str();
  return true;
}

template <typename T>
bool write_binary_value(
  std::ofstream &output,
  const T &value
)
{
  output.write(reinterpret_cast<const char *>(&value), static_cast<std::streamsize>(sizeof(T)));
  return static_cast<bool>(output);
}

template <typename T>
bool write_binary_values(
  std::ofstream &output,
  const T *values,
  std::size_t count
)
{
  if(count == 0U) {
    return true;
  }

  output.write(
    reinterpret_cast<const char *>(values),
    static_cast<std::streamsize>(sizeof(T) * count)
  );
  return static_cast<bool>(output);
}

bool write_msh22_binary_fixture(const std::filesystem::path &path)
{
  std::ofstream output(path, std::ios::trunc | std::ios::binary);
  if(!output) {
    std::fprintf(stderr, "mesh_io_smoke: failed to open %s for writing\n", path.string().c_str());
    return false;
  }

  output << "$MeshFormat\n2.2 1 " << sizeof(double) << "\n";
  const int one = 1;
  if(!write_binary_value(output, one)) {
    return false;
  }
  output << "\n$EndMeshFormat\n";
  output << "$PhysicalNames\n3\n";
  output << "1 11 \"feature_edges\"\n";
  output << "2 21 \"wall\"\n";
  output << "3 31 \"volume\"\n";
  output << "$EndPhysicalNames\n";

  output << "$Nodes\n4\n";
  const std::array<std::array<double, 3>, 4> nodes = {{
    {0.0, 0.0, 0.0},
    {1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0},
  }};
  for(int node_id = 1; node_id <= 4; ++node_id) {
    if(!write_binary_value(output, node_id) ||
       !write_binary_values(output, nodes[static_cast<std::size_t>(node_id - 1)].data(), 3U)) {
      return false;
    }
  }
  output << "\n$EndNodes\n";

  output << "$Elements\n11\n";
  const std::array<int, 3> edge_header = {1, 6, 2};
  const std::array<int, 3> face_header = {2, 4, 2};
  const std::array<int, 3> cell_header = {4, 1, 2};
  if(!write_binary_values(output, edge_header.data(), edge_header.size())) {
    return false;
  }
  const std::array<std::array<int, 5>, 6> edge_records = {{
    {1, 11, 11, 1, 2},
    {2, 11, 11, 2, 3},
    {3, 11, 11, 3, 1},
    {4, 11, 11, 1, 4},
    {5, 11, 11, 2, 4},
    {6, 11, 11, 3, 4},
  }};
  for(const auto &record : edge_records) {
    if(!write_binary_values(output, record.data(), record.size())) {
      return false;
    }
  }
  if(!write_binary_values(output, face_header.data(), face_header.size())) {
    return false;
  }
  const std::array<std::array<int, 6>, 4> face_records = {{
    {7, 21, 21, 1, 3, 2},
    {8, 21, 21, 1, 2, 4},
    {9, 21, 21, 2, 3, 4},
    {10, 21, 21, 3, 1, 4},
  }};
  for(const auto &record : face_records) {
    if(!write_binary_values(output, record.data(), record.size())) {
      return false;
    }
  }
  if(!write_binary_values(output, cell_header.data(), cell_header.size())) {
    return false;
  }
  const std::array<int, 7> cell_record = {11, 31, 31, 1, 2, 3, 4};
  if(!write_binary_values(output, cell_record.data(), cell_record.size())) {
    return false;
  }
  output << "\n$EndElements\n";

  return static_cast<bool>(output);
}

template <typename SizeType>
bool write_msh41_binary_fixture_impl(const std::filesystem::path &path)
{
  std::ofstream output(path, std::ios::trunc | std::ios::binary);
  if(!output) {
    std::fprintf(stderr, "mesh_io_smoke: failed to open %s for writing\n", path.string().c_str());
    return false;
  }

  output << "$MeshFormat\n4.1 1 " << sizeof(SizeType) << "\n";
  const int one = 1;
  if(!write_binary_value(output, one)) {
    return false;
  }
  output << "\n$EndMeshFormat\n";
  output << "$PhysicalNames\n3\n";
  output << "1 11 \"feature_edges\"\n";
  output << "2 21 \"wall\"\n";
  output << "3 31 \"volume\"\n";
  output << "$EndPhysicalNames\n";

  output << "$Entities\n";
  const std::array<SizeType, 4> entity_counts = {0U, 1U, 1U, 1U};
  if(!write_binary_values(output, entity_counts.data(), entity_counts.size())) {
    return false;
  }
  const std::array<double, 6> bounds = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
  for(const int entity_tag : {11, 21, 31}) {
    if(!write_binary_value(output, entity_tag) ||
       !write_binary_values(output, bounds.data(), bounds.size())) {
      return false;
    }
    const SizeType physical_count = 1U;
    const SizeType boundary_count = 0U;
    if(!write_binary_value(output, physical_count) ||
       !write_binary_value(output, entity_tag) ||
       !write_binary_value(output, boundary_count)) {
      return false;
    }
  }
  output << "\n$EndEntities\n";

  output << "$Nodes\n";
  const std::array<SizeType, 4> node_header = {1U, 4U, 1U, 4U};
  if(!write_binary_values(output, node_header.data(), node_header.size())) {
    return false;
  }
  const std::array<int, 3> node_block_header = {3, 31, 0};
  const SizeType block_size = 4U;
  const std::array<SizeType, 4> node_tags = {1U, 2U, 3U, 4U};
  const std::array<std::array<double, 3>, 4> nodes = {{
    {0.0, 0.0, 0.0},
    {1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0},
  }};
  if(!write_binary_values(output, node_block_header.data(), node_block_header.size()) ||
     !write_binary_value(output, block_size) ||
     !write_binary_values(output, node_tags.data(), node_tags.size())) {
    return false;
  }
  for(const auto &node : nodes) {
    if(!write_binary_values(output, node.data(), node.size())) {
      return false;
    }
  }
  output << "\n$EndNodes\n";

  output << "$Elements\n";
  const std::array<SizeType, 4> element_header = {3U, 11U, 1U, 11U};
  if(!write_binary_values(output, element_header.data(), element_header.size())) {
    return false;
  }
  const std::array<int, 3> edge_block_header = {1, 11, 1};
  const std::array<int, 3> face_block_header = {2, 21, 2};
  const std::array<int, 3> cell_block_header = {3, 31, 4};
  if(!write_binary_values(output, edge_block_header.data(), edge_block_header.size()) ||
     !write_binary_value(output, static_cast<SizeType>(6U))) {
    return false;
  }
  const std::array<std::array<SizeType, 3>, 6> edge_records = {{
    {1U, 1U, 2U},
    {2U, 2U, 3U},
    {3U, 3U, 1U},
    {4U, 1U, 4U},
    {5U, 2U, 4U},
    {6U, 3U, 4U},
  }};
  for(const auto &record : edge_records) {
    if(!write_binary_values(output, record.data(), record.size())) {
      return false;
    }
  }
  if(!write_binary_values(output, face_block_header.data(), face_block_header.size()) ||
     !write_binary_value(output, static_cast<SizeType>(4U))) {
    return false;
  }
  const std::array<std::array<SizeType, 4>, 4> face_records = {{
    {7U, 1U, 3U, 2U},
    {8U, 1U, 2U, 4U},
    {9U, 2U, 3U, 4U},
    {10U, 3U, 1U, 4U},
  }};
  for(const auto &record : face_records) {
    if(!write_binary_values(output, record.data(), record.size())) {
      return false;
    }
  }
  if(!write_binary_values(output, cell_block_header.data(), cell_block_header.size()) ||
     !write_binary_value(output, static_cast<SizeType>(1U))) {
    return false;
  }
  const std::array<SizeType, 5> cell_record = {11U, 1U, 2U, 3U, 4U};
  if(!write_binary_values(output, cell_record.data(), cell_record.size())) {
    return false;
  }
  output << "\n$EndElements\n";

  return static_cast<bool>(output);
}

bool write_msh41_binary_fixture(const std::filesystem::path &path)
{
  return write_msh41_binary_fixture_impl<std::size_t>(path);
}

bool write_msh41_binary_fixture_32bit(const std::filesystem::path &path)
{
  return write_msh41_binary_fixture_impl<std::uint32_t>(path);
}

template <typename SizeType>
bool write_msh41_partitioned_binary_fixture_impl(const std::filesystem::path &path)
{
  std::ofstream output(path, std::ios::trunc | std::ios::binary);
  if(!output) {
    std::fprintf(stderr, "mesh_io_smoke: failed to open %s for writing\n", path.string().c_str());
    return false;
  }

  output << "$MeshFormat\n4.1 1 " << sizeof(SizeType) << "\n";
  const int one = 1;
  if(!write_binary_value(output, one)) {
    return false;
  }
  output << "\n$EndMeshFormat\n";
  output << "$PhysicalNames\n3\n";
  output << "1 11 \"feature_edges\"\n";
  output << "2 21 \"wall\"\n";
  output << "3 31 \"volume\"\n";
  output << "$EndPhysicalNames\n";

  output << "$PartitionedEntities\n";
  const SizeType partition_count = 2U;
  const SizeType ghost_count = 1U;
  const std::array<int, 2> ghost_record = {999, 2};
  const std::array<SizeType, 4> entity_counts = {0U, 1U, 1U, 1U};
  if(!write_binary_value(output, partition_count) ||
     !write_binary_value(output, ghost_count) ||
     !write_binary_values(output, ghost_record.data(), ghost_record.size()) ||
     !write_binary_values(output, entity_counts.data(), entity_counts.size())) {
    return false;
  }

  const std::array<double, 6> bounds = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
  const SizeType one_count = 1U;
  const SizeType zero_count = 0U;
  if(!write_binary_values(output, std::array<int, 3> {111, 1, 11}.data(), 3U) ||
     !write_binary_value(output, one_count) ||
     !write_binary_value(output, 1) ||
     !write_binary_values(output, bounds.data(), bounds.size()) ||
     !write_binary_value(output, one_count) ||
     !write_binary_value(output, 11) ||
     !write_binary_value(output, zero_count)) {
    return false;
  }
  if(!write_binary_values(output, std::array<int, 3> {121, 2, 21}.data(), 3U) ||
     !write_binary_value(output, one_count) ||
     !write_binary_value(output, 1) ||
     !write_binary_values(output, bounds.data(), bounds.size()) ||
     !write_binary_value(output, one_count) ||
     !write_binary_value(output, 21) ||
     !write_binary_value(output, zero_count)) {
    return false;
  }
  if(!write_binary_values(output, std::array<int, 3> {131, 3, 31}.data(), 3U) ||
     !write_binary_value(output, one_count) ||
     !write_binary_value(output, 2) ||
     !write_binary_values(output, bounds.data(), bounds.size()) ||
     !write_binary_value(output, one_count) ||
     !write_binary_value(output, 31) ||
     !write_binary_value(output, zero_count)) {
    return false;
  }
  output << "\n$EndPartitionedEntities\n";

  output << "$Nodes\n";
  const std::array<SizeType, 4> node_header = {1U, 4U, 1U, 4U};
  const std::array<int, 3> node_block_header = {3, 131, 1};
  const SizeType node_block_size = 4U;
  const std::array<SizeType, 4> node_tags = {1U, 2U, 3U, 4U};
  const std::array<std::array<double, 6>, 4> nodes = {{
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {1.0, 0.0, 0.0, 1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0, 0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0, 0.0, 0.0, 1.0},
  }};
  if(!write_binary_values(output, node_header.data(), node_header.size()) ||
     !write_binary_values(output, node_block_header.data(), node_block_header.size()) ||
     !write_binary_value(output, node_block_size) ||
     !write_binary_values(output, node_tags.data(), node_tags.size())) {
    return false;
  }
  for(const auto &node : nodes) {
    if(!write_binary_values(output, node.data(), node.size())) {
      return false;
    }
  }
  output << "\n$EndNodes\n";

  output << "$Elements\n";
  const std::array<SizeType, 4> element_header = {3U, 11U, 1U, 11U};
  const std::array<int, 3> edge_block_header = {1, 111, 1};
  const std::array<int, 3> face_block_header = {2, 121, 2};
  const std::array<int, 3> cell_block_header = {3, 131, 4};
  if(!write_binary_values(output, element_header.data(), element_header.size()) ||
     !write_binary_values(output, edge_block_header.data(), edge_block_header.size()) ||
     !write_binary_value(output, static_cast<SizeType>(6U))) {
    return false;
  }
  const std::array<std::array<SizeType, 3>, 6> edge_records = {{
    {1U, 1U, 2U},
    {2U, 2U, 3U},
    {3U, 3U, 1U},
    {4U, 1U, 4U},
    {5U, 2U, 4U},
    {6U, 3U, 4U},
  }};
  for(const auto &record : edge_records) {
    if(!write_binary_values(output, record.data(), record.size())) {
      return false;
    }
  }
  if(!write_binary_values(output, face_block_header.data(), face_block_header.size()) ||
     !write_binary_value(output, static_cast<SizeType>(4U))) {
    return false;
  }
  const std::array<std::array<SizeType, 4>, 4> face_records = {{
    {7U, 1U, 3U, 2U},
    {8U, 1U, 2U, 4U},
    {9U, 2U, 3U, 4U},
    {10U, 3U, 1U, 4U},
  }};
  for(const auto &record : face_records) {
    if(!write_binary_values(output, record.data(), record.size())) {
      return false;
    }
  }
  if(!write_binary_values(output, cell_block_header.data(), cell_block_header.size()) ||
     !write_binary_value(output, static_cast<SizeType>(1U))) {
    return false;
  }
  const std::array<SizeType, 5> cell_record = {11U, 1U, 2U, 3U, 4U};
  if(!write_binary_values(output, cell_record.data(), cell_record.size())) {
    return false;
  }
  output << "\n$EndElements\n";

  output << "$Periodic\n";
  const SizeType periodic_link_count = 1U;
  const std::array<int, 3> periodic_header = {2, 121, 121};
  const SizeType affine_count = 0U;
  const SizeType corresponding_node_count = 2U;
  const std::array<SizeType, 4> corresponding_nodes = {1U, 1U, 2U, 2U};
  if(!write_binary_value(output, periodic_link_count) ||
     !write_binary_values(output, periodic_header.data(), periodic_header.size()) ||
     !write_binary_value(output, affine_count) ||
     !write_binary_value(output, corresponding_node_count) ||
     !write_binary_values(output, corresponding_nodes.data(), corresponding_nodes.size())) {
    return false;
  }
  output << "\n$EndPeriodic\n";

  output << "$GhostElements\n";
  const SizeType ghost_element_count = 1U;
  const SizeType ghost_element_tag = 11U;
  const int partition_tag = 7;
  const SizeType ghost_partition_count = 1U;
  const int ghost_partition_tag = 9;
  if(!write_binary_value(output, ghost_element_count) ||
     !write_binary_value(output, ghost_element_tag) ||
     !write_binary_value(output, partition_tag) ||
     !write_binary_value(output, ghost_partition_count) ||
     !write_binary_value(output, ghost_partition_tag)) {
    return false;
  }
  output << "\n$EndGhostElements\n";

  output << "$Parametrizations\n";
  const std::array<SizeType, 2> parametrization_header = {1U, 1U};
  const SizeType curve_param_node_count = 2U;
  const std::array<std::array<double, 4>, 2> curve_param_nodes = {{
    {0.0, 0.0, 0.0, 0.0},
    {1.0, 0.0, 0.0, 1.0},
  }};
  const SizeType surface_param_node_count = 3U;
  const SizeType surface_param_triangle_count = 1U;
  const std::array<std::array<double, 11>, 3> surface_param_nodes = {{
    {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0},
    {1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0},
    {0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0},
  }};
  const std::array<int, 3> surface_triangle = {1, 2, 3};
  if(!write_binary_values(output, parametrization_header.data(), parametrization_header.size()) ||
     !write_binary_value(output, 111) ||
     !write_binary_value(output, curve_param_node_count)) {
    return false;
  }
  for(const auto &node : curve_param_nodes) {
    if(!write_binary_values(output, node.data(), node.size())) {
      return false;
    }
  }
  if(!write_binary_value(output, 121) ||
     !write_binary_value(output, surface_param_node_count) ||
     !write_binary_value(output, surface_param_triangle_count)) {
    return false;
  }
  for(const auto &node : surface_param_nodes) {
    if(!write_binary_values(output, node.data(), node.size())) {
      return false;
    }
  }
  if(!write_binary_values(output, surface_triangle.data(), surface_triangle.size())) {
    return false;
  }
  output << "\n$EndParametrizations\n";

  return static_cast<bool>(output);
}

bool write_msh41_partitioned_binary_fixture(const std::filesystem::path &path)
{
  return write_msh41_partitioned_binary_fixture_impl<std::uint32_t>(path);
}

void append_small_field_card(
  std::ostringstream &stream,
  const char *card,
  const std::vector<std::string> &fields
)
{
  stream << std::left << std::setw(8) << card << std::right;
  for(const auto &field : fields) {
    stream << std::setw(8) << field;
  }
  stream << '\n';
}

void append_large_field_card(
  std::ostringstream &stream,
  const char *card,
  const std::vector<std::string> &fields
)
{
  stream << std::left << std::setw(8) << (std::string(card) + '*');
  for(std::size_t index = 0U; index < 4U && index < fields.size(); ++index) {
    stream << std::setw(16) << fields[index];
  }
  stream << '\n';

  for(std::size_t offset = 4U; offset < fields.size(); offset += 4U) {
    stream << std::left << std::setw(8) << "*";
    for(std::size_t index = 0U; index < 4U && offset + index < fields.size(); ++index) {
      stream << std::setw(16) << fields[offset + index];
    }
    stream << '\n';
  }
}

std::string make_nastran_small_field_sample()
{
  std::ostringstream stream;
  stream << "BEGIN BULK\n";
  append_small_field_card(stream, "MAT1", {"31", "210000.", "", "0.3"});
  append_small_field_card(stream, "PROD", {"11", "31", "1.0"});
  append_small_field_card(stream, "PSHELL", {"21", "31", "1.0"});
  append_small_field_card(stream, "PSOLID", {"31", "31"});
  append_small_field_card(stream, "GRID", {"1", "", "0.0", "0.0", "0.0", "0", "0", "0"});
  append_small_field_card(stream, "GRID", {"2", "", "1.0", "0.0", "0.0", "0", "0", "0"});
  append_small_field_card(stream, "GRID", {"3", "", "0.0", "1.0", "0.0", "0", "0", "0"});
  append_small_field_card(stream, "GRID", {"4", "", "0.0", "0.0", "1.0", "0", "0", "0"});
  append_small_field_card(stream, "CROD", {"1", "11", "1", "2"});
  append_small_field_card(stream, "CROD", {"2", "11", "2", "3"});
  append_small_field_card(stream, "CROD", {"3", "11", "3", "1"});
  append_small_field_card(stream, "CROD", {"4", "11", "1", "4"});
  append_small_field_card(stream, "CROD", {"5", "11", "2", "4"});
  append_small_field_card(stream, "CROD", {"6", "11", "3", "4"});
  append_small_field_card(stream, "CTRIA3", {"7", "21", "1", "3", "2"});
  append_small_field_card(stream, "CTRIA3", {"8", "21", "1", "2", "4"});
  append_small_field_card(stream, "CTRIA3", {"9", "21", "2", "3", "4"});
  append_small_field_card(stream, "CTRIA3", {"10", "21", "3", "1", "4"});
  append_small_field_card(stream, "CTETRA", {"11", "31", "1", "2", "3", "4"});
  stream << "ENDDATA\n";
  return stream.str();
}

std::string make_nastran_small_field_broad_sample()
{
  std::ostringstream stream;
  stream << "BEGIN BULK\n";
  append_small_field_card(stream, "MAT1", {"301", "210000.", "", "0.3"});
  append_small_field_card(stream, "PROD", {"101", "301", "1.0"});
  append_small_field_card(stream, "PSHELL", {"201", "301", "1.0"});
  append_small_field_card(stream, "GRID", {"1", "", "0.0", "0.0", "0.0", "0", "0", "0"});
  append_small_field_card(stream, "GRID", {"2", "", "1.0", "0.0", "0.0", "0", "0", "0"});
  append_small_field_card(stream, "GRID", {"3", "", "1.0", "1.0", "0.0", "0", "0", "0"});
  append_small_field_card(stream, "GRID", {"4", "", "0.0", "1.0", "0.0", "0", "0", "0"});
  append_small_field_card(stream, "GRID", {"5", "", "0.5", "0.5", "1.0", "0", "0", "0"});
  append_small_field_card(stream, "CBAR", {"1", "101", "1", "2", "0.", "0.", "1."});
  append_small_field_card(stream, "CBEAM", {"2", "101", "2", "3", "0.", "0.", "1."});
  append_small_field_card(stream, "CROD", {"3", "101", "3", "4"});
  append_small_field_card(stream, "CQUAD4", {"10", "201", "1", "2", "3", "4", "0."});
  append_small_field_card(stream, "CTRIA3", {"11", "201", "1", "3", "5", "0."});
  stream << "ENDDATA\n";
  return stream.str();
}

std::string make_nastran_large_field_broad_sample()
{
  std::ostringstream stream;
  stream << "BEGIN BULK\n";
  append_large_field_card(stream, "MAT1", {"301", "210000.", "", "0.3"});
  append_large_field_card(stream, "PROD", {"101", "301", "1.0"});
  append_large_field_card(stream, "PSHELL", {"201", "301", "1.0"});
  append_large_field_card(stream, "GRID", {"1", "", "0.0", "0.0", "0.0", "0", "0", "0"});
  append_large_field_card(stream, "GRID", {"2", "", "1.0", "0.0", "0.0", "0", "0", "0"});
  append_large_field_card(stream, "GRID", {"3", "", "1.0", "1.0", "0.0", "0", "0", "0"});
  append_large_field_card(stream, "GRID", {"4", "", "0.0", "1.0", "0.0", "0", "0", "0"});
  append_large_field_card(stream, "GRID", {"5", "", "0.5", "0.5", "1.0", "0", "0", "0"});
  append_large_field_card(stream, "CBAR", {"1", "101", "1", "2", "0.", "0.", "1."});
  append_large_field_card(stream, "CBEAM", {"2", "101", "2", "3", "0.", "0.", "1."});
  append_large_field_card(stream, "CROD", {"3", "101", "3", "4"});
  append_large_field_card(stream, "CQUAD4", {"10", "201", "1", "2", "3", "4", "0."});
  append_large_field_card(stream, "CTRIA3", {"11", "201", "1", "3", "5", "0."});
  stream << "ENDDATA\n";
  return stream.str();
}

const sqmesh::mesh::EntityGroup *find_entity_group(
  const sqmesh::mesh::Domain &domain,
  sqmesh::mesh::EntityOrder order,
  const char *name
)
{
  for(const auto &entity_group : domain.entity_groups()) {
    if(entity_group.order() == order && entity_group.name() == std::string_view(name)) {
      return &entity_group;
    }
  }
  return nullptr;
}

bool expect_optional_real(
  double actual,
  double expected,
  std::string_view message
)
{
  if(!std::isfinite(expected)) {
    return expect(!std::isfinite(actual), std::string(message).c_str());
  }

  return expect(
    std::isfinite(actual) && std::abs(actual - expected) <= 1.0e-12,
    std::string(message).c_str()
  );
}

bool verify_nastran_entity_group_import_info(
  const sqmesh::mesh::EntityGroup &entity_group,
  sqmesh::mesh::NastranEntityGroupSourceCard expected_source_card,
  sqmesh::mesh::NastranPropertyCard expected_property_card,
  sqmesh::mesh::NastranMaterialCard expected_material_card,
  std::uint32_t expected_material_id,
  double expected_rod_area,
  double expected_shell_thickness,
  double expected_youngs_modulus,
  double expected_shear_modulus,
  double expected_poisson_ratio,
  std::string_view label
)
{
  const auto &import_info = entity_group.import_info();
  const auto &nastran = import_info.nastran;
  if(!expect(
       import_info.format == sqmesh::mesh::EntityGroupImportFormat::nastran,
       (std::string(label) + " should publish NASTRAN import metadata").c_str()
     )) {
    return false;
  }
  if(!expect(
       nastran.source_card == expected_source_card,
       (std::string(label) + " should preserve the expected NASTRAN source-card identity").c_str()
     )) {
    return false;
  }
  if(!expect(
       nastran.property_card == expected_property_card,
       (std::string(label) + " should preserve the expected NASTRAN property-card identity").c_str()
     )) {
    return false;
  }
  if(!expect(
       nastran.material_card == expected_material_card,
       (std::string(label) + " should preserve the expected NASTRAN material-card identity").c_str()
     )) {
    return false;
  }
  if(!expect(
       nastran.material_id == expected_material_id,
       (std::string(label) + " should preserve the expected NASTRAN material id").c_str()
     )) {
    return false;
  }
  if(!expect_optional_real(
       nastran.rod_area,
       expected_rod_area,
       std::string(label) + " should preserve the expected NASTRAN rod area"
     )) {
    return false;
  }
  if(!expect_optional_real(
       nastran.shell_thickness,
       expected_shell_thickness,
       std::string(label) + " should preserve the expected NASTRAN shell thickness"
     )) {
    return false;
  }
  if(!expect_optional_real(
       nastran.youngs_modulus,
       expected_youngs_modulus,
       std::string(label) + " should preserve the expected NASTRAN Young's modulus"
     )) {
    return false;
  }
  if(!expect_optional_real(
       nastran.shear_modulus,
       expected_shear_modulus,
       std::string(label) + " should preserve the expected NASTRAN shear modulus"
     )) {
    return false;
  }
  if(!expect_optional_real(
       nastran.poisson_ratio,
       expected_poisson_ratio,
       std::string(label) + " should preserve the expected NASTRAN Poisson ratio"
     )) {
    return false;
  }
  return true;
}

const sqmesh::mesh::EntityGroup *find_entity_group_with_source_entity_tag(
  const sqmesh::mesh::Domain &domain,
  sqmesh::mesh::EntityOrder order,
  const char *name,
  std::uint32_t source_entity_tag
)
{
  for(const auto &entity_group : domain.entity_groups()) {
    if(entity_group.order() == order && entity_group.name() == std::string_view(name) &&
       entity_group.source_entity_tag() == source_entity_tag) {
      return &entity_group;
    }
  }
  return nullptr;
}

bool verify_tetra_domain(
  const sqmesh::mesh::Domain &domain,
  const char *expected_domain_name,
  const char *edge_name,
  std::uint32_t edge_zone_id,
  const char *face_name,
  std::uint32_t face_zone_id,
  const char *cell_name,
  std::uint32_t cell_zone_id
)
{
  const auto summary = domain.summary();
  if(!expect(summary.node_count == 4U, "tetra import should preserve four nodes")) {
    return false;
  }
  if(!expect(summary.edge_count == 6U, "tetra import should preserve six line elements")) {
    return false;
  }
  if(!expect(summary.face_count == 4U, "tetra import should preserve four triangle faces")) {
    return false;
  }
  if(!expect(summary.cell_count == 1U, "tetra import should preserve one tetra cell")) {
    return false;
  }
  if(!expect(domain.entity_group_count() == 4U, "tetra import should build node, edge, face, and cell entity_groups")) {
    return false;
  }
  if(expected_domain_name != nullptr &&
     !expect(
       domain.name() == std::string_view(expected_domain_name),
       "tetra import should preserve the expected domain name"
     )) {
    return false;
  }

  const auto *edge_entity_group = find_entity_group(
    domain,
    sqmesh::mesh::EntityOrder::edge,
    edge_name
  );
  const auto *face_entity_group = find_entity_group(
    domain,
    sqmesh::mesh::EntityOrder::face,
    face_name
  );
  const auto *cell_entity_group = find_entity_group(
    domain,
    sqmesh::mesh::EntityOrder::cell,
    cell_name
  );

  if(!expect(edge_entity_group != nullptr, "tetra import should preserve the line section label")) {
    return false;
  }
  if(!expect(face_entity_group != nullptr, "tetra import should preserve the face section label")) {
    return false;
  }
  if(!expect(cell_entity_group != nullptr, "tetra import should preserve the cell section label")) {
    return false;
  }
  if(!expect(edge_entity_group->zone_id() == edge_zone_id, "tetra line entity_group should preserve the section id")) {
    return false;
  }
  if(!expect(face_entity_group->zone_id() == face_zone_id, "tetra face entity_group should preserve the section id")) {
    return false;
  }
  if(!expect(cell_entity_group->zone_id() == cell_zone_id, "tetra cell entity_group should preserve the section id")) {
    return false;
  }

  const auto first_edge = sqmesh::mesh::EntityRef {edge_entity_group->id(), 0U};
  if(!expect(
       sqmesh::mesh::is_valid(domain.adjacent_face(first_edge, sqmesh::mesh::FaceSide::left)),
       "tetra line elements should map back to adjacent face refs"
     )) {
    return false;
  }
  if(!expect(
       sqmesh::mesh::is_valid(domain.adjacent_face(first_edge, sqmesh::mesh::FaceSide::right)),
       "closed tetra boundary edges should retain their second adjacent face"
     )) {
    return false;
  }

  const auto first_face = sqmesh::mesh::EntityRef {face_entity_group->id(), 0U};
  if(!expect(
       domain.adjacent_cell(first_face, sqmesh::mesh::FaceSide::left).entity_group == cell_entity_group->id(),
       "tetra boundary faces should keep their owning tetra cell adjacency"
     )) {
    return false;
  }
  if(!expect(
       !sqmesh::mesh::is_valid(domain.adjacent_cell(first_face, sqmesh::mesh::FaceSide::right)),
       "tetra boundary faces should not expose a right cell on a single tetra"
     )) {
    return false;
  }

  return true;
}

bool verify_msh_domain(const sqmesh::mesh::Domain &domain)
{
  return verify_tetra_domain(
    domain,
    nullptr,
    "feature_edges",
    11U,
    "wall",
    21U,
    "volume",
    31U
  );
}

bool verify_msh_interface_domain(const sqmesh::mesh::Domain &domain)
{
  const auto summary = domain.summary();
  if(!expect(summary.node_count == 5U, "multi-region MSH import should preserve five nodes")) {
    return false;
  }
  if(!expect(summary.edge_count == 0U, "multi-region MSH interface fixture should not invent line elements")) {
    return false;
  }
  if(!expect(summary.face_count == 7U, "multi-region MSH import should preserve seven unique triangle faces")) {
    return false;
  }
  if(!expect(summary.cell_count == 2U, "multi-region MSH import should preserve two tetra cells")) {
    return false;
  }
  if(!expect(domain.entity_group_count() == 6U, "multi-region MSH import should build node, interface, boundary, and region entity_groups")) {
    return false;
  }

  const auto *interface_entity_group = find_entity_group(
    domain,
    sqmesh::mesh::EntityOrder::face,
    "fluid_solid_interface"
  );
  const auto *fluid_boundary_entity_group = find_entity_group(
    domain,
    sqmesh::mesh::EntityOrder::face,
    "z101__boundary_faces"
  );
  const auto *solid_boundary_entity_group = find_entity_group(
    domain,
    sqmesh::mesh::EntityOrder::face,
    "z202__boundary_faces"
  );
  const auto *fluid_region_entity_group = find_entity_group(
    domain,
    sqmesh::mesh::EntityOrder::cell,
    "fluid_region"
  );
  const auto *solid_region_entity_group = find_entity_group(
    domain,
    sqmesh::mesh::EntityOrder::cell,
    "solid_region"
  );

  if(!expect(interface_entity_group != nullptr, "multi-region MSH import should preserve the explicit interface face label")) {
    return false;
  }
  if(!expect(fluid_boundary_entity_group != nullptr && solid_boundary_entity_group != nullptr,
             "multi-region MSH import should synthesize region-qualified boundary face names")) {
    return false;
  }
  if(!expect(fluid_region_entity_group != nullptr && solid_region_entity_group != nullptr,
             "multi-region MSH import should preserve both region cell entity_group names")) {
    return false;
  }

  if(!expect(interface_entity_group->semantic() == sqmesh::mesh::EntityGroupSemantic::interface,
             "shared cross-region MSH faces should become explicit interface entity_groups")) {
    return false;
  }
  if(!expect(!interface_entity_group->is_boundary(),
             "interface face entity_groups should remain distinct from boundary-only face entity_groups")) {
    return false;
  }
  if(!expect(
       interface_entity_group->zone_id() == 42U &&
         interface_entity_group->primary_region_zone_id() == 101U &&
         interface_entity_group->secondary_region_zone_id() == 202U,
       "interface entity_groups should preserve their source zone id and both adjacent region ids"
     )) {
    return false;
  }
  if(!expect(interface_entity_group->faces().size() == 1U,
             "the interface entity_group should contain the single shared triangle")) {
    return false;
  }

  if(!expect(fluid_boundary_entity_group->semantic() == sqmesh::mesh::EntityGroupSemantic::boundary &&
               fluid_boundary_entity_group->primary_region_zone_id() == 101U &&
               fluid_boundary_entity_group->faces().size() == 3U,
             "fluid boundary faces should remain explicit boundary semantics bound to region 101")) {
    return false;
  }
  if(!expect(solid_boundary_entity_group->semantic() == sqmesh::mesh::EntityGroupSemantic::boundary &&
               solid_boundary_entity_group->primary_region_zone_id() == 202U &&
               solid_boundary_entity_group->faces().size() == 3U,
             "solid boundary faces should remain explicit boundary semantics bound to region 202")) {
    return false;
  }
  if(!expect(fluid_region_entity_group->semantic() == sqmesh::mesh::EntityGroupSemantic::region &&
               fluid_region_entity_group->zone_id() == 101U &&
               fluid_region_entity_group->primary_region_zone_id() == 101U,
             "the first cell entity_group should remain an explicit region entity_group for zone 101")) {
    return false;
  }
  if(!expect(solid_region_entity_group->semantic() == sqmesh::mesh::EntityGroupSemantic::region &&
               solid_region_entity_group->zone_id() == 202U &&
               solid_region_entity_group->primary_region_zone_id() == 202U,
             "the second cell entity_group should remain an explicit region entity_group for zone 202")) {
    return false;
  }

  const auto interface_face = sqmesh::mesh::EntityRef {interface_entity_group->id(), 0U};
  if(!expect(
       sqmesh::mesh::is_interface_entity(domain.face(interface_face).header.flags),
       "the shared MSH interface face should carry the interface entity flag"
     )) {
    return false;
  }
  if(!expect(
       domain.adjacent_cell(interface_face, sqmesh::mesh::FaceSide::left).entity_group ==
         fluid_region_entity_group->id() &&
         domain.adjacent_cell(interface_face, sqmesh::mesh::FaceSide::right).entity_group ==
           solid_region_entity_group->id(),
       "the interface face should keep left/right adjacency to the two region entity_groups"
     )) {
    return false;
  }

  return true;
}

bool verify_msh41_shared_physical_entity_domain(const sqmesh::mesh::Domain &domain)
{
  const auto summary = domain.summary();
  if(!expect(
       summary.node_count == 6U,
       "shared-physical MSH import should preserve six surface vertices"
     )) {
    return false;
  }
  if(!expect(
       summary.edge_count == 0U,
       "shared-physical MSH import should not invent line elements"
     )) {
    return false;
  }
  if(!expect(
       summary.face_count == 2U,
       "shared-physical MSH import should preserve two triangle faces"
     )) {
    return false;
  }
  if(!expect(
       summary.cell_count == 0U,
       "shared-physical MSH import should remain surface-only"
     )) {
    return false;
  }
  if(!expect(
       domain.entity_group_count() == 3U,
       "shared-physical MSH import should build one node entity_group plus one face entity_group per source entity"
     )) {
    return false;
  }

  const auto *first_wall_entity_group = find_entity_group_with_source_entity_tag(
    domain,
    sqmesh::mesh::EntityOrder::face,
    "wall",
    101U
  );
  const auto *second_wall_entity_group = find_entity_group_with_source_entity_tag(
    domain,
    sqmesh::mesh::EntityOrder::face,
    "wall",
    202U
  );
  if(!expect(
       first_wall_entity_group != nullptr && second_wall_entity_group != nullptr,
       "shared-physical MSH import should keep distinct face entity_groups for each source entity tag"
     )) {
    return false;
  }
  if(!expect(
       first_wall_entity_group->id() != second_wall_entity_group->id(),
       "shared-physical MSH import should not merge distinct source entities into one entity_group"
     )) {
    return false;
  }
  if(!expect(
       first_wall_entity_group->zone_id() == 21U && second_wall_entity_group->zone_id() == 21U,
       "shared-physical MSH face entity_groups should still preserve the shared physical tag"
     )) {
    return false;
  }
  if(!expect(
       first_wall_entity_group->semantic() == sqmesh::mesh::EntityGroupSemantic::boundary &&
         second_wall_entity_group->semantic() == sqmesh::mesh::EntityGroupSemantic::boundary &&
         first_wall_entity_group->is_boundary() && second_wall_entity_group->is_boundary(),
       "shared-physical MSH surface entities should remain explicit boundary face entity_groups"
     )) {
    return false;
  }
  if(!expect(
       first_wall_entity_group->faces().size() == 1U && second_wall_entity_group->faces().size() == 1U,
       "shared-physical MSH face entity_groups should each contain their own element block"
     )) {
    return false;
  }

  const auto first_face = sqmesh::mesh::EntityRef {first_wall_entity_group->id(), 0U};
  const auto second_face = sqmesh::mesh::EntityRef {second_wall_entity_group->id(), 0U};
  if(!expect(
       !sqmesh::mesh::is_valid(domain.adjacent_cell(first_face, sqmesh::mesh::FaceSide::left)) &&
         !sqmesh::mesh::is_valid(domain.adjacent_cell(second_face, sqmesh::mesh::FaceSide::left)),
       "shared-physical MSH surface faces should stay surface-only after import"
     )) {
    return false;
  }

  return true;
}

bool verify_obj_domain(const sqmesh::mesh::Domain &domain)
{
  const auto summary = domain.summary();
  if(!expect(summary.node_count == 3U, "OBJ import should preserve three vertices")) {
    return false;
  }
  if(!expect(summary.edge_count == 3U, "OBJ import should preserve three line records")) {
    return false;
  }
  if(!expect(summary.face_count == 1U, "OBJ import should preserve one triangular face")) {
    return false;
  }
  if(!expect(summary.cell_count == 0U, "OBJ import should remain surface-only")) {
    return false;
  }
  if(!expect(domain.entity_group_count() == 3U, "OBJ import should build node, edge, and face entity_groups")) {
    return false;
  }
  if(!expect(domain.name() == std::string_view("tri_patch"), "OBJ import should keep the object name")) {
    return false;
  }

  const auto *edge_entity_group = find_entity_group(
    domain,
    sqmesh::mesh::EntityOrder::edge,
    "feature_loop"
  );
  const auto *face_entity_group = find_entity_group(
    domain,
    sqmesh::mesh::EntityOrder::face,
    "surface_patch"
  );

  if(!expect(edge_entity_group != nullptr, "OBJ import should preserve the edge group label")) {
    return false;
  }
  if(!expect(face_entity_group != nullptr, "OBJ import should preserve the face group label")) {
    return false;
  }
  if(!expect(edge_entity_group->is_boundary(), "OBJ open polyline edges should remain boundary edges")) {
    return false;
  }
  if(!expect(face_entity_group->is_boundary(), "OBJ surface faces should remain boundary faces")) {
    return false;
  }

  const auto first_edge = sqmesh::mesh::EntityRef {edge_entity_group->id(), 0U};
  if(!expect(
       sqmesh::mesh::is_valid(domain.adjacent_face(first_edge, sqmesh::mesh::FaceSide::left)),
       "OBJ line edges should retain their adjacent face mapping"
     )) {
    return false;
  }
  if(!expect(
       !sqmesh::mesh::is_valid(domain.adjacent_face(first_edge, sqmesh::mesh::FaceSide::right)),
       "OBJ boundary edges should have no right-adjacent face"
     )) {
    return false;
  }

  return true;
}

bool verify_grouped_obj_domain(const sqmesh::mesh::Domain &domain)
{
  const auto summary = domain.summary();
  if(!expect(summary.node_count == 5U, "grouped OBJ import should preserve five vertices")) {
    return false;
  }
  if(!expect(summary.edge_count == 4U, "grouped OBJ import should expand one polyline into four edges")) {
    return false;
  }
  if(!expect(summary.face_count == 4U, "grouped OBJ import should triangulate grouped polygon faces")) {
    return false;
  }
  if(!expect(summary.cell_count == 0U, "grouped OBJ import should remain surface-only")) {
    return false;
  }
  if(!expect(domain.entity_group_count() == 5U, "grouped OBJ import should build one edge entity_group plus three face entity_groups")) {
    return false;
  }
  if(!expect(
       domain.name() == std::string_view("grouped_patch"),
       "grouped OBJ import should keep the first object name as the domain name"
     )) {
    return false;
  }

  const auto *edge_entity_group = find_entity_group(
    domain,
    sqmesh::mesh::EntityOrder::edge,
    "boundary_loop"
  );
  const auto *panel_entity_group = find_entity_group(
    domain,
    sqmesh::mesh::EntityOrder::face,
    "panel_A seam"
  );
  const auto *cap_entity_group = find_entity_group(
    domain,
    sqmesh::mesh::EntityOrder::face,
    "cap"
  );
  const auto *object_entity_group = find_entity_group(
    domain,
    sqmesh::mesh::EntityOrder::face,
    "loose_object"
  );

  if(!expect(edge_entity_group != nullptr, "grouped OBJ import should preserve the polyline group label")) {
    return false;
  }
  if(!expect(panel_entity_group != nullptr, "grouped OBJ import should preserve the multi-token face group label")) {
    return false;
  }
  if(!expect(cap_entity_group != nullptr, "grouped OBJ import should preserve the cap face group label")) {
    return false;
  }
  if(!expect(object_entity_group != nullptr, "grouped OBJ import should fall back to the object name when no group is active")) {
    return false;
  }
  if(!expect(panel_entity_group->faces().size() == 2U, "grouped OBJ import should triangulate the quad into two grouped faces")) {
    return false;
  }
  if(!expect(cap_entity_group->faces().size() == 1U, "grouped OBJ import should keep the cap triangle in its own group")) {
    return false;
  }
  if(!expect(object_entity_group->faces().size() == 1U, "grouped OBJ import should keep the object-fallback triangle in its own entity_group")) {
    return false;
  }
  if(!expect(edge_entity_group->is_boundary(), "grouped OBJ polyline edges should remain boundary edges")) {
    return false;
  }
  if(!expect(
       panel_entity_group->is_boundary() && cap_entity_group->is_boundary() && object_entity_group->is_boundary(),
       "grouped OBJ face entity_groups should remain boundary faces"
     )) {
    return false;
  }

  const auto first_edge = sqmesh::mesh::EntityRef {edge_entity_group->id(), 0U};
  if(!expect(
       sqmesh::mesh::is_valid(domain.adjacent_face(first_edge, sqmesh::mesh::FaceSide::left)),
       "grouped OBJ polyline edges should retain their adjacent face mapping"
     )) {
    return false;
  }
  if(!expect(
       !sqmesh::mesh::is_valid(domain.adjacent_face(first_edge, sqmesh::mesh::FaceSide::right)),
       "grouped OBJ boundary polyline edges should not expose a right-adjacent face"
     )) {
    return false;
  }

  return true;
}

using Point2 = std::array<double, 2>;

constexpr std::array<Point2, 6> kExpectedConcaveObjPolygon = {{
  {3.0, 1.0},
  {1.0, 1.0},
  {1.0, 3.0},
  {0.0, 3.0},
  {0.0, 0.0},
  {3.0, 0.0},
}};

template <std::size_t N>
[[nodiscard]] double polygon_area_2d(const std::array<Point2, N> &polygon)
{
  double area = 0.0;
  for(std::size_t index = 0U; index < polygon.size(); ++index) {
    const auto &current = polygon[index];
    const auto &next = polygon[(index + 1U) % polygon.size()];
    area += current[0] * next[1] - next[0] * current[1];
  }
  return std::abs(0.5 * area);
}

[[nodiscard]] double triangle_area_2d(const std::array<Point2, 3> &triangle)
{
  return std::abs(
    0.5 *
    ((triangle[1][0] - triangle[0][0]) * (triangle[2][1] - triangle[0][1]) -
     (triangle[2][0] - triangle[0][0]) * (triangle[1][1] - triangle[0][1]))
  );
}

template <std::size_t N>
[[nodiscard]] bool point_in_polygon_2d(
  const std::array<Point2, N> &polygon,
  const Point2 &point
)
{
  bool inside = false;
  for(std::size_t index = 0U, previous = polygon.size() - 1U;
      index < polygon.size();
      previous = index++) {
    const auto &current = polygon[index];
    const auto &prior = polygon[previous];
    const bool crosses =
      ((current[1] > point[1]) != (prior[1] > point[1])) &&
      (point[0] <
       (prior[0] - current[0]) * (point[1] - current[1]) / (prior[1] - current[1]) +
         current[0]);
    if(crosses) {
      inside = !inside;
    }
  }
  return inside;
}

bool verify_concave_obj_domain(const sqmesh::mesh::Domain &domain)
{
  const auto summary = domain.summary();
  if(!expect(summary.node_count == 6U, "concave OBJ import should preserve six vertices")) {
    return false;
  }
  if(!expect(summary.edge_count == 0U, "concave OBJ import should remain face-only when no polylines are present")) {
    return false;
  }
  if(!expect(summary.face_count == 4U, "concave OBJ import should triangulate one six-vertex concave polygon into four faces")) {
    return false;
  }
  if(!expect(summary.cell_count == 0U, "concave OBJ import should remain surface-only")) {
    return false;
  }
  if(!expect(domain.entity_group_count() == 2U, "concave OBJ import should build node and face entity_groups")) {
    return false;
  }
  if(!expect(domain.name() == std::string_view("concave_patch"), "concave OBJ import should keep the object name")) {
    return false;
  }

  const auto *face_entity_group = find_entity_group(
    domain,
    sqmesh::mesh::EntityOrder::face,
    "concave_face"
  );
  if(!expect(face_entity_group != nullptr, "concave OBJ import should preserve the face group label")) {
    return false;
  }
  if(!expect(face_entity_group->is_boundary(), "concave OBJ face entity_groups should remain boundary faces")) {
    return false;
  }
  if(!expect(face_entity_group->faces().size() == 4U, "concave OBJ import should keep all ear-clipped triangles in one face entity_group")) {
    return false;
  }

  std::array<bool, kExpectedConcaveObjPolygon.size()> seen_vertices = {
    false,
    false,
    false,
    false,
    false,
    false,
  };
  double total_triangle_area = 0.0;
  for(std::uint32_t face_index = 0U;
      face_index < static_cast<std::uint32_t>(face_entity_group->faces().size());
      ++face_index) {
    std::array<Point2, 3> triangle {};
    std::size_t node_count = 0U;
    for(const auto node_ref : domain.face_nodes({face_entity_group->id(), face_index})) {
      if(!expect(node_ref.index < seen_vertices.size(), "concave OBJ triangles should reference only imported vertices")) {
        return false;
      }
      seen_vertices[node_ref.index] = true;
      const auto &node = domain.node(node_ref);
      if(!expect(
           std::abs(node.coordinates[2]) <= 1.0e-12,
           "concave OBJ triangles should stay on the original planar z=0 slice"
         )) {
        return false;
      }
      triangle[node_count] = {node.coordinates[0], node.coordinates[1]};
      ++node_count;
    }
    if(!expect(node_count == 3U, "concave OBJ triangles should stay first-order triangles")) {
      return false;
    }

    total_triangle_area += triangle_area_2d(triangle);
    const Point2 centroid = {
      (triangle[0][0] + triangle[1][0] + triangle[2][0]) / 3.0,
      (triangle[0][1] + triangle[1][1] + triangle[2][1]) / 3.0,
    };
    if(!expect(
         point_in_polygon_2d(kExpectedConcaveObjPolygon, centroid),
         "concave OBJ ear clipping should keep each triangle centroid inside the source polygon"
       )) {
      return false;
    }
  }

  if(!expect(
       std::count(seen_vertices.begin(), seen_vertices.end(), true) == 6,
       "concave OBJ triangulation should cover all polygon vertices"
     )) {
    return false;
  }
  if(!expect(
       std::abs(total_triangle_area - polygon_area_2d(kExpectedConcaveObjPolygon)) <= 1.0e-12,
       "concave OBJ ear clipping should preserve the source polygon area"
     )) {
    return false;
  }

  return true;
}

bool verify_obj_attribute_review(
  const sqmesh::mesh::testing::ObjImportReview &review
)
{
  if(!expect(
       review.name == std::string("obj_attr_patch"),
       "OBJ parser review should preserve the object name"
     )) {
    return false;
  }
  if(!expect(
       review.texture_coordinates.size() == 5U,
       "OBJ parser review should preserve five texture-coordinate records"
     )) {
    return false;
  }
  if(!expect(
       review.normals.size() == 2U,
       "OBJ parser review should preserve two normal records"
     )) {
    return false;
  }
  if(!expect(
       review.parameter_vertices.size() == 3U,
       "OBJ parser review should preserve three parameter-vertex records"
     )) {
    return false;
  }
  if(!expect(
       review.material_libraries.size() == 1U,
       "OBJ parser review should preserve one mtllib record"
     )) {
    return false;
  }
  if(!expect(
       review.material_libraries.front().entries.size() == 2U &&
         review.material_libraries.front().entries[0] == std::string("base.mtl") &&
         review.material_libraries.front().entries[1] == std::string("accent.mtl"),
       "OBJ parser review should preserve the mtllib entries in order"
     )) {
    return false;
  }
  if(!expect(
       review.texture_coordinates[0].component_count == 1U &&
         review.texture_coordinates[1].component_count == 2U &&
         review.texture_coordinates[2].component_count == 3U,
       "OBJ parser review should preserve bounded texture-coordinate arity"
     )) {
    return false;
  }
  if(!expect(
       review.parameter_vertices[0].component_count == 1U &&
         review.parameter_vertices[1].component_count == 2U &&
         review.parameter_vertices[2].component_count == 3U,
       "OBJ parser review should preserve bounded parameter-vertex arity"
     )) {
    return false;
  }
  if(!expect(
       review.edges.size() == 2U,
       "OBJ parser review should preserve both polyline segments"
     )) {
    return false;
  }
  if(!expect(
       review.edges[0].material_name == std::string("wire") &&
         review.edges[1].material_name == std::string("wire"),
       "OBJ parser review should preserve usemtl on parsed line segments"
     )) {
    return false;
  }
  if(!expect(
       review.faces.size() == 4U,
       "OBJ parser review should preserve four parsed face records"
     )) {
    return false;
  }

  const auto &vertex_only_face = review.faces[0];
  if(!expect(
       vertex_only_face.material_name == std::string("plain"),
       "OBJ parser review should preserve usemtl on vertex-only faces"
     )) {
    return false;
  }
  if(!expect(
       vertex_only_face.corners[0].vertex_index == 0U &&
         !vertex_only_face.corners[0].texture_index.present &&
         !vertex_only_face.corners[0].normal_index.present &&
         vertex_only_face.corners[1].vertex_index == 1U &&
         vertex_only_face.corners[2].vertex_index == 2U,
       "OBJ parser review should preserve plain `v` face tokens without inventing texture or normal bindings"
     )) {
    return false;
  }

  const auto &texture_face = review.faces[1];
  if(!expect(
       texture_face.material_name == std::string("tex_only"),
       "OBJ parser review should preserve usemtl on `v/vt` faces"
     )) {
    return false;
  }
  if(!expect(
       texture_face.corners[0].vertex_index == 0U &&
         texture_face.corners[0].texture_index.present &&
         texture_face.corners[0].texture_index.value == 0U &&
         !texture_face.corners[0].normal_index.present &&
         texture_face.corners[1].vertex_index == 2U &&
         texture_face.corners[1].texture_index.present &&
         texture_face.corners[1].texture_index.value == 1U &&
         texture_face.corners[2].vertex_index == 3U &&
         texture_face.corners[2].texture_index.present &&
         texture_face.corners[2].texture_index.value == 3U,
       "OBJ parser review should preserve `v/vt` face-token texture indices"
     )) {
    return false;
  }

  const auto &normal_face = review.faces[2];
  if(!expect(
       normal_face.material_name == std::string("normal_only"),
       "OBJ parser review should preserve usemtl on `v//vn` faces"
     )) {
    return false;
  }
  if(!expect(
       normal_face.corners[0].vertex_index == 0U &&
         !normal_face.corners[0].texture_index.present &&
         normal_face.corners[0].normal_index.present &&
         normal_face.corners[0].normal_index.value == 0U &&
         normal_face.corners[1].vertex_index == 3U &&
         normal_face.corners[1].normal_index.present &&
         normal_face.corners[2].vertex_index == 4U &&
         normal_face.corners[2].normal_index.present,
       "OBJ parser review should preserve `v//vn` face-token normal indices"
     )) {
    return false;
  }

  const auto &full_face = review.faces[3];
  if(!expect(
       full_face.material_name == std::string("full_negative"),
       "OBJ parser review should preserve usemtl on `v/vt/vn` faces"
     )) {
    return false;
  }
  if(!expect(
       full_face.corners[0].vertex_index == 4U &&
         full_face.corners[0].texture_index.present &&
         full_face.corners[0].texture_index.value == 2U &&
         full_face.corners[0].normal_index.present &&
         full_face.corners[0].normal_index.value == 1U &&
         full_face.corners[1].vertex_index == 5U &&
         full_face.corners[1].texture_index.present &&
         full_face.corners[1].texture_index.value == 3U &&
         full_face.corners[1].normal_index.present &&
         full_face.corners[1].normal_index.value == 1U &&
         full_face.corners[2].vertex_index == 6U &&
         full_face.corners[2].texture_index.present &&
         full_face.corners[2].texture_index.value == 4U &&
         full_face.corners[2].normal_index.present &&
         full_face.corners[2].normal_index.value == 1U,
       "OBJ parser review should resolve negative `v/vt/vn` face-token indices safely"
     )) {
    return false;
  }

  return true;
}

bool verify_obj_attribute_export_truthfulness(const std::string &text)
{
  if(!expect(
       text.find("o obj_attr_patch") != std::string::npos,
       "OBJ attribute-channel export should keep the object name on the bounded subset"
     )) {
    return false;
  }
  if(!expect(
       text.find("g feature_edges") != std::string::npos &&
         text.find("g faces") != std::string::npos,
       "OBJ attribute-channel export should keep grouped edge and face labels"
     )) {
    return false;
  }
  if(!expect(
       text.find("vt ") == std::string::npos &&
         text.find("vn ") == std::string::npos &&
         text.find("vp ") == std::string::npos &&
         text.find("mtllib") == std::string::npos &&
         text.find("usemtl") == std::string::npos,
       "OBJ export should stay honest about not re-emitting parser-only attribute/material metadata"
     )) {
    return false;
  }

  std::istringstream stream(text);
  std::string line;
  while(std::getline(stream, line)) {
    if((line.rfind("f ", 0U) == 0U || line.rfind("l ", 0U) == 0U) &&
       line.find('/') != std::string::npos) {
      return expect(
        false,
        "OBJ export should emit vertex-only line and face references on the supported subset"
      );
    }
  }

  return true;
}

bool verify_nastran_domain(const sqmesh::mesh::Domain &domain)
{
  if(!verify_tetra_domain(
       domain,
       nullptr,
       "line_pid_11",
       11U,
       "surface_pid_21",
       21U,
       "volume_pid_31",
       31U
     )) {
    return false;
  }

  const auto *edge_entity_group = find_entity_group(
    domain,
    sqmesh::mesh::EntityOrder::edge,
    "line_pid_11"
  );
  const auto *face_entity_group = find_entity_group(
    domain,
    sqmesh::mesh::EntityOrder::face,
    "surface_pid_21"
  );
  const auto *cell_entity_group = find_entity_group(
    domain,
    sqmesh::mesh::EntityOrder::cell,
    "volume_pid_31"
  );
  if(!expect(
       edge_entity_group != nullptr && face_entity_group != nullptr && cell_entity_group != nullptr,
       "bounded NASTRAN tetra fixture should expose deterministic line, face, and cell entity_groups"
     )) {
    return false;
  }
  if(!verify_nastran_entity_group_import_info(
       *edge_entity_group,
       sqmesh::mesh::NastranEntityGroupSourceCard::crod,
       sqmesh::mesh::NastranPropertyCard::prod,
       sqmesh::mesh::NastranMaterialCard::mat1,
       31U,
       1.0,
       std::numeric_limits<double>::quiet_NaN(),
       210000.0,
       std::numeric_limits<double>::quiet_NaN(),
       0.3,
       "bounded NASTRAN line entity_group"
     )) {
    return false;
  }
  if(!verify_nastran_entity_group_import_info(
       *face_entity_group,
       sqmesh::mesh::NastranEntityGroupSourceCard::ctria3,
       sqmesh::mesh::NastranPropertyCard::pshell,
       sqmesh::mesh::NastranMaterialCard::mat1,
       31U,
       std::numeric_limits<double>::quiet_NaN(),
       1.0,
       210000.0,
       std::numeric_limits<double>::quiet_NaN(),
       0.3,
       "bounded NASTRAN face entity_group"
     )) {
    return false;
  }
  if(!verify_nastran_entity_group_import_info(
       *cell_entity_group,
       sqmesh::mesh::NastranEntityGroupSourceCard::ctetra,
       sqmesh::mesh::NastranPropertyCard::psolid,
       sqmesh::mesh::NastranMaterialCard::mat1,
       31U,
       std::numeric_limits<double>::quiet_NaN(),
       std::numeric_limits<double>::quiet_NaN(),
       210000.0,
       std::numeric_limits<double>::quiet_NaN(),
       0.3,
       "bounded NASTRAN cell entity_group"
     )) {
    return false;
  }
  return true;
}

bool verify_broadened_nastran_surface_domain(
  const sqmesh::mesh::Domain &domain,
  sqmesh::mesh::NastranEntityGroupSourceCard expected_edge_source_card =
    sqmesh::mesh::NastranEntityGroupSourceCard::mixed
)
{
  const auto summary = domain.summary();
  if(!expect(summary.node_count == 5U, "broadened NASTRAN import should preserve five GRID nodes")) {
    return false;
  }
  if(!expect(summary.edge_count == 3U, "broadened NASTRAN import should preserve three line elements")) {
    return false;
  }
  if(!expect(summary.face_count == 3U, "broadened NASTRAN import should triangulate one CQUAD4 plus one CTRIA3 into three faces")) {
    return false;
  }
  if(!expect(summary.cell_count == 0U, "broadened NASTRAN surface fixture should remain surface-only")) {
    return false;
  }
  if(!expect(
       domain.entity_group_count() == 4U,
       "broadened NASTRAN surface fixture should keep distinct triangle and CQUAD4 face entity_groups"
     )) {
    return false;
  }

  const auto *edge_entity_group = find_entity_group(
    domain,
    sqmesh::mesh::EntityOrder::edge,
    "line_pid_101"
  );
  const auto *face_entity_group = find_entity_group(
    domain,
    sqmesh::mesh::EntityOrder::face,
    "surface_pid_201"
  );
  const auto *quad_face_entity_group = find_entity_group(
    domain,
    sqmesh::mesh::EntityOrder::face,
    "cquad4_pid_201"
  );
  if(!expect(edge_entity_group != nullptr, "broadened NASTRAN import should map line PIDs onto deterministic entity_group names")) {
    return false;
  }
  if(!expect(
       face_entity_group != nullptr && quad_face_entity_group != nullptr,
       "broadened NASTRAN import should keep CTRIA3 and CQUAD4 shell families in separate deterministic face entity_groups"
     )) {
    return false;
  }
  if(!verify_nastran_entity_group_import_info(
       *edge_entity_group,
       expected_edge_source_card,
       sqmesh::mesh::NastranPropertyCard::prod,
       sqmesh::mesh::NastranMaterialCard::mat1,
       301U,
       1.0,
       std::numeric_limits<double>::quiet_NaN(),
       210000.0,
       std::numeric_limits<double>::quiet_NaN(),
       0.3,
       "broadened NASTRAN line entity_group"
     )) {
    return false;
  }
  if(!verify_nastran_entity_group_import_info(
       *face_entity_group,
       sqmesh::mesh::NastranEntityGroupSourceCard::ctria3,
       sqmesh::mesh::NastranPropertyCard::pshell,
       sqmesh::mesh::NastranMaterialCard::mat1,
       301U,
       std::numeric_limits<double>::quiet_NaN(),
       1.0,
       210000.0,
       std::numeric_limits<double>::quiet_NaN(),
       0.3,
       "broadened NASTRAN CTRIA3 face entity_group"
     )) {
    return false;
  }
  if(!verify_nastran_entity_group_import_info(
       *quad_face_entity_group,
       sqmesh::mesh::NastranEntityGroupSourceCard::cquad4,
       sqmesh::mesh::NastranPropertyCard::pshell,
       sqmesh::mesh::NastranMaterialCard::mat1,
       301U,
       std::numeric_limits<double>::quiet_NaN(),
       1.0,
       210000.0,
       std::numeric_limits<double>::quiet_NaN(),
       0.3,
       "broadened NASTRAN CQUAD4 face entity_group"
     )) {
    return false;
  }
  if(!expect(edge_entity_group->zone_id() == 101U, "broadened NASTRAN edge entity_group should preserve the PROD PID")) {
    return false;
  }
  if(!expect(
       face_entity_group->zone_id() == 201U && quad_face_entity_group->zone_id() == 201U,
       "broadened NASTRAN face entity_groups should preserve the PSHELL PID across both shell families"
     )) {
    return false;
  }
  if(!expect(edge_entity_group->edges().size() == 3U, "broadened NASTRAN import should keep CBAR, CBEAM, and CROD records in one line entity_group")) {
    return false;
  }
  if(!expect(
       face_entity_group->faces().size() == 1U &&
         quad_face_entity_group->faces().size() == 2U,
       "broadened NASTRAN import should keep CTRIA3 faces separate from the paired triangles that came from CQUAD4"
     )) {
    return false;
  }
  const auto triangle_face_ref = sqmesh::mesh::EntityRef {face_entity_group->id(), 0U};
  const auto quad_face_ref_0 = sqmesh::mesh::EntityRef {quad_face_entity_group->id(), 0U};
  const auto quad_face_ref_1 = sqmesh::mesh::EntityRef {quad_face_entity_group->id(), 1U};
  if(!expect(
       domain.face_source_entity_tag(triangle_face_ref) != sqmesh::mesh::invalid_index &&
         domain.face_source_entity_tag(quad_face_ref_0) != sqmesh::mesh::invalid_index &&
         domain.face_source_entity_tag(quad_face_ref_0) ==
           domain.face_source_entity_tag(quad_face_ref_1) &&
         domain.face_source_entity_tag(triangle_face_ref) !=
           domain.face_source_entity_tag(quad_face_ref_0),
       "broadened NASTRAN import should preserve per-face source tags so exported CQUAD4 pairing stays evidence-based"
     )) {
    return false;
  }
  if(!expect(edge_entity_group->is_boundary(), "broadened NASTRAN line entity_group should remain boundary-only for the surface fixture")) {
    return false;
  }
  if(!expect(
       face_entity_group->is_boundary() && quad_face_entity_group->is_boundary(),
       "broadened NASTRAN face entity_groups should remain boundary faces"
     )) {
    return false;
  }

  const auto first_edge = sqmesh::mesh::EntityRef {edge_entity_group->id(), 0U};
  if(!expect(
       sqmesh::mesh::is_valid(domain.adjacent_face(first_edge, sqmesh::mesh::FaceSide::left)),
       "broadened NASTRAN surface edges should retain their adjacent face mapping"
     )) {
    return false;
  }

  return true;
}

sqmesh::mesh::Domain make_msh_physical_name_conflict_domain()
{
  sqmesh::mesh::Domain domain("msh_physical_name_conflict");
  const auto node_entity_group =
    domain.create_entity_group({sqmesh::mesh::EntityOrder::node, "nodes"});

  sqmesh::mesh::EntityGroupDefinition left_definition;
  left_definition.order = sqmesh::mesh::EntityOrder::face;
  left_definition.name = "left_wall";
  left_definition.zone_id = 21U;
  left_definition.boundary = true;
  left_definition.default_kind = sqmesh::mesh::EntityKind::face_triangle;
  left_definition.semantic = sqmesh::mesh::EntityGroupSemantic::boundary;
  const auto left_entity_group = domain.create_entity_group(std::move(left_definition));

  sqmesh::mesh::EntityGroupDefinition right_definition;
  right_definition.order = sqmesh::mesh::EntityOrder::face;
  right_definition.name = "right_wall";
  right_definition.zone_id = 21U;
  right_definition.boundary = true;
  right_definition.default_kind = sqmesh::mesh::EntityKind::face_triangle;
  right_definition.semantic = sqmesh::mesh::EntityGroupSemantic::boundary;
  const auto right_entity_group = domain.create_entity_group(std::move(right_definition));

  const auto n0 = domain.add_node(node_entity_group, {0.0, 0.0, 0.0});
  const auto n1 = domain.add_node(node_entity_group, {1.0, 0.0, 0.0});
  const auto n2 = domain.add_node(node_entity_group, {0.0, 1.0, 0.0});
  const auto n3 = domain.add_node(node_entity_group, {2.0, 0.0, 0.0});
  const auto n4 = domain.add_node(node_entity_group, {2.0, 1.0, 0.0});

  static_cast<void>(domain.add_triangle_face(left_entity_group, {n0, n1, n2}));
  static_cast<void>(domain.add_triangle_face(right_entity_group, {n1, n3, n4}));
  return domain;
}

sqmesh::mesh::Domain make_nastran_unpaired_cquad4_domain()
{
  sqmesh::mesh::Domain domain("nastran_unpaired_cquad4");
  const auto node_entity_group =
    domain.create_entity_group({sqmesh::mesh::EntityOrder::node, "nodes"});

  sqmesh::mesh::EntityGroupDefinition face_definition;
  face_definition.order = sqmesh::mesh::EntityOrder::face;
  face_definition.name = "cquad4_pid_201";
  face_definition.zone_id = 201U;
  face_definition.boundary = true;
  face_definition.default_kind = sqmesh::mesh::EntityKind::face_triangle;
  face_definition.semantic = sqmesh::mesh::EntityGroupSemantic::boundary;
  face_definition.import_info.format = sqmesh::mesh::EntityGroupImportFormat::nastran;
  face_definition.import_info.nastran.source_card =
    sqmesh::mesh::NastranEntityGroupSourceCard::cquad4;
  const auto face_entity_group = domain.create_entity_group(std::move(face_definition));

  const auto n0 = domain.add_node(node_entity_group, {0.0, 0.0, 0.0});
  const auto n1 = domain.add_node(node_entity_group, {1.0, 0.0, 0.0});
  const auto n2 = domain.add_node(node_entity_group, {1.0, 1.0, 0.0});
  const auto n3 = domain.add_node(node_entity_group, {0.0, 1.0, 0.0});

  static_cast<void>(domain.add_triangle_face(face_entity_group, {n0, n1, n2}));
  static_cast<void>(domain.add_triangle_face(face_entity_group, {n0, n2, n3}));
  return domain;
}

#ifdef sqmesh_TEST_CGNS_ENABLED

sqmesh::mesh::EntityGroupImportInfo make_cgns_entity_group_import_info(
  std::uint32_t zone_index,
  std::string_view zone_name,
  std::string_view local_name,
  std::int32_t bc_type_value = -1,
  std::uint32_t base_index = 1U,
  std::string_view base_name = "sqmesh_base"
)
{
  sqmesh::mesh::EntityGroupImportInfo import_info;
  import_info.format = sqmesh::mesh::EntityGroupImportFormat::cgns;
  import_info.cgns.base_index = base_index;
  import_info.cgns.base_name = std::string(base_name);
  import_info.cgns.zone_index = zone_index;
  import_info.cgns.zone_name = std::string(zone_name);
  import_info.cgns.local_name = std::string(local_name);
  import_info.cgns.bc_type_value = bc_type_value;
  return import_info;
}

sqmesh::mesh::Domain make_cgns_multi_zone_interface_domain(
  bool duplicate_secondary_interface_nodes
)
{
  sqmesh::mesh::Domain domain(
    duplicate_secondary_interface_nodes
      ? "manual_cgns_nonshared_interface"
      : "manual_cgns_shared_interface"
  );
  const auto node_entity_group =
    domain.create_entity_group({sqmesh::mesh::EntityOrder::node, "nodes"});

  sqmesh::mesh::EntityGroupDefinition zone_a_volume_definition;
  zone_a_volume_definition.order = sqmesh::mesh::EntityOrder::cell;
  zone_a_volume_definition.name = "z1_zone_a__volume";
  zone_a_volume_definition.zone_id = 31U;
  zone_a_volume_definition.default_kind = sqmesh::mesh::EntityKind::cell_tetra;
  zone_a_volume_definition.semantic = sqmesh::mesh::EntityGroupSemantic::region;
  zone_a_volume_definition.primary_region_zone_id = 1U;
  zone_a_volume_definition.import_info =
    make_cgns_entity_group_import_info(1U, "zone_a", "volume");
  const auto zone_a_volume_entity_group =
    domain.create_entity_group(std::move(zone_a_volume_definition));

  sqmesh::mesh::EntityGroupDefinition zone_a_boundary_definition;
  zone_a_boundary_definition.order = sqmesh::mesh::EntityOrder::face;
  zone_a_boundary_definition.name = "z1_zone_a__wall";
  zone_a_boundary_definition.zone_id = 21U;
  zone_a_boundary_definition.boundary = true;
  zone_a_boundary_definition.default_kind =
    sqmesh::mesh::EntityKind::face_triangle;
  zone_a_boundary_definition.semantic = sqmesh::mesh::EntityGroupSemantic::boundary;
  zone_a_boundary_definition.primary_region_zone_id = 1U;
  zone_a_boundary_definition.import_info = make_cgns_entity_group_import_info(
    1U,
    "zone_a",
    "wall",
    static_cast<std::int32_t>(CGNS_ENUMV(BCWall))
  );
  const auto zone_a_boundary_entity_group =
    domain.create_entity_group(std::move(zone_a_boundary_definition));

  sqmesh::mesh::EntityGroupDefinition zone_b_volume_definition;
  zone_b_volume_definition.order = sqmesh::mesh::EntityOrder::cell;
  zone_b_volume_definition.name = "z2_zone_b__volume";
  zone_b_volume_definition.zone_id = 31U;
  zone_b_volume_definition.default_kind = sqmesh::mesh::EntityKind::cell_tetra;
  zone_b_volume_definition.semantic = sqmesh::mesh::EntityGroupSemantic::region;
  zone_b_volume_definition.primary_region_zone_id = 2U;
  zone_b_volume_definition.import_info =
    make_cgns_entity_group_import_info(2U, "zone_b", "volume");
  const auto zone_b_volume_entity_group =
    domain.create_entity_group(std::move(zone_b_volume_definition));

  sqmesh::mesh::EntityGroupDefinition zone_b_boundary_definition;
  zone_b_boundary_definition.order = sqmesh::mesh::EntityOrder::face;
  zone_b_boundary_definition.name = "z2_zone_b__wall";
  zone_b_boundary_definition.zone_id = 21U;
  zone_b_boundary_definition.boundary = true;
  zone_b_boundary_definition.default_kind =
    sqmesh::mesh::EntityKind::face_triangle;
  zone_b_boundary_definition.semantic = sqmesh::mesh::EntityGroupSemantic::boundary;
  zone_b_boundary_definition.primary_region_zone_id = 2U;
  zone_b_boundary_definition.import_info = make_cgns_entity_group_import_info(
    2U,
    "zone_b",
    "wall",
    static_cast<std::int32_t>(CGNS_ENUMV(BCWallInviscid))
  );
  const auto zone_b_boundary_entity_group =
    domain.create_entity_group(std::move(zone_b_boundary_definition));

  sqmesh::mesh::EntityGroupDefinition interface_definition;
  interface_definition.order = sqmesh::mesh::EntityOrder::face;
  interface_definition.name = "zone_interface";
  interface_definition.default_kind = sqmesh::mesh::EntityKind::face_triangle;
  interface_definition.semantic = sqmesh::mesh::EntityGroupSemantic::interface;
  interface_definition.primary_region_zone_id = 1U;
  interface_definition.secondary_region_zone_id = 2U;
  const auto interface_entity_group = domain.create_entity_group(std::move(interface_definition));

  const auto n0 = domain.add_node(node_entity_group, {0.0, 0.0, 0.0});
  const auto n1 = domain.add_node(node_entity_group, {1.0, 0.0, 0.0});
  const auto n2 = domain.add_node(node_entity_group, {0.0, 1.0, 0.0});
  const auto n3 = domain.add_node(node_entity_group, {0.0, 0.0, 1.0});
  const auto n4 = domain.add_node(node_entity_group, {0.0, 0.0, -1.0});

  const auto nb0 =
    duplicate_secondary_interface_nodes
      ? domain.add_node(node_entity_group, {0.0, 0.0, 0.0})
      : n0;
  const auto nb1 =
    duplicate_secondary_interface_nodes
      ? domain.add_node(node_entity_group, {1.0, 0.0, 0.0})
      : n1;
  const auto nb2 =
    duplicate_secondary_interface_nodes
      ? domain.add_node(node_entity_group, {0.0, 1.0, 0.0})
      : n2;

  const auto zone_a_boundary_face_0 =
    domain.add_triangle_face(zone_a_boundary_entity_group, {n0, n1, n3});
  const auto zone_a_boundary_face_1 =
    domain.add_triangle_face(zone_a_boundary_entity_group, {n1, n2, n3});
  const auto zone_a_boundary_face_2 =
    domain.add_triangle_face(zone_a_boundary_entity_group, {n2, n0, n3});
  const auto zone_b_boundary_face_0 =
    domain.add_triangle_face(zone_b_boundary_entity_group, {nb0, n4, nb1});
  const auto zone_b_boundary_face_1 =
    domain.add_triangle_face(zone_b_boundary_entity_group, {nb1, n4, nb2});
  const auto zone_b_boundary_face_2 =
    domain.add_triangle_face(zone_b_boundary_entity_group, {nb2, n4, nb0});
  const auto interface_face =
    domain.add_triangle_face(interface_entity_group, {n0, n1, n2});

  const auto zone_a_cell = domain.add_tetra_cell(
    zone_a_volume_entity_group,
    {n0, n1, n2, n3},
    {
      interface_face,
      zone_a_boundary_face_0,
      zone_a_boundary_face_1,
      zone_a_boundary_face_2,
    }
  );
  const auto zone_b_cell = domain.add_tetra_cell(
    zone_b_volume_entity_group,
    {nb0, nb1, nb2, n4},
    {
      interface_face,
      zone_b_boundary_face_0,
      zone_b_boundary_face_1,
      zone_b_boundary_face_2,
    }
  );

  domain.set_face_cells(zone_a_boundary_face_0, zone_a_cell);
  domain.set_face_cells(zone_a_boundary_face_1, zone_a_cell);
  domain.set_face_cells(zone_a_boundary_face_2, zone_a_cell);
  domain.set_face_cells(zone_b_boundary_face_0, zone_b_cell);
  domain.set_face_cells(zone_b_boundary_face_1, zone_b_cell);
  domain.set_face_cells(zone_b_boundary_face_2, zone_b_cell);
  domain.set_face_cells(interface_face, zone_a_cell, zone_b_cell);

  return domain;
}

sqmesh::mesh::Domain make_cgns_supported_multi_zone_interface_domain()
{
  return make_cgns_multi_zone_interface_domain(false);
}

sqmesh::mesh::Domain make_cgns_unsupported_multi_zone_interface_domain()
{
  return make_cgns_multi_zone_interface_domain(true);
}

sqmesh::mesh::Domain make_cgns_surface_multi_zone_interface_domain()
{
  sqmesh::mesh::Domain domain("manual_cgns_surface_shared_interface");
  const auto node_entity_group =
    domain.create_entity_group({sqmesh::mesh::EntityOrder::node, "nodes"});

  sqmesh::mesh::EntityGroupDefinition zone_a_face_definition;
  zone_a_face_definition.order = sqmesh::mesh::EntityOrder::face;
  zone_a_face_definition.name = "z1_zone_a__surface_patch";
  zone_a_face_definition.zone_id = 21U;
  zone_a_face_definition.boundary = true;
  zone_a_face_definition.default_kind =
    sqmesh::mesh::EntityKind::face_triangle;
  zone_a_face_definition.semantic = sqmesh::mesh::EntityGroupSemantic::boundary;
  zone_a_face_definition.primary_region_zone_id = 1U;
  zone_a_face_definition.import_info = make_cgns_entity_group_import_info(
    1U,
    "zone_a",
    "surface_patch",
    static_cast<std::int32_t>(CGNS_ENUMV(BCWall))
  );
  const auto zone_a_face_entity_group =
    domain.create_entity_group(std::move(zone_a_face_definition));

  sqmesh::mesh::EntityGroupDefinition zone_a_edge_definition;
  zone_a_edge_definition.order = sqmesh::mesh::EntityOrder::edge;
  zone_a_edge_definition.name = "z1_zone_a__outer_loop";
  zone_a_edge_definition.zone_id = 11U;
  zone_a_edge_definition.boundary = true;
  zone_a_edge_definition.default_kind = sqmesh::mesh::EntityKind::edge_line;
  zone_a_edge_definition.semantic = sqmesh::mesh::EntityGroupSemantic::boundary;
  zone_a_edge_definition.primary_region_zone_id = 1U;
  zone_a_edge_definition.import_info = make_cgns_entity_group_import_info(
    1U,
    "zone_a",
    "outer_loop",
    static_cast<std::int32_t>(CGNS_ENUMV(BCSymmetryPlane))
  );
  const auto zone_a_edge_entity_group =
    domain.create_entity_group(std::move(zone_a_edge_definition));

  sqmesh::mesh::EntityGroupDefinition zone_b_face_definition;
  zone_b_face_definition.order = sqmesh::mesh::EntityOrder::face;
  zone_b_face_definition.name = "z2_zone_b__surface_patch";
  zone_b_face_definition.zone_id = 21U;
  zone_b_face_definition.boundary = true;
  zone_b_face_definition.default_kind =
    sqmesh::mesh::EntityKind::face_triangle;
  zone_b_face_definition.semantic = sqmesh::mesh::EntityGroupSemantic::boundary;
  zone_b_face_definition.primary_region_zone_id = 2U;
  zone_b_face_definition.import_info = make_cgns_entity_group_import_info(
    2U,
    "zone_b",
    "surface_patch",
    static_cast<std::int32_t>(CGNS_ENUMV(BCWallInviscid))
  );
  const auto zone_b_face_entity_group =
    domain.create_entity_group(std::move(zone_b_face_definition));

  sqmesh::mesh::EntityGroupDefinition zone_b_edge_definition;
  zone_b_edge_definition.order = sqmesh::mesh::EntityOrder::edge;
  zone_b_edge_definition.name = "z2_zone_b__outer_loop";
  zone_b_edge_definition.zone_id = 11U;
  zone_b_edge_definition.boundary = true;
  zone_b_edge_definition.default_kind = sqmesh::mesh::EntityKind::edge_line;
  zone_b_edge_definition.semantic = sqmesh::mesh::EntityGroupSemantic::boundary;
  zone_b_edge_definition.primary_region_zone_id = 2U;
  zone_b_edge_definition.import_info = make_cgns_entity_group_import_info(
    2U,
    "zone_b",
    "outer_loop",
    static_cast<std::int32_t>(CGNS_ENUMV(BCOutflow))
  );
  const auto zone_b_edge_entity_group =
    domain.create_entity_group(std::move(zone_b_edge_definition));

  sqmesh::mesh::EntityGroupDefinition interface_definition;
  interface_definition.order = sqmesh::mesh::EntityOrder::edge;
  interface_definition.name = "zone_interface";
  interface_definition.default_kind = sqmesh::mesh::EntityKind::edge_line;
  interface_definition.semantic = sqmesh::mesh::EntityGroupSemantic::interface;
  interface_definition.primary_region_zone_id = 1U;
  interface_definition.secondary_region_zone_id = 2U;
  const auto interface_entity_group =
    domain.create_entity_group(std::move(interface_definition));

  const auto n0 = domain.add_node(node_entity_group, {0.0, 0.0, 0.0});
  const auto n1 = domain.add_node(node_entity_group, {1.0, 0.0, 0.0});
  const auto n2 = domain.add_node(node_entity_group, {0.0, 1.0, 0.0});
  const auto n3 = domain.add_node(node_entity_group, {1.0, 1.0, 0.0});

  const auto zone_a_face =
    domain.add_triangle_face(zone_a_face_entity_group, {n0, n1, n2});
  const auto zone_b_face =
    domain.add_triangle_face(zone_b_face_entity_group, {n0, n1, n3});

  const auto zone_a_edge_0 =
    domain.add_edge(zone_a_edge_entity_group, {n1, n2});
  const auto zone_a_edge_1 =
    domain.add_edge(zone_a_edge_entity_group, {n2, n0});
  const auto zone_b_edge_0 =
    domain.add_edge(zone_b_edge_entity_group, {n1, n3});
  const auto zone_b_edge_1 =
    domain.add_edge(zone_b_edge_entity_group, {n3, n0});
  const auto interface_edge =
    domain.add_edge(interface_entity_group, {n0, n1});

  domain.set_edge_faces(zone_a_edge_0, zone_a_face);
  domain.set_edge_faces(zone_a_edge_1, zone_a_face);
  domain.set_edge_faces(zone_b_edge_0, zone_b_face);
  domain.set_edge_faces(zone_b_edge_1, zone_b_face);
  domain.set_edge_faces(interface_edge, zone_a_face, zone_b_face);

  return domain;
}

sqmesh::mesh::Domain make_cgns_multi_face_interface_domain()
{
  sqmesh::mesh::Domain domain("manual_cgns_shared_patch_interface");
  const auto node_entity_group =
    domain.create_entity_group({sqmesh::mesh::EntityOrder::node, "nodes"});

  sqmesh::mesh::EntityGroupDefinition zone_a_volume_definition;
  zone_a_volume_definition.order = sqmesh::mesh::EntityOrder::cell;
  zone_a_volume_definition.name = "z1_zone_a__volume";
  zone_a_volume_definition.zone_id = 31U;
  zone_a_volume_definition.default_kind = sqmesh::mesh::EntityKind::cell_tetra;
  zone_a_volume_definition.semantic = sqmesh::mesh::EntityGroupSemantic::region;
  zone_a_volume_definition.primary_region_zone_id = 1U;
  zone_a_volume_definition.import_info =
    make_cgns_entity_group_import_info(1U, "zone_a", "volume");
  const auto zone_a_volume_entity_group =
    domain.create_entity_group(std::move(zone_a_volume_definition));

  sqmesh::mesh::EntityGroupDefinition zone_a_boundary_definition;
  zone_a_boundary_definition.order = sqmesh::mesh::EntityOrder::face;
  zone_a_boundary_definition.name = "z1_zone_a__wall";
  zone_a_boundary_definition.zone_id = 21U;
  zone_a_boundary_definition.boundary = true;
  zone_a_boundary_definition.default_kind =
    sqmesh::mesh::EntityKind::face_triangle;
  zone_a_boundary_definition.semantic = sqmesh::mesh::EntityGroupSemantic::boundary;
  zone_a_boundary_definition.primary_region_zone_id = 1U;
  zone_a_boundary_definition.import_info = make_cgns_entity_group_import_info(
    1U,
    "zone_a",
    "wall",
    static_cast<std::int32_t>(CGNS_ENUMV(BCWall))
  );
  const auto zone_a_boundary_entity_group =
    domain.create_entity_group(std::move(zone_a_boundary_definition));

  sqmesh::mesh::EntityGroupDefinition zone_b_volume_definition;
  zone_b_volume_definition.order = sqmesh::mesh::EntityOrder::cell;
  zone_b_volume_definition.name = "z2_zone_b__volume";
  zone_b_volume_definition.zone_id = 31U;
  zone_b_volume_definition.default_kind = sqmesh::mesh::EntityKind::cell_tetra;
  zone_b_volume_definition.semantic = sqmesh::mesh::EntityGroupSemantic::region;
  zone_b_volume_definition.primary_region_zone_id = 2U;
  zone_b_volume_definition.import_info =
    make_cgns_entity_group_import_info(2U, "zone_b", "volume");
  const auto zone_b_volume_entity_group =
    domain.create_entity_group(std::move(zone_b_volume_definition));

  sqmesh::mesh::EntityGroupDefinition zone_b_boundary_definition;
  zone_b_boundary_definition.order = sqmesh::mesh::EntityOrder::face;
  zone_b_boundary_definition.name = "z2_zone_b__wall";
  zone_b_boundary_definition.zone_id = 21U;
  zone_b_boundary_definition.boundary = true;
  zone_b_boundary_definition.default_kind =
    sqmesh::mesh::EntityKind::face_triangle;
  zone_b_boundary_definition.semantic = sqmesh::mesh::EntityGroupSemantic::boundary;
  zone_b_boundary_definition.primary_region_zone_id = 2U;
  zone_b_boundary_definition.import_info = make_cgns_entity_group_import_info(
    2U,
    "zone_b",
    "wall",
    static_cast<std::int32_t>(CGNS_ENUMV(BCWallInviscid))
  );
  const auto zone_b_boundary_entity_group =
    domain.create_entity_group(std::move(zone_b_boundary_definition));

  sqmesh::mesh::EntityGroupDefinition interface_definition;
  interface_definition.order = sqmesh::mesh::EntityOrder::face;
  interface_definition.name = "zone_interface";
  interface_definition.default_kind = sqmesh::mesh::EntityKind::face_triangle;
  interface_definition.semantic = sqmesh::mesh::EntityGroupSemantic::interface;
  interface_definition.primary_region_zone_id = 1U;
  interface_definition.secondary_region_zone_id = 2U;
  const auto interface_entity_group = domain.create_entity_group(std::move(interface_definition));

  const auto n0 = domain.add_node(node_entity_group, {0.0, 0.0, 0.0});
  const auto n1 = domain.add_node(node_entity_group, {1.0, 0.0, 0.0});
  const auto n2 = domain.add_node(node_entity_group, {1.0, 1.0, 0.0});
  const auto n3 = domain.add_node(node_entity_group, {0.0, 1.0, 0.0});
  const auto n4 = domain.add_node(node_entity_group, {0.0, 0.0, 1.0});
  const auto n5 = domain.add_node(node_entity_group, {0.0, 0.0, -1.0});
  const auto n6 = domain.add_node(node_entity_group, {0.0, 1.0, 1.0});
  const auto n7 = domain.add_node(node_entity_group, {0.0, 1.0, -1.0});

  const auto zone_a_boundary_face_0 =
    domain.add_triangle_face(zone_a_boundary_entity_group, {n0, n1, n4});
  const auto zone_a_boundary_face_1 =
    domain.add_triangle_face(zone_a_boundary_entity_group, {n1, n2, n4});
  const auto zone_a_boundary_face_2 =
    domain.add_triangle_face(zone_a_boundary_entity_group, {n2, n0, n4});
  const auto zone_a_boundary_face_3 =
    domain.add_triangle_face(zone_a_boundary_entity_group, {n0, n2, n6});
  const auto zone_a_boundary_face_4 =
    domain.add_triangle_face(zone_a_boundary_entity_group, {n2, n3, n6});
  const auto zone_a_boundary_face_5 =
    domain.add_triangle_face(zone_a_boundary_entity_group, {n3, n0, n6});
  const auto zone_b_boundary_face_0 =
    domain.add_triangle_face(zone_b_boundary_entity_group, {n0, n5, n1});
  const auto zone_b_boundary_face_1 =
    domain.add_triangle_face(zone_b_boundary_entity_group, {n1, n5, n2});
  const auto zone_b_boundary_face_2 =
    domain.add_triangle_face(zone_b_boundary_entity_group, {n2, n5, n0});
  const auto zone_b_boundary_face_3 =
    domain.add_triangle_face(zone_b_boundary_entity_group, {n0, n7, n2});
  const auto zone_b_boundary_face_4 =
    domain.add_triangle_face(zone_b_boundary_entity_group, {n2, n7, n3});
  const auto zone_b_boundary_face_5 =
    domain.add_triangle_face(zone_b_boundary_entity_group, {n3, n7, n0});
  const auto interface_face_0 =
    domain.add_triangle_face(interface_entity_group, {n0, n1, n2});
  const auto interface_face_1 =
    domain.add_triangle_face(interface_entity_group, {n0, n2, n3});

  const auto zone_a_cell_0 = domain.add_tetra_cell(
    zone_a_volume_entity_group,
    {n0, n1, n2, n4},
    {
      interface_face_0,
      zone_a_boundary_face_0,
      zone_a_boundary_face_1,
      zone_a_boundary_face_2,
    }
  );
  const auto zone_a_cell_1 = domain.add_tetra_cell(
    zone_a_volume_entity_group,
    {n0, n2, n3, n6},
    {
      interface_face_1,
      zone_a_boundary_face_3,
      zone_a_boundary_face_4,
      zone_a_boundary_face_5,
    }
  );
  const auto zone_b_cell_0 = domain.add_tetra_cell(
    zone_b_volume_entity_group,
    {n0, n1, n2, n5},
    {
      interface_face_0,
      zone_b_boundary_face_0,
      zone_b_boundary_face_1,
      zone_b_boundary_face_2,
    }
  );
  const auto zone_b_cell_1 = domain.add_tetra_cell(
    zone_b_volume_entity_group,
    {n0, n2, n3, n7},
    {
      interface_face_1,
      zone_b_boundary_face_3,
      zone_b_boundary_face_4,
      zone_b_boundary_face_5,
    }
  );

  domain.set_face_cells(zone_a_boundary_face_0, zone_a_cell_0);
  domain.set_face_cells(zone_a_boundary_face_1, zone_a_cell_0);
  domain.set_face_cells(zone_a_boundary_face_2, zone_a_cell_0);
  domain.set_face_cells(zone_a_boundary_face_3, zone_a_cell_1);
  domain.set_face_cells(zone_a_boundary_face_4, zone_a_cell_1);
  domain.set_face_cells(zone_a_boundary_face_5, zone_a_cell_1);
  domain.set_face_cells(zone_b_boundary_face_0, zone_b_cell_0);
  domain.set_face_cells(zone_b_boundary_face_1, zone_b_cell_0);
  domain.set_face_cells(zone_b_boundary_face_2, zone_b_cell_0);
  domain.set_face_cells(zone_b_boundary_face_3, zone_b_cell_1);
  domain.set_face_cells(zone_b_boundary_face_4, zone_b_cell_1);
  domain.set_face_cells(zone_b_boundary_face_5, zone_b_cell_1);
  domain.set_face_cells(interface_face_0, zone_a_cell_0, zone_b_cell_0);
  domain.set_face_cells(interface_face_1, zone_a_cell_1, zone_b_cell_1);

  return domain;
}

sqmesh::mesh::Domain make_cgns_surface_multi_edge_interface_domain()
{
  sqmesh::mesh::Domain domain("manual_cgns_surface_polyline_interface");
  const auto node_entity_group =
    domain.create_entity_group({sqmesh::mesh::EntityOrder::node, "nodes"});

  sqmesh::mesh::EntityGroupDefinition zone_a_face_definition;
  zone_a_face_definition.order = sqmesh::mesh::EntityOrder::face;
  zone_a_face_definition.name = "z1_zone_a__surface_patch";
  zone_a_face_definition.zone_id = 21U;
  zone_a_face_definition.boundary = true;
  zone_a_face_definition.default_kind =
    sqmesh::mesh::EntityKind::face_triangle;
  zone_a_face_definition.semantic = sqmesh::mesh::EntityGroupSemantic::boundary;
  zone_a_face_definition.primary_region_zone_id = 1U;
  zone_a_face_definition.import_info = make_cgns_entity_group_import_info(
    1U,
    "zone_a",
    "surface_patch",
    static_cast<std::int32_t>(CGNS_ENUMV(BCWall))
  );
  const auto zone_a_face_entity_group =
    domain.create_entity_group(std::move(zone_a_face_definition));

  sqmesh::mesh::EntityGroupDefinition zone_a_edge_definition;
  zone_a_edge_definition.order = sqmesh::mesh::EntityOrder::edge;
  zone_a_edge_definition.name = "z1_zone_a__outer_loop";
  zone_a_edge_definition.zone_id = 11U;
  zone_a_edge_definition.boundary = true;
  zone_a_edge_definition.default_kind = sqmesh::mesh::EntityKind::edge_line;
  zone_a_edge_definition.semantic = sqmesh::mesh::EntityGroupSemantic::boundary;
  zone_a_edge_definition.primary_region_zone_id = 1U;
  zone_a_edge_definition.import_info = make_cgns_entity_group_import_info(
    1U,
    "zone_a",
    "outer_loop",
    static_cast<std::int32_t>(CGNS_ENUMV(BCSymmetryPlane))
  );
  const auto zone_a_edge_entity_group =
    domain.create_entity_group(std::move(zone_a_edge_definition));

  sqmesh::mesh::EntityGroupDefinition zone_b_face_definition;
  zone_b_face_definition.order = sqmesh::mesh::EntityOrder::face;
  zone_b_face_definition.name = "z2_zone_b__surface_patch";
  zone_b_face_definition.zone_id = 21U;
  zone_b_face_definition.boundary = true;
  zone_b_face_definition.default_kind =
    sqmesh::mesh::EntityKind::face_triangle;
  zone_b_face_definition.semantic = sqmesh::mesh::EntityGroupSemantic::boundary;
  zone_b_face_definition.primary_region_zone_id = 2U;
  zone_b_face_definition.import_info = make_cgns_entity_group_import_info(
    2U,
    "zone_b",
    "surface_patch",
    static_cast<std::int32_t>(CGNS_ENUMV(BCWallInviscid))
  );
  const auto zone_b_face_entity_group =
    domain.create_entity_group(std::move(zone_b_face_definition));

  sqmesh::mesh::EntityGroupDefinition zone_b_edge_definition;
  zone_b_edge_definition.order = sqmesh::mesh::EntityOrder::edge;
  zone_b_edge_definition.name = "z2_zone_b__outer_loop";
  zone_b_edge_definition.zone_id = 11U;
  zone_b_edge_definition.boundary = true;
  zone_b_edge_definition.default_kind = sqmesh::mesh::EntityKind::edge_line;
  zone_b_edge_definition.semantic = sqmesh::mesh::EntityGroupSemantic::boundary;
  zone_b_edge_definition.primary_region_zone_id = 2U;
  zone_b_edge_definition.import_info = make_cgns_entity_group_import_info(
    2U,
    "zone_b",
    "outer_loop",
    static_cast<std::int32_t>(CGNS_ENUMV(BCOutflow))
  );
  const auto zone_b_edge_entity_group =
    domain.create_entity_group(std::move(zone_b_edge_definition));

  sqmesh::mesh::EntityGroupDefinition interface_definition;
  interface_definition.order = sqmesh::mesh::EntityOrder::edge;
  interface_definition.name = "zone_interface";
  interface_definition.default_kind = sqmesh::mesh::EntityKind::edge_line;
  interface_definition.semantic = sqmesh::mesh::EntityGroupSemantic::interface;
  interface_definition.primary_region_zone_id = 1U;
  interface_definition.secondary_region_zone_id = 2U;
  const auto interface_entity_group =
    domain.create_entity_group(std::move(interface_definition));

  const auto n0 = domain.add_node(node_entity_group, {0.0, 0.0, 0.0});
  const auto n1 = domain.add_node(node_entity_group, {1.0, 0.0, 0.0});
  const auto n2 = domain.add_node(node_entity_group, {2.0, 0.0, 0.0});
  const auto n3 = domain.add_node(node_entity_group, {0.0, 1.0, 0.0});
  const auto n4 = domain.add_node(node_entity_group, {2.0, 1.0, 0.0});
  const auto n5 = domain.add_node(node_entity_group, {0.0, -1.0, 0.0});
  const auto n6 = domain.add_node(node_entity_group, {2.0, -1.0, 0.0});

  const auto zone_a_face_0 =
    domain.add_triangle_face(zone_a_face_entity_group, {n0, n1, n3});
  const auto zone_a_face_1 =
    domain.add_triangle_face(zone_a_face_entity_group, {n1, n2, n4});
  const auto zone_b_face_0 =
    domain.add_triangle_face(zone_b_face_entity_group, {n0, n5, n1});
  const auto zone_b_face_1 =
    domain.add_triangle_face(zone_b_face_entity_group, {n1, n6, n2});

  const auto zone_a_edge_0 =
    domain.add_edge(zone_a_edge_entity_group, {n1, n3});
  const auto zone_a_edge_1 =
    domain.add_edge(zone_a_edge_entity_group, {n3, n0});
  const auto zone_a_edge_2 =
    domain.add_edge(zone_a_edge_entity_group, {n2, n4});
  const auto zone_a_edge_3 =
    domain.add_edge(zone_a_edge_entity_group, {n4, n1});
  const auto zone_b_edge_0 =
    domain.add_edge(zone_b_edge_entity_group, {n0, n5});
  const auto zone_b_edge_1 =
    domain.add_edge(zone_b_edge_entity_group, {n5, n1});
  const auto zone_b_edge_2 =
    domain.add_edge(zone_b_edge_entity_group, {n1, n6});
  const auto zone_b_edge_3 =
    domain.add_edge(zone_b_edge_entity_group, {n6, n2});
  const auto interface_edge_0 =
    domain.add_edge(interface_entity_group, {n0, n1});
  const auto interface_edge_1 =
    domain.add_edge(interface_entity_group, {n1, n2});

  domain.set_edge_faces(zone_a_edge_0, zone_a_face_0);
  domain.set_edge_faces(zone_a_edge_1, zone_a_face_0);
  domain.set_edge_faces(zone_a_edge_2, zone_a_face_1);
  domain.set_edge_faces(zone_a_edge_3, zone_a_face_1);
  domain.set_edge_faces(zone_b_edge_0, zone_b_face_0);
  domain.set_edge_faces(zone_b_edge_1, zone_b_face_0);
  domain.set_edge_faces(zone_b_edge_2, zone_b_face_1);
  domain.set_edge_faces(zone_b_edge_3, zone_b_face_1);
  domain.set_edge_faces(interface_edge_0, zone_a_face_0, zone_b_face_0);
  domain.set_edge_faces(interface_edge_1, zone_a_face_1, zone_b_face_1);

  return domain;
}

#endif

#ifdef sqmesh_TEST_CGNS_ENABLED

struct CgnsSectionExpectation final {
  const char *name = nullptr;
  CGNS_ENUMT(ElementType_t) type = CGNS_ENUMV(ElementTypeNull);
  const char *zone_id = nullptr;
  cgsize_t start = 0;
  cgsize_t end = 0;
};

struct CgnsBoundaryExpectation final {
  const char *name = nullptr;
  CGNS_ENUMT(BCType_t) type = CGNS_ENUMV(BCTypeNull);
  CGNS_ENUMT(PointSetType_t) point_set_type = CGNS_ENUMV(PointSetTypeNull);
  CGNS_ENUMT(GridLocation_t) location = CGNS_ENUMV(GridLocationNull);
  std::vector<cgsize_t> points {};
};

struct CgnsConnectivityExpectation final {
  const char *name = nullptr;
  CGNS_ENUMT(GridLocation_t) location = CGNS_ENUMV(GridLocationNull);
  CGNS_ENUMT(GridConnectivityType_t) type =
    CGNS_ENUMV(GridConnectivityTypeNull);
  CGNS_ENUMT(PointSetType_t) point_set_type = CGNS_ENUMV(PointSetTypeNull);
  const char *donor_name = nullptr;
  CGNS_ENUMT(ZoneType_t) donor_zone_type = CGNS_ENUMV(ZoneTypeNull);
  CGNS_ENUMT(PointSetType_t) donor_point_set_type =
    CGNS_ENUMV(PointSetTypeNull);
  CGNS_ENUMT(DataType_t) donor_data_type = CGNS_ENUMV(DataTypeNull);
  std::vector<cgsize_t> points {};
  std::vector<cgsize_t> donor_points {};
};

struct CgnsZoneExpectation final {
  const char *name = nullptr;
  std::vector<CgnsSectionExpectation> sections {};
  std::vector<CgnsBoundaryExpectation> boundaries {};
  std::vector<CgnsConnectivityExpectation> connections {};
};

bool check_cgns(int status, const char *message)
{
  if(status == CG_OK) {
    return true;
  }

  std::fprintf(
    stderr,
    "mesh_io_smoke: %s (%s)\n",
    message,
    cg_get_error()
  );
  return false;
}

bool verify_cgns_entity_group_import_info(
  const sqmesh::mesh::EntityGroup &entity_group,
  std::uint32_t expected_base_index,
  std::string_view expected_base_name,
  std::uint32_t expected_zone_index,
  std::string_view expected_zone_name,
  std::string_view expected_local_name,
  std::int32_t expected_bc_type_value,
  std::string_view label
)
{
  const auto &import_info = entity_group.import_info();
  if(!expect(
       import_info.format == sqmesh::mesh::EntityGroupImportFormat::cgns,
       (std::string(label) + " should publish CGNS import metadata").c_str()
     )) {
    return false;
  }
  if(expected_base_index != sqmesh::mesh::invalid_index ||
     !expected_base_name.empty()) {
    if(!expect(
         import_info.cgns.base_index == expected_base_index,
         (std::string(label) + " should preserve the expected CGNS base index").c_str()
       )) {
      return false;
    }
    if(!expect(
         import_info.cgns.base_name == expected_base_name,
         (std::string(label) + " should preserve the expected CGNS base name").c_str()
       )) {
      return false;
    }
  }
  if(!expect(
       import_info.cgns.zone_index == expected_zone_index,
       (std::string(label) + " should preserve the expected CGNS zone index").c_str()
     )) {
    return false;
  }
  if(!expect(
       import_info.cgns.zone_name == expected_zone_name,
       (std::string(label) + " should preserve the expected CGNS zone name").c_str()
     )) {
    return false;
  }
  if(!expect(
       import_info.cgns.local_name == expected_local_name,
       (std::string(label) + " should preserve the expected CGNS local name").c_str()
     )) {
    return false;
  }
  if(!expect(
       import_info.cgns.bc_type_value == expected_bc_type_value,
       (std::string(label) + " should preserve the expected raw CGNS BCType_t value").c_str()
     )) {
    return false;
  }
  return true;
}

bool write_cgns_section_zone_id_descriptor(
  int file_index,
  int base_index,
  int zone_index,
  int section_index,
  const char *zone_id
)
{
  return check_cgns(
           cg_goto(file_index, base_index, "Zone_t", zone_index, "Elements_t", section_index, "end"),
           "failed to access CGNS section metadata"
         ) &&
         check_cgns(
           cg_descriptor_write("SQMeshZoneId", zone_id),
           "failed to write CGNS section zone-id descriptor"
         );
}

bool write_cgns_boundary_condition(
  int file_index,
  int base_index,
  int zone_index,
  const char *name,
  CGNS_ENUMT(BCType_t) boundary_type,
  CGNS_ENUMT(PointSetType_t) point_set_type,
  const std::vector<cgsize_t> &points,
  CGNS_ENUMT(GridLocation_t) location,
  bool write_location
)
{
  int boundary_index = 0;
  if(!check_cgns(
       cg_boco_write(
         file_index,
         base_index,
         zone_index,
         name,
         boundary_type,
         point_set_type,
         static_cast<cgsize_t>(points.size()),
         points.data(),
         &boundary_index
       ),
       "failed to write CGNS boundary condition"
     )) {
    return false;
  }

  if(!write_location) {
    return true;
  }

  return check_cgns(
    cg_boco_gridlocation_write(file_index, base_index, zone_index, boundary_index, location),
    "failed to write CGNS boundary-condition grid location"
  );
}

bool write_supported_cgns_fixture(const std::filesystem::path &path)
{
  int file_index = 0;
  if(!check_cgns(cg_open(path.string().c_str(), CG_MODE_WRITE, &file_index), "failed to open CGNS fixture")) {
    return false;
  }

  int base_index = 0;
  if(!check_cgns(cg_base_write(file_index, "sqmesh_base", 3, 3, &base_index), "failed to write CGNS base")) {
    static_cast<void>(cg_close(file_index));
    return false;
  }

  const cgsize_t zone_size[3] = {4, 1, 0};
  int zone_index = 0;
  if(!check_cgns(
       cg_zone_write(
         file_index,
         base_index,
         "sqmesh_cgns_zone",
         zone_size,
         CGNS_ENUMV(Unstructured),
         &zone_index
       ),
       "failed to write CGNS zone"
     )) {
    static_cast<void>(cg_close(file_index));
    return false;
  }

  const double x[4] = {0.0, 1.0, 0.0, 0.0};
  const double y[4] = {0.0, 0.0, 1.0, 0.0};
  const double z[4] = {0.0, 0.0, 0.0, 1.0};
  int coordinate_index = 0;
  if(!check_cgns(
       cg_coord_write(file_index, base_index, zone_index, CGNS_ENUMV(RealDouble), "CoordinateX", x, &coordinate_index),
       "failed to write CoordinateX"
     ) ||
     !check_cgns(
       cg_coord_write(file_index, base_index, zone_index, CGNS_ENUMV(RealDouble), "CoordinateY", y, &coordinate_index),
       "failed to write CoordinateY"
     ) ||
     !check_cgns(
       cg_coord_write(file_index, base_index, zone_index, CGNS_ENUMV(RealDouble), "CoordinateZ", z, &coordinate_index),
       "failed to write CoordinateZ"
     )) {
    static_cast<void>(cg_close(file_index));
    return false;
  }

  const cgsize_t edges[12] = {
    1, 2,
    2, 3,
    3, 1,
    1, 4,
    2, 4,
    3, 4,
  };
  const cgsize_t faces[12] = {
    1, 3, 2,
    1, 2, 4,
    2, 3, 4,
    3, 1, 4,
  };
  const cgsize_t cells[4] = {1, 2, 3, 4};

  int section_index = 0;
  if(!check_cgns(
       cg_section_write(file_index, base_index, zone_index, "edge_section", CGNS_ENUMV(BAR_2), 1, 6, 0, edges, &section_index),
       "failed to write CGNS edge section"
     ) ||
     !write_cgns_section_zone_id_descriptor(file_index, base_index, zone_index, section_index, "11")) {
    static_cast<void>(cg_close(file_index));
    return false;
  }

  if(!check_cgns(
       cg_section_write(file_index, base_index, zone_index, "surface_section", CGNS_ENUMV(TRI_3), 7, 10, 0, faces, &section_index),
       "failed to write CGNS face section"
     ) ||
     !write_cgns_section_zone_id_descriptor(file_index, base_index, zone_index, section_index, "21")) {
    static_cast<void>(cg_close(file_index));
    return false;
  }

  if(!check_cgns(
       cg_section_write(file_index, base_index, zone_index, "volume", CGNS_ENUMV(TETRA_4), 11, 11, 0, cells, &section_index),
       "failed to write CGNS cell section"
     ) ||
     !write_cgns_section_zone_id_descriptor(file_index, base_index, zone_index, section_index, "31") ||
     !write_cgns_boundary_condition(
       file_index,
       base_index,
       zone_index,
       "feature_edges",
       CGNS_ENUMV(BCSymmetryPlane),
       CGNS_ENUMV(PointRange),
       {1, 6},
       CGNS_ENUMV(EdgeCenter),
       true
     ) ||
     !write_cgns_boundary_condition(
       file_index,
       base_index,
       zone_index,
       "wall",
       CGNS_ENUMV(BCWall),
       CGNS_ENUMV(PointList),
       {7, 8, 9, 10},
       CGNS_ENUMV(FaceCenter),
       true
     )) {
    static_cast<void>(cg_close(file_index));
    return false;
  }

  return check_cgns(cg_close(file_index), "failed to close CGNS fixture");
}

bool write_mixed_surface_cgns_fixture(const std::filesystem::path &path)
{
  int file_index = 0;
  if(!check_cgns(cg_open(path.string().c_str(), CG_MODE_WRITE, &file_index), "failed to open mixed CGNS fixture")) {
    return false;
  }

  int base_index = 0;
  if(!check_cgns(cg_base_write(file_index, "sqmesh_base", 2, 3, &base_index), "failed to write mixed CGNS base")) {
    static_cast<void>(cg_close(file_index));
    return false;
  }

  const cgsize_t zone_size[3] = {5, 2, 0};
  int zone_index = 0;
  if(!check_cgns(
       cg_zone_write(
         file_index,
         base_index,
         "sqmesh_cgns_mixed_zone",
         zone_size,
         CGNS_ENUMV(Unstructured),
         &zone_index
       ),
       "failed to write mixed CGNS zone"
     )) {
    static_cast<void>(cg_close(file_index));
    return false;
  }

  const double x[5] = {0.0, 1.0, 1.0, 0.0, -0.25};
  const double y[5] = {0.0, 0.0, 1.0, 1.0, 0.5};
  const double z[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  int coordinate_index = 0;
  if(!check_cgns(
       cg_coord_write(file_index, base_index, zone_index, CGNS_ENUMV(RealDouble), "CoordinateX", x, &coordinate_index),
       "failed to write mixed CoordinateX"
     ) ||
     !check_cgns(
       cg_coord_write(file_index, base_index, zone_index, CGNS_ENUMV(RealDouble), "CoordinateY", y, &coordinate_index),
       "failed to write mixed CoordinateY"
     ) ||
     !check_cgns(
       cg_coord_write(file_index, base_index, zone_index, CGNS_ENUMV(RealDouble), "CoordinateZ", z, &coordinate_index),
       "failed to write mixed CoordinateZ"
     )) {
    static_cast<void>(cg_close(file_index));
    return false;
  }

  const cgsize_t mixed_elements[24] = {
    CGNS_ENUMV(BAR_2), 1, 2,
    CGNS_ENUMV(BAR_2), 2, 3,
    CGNS_ENUMV(BAR_2), 3, 4,
    CGNS_ENUMV(BAR_2), 4, 5,
    CGNS_ENUMV(BAR_2), 5, 1,
    CGNS_ENUMV(QUAD_4), 1, 2, 3, 4,
    CGNS_ENUMV(TRI_3), 1, 4, 5,
  };
  const cgsize_t mixed_offsets[8] = {0, 3, 6, 9, 12, 15, 20, 24};
  int section_index = 0;
  if(!check_cgns(
       cg_poly_section_write(
         file_index,
         base_index,
         zone_index,
         "mixed_surface_section",
         CGNS_ENUMV(MIXED),
         1,
         7,
         0,
         mixed_elements,
         mixed_offsets,
         &section_index
       ),
       "failed to write mixed CGNS section"
     ) ||
     !write_cgns_section_zone_id_descriptor(file_index, base_index, zone_index, section_index, "41") ||
     !write_cgns_boundary_condition(
       file_index,
       base_index,
       zone_index,
       "outer_loop",
       CGNS_ENUMV(BCSymmetryPlane),
       CGNS_ENUMV(PointRange),
       {1, 5},
       CGNS_ENUMV(EdgeCenter),
       true
     ) ||
     !write_cgns_boundary_condition(
       file_index,
       base_index,
       zone_index,
       "surface_patch",
       CGNS_ENUMV(BCInflow),
       CGNS_ENUMV(PointList),
       {6, 7},
       CGNS_ENUMV(CellCenter),
       true
     )) {
    static_cast<void>(cg_close(file_index));
    return false;
  }

  return check_cgns(cg_close(file_index), "failed to close mixed CGNS fixture");
}

bool write_multi_zone_cgns_fixture(const std::filesystem::path &path)
{
  int file_index = 0;
  if(!check_cgns(cg_open(path.string().c_str(), CG_MODE_WRITE, &file_index), "failed to open multi-zone CGNS fixture")) {
    return false;
  }

  int base_index = 0;
  if(!check_cgns(cg_base_write(file_index, "sqmesh_base", 3, 3, &base_index), "failed to write multi-zone CGNS base")) {
    static_cast<void>(cg_close(file_index));
    return false;
  }

  const cgsize_t zone_size[3] = {4, 1, 0};
  const cgsize_t edges[12] = {
    1, 2,
    2, 3,
    3, 1,
    1, 4,
    2, 4,
    3, 4,
  };
  const cgsize_t faces[12] = {
    1, 3, 2,
    1, 2, 4,
    2, 3, 4,
    3, 1, 4,
  };
  const cgsize_t cells[4] = {1, 2, 3, 4};

  const auto write_zone = [&](
                            const char *zone_name,
                            double x_offset,
                            CGNS_ENUMT(BCType_t) edge_boundary_type,
                            CGNS_ENUMT(BCType_t) face_boundary_type
                          ) {
    int zone_index = 0;
    if(!check_cgns(
         cg_zone_write(
           file_index,
           base_index,
           zone_name,
           zone_size,
           CGNS_ENUMV(Unstructured),
           &zone_index
         ),
         "failed to write multi-zone CGNS zone"
       )) {
      return false;
    }

    const double x[4] = {x_offset + 0.0, x_offset + 1.0, x_offset + 0.0, x_offset + 0.0};
    const double y[4] = {0.0, 0.0, 1.0, 0.0};
    const double z[4] = {0.0, 0.0, 0.0, 1.0};
    int coordinate_index = 0;
    if(!check_cgns(
         cg_coord_write(file_index, base_index, zone_index, CGNS_ENUMV(RealDouble), "CoordinateX", x, &coordinate_index),
         "failed to write multi-zone CoordinateX"
       ) ||
       !check_cgns(
         cg_coord_write(file_index, base_index, zone_index, CGNS_ENUMV(RealDouble), "CoordinateY", y, &coordinate_index),
         "failed to write multi-zone CoordinateY"
       ) ||
       !check_cgns(
         cg_coord_write(file_index, base_index, zone_index, CGNS_ENUMV(RealDouble), "CoordinateZ", z, &coordinate_index),
         "failed to write multi-zone CoordinateZ"
       )) {
      return false;
    }

    int section_index = 0;
    if(!check_cgns(
         cg_section_write(file_index, base_index, zone_index, "edge_section", CGNS_ENUMV(BAR_2), 1, 6, 0, edges, &section_index),
         "failed to write multi-zone CGNS edge section"
       ) ||
       !write_cgns_section_zone_id_descriptor(file_index, base_index, zone_index, section_index, "11")) {
      return false;
    }

    if(!check_cgns(
         cg_section_write(file_index, base_index, zone_index, "surface_section", CGNS_ENUMV(TRI_3), 7, 10, 0, faces, &section_index),
         "failed to write multi-zone CGNS face section"
       ) ||
       !write_cgns_section_zone_id_descriptor(file_index, base_index, zone_index, section_index, "21")) {
      return false;
    }

    return check_cgns(
             cg_section_write(file_index, base_index, zone_index, "volume", CGNS_ENUMV(TETRA_4), 11, 11, 0, cells, &section_index),
             "failed to write multi-zone CGNS cell section"
           ) &&
           write_cgns_section_zone_id_descriptor(file_index, base_index, zone_index, section_index, "31") &&
           write_cgns_boundary_condition(
             file_index,
             base_index,
             zone_index,
             "feature_edges",
             edge_boundary_type,
             CGNS_ENUMV(PointRange),
             {1, 6},
             CGNS_ENUMV(EdgeCenter),
             true
           ) &&
           write_cgns_boundary_condition(
             file_index,
             base_index,
             zone_index,
             "wall",
             face_boundary_type,
             CGNS_ENUMV(PointList),
             {7, 8, 9, 10},
             CGNS_ENUMV(FaceCenter),
             true
           );
  };

  if(!write_zone("zone_a", 0.0, CGNS_ENUMV(BCSymmetryPlane), CGNS_ENUMV(BCWall)) ||
     !write_zone("zone_b", 2.0, CGNS_ENUMV(BCOutflow), CGNS_ENUMV(BCWallInviscid))) {
    static_cast<void>(cg_close(file_index));
    return false;
  }

  return check_cgns(cg_close(file_index), "failed to close multi-zone CGNS fixture");
}

bool write_multi_base_cgns_fixture(const std::filesystem::path &path)
{
  int file_index = 0;
  if(!check_cgns(cg_open(path.string().c_str(), CG_MODE_WRITE, &file_index), "failed to open multi-base CGNS fixture")) {
    return false;
  }

  const cgsize_t zone_size[3] = {4, 1, 0};
  const cgsize_t edges[12] = {
    1, 2,
    2, 3,
    3, 1,
    1, 4,
    2, 4,
    3, 4,
  };
  const cgsize_t faces[12] = {
    1, 3, 2,
    1, 2, 4,
    2, 3, 4,
    3, 1, 4,
  };
  const cgsize_t cells[4] = {1, 2, 3, 4};

  const auto write_base = [&](
                            const char *base_name,
                            double x_offset,
                            CGNS_ENUMT(BCType_t) edge_boundary_type,
                            CGNS_ENUMT(BCType_t) face_boundary_type
                          ) {
    int base_index = 0;
    if(!check_cgns(
         cg_base_write(file_index, base_name, 3, 3, &base_index),
         "failed to write multi-base CGNS base"
       )) {
      return false;
    }

    int zone_index = 0;
    if(!check_cgns(
         cg_zone_write(
           file_index,
           base_index,
           "shared_zone",
           zone_size,
           CGNS_ENUMV(Unstructured),
           &zone_index
         ),
         "failed to write multi-base CGNS zone"
       )) {
      return false;
    }

    const double x[4] = {x_offset + 0.0, x_offset + 1.0, x_offset + 0.0, x_offset + 0.0};
    const double y[4] = {0.0, 0.0, 1.0, 0.0};
    const double z[4] = {0.0, 0.0, 0.0, 1.0};
    int coordinate_index = 0;
    if(!check_cgns(
         cg_coord_write(file_index, base_index, zone_index, CGNS_ENUMV(RealDouble), "CoordinateX", x, &coordinate_index),
         "failed to write multi-base CoordinateX"
       ) ||
       !check_cgns(
         cg_coord_write(file_index, base_index, zone_index, CGNS_ENUMV(RealDouble), "CoordinateY", y, &coordinate_index),
         "failed to write multi-base CoordinateY"
       ) ||
       !check_cgns(
         cg_coord_write(file_index, base_index, zone_index, CGNS_ENUMV(RealDouble), "CoordinateZ", z, &coordinate_index),
         "failed to write multi-base CoordinateZ"
       )) {
      return false;
    }

    int section_index = 0;
    if(!check_cgns(
         cg_section_write(file_index, base_index, zone_index, "edge_section", CGNS_ENUMV(BAR_2), 1, 6, 0, edges, &section_index),
         "failed to write multi-base CGNS edge section"
       ) ||
       !write_cgns_section_zone_id_descriptor(file_index, base_index, zone_index, section_index, "11")) {
      return false;
    }

    if(!check_cgns(
         cg_section_write(file_index, base_index, zone_index, "surface_section", CGNS_ENUMV(TRI_3), 7, 10, 0, faces, &section_index),
         "failed to write multi-base CGNS face section"
       ) ||
       !write_cgns_section_zone_id_descriptor(file_index, base_index, zone_index, section_index, "21")) {
      return false;
    }

    return check_cgns(
             cg_section_write(file_index, base_index, zone_index, "volume", CGNS_ENUMV(TETRA_4), 11, 11, 0, cells, &section_index),
             "failed to write multi-base CGNS cell section"
           ) &&
           write_cgns_section_zone_id_descriptor(file_index, base_index, zone_index, section_index, "31") &&
           write_cgns_boundary_condition(
             file_index,
             base_index,
             zone_index,
             "feature_edges",
             edge_boundary_type,
             CGNS_ENUMV(PointRange),
             {1, 6},
             CGNS_ENUMV(EdgeCenter),
             true
           ) &&
           write_cgns_boundary_condition(
             file_index,
             base_index,
             zone_index,
             "wall",
             face_boundary_type,
             CGNS_ENUMV(PointList),
             {7, 8, 9, 10},
             CGNS_ENUMV(FaceCenter),
             true
           );
  };

  if(!write_base("sqmesh_base_a", 0.0, CGNS_ENUMV(BCSymmetryPlane), CGNS_ENUMV(BCWall)) ||
     !write_base("sqmesh_base_b", 2.0, CGNS_ENUMV(BCOutflow), CGNS_ENUMV(BCWallInviscid))) {
    static_cast<void>(cg_close(file_index));
    return false;
  }

  return check_cgns(cg_close(file_index), "failed to close multi-base CGNS fixture");
}

bool write_unsupported_cgns_physical_dim_fixture(const std::filesystem::path &path)
{
  int file_index = 0;
  if(!check_cgns(
       cg_open(path.string().c_str(), CG_MODE_WRITE, &file_index),
       "failed to open unsupported CGNS physical-dimension fixture"
     )) {
    return false;
  }

  const auto write_surface_base =
    [&](const char *base_name, int physical_dim, double x_offset) {
      int base_index = 0;
      if(!check_cgns(
           cg_base_write(file_index, base_name, 2, physical_dim, &base_index),
           "failed to write unsupported CGNS physical-dimension base"
         )) {
        return false;
      }

      const cgsize_t zone_size[3] = {3, 1, 0};
      int zone_index = 0;
      if(!check_cgns(
           cg_zone_write(
             file_index,
             base_index,
             "surface_zone",
             zone_size,
             CGNS_ENUMV(Unstructured),
             &zone_index
           ),
           "failed to write unsupported CGNS physical-dimension zone"
         )) {
        return false;
      }

      const double x[3] = {x_offset + 0.0, x_offset + 1.0, x_offset + 0.0};
      const double y[3] = {0.0, 0.0, 1.0};
      const double z[3] = {0.0, 0.0, 0.0};
      int coordinate_index = 0;
      if(!check_cgns(
           cg_coord_write(
             file_index,
             base_index,
             zone_index,
             CGNS_ENUMV(RealDouble),
             "CoordinateX",
             x,
             &coordinate_index
           ),
           "failed to write unsupported CGNS physical-dimension CoordinateX"
         ) ||
         !check_cgns(
           cg_coord_write(
             file_index,
             base_index,
             zone_index,
             CGNS_ENUMV(RealDouble),
             "CoordinateY",
             y,
             &coordinate_index
           ),
           "failed to write unsupported CGNS physical-dimension CoordinateY"
         )) {
        return false;
      }
      if(physical_dim == 3 &&
         !check_cgns(
           cg_coord_write(
             file_index,
             base_index,
             zone_index,
             CGNS_ENUMV(RealDouble),
             "CoordinateZ",
             z,
             &coordinate_index
           ),
           "failed to write unsupported CGNS physical-dimension CoordinateZ"
         )) {
        return false;
      }

      const cgsize_t surface_faces[3] = {1, 2, 3};
      int section_index = 0;
      return check_cgns(
        cg_section_write(
          file_index,
          base_index,
          zone_index,
          "surface",
          CGNS_ENUMV(TRI_3),
          1,
          1,
          0,
          surface_faces,
          &section_index
        ),
        "failed to write unsupported CGNS physical-dimension section"
      );
    };

  if(!write_surface_base("sqmesh_base_a", 2, 0.0) ||
     !write_surface_base("sqmesh_base_b", 3, 2.0)) {
    static_cast<void>(cg_close(file_index));
    return false;
  }

  return check_cgns(
    cg_close(file_index),
    "failed to close unsupported CGNS physical-dimension fixture"
  );
}

bool write_unsupported_cgns_fixture(const std::filesystem::path &path)
{
  int file_index = 0;
  if(!check_cgns(cg_open(path.string().c_str(), CG_MODE_WRITE, &file_index), "failed to open unsupported CGNS fixture")) {
    return false;
  }

  int base_index = 0;
  if(!check_cgns(
       cg_base_write(file_index, "sqmesh_base_a", 3, 3, &base_index),
       "failed to write unsupported CGNS base_a"
     )) {
    static_cast<void>(cg_close(file_index));
    return false;
  }

  const cgsize_t volume_zone_size[3] = {4, 1, 0};
  int zone_index = 0;
  if(!check_cgns(
       cg_zone_write(
         file_index,
         base_index,
         "volume_zone",
         volume_zone_size,
         CGNS_ENUMV(Unstructured),
         &zone_index
       ),
       "failed to write unsupported CGNS volume zone"
     )) {
    static_cast<void>(cg_close(file_index));
    return false;
  }

  const double volume_x[4] = {0.0, 1.0, 0.0, 0.0};
  const double volume_y[4] = {0.0, 0.0, 1.0, 0.0};
  const double volume_z[4] = {0.0, 0.0, 0.0, 1.0};
  int coordinate_index = 0;
  if(!check_cgns(
       cg_coord_write(file_index, base_index, zone_index, CGNS_ENUMV(RealDouble), "CoordinateX", volume_x, &coordinate_index),
       "failed to write unsupported CGNS volume CoordinateX"
     ) ||
     !check_cgns(
       cg_coord_write(file_index, base_index, zone_index, CGNS_ENUMV(RealDouble), "CoordinateY", volume_y, &coordinate_index),
       "failed to write unsupported CGNS volume CoordinateY"
     ) ||
     !check_cgns(
       cg_coord_write(file_index, base_index, zone_index, CGNS_ENUMV(RealDouble), "CoordinateZ", volume_z, &coordinate_index),
       "failed to write unsupported CGNS volume CoordinateZ"
     )) {
    static_cast<void>(cg_close(file_index));
    return false;
  }

  int section_index = 0;
  const cgsize_t volume_cells[4] = {1, 2, 3, 4};
  if(!check_cgns(
       cg_section_write(file_index, base_index, zone_index, "volume", CGNS_ENUMV(TETRA_4), 1, 1, 0, volume_cells, &section_index),
       "failed to write unsupported CGNS volume section"
     )) {
    static_cast<void>(cg_close(file_index));
    return false;
  }

  if(!check_cgns(
       cg_base_write(file_index, "sqmesh_base_b", 2, 3, &base_index),
       "failed to write unsupported CGNS base_b"
     )) {
    static_cast<void>(cg_close(file_index));
    return false;
  }

  const cgsize_t surface_zone_size[3] = {3, 1, 0};
  if(!check_cgns(
       cg_zone_write(
         file_index,
         base_index,
         "surface_zone",
         surface_zone_size,
         CGNS_ENUMV(Unstructured),
         &zone_index
       ),
       "failed to write unsupported CGNS surface zone"
     )) {
    static_cast<void>(cg_close(file_index));
    return false;
  }

  const double surface_x[3] = {2.0, 3.0, 2.0};
  const double surface_y[3] = {0.0, 0.0, 1.0};
  const double surface_z[3] = {0.0, 0.0, 0.0};
  if(!check_cgns(
       cg_coord_write(file_index, base_index, zone_index, CGNS_ENUMV(RealDouble), "CoordinateX", surface_x, &coordinate_index),
       "failed to write unsupported CGNS surface CoordinateX"
     ) ||
     !check_cgns(
       cg_coord_write(file_index, base_index, zone_index, CGNS_ENUMV(RealDouble), "CoordinateY", surface_y, &coordinate_index),
       "failed to write unsupported CGNS surface CoordinateY"
     ) ||
     !check_cgns(
       cg_coord_write(file_index, base_index, zone_index, CGNS_ENUMV(RealDouble), "CoordinateZ", surface_z, &coordinate_index),
       "failed to write unsupported CGNS surface CoordinateZ"
     )) {
    static_cast<void>(cg_close(file_index));
    return false;
  }

  const cgsize_t surface_faces[3] = {1, 2, 3};
  if(!check_cgns(
       cg_section_write(file_index, base_index, zone_index, "surface", CGNS_ENUMV(TRI_3), 1, 1, 0, surface_faces, &section_index),
       "failed to write unsupported CGNS surface section"
     )) {
    static_cast<void>(cg_close(file_index));
    return false;
  }

  return check_cgns(cg_close(file_index), "failed to close unsupported CGNS fixture");
}

bool verify_exported_cgns_file(
  const std::filesystem::path &path,
  const char *expected_zone_name,
  const std::vector<CgnsSectionExpectation> &expected_sections,
  const std::vector<CgnsBoundaryExpectation> &expected_boundaries
)
{
  const auto verify_exported_cgns_zone =
    [&](int file_index,
        int zone_index,
        const CgnsZoneExpectation &expected_zone) {
      char zone_name[33] = {};
      cgsize_t zone_size[3] = {0, 0, 0};
      if(!check_cgns(
           cg_zone_read(file_index, 1, zone_index, zone_name, zone_size),
           "failed to read exported CGNS zone"
         )) {
        return false;
      }
      if(!expect(
           std::string_view(zone_name) == std::string_view(expected_zone.name),
           "export_cgns should preserve the CGNS zone name"
         )) {
        return false;
      }

      int section_count = 0;
      if(!check_cgns(
           cg_nsections(file_index, 1, zone_index, &section_count),
           "failed to read exported CGNS sections"
         )) {
        return false;
      }
      if(!expect(
           section_count == static_cast<int>(expected_zone.sections.size()),
           "export_cgns should write the expected number of supported sections"
         )) {
        return false;
      }

      for(int section = 1; section <= section_count; ++section) {
        const auto &expected_section =
          expected_zone.sections[static_cast<std::size_t>(section - 1)];
        char section_name[33] = {};
        CGNS_ENUMT(ElementType_t) type = CGNS_ENUMV(ElementTypeNull);
        cgsize_t start = 0;
        cgsize_t end = 0;
        int boundary_count = 0;
        int parent_flag = 0;
        if(!check_cgns(
             cg_section_read(
               file_index,
               1,
               zone_index,
               section,
               section_name,
               &type,
               &start,
               &end,
               &boundary_count,
               &parent_flag
             ),
             "failed to read exported CGNS section"
           )) {
          return false;
        }
        if(!expect(type == expected_section.type, "export_cgns should preserve section element types")) {
          return false;
        }
        if(!expect(
             std::string_view(section_name) == std::string_view(expected_section.name),
             "export_cgns should preserve section names"
           )) {
          return false;
        }
        if(!expect(
             start == expected_section.start && end == expected_section.end,
             "export_cgns should preserve section element ranges"
           )) {
          return false;
        }
        if(!check_cgns(
             cg_goto(file_index, 1, "Zone_t", zone_index, "Elements_t", section, "end"),
             "failed to access exported CGNS section metadata"
           )) {
          return false;
        }
        int descriptor_count = 0;
        if(!check_cgns(cg_ndescriptors(&descriptor_count), "failed to read exported CGNS descriptor count")) {
          return false;
        }
        bool found_zone_id = false;
        for(int descriptor_index = 1; descriptor_index <= descriptor_count; ++descriptor_index) {
          char descriptor_name[33] = {};
          char *descriptor_text = nullptr;
          if(!check_cgns(
               cg_descriptor_read(descriptor_index, descriptor_name, &descriptor_text),
               "failed to read exported CGNS descriptor"
             )) {
            return false;
          }

          const std::string descriptor_value =
            descriptor_text == nullptr ? std::string() : std::string(descriptor_text);
          if(descriptor_text != nullptr) {
            static_cast<void>(cg_free(descriptor_text));
          }

          if(std::string_view(descriptor_name) != std::string_view("SQMeshZoneId")) {
            continue;
          }

          found_zone_id = true;
          if(!expect(
               descriptor_value == expected_section.zone_id,
               "export_cgns should preserve section zone-id descriptors"
             )) {
            return false;
          }
        }
        if(!expect(found_zone_id, "export_cgns should write an SQMeshZoneId descriptor for each section")) {
          return false;
        }
      }

      int boundary_count = 0;
      if(!check_cgns(
           cg_nbocos(file_index, 1, zone_index, &boundary_count),
           "failed to read exported CGNS boundary-condition count"
         )) {
        return false;
      }
      if(!expect(
           boundary_count == static_cast<int>(expected_zone.boundaries.size()),
           "export_cgns should write the expected number of boundary-condition nodes"
         )) {
        std::fprintf(
          stderr,
          "mesh_io_smoke: exported CGNS boundary count actual=%d expected=%zu\n",
          boundary_count,
          expected_zone.boundaries.size()
        );
        return false;
      }

      for(int boundary = 1; boundary <= boundary_count; ++boundary) {
        const auto &expected_boundary =
          expected_zone.boundaries[static_cast<std::size_t>(boundary - 1)];
        char boundary_name[33] = {};
        CGNS_ENUMT(BCType_t) boundary_type = CGNS_ENUMV(BCTypeNull);
        CGNS_ENUMT(PointSetType_t) point_set_type = CGNS_ENUMV(PointSetTypeNull);
        cgsize_t point_count = 0;
        int normal_index[3] = {0, 0, 0};
        cgsize_t normal_list_size = 0;
        CGNS_ENUMT(DataType_t) normal_data_type = CGNS_ENUMV(DataTypeNull);
        int dataset_count = 0;
        if(!check_cgns(
             cg_boco_info(
               file_index,
               1,
               zone_index,
               boundary,
               boundary_name,
               &boundary_type,
               &point_set_type,
               &point_count,
               normal_index,
               &normal_list_size,
               &normal_data_type,
               &dataset_count
             ),
             "failed to read exported CGNS boundary-condition metadata"
           )) {
          return false;
        }
        static_cast<void>(normal_index);
        static_cast<void>(normal_list_size);
        static_cast<void>(normal_data_type);
        static_cast<void>(dataset_count);

        if(!expect(
             std::string_view(boundary_name) == std::string_view(expected_boundary.name),
             "export_cgns should preserve boundary-condition names"
           ) ||
           !expect(
             boundary_type == expected_boundary.type,
             "export_cgns should preserve boundary-condition types"
           ) ||
           !expect(
             point_set_type == expected_boundary.point_set_type,
             "export_cgns should preserve boundary-condition point-set types"
           )) {
          return false;
        }

        CGNS_ENUMT(GridLocation_t) location = CGNS_ENUMV(GridLocationNull);
        if(!check_cgns(
             cg_boco_gridlocation_read(file_index, 1, zone_index, boundary, &location),
             "failed to read exported CGNS boundary-condition grid location"
           )) {
          return false;
        }
        if(!expect(location == expected_boundary.location, "export_cgns should preserve boundary-condition grid locations")) {
          return false;
        }

        if(!expect(
             point_count == static_cast<cgsize_t>(expected_boundary.points.size()),
             "export_cgns should preserve boundary-condition point counts"
           )) {
          return false;
        }

        std::vector<cgsize_t> points(expected_boundary.points.size(), 0);
        if(!check_cgns(
             cg_boco_read(file_index, 1, zone_index, boundary, points.data(), nullptr),
             "failed to read exported CGNS boundary-condition points"
           )) {
          return false;
        }
        if(!expect(points == expected_boundary.points, "export_cgns should preserve boundary-condition point payloads")) {
          return false;
        }
      }

      int connection_count = 0;
      if(!check_cgns(
           cg_nconns(file_index, 1, zone_index, &connection_count),
           "failed to read exported CGNS interface-connectivity count"
         )) {
        return false;
      }
      if(!expect(
           connection_count == static_cast<int>(expected_zone.connections.size()),
           "export_cgns should write the expected number of interface-connectivity nodes"
         )) {
        return false;
      }

      for(int connection = 1; connection <= connection_count; ++connection) {
        const auto &expected_connection =
          expected_zone.connections[static_cast<std::size_t>(connection - 1)];
        char connection_name[33] = {};
        CGNS_ENUMT(GridLocation_t) location = CGNS_ENUMV(GridLocationNull);
        CGNS_ENUMT(GridConnectivityType_t) connection_type =
          CGNS_ENUMV(GridConnectivityTypeNull);
        CGNS_ENUMT(PointSetType_t) point_set_type = CGNS_ENUMV(PointSetTypeNull);
        cgsize_t point_count = 0;
        char donor_name[33] = {};
        CGNS_ENUMT(ZoneType_t) donor_zone_type = CGNS_ENUMV(ZoneTypeNull);
        CGNS_ENUMT(PointSetType_t) donor_point_set_type =
          CGNS_ENUMV(PointSetTypeNull);
        CGNS_ENUMT(DataType_t) donor_data_type = CGNS_ENUMV(DataTypeNull);
        cgsize_t donor_point_count = 0;
        if(!check_cgns(
             cg_conn_info(
               file_index,
               1,
               zone_index,
               connection,
               connection_name,
               &location,
               &connection_type,
               &point_set_type,
               &point_count,
               donor_name,
               &donor_zone_type,
               &donor_point_set_type,
               &donor_data_type,
               &donor_point_count
             ),
             "failed to read exported CGNS interface-connectivity metadata"
           )) {
          return false;
        }

        if(!expect(
             std::string_view(connection_name) ==
               std::string_view(expected_connection.name),
             "export_cgns should preserve interface-connectivity names"
           ) ||
           !expect(
             location == expected_connection.location,
             "export_cgns should preserve interface-connectivity locations"
           ) ||
           !expect(
             connection_type == expected_connection.type,
             "export_cgns should preserve interface-connectivity types"
           ) ||
           !expect(
             point_set_type == expected_connection.point_set_type,
             "export_cgns should preserve interface-connectivity point-set types"
           ) ||
           !expect(
             std::string_view(donor_name) ==
               std::string_view(expected_connection.donor_name),
             "export_cgns should preserve interface-connectivity donor zone names"
           ) ||
           !expect(
             donor_zone_type == expected_connection.donor_zone_type,
             "export_cgns should preserve interface-connectivity donor zone types"
           ) ||
           !expect(
             donor_point_set_type == expected_connection.donor_point_set_type,
             "export_cgns should preserve interface-connectivity donor point-set types"
           ) ||
           !expect(
             donor_data_type == expected_connection.donor_data_type,
             "export_cgns should preserve interface-connectivity donor data types"
           )) {
          return false;
        }

        if(!expect(
             point_count ==
               static_cast<cgsize_t>(expected_connection.points.size()),
             "export_cgns should preserve interface-connectivity point counts"
           ) ||
           !expect(
             donor_point_count ==
               static_cast<cgsize_t>(expected_connection.donor_points.size()),
             "export_cgns should preserve interface-connectivity donor point counts"
           )) {
          return false;
        }

        std::vector<cgsize_t> points(expected_connection.points.size(), 0);
        std::vector<cgsize_t> donor_points(
          expected_connection.donor_points.size(),
          0
        );
        if(!check_cgns(
             cg_conn_read(
               file_index,
               1,
               zone_index,
               connection,
               points.data(),
               donor_data_type,
               donor_points.data()
             ),
             "failed to read exported CGNS interface-connectivity point lists"
           )) {
          return false;
        }
        if(!expect(
             points == expected_connection.points,
             "export_cgns should preserve interface-connectivity point payloads"
           ) ||
           !expect(
             donor_points == expected_connection.donor_points,
             "export_cgns should preserve interface-connectivity donor payloads"
           )) {
          return false;
        }
      }

      return true;
    };

  int file_index = 0;
  if(!check_cgns(cg_open(path.string().c_str(), CG_MODE_READ, &file_index), "failed to reopen exported CGNS file")) {
    return false;
  }

  int base_count = 0;
  if(!check_cgns(cg_nbases(file_index, &base_count), "failed to read exported CGNS base count")) {
    static_cast<void>(cg_close(file_index));
    return false;
  }
  if(!expect(base_count == 1, "export_cgns should write exactly one base")) {
    static_cast<void>(cg_close(file_index));
    return false;
  }

  int zone_count = 0;
  if(!check_cgns(cg_nzones(file_index, 1, &zone_count), "failed to read exported CGNS zone count")) {
    static_cast<void>(cg_close(file_index));
    return false;
  }
  if(!expect(zone_count == 1, "export_cgns should write exactly one zone")) {
    static_cast<void>(cg_close(file_index));
    return false;
  }

  const bool verified = verify_exported_cgns_zone(
    file_index,
    1,
    {expected_zone_name, expected_sections, expected_boundaries}
  );
  const bool closed = check_cgns(cg_close(file_index), "failed to close exported CGNS file");
  return verified && closed;
}

bool verify_exported_multi_zone_cgns_file(
  const std::filesystem::path &path,
  const std::vector<CgnsZoneExpectation> &expected_zones
)
{
  auto verify_exported_cgns_zone =
    [&](int file_index,
        int zone_index,
        const CgnsZoneExpectation &expected_zone) {
      char zone_name[33] = {};
      cgsize_t zone_size[3] = {0, 0, 0};
      if(!check_cgns(
           cg_zone_read(file_index, 1, zone_index, zone_name, zone_size),
           "failed to read exported CGNS zone"
         )) {
        return false;
      }
      if(!expect(
           std::string_view(zone_name) == std::string_view(expected_zone.name),
           "export_cgns should preserve the CGNS zone name"
         )) {
        return false;
      }

      int section_count = 0;
      if(!check_cgns(
           cg_nsections(file_index, 1, zone_index, &section_count),
           "failed to read exported CGNS sections"
         )) {
        return false;
      }
      if(!expect(
           section_count == static_cast<int>(expected_zone.sections.size()),
           "export_cgns should write the expected number of supported sections"
         )) {
        return false;
      }

      for(int section = 1; section <= section_count; ++section) {
        const auto &expected_section =
          expected_zone.sections[static_cast<std::size_t>(section - 1)];
        char section_name[33] = {};
        CGNS_ENUMT(ElementType_t) type = CGNS_ENUMV(ElementTypeNull);
        cgsize_t start = 0;
        cgsize_t end = 0;
        int boundary_count = 0;
        int parent_flag = 0;
        if(!check_cgns(
             cg_section_read(
               file_index,
               1,
               zone_index,
               section,
               section_name,
               &type,
               &start,
               &end,
               &boundary_count,
               &parent_flag
             ),
             "failed to read exported CGNS section"
           )) {
          return false;
        }
        if(!expect(type == expected_section.type, "export_cgns should preserve section element types")) {
          return false;
        }
        if(!expect(
             std::string_view(section_name) == std::string_view(expected_section.name),
             "export_cgns should preserve section names"
           )) {
          return false;
        }
        if(!expect(
             start == expected_section.start && end == expected_section.end,
             "export_cgns should preserve section element ranges"
           )) {
          return false;
        }
        if(!check_cgns(
             cg_goto(file_index, 1, "Zone_t", zone_index, "Elements_t", section, "end"),
             "failed to access exported CGNS section metadata"
           )) {
          return false;
        }
        int descriptor_count = 0;
        if(!check_cgns(cg_ndescriptors(&descriptor_count), "failed to read exported CGNS descriptor count")) {
          return false;
        }
        bool found_zone_id = false;
        for(int descriptor_index = 1; descriptor_index <= descriptor_count; ++descriptor_index) {
          char descriptor_name[33] = {};
          char *descriptor_text = nullptr;
          if(!check_cgns(
               cg_descriptor_read(descriptor_index, descriptor_name, &descriptor_text),
               "failed to read exported CGNS descriptor"
             )) {
            return false;
          }

          const std::string descriptor_value =
            descriptor_text == nullptr ? std::string() : std::string(descriptor_text);
          if(descriptor_text != nullptr) {
            static_cast<void>(cg_free(descriptor_text));
          }

          if(std::string_view(descriptor_name) != std::string_view("SQMeshZoneId")) {
            continue;
          }

          found_zone_id = true;
          if(!expect(
               descriptor_value == expected_section.zone_id,
               "export_cgns should preserve section zone-id descriptors"
             )) {
            return false;
          }
        }
        if(!expect(found_zone_id, "export_cgns should write an SQMeshZoneId descriptor for each section")) {
          return false;
        }
      }

      int boundary_count = 0;
      if(!check_cgns(
           cg_nbocos(file_index, 1, zone_index, &boundary_count),
           "failed to read exported CGNS boundary-condition count"
         )) {
        return false;
      }
      if(!expect(
           boundary_count == static_cast<int>(expected_zone.boundaries.size()),
           "export_cgns should write the expected number of boundary-condition nodes"
         )) {
        return false;
      }

      for(int boundary = 1; boundary <= boundary_count; ++boundary) {
        const auto &expected_boundary =
          expected_zone.boundaries[static_cast<std::size_t>(boundary - 1)];
        char boundary_name[33] = {};
        CGNS_ENUMT(BCType_t) boundary_type = CGNS_ENUMV(BCTypeNull);
        CGNS_ENUMT(PointSetType_t) point_set_type = CGNS_ENUMV(PointSetTypeNull);
        cgsize_t point_count = 0;
        int normal_index[3] = {0, 0, 0};
        cgsize_t normal_list_size = 0;
        CGNS_ENUMT(DataType_t) normal_data_type = CGNS_ENUMV(DataTypeNull);
        int dataset_count = 0;
        if(!check_cgns(
             cg_boco_info(
               file_index,
               1,
               zone_index,
               boundary,
               boundary_name,
               &boundary_type,
               &point_set_type,
               &point_count,
               normal_index,
               &normal_list_size,
               &normal_data_type,
               &dataset_count
             ),
             "failed to read exported CGNS boundary-condition metadata"
           )) {
          return false;
        }
        static_cast<void>(normal_index);
        static_cast<void>(normal_list_size);
        static_cast<void>(normal_data_type);
        static_cast<void>(dataset_count);

        if(!expect(
             std::string_view(boundary_name) == std::string_view(expected_boundary.name),
             "export_cgns should preserve boundary-condition names"
           ) ||
           !expect(
             boundary_type == expected_boundary.type,
             "export_cgns should preserve boundary-condition types"
           ) ||
           !expect(
             point_set_type == expected_boundary.point_set_type,
             "export_cgns should preserve boundary-condition point-set types"
           )) {
          return false;
        }

        CGNS_ENUMT(GridLocation_t) location = CGNS_ENUMV(GridLocationNull);
        if(!check_cgns(
             cg_boco_gridlocation_read(file_index, 1, zone_index, boundary, &location),
             "failed to read exported CGNS boundary-condition grid location"
           )) {
          return false;
        }
        if(!expect(location == expected_boundary.location, "export_cgns should preserve boundary-condition grid locations")) {
          return false;
        }

        if(!expect(
             point_count == static_cast<cgsize_t>(expected_boundary.points.size()),
             "export_cgns should preserve boundary-condition point counts"
           )) {
          return false;
        }

        std::vector<cgsize_t> points(expected_boundary.points.size(), 0);
        if(!check_cgns(
             cg_boco_read(file_index, 1, zone_index, boundary, points.data(), nullptr),
             "failed to read exported CGNS boundary-condition points"
           )) {
          return false;
        }
        if(!expect(points == expected_boundary.points, "export_cgns should preserve boundary-condition point payloads")) {
          return false;
        }
      }

      return true;
    };

  int file_index = 0;
  if(!check_cgns(cg_open(path.string().c_str(), CG_MODE_READ, &file_index), "failed to reopen exported CGNS file")) {
    return false;
  }

  int base_count = 0;
  if(!check_cgns(cg_nbases(file_index, &base_count), "failed to read exported CGNS base count")) {
    static_cast<void>(cg_close(file_index));
    return false;
  }
  if(!expect(base_count == 1, "export_cgns should write exactly one base")) {
    static_cast<void>(cg_close(file_index));
    return false;
  }

  int zone_count = 0;
  if(!check_cgns(cg_nzones(file_index, 1, &zone_count), "failed to read exported CGNS zone count")) {
    static_cast<void>(cg_close(file_index));
    return false;
  }
  if(!expect(
       zone_count == static_cast<int>(expected_zones.size()),
       "export_cgns should write the expected number of CGNS zones"
     )) {
    static_cast<void>(cg_close(file_index));
    return false;
  }

  for(int zone = 1; zone <= zone_count; ++zone) {
    if(!verify_exported_cgns_zone(
         file_index,
         zone,
         expected_zones[static_cast<std::size_t>(zone - 1)]
       )) {
      static_cast<void>(cg_close(file_index));
      return false;
    }
  }

  return check_cgns(cg_close(file_index), "failed to close exported CGNS file");
}

bool verify_cgns_domain(
  const sqmesh::mesh::Domain &domain,
  const char *expected_base_name
)
{
  const auto expected_base_index =
    std::string_view(expected_base_name).empty()
      ? sqmesh::mesh::invalid_index
      : 1U;
  if(!verify_tetra_domain(
       domain,
       "sqmesh_cgns_zone",
       "feature_edges",
       11U,
       "wall",
       21U,
       "volume",
       31U
     )) {
    return false;
  }

  const auto *edge_entity_group = find_entity_group(
    domain,
    sqmesh::mesh::EntityOrder::edge,
    "feature_edges"
  );
  const auto *face_entity_group = find_entity_group(
    domain,
    sqmesh::mesh::EntityOrder::face,
    "wall"
  );
  const auto *cell_entity_group = find_entity_group(
    domain,
    sqmesh::mesh::EntityOrder::cell,
    "volume"
  );
  if(!expect(
       edge_entity_group != nullptr && face_entity_group != nullptr && cell_entity_group != nullptr,
       "CGNS tetra fixture should expose deterministic edge, face, and cell entity_groups"
     )) {
    return false;
  }
  if(!verify_cgns_entity_group_import_info(
       *edge_entity_group,
       expected_base_index,
       expected_base_name,
       1U,
       "sqmesh_cgns_zone",
       "feature_edges",
       static_cast<std::int32_t>(CGNS_ENUMV(BCSymmetryPlane)),
       "CGNS tetra edge entity_group"
     )) {
    return false;
  }
  if(!verify_cgns_entity_group_import_info(
       *face_entity_group,
       expected_base_index,
       expected_base_name,
       1U,
       "sqmesh_cgns_zone",
       "wall",
       static_cast<std::int32_t>(CGNS_ENUMV(BCWall)),
       "CGNS tetra face entity_group"
     )) {
    return false;
  }
  if(!verify_cgns_entity_group_import_info(
       *cell_entity_group,
       expected_base_index,
       expected_base_name,
       1U,
       "sqmesh_cgns_zone",
       "volume",
       -1,
       "CGNS tetra cell entity_group"
     )) {
    return false;
  }
  return true;
}

bool verify_cgns_mixed_surface_domain(
  const sqmesh::mesh::Domain &domain,
  const char *expected_base_name
)
{
  const auto expected_base_index =
    std::string_view(expected_base_name).empty()
      ? sqmesh::mesh::invalid_index
      : 1U;
  const auto summary = domain.summary();
  if(!expect(summary.node_count == 5U, "mixed CGNS import should preserve five nodes")) {
    return false;
  }
  if(!expect(summary.edge_count == 5U, "mixed CGNS import should preserve five line elements")) {
    return false;
  }
  if(!expect(summary.face_count == 3U, "mixed CGNS import should triangulate one quad plus one triangle into three faces")) {
    return false;
  }
  if(!expect(summary.cell_count == 0U, "mixed CGNS surface import should remain surface-only")) {
    return false;
  }
  if(!expect(domain.entity_group_count() == 3U, "mixed CGNS surface import should build node, edge, and face entity_groups")) {
    return false;
  }
  if(!expect(
       domain.name() == std::string_view("sqmesh_cgns_mixed_zone"),
       "mixed CGNS import should preserve the zone name"
     )) {
    return false;
  }

  const auto *edge_entity_group = find_entity_group(
    domain,
    sqmesh::mesh::EntityOrder::edge,
    "outer_loop"
  );
  const auto *face_entity_group = find_entity_group(
    domain,
    sqmesh::mesh::EntityOrder::face,
    "surface_patch"
  );
  if(!expect(edge_entity_group != nullptr, "mixed CGNS import should preserve the edge boundary-condition name")) {
    return false;
  }
  if(!expect(face_entity_group != nullptr, "mixed CGNS import should preserve the face boundary-condition name")) {
    return false;
  }
  if(!expect(edge_entity_group->zone_id() == 41U, "mixed CGNS edge entity_group should preserve the section descriptor")) {
    return false;
  }
  if(!expect(face_entity_group->zone_id() == 41U, "mixed CGNS face entity_group should preserve the section descriptor")) {
    return false;
  }
  if(!expect(edge_entity_group->is_boundary(), "mixed CGNS boundary edges should remain boundary edges")) {
    return false;
  }
  if(!expect(face_entity_group->is_boundary(), "mixed CGNS surface faces should remain boundary faces")) {
    return false;
  }
  if(!expect(edge_entity_group->edges().size() == 5U, "mixed CGNS boundary edge entity_group should contain the five edge records")) {
    return false;
  }
  if(!expect(face_entity_group->faces().size() == 3U, "mixed CGNS face entity_group should contain the triangulated surface faces")) {
    return false;
  }
  if(!verify_cgns_entity_group_import_info(
       *edge_entity_group,
       expected_base_index,
       expected_base_name,
       1U,
       "sqmesh_cgns_mixed_zone",
       "outer_loop",
       static_cast<std::int32_t>(CGNS_ENUMV(BCSymmetryPlane)),
       "mixed CGNS edge entity_group"
     )) {
    return false;
  }
  if(!verify_cgns_entity_group_import_info(
       *face_entity_group,
       expected_base_index,
       expected_base_name,
       1U,
       "sqmesh_cgns_mixed_zone",
       "surface_patch",
       static_cast<std::int32_t>(CGNS_ENUMV(BCInflow)),
       "mixed CGNS face entity_group"
     )) {
    return false;
  }

  const auto first_edge = sqmesh::mesh::EntityRef {edge_entity_group->id(), 0U};
  if(!expect(
       sqmesh::mesh::is_valid(domain.adjacent_face(first_edge, sqmesh::mesh::FaceSide::left)),
       "mixed CGNS boundary edges should map back to an adjacent surface face"
     )) {
    return false;
  }
  if(!expect(
       !sqmesh::mesh::is_valid(domain.adjacent_face(first_edge, sqmesh::mesh::FaceSide::right)),
       "mixed CGNS boundary edges should not expose a right-adjacent face"
     )) {
    return false;
  }

  return true;
}

bool verify_cgns_multi_zone_domain(
  const sqmesh::mesh::Domain &domain,
  const char *expected_domain_name,
  const char *expected_base_name
)
{
  const auto expected_base_index =
    std::string_view(expected_base_name).empty()
      ? sqmesh::mesh::invalid_index
      : 1U;
  const auto summary = domain.summary();
  if(!expect(summary.node_count == 8U, "multi-zone CGNS import should preserve eight nodes across two tetra zones")) {
    return false;
  }
  if(!expect(summary.edge_count == 12U, "multi-zone CGNS import should preserve twelve line elements")) {
    return false;
  }
  if(!expect(summary.face_count == 8U, "multi-zone CGNS import should preserve eight triangle faces")) {
    return false;
  }
  if(!expect(summary.cell_count == 2U, "multi-zone CGNS import should preserve two tetra cells")) {
    return false;
  }
  if(!expect(domain.entity_group_count() == 7U, "multi-zone CGNS import should keep separate edge, face, and cell entity_groups per zone")) {
    return false;
  }
  if(!expect(
       domain.name() == std::string_view(expected_domain_name),
       "multi-zone CGNS import should preserve the expected flattened domain name"
     )) {
    return false;
  }

  const auto verify_zone =
    [&](std::uint32_t zone_index,
        const char *zone_name,
        std::int32_t edge_bc_type_value,
        std::int32_t face_bc_type_value) {
    const auto prefix = std::string("z") + std::to_string(zone_index) + "_" + zone_name;
    const auto edge_name = std::string(prefix) + "__feature_edges";
    const auto face_name = std::string(prefix) + "__wall";
    const auto cell_name = std::string(prefix) + "__volume";

    const auto *edge_entity_group = find_entity_group(
      domain,
      sqmesh::mesh::EntityOrder::edge,
      edge_name.c_str()
    );
    const auto *face_entity_group = find_entity_group(
      domain,
      sqmesh::mesh::EntityOrder::face,
      face_name.c_str()
    );
    const auto *cell_entity_group = find_entity_group(
      domain,
      sqmesh::mesh::EntityOrder::cell,
      cell_name.c_str()
    );

    if(!expect(edge_entity_group != nullptr, "multi-zone CGNS import should preserve zone-qualified edge boundary names")) {
      return false;
    }
    if(!expect(face_entity_group != nullptr, "multi-zone CGNS import should preserve zone-qualified face boundary names")) {
      return false;
    }
    if(!expect(cell_entity_group != nullptr, "multi-zone CGNS import should preserve zone-qualified cell section names")) {
      return false;
    }
    if(!verify_cgns_entity_group_import_info(
         *edge_entity_group,
         expected_base_index,
         expected_base_name,
         zone_index,
         zone_name,
         "feature_edges",
         edge_bc_type_value,
         "multi-zone CGNS edge entity_group"
       )) {
      return false;
    }
    if(!verify_cgns_entity_group_import_info(
         *face_entity_group,
         expected_base_index,
         expected_base_name,
         zone_index,
         zone_name,
         "wall",
         face_bc_type_value,
         "multi-zone CGNS face entity_group"
       )) {
      return false;
    }
    if(!verify_cgns_entity_group_import_info(
         *cell_entity_group,
         expected_base_index,
         expected_base_name,
         zone_index,
         zone_name,
         "volume",
         -1,
         "multi-zone CGNS cell entity_group"
       )) {
      return false;
    }
    if(!expect(edge_entity_group->zone_id() == 11U, "multi-zone CGNS edge entity_groups should preserve section descriptors")) {
      return false;
    }
    if(!expect(face_entity_group->zone_id() == 21U, "multi-zone CGNS face entity_groups should preserve section descriptors")) {
      return false;
    }
    if(!expect(cell_entity_group->zone_id() == 31U, "multi-zone CGNS cell entity_groups should preserve section descriptors")) {
      return false;
    }
    if(!expect(!edge_entity_group->is_boundary(), "multi-zone CGNS edge entity_groups should still keep their two-face tetra adjacency")) {
      return false;
    }
    if(!expect(face_entity_group->is_boundary(), "multi-zone CGNS face entity_groups should remain boundary-only")) {
      return false;
    }
    if(!expect(!cell_entity_group->is_boundary(), "multi-zone CGNS cell entity_groups should stay interior")) {
      return false;
    }
    if(!expect(edge_entity_group->edges().size() == 6U, "each multi-zone CGNS edge entity_group should contain six BAR_2 records")) {
      return false;
    }
    if(!expect(face_entity_group->faces().size() == 4U, "each multi-zone CGNS face entity_group should contain four TRI_3 records")) {
      return false;
    }
    if(!expect(cell_entity_group->cells().size() == 1U, "each multi-zone CGNS cell entity_group should contain one TETRA_4 record")) {
      return false;
    }

    const auto first_edge = sqmesh::mesh::EntityRef {edge_entity_group->id(), 0U};
    if(!expect(
         sqmesh::mesh::is_valid(domain.adjacent_face(first_edge, sqmesh::mesh::FaceSide::left)),
         "multi-zone CGNS boundary edges should still map to adjacent faces"
       )) {
      return false;
    }
    if(!expect(
         sqmesh::mesh::is_valid(domain.adjacent_face(first_edge, sqmesh::mesh::FaceSide::right)),
         "multi-zone CGNS edge entity_groups should retain their second adjacent face inside each tetra zone"
       )) {
      return false;
    }

    return true;
  };

  return verify_zone(
           1U,
           "zone_a",
           static_cast<std::int32_t>(CGNS_ENUMV(BCSymmetryPlane)),
           static_cast<std::int32_t>(CGNS_ENUMV(BCWall))
         ) &&
         verify_zone(
           2U,
           "zone_b",
           static_cast<std::int32_t>(CGNS_ENUMV(BCOutflow)),
           static_cast<std::int32_t>(CGNS_ENUMV(BCWallInviscid))
         );
}

bool verify_cgns_multi_base_domain(
  const sqmesh::mesh::Domain &domain,
  const char *expected_domain_name
)
{
  const auto summary = domain.summary();
  if(!expect(summary.node_count == 8U, "multi-base CGNS import should preserve eight nodes across two tetra zones")) {
    return false;
  }
  if(!expect(summary.edge_count == 12U, "multi-base CGNS import should preserve twelve line elements")) {
    return false;
  }
  if(!expect(summary.face_count == 8U, "multi-base CGNS import should preserve eight triangle faces")) {
    return false;
  }
  if(!expect(summary.cell_count == 2U, "multi-base CGNS import should preserve two tetra cells")) {
    return false;
  }
  if(!expect(domain.entity_group_count() == 7U, "multi-base CGNS import should keep separate edge, face, and cell entity_groups per base")) {
    return false;
  }
  if(!expect(
       domain.name() == std::string_view(expected_domain_name),
       "multi-base CGNS import should preserve the expected flattened domain name"
     )) {
    return false;
  }

  const auto verify_base =
    [&](std::uint32_t base_index,
        const char *base_name,
        std::int32_t edge_bc_type_value,
        std::int32_t face_bc_type_value) {
    const auto prefix =
      std::string("b") + std::to_string(base_index) + "_" + base_name +
      "__z1_shared_zone";
    const auto edge_name = prefix + "__feature_edges";
    const auto face_name = prefix + "__wall";
    const auto cell_name = prefix + "__volume";

    const auto *edge_entity_group = find_entity_group(
      domain,
      sqmesh::mesh::EntityOrder::edge,
      edge_name.c_str()
    );
    const auto *face_entity_group = find_entity_group(
      domain,
      sqmesh::mesh::EntityOrder::face,
      face_name.c_str()
    );
    const auto *cell_entity_group = find_entity_group(
      domain,
      sqmesh::mesh::EntityOrder::cell,
      cell_name.c_str()
    );

    if(!expect(edge_entity_group != nullptr, "multi-base CGNS import should preserve base-qualified edge boundary names")) {
      return false;
    }
    if(!expect(face_entity_group != nullptr, "multi-base CGNS import should preserve base-qualified face boundary names")) {
      return false;
    }
    if(!expect(cell_entity_group != nullptr, "multi-base CGNS import should preserve base-qualified cell section names")) {
      return false;
    }
    if(!verify_cgns_entity_group_import_info(
         *edge_entity_group,
         base_index,
         base_name,
         1U,
         "shared_zone",
         "feature_edges",
         edge_bc_type_value,
         "multi-base CGNS edge entity_group"
       )) {
      return false;
    }
    if(!verify_cgns_entity_group_import_info(
         *face_entity_group,
         base_index,
         base_name,
         1U,
         "shared_zone",
         "wall",
         face_bc_type_value,
         "multi-base CGNS face entity_group"
       )) {
      return false;
    }
    if(!verify_cgns_entity_group_import_info(
         *cell_entity_group,
         base_index,
         base_name,
         1U,
         "shared_zone",
         "volume",
         -1,
         "multi-base CGNS cell entity_group"
       )) {
      return false;
    }
    if(!expect(edge_entity_group->zone_id() == 11U, "multi-base CGNS edge entity_groups should preserve section descriptors")) {
      return false;
    }
    if(!expect(face_entity_group->zone_id() == 21U, "multi-base CGNS face entity_groups should preserve section descriptors")) {
      return false;
    }
    if(!expect(cell_entity_group->zone_id() == 31U, "multi-base CGNS cell entity_groups should preserve section descriptors")) {
      return false;
    }
    if(!expect(!edge_entity_group->is_boundary(), "multi-base CGNS edge entity_groups should still keep their two-face tetra adjacency")) {
      return false;
    }
    if(!expect(face_entity_group->is_boundary(), "multi-base CGNS face entity_groups should remain boundary-only")) {
      return false;
    }
    if(!expect(!cell_entity_group->is_boundary(), "multi-base CGNS cell entity_groups should stay interior")) {
      return false;
    }
    if(!expect(edge_entity_group->edges().size() == 6U, "each multi-base CGNS edge entity_group should contain six BAR_2 records")) {
      return false;
    }
    if(!expect(face_entity_group->faces().size() == 4U, "each multi-base CGNS face entity_group should contain four TRI_3 records")) {
      return false;
    }
    if(!expect(cell_entity_group->cells().size() == 1U, "each multi-base CGNS cell entity_group should contain one TETRA_4 record")) {
      return false;
    }

    const auto first_edge = sqmesh::mesh::EntityRef {edge_entity_group->id(), 0U};
    if(!expect(
         sqmesh::mesh::is_valid(domain.adjacent_face(first_edge, sqmesh::mesh::FaceSide::left)),
         "multi-base CGNS boundary edges should still map to adjacent faces"
       )) {
      return false;
    }
    if(!expect(
         sqmesh::mesh::is_valid(domain.adjacent_face(first_edge, sqmesh::mesh::FaceSide::right)),
         "multi-base CGNS edge entity_groups should retain their second adjacent face inside each tetra zone"
       )) {
      return false;
    }

    return true;
  };

  return verify_base(
           1U,
           "sqmesh_base_a",
           static_cast<std::int32_t>(CGNS_ENUMV(BCSymmetryPlane)),
           static_cast<std::int32_t>(CGNS_ENUMV(BCWall))
         ) &&
         verify_base(
           2U,
           "sqmesh_base_b",
           static_cast<std::int32_t>(CGNS_ENUMV(BCOutflow)),
           static_cast<std::int32_t>(CGNS_ENUMV(BCWallInviscid))
         );
}

#endif

} // namespace

int main()
{
  constexpr const char *kMshSample = R"MSH($MeshFormat
2.2 0 8
$EndMeshFormat
$PhysicalNames
3
1 11 "feature_edges"
2 21 "wall"
3 31 "volume"
$EndPhysicalNames
$Nodes
4
1 0 0 0
2 1 0 0
3 0 1 0
4 0 0 1
$EndNodes
$Elements
11
1 1 2 11 11 1 2
2 1 2 11 11 2 3
3 1 2 11 11 3 1
4 1 2 11 11 1 4
5 1 2 11 11 2 4
6 1 2 11 11 3 4
7 2 2 21 21 1 3 2
8 2 2 21 21 1 2 4
9 2 2 21 21 2 3 4
10 2 2 21 21 3 1 4
11 4 2 31 31 1 2 3 4
$EndElements
)MSH";

  constexpr const char *kMshInterfaceSample = R"MSH($MeshFormat
2.2 0 8
$EndMeshFormat
$PhysicalNames
3
2 42 "fluid_solid_interface"
3 101 "fluid_region"
3 202 "solid_region"
$EndPhysicalNames
$Nodes
5
1 0 0 0
2 1 0 0
3 0 1 0
4 0 0 1
5 1 1 1
$EndNodes
$Elements
3
1 2 2 42 42 2 3 4
2 4 2 101 101 1 2 3 4
3 4 2 202 202 5 2 3 4
$EndElements
)MSH";

  constexpr const char *kMsh41Sample = R"MSH($MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
3
1 11 "feature_edges"
2 21 "wall"
3 31 "volume"
$EndPhysicalNames
$Entities
0 1 1 1
11 0 0 0 1 1 1 1 11 0
21 0 0 0 1 1 1 1 21 0
31 0 0 0 1 1 1 1 31 0
$EndEntities
$Nodes
1 4 1 4
3 31 0 4
1
2
3
4
0 0 0
1 0 0
0 1 0
0 0 1
$EndNodes
$Elements
3 11 1 11
1 11 1 6
1 1 2
2 2 3
3 3 1
4 1 4
5 2 4
6 3 4
2 21 2 4
7 1 3 2
8 1 2 4
9 2 3 4
10 3 1 4
3 31 4 1
11 1 2 3 4
$EndElements
)MSH";

  constexpr const char *kMsh41PartitionedSample = R"MSH($MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
3
1 11 "feature_edges"
2 21 "wall"
3 31 "volume"
$EndPhysicalNames
$PartitionedEntities
2
1
999 2
0 1 1 1
111 1 11 1 1 0 0 0 1 1 1 1 11 0
121 2 21 1 1 0 0 0 1 1 1 1 21 0
131 3 31 1 2 0 0 0 1 1 1 1 31 0
$EndPartitionedEntities
$Nodes
1 4 1 4
3 131 1 4
1
2
3
4
0 0 0 0 0 0
1 0 0 1 0 0
0 1 0 0 1 0
0 0 1 0 0 1
$EndNodes
$Elements
3 11 1 11
1 111 1 6
1 1 2
2 2 3
3 3 1
4 1 4
5 2 4
6 3 4
2 121 2 4
7 1 3 2
8 1 2 4
9 2 3 4
10 3 1 4
3 131 4 1
11 1 2 3 4
$EndElements
$Periodic
1
2 121 121
0
2
1 1
2 2
$EndPeriodic
$GhostElements
1
11 7 1 9
$EndGhostElements
$Parametrizations
1 1
111 2
0 0 0 0
1 0 0 1
121 3 1
0 0 0 0 0 1 0 0 0 1 0
1 0 0 1 0 1 0 0 0 1 0
0 1 0 0 1 1 0 0 0 1 0
1 2 3
$EndParametrizations
)MSH";

  constexpr const char *kMsh41SharedPhysicalEntitiesSample = R"MSH($MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
1
2 21 "wall"
$EndPhysicalNames
$Entities
0 0 2 0
101 0 0 0 1 1 0 1 21 0
202 2 0 0 3 1 0 1 21 0
$EndEntities
$Nodes
2 6 1 6
2 101 0 3
1
2
3
0 0 0
1 0 0
0 1 0
2 202 0 3
4
5
6
2 0 0
3 0 0
2 1 0
$EndNodes
$Elements
2 2 1 2
2 101 2 1
1 1 2 3
2 202 2 1
2 4 5 6
$EndElements
)MSH";

  constexpr const char *kObjSample = R"OBJ(# SQMesh OBJ sample
o tri_patch
v 0 0 0
v 1 0 0
v 0 1 0
g feature_loop
l 1 2
l 2 3
l 3 1
g surface_patch
f 1 3 2
)OBJ";

  constexpr const char *kObjGroupedSample = R"OBJ(# SQMesh OBJ grouped sample
o grouped_patch
v 0 0 0
v 1 0 0
v 1 1 0
v 0 1 0
v 0.5 0.5 1
vt 0 0
vn 0 0 1
g boundary_loop
l -5 -4 -3 -2 -5
g panel_A seam
f 1/1/1 2/1/1 3/1/1 4/1/1 # quad with texture/normal refs
g cap
f -5 -3 -1
g
o loose_object
f 2 4 5
)OBJ";

  constexpr const char *kObjAttributeSample = R"OBJ(# SQMesh OBJ attribute-channel sample
mtllib base.mtl accent.mtl
o obj_attr_patch
g feature_edges
v 0 0 0
v 1 0 0
v 1 1 0
v 0 1 0
v 0 0 1
v 1 0 1
v 1 1 1
vt 0.0
vt 1.0 0.0
vt 1.0 1.0 0.0
vt 0.0 1.0
vt 0.5 0.5
vn 0 0 1
vn 0 0 -1
vp 0.25
vp 0.5 0.75
vp 0.1 0.2 1.0
usemtl wire
l 1 2 3
g faces
usemtl plain
f 1 2 3
usemtl tex_only
f 1/1 3/2 4/4
usemtl normal_only
f 1//1 4//1 5//1
usemtl full_negative
f -3/-3/-1 -2/-2/-1 -1/-1/-1
)OBJ";

  constexpr const char *kObjConcaveSample = R"OBJ(# SQMesh OBJ concave polygon sample
o concave_patch
v 3 1 0
v 1 1 0
v 1 3 0
v 0 3 0
v 0 0 0
v 3 0 0
g concave_face
f 1 2 3 4 5 6
)OBJ";

  constexpr const char *kObjCollinearPolygonSample = R"OBJ(# SQMesh OBJ consecutive-collinear polygon sample
o collinear_patch
v 0 0 0
v 1 0 0
v 2 0 0
v 2 1 0
v 0 1 0
f 1 2 3 4 5
)OBJ";

  constexpr const char *kObjUnsupportedPolygonSample = R"OBJ(# SQMesh OBJ unsupported polygon sample
o bow_tie_patch
v 0 0 0
v 2 0 0
v 0 2 0
v 2 1 0
v 0 1 0
f 1 2 3 4 5
)OBJ";

  constexpr const char *kNastranSample = R"NAS($ SQMesh NASTRAN sample
BEGIN BULK
MAT1,31,210000.,,0.3
PROD,11,31,1.0
PSHELL,21,31,1.0
PSOLID,31,31
GRID,1,,0.0,0.0,0.0,0,0,0
GRID,2,,1.0,0.0,0.0,0,0,0
GRID,3,,0.0,1.0,0.0,0,0,0
GRID,4,,0.0,0.0,1.0,0,0,0
CROD,1,11,1,2
CROD,2,11,2,3
CROD,3,11,3,1
CROD,4,11,1,4
CROD,5,11,2,4
CROD,6,11,3,4
CTRIA3,7,21,1,3,2
CTRIA3,8,21,1,2,4
CTRIA3,9,21,2,3,4
CTRIA3,10,21,3,1,4
CTETRA,11,31,1,2,3,4
ENDDATA
)NAS";

  constexpr const char *kNastranBroadSample = R"NAS($ SQMesh broadened NASTRAN sample
BEGIN BULK
MAT1,301,210000.,,0.3
PROD,101,301,1.0
PSHELL,201,301,1.0
GRID,1,,0.0,0.0,0.0,0,0,0
GRID,2,,1.0,0.0,0.0,0,0,0
GRID,3,,1.0,1.0,0.0,0,0,0
GRID,4,,0.0,1.0,0.0,0,0,0
GRID,5,,0.5,0.5,1.0,0,0,0
CBAR,1,101,1,2,0.0,0.0,1.0
CBEAM,2,101,2,3,0.0,0.0,1.0
CROD,3,101,3,4
CQUAD4,10,201,1,2,3,4,0.0
CTRIA3,11,201,1,3,5,0.0
ENDDATA
)NAS";

  constexpr const char *kNastranUnsupportedSample = R"NAS(BEGIN BULK
GRID,1,,0.0,0.0,0.0
GRID,2,,1.0,0.0,0.0
GRID,3,,0.0,1.0,0.0
GRID,4,,0.0,0.0,1.0
GRID,5,,1.0,1.0,0.0
GRID,6,,0.0,0.0,1.0
GRID,7,,1.0,0.0,1.0
GRID,8,,1.0,1.0,1.0
CHEXA,1,41,1,2,5,3,4,7,8,6
ENDDATA
)NAS";

  const auto workdir = std::filesystem::current_path();
  const auto msh_input = workdir / "sqmesh_mesh_io_input.msh";
  const auto msh_output = workdir / "sqmesh_mesh_io_roundtrip.msh";
  const auto msh_interface_input = workdir / "sqmesh_mesh_io_interface_input.msh";
  const auto msh_interface_output = workdir / "sqmesh_mesh_io_interface_roundtrip.msh";
  const auto msh_binary_input = workdir / "sqmesh_mesh_io_input_binary.msh";
  const auto msh_binary_output = workdir / "sqmesh_mesh_io_roundtrip_binary.msh";
  const auto msh41_input = workdir / "sqmesh_mesh_io_input_v41.msh";
  const auto msh41_output = workdir / "sqmesh_mesh_io_roundtrip_v41.msh";
  const auto msh41_binary_input = workdir / "sqmesh_mesh_io_input_v41_binary.msh";
  const auto msh41_binary_output = workdir / "sqmesh_mesh_io_roundtrip_v41_binary.msh";
  const auto msh41_binary_32bit_input = workdir / "sqmesh_mesh_io_input_v41_binary_32bit.msh";
  const auto msh41_binary_32bit_output =
    workdir / "sqmesh_mesh_io_roundtrip_v41_binary_32bit.msh";
  const auto msh41_partitioned_input = workdir / "sqmesh_mesh_io_input_v41_partitioned.msh";
  const auto msh41_partitioned_output = workdir / "sqmesh_mesh_io_roundtrip_v41_partitioned.msh";
  const auto msh41_partitioned_binary_input =
    workdir / "sqmesh_mesh_io_input_v41_partitioned_binary.msh";
  const auto msh41_partitioned_binary_output =
    workdir / "sqmesh_mesh_io_roundtrip_v41_partitioned_binary.msh";
  const auto msh41_shared_physical_input =
    workdir / "sqmesh_mesh_io_input_v41_shared_physical.msh";
  const auto msh41_shared_physical_output =
    workdir / "sqmesh_mesh_io_roundtrip_v41_shared_physical.msh";
  const auto msh41_shared_physical_binary_output =
    workdir / "sqmesh_mesh_io_roundtrip_v41_shared_physical_binary.msh";
  const auto msh_conflicting_names_output =
    workdir / "sqmesh_mesh_io_conflicting_physical_names.msh";
  const auto obj_input = workdir / "sqmesh_mesh_io_input.obj";
  const auto obj_output = workdir / "sqmesh_mesh_io_roundtrip.obj";
  const auto obj_attribute_input = workdir / "sqmesh_mesh_io_attribute_input.obj";
  const auto obj_attribute_output = workdir / "sqmesh_mesh_io_attribute_roundtrip.obj";
  const auto obj_concave_input = workdir / "sqmesh_mesh_io_concave_input.obj";
  const auto obj_concave_output = workdir / "sqmesh_mesh_io_concave_roundtrip.obj";
  const auto obj_collinear_polygon_input = workdir / "sqmesh_mesh_io_collinear_polygon.obj";
  const auto obj_grouped_input = workdir / "sqmesh_mesh_io_grouped_input.obj";
  const auto obj_grouped_output = workdir / "sqmesh_mesh_io_grouped_roundtrip.obj";
  const auto obj_unsupported_polygon_input = workdir / "sqmesh_mesh_io_unsupported_polygon.obj";
  const auto nastran_input = workdir / "sqmesh_mesh_io_input.bdf";
  const auto nastran_output = workdir / "sqmesh_mesh_io_roundtrip.bdf";
  const auto nastran_fixed_input = workdir / "sqmesh_mesh_io_small_field.dat";
  const auto nastran_broad_input = workdir / "sqmesh_mesh_io_broadened_input.bdf";
  const auto nastran_broad_output = workdir / "sqmesh_mesh_io_broadened_roundtrip.bdf";
  const auto nastran_broad_fixed_input = workdir / "sqmesh_mesh_io_broadened_small_field.dat";
  const auto nastran_broad_fixed_output = workdir / "sqmesh_mesh_io_broadened_roundtrip_small_field.dat";
  const auto nastran_broad_large_input = workdir / "sqmesh_mesh_io_broadened_large_field.dat";
  const auto nastran_broad_large_output = workdir / "sqmesh_mesh_io_broadened_roundtrip_large_field.dat";
  const auto nastran_unpaired_quad_output =
    workdir / "sqmesh_mesh_io_unpaired_cquad4_roundtrip.bdf";
  const auto nastran_unsupported = workdir / "sqmesh_mesh_io_unsupported.bdf";
  const auto cgns_input = workdir / "sqmesh_mesh_io_input.cgns";
  const auto cgns_output = workdir / "sqmesh_mesh_io_roundtrip.cgns";
  const auto cgns_mixed_input = workdir / "sqmesh_mesh_io_mixed_input.cgns";
  const auto cgns_mixed_output = workdir / "sqmesh_mesh_io_mixed_roundtrip.cgns";
  const auto cgns_multizone_input = workdir / "sqmesh_mesh_io_multizone_input.cgns";
  const auto cgns_multizone_output = workdir / "sqmesh_mesh_io_multizone_roundtrip.cgns";
  const auto cgns_multizone_preserved_output =
    workdir / "sqmesh_mesh_io_multizone_preserved_roundtrip.cgns";
  const auto cgns_multibase_input = workdir / "sqmesh_mesh_io_multibase_input.cgns";
  const auto cgns_multibase_output =
    workdir / "sqmesh_mesh_io_multibase_roundtrip.cgns";
  const auto cgns_multibase_preserved_output =
    workdir / "sqmesh_mesh_io_multibase_preserved_roundtrip.cgns";
  const auto cgns_multizone_interface_output =
    workdir / "sqmesh_mesh_io_multizone_interface_roundtrip.cgns";
  const auto cgns_surface_multizone_interface_output =
    workdir / "sqmesh_mesh_io_surface_multizone_interface_roundtrip.cgns";
  const auto cgns_multizone_patch_interface_output =
    workdir / "sqmesh_mesh_io_multizone_patch_interface_roundtrip.cgns";
  const auto cgns_surface_multizone_polyline_interface_output =
    workdir / "sqmesh_mesh_io_surface_multizone_polyline_interface_roundtrip.cgns";
  const auto cgns_multizone_interface_unsupported_output =
    workdir / "sqmesh_mesh_io_multizone_interface_unsupported_roundtrip.cgns";
  const auto cgns_unsupported_physical_dim =
    workdir / "sqmesh_mesh_io_unsupported_physical_dim.cgns";
  const auto cgns_unsupported = workdir / "sqmesh_mesh_io_unsupported.cgns";

  std::filesystem::remove(msh_input);
  std::filesystem::remove(msh_output);
  std::filesystem::remove(msh_interface_input);
  std::filesystem::remove(msh_interface_output);
  std::filesystem::remove(msh_binary_input);
  std::filesystem::remove(msh_binary_output);
  std::filesystem::remove(msh41_input);
  std::filesystem::remove(msh41_output);
  std::filesystem::remove(msh41_binary_input);
  std::filesystem::remove(msh41_binary_output);
  std::filesystem::remove(msh41_binary_32bit_input);
  std::filesystem::remove(msh41_binary_32bit_output);
  std::filesystem::remove(msh41_partitioned_input);
  std::filesystem::remove(msh41_partitioned_output);
  std::filesystem::remove(msh41_partitioned_binary_input);
  std::filesystem::remove(msh41_partitioned_binary_output);
  std::filesystem::remove(msh41_shared_physical_input);
  std::filesystem::remove(msh41_shared_physical_output);
  std::filesystem::remove(msh41_shared_physical_binary_output);
  std::filesystem::remove(msh_conflicting_names_output);
  std::filesystem::remove(obj_input);
  std::filesystem::remove(obj_output);
  std::filesystem::remove(obj_attribute_input);
  std::filesystem::remove(obj_attribute_output);
  std::filesystem::remove(obj_concave_input);
  std::filesystem::remove(obj_concave_output);
  std::filesystem::remove(obj_collinear_polygon_input);
  std::filesystem::remove(obj_grouped_input);
  std::filesystem::remove(obj_grouped_output);
  std::filesystem::remove(obj_unsupported_polygon_input);
  std::filesystem::remove(nastran_input);
  std::filesystem::remove(nastran_output);
  std::filesystem::remove(nastran_fixed_input);
  std::filesystem::remove(nastran_broad_input);
  std::filesystem::remove(nastran_broad_output);
  std::filesystem::remove(nastran_broad_fixed_input);
  std::filesystem::remove(nastran_broad_fixed_output);
  std::filesystem::remove(nastran_broad_large_input);
  std::filesystem::remove(nastran_broad_large_output);
  std::filesystem::remove(nastran_unpaired_quad_output);
  std::filesystem::remove(nastran_unsupported);
  std::filesystem::remove(cgns_input);
  std::filesystem::remove(cgns_output);
  std::filesystem::remove(cgns_mixed_input);
  std::filesystem::remove(cgns_mixed_output);
  std::filesystem::remove(cgns_multizone_input);
  std::filesystem::remove(cgns_multizone_output);
  std::filesystem::remove(cgns_multizone_preserved_output);
  std::filesystem::remove(cgns_multibase_input);
  std::filesystem::remove(cgns_multibase_output);
  std::filesystem::remove(cgns_multibase_preserved_output);
  std::filesystem::remove(cgns_multizone_interface_output);
  std::filesystem::remove(cgns_surface_multizone_interface_output);
  std::filesystem::remove(cgns_multizone_patch_interface_output);
  std::filesystem::remove(cgns_surface_multizone_polyline_interface_output);
  std::filesystem::remove(cgns_multizone_interface_unsupported_output);
  std::filesystem::remove(cgns_unsupported_physical_dim);
  std::filesystem::remove(cgns_unsupported);

  if(!write_text_file(msh_input, kMshSample) ||
     !write_text_file(msh_interface_input, kMshInterfaceSample) ||
     !write_msh22_binary_fixture(msh_binary_input) ||
     !write_text_file(msh41_input, kMsh41Sample) ||
     !write_msh41_binary_fixture(msh41_binary_input) ||
     !write_msh41_binary_fixture_32bit(msh41_binary_32bit_input) ||
     !write_text_file(msh41_partitioned_input, kMsh41PartitionedSample) ||
     !write_msh41_partitioned_binary_fixture(msh41_partitioned_binary_input) ||
     !write_text_file(msh41_shared_physical_input, kMsh41SharedPhysicalEntitiesSample) ||
     !write_text_file(obj_input, kObjSample) ||
     !write_text_file(obj_attribute_input, kObjAttributeSample) ||
     !write_text_file(obj_concave_input, kObjConcaveSample) ||
     !write_text_file(obj_collinear_polygon_input, kObjCollinearPolygonSample) ||
     !write_text_file(obj_grouped_input, kObjGroupedSample) ||
     !write_text_file(obj_unsupported_polygon_input, kObjUnsupportedPolygonSample) ||
     !write_text_file(nastran_input, kNastranSample) ||
     !write_text_file(nastran_fixed_input, make_nastran_small_field_sample()) ||
     !write_text_file(nastran_broad_input, kNastranBroadSample) ||
     !write_text_file(nastran_broad_fixed_input, make_nastran_small_field_broad_sample()) ||
     !write_text_file(nastran_broad_large_input, make_nastran_large_field_broad_sample()) ||
     !write_text_file(nastran_unsupported, kNastranUnsupportedSample)) {
    return EXIT_FAILURE;
  }

  sqmesh::base::ContextHandle context = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::base::initialize(context) == sqmesh::base::StatusCode::ok,
       "initialize should create a runtime context"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle msh_mesh = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_msh(msh_input.string(), msh_mesh, {}, context) ==
         sqmesh::base::StatusCode::ok,
       "import_msh should load the supported ASCII Gmsh 2.2 subset"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       msh_mesh != sqmesh::invalid_handle,
       "import_msh should return a runtime mesh handle"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshSummary msh_summary;
  sqmesh::mesh::Domain msh_domain;
  if(!expect(
       sqmesh::mesh::mesh_summary(msh_mesh, msh_summary, context) ==
         sqmesh::base::StatusCode::ok,
       "mesh_summary should resolve imported MSH meshes"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::mesh::domain_snapshot(msh_mesh, msh_domain, context) ==
         sqmesh::base::StatusCode::ok,
       "domain_snapshot should materialize imported MSH data into Domain/EntityGroup storage"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       msh_summary.node_count == 4U &&
         msh_summary.edge_count == 6U &&
         msh_summary.face_count == 4U &&
         msh_summary.cell_count == 1U,
       "mesh_summary should report the imported MSH counts"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_msh_domain(msh_domain)) {
    return EXIT_FAILURE;
  }

  if(!expect(
       sqmesh::mesh::export_msh(msh_mesh, msh_output.string(), {}, context) ==
         sqmesh::base::StatusCode::ok,
       "export_msh should write the supported ASCII Gmsh 2.2 subset"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle msh_roundtrip = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_msh(msh_output.string(), msh_roundtrip, {}, context) ==
         sqmesh::base::StatusCode::ok,
       "exported MSH files should be re-readable by import_msh"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::Domain msh_roundtrip_domain;
  if(!expect(
       sqmesh::mesh::domain_snapshot(msh_roundtrip, msh_roundtrip_domain, context) ==
         sqmesh::base::StatusCode::ok,
       "domain_snapshot should resolve the round-tripped MSH mesh"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_msh_domain(msh_roundtrip_domain)) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle msh_interface_mesh = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_msh(msh_interface_input.string(), msh_interface_mesh, {}, context) ==
         sqmesh::base::StatusCode::ok,
       "import_msh should load an explicit two-region shared-face interface case"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::Domain msh_interface_domain;
  if(!expect(
       sqmesh::mesh::domain_snapshot(msh_interface_mesh, msh_interface_domain, context) ==
         sqmesh::base::StatusCode::ok,
       "domain_snapshot should materialize multi-region MSH interface data into Domain/EntityGroup storage"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_msh_interface_domain(msh_interface_domain)) {
    return EXIT_FAILURE;
  }

  if(!expect(
       sqmesh::mesh::export_msh(msh_interface_mesh, msh_interface_output.string(), {}, context) ==
         sqmesh::base::StatusCode::ok,
       "export_msh should preserve the supported multi-region interface subset"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle msh_interface_roundtrip = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_msh(
         msh_interface_output.string(),
         msh_interface_roundtrip,
         {},
         context
       ) == sqmesh::base::StatusCode::ok,
       "exported multi-region MSH files should remain re-readable by import_msh"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::Domain msh_interface_roundtrip_domain;
  if(!expect(
       sqmesh::mesh::domain_snapshot(
         msh_interface_roundtrip,
         msh_interface_roundtrip_domain,
         context
       ) == sqmesh::base::StatusCode::ok,
       "domain_snapshot should resolve the round-tripped multi-region MSH mesh"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_msh_interface_domain(msh_interface_roundtrip_domain)) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle msh_binary_mesh = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_msh(msh_binary_input.string(), msh_binary_mesh, {}, context) ==
         sqmesh::base::StatusCode::ok,
       "import_msh should load the supported binary Gmsh 2.2 subset"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::Domain msh_binary_domain;
  if(!expect(
       sqmesh::mesh::domain_snapshot(msh_binary_mesh, msh_binary_domain, context) ==
         sqmesh::base::StatusCode::ok,
       "domain_snapshot should materialize imported binary Gmsh 2.2 data"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_msh_domain(msh_binary_domain)) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MshExportOptions msh_binary_export_options {};
  msh_binary_export_options.format_version = sqmesh::mesh::MshFormatVersion::gmsh22_binary;
  if(!expect(
       sqmesh::mesh::export_msh(
         msh_binary_mesh,
         msh_binary_output.string(),
         msh_binary_export_options,
         context
       ) == sqmesh::base::StatusCode::ok,
       "export_msh should write the supported binary Gmsh 2.2 subset"
     )) {
    return EXIT_FAILURE;
  }

  std::string msh_binary_output_bytes;
  if(!expect(
       read_binary_file(msh_binary_output, msh_binary_output_bytes),
       "mesh_io_smoke should be able to read the exported binary Gmsh 2.2 file"
     )) {
    return EXIT_FAILURE;
  }
  const auto msh22_binary_header =
    std::string("$MeshFormat\n2.2 1 ") + std::to_string(sizeof(double)) + "\n";
  if(!expect(
       msh_binary_output_bytes.find(msh22_binary_header) != std::string::npos &&
         msh_binary_output_bytes.find("$PhysicalNames\n") != std::string::npos,
       "exported binary Gmsh 2.2 output should advertise binary MeshFormat and PhysicalNames data"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle msh_binary_roundtrip = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_msh(msh_binary_output.string(), msh_binary_roundtrip, {}, context) ==
         sqmesh::base::StatusCode::ok,
       "exported binary Gmsh 2.2 files should be re-readable by import_msh"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::Domain msh_binary_roundtrip_domain;
  if(!expect(
       sqmesh::mesh::domain_snapshot(msh_binary_roundtrip, msh_binary_roundtrip_domain, context) ==
         sqmesh::base::StatusCode::ok,
       "domain_snapshot should resolve the round-tripped binary Gmsh 2.2 mesh"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_msh_domain(msh_binary_roundtrip_domain)) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle msh41_mesh = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_msh(msh41_input.string(), msh41_mesh, {}, context) ==
         sqmesh::base::StatusCode::ok,
       "import_msh should load the supported ASCII Gmsh 4.1 subset"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshSummary msh41_summary;
  sqmesh::mesh::Domain msh41_domain;
  if(!expect(
       sqmesh::mesh::mesh_summary(msh41_mesh, msh41_summary, context) ==
         sqmesh::base::StatusCode::ok,
       "mesh_summary should resolve imported ASCII Gmsh 4.1 meshes"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::mesh::domain_snapshot(msh41_mesh, msh41_domain, context) ==
         sqmesh::base::StatusCode::ok,
       "domain_snapshot should materialize imported ASCII Gmsh 4.1 data into Domain/EntityGroup storage"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       msh41_summary.node_count == 4U &&
         msh41_summary.edge_count == 6U &&
         msh41_summary.face_count == 4U &&
         msh41_summary.cell_count == 1U,
       "mesh_summary should report the imported ASCII Gmsh 4.1 counts"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_msh_domain(msh41_domain)) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MshExportOptions msh41_export_options {};
  msh41_export_options.format_version = sqmesh::mesh::MshFormatVersion::gmsh41_ascii;
  sqmesh::mesh::MshExportOptions msh41_binary_export_options {};
  msh41_binary_export_options.format_version =
    sqmesh::mesh::MshFormatVersion::gmsh41_binary;
  if(!expect(
       sqmesh::mesh::export_msh(msh41_mesh, msh41_output.string(), msh41_export_options, context) ==
         sqmesh::base::StatusCode::ok,
       "export_msh should write the supported ASCII Gmsh 4.1 subset"
     )) {
    return EXIT_FAILURE;
  }

  std::string msh41_output_text;
  if(!expect(
       read_text_file(msh41_output, msh41_output_text),
       "mesh_io_smoke should be able to read the exported Gmsh 4.1 file"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       msh41_output_text.find("4.1 0 8") != std::string::npos &&
         msh41_output_text.find("$Entities") != std::string::npos &&
         msh41_output_text.find("\"feature_edges\"") != std::string::npos,
       "exported Gmsh 4.1 output should contain MeshFormat 4.1, Entities, and PhysicalNames data"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle msh41_roundtrip = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_msh(msh41_output.string(), msh41_roundtrip, {}, context) ==
         sqmesh::base::StatusCode::ok,
       "exported ASCII Gmsh 4.1 files should be re-readable by import_msh"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::Domain msh41_roundtrip_domain;
  if(!expect(
       sqmesh::mesh::domain_snapshot(msh41_roundtrip, msh41_roundtrip_domain, context) ==
         sqmesh::base::StatusCode::ok,
       "domain_snapshot should resolve the round-tripped ASCII Gmsh 4.1 mesh"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_msh_domain(msh41_roundtrip_domain)) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle msh41_shared_physical_mesh = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_msh(
         msh41_shared_physical_input.string(),
         msh41_shared_physical_mesh,
         {},
         context
       ) == sqmesh::base::StatusCode::ok,
       "import_msh should preserve distinct Gmsh 4.1 entity blocks that share one physical group"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::Domain msh41_shared_physical_domain;
  if(!expect(
       sqmesh::mesh::domain_snapshot(
         msh41_shared_physical_mesh,
         msh41_shared_physical_domain,
         context
       ) == sqmesh::base::StatusCode::ok,
       "domain_snapshot should materialize shared-physical Gmsh 4.1 entity-block data"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_msh41_shared_physical_entity_domain(msh41_shared_physical_domain)) {
    return EXIT_FAILURE;
  }

  if(!expect(
       sqmesh::mesh::export_msh(
         msh41_shared_physical_mesh,
         msh41_shared_physical_output.string(),
         msh41_export_options,
         context
       ) == sqmesh::base::StatusCode::ok,
       "export_msh should keep separate Gmsh 4.1 entity blocks when they share one physical group"
     )) {
    return EXIT_FAILURE;
  }

  std::string msh41_shared_physical_output_text;
  if(!expect(
       read_text_file(
         msh41_shared_physical_output,
         msh41_shared_physical_output_text
       ),
       "mesh_io_smoke should be able to read the exported shared-physical Gmsh 4.1 file"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       msh41_shared_physical_output_text.find("$PhysicalNames\n1\n2 21 \"wall\"\n") !=
           std::string::npos &&
         msh41_shared_physical_output_text.find("$Entities\n0 0 2 0\n") != std::string::npos &&
         msh41_shared_physical_output_text.find("\n2 101 2 1\n") != std::string::npos &&
         msh41_shared_physical_output_text.find("\n2 202 2 1\n") != std::string::npos,
       "exported shared-physical Gmsh 4.1 output should retain one physical name with two distinct face entity blocks"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle msh41_shared_physical_roundtrip = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_msh(
         msh41_shared_physical_output.string(),
         msh41_shared_physical_roundtrip,
         {},
         context
       ) == sqmesh::base::StatusCode::ok,
       "exported shared-physical Gmsh 4.1 files should be re-readable by import_msh"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::Domain msh41_shared_physical_roundtrip_domain;
  if(!expect(
       sqmesh::mesh::domain_snapshot(
         msh41_shared_physical_roundtrip,
         msh41_shared_physical_roundtrip_domain,
         context
       ) == sqmesh::base::StatusCode::ok,
       "domain_snapshot should resolve the round-tripped shared-physical Gmsh 4.1 mesh"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_msh41_shared_physical_entity_domain(msh41_shared_physical_roundtrip_domain)) {
    return EXIT_FAILURE;
  }

  if(!expect(
       sqmesh::mesh::export_msh(
         msh41_shared_physical_mesh,
         msh41_shared_physical_binary_output.string(),
         msh41_binary_export_options,
         context
       ) == sqmesh::base::StatusCode::ok,
       "export_msh should also keep shared-physical Gmsh 4.1 entity blocks on the binary export path"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle msh41_shared_physical_binary_roundtrip =
    sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_msh(
         msh41_shared_physical_binary_output.string(),
         msh41_shared_physical_binary_roundtrip,
         {},
         context
       ) == sqmesh::base::StatusCode::ok,
       "binary-exported shared-physical Gmsh 4.1 files should be re-readable by import_msh"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::Domain msh41_shared_physical_binary_roundtrip_domain;
  if(!expect(
       sqmesh::mesh::domain_snapshot(
         msh41_shared_physical_binary_roundtrip,
         msh41_shared_physical_binary_roundtrip_domain,
         context
       ) == sqmesh::base::StatusCode::ok,
       "domain_snapshot should resolve the binary round-tripped shared-physical Gmsh 4.1 mesh"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_msh41_shared_physical_entity_domain(
       msh41_shared_physical_binary_roundtrip_domain
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle msh_conflicting_names_mesh = sqmesh::invalid_handle;
  if(!expect(
       store_manual_mesh(
         make_msh_physical_name_conflict_domain(),
         "manual msh physical-name conflict",
         msh_conflicting_names_mesh,
         context
       ),
       "mesh_io_smoke should be able to publish a manual mesh with conflicting MSH physical names"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::mesh::export_msh(
         msh_conflicting_names_mesh,
         msh_conflicting_names_output.string(),
         msh41_export_options,
         context
       ) == sqmesh::base::StatusCode::unsupported,
       "export_msh should fail truthfully when two entity_groups map one physical tag to different names"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle msh41_binary_mesh = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_msh(msh41_binary_input.string(), msh41_binary_mesh, {}, context) ==
         sqmesh::base::StatusCode::ok,
       "import_msh should load the supported binary Gmsh 4.1 subset"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::Domain msh41_binary_domain;
  if(!expect(
       sqmesh::mesh::domain_snapshot(msh41_binary_mesh, msh41_binary_domain, context) ==
         sqmesh::base::StatusCode::ok,
       "domain_snapshot should materialize imported binary Gmsh 4.1 data"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_msh_domain(msh41_binary_domain)) {
    return EXIT_FAILURE;
  }

  if(!expect(
       sqmesh::mesh::export_msh(
         msh41_binary_mesh,
         msh41_binary_output.string(),
         msh41_binary_export_options,
         context
       ) == sqmesh::base::StatusCode::ok,
       "export_msh should write the supported binary Gmsh 4.1 subset"
     )) {
    return EXIT_FAILURE;
  }

  std::string msh41_binary_output_bytes;
  if(!expect(
       read_binary_file(msh41_binary_output, msh41_binary_output_bytes),
       "mesh_io_smoke should be able to read the exported binary Gmsh 4.1 file"
     )) {
    return EXIT_FAILURE;
  }
  const auto msh41_binary_header =
    std::string("$MeshFormat\n4.1 1 ") + std::to_string(sizeof(std::size_t)) + "\n";
  if(!expect(
       msh41_binary_output_bytes.find(msh41_binary_header) != std::string::npos &&
         msh41_binary_output_bytes.find("$Entities\n") != std::string::npos &&
         msh41_binary_output_bytes.find("\"feature_edges\"") != std::string::npos,
       "exported binary Gmsh 4.1 output should contain binary MeshFormat, Entities, and PhysicalNames data"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle msh41_binary_roundtrip = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_msh(
         msh41_binary_output.string(),
         msh41_binary_roundtrip,
         {},
         context
       ) == sqmesh::base::StatusCode::ok,
       "exported binary Gmsh 4.1 files should be re-readable by import_msh"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::Domain msh41_binary_roundtrip_domain;
  if(!expect(
       sqmesh::mesh::domain_snapshot(
         msh41_binary_roundtrip,
         msh41_binary_roundtrip_domain,
         context
       ) == sqmesh::base::StatusCode::ok,
       "domain_snapshot should resolve the round-tripped binary Gmsh 4.1 mesh"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_msh_domain(msh41_binary_roundtrip_domain)) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle msh41_binary_32bit_mesh = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_msh(
         msh41_binary_32bit_input.string(),
         msh41_binary_32bit_mesh,
         {},
         context
       ) == sqmesh::base::StatusCode::ok,
       "import_msh should load binary Gmsh 4.1 fixtures with 4-byte MeshFormat data sizes"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::Domain msh41_binary_32bit_domain;
  if(!expect(
       sqmesh::mesh::domain_snapshot(
         msh41_binary_32bit_mesh,
         msh41_binary_32bit_domain,
         context
       ) == sqmesh::base::StatusCode::ok,
       "domain_snapshot should materialize 4-byte-tag binary Gmsh 4.1 data"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_msh_domain(msh41_binary_32bit_domain)) {
    return EXIT_FAILURE;
  }

  if(!expect(
       sqmesh::mesh::export_msh(
         msh41_binary_32bit_mesh,
         msh41_binary_32bit_output.string(),
         msh41_binary_export_options,
         context
       ) == sqmesh::base::StatusCode::ok,
       "export_msh should rewrite 4-byte-tag binary Gmsh 4.1 input through the supported binary export path"
     )) {
    return EXIT_FAILURE;
  }

  std::string msh41_binary_32bit_output_bytes;
  if(!expect(
       read_binary_file(msh41_binary_32bit_output, msh41_binary_32bit_output_bytes),
       "mesh_io_smoke should be able to read the rewritten 4-byte-tag binary Gmsh 4.1 file"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       msh41_binary_32bit_output_bytes.find(msh41_binary_header) != std::string::npos &&
         msh41_binary_32bit_output_bytes.find("$Entities\n") != std::string::npos,
       "rewritten 4-byte-tag binary Gmsh 4.1 output should advertise the supported binary MeshFormat boundary"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle msh41_binary_32bit_roundtrip = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_msh(
         msh41_binary_32bit_output.string(),
         msh41_binary_32bit_roundtrip,
         {},
         context
       ) == sqmesh::base::StatusCode::ok,
       "rewritten 4-byte-tag binary Gmsh 4.1 files should be re-readable by import_msh"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::Domain msh41_binary_32bit_roundtrip_domain;
  if(!expect(
       sqmesh::mesh::domain_snapshot(
         msh41_binary_32bit_roundtrip,
         msh41_binary_32bit_roundtrip_domain,
         context
       ) == sqmesh::base::StatusCode::ok,
       "domain_snapshot should resolve the rewritten 4-byte-tag binary Gmsh 4.1 mesh"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_msh_domain(msh41_binary_32bit_roundtrip_domain)) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle msh41_partitioned_mesh = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_msh(
         msh41_partitioned_input.string(),
         msh41_partitioned_mesh,
         {},
         context
       ) == sqmesh::base::StatusCode::ok,
       "import_msh should load ASCII Gmsh 4.1 files with partitioned entities, optional sections, and parametric volume-node payloads"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::Domain msh41_partitioned_domain;
  if(!expect(
       sqmesh::mesh::domain_snapshot(
         msh41_partitioned_mesh,
         msh41_partitioned_domain,
         context
       ) == sqmesh::base::StatusCode::ok,
       "domain_snapshot should materialize broadened ASCII Gmsh 4.1 partitioned data"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_msh_domain(msh41_partitioned_domain)) {
    return EXIT_FAILURE;
  }

  if(!expect(
       sqmesh::mesh::export_msh(
         msh41_partitioned_mesh,
         msh41_partitioned_output.string(),
         msh41_export_options,
         context
       ) == sqmesh::base::StatusCode::ok,
       "export_msh should rewrite broadened ASCII Gmsh 4.1 imports through the supported core 4.1 export path"
     )) {
    return EXIT_FAILURE;
  }

  std::string msh41_partitioned_output_text;
  if(!expect(
       read_text_file(msh41_partitioned_output, msh41_partitioned_output_text),
       "mesh_io_smoke should be able to read the rewritten broadened ASCII Gmsh 4.1 file"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       msh41_partitioned_output_text.find("$Entities") != std::string::npos &&
         msh41_partitioned_output_text.find("$PartitionedEntities") == std::string::npos &&
         msh41_partitioned_output_text.find("$Periodic") == std::string::npos &&
         msh41_partitioned_output_text.find("$GhostElements") == std::string::npos &&
         msh41_partitioned_output_text.find("$Parametrizations") == std::string::npos,
       "rewritten broadened ASCII Gmsh 4.1 output should stay within the documented core export boundary"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle msh41_partitioned_roundtrip = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_msh(
         msh41_partitioned_output.string(),
         msh41_partitioned_roundtrip,
         {},
         context
       ) == sqmesh::base::StatusCode::ok,
       "rewritten broadened ASCII Gmsh 4.1 files should be re-readable by import_msh"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::Domain msh41_partitioned_roundtrip_domain;
  if(!expect(
       sqmesh::mesh::domain_snapshot(
         msh41_partitioned_roundtrip,
         msh41_partitioned_roundtrip_domain,
         context
       ) == sqmesh::base::StatusCode::ok,
       "domain_snapshot should resolve the rewritten broadened ASCII Gmsh 4.1 mesh"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_msh_domain(msh41_partitioned_roundtrip_domain)) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle msh41_partitioned_binary_mesh = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_msh(
         msh41_partitioned_binary_input.string(),
         msh41_partitioned_binary_mesh,
         {},
         context
       ) == sqmesh::base::StatusCode::ok,
       "import_msh should load binary Gmsh 4.1 files with partitioned entities and optional 4.1 sections"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::Domain msh41_partitioned_binary_domain;
  if(!expect(
       sqmesh::mesh::domain_snapshot(
         msh41_partitioned_binary_mesh,
         msh41_partitioned_binary_domain,
         context
       ) == sqmesh::base::StatusCode::ok,
       "domain_snapshot should materialize broadened binary Gmsh 4.1 partitioned data"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_msh_domain(msh41_partitioned_binary_domain)) {
    return EXIT_FAILURE;
  }

  if(!expect(
       sqmesh::mesh::export_msh(
         msh41_partitioned_binary_mesh,
         msh41_partitioned_binary_output.string(),
         msh41_binary_export_options,
         context
       ) == sqmesh::base::StatusCode::ok,
       "export_msh should rewrite broadened binary Gmsh 4.1 imports through the supported core binary export path"
     )) {
    return EXIT_FAILURE;
  }

  std::string msh41_partitioned_binary_output_bytes;
  if(!expect(
       read_binary_file(
         msh41_partitioned_binary_output,
         msh41_partitioned_binary_output_bytes
       ),
       "mesh_io_smoke should be able to read the rewritten broadened binary Gmsh 4.1 file"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       msh41_partitioned_binary_output_bytes.find(msh41_binary_header) != std::string::npos &&
         msh41_partitioned_binary_output_bytes.find("$Entities\n") != std::string::npos &&
         msh41_partitioned_binary_output_bytes.find("$PartitionedEntities\n") == std::string::npos &&
         msh41_partitioned_binary_output_bytes.find("$Periodic\n") == std::string::npos &&
         msh41_partitioned_binary_output_bytes.find("$GhostElements\n") == std::string::npos &&
         msh41_partitioned_binary_output_bytes.find("$Parametrizations\n") == std::string::npos,
       "rewritten broadened binary Gmsh 4.1 output should stay within the documented core export boundary"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle msh41_partitioned_binary_roundtrip = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_msh(
         msh41_partitioned_binary_output.string(),
         msh41_partitioned_binary_roundtrip,
         {},
         context
       ) == sqmesh::base::StatusCode::ok,
       "rewritten broadened binary Gmsh 4.1 files should be re-readable by import_msh"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::Domain msh41_partitioned_binary_roundtrip_domain;
  if(!expect(
       sqmesh::mesh::domain_snapshot(
         msh41_partitioned_binary_roundtrip,
         msh41_partitioned_binary_roundtrip_domain,
         context
       ) == sqmesh::base::StatusCode::ok,
       "domain_snapshot should resolve the rewritten broadened binary Gmsh 4.1 mesh"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_msh_domain(msh41_partitioned_binary_roundtrip_domain)) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle obj_mesh = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_obj(obj_input.string(), obj_mesh, {}, context) ==
         sqmesh::base::StatusCode::ok,
       "import_obj should load the supported OBJ subset"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshSummary obj_summary;
  sqmesh::mesh::Domain obj_domain;
  if(!expect(
       sqmesh::mesh::mesh_summary(obj_mesh, obj_summary, context) ==
         sqmesh::base::StatusCode::ok,
       "mesh_summary should resolve imported OBJ meshes"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::mesh::domain_snapshot(obj_mesh, obj_domain, context) ==
         sqmesh::base::StatusCode::ok,
       "domain_snapshot should materialize imported OBJ data into Domain/EntityGroup storage"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       obj_summary.node_count == 3U &&
         obj_summary.edge_count == 3U &&
         obj_summary.face_count == 1U &&
         obj_summary.cell_count == 0U,
       "mesh_summary should report the imported OBJ counts"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_obj_domain(obj_domain)) {
    return EXIT_FAILURE;
  }

  if(!expect(
       sqmesh::mesh::export_obj(obj_mesh, obj_output.string(), {}, context) ==
         sqmesh::base::StatusCode::ok,
       "export_obj should write the supported surface OBJ subset"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle obj_roundtrip = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_obj(obj_output.string(), obj_roundtrip, {}, context) ==
         sqmesh::base::StatusCode::ok,
       "exported OBJ files should be re-readable by import_obj"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::Domain obj_roundtrip_domain;
  if(!expect(
       sqmesh::mesh::domain_snapshot(obj_roundtrip, obj_roundtrip_domain, context) ==
         sqmesh::base::StatusCode::ok,
       "domain_snapshot should resolve the round-tripped OBJ mesh"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_obj_domain(obj_roundtrip_domain)) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::testing::ObjImportReview obj_attribute_review;
  if(!expect(
       sqmesh::mesh::testing::inspect_obj_import_review(
         obj_attribute_input.string(),
         obj_attribute_review
       ) == sqmesh::base::StatusCode::ok,
       "OBJ parser review seam should parse bounded attribute-channel records"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_obj_attribute_review(obj_attribute_review)) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle obj_attribute_mesh = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_obj(obj_attribute_input.string(), obj_attribute_mesh, {}, context) ==
         sqmesh::base::StatusCode::ok,
       "import_obj should accept bounded OBJ face-token forms with attribute references"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshSummary obj_attribute_summary;
  sqmesh::mesh::Domain obj_attribute_domain;
  if(!expect(
       sqmesh::mesh::mesh_summary(obj_attribute_mesh, obj_attribute_summary, context) ==
         sqmesh::base::StatusCode::ok,
       "mesh_summary should resolve OBJ meshes imported from attribute-channel fixtures"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::mesh::domain_snapshot(obj_attribute_mesh, obj_attribute_domain, context) ==
         sqmesh::base::StatusCode::ok,
       "domain_snapshot should materialize OBJ meshes imported from attribute-channel fixtures"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       obj_attribute_summary.node_count == 7U &&
         obj_attribute_summary.edge_count == 2U &&
         obj_attribute_summary.face_count == 4U &&
         obj_attribute_summary.cell_count == 0U,
       "mesh_summary should report the bounded OBJ attribute-channel fixture counts"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       obj_attribute_domain.name() == std::string_view("obj_attr_patch") &&
         obj_attribute_domain.entity_group_count() == 3U,
       "OBJ attribute-channel import should keep the object name and surface-only entity_group layout"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       find_entity_group(obj_attribute_domain, sqmesh::mesh::EntityOrder::edge, "feature_edges") != nullptr &&
         find_entity_group(obj_attribute_domain, sqmesh::mesh::EntityOrder::face, "faces") != nullptr,
       "OBJ attribute-channel import should keep grouped edge and face entity_group names on the bounded subset"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::mesh::export_obj(obj_attribute_mesh, obj_attribute_output.string(), {}, context) ==
         sqmesh::base::StatusCode::ok,
       "export_obj should keep OBJ attribute-channel meshes on the truthful bounded surface subset"
     )) {
    return EXIT_FAILURE;
  }

  std::string obj_attribute_output_text;
  if(!expect(
       read_text_file(obj_attribute_output, obj_attribute_output_text),
       "mesh_io_smoke should be able to read the exported OBJ attribute-channel file"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_obj_attribute_export_truthfulness(obj_attribute_output_text)) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle obj_attribute_roundtrip = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_obj(obj_attribute_output.string(), obj_attribute_roundtrip, {}, context) ==
         sqmesh::base::StatusCode::ok,
       "attribute-channel OBJ exports should stay re-readable after dropping parser-only metadata"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshSummary obj_attribute_roundtrip_summary;
  sqmesh::mesh::Domain obj_attribute_roundtrip_domain;
  if(!expect(
       sqmesh::mesh::mesh_summary(
         obj_attribute_roundtrip,
         obj_attribute_roundtrip_summary,
         context
       ) == sqmesh::base::StatusCode::ok,
       "mesh_summary should resolve round-tripped OBJ meshes exported from the attribute-channel fixture"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::mesh::domain_snapshot(
         obj_attribute_roundtrip,
         obj_attribute_roundtrip_domain,
         context
       ) == sqmesh::base::StatusCode::ok,
       "domain_snapshot should materialize round-tripped OBJ attribute-channel meshes"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       obj_attribute_roundtrip_summary.node_count == 7U &&
         obj_attribute_roundtrip_summary.edge_count == 2U &&
         obj_attribute_roundtrip_summary.face_count == 4U &&
         obj_attribute_roundtrip_summary.cell_count == 0U,
       "mesh_summary should report the round-tripped bounded OBJ attribute-channel counts"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       obj_attribute_roundtrip_domain.name() == std::string_view("obj_attr_patch") &&
         obj_attribute_roundtrip_domain.entity_group_count() == 3U,
       "round-tripped OBJ attribute-channel export should keep the object name and surface-only entity_group layout"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       find_entity_group(obj_attribute_roundtrip_domain, sqmesh::mesh::EntityOrder::edge, "feature_edges") != nullptr &&
         find_entity_group(obj_attribute_roundtrip_domain, sqmesh::mesh::EntityOrder::face, "faces") != nullptr,
       "round-tripped OBJ attribute-channel export should keep grouped edge and face entity_group names"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle obj_concave_mesh = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_obj(obj_concave_input.string(), obj_concave_mesh, {}, context) ==
         sqmesh::base::StatusCode::ok,
       "import_obj should load bounded simple planar concave OBJ polygons through ear clipping"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshSummary obj_concave_summary;
  sqmesh::mesh::Domain obj_concave_domain;
  if(!expect(
       sqmesh::mesh::mesh_summary(obj_concave_mesh, obj_concave_summary, context) ==
         sqmesh::base::StatusCode::ok,
       "mesh_summary should resolve concave OBJ meshes"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::mesh::domain_snapshot(obj_concave_mesh, obj_concave_domain, context) ==
         sqmesh::base::StatusCode::ok,
       "domain_snapshot should materialize concave OBJ data into Domain/EntityGroup storage"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       obj_concave_summary.node_count == 6U &&
         obj_concave_summary.edge_count == 0U &&
         obj_concave_summary.face_count == 4U &&
         obj_concave_summary.cell_count == 0U,
       "mesh_summary should report the concave OBJ counts"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_concave_obj_domain(obj_concave_domain)) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::mesh::export_obj(obj_concave_mesh, obj_concave_output.string(), {}, context) ==
         sqmesh::base::StatusCode::ok,
       "export_obj should emit the ear-clipped concave OBJ mesh as a bounded triangle-only surface"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle obj_concave_roundtrip = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_obj(obj_concave_output.string(), obj_concave_roundtrip, {}, context) ==
         sqmesh::base::StatusCode::ok,
       "ear-clipped concave OBJ exports should stay re-readable after triangle-only export"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshSummary obj_concave_roundtrip_summary;
  sqmesh::mesh::Domain obj_concave_roundtrip_domain;
  if(!expect(
       sqmesh::mesh::mesh_summary(
         obj_concave_roundtrip,
         obj_concave_roundtrip_summary,
         context
       ) == sqmesh::base::StatusCode::ok,
       "mesh_summary should resolve round-tripped concave OBJ meshes"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::mesh::domain_snapshot(
         obj_concave_roundtrip,
         obj_concave_roundtrip_domain,
         context
       ) == sqmesh::base::StatusCode::ok,
       "domain_snapshot should materialize round-tripped concave OBJ meshes"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       obj_concave_roundtrip_summary.node_count == 6U &&
         obj_concave_roundtrip_summary.edge_count == 0U &&
         obj_concave_roundtrip_summary.face_count == 4U &&
         obj_concave_roundtrip_summary.cell_count == 0U,
       "mesh_summary should report the round-tripped concave OBJ counts"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_concave_obj_domain(obj_concave_roundtrip_domain)) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle obj_grouped_mesh = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_obj(obj_grouped_input.string(), obj_grouped_mesh, {}, context) ==
         sqmesh::base::StatusCode::ok,
       "import_obj should load grouped OBJ surfaces with polyline edges and polygon faces"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshSummary obj_grouped_summary;
  sqmesh::mesh::Domain obj_grouped_domain;
  if(!expect(
       sqmesh::mesh::mesh_summary(obj_grouped_mesh, obj_grouped_summary, context) ==
         sqmesh::base::StatusCode::ok,
       "mesh_summary should resolve grouped OBJ meshes"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::mesh::domain_snapshot(obj_grouped_mesh, obj_grouped_domain, context) ==
         sqmesh::base::StatusCode::ok,
       "domain_snapshot should materialize grouped OBJ data into Domain/EntityGroup storage"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       obj_grouped_summary.node_count == 5U &&
         obj_grouped_summary.edge_count == 4U &&
         obj_grouped_summary.face_count == 4U &&
         obj_grouped_summary.cell_count == 0U,
       "mesh_summary should report the grouped OBJ counts"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_grouped_obj_domain(obj_grouped_domain)) {
    return EXIT_FAILURE;
  }

  if(!expect(
       sqmesh::mesh::export_obj(obj_grouped_mesh, obj_grouped_output.string(), {}, context) ==
         sqmesh::base::StatusCode::ok,
       "export_obj should preserve grouped OBJ entity_group labels on the supported broadened subset"
     )) {
    return EXIT_FAILURE;
  }

  std::string obj_grouped_output_text;
  if(!expect(
       read_text_file(obj_grouped_output, obj_grouped_output_text),
       "mesh_io_smoke should be able to read the exported grouped OBJ file"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       obj_grouped_output_text.find("g panel_A seam") != std::string::npos &&
         obj_grouped_output_text.find("g cap") != std::string::npos &&
         obj_grouped_output_text.find("g loose_object") != std::string::npos,
       "exported grouped OBJ output should retain the supported group labels"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle obj_grouped_roundtrip = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_obj(obj_grouped_output.string(), obj_grouped_roundtrip, {}, context) ==
         sqmesh::base::StatusCode::ok,
       "exported grouped OBJ files should be re-readable by import_obj"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::Domain obj_grouped_roundtrip_domain;
  if(!expect(
       sqmesh::mesh::domain_snapshot(obj_grouped_roundtrip, obj_grouped_roundtrip_domain, context) ==
         sqmesh::base::StatusCode::ok,
       "domain_snapshot should resolve the round-tripped grouped OBJ mesh"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_grouped_obj_domain(obj_grouped_roundtrip_domain)) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle obj_unsupported_polygon_mesh = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_obj(
         obj_unsupported_polygon_input.string(),
         obj_unsupported_polygon_mesh,
         {},
         context
       ) == sqmesh::base::StatusCode::unsupported,
       "import_obj should reject self-intersecting OBJ polygons cleanly"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       std::string(sqmesh::base::last_error_message()).find("self-intersecting") != std::string::npos,
       "self-intersecting OBJ rejection should report a reviewer-readable diagnostic"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle obj_collinear_polygon_mesh = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_obj(
         obj_collinear_polygon_input.string(),
         obj_collinear_polygon_mesh,
         {},
         context
       ) == sqmesh::base::StatusCode::unsupported,
       "import_obj should reject OBJ polygons with three consecutive collinear vertices on the bounded subset"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       std::string(sqmesh::base::last_error_message()).find("non-collinear vertices") != std::string::npos,
       "consecutive-collinear OBJ rejection should report the bounded non-collinear constraint"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle nastran_mesh = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_nastran(nastran_input.string(), nastran_mesh, {}, context) ==
         sqmesh::base::StatusCode::ok,
       "import_nastran should load the supported free-field bulk-data subset"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshSummary nastran_summary;
  sqmesh::mesh::Domain nastran_domain;
  if(!expect(
       sqmesh::mesh::mesh_summary(nastran_mesh, nastran_summary, context) ==
         sqmesh::base::StatusCode::ok,
       "mesh_summary should resolve imported NASTRAN meshes"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::mesh::domain_snapshot(nastran_mesh, nastran_domain, context) ==
         sqmesh::base::StatusCode::ok,
       "domain_snapshot should materialize imported NASTRAN data into Domain/EntityGroup storage"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       nastran_summary.node_count == 4U &&
         nastran_summary.edge_count == 6U &&
         nastran_summary.face_count == 4U &&
         nastran_summary.cell_count == 1U,
       "mesh_summary should report the imported NASTRAN counts"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_nastran_domain(nastran_domain)) {
    return EXIT_FAILURE;
  }

  if(!expect(
       sqmesh::mesh::export_nastran(nastran_mesh, nastran_output.string(), {}, context) ==
         sqmesh::base::StatusCode::ok,
       "export_nastran should write the supported free-field bulk-data subset"
     )) {
    return EXIT_FAILURE;
  }

  std::string nastran_output_text;
  if(!expect(
       read_text_file(nastran_output, nastran_output_text),
       "mesh_io_smoke should be able to read the exported bounded NASTRAN file"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       nastran_output_text.find("MAT1,31,210000") != std::string::npos &&
         nastran_output_text.find("PROD,11,31,1") != std::string::npos &&
         nastran_output_text.find("PSHELL,21,31,1") != std::string::npos &&
         nastran_output_text.find("PSOLID,31,31") != std::string::npos,
       "exported bounded NASTRAN output should preserve MAT1, PROD, PSHELL, and PSOLID metadata"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle nastran_roundtrip = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_nastran(nastran_output.string(), nastran_roundtrip, {}, context) ==
         sqmesh::base::StatusCode::ok,
       "exported NASTRAN files should be re-readable by import_nastran"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::Domain nastran_roundtrip_domain;
  if(!expect(
       sqmesh::mesh::domain_snapshot(nastran_roundtrip, nastran_roundtrip_domain, context) ==
         sqmesh::base::StatusCode::ok,
       "domain_snapshot should resolve the round-tripped NASTRAN mesh"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_nastran_domain(nastran_roundtrip_domain)) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle nastran_small_field_mesh = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_nastran(
         nastran_fixed_input.string(),
         nastran_small_field_mesh,
         {},
         context
       ) == sqmesh::base::StatusCode::ok,
       "import_nastran should also accept the supported 8-character small-field subset"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::Domain nastran_small_field_domain;
  if(!expect(
       sqmesh::mesh::domain_snapshot(
         nastran_small_field_mesh,
         nastran_small_field_domain,
         context
       ) == sqmesh::base::StatusCode::ok,
       "domain_snapshot should resolve supported small-field NASTRAN imports"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_nastran_domain(nastran_small_field_domain)) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle nastran_broad_mesh = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_nastran(nastran_broad_input.string(), nastran_broad_mesh, {}, context) ==
         sqmesh::base::StatusCode::ok,
       "import_nastran should load the broadened free-field NASTRAN subset with beam, shell, and property cards"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::Domain nastran_broad_domain;
  if(!expect(
       sqmesh::mesh::domain_snapshot(nastran_broad_mesh, nastran_broad_domain, context) ==
         sqmesh::base::StatusCode::ok,
       "domain_snapshot should resolve the broadened free-field NASTRAN fixture"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_broadened_nastran_surface_domain(nastran_broad_domain)) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle nastran_broad_small_field_mesh = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_nastran(
         nastran_broad_fixed_input.string(),
         nastran_broad_small_field_mesh,
         {},
         context
       ) == sqmesh::base::StatusCode::ok,
       "import_nastran should load the broadened 8-character small-field NASTRAN subset"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::Domain nastran_broad_small_field_domain;
  if(!expect(
       sqmesh::mesh::domain_snapshot(
         nastran_broad_small_field_mesh,
         nastran_broad_small_field_domain,
         context
       ) == sqmesh::base::StatusCode::ok,
       "domain_snapshot should resolve the broadened small-field NASTRAN fixture"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_broadened_nastran_surface_domain(nastran_broad_small_field_domain)) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle nastran_broad_large_field_mesh = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_nastran(
         nastran_broad_large_input.string(),
         nastran_broad_large_field_mesh,
         {},
         context
       ) == sqmesh::base::StatusCode::ok,
       "import_nastran should load the broadened 16-character large-field NASTRAN subset"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::Domain nastran_broad_large_field_domain;
  if(!expect(
       sqmesh::mesh::domain_snapshot(
         nastran_broad_large_field_mesh,
         nastran_broad_large_field_domain,
         context
       ) == sqmesh::base::StatusCode::ok,
       "domain_snapshot should resolve the broadened large-field NASTRAN fixture"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_broadened_nastran_surface_domain(nastran_broad_large_field_domain)) {
    return EXIT_FAILURE;
  }

  if(!expect(
       sqmesh::mesh::export_nastran(nastran_broad_mesh, nastran_broad_output.string(), {}, context) ==
         sqmesh::base::StatusCode::ok,
       "export_nastran should preserve the broadened supported shell-family subset as free-field bulk data"
     )) {
    return EXIT_FAILURE;
  }

  std::string nastran_broad_output_text;
  if(!expect(
       read_text_file(nastran_broad_output, nastran_broad_output_text),
       "mesh_io_smoke should be able to read the exported free-field NASTRAN file"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       nastran_broad_output_text.find("CQUAD4,") != std::string::npos &&
         nastran_broad_output_text.find("CTRIA3,") != std::string::npos &&
         nastran_broad_output_text.find("MAT1,301,210000") != std::string::npos &&
         nastran_broad_output_text.find("PROD,101,301,1") != std::string::npos &&
         nastran_broad_output_text.find("PSHELL,201,301,1") != std::string::npos,
       "exported free-field NASTRAN output should preserve bounded property/material metadata alongside both shell families when the imported quad pairing is intact"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle nastran_broad_roundtrip = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_nastran(nastran_broad_output.string(), nastran_broad_roundtrip, {}, context) ==
         sqmesh::base::StatusCode::ok,
       "exported broadened free-field NASTRAN files should be re-readable by import_nastran"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::Domain nastran_broad_roundtrip_domain;
  if(!expect(
       sqmesh::mesh::domain_snapshot(
         nastran_broad_roundtrip,
         nastran_broad_roundtrip_domain,
         context
       ) == sqmesh::base::StatusCode::ok,
       "domain_snapshot should resolve the broadened free-field NASTRAN round-trip mesh"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_broadened_nastran_surface_domain(
       nastran_broad_roundtrip_domain,
       sqmesh::mesh::NastranEntityGroupSourceCard::crod
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::NastranExportOptions nastran_small_export_options {};
  nastran_small_export_options.field_format = sqmesh::mesh::NastranFieldFormat::small_fixed;
  if(!expect(
       sqmesh::mesh::export_nastran(
         nastran_broad_mesh,
         nastran_broad_fixed_output.string(),
         nastran_small_export_options,
         context
       ) == sqmesh::base::StatusCode::ok,
       "export_nastran should also write the supported small-field NASTRAN subset"
     )) {
    return EXIT_FAILURE;
  }

  std::string nastran_broad_fixed_output_text;
  if(!expect(
       read_text_file(nastran_broad_fixed_output, nastran_broad_fixed_output_text),
       "mesh_io_smoke should be able to read the exported small-field NASTRAN file"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       nastran_broad_fixed_output_text.find("GRID    ") != std::string::npos &&
         nastran_broad_fixed_output_text.find("CQUAD4") != std::string::npos &&
         nastran_broad_fixed_output_text.find("MAT1") != std::string::npos &&
         nastran_broad_fixed_output_text.find("PROD") != std::string::npos &&
         nastran_broad_fixed_output_text.find("PSHELL") != std::string::npos &&
         nastran_broad_fixed_output_text.find("GRID,") == std::string::npos &&
         nastran_broad_fixed_output_text.find("CQUAD4,") == std::string::npos,
       "exported small-field NASTRAN output should use fixed-width cards and preserve bounded property/material metadata alongside CQUAD4"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle nastran_broad_fixed_roundtrip = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_nastran(
         nastran_broad_fixed_output.string(),
         nastran_broad_fixed_roundtrip,
         {},
         context
       ) == sqmesh::base::StatusCode::ok,
       "exported broadened small-field NASTRAN files should be re-readable by import_nastran"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::Domain nastran_broad_fixed_roundtrip_domain;
  if(!expect(
       sqmesh::mesh::domain_snapshot(
         nastran_broad_fixed_roundtrip,
         nastran_broad_fixed_roundtrip_domain,
         context
       ) == sqmesh::base::StatusCode::ok,
       "domain_snapshot should resolve the broadened small-field NASTRAN round-trip mesh"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_broadened_nastran_surface_domain(
       nastran_broad_fixed_roundtrip_domain,
       sqmesh::mesh::NastranEntityGroupSourceCard::crod
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::NastranExportOptions nastran_large_export_options {};
  nastran_large_export_options.field_format = sqmesh::mesh::NastranFieldFormat::large_fixed;
  if(!expect(
       sqmesh::mesh::export_nastran(
         nastran_broad_mesh,
         nastran_broad_large_output.string(),
         nastran_large_export_options,
         context
       ) == sqmesh::base::StatusCode::ok,
       "export_nastran should also write the supported large-field NASTRAN subset"
     )) {
    return EXIT_FAILURE;
  }

  std::string nastran_broad_large_output_text;
  if(!expect(
       read_text_file(nastran_broad_large_output, nastran_broad_large_output_text),
       "mesh_io_smoke should be able to read the exported large-field NASTRAN file"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       nastran_broad_large_output_text.find("GRID*") != std::string::npos &&
         nastran_broad_large_output_text.find("CQUAD4*") != std::string::npos &&
         nastran_broad_large_output_text.find("MAT1*") != std::string::npos &&
         nastran_broad_large_output_text.find("PROD*") != std::string::npos &&
         nastran_broad_large_output_text.find("PSHELL*") != std::string::npos &&
         nastran_broad_large_output_text.find("*       ") != std::string::npos &&
         nastran_broad_large_output_text.find("GRID,") == std::string::npos &&
         nastran_broad_large_output_text.find("CQUAD4,") == std::string::npos,
       "exported large-field NASTRAN output should use canonical star continuations and preserve bounded property/material metadata alongside CQUAD4"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle nastran_broad_large_roundtrip = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_nastran(
         nastran_broad_large_output.string(),
         nastran_broad_large_roundtrip,
         {},
         context
       ) == sqmesh::base::StatusCode::ok,
       "exported broadened large-field NASTRAN files should be re-readable by import_nastran"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::Domain nastran_broad_large_roundtrip_domain;
  if(!expect(
       sqmesh::mesh::domain_snapshot(
         nastran_broad_large_roundtrip,
         nastran_broad_large_roundtrip_domain,
         context
       ) == sqmesh::base::StatusCode::ok,
       "domain_snapshot should resolve the broadened large-field NASTRAN round-trip mesh"
     )) {
    return EXIT_FAILURE;
  }
  if(!verify_broadened_nastran_surface_domain(
       nastran_broad_large_roundtrip_domain,
       sqmesh::mesh::NastranEntityGroupSourceCard::crod
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle nastran_unpaired_quad_mesh = sqmesh::invalid_handle;
  if(!expect(
       store_manual_mesh(
         make_nastran_unpaired_cquad4_domain(),
         "manual nastran unpaired cquad4",
         nastran_unpaired_quad_mesh,
         context
       ),
       "mesh_io_smoke should be able to publish a manual CQUAD4-entity_group mesh without per-face pairing metadata"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::mesh::export_nastran(
         nastran_unpaired_quad_mesh,
         nastran_unpaired_quad_output.string(),
         {},
         context
       ) == sqmesh::base::StatusCode::ok,
       "export_nastran should still write triangle faces when a CQUAD4 entity_group lacks trustworthy pairing metadata"
     )) {
    return EXIT_FAILURE;
  }

  std::string nastran_unpaired_quad_output_text;
  if(!expect(
       read_text_file(nastran_unpaired_quad_output, nastran_unpaired_quad_output_text),
       "mesh_io_smoke should be able to read the exported unpaired CQUAD4 regression file"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       nastran_unpaired_quad_output_text.find("CQUAD4,") == std::string::npos &&
         nastran_unpaired_quad_output_text.find("CTRIA3,") != std::string::npos,
       "export_nastran should fall back to CTRIA3 when a CQUAD4 entity_group has no matching per-face source tags"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle nastran_unsupported_mesh = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::import_nastran(
         nastran_unsupported.string(),
         nastran_unsupported_mesh,
         {},
         context
       ) == sqmesh::base::StatusCode::unsupported,
       "import_nastran should reject unsupported bulk cards cleanly"
     )) {
    return EXIT_FAILURE;
  }

  if(!sqmesh::mesh::cgns_io_available()) {
    sqmesh::mesh::MeshHandle cgns_mesh = sqmesh::invalid_handle;
    if(!expect(
         sqmesh::mesh::import_cgns("missing.cgns", cgns_mesh, {}, context) ==
           sqmesh::base::StatusCode::unsupported,
         "import_cgns should report unsupported when CGNS IO is disabled"
       )) {
      return EXIT_FAILURE;
    }
    if(!expect(
         sqmesh::mesh::export_cgns(msh_mesh, cgns_output.string(), {}, context) ==
           sqmesh::base::StatusCode::unsupported,
         "export_cgns should report unsupported when CGNS IO is disabled"
       )) {
      return EXIT_FAILURE;
    }
  }
  else {
#ifdef sqmesh_TEST_CGNS_ENABLED
    const std::vector<CgnsSectionExpectation> expected_tetra_sections = {
      {"feature_edges", CGNS_ENUMV(BAR_2), "11", 1, 6},
      {"wall", CGNS_ENUMV(TRI_3), "21", 7, 10},
      {"volume", CGNS_ENUMV(TETRA_4), "31", 11, 11},
    };
    const std::vector<CgnsBoundaryExpectation> expected_tetra_boundaries = {
      {"feature_edges", CGNS_ENUMV(BCSymmetryPlane), CGNS_ENUMV(PointRange), CGNS_ENUMV(EdgeCenter), {1, 6}},
      {"wall", CGNS_ENUMV(BCWall), CGNS_ENUMV(PointRange), CGNS_ENUMV(FaceCenter), {7, 10}},
    };
    const std::vector<CgnsSectionExpectation> expected_mixed_sections = {
      {"outer_loop", CGNS_ENUMV(BAR_2), "41", 1, 5},
      {"surface_patch", CGNS_ENUMV(TRI_3), "41", 6, 8},
    };
    const std::vector<CgnsBoundaryExpectation> expected_mixed_boundaries = {
      {"outer_loop", CGNS_ENUMV(BCSymmetryPlane), CGNS_ENUMV(PointRange), CGNS_ENUMV(EdgeCenter), {1, 5}},
      {"surface_patch", CGNS_ENUMV(BCInflow), CGNS_ENUMV(PointRange), CGNS_ENUMV(CellCenter), {6, 8}},
    };
    const std::vector<CgnsSectionExpectation> expected_multizone_sections = {
      {"z1_zone_a__feature_edges", CGNS_ENUMV(BAR_2), "11", 1, 6},
      {"z2_zone_b__feature_edges", CGNS_ENUMV(BAR_2), "11", 7, 12},
      {"z1_zone_a__wall", CGNS_ENUMV(TRI_3), "21", 13, 16},
      {"z2_zone_b__wall", CGNS_ENUMV(TRI_3), "21", 17, 20},
      {"z1_zone_a__volume", CGNS_ENUMV(TETRA_4), "31", 21, 21},
      {"z2_zone_b__volume", CGNS_ENUMV(TETRA_4), "31", 22, 22},
    };
    const std::vector<CgnsBoundaryExpectation> expected_multizone_boundaries = {
      {"z1_zone_a__feature_edges", CGNS_ENUMV(BCSymmetryPlane), CGNS_ENUMV(PointRange), CGNS_ENUMV(EdgeCenter), {1, 6}},
      {"z2_zone_b__feature_edges", CGNS_ENUMV(BCOutflow), CGNS_ENUMV(PointRange), CGNS_ENUMV(EdgeCenter), {7, 12}},
      {"z1_zone_a__wall", CGNS_ENUMV(BCWall), CGNS_ENUMV(PointRange), CGNS_ENUMV(FaceCenter), {13, 16}},
      {"z2_zone_b__wall", CGNS_ENUMV(BCWallInviscid), CGNS_ENUMV(PointRange), CGNS_ENUMV(FaceCenter), {17, 20}},
    };
    const std::vector<CgnsSectionExpectation> expected_multibase_sections = {
      {"b1_sqmesh_base_a__z1_shared_zone", CGNS_ENUMV(BAR_2), "11", 1, 6},
      {"b2_sqmesh_base_b__z1_shared_zone", CGNS_ENUMV(BAR_2), "11", 7, 12},
      {"b1_sqmesh_base_a__z1_shared_zo_1", CGNS_ENUMV(TRI_3), "21", 13, 16},
      {"b2_sqmesh_base_b__z1_shared_zo_1", CGNS_ENUMV(TRI_3), "21", 17, 20},
      {"b1_sqmesh_base_a__z1_shared_zo_2", CGNS_ENUMV(TETRA_4), "31", 21, 21},
      {"b2_sqmesh_base_b__z1_shared_zo_2", CGNS_ENUMV(TETRA_4), "31", 22, 22},
    };
    const std::vector<CgnsBoundaryExpectation> expected_multibase_boundaries = {
      {"b1_sqmesh_base_a__z1_shared_zone", CGNS_ENUMV(BCSymmetryPlane), CGNS_ENUMV(PointRange), CGNS_ENUMV(EdgeCenter), {1, 6}},
      {"b2_sqmesh_base_b__z1_shared_zone", CGNS_ENUMV(BCOutflow), CGNS_ENUMV(PointRange), CGNS_ENUMV(EdgeCenter), {7, 12}},
      {"b1_sqmesh_base_a__z1_shared_zo_1", CGNS_ENUMV(BCWall), CGNS_ENUMV(PointRange), CGNS_ENUMV(FaceCenter), {13, 16}},
      {"b2_sqmesh_base_b__z1_shared_zo_1", CGNS_ENUMV(BCWallInviscid), CGNS_ENUMV(PointRange), CGNS_ENUMV(FaceCenter), {17, 20}},
    };
    const std::vector<CgnsZoneExpectation> expected_preserved_multizone_zones = {
      {
        "zone_a",
        {
          {"feature_edges", CGNS_ENUMV(BAR_2), "11", 1, 6},
          {"wall", CGNS_ENUMV(TRI_3), "21", 7, 10},
          {"volume", CGNS_ENUMV(TETRA_4), "31", 11, 11},
        },
        {
          {"feature_edges", CGNS_ENUMV(BCSymmetryPlane), CGNS_ENUMV(PointRange), CGNS_ENUMV(EdgeCenter), {1, 6}},
          {"wall", CGNS_ENUMV(BCWall), CGNS_ENUMV(PointRange), CGNS_ENUMV(FaceCenter), {7, 10}},
        },
      },
      {
        "zone_b",
        {
          {"feature_edges", CGNS_ENUMV(BAR_2), "11", 1, 6},
          {"wall", CGNS_ENUMV(TRI_3), "21", 7, 10},
          {"volume", CGNS_ENUMV(TETRA_4), "31", 11, 11},
        },
        {
          {"feature_edges", CGNS_ENUMV(BCOutflow), CGNS_ENUMV(PointRange), CGNS_ENUMV(EdgeCenter), {1, 6}},
          {"wall", CGNS_ENUMV(BCWallInviscid), CGNS_ENUMV(PointRange), CGNS_ENUMV(FaceCenter), {7, 10}},
        },
      },
    };
    const std::vector<CgnsZoneExpectation> expected_surface_interface_multizone_zones = {
      {
        "zone_a",
        {
          {"outer_loop", CGNS_ENUMV(BAR_2), "11", 1, 2},
          {"surface_patch", CGNS_ENUMV(TRI_3), "21", 3, 3},
        },
        {
          {"outer_loop", CGNS_ENUMV(BCSymmetryPlane), CGNS_ENUMV(PointRange), CGNS_ENUMV(EdgeCenter), {1, 2}},
          {"surface_patch", CGNS_ENUMV(BCWall), CGNS_ENUMV(PointRange), CGNS_ENUMV(CellCenter), {3, 3}},
        },
        {
          {
            "zone_interface",
            CGNS_ENUMV(Vertex),
            CGNS_ENUMV(Abutting1to1),
            CGNS_ENUMV(PointList),
            "zone_b",
            CGNS_ENUMV(Unstructured),
            CGNS_ENUMV(PointListDonor),
            CGNS_ENUMV(DataTypeNull),
            {1, 2},
            {1, 2},
          },
        },
      },
      {
        "zone_b",
        {
          {"outer_loop", CGNS_ENUMV(BAR_2), "11", 1, 2},
          {"surface_patch", CGNS_ENUMV(TRI_3), "21", 3, 3},
        },
        {
          {"outer_loop", CGNS_ENUMV(BCOutflow), CGNS_ENUMV(PointRange), CGNS_ENUMV(EdgeCenter), {1, 2}},
          {"surface_patch", CGNS_ENUMV(BCWallInviscid), CGNS_ENUMV(PointRange), CGNS_ENUMV(CellCenter), {3, 3}},
        },
        {
          {
            "zone_interface",
            CGNS_ENUMV(Vertex),
            CGNS_ENUMV(Abutting1to1),
            CGNS_ENUMV(PointList),
            "zone_a",
            CGNS_ENUMV(Unstructured),
            CGNS_ENUMV(PointListDonor),
            CGNS_ENUMV(DataTypeNull),
            {1, 2},
            {1, 2},
          },
        },
      },
    };
    const std::vector<CgnsZoneExpectation> expected_surface_polyline_interface_multizone_zones = {
      {
        "zone_a",
        {
          {"outer_loop", CGNS_ENUMV(BAR_2), "11", 1, 4},
          {"surface_patch", CGNS_ENUMV(TRI_3), "21", 5, 6},
        },
        {
          {"outer_loop", CGNS_ENUMV(BCSymmetryPlane), CGNS_ENUMV(PointRange), CGNS_ENUMV(EdgeCenter), {1, 4}},
          {"surface_patch", CGNS_ENUMV(BCWall), CGNS_ENUMV(PointRange), CGNS_ENUMV(CellCenter), {5, 6}},
        },
        {
          {
            "zone_interface",
            CGNS_ENUMV(Vertex),
            CGNS_ENUMV(Abutting1to1),
            CGNS_ENUMV(PointList),
            "zone_b",
            CGNS_ENUMV(Unstructured),
            CGNS_ENUMV(PointListDonor),
            CGNS_ENUMV(DataTypeNull),
            {1, 2, 4},
            {1, 3, 5},
          },
        },
      },
      {
        "zone_b",
        {
          {"outer_loop", CGNS_ENUMV(BAR_2), "11", 1, 4},
          {"surface_patch", CGNS_ENUMV(TRI_3), "21", 5, 6},
        },
        {
          {"outer_loop", CGNS_ENUMV(BCOutflow), CGNS_ENUMV(PointRange), CGNS_ENUMV(EdgeCenter), {1, 4}},
          {"surface_patch", CGNS_ENUMV(BCWallInviscid), CGNS_ENUMV(PointRange), CGNS_ENUMV(CellCenter), {5, 6}},
        },
        {
          {
            "zone_interface",
            CGNS_ENUMV(Vertex),
            CGNS_ENUMV(Abutting1to1),
            CGNS_ENUMV(PointList),
            "zone_a",
            CGNS_ENUMV(Unstructured),
            CGNS_ENUMV(PointListDonor),
            CGNS_ENUMV(DataTypeNull),
            {1, 3, 5},
            {1, 2, 4},
          },
        },
      },
    };
    const std::vector<CgnsZoneExpectation> expected_interface_multizone_zones = {
      {
        "zone_a",
        {
          {"wall", CGNS_ENUMV(TRI_3), "21", 1, 3},
          {"volume", CGNS_ENUMV(TETRA_4), "31", 4, 4},
        },
        {
          {"wall", CGNS_ENUMV(BCWall), CGNS_ENUMV(PointRange), CGNS_ENUMV(FaceCenter), {1, 3}},
        },
        {
          {
            "zone_interface",
            CGNS_ENUMV(Vertex),
            CGNS_ENUMV(Abutting1to1),
            CGNS_ENUMV(PointList),
            "zone_b",
            CGNS_ENUMV(Unstructured),
            CGNS_ENUMV(PointListDonor),
            CGNS_ENUMV(DataTypeNull),
            {1, 2, 3},
            {1, 2, 3},
          },
        },
      },
      {
        "zone_b",
        {
          {"wall", CGNS_ENUMV(TRI_3), "21", 1, 3},
          {"volume", CGNS_ENUMV(TETRA_4), "31", 4, 4},
        },
        {
          {"wall", CGNS_ENUMV(BCWallInviscid), CGNS_ENUMV(PointRange), CGNS_ENUMV(FaceCenter), {1, 3}},
        },
        {
          {
            "zone_interface",
            CGNS_ENUMV(Vertex),
            CGNS_ENUMV(Abutting1to1),
            CGNS_ENUMV(PointList),
            "zone_a",
            CGNS_ENUMV(Unstructured),
            CGNS_ENUMV(PointListDonor),
            CGNS_ENUMV(DataTypeNull),
            {1, 2, 3},
            {1, 2, 3},
          },
        },
      },
    };
    const std::vector<CgnsZoneExpectation> expected_patch_interface_multizone_zones = {
      {
        "zone_a",
        {
          {"wall", CGNS_ENUMV(TRI_3), "21", 1, 6},
          {"volume", CGNS_ENUMV(TETRA_4), "31", 7, 8},
        },
        {
          {"wall", CGNS_ENUMV(BCWall), CGNS_ENUMV(PointRange), CGNS_ENUMV(FaceCenter), {1, 6}},
        },
        {
          {
            "zone_interface",
            CGNS_ENUMV(Vertex),
            CGNS_ENUMV(Abutting1to1),
            CGNS_ENUMV(PointList),
            "zone_b",
            CGNS_ENUMV(Unstructured),
            CGNS_ENUMV(PointListDonor),
            CGNS_ENUMV(DataTypeNull),
            {1, 2, 3, 5},
            {1, 2, 3, 5},
          },
        },
      },
      {
        "zone_b",
        {
          {"wall", CGNS_ENUMV(TRI_3), "21", 1, 6},
          {"volume", CGNS_ENUMV(TETRA_4), "31", 7, 8},
        },
        {
          {"wall", CGNS_ENUMV(BCWallInviscid), CGNS_ENUMV(PointRange), CGNS_ENUMV(FaceCenter), {1, 6}},
        },
        {
          {
            "zone_interface",
            CGNS_ENUMV(Vertex),
            CGNS_ENUMV(Abutting1to1),
            CGNS_ENUMV(PointList),
            "zone_a",
            CGNS_ENUMV(Unstructured),
            CGNS_ENUMV(PointListDonor),
            CGNS_ENUMV(DataTypeNull),
            {1, 2, 3, 5},
            {1, 2, 3, 5},
          },
        },
      },
    };

    if(!expect(
         write_supported_cgns_fixture(cgns_input),
         "the CGNS smoke fixture should be writable with cgnslib"
       ) ||
       !expect(
         write_mixed_surface_cgns_fixture(cgns_mixed_input),
         "the mixed CGNS smoke fixture should be writable with cgnslib"
       ) ||
       !expect(
         write_multi_zone_cgns_fixture(cgns_multizone_input),
         "the multi-zone CGNS smoke fixture should be writable with cgnslib"
       ) ||
       !expect(
         write_multi_base_cgns_fixture(cgns_multibase_input),
         "the multi-base CGNS smoke fixture should be writable with cgnslib"
       ) ||
       !expect(
         write_unsupported_cgns_physical_dim_fixture(cgns_unsupported_physical_dim),
         "the unsupported physical-dimension CGNS fixture should be writable with cgnslib"
       ) ||
       !expect(
         write_unsupported_cgns_fixture(cgns_unsupported),
         "the unsupported CGNS smoke fixture should be writable with cgnslib"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::MeshHandle cgns_mesh = sqmesh::invalid_handle;
    if(!expect(
         sqmesh::mesh::import_cgns(cgns_input.string(), cgns_mesh, {}, context) ==
           sqmesh::base::StatusCode::ok,
         "import_cgns should load the supported single-base single-zone subset"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::Domain cgns_domain;
    if(!expect(
         sqmesh::mesh::domain_snapshot(cgns_mesh, cgns_domain, context) ==
           sqmesh::base::StatusCode::ok,
         "domain_snapshot should materialize imported CGNS data into Domain/EntityGroup storage"
       )) {
      return EXIT_FAILURE;
    }
    if(!verify_cgns_domain(cgns_domain, "")) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::CgnsExportOptions cgns_export_options;
    cgns_export_options.base_name = "sqmesh_base";
    cgns_export_options.zone_name = "sqmesh_cgns_zone";
    if(!expect(
         sqmesh::mesh::export_cgns(cgns_mesh, cgns_output.string(), cgns_export_options, context) ==
           sqmesh::base::StatusCode::ok,
         "export_cgns should write the supported CGNS subset"
       )) {
      return EXIT_FAILURE;
    }
    if(!expect(
         verify_exported_cgns_file(
           cgns_output,
           "sqmesh_cgns_zone",
           expected_tetra_sections,
           expected_tetra_boundaries
         ),
         "export_cgns should produce the supported tetrahedral CGNS subset with boundary conditions"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::MeshHandle cgns_roundtrip = sqmesh::invalid_handle;
    if(!expect(
         sqmesh::mesh::import_cgns(cgns_output.string(), cgns_roundtrip, {}, context) ==
           sqmesh::base::StatusCode::ok,
         "exported CGNS files should be re-readable by import_cgns"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::Domain cgns_roundtrip_domain;
    if(!expect(
         sqmesh::mesh::domain_snapshot(cgns_roundtrip, cgns_roundtrip_domain, context) ==
           sqmesh::base::StatusCode::ok,
         "domain_snapshot should resolve the round-tripped CGNS mesh"
       )) {
      return EXIT_FAILURE;
    }
    if(!verify_cgns_domain(cgns_roundtrip_domain, "")) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::MeshHandle cgns_mixed_mesh = sqmesh::invalid_handle;
    if(!expect(
         sqmesh::mesh::import_cgns(cgns_mixed_input.string(), cgns_mixed_mesh, {}, context) ==
           sqmesh::base::StatusCode::ok,
         "import_cgns should load the broadened MIXED/QUAD_4 surface subset"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::Domain cgns_mixed_domain;
    if(!expect(
         sqmesh::mesh::domain_snapshot(cgns_mixed_mesh, cgns_mixed_domain, context) ==
           sqmesh::base::StatusCode::ok,
         "domain_snapshot should materialize broadened CGNS surface data into Domain/EntityGroup storage"
       )) {
      return EXIT_FAILURE;
    }
    if(!verify_cgns_mixed_surface_domain(cgns_mixed_domain, "")) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::CgnsExportOptions cgns_mixed_export_options;
    cgns_mixed_export_options.base_name = "sqmesh_base";
    cgns_mixed_export_options.zone_name = "sqmesh_cgns_mixed_zone";
    if(!expect(
         sqmesh::mesh::export_cgns(
           cgns_mixed_mesh,
           cgns_mixed_output.string(),
           cgns_mixed_export_options,
           context
         ) == sqmesh::base::StatusCode::ok,
         "export_cgns should write the supported triangulated surface CGNS subset after MIXED/QUAD_4 import"
       )) {
      return EXIT_FAILURE;
    }
    if(!expect(
         verify_exported_cgns_file(
           cgns_mixed_output,
           "sqmesh_cgns_mixed_zone",
           expected_mixed_sections,
           expected_mixed_boundaries
         ),
         "export_cgns should preserve broadened CGNS boundary metadata on supported surface output"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::MeshHandle cgns_mixed_roundtrip = sqmesh::invalid_handle;
    if(!expect(
         sqmesh::mesh::import_cgns(cgns_mixed_output.string(), cgns_mixed_roundtrip, {}, context) ==
           sqmesh::base::StatusCode::ok,
         "exported broadened CGNS surface files should be re-readable by import_cgns"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::Domain cgns_mixed_roundtrip_domain;
    if(!expect(
         sqmesh::mesh::domain_snapshot(cgns_mixed_roundtrip, cgns_mixed_roundtrip_domain, context) ==
           sqmesh::base::StatusCode::ok,
         "domain_snapshot should resolve the round-tripped broadened CGNS surface mesh"
       )) {
      return EXIT_FAILURE;
    }
    if(!verify_cgns_mixed_surface_domain(
         cgns_mixed_roundtrip_domain,
         ""
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::MeshHandle cgns_multizone_mesh = sqmesh::invalid_handle;
    if(!expect(
         sqmesh::mesh::import_cgns(cgns_multizone_input.string(), cgns_multizone_mesh, {}, context) ==
           sqmesh::base::StatusCode::ok,
         "import_cgns should flatten the supported single-base multi-zone subset into zone-qualified entity_groups"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::Domain cgns_multizone_domain;
    if(!expect(
         sqmesh::mesh::domain_snapshot(cgns_multizone_mesh, cgns_multizone_domain, context) ==
           sqmesh::base::StatusCode::ok,
         "domain_snapshot should materialize flattened multi-zone CGNS data into Domain/EntityGroup storage"
       )) {
      return EXIT_FAILURE;
    }
    if(!verify_cgns_multi_zone_domain(
         cgns_multizone_domain,
         "sqmesh_base",
         ""
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::CgnsExportOptions cgns_multizone_export_options;
    cgns_multizone_export_options.base_name = "sqmesh_base";
    cgns_multizone_export_options.zone_name = "sqmesh_cgns_flattened_multi_zone";
    if(!expect(
         sqmesh::mesh::export_cgns(
           cgns_multizone_mesh,
           cgns_multizone_output.string(),
           cgns_multizone_export_options,
           context
         ) == sqmesh::base::StatusCode::ok,
         "export_cgns should rewrite flattened multi-zone input through the documented single-zone export path"
       )) {
      return EXIT_FAILURE;
    }
    if(!expect(
         verify_exported_cgns_file(
           cgns_multizone_output,
           "sqmesh_cgns_flattened_multi_zone",
           expected_multizone_sections,
           expected_multizone_boundaries
         ),
         "export_cgns should preserve the flattened multi-zone section and boundary partitioning inside one written zone"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::MeshHandle cgns_multizone_roundtrip = sqmesh::invalid_handle;
    if(!expect(
         sqmesh::mesh::import_cgns(cgns_multizone_output.string(), cgns_multizone_roundtrip, {}, context) ==
           sqmesh::base::StatusCode::ok,
         "the flattened single-zone CGNS export from a multi-zone import should be re-readable by import_cgns"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::Domain cgns_multizone_roundtrip_domain;
    if(!expect(
         sqmesh::mesh::domain_snapshot(cgns_multizone_roundtrip, cgns_multizone_roundtrip_domain, context) ==
           sqmesh::base::StatusCode::ok,
         "domain_snapshot should resolve the flattened multi-zone CGNS round-trip mesh"
       )) {
      return EXIT_FAILURE;
    }
    if(!verify_cgns_multi_zone_domain(
         cgns_multizone_roundtrip_domain,
         "sqmesh_cgns_flattened_multi_zone",
         ""
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::CgnsExportOptions cgns_preserved_multizone_export_options;
    cgns_preserved_multizone_export_options.base_name = "sqmesh_base";
    cgns_preserved_multizone_export_options.write_multi_zone = true;
    if(!expect(
         sqmesh::mesh::export_cgns(
           cgns_multizone_mesh,
           cgns_multizone_preserved_output.string(),
           cgns_preserved_multizone_export_options,
           context
         ) == sqmesh::base::StatusCode::ok,
         "export_cgns should preserve the documented multi-zone CGNS subset when write_multi_zone is enabled"
       )) {
      return EXIT_FAILURE;
    }
    if(!expect(
         verify_exported_multi_zone_cgns_file(
           cgns_multizone_preserved_output,
           expected_preserved_multizone_zones
         ),
         "export_cgns should preserve multiple supported CGNS zones with local section and boundary names"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::MeshHandle cgns_preserved_multizone_roundtrip = sqmesh::invalid_handle;
    if(!expect(
         sqmesh::mesh::import_cgns(
           cgns_multizone_preserved_output.string(),
           cgns_preserved_multizone_roundtrip,
           {},
           context
         ) == sqmesh::base::StatusCode::ok,
         "the multi-zone CGNS export should be re-readable by import_cgns"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::Domain cgns_preserved_multizone_roundtrip_domain;
    if(!expect(
         sqmesh::mesh::domain_snapshot(
           cgns_preserved_multizone_roundtrip,
           cgns_preserved_multizone_roundtrip_domain,
           context
         ) == sqmesh::base::StatusCode::ok,
         "domain_snapshot should resolve the preserved multi-zone CGNS round-trip mesh"
       )) {
      return EXIT_FAILURE;
    }
    if(!verify_cgns_multi_zone_domain(
         cgns_preserved_multizone_roundtrip_domain,
         "sqmesh_base",
         ""
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::MeshHandle cgns_multibase_mesh = sqmesh::invalid_handle;
    if(!expect(
         sqmesh::mesh::import_cgns(cgns_multibase_input.string(), cgns_multibase_mesh, {}, context) ==
           sqmesh::base::StatusCode::ok,
         "import_cgns should flatten the supported multi-base CGNS subset into base-qualified entity_groups"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::Domain cgns_multibase_domain;
    if(!expect(
         sqmesh::mesh::domain_snapshot(cgns_multibase_mesh, cgns_multibase_domain, context) ==
           sqmesh::base::StatusCode::ok,
         "domain_snapshot should materialize flattened multi-base CGNS data into Domain/EntityGroup storage"
       )) {
      return EXIT_FAILURE;
    }
    if(!verify_cgns_multi_base_domain(
         cgns_multibase_domain,
         "sqmesh_mesh_io_multibase_input"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::CgnsExportOptions cgns_multibase_flattened_export_options;
    cgns_multibase_flattened_export_options.base_name = "sqmesh_base";
    cgns_multibase_flattened_export_options.zone_name =
      "sqmesh_cgns_flattened_multi_base";
    if(!expect(
         sqmesh::mesh::export_cgns(
           cgns_multibase_mesh,
           cgns_multibase_output.string(),
           cgns_multibase_flattened_export_options,
           context
         ) == sqmesh::base::StatusCode::ok,
         "export_cgns should keep the default single-zone rewrite path for supported multi-base imports"
       )) {
      return EXIT_FAILURE;
    }
    if(!expect(
         verify_exported_cgns_file(
           cgns_multibase_output,
           "sqmesh_cgns_flattened_multi_base",
           expected_multibase_sections,
           expected_multibase_boundaries
         ),
         "export_cgns should flatten supported multi-base imports into one written zone with reviewable base-qualified sections and boundaries"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::CgnsExportOptions cgns_multibase_export_options;
    cgns_multibase_export_options.base_name = "sqmesh_base";
    cgns_multibase_export_options.write_multi_zone = true;
    if(!expect(
         sqmesh::mesh::export_cgns(
           cgns_multibase_mesh,
           cgns_multibase_preserved_output.string(),
           cgns_multibase_export_options,
           context
         ) == sqmesh::base::StatusCode::unsupported,
         "export_cgns should reject preserved multi-zone export when imported CGNS provenance spans more than one base"
       )) {
      return EXIT_FAILURE;
    }
    if(!expect(
         std::string(sqmesh::base::last_error_message()).find("multiple imported bases") !=
           std::string::npos,
         "multi-base CGNS export rejection should report that preserved multi-base export is still unsupported"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::MeshHandle cgns_interface_mesh = sqmesh::invalid_handle;
    if(!expect(
         store_manual_mesh(
           make_cgns_supported_multi_zone_interface_domain(),
           "manual cgns multi-zone interface",
           cgns_interface_mesh,
           context
         ),
         "mesh_io_smoke should be able to publish a manual shared-node CGNS multi-zone interface mesh"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::CgnsExportOptions cgns_interface_export_options;
    cgns_interface_export_options.base_name = "sqmesh_base";
    cgns_interface_export_options.write_multi_zone = true;
    if(!expect(
         sqmesh::mesh::export_cgns(
           cgns_interface_mesh,
           cgns_multizone_interface_output.string(),
           cgns_interface_export_options,
           context
         ) == sqmesh::base::StatusCode::ok,
         "export_cgns should reconstruct the bounded shared-node multi-zone interface subset when write_multi_zone is enabled"
       )) {
      return EXIT_FAILURE;
    }
    if(!expect(
         verify_exported_multi_zone_cgns_file(
           cgns_multizone_interface_output,
           expected_interface_multizone_zones
         ),
         "export_cgns should write reviewable interface connectivity while preserving supported boundary BC_t passthrough"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::MeshHandle cgns_patch_interface_mesh =
      sqmesh::invalid_handle;
    if(!expect(
         store_manual_mesh(
           make_cgns_multi_face_interface_domain(),
           "manual cgns multi-face interface",
           cgns_patch_interface_mesh,
           context
         ),
         "mesh_io_smoke should be able to publish a manual shared-node multi-face CGNS interface mesh"
       )) {
      return EXIT_FAILURE;
    }
    if(!expect(
         sqmesh::mesh::export_cgns(
           cgns_patch_interface_mesh,
           cgns_multizone_patch_interface_output.string(),
           cgns_interface_export_options,
           context
         ) == sqmesh::base::StatusCode::ok,
         "export_cgns should reconstruct shared-node CGNS patch interfaces that span multiple interface faces"
       )) {
      return EXIT_FAILURE;
    }
    if(!expect(
         verify_exported_multi_zone_cgns_file(
           cgns_multizone_patch_interface_output,
           expected_patch_interface_multizone_zones
         ),
         "export_cgns should write one vertex connectivity record per zone even when the reconstructed interface spans multiple faces"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::MeshHandle cgns_unsupported_interface_mesh =
      sqmesh::invalid_handle;
    if(!expect(
         store_manual_mesh(
           make_cgns_unsupported_multi_zone_interface_domain(),
           "manual cgns unsupported multi-zone interface",
           cgns_unsupported_interface_mesh,
           context
         ),
         "mesh_io_smoke should be able to publish a manual non-shared-node CGNS interface mesh"
       )) {
      return EXIT_FAILURE;
    }
    if(!expect(
         sqmesh::mesh::export_cgns(
           cgns_unsupported_interface_mesh,
           cgns_multizone_interface_unsupported_output.string(),
           cgns_interface_export_options,
           context
         ) == sqmesh::base::StatusCode::unsupported,
         "export_cgns should keep non-shared-node multi-zone interfaces outside the bounded supported subset"
       )) {
      return EXIT_FAILURE;
    }
    if(!expect(
         std::string(sqmesh::base::last_error_message()).find("shared-node interfaces") !=
           std::string::npos,
         "unsupported multi-zone interface export should report the bounded shared-node diagnostic"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::MeshHandle cgns_surface_interface_mesh =
      sqmesh::invalid_handle;
    if(!expect(
         store_manual_mesh(
           make_cgns_surface_multi_zone_interface_domain(),
           "manual cgns surface multi-zone interface",
           cgns_surface_interface_mesh,
           context
         ),
         "mesh_io_smoke should be able to publish a manual shared-edge CGNS surface interface mesh"
       )) {
      return EXIT_FAILURE;
    }
    if(!expect(
         sqmesh::mesh::export_cgns(
           cgns_surface_interface_mesh,
           cgns_surface_multizone_interface_output.string(),
           cgns_interface_export_options,
           context
         ) == sqmesh::base::StatusCode::ok,
         "export_cgns should reconstruct the bounded shared-node surface-zone interface subset when write_multi_zone is enabled"
       )) {
      return EXIT_FAILURE;
    }
    if(!expect(
         verify_exported_multi_zone_cgns_file(
           cgns_surface_multizone_interface_output,
           expected_surface_interface_multizone_zones
         ),
         "export_cgns should write reviewable edge-based interface connectivity for supported surface zones"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::MeshHandle cgns_surface_polyline_interface_mesh =
      sqmesh::invalid_handle;
    if(!expect(
         store_manual_mesh(
           make_cgns_surface_multi_edge_interface_domain(),
           "manual cgns surface multi-edge interface",
           cgns_surface_polyline_interface_mesh,
           context
         ),
         "mesh_io_smoke should be able to publish a manual shared-edge polyline CGNS interface mesh"
       )) {
      return EXIT_FAILURE;
    }
    if(!expect(
         sqmesh::mesh::export_cgns(
           cgns_surface_polyline_interface_mesh,
           cgns_surface_multizone_polyline_interface_output.string(),
           cgns_interface_export_options,
           context
         ) == sqmesh::base::StatusCode::ok,
         "export_cgns should reconstruct shared-edge CGNS surface interfaces that span multiple interface edges"
       )) {
      return EXIT_FAILURE;
    }
    if(!expect(
         verify_exported_multi_zone_cgns_file(
           cgns_surface_multizone_polyline_interface_output,
           expected_surface_polyline_interface_multizone_zones
         ),
         "export_cgns should keep surface-zone interface point lists stable when the shared polyline spans multiple edges"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::MeshHandle unsupported_cgns_mesh = sqmesh::invalid_handle;
    if(!expect(
         sqmesh::mesh::import_cgns(cgns_unsupported.string(), unsupported_cgns_mesh, {}, context) ==
           sqmesh::base::StatusCode::unsupported,
         "import_cgns should reject unsupported multi-base files cleanly"
       )) {
      return EXIT_FAILURE;
    }
    if(!expect(
         std::string(sqmesh::base::last_error_message()).find(
           "share one cell dimension and one physical coordinate dimension"
         ) != std::string::npos,
         "unsupported multi-base import should report the bounded same-dimension requirement"
       )) {
      return EXIT_FAILURE;
    }

    sqmesh::mesh::MeshHandle unsupported_cgns_physical_dim_mesh =
      sqmesh::invalid_handle;
    if(!expect(
         sqmesh::mesh::import_cgns(
           cgns_unsupported_physical_dim.string(),
           unsupported_cgns_physical_dim_mesh,
           {},
           context
         ) == sqmesh::base::StatusCode::unsupported,
         "import_cgns should reject multi-base files whose bases disagree only on physical coordinate dimension"
       )) {
      return EXIT_FAILURE;
    }
    if(!expect(
         std::string(sqmesh::base::last_error_message()).find(
           "share one cell dimension and one physical coordinate dimension"
         ) != std::string::npos,
         "unsupported multi-base physical-dimension import should report the bounded same-dimension requirement"
       )) {
      return EXIT_FAILURE;
    }
#else
    std::fprintf(
      stderr,
      "mesh_io_smoke: CGNS IO is enabled but the test was built without CGNS support.\n"
    );
    return EXIT_FAILURE;
#endif
  }

  sqmesh::geo::ModelHandle model = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::geo::create_placeholder_model(model, context) ==
         sqmesh::base::StatusCode::ok,
       "create_placeholder_model should seed a source model for dummy meshing"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::mesh::ParameterDictionary parameters;
  parameters.set_number("minimum_length", 0.25);
  parameters.set_number("growth_rate", 1.2);
  parameters.set_text("element_type", "tetra");

  sqmesh::mesh::MeshHandle volume_mesh = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::create_volume_mesh(
         model,
         "Dummy Mesher",
         parameters,
         volume_mesh,
         context
       ) == sqmesh::base::StatusCode::ok,
       "create_volume_mesh should produce a tetra mesh for the OBJ export negative check"
     )) {
    return EXIT_FAILURE;
  }

  if(!expect(
       sqmesh::mesh::export_obj(volume_mesh, obj_output.string(), {}, context) ==
         sqmesh::base::StatusCode::unsupported,
       "export_obj should stay honest about only supporting surface meshes"
     )) {
    return EXIT_FAILURE;
  }

  if(!expect(
       sqmesh::base::shutdown_all() == sqmesh::base::StatusCode::ok,
       "shutdown_all should release the runtime after mesh IO checks"
     )) {
    return EXIT_FAILURE;
  }

  std::filesystem::remove(msh_input);
  std::filesystem::remove(msh_output);
  std::filesystem::remove(msh_binary_input);
  std::filesystem::remove(msh_binary_output);
  std::filesystem::remove(msh41_input);
  std::filesystem::remove(msh41_output);
  std::filesystem::remove(msh41_binary_input);
  std::filesystem::remove(msh41_binary_output);
  std::filesystem::remove(msh41_binary_32bit_input);
  std::filesystem::remove(msh41_binary_32bit_output);
  std::filesystem::remove(msh41_partitioned_input);
  std::filesystem::remove(msh41_partitioned_output);
  std::filesystem::remove(msh41_partitioned_binary_input);
  std::filesystem::remove(msh41_partitioned_binary_output);
  std::filesystem::remove(msh41_shared_physical_input);
  std::filesystem::remove(msh41_shared_physical_output);
  std::filesystem::remove(msh41_shared_physical_binary_output);
  std::filesystem::remove(obj_input);
  std::filesystem::remove(obj_output);
  std::filesystem::remove(obj_attribute_input);
  std::filesystem::remove(obj_attribute_output);
  std::filesystem::remove(obj_concave_input);
  std::filesystem::remove(obj_concave_output);
  std::filesystem::remove(obj_collinear_polygon_input);
  std::filesystem::remove(obj_grouped_input);
  std::filesystem::remove(obj_grouped_output);
  std::filesystem::remove(obj_unsupported_polygon_input);
  std::filesystem::remove(nastran_input);
  std::filesystem::remove(nastran_output);
  std::filesystem::remove(nastran_fixed_input);
  std::filesystem::remove(nastran_broad_input);
  std::filesystem::remove(nastran_broad_output);
  std::filesystem::remove(nastran_broad_fixed_input);
  std::filesystem::remove(nastran_broad_fixed_output);
  std::filesystem::remove(nastran_broad_large_input);
  std::filesystem::remove(nastran_broad_large_output);
  std::filesystem::remove(nastran_unsupported);
  std::filesystem::remove(cgns_input);
  std::filesystem::remove(cgns_output);
  std::filesystem::remove(cgns_mixed_input);
  std::filesystem::remove(cgns_mixed_output);
  std::filesystem::remove(cgns_multizone_input);
  std::filesystem::remove(cgns_multizone_output);
  std::filesystem::remove(cgns_multizone_preserved_output);
  std::filesystem::remove(cgns_multibase_input);
  std::filesystem::remove(cgns_multibase_output);
  std::filesystem::remove(cgns_multibase_preserved_output);
  std::filesystem::remove(cgns_multizone_interface_output);
  std::filesystem::remove(cgns_surface_multizone_interface_output);
  std::filesystem::remove(cgns_multizone_patch_interface_output);
  std::filesystem::remove(cgns_surface_multizone_polyline_interface_output);
  std::filesystem::remove(cgns_multizone_interface_unsupported_output);
  std::filesystem::remove(cgns_unsupported_physical_dim);
  std::filesystem::remove(cgns_unsupported);
  return EXIT_SUCCESS;
}
