// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "bindings.hpp"
#include "exceptions.hpp"

#include "sqmesh/mesh/api.hpp"

#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <cstdint>
#include <string>

namespace {

std::size_t hash_entity_ref(sqmesh::mesh::EntityRef entity) noexcept
{
  return (static_cast<std::size_t>(entity.entity_group) << 32U) ^
         static_cast<std::size_t>(entity.index);
}

py::object cast_parameter_value(const sqmesh::mesh::ParameterValue &value)
{
  if(value.is_boolean()) {
    return py::bool_(value.boolean());
  }
  if(value.is_integer()) {
    return py::int_(value.integer());
  }
  if(value.is_text()) {
    return py::str(value.text());
  }
  if(value.is_number()) {
    return py::float_(value.number());
  }
  return py::none();
}

} // namespace

void bind_mesh(py::module_ &module)
{
  namespace mesh = sqmesh::mesh;
  namespace base = sqmesh::base;
  namespace geo = sqmesh::geo;

  auto mesh_module = module.def_submodule("mesh", "SQMesh mesh bindings.");
  mesh_module.attr("invalid_index") = py::int_(mesh::invalid_index);
  mesh_module.attr("entity_order_mask") = py::int_(mesh::entity_order_mask);
  mesh_module.attr("entity_boundary_mask") = py::int_(mesh::entity_boundary_mask);
  mesh_module.attr("entity_interface_mask") = py::int_(mesh::entity_interface_mask);
  mesh_module.attr("entity_meshed_mask") = py::int_(mesh::entity_meshed_mask);
  mesh_module.attr("entity_live_mask") = py::int_(mesh::entity_live_mask);

  py::enum_<mesh::EntityOrder>(mesh_module, "EntityOrder")
    .value("node", mesh::EntityOrder::node)
    .value("edge", mesh::EntityOrder::edge)
    .value("face", mesh::EntityOrder::face)
    .value("cell", mesh::EntityOrder::cell);

  py::enum_<mesh::FaceSide>(mesh_module, "FaceSide")
    .value("left", mesh::FaceSide::left)
    .value("right", mesh::FaceSide::right);

  py::enum_<mesh::EntityKind>(mesh_module, "EntityKind")
    .value("invalid", mesh::EntityKind::invalid)
    .value("node_point", mesh::EntityKind::node_point)
    .value("edge_line", mesh::EntityKind::edge_line)
    .value("face_triangle", mesh::EntityKind::face_triangle)
    .value("face_quad", mesh::EntityKind::face_quad)
    .value("cell_tetra", mesh::EntityKind::cell_tetra)
    .value("cell_hexa", mesh::EntityKind::cell_hexa);

  py::enum_<mesh::ParameterType>(mesh_module, "ParameterType")
    .value("empty", mesh::ParameterType::empty)
    .value("integer", mesh::ParameterType::integer)
    .value("number", mesh::ParameterType::number)
    .value("boolean", mesh::ParameterType::boolean)
    .value("text", mesh::ParameterType::text);

  py::enum_<mesh::QualityStatus>(mesh_module, "QualityStatus")
    .value("unsupported", mesh::QualityStatus::unsupported)
    .value("valid", mesh::QualityStatus::valid)
    .value("degenerate", mesh::QualityStatus::degenerate)
    .value("inverted", mesh::QualityStatus::inverted);

  py::enum_<mesh::EntityGroupSemantic>(mesh_module, "EntityGroupSemantic")
    .value("unspecified", mesh::EntityGroupSemantic::unspecified)
    .value("node", mesh::EntityGroupSemantic::node)
    .value("interior", mesh::EntityGroupSemantic::interior)
    .value("boundary", mesh::EntityGroupSemantic::boundary)
    .value("interface", mesh::EntityGroupSemantic::interface)
    .value("region", mesh::EntityGroupSemantic::region);

  py::enum_<mesh::EntityGroupRole>(mesh_module, "EntityGroupRole")
    .value("computational", mesh::EntityGroupRole::computational)
    .value("geometric_proxy", mesh::EntityGroupRole::geometric_proxy)
    .value("annotation", mesh::EntityGroupRole::annotation);

  py::enum_<mesh::EntityGroupImportFormat>(mesh_module, "EntityGroupImportFormat")
    .value("none", mesh::EntityGroupImportFormat::none)
    .value("cgns", mesh::EntityGroupImportFormat::cgns)
    .value("nastran", mesh::EntityGroupImportFormat::nastran);

  py::enum_<mesh::NastranEntityGroupSourceCard>(mesh_module, "NastranEntityGroupSourceCard")
    .value("unspecified", mesh::NastranEntityGroupSourceCard::unspecified)
    .value("cbar", mesh::NastranEntityGroupSourceCard::cbar)
    .value("cbeam", mesh::NastranEntityGroupSourceCard::cbeam)
    .value("crod", mesh::NastranEntityGroupSourceCard::crod)
    .value("ctria3", mesh::NastranEntityGroupSourceCard::ctria3)
    .value("cquad4", mesh::NastranEntityGroupSourceCard::cquad4)
    .value("ctetra", mesh::NastranEntityGroupSourceCard::ctetra)
    .value("mixed", mesh::NastranEntityGroupSourceCard::mixed);

  py::enum_<mesh::NastranPropertyCard>(mesh_module, "NastranPropertyCard")
    .value("unspecified", mesh::NastranPropertyCard::unspecified)
    .value("prod", mesh::NastranPropertyCard::prod)
    .value("pshell", mesh::NastranPropertyCard::pshell)
    .value("psolid", mesh::NastranPropertyCard::psolid);

  py::enum_<mesh::NastranMaterialCard>(mesh_module, "NastranMaterialCard")
    .value("unspecified", mesh::NastranMaterialCard::unspecified)
    .value("mat1", mesh::NastranMaterialCard::mat1);

  py::enum_<mesh::MshFormatVersion>(mesh_module, "MshFormatVersion")
    .value("gmsh22_ascii", mesh::MshFormatVersion::gmsh22_ascii)
    .value("gmsh22_binary", mesh::MshFormatVersion::gmsh22_binary)
    .value("gmsh41_ascii", mesh::MshFormatVersion::gmsh41_ascii)
    .value("gmsh41_binary", mesh::MshFormatVersion::gmsh41_binary);

  py::enum_<mesh::NastranFieldFormat>(mesh_module, "NastranFieldFormat")
    .value("free_field", mesh::NastranFieldFormat::free_field)
    .value("small_fixed", mesh::NastranFieldFormat::small_fixed)
    .value("large_fixed", mesh::NastranFieldFormat::large_fixed);

  py::class_<mesh::EntityRef>(mesh_module, "EntityRef")
    .def(py::init<>())
    .def(
      py::init<mesh::EntityGroupIndex, std::uint32_t>(),
      py::arg("entity_group"),
      py::arg("index")
    )
    .def_readwrite("entity_group", &mesh::EntityRef::entity_group)
    .def_readwrite("index", &mesh::EntityRef::index)
    .def("is_valid", &mesh::is_valid)
    .def("__bool__", &mesh::is_valid)
    .def(
      "__repr__",
      [](mesh::EntityRef self) {
        return "<sqmesh.mesh.EntityRef entity_group=" +
               std::to_string(self.entity_group) +
               " index=" + std::to_string(self.index) + ">";
      }
    )
    .def("__hash__", &hash_entity_ref)
    .def(py::self == py::self)
    .def(py::self != py::self);

  py::class_<mesh::ConnectivitySpan>(mesh_module, "ConnectivitySpan")
    .def(py::init<>())
    .def_readwrite("offset", &mesh::ConnectivitySpan::offset)
    .def_readwrite("count", &mesh::ConnectivitySpan::count);

  py::class_<mesh::EntityHeader>(mesh_module, "EntityHeader")
    .def(py::init<>())
    .def_readwrite("id", &mesh::EntityHeader::id)
    .def_readwrite("flags", &mesh::EntityHeader::flags)
    .def_readwrite("index", &mesh::EntityHeader::index)
    .def_readwrite("kind", &mesh::EntityHeader::kind)
    .def_readwrite("reserved", &mesh::EntityHeader::reserved)
    .def_readwrite("entity_group", &mesh::EntityHeader::entity_group);

  py::class_<mesh::Node>(mesh_module, "Node")
    .def(py::init<>())
    .def_readwrite("header", &mesh::Node::header)
    .def_readwrite("coordinates", &mesh::Node::coordinates);

  py::class_<mesh::Edge>(mesh_module, "Edge")
    .def(py::init<>())
    .def_readwrite("header", &mesh::Edge::header)
    .def_readwrite("node_span", &mesh::Edge::node_span)
    .def_readwrite("left_face", &mesh::Edge::left_face)
    .def_readwrite("right_face", &mesh::Edge::right_face);

  py::class_<mesh::Face>(mesh_module, "Face")
    .def(py::init<>())
    .def_readwrite("header", &mesh::Face::header)
    .def_readwrite("node_span", &mesh::Face::node_span)
    .def_readwrite("left_cell", &mesh::Face::left_cell)
    .def_readwrite("right_cell", &mesh::Face::right_cell);

  py::class_<mesh::Cell>(mesh_module, "Cell")
    .def(py::init<>())
    .def_readwrite("header", &mesh::Cell::header)
    .def_readwrite("node_span", &mesh::Cell::node_span)
    .def_readwrite("face_span", &mesh::Cell::face_span);

  py::class_<mesh::MeshCoreLayout>(mesh_module, "MeshCoreLayout")
    .def(py::init<>())
    .def_readwrite("entity_ref_size", &mesh::MeshCoreLayout::entity_ref_size)
    .def_readwrite(
      "entity_ref_alignment",
      &mesh::MeshCoreLayout::entity_ref_alignment
    )
    .def_readwrite(
      "connectivity_span_size",
      &mesh::MeshCoreLayout::connectivity_span_size
    )
    .def_readwrite(
      "connectivity_span_alignment",
      &mesh::MeshCoreLayout::connectivity_span_alignment
    )
    .def_readwrite(
      "entity_header_size",
      &mesh::MeshCoreLayout::entity_header_size
    )
    .def_readwrite(
      "entity_header_alignment",
      &mesh::MeshCoreLayout::entity_header_alignment
    )
    .def_readwrite("node_size", &mesh::MeshCoreLayout::node_size)
    .def_readwrite("node_alignment", &mesh::MeshCoreLayout::node_alignment)
    .def_readwrite("edge_size", &mesh::MeshCoreLayout::edge_size)
    .def_readwrite("edge_alignment", &mesh::MeshCoreLayout::edge_alignment)
    .def_readwrite("face_size", &mesh::MeshCoreLayout::face_size)
    .def_readwrite("face_alignment", &mesh::MeshCoreLayout::face_alignment)
    .def_readwrite("cell_size", &mesh::MeshCoreLayout::cell_size)
    .def_readwrite("cell_alignment", &mesh::MeshCoreLayout::cell_alignment)
    .def_readwrite(
      "node_header_offset",
      &mesh::MeshCoreLayout::node_header_offset
    )
    .def_readwrite(
      "edge_header_offset",
      &mesh::MeshCoreLayout::edge_header_offset
    )
    .def_readwrite(
      "edge_node_span_offset",
      &mesh::MeshCoreLayout::edge_node_span_offset
    )
    .def_readwrite(
      "edge_left_face_offset",
      &mesh::MeshCoreLayout::edge_left_face_offset
    )
    .def_readwrite(
      "edge_right_face_offset",
      &mesh::MeshCoreLayout::edge_right_face_offset
    )
    .def_readwrite(
      "face_header_offset",
      &mesh::MeshCoreLayout::face_header_offset
    )
    .def_readwrite(
      "face_node_span_offset",
      &mesh::MeshCoreLayout::face_node_span_offset
    )
    .def_readwrite(
      "face_left_cell_offset",
      &mesh::MeshCoreLayout::face_left_cell_offset
    )
    .def_readwrite(
      "face_right_cell_offset",
      &mesh::MeshCoreLayout::face_right_cell_offset
    )
    .def_readwrite(
      "cell_header_offset",
      &mesh::MeshCoreLayout::cell_header_offset
    )
    .def_readwrite(
      "cell_node_span_offset",
      &mesh::MeshCoreLayout::cell_node_span_offset
    )
    .def_readwrite(
      "cell_face_span_offset",
      &mesh::MeshCoreLayout::cell_face_span_offset
    );

  py::class_<mesh::MeshSummary>(mesh_module, "MeshSummary")
    .def(py::init<>())
    .def_readwrite("node_count", &mesh::MeshSummary::node_count)
    .def_readwrite("edge_count", &mesh::MeshSummary::edge_count)
    .def_readwrite("face_count", &mesh::MeshSummary::face_count)
    .def_readwrite("cell_count", &mesh::MeshSummary::cell_count)
    .def_readwrite(
      "source_topology_revision",
      &mesh::MeshSummary::source_topology_revision
    );

  py::class_<mesh::ParameterValue>(mesh_module, "ParameterValue")
    .def(py::init<>())
    .def(py::init<bool>(), py::arg("value"))
    .def(py::init<std::int64_t>(), py::arg("value"))
    .def(py::init<double>(), py::arg("value"))
    .def(py::init<std::string>(), py::arg("value"))
    .def("type", &mesh::ParameterValue::type)
    .def("empty", &mesh::ParameterValue::empty)
    .def("is_integer", &mesh::ParameterValue::is_integer)
    .def("is_number", &mesh::ParameterValue::is_number)
    .def("is_boolean", &mesh::ParameterValue::is_boolean)
    .def("is_text", &mesh::ParameterValue::is_text)
    .def(
      "integer",
      &mesh::ParameterValue::integer,
      py::arg("fallback") = std::int64_t(0)
    )
    .def("number", &mesh::ParameterValue::number, py::arg("fallback") = 0.0)
    .def(
      "boolean",
      &mesh::ParameterValue::boolean,
      py::arg("fallback") = false
    )
    .def("text", [](const mesh::ParameterValue &self) {
      return std::string(self.text());
    })
    .def("to_python", &cast_parameter_value)
    .def("__repr__", [](const mesh::ParameterValue &self) {
      auto value = cast_parameter_value(self);
      return "<sqmesh.mesh.ParameterValue " + py::repr(value).cast<std::string>() + ">";
    });

  py::class_<mesh::ParameterDictionary>(mesh_module, "ParameterDictionary")
    .def(py::init<>())
    .def("set", &mesh::ParameterDictionary::set, py::arg("key"), py::arg("value"))
    .def(
      "set_integer",
      &mesh::ParameterDictionary::set_integer,
      py::arg("key"),
      py::arg("value")
    )
    .def(
      "set_number",
      &mesh::ParameterDictionary::set_number,
      py::arg("key"),
      py::arg("value")
    )
    .def(
      "set_boolean",
      &mesh::ParameterDictionary::set_boolean,
      py::arg("key"),
      py::arg("value")
    )
    .def(
      "set_text",
      &mesh::ParameterDictionary::set_text,
      py::arg("key"),
      py::arg("value")
    )
    .def("contains", &mesh::ParameterDictionary::contains, py::arg("key"))
    .def("__contains__", &mesh::ParameterDictionary::contains, py::arg("key"))
    .def("size", &mesh::ParameterDictionary::size)
    .def("__len__", &mesh::ParameterDictionary::size)
    .def(
      "find",
      [](const mesh::ParameterDictionary &self, std::string_view key) -> py::object {
        const auto *value = self.find(key);
        if(value == nullptr) {
          return py::none();
        }
        return py::cast(*value);
      },
      py::arg("key")
    )
    .def(
      "get",
      [](const mesh::ParameterDictionary &self, std::string_view key) -> py::object {
        const auto *value = self.find(key);
        if(value == nullptr) {
          return py::none();
        }
        return cast_parameter_value(*value);
      },
      py::arg("key")
    )
    .def(
      "try_get_integer",
      [](const mesh::ParameterDictionary &self, std::string_view key) -> py::object {
        std::int64_t value = 0;
        if(!self.try_get_integer(key, value)) {
          return py::none();
        }
        return py::int_(value);
      },
      py::arg("key")
    )
    .def(
      "try_get_number",
      [](const mesh::ParameterDictionary &self, std::string_view key) -> py::object {
        double value = 0.0;
        if(!self.try_get_number(key, value)) {
          return py::none();
        }
        return py::float_(value);
      },
      py::arg("key")
    )
    .def(
      "try_get_boolean",
      [](const mesh::ParameterDictionary &self, std::string_view key) -> py::object {
        bool value = false;
        if(!self.try_get_boolean(key, value)) {
          return py::none();
        }
        return py::bool_(value);
      },
      py::arg("key")
    )
    .def(
      "try_get_text",
      [](const mesh::ParameterDictionary &self, std::string_view key) -> py::object {
        std::string_view value;
        if(!self.try_get_text(key, value)) {
          return py::none();
        }
        return py::str(value);
      },
      py::arg("key")
    );

  py::class_<mesh::MeshSizeControl>(mesh_module, "MeshSizeControl")
    .def(py::init<>())
    .def_readwrite("entity", &mesh::MeshSizeControl::entity)
    .def_readwrite("target_size", &mesh::MeshSizeControl::target_size);

  py::class_<mesh::MeshSizeControls>(mesh_module, "MeshSizeControls")
    .def(py::init<>())
    .def_readwrite("topology_revision", &mesh::MeshSizeControls::topology_revision)
    .def_readwrite("local_sizes", &mesh::MeshSizeControls::local_sizes)
    .def("add_local_size", &mesh::MeshSizeControls::add_local_size)
    .def("empty", &mesh::MeshSizeControls::empty);

  py::class_<mesh::MeshingOptions>(mesh_module, "MeshingOptions")
    .def(py::init<>())
    .def_readwrite("parameters", &mesh::MeshingOptions::parameters)
    .def_readwrite("size_controls", &mesh::MeshingOptions::size_controls);

  py::class_<mesh::DomainStatistics>(mesh_module, "DomainStatistics")
    .def(py::init<>())
    .def_readwrite("summary", &mesh::DomainStatistics::summary)
    .def_readwrite("entity_group_count", &mesh::DomainStatistics::entity_group_count)
    .def_readwrite(
      "boundary_edge_count",
      &mesh::DomainStatistics::boundary_edge_count
    )
    .def_readwrite(
      "interior_edge_count",
      &mesh::DomainStatistics::interior_edge_count
    )
    .def_readwrite(
      "boundary_face_count",
      &mesh::DomainStatistics::boundary_face_count
    )
    .def_readwrite(
      "interior_face_count",
      &mesh::DomainStatistics::interior_face_count
    )
    .def_readwrite("line_edge_count", &mesh::DomainStatistics::line_edge_count)
    .def_readwrite(
      "triangle_face_count",
      &mesh::DomainStatistics::triangle_face_count
    )
    .def_readwrite(
      "tetra_cell_count",
      &mesh::DomainStatistics::tetra_cell_count
    );

  py::class_<mesh::QualityMetricSummary>(mesh_module, "QualityMetricSummary")
    .def(py::init<>())
    .def_readwrite("count", &mesh::QualityMetricSummary::count)
    .def_readwrite("minimum", &mesh::QualityMetricSummary::minimum)
    .def_readwrite("maximum", &mesh::QualityMetricSummary::maximum)
    .def_readwrite("average", &mesh::QualityMetricSummary::average);

  py::class_<mesh::ElementQuality>(mesh_module, "ElementQuality")
    .def(py::init<>())
    .def_readwrite("entity", &mesh::ElementQuality::entity)
    .def_readwrite("kind", &mesh::ElementQuality::kind)
    .def_readwrite("status", &mesh::ElementQuality::status)
    .def_readwrite("supported", &mesh::ElementQuality::supported)
    .def_readwrite("degenerate", &mesh::ElementQuality::degenerate)
    .def_readwrite("inverted", &mesh::ElementQuality::inverted)
    .def_readwrite(
      "jacobian_is_signed",
      &mesh::ElementQuality::jacobian_is_signed
    )
    .def_readwrite("jacobian", &mesh::ElementQuality::jacobian)
    .def_readwrite("skewness", &mesh::ElementQuality::skewness)
    .def_readwrite("aspect_ratio", &mesh::ElementQuality::aspect_ratio)
    .def_readwrite("radius_ratio", &mesh::ElementQuality::radius_ratio)
    .def_readwrite("min_angle", &mesh::ElementQuality::min_angle)
    .def_readwrite("max_angle", &mesh::ElementQuality::max_angle);

  py::class_<mesh::KindQualitySummary>(mesh_module, "KindQualitySummary")
    .def(py::init<>())
    .def_readwrite("kind", &mesh::KindQualitySummary::kind)
    .def_readwrite(
      "supported_element_count",
      &mesh::KindQualitySummary::supported_element_count
    )
    .def_readwrite(
      "valid_element_count",
      &mesh::KindQualitySummary::valid_element_count
    )
    .def_readwrite(
      "degenerate_element_count",
      &mesh::KindQualitySummary::degenerate_element_count
    )
    .def_readwrite(
      "inverted_element_count",
      &mesh::KindQualitySummary::inverted_element_count
    )
    .def_readwrite("jacobian", &mesh::KindQualitySummary::jacobian)
    .def_readwrite("skewness", &mesh::KindQualitySummary::skewness)
    .def_readwrite("aspect_ratio", &mesh::KindQualitySummary::aspect_ratio)
    .def_readwrite("radius_ratio", &mesh::KindQualitySummary::radius_ratio)
    .def_readwrite("min_angle", &mesh::KindQualitySummary::min_angle)
    .def_readwrite("max_angle", &mesh::KindQualitySummary::max_angle);

  py::class_<mesh::MeshQualityReport>(mesh_module, "MeshQualityReport")
    .def(py::init<>())
    .def_readwrite("summary", &mesh::MeshQualityReport::summary)
    .def_readwrite(
      "supported_element_count",
      &mesh::MeshQualityReport::supported_element_count
    )
    .def_readwrite(
      "unsupported_element_count",
      &mesh::MeshQualityReport::unsupported_element_count
    )
    .def_readwrite(
      "valid_element_count",
      &mesh::MeshQualityReport::valid_element_count
    )
    .def_readwrite(
      "degenerate_element_count",
      &mesh::MeshQualityReport::degenerate_element_count
    )
    .def_readwrite(
      "inverted_element_count",
      &mesh::MeshQualityReport::inverted_element_count
    )
    .def_readwrite("kinds", &mesh::MeshQualityReport::kinds)
    .def_readwrite("elements", &mesh::MeshQualityReport::elements);

  py::class_<mesh::CgnsEntityGroupImportInfo>(mesh_module, "CgnsEntityGroupImportInfo")
    .def(py::init<>())
    .def_readwrite("base_index", &mesh::CgnsEntityGroupImportInfo::base_index)
    .def_readwrite("base_name", &mesh::CgnsEntityGroupImportInfo::base_name)
    .def_readwrite("zone_index", &mesh::CgnsEntityGroupImportInfo::zone_index)
    .def_readwrite("zone_name", &mesh::CgnsEntityGroupImportInfo::zone_name)
    .def_readwrite("local_name", &mesh::CgnsEntityGroupImportInfo::local_name)
    .def_readwrite(
      "bc_type_value",
      &mesh::CgnsEntityGroupImportInfo::bc_type_value
    );

  py::class_<mesh::NastranEntityGroupImportInfo>(mesh_module, "NastranEntityGroupImportInfo")
    .def(py::init<>())
    .def_readwrite(
      "source_card",
      &mesh::NastranEntityGroupImportInfo::source_card
    )
    .def_readwrite(
      "property_card",
      &mesh::NastranEntityGroupImportInfo::property_card
    )
    .def_readwrite(
      "material_card",
      &mesh::NastranEntityGroupImportInfo::material_card
    )
    .def_readwrite("material_id", &mesh::NastranEntityGroupImportInfo::material_id)
    .def_readwrite("rod_area", &mesh::NastranEntityGroupImportInfo::rod_area)
    .def_readwrite(
      "shell_thickness",
      &mesh::NastranEntityGroupImportInfo::shell_thickness
    )
    .def_readwrite(
      "youngs_modulus",
      &mesh::NastranEntityGroupImportInfo::youngs_modulus
    )
    .def_readwrite(
      "shear_modulus",
      &mesh::NastranEntityGroupImportInfo::shear_modulus
    )
    .def_readwrite(
      "poisson_ratio",
      &mesh::NastranEntityGroupImportInfo::poisson_ratio
    );

  py::class_<mesh::EntityGroupImportInfo>(mesh_module, "EntityGroupImportInfo")
    .def(py::init<>())
    .def_readwrite("format", &mesh::EntityGroupImportInfo::format)
    .def_readwrite("cgns", &mesh::EntityGroupImportInfo::cgns)
    .def_readwrite("nastran", &mesh::EntityGroupImportInfo::nastran);

  py::class_<mesh::EntityGroupInfo>(mesh_module, "EntityGroupInfo")
    .def(py::init<>())
    .def_readwrite("id", &mesh::EntityGroupInfo::id)
    .def_readwrite("zone_id", &mesh::EntityGroupInfo::zone_id)
    .def_readwrite(
      "source_entity_tag",
      &mesh::EntityGroupInfo::source_entity_tag
    )
    .def_readwrite("order", &mesh::EntityGroupInfo::order)
    .def_readwrite("boundary", &mesh::EntityGroupInfo::boundary)
    .def_readwrite("default_kind", &mesh::EntityGroupInfo::default_kind)
    .def_readwrite("semantic", &mesh::EntityGroupInfo::semantic)
    .def_readwrite("role", &mesh::EntityGroupInfo::role)
    .def_readwrite(
      "primary_region_zone_id",
      &mesh::EntityGroupInfo::primary_region_zone_id
    )
    .def_readwrite(
      "secondary_region_zone_id",
      &mesh::EntityGroupInfo::secondary_region_zone_id
    )
    .def_readwrite("name", &mesh::EntityGroupInfo::name)
    .def_readwrite("import_info", &mesh::EntityGroupInfo::import_info);

  py::class_<mesh::EntityGroup>(mesh_module, "EntityGroup")
    .def("info", [](const mesh::EntityGroup &self) { return self.info(); })
    .def("import_info", [](const mesh::EntityGroup &self) { return self.import_info(); })
    .def("id", &mesh::EntityGroup::id)
    .def("zone_id", &mesh::EntityGroup::zone_id)
    .def("source_entity_tag", &mesh::EntityGroup::source_entity_tag)
    .def("order", &mesh::EntityGroup::order)
    .def("is_boundary", &mesh::EntityGroup::is_boundary)
    .def("semantic", &mesh::EntityGroup::semantic)
    .def("role", &mesh::EntityGroup::role)
    .def("is_interface", &mesh::EntityGroup::is_interface)
    .def(
      "primary_region_zone_id",
      &mesh::EntityGroup::primary_region_zone_id
    )
    .def(
      "secondary_region_zone_id",
      &mesh::EntityGroup::secondary_region_zone_id
    )
    .def("default_kind", &mesh::EntityGroup::default_kind)
    .def("name", [](const mesh::EntityGroup &self) {
      return std::string(self.name());
    })
    .def("entity_count", &mesh::EntityGroup::entity_count)
    .def("nodes", [](const mesh::EntityGroup &self) { return self.nodes(); })
    .def("edges", [](const mesh::EntityGroup &self) { return self.edges(); })
    .def("faces", [](const mesh::EntityGroup &self) { return self.faces(); })
    .def("cells", [](const mesh::EntityGroup &self) { return self.cells(); })
    .def("node_channel", [](const mesh::EntityGroup &self) { return self.node_channel(); })
    .def("face_channel", [](const mesh::EntityGroup &self) { return self.face_channel(); });

  py::class_<mesh::Domain>(mesh_module, "Domain")
    .def("name", [](const mesh::Domain &self) {
      return std::string(self.name());
    })
    .def("entity_groups", [](const mesh::Domain &self) { return self.entity_groups(); })
    .def("entity_group", [](const mesh::Domain &self, mesh::EntityGroupIndex index) {
      return self.entity_group(index);
    })
    .def(
      "entity_group_count",
      [](const mesh::Domain &self) { return self.entity_group_count(); }
    )
    .def(
      "entity_group_count_by_order",
      [](const mesh::Domain &self, mesh::EntityOrder order) {
        return self.entity_group_count(order);
      },
      py::arg("order")
    )
    .def(
      "entity_group_count_by_role",
      [](const mesh::Domain &self, mesh::EntityGroupRole role) {
        return self.entity_group_count(role);
      },
      py::arg("role")
    )
    .def("source_topology_revision", &mesh::Domain::source_topology_revision)
    .def("node_count", &mesh::Domain::node_count)
    .def("edge_count", &mesh::Domain::edge_count)
    .def("face_count", &mesh::Domain::face_count)
    .def("cell_count", &mesh::Domain::cell_count)
    .def("summary", &mesh::Domain::summary)
    .def("statistics", &mesh::Domain::statistics)
    .def("element_quality", &mesh::Domain::element_quality, py::arg("entity_ref"))
    .def("quality_report", &mesh::Domain::quality_report)
    .def("node", &mesh::Domain::node, py::arg("node_ref"))
    .def("edge", &mesh::Domain::edge, py::arg("edge_ref"))
    .def("face", &mesh::Domain::face, py::arg("face_ref"))
    .def("cell", &mesh::Domain::cell, py::arg("cell_ref"))
    .def(
      "edge_nodes",
      [](const mesh::Domain &self, mesh::EntityRef edge_ref) {
        const auto range = self.edge_nodes(edge_ref);
        return std::vector<mesh::EntityRef>(range.begin(), range.end());
      },
      py::arg("edge_ref")
    )
    .def(
      "face_nodes",
      [](const mesh::Domain &self, mesh::EntityRef face_ref) {
        const auto range = self.face_nodes(face_ref);
        return std::vector<mesh::EntityRef>(range.begin(), range.end());
      },
      py::arg("face_ref")
    )
    .def(
      "cell_nodes",
      [](const mesh::Domain &self, mesh::EntityRef cell_ref) {
        const auto range = self.cell_nodes(cell_ref);
        return std::vector<mesh::EntityRef>(range.begin(), range.end());
      },
      py::arg("cell_ref")
    )
    .def(
      "cell_faces",
      [](const mesh::Domain &self, mesh::EntityRef cell_ref) {
        const auto range = self.cell_faces(cell_ref);
        return std::vector<mesh::EntityRef>(range.begin(), range.end());
      },
      py::arg("cell_ref")
    )
    .def(
      "adjacent_face",
      &mesh::Domain::adjacent_face,
      py::arg("edge_ref"),
      py::arg("side")
    )
    .def(
      "adjacent_cell",
      &mesh::Domain::adjacent_cell,
      py::arg("face_ref"),
      py::arg("side")
    )
    .def(
      "edge_topology_owner",
      &mesh::Domain::edge_topology_owner,
      py::arg("edge_ref")
    )
    .def(
      "face_topology_owner",
      &mesh::Domain::face_topology_owner,
      py::arg("face_ref")
    )
    .def(
      "face_source_entity_tag",
      &mesh::Domain::face_source_entity_tag,
      py::arg("face_ref")
    );

  py::class_<mesh::MshImportOptions>(mesh_module, "MshImportOptions")
    .def(py::init<>())
    .def_readwrite(
      "read_physical_names",
      &mesh::MshImportOptions::read_physical_names
    );

  py::class_<mesh::MshExportOptions>(mesh_module, "MshExportOptions")
    .def(py::init<>())
    .def_readwrite(
      "write_physical_names",
      &mesh::MshExportOptions::write_physical_names
    )
    .def_readwrite(
      "format_version",
      &mesh::MshExportOptions::format_version
    );

  py::class_<mesh::ObjImportOptions>(mesh_module, "ObjImportOptions")
    .def(py::init<>())
    .def_readwrite("read_groups", &mesh::ObjImportOptions::read_groups);

  py::class_<mesh::ObjExportOptions>(mesh_module, "ObjExportOptions")
    .def(py::init<>())
    .def_readwrite("write_groups", &mesh::ObjExportOptions::write_groups)
    .def_readwrite(
      "write_object_name",
      &mesh::ObjExportOptions::write_object_name
    );

  py::class_<mesh::CgnsImportOptions>(mesh_module, "CgnsImportOptions")
    .def(py::init<>())
    .def_readwrite(
      "read_section_names",
      &mesh::CgnsImportOptions::read_section_names
    );

  py::class_<mesh::CgnsExportOptions>(mesh_module, "CgnsExportOptions")
    .def(py::init<>())
    .def_readwrite("base_name", &mesh::CgnsExportOptions::base_name)
    .def_readwrite("zone_name", &mesh::CgnsExportOptions::zone_name)
    .def_readwrite(
      "write_multi_zone",
      &mesh::CgnsExportOptions::write_multi_zone
    );

  py::class_<mesh::NastranImportOptions>(mesh_module, "NastranImportOptions")
    .def(py::init<>())
    .def_readwrite(
      "allow_free_fields",
      &mesh::NastranImportOptions::allow_free_fields
    )
    .def_readwrite(
      "allow_small_fixed_fields",
      &mesh::NastranImportOptions::allow_small_fixed_fields
    )
    .def_readwrite(
      "allow_large_fixed_fields",
      &mesh::NastranImportOptions::allow_large_fixed_fields
    );

  py::class_<mesh::NastranExportOptions>(mesh_module, "NastranExportOptions")
    .def(py::init<>())
    .def_readwrite("field_format", &mesh::NastranExportOptions::field_format)
    .def_readwrite(
      "write_begin_bulk",
      &mesh::NastranExportOptions::write_begin_bulk
    )
    .def_readwrite("write_enddata", &mesh::NastranExportOptions::write_enddata);

  mesh_module.def("module_name", []() {
    return std::string(mesh::module_name());
  });
  mesh_module.def("make_dummy_tetra_domain", &mesh::make_dummy_tetra_domain);
  mesh_module.def(
    "encode_entity_kind",
    &mesh::encode_entity_kind,
    py::arg("order"),
    py::arg("arity")
  );
  mesh_module.def(
    "entity_order_from_kind",
    py::overload_cast<mesh::EntityKind>(&mesh::entity_order),
    py::arg("kind")
  );
  mesh_module.def("entity_arity", &mesh::entity_arity, py::arg("kind"));
  mesh_module.def(
    "is_boundary_entity",
    &mesh::is_boundary_entity,
    py::arg("flags")
  );
  mesh_module.def(
    "is_interface_entity",
    &mesh::is_interface_entity,
    py::arg("flags")
  );
  mesh_module.def(
    "entity_order_from_flags",
    py::overload_cast<std::uint32_t>(&mesh::entity_order),
    py::arg("flags")
  );
  mesh_module.def("mesh_core_layout", &mesh::mesh_core_layout);
  mesh_module.def("expected_mesh_core_layout", &mesh::expected_mesh_core_layout);
  mesh_module.def(
    "same_mesh_core_layout",
    &mesh::same_mesh_core_layout,
    py::arg("lhs"),
    py::arg("rhs")
  );
  mesh_module.def(
    "matches_mesh_core_layout_baseline",
    &mesh::matches_mesh_core_layout_baseline,
    py::arg("layout")
  );
  mesh_module.def(
    "layout_benchmark_bytes_per_cell",
    &mesh::layout_benchmark_bytes_per_cell,
    py::arg("layout")
  );
  mesh_module.def(
    "expected_layout_benchmark_bytes_per_cell",
    &mesh::expected_layout_benchmark_bytes_per_cell
  );
  mesh_module.def(
    "is_algorithm_registered",
    &mesh::is_algorithm_registered,
    py::arg("algorithm_name")
  );
  mesh_module.def("cgns_io_available", &mesh::cgns_io_available);
  mesh_module.def(
    "import_msh",
    [](const std::string &path,
       const mesh::MshImportOptions &options,
       base::ContextHandle context_handle) {
      mesh::MeshHandle mesh_handle = sqmesh::invalid_handle;
      throw_on_status(mesh::import_msh(path, mesh_handle, options, context_handle));
      return mesh_handle;
    },
    py::arg("path"),
    py::arg("options") = mesh::MshImportOptions {},
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  mesh_module.def(
    "export_msh",
    [](mesh::MeshHandle mesh_handle,
       const std::string &path,
       const mesh::MshExportOptions &options,
       base::ContextHandle context_handle) {
      throw_on_status(mesh::export_msh(mesh_handle, path, options, context_handle));
    },
    py::arg("mesh_handle"),
    py::arg("path"),
    py::arg("options") = mesh::MshExportOptions {},
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  mesh_module.def(
    "import_obj",
    [](const std::string &path,
       const mesh::ObjImportOptions &options,
       base::ContextHandle context_handle) {
      mesh::MeshHandle mesh_handle = sqmesh::invalid_handle;
      throw_on_status(mesh::import_obj(path, mesh_handle, options, context_handle));
      return mesh_handle;
    },
    py::arg("path"),
    py::arg("options") = mesh::ObjImportOptions {},
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  mesh_module.def(
    "export_obj",
    [](mesh::MeshHandle mesh_handle,
       const std::string &path,
       const mesh::ObjExportOptions &options,
       base::ContextHandle context_handle) {
      throw_on_status(mesh::export_obj(mesh_handle, path, options, context_handle));
    },
    py::arg("mesh_handle"),
    py::arg("path"),
    py::arg("options") = mesh::ObjExportOptions {},
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  mesh_module.def(
    "import_cgns",
    [](const std::string &path,
       const mesh::CgnsImportOptions &options,
       base::ContextHandle context_handle) {
      mesh::MeshHandle mesh_handle = sqmesh::invalid_handle;
      throw_on_status(mesh::import_cgns(path, mesh_handle, options, context_handle));
      return mesh_handle;
    },
    py::arg("path"),
    py::arg("options") = mesh::CgnsImportOptions {},
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  mesh_module.def(
    "export_cgns",
    [](mesh::MeshHandle mesh_handle,
       const std::string &path,
       const mesh::CgnsExportOptions &options,
       base::ContextHandle context_handle) {
      throw_on_status(mesh::export_cgns(mesh_handle, path, options, context_handle));
    },
    py::arg("mesh_handle"),
    py::arg("path"),
    py::arg("options") = mesh::CgnsExportOptions {},
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  mesh_module.def(
    "import_nastran",
    [](const std::string &path,
       const mesh::NastranImportOptions &options,
       base::ContextHandle context_handle) {
      mesh::MeshHandle mesh_handle = sqmesh::invalid_handle;
      throw_on_status(
        mesh::import_nastran(path, mesh_handle, options, context_handle)
      );
      return mesh_handle;
    },
    py::arg("path"),
    py::arg("options") = mesh::NastranImportOptions {},
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  mesh_module.def(
    "export_nastran",
    [](mesh::MeshHandle mesh_handle,
       const std::string &path,
       const mesh::NastranExportOptions &options,
       base::ContextHandle context_handle) {
      throw_on_status(
        mesh::export_nastran(mesh_handle, path, options, context_handle)
      );
    },
    py::arg("mesh_handle"),
    py::arg("path"),
    py::arg("options") = mesh::NastranExportOptions {},
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  mesh_module.def(
    "create_surface_mesh",
    [](geo::ModelHandle model_handle,
       const std::string &algorithm_name,
       const mesh::MeshingOptions &options,
       base::ContextHandle context_handle) {
      mesh::MeshHandle mesh_handle = sqmesh::invalid_handle;
      throw_on_status(
        mesh::create_surface_mesh(
          model_handle,
          algorithm_name,
          options,
          mesh_handle,
          context_handle
        )
      );
      return mesh_handle;
    },
    py::arg("model_handle"),
    py::arg("algorithm_name"),
    py::arg("options"),
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  mesh_module.def(
    "create_surface_mesh",
    [](geo::ModelHandle model_handle,
       const std::string &algorithm_name,
       const mesh::ParameterDictionary &parameters,
       base::ContextHandle context_handle) {
      mesh::MeshHandle mesh_handle = sqmesh::invalid_handle;
      throw_on_status(
        mesh::create_surface_mesh(
          model_handle,
          algorithm_name,
          parameters,
          mesh_handle,
          context_handle
        )
      );
      return mesh_handle;
    },
    py::arg("model_handle"),
    py::arg("algorithm_name"),
    py::arg("parameters"),
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  mesh_module.def(
    "create_volume_mesh",
    [](geo::ModelHandle model_handle,
       const std::string &algorithm_name,
       const mesh::MeshingOptions &options,
       base::ContextHandle context_handle) {
      mesh::MeshHandle mesh_handle = sqmesh::invalid_handle;
      throw_on_status(
        mesh::create_volume_mesh(
          model_handle,
          algorithm_name,
          options,
          mesh_handle,
          context_handle
        )
      );
      return mesh_handle;
    },
    py::arg("model_handle"),
    py::arg("algorithm_name"),
    py::arg("options"),
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  mesh_module.def(
    "create_volume_mesh",
    [](geo::ModelHandle model_handle,
       const std::string &algorithm_name,
       const mesh::ParameterDictionary &parameters,
       base::ContextHandle context_handle) {
      mesh::MeshHandle mesh_handle = sqmesh::invalid_handle;
      throw_on_status(
        mesh::create_volume_mesh(
          model_handle,
          algorithm_name,
          parameters,
          mesh_handle,
          context_handle
        )
      );
      return mesh_handle;
    },
    py::arg("model_handle"),
    py::arg("algorithm_name"),
    py::arg("parameters"),
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  mesh_module.def(
    "mesh_summary",
    [](mesh::MeshHandle mesh_handle, base::ContextHandle context_handle) {
      mesh::MeshSummary summary;
      throw_on_status(mesh::mesh_summary(mesh_handle, summary, context_handle));
      return summary;
    },
    py::arg("mesh_handle"),
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  mesh_module.def(
    "mesh_quality_report",
    [](mesh::MeshHandle mesh_handle, base::ContextHandle context_handle) {
      mesh::MeshQualityReport report;
      throw_on_status(
        mesh::mesh_quality_report(mesh_handle, report, context_handle)
      );
      return report;
    },
    py::arg("mesh_handle"),
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  mesh_module.def(
    "domain_snapshot",
    [](mesh::MeshHandle mesh_handle, base::ContextHandle context_handle) {
      mesh::Domain domain;
      throw_on_status(mesh::domain_snapshot(mesh_handle, domain, context_handle));
      return domain;
    },
    py::arg("mesh_handle"),
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  mesh_module.def(
    "nodes_count",
    [](mesh::MeshHandle mesh_handle, base::ContextHandle context_handle) {
      std::size_t count = 0U;
      throw_on_status(mesh::nodes_count(mesh_handle, count, context_handle));
      return count;
    },
    py::arg("mesh_handle"),
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  mesh_module.def(
    "cells_count",
    [](mesh::MeshHandle mesh_handle, base::ContextHandle context_handle) {
      std::size_t count = 0U;
      throw_on_status(mesh::cells_count(mesh_handle, count, context_handle));
      return count;
    },
    py::arg("mesh_handle"),
    py::arg("context_handle") = sqmesh::invalid_handle
  );
}
