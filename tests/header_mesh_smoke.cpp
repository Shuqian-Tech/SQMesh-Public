// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "sqmesh/mesh/api.hpp"

#include <cstdlib>

int main()
{
  const sqmesh::mesh::MshImportOptions msh_import_options {};
  const sqmesh::mesh::MshExportOptions msh_export_options {};
  const sqmesh::mesh::ObjImportOptions obj_import_options {};
  const sqmesh::mesh::ObjExportOptions obj_export_options {};
  const sqmesh::mesh::CgnsImportOptions cgns_import_options {};
  const sqmesh::mesh::CgnsExportOptions cgns_export_options {};
  const sqmesh::mesh::MeshSizeControl mesh_size_control {};
  const sqmesh::mesh::MeshSizeControls mesh_size_controls {};
  const sqmesh::mesh::MeshingOptions meshing_options {};
  const sqmesh::mesh::AutoCfdSpacingOptions spacing_options {};
  const sqmesh::mesh::AutoCfdSurfaceDefaults surface_defaults {};
  sqmesh::mesh::MeshingOptions configured_options {};
  sqmesh::mesh::apply_auto_cfd_spacing(configured_options, spacing_options);
  sqmesh::mesh::apply_auto_cfd_surface_defaults(configured_options, surface_defaults);
  const auto domain = sqmesh::mesh::make_dummy_tetra_domain();
  const auto summary = domain.summary();
  const auto quality_report = domain.quality_report();
  return summary.node_count == 4U && summary.edge_count == 0U &&
           summary.cell_count == 1U &&
           quality_report.supported_element_count == 5U &&
           quality_report.inverted_element_count == 0U &&
           msh_import_options.read_physical_names &&
           msh_export_options.write_physical_names &&
           msh_export_options.format_version ==
             sqmesh::mesh::MshFormatVersion::gmsh22_ascii &&
           obj_import_options.read_groups &&
           obj_export_options.write_groups &&
           obj_export_options.write_object_name &&
           cgns_import_options.read_section_names &&
           cgns_export_options.base_name.empty() &&
           cgns_export_options.zone_name.empty() &&
           !sqmesh::geo::is_valid(mesh_size_control.entity) &&
           mesh_size_control.target_size == 0.0 &&
           mesh_size_controls.topology_revision == 0U &&
           mesh_size_controls.empty() &&
           spacing_options.growth_rate == 1.2 &&
           spacing_options.feature_angle == 15.0 &&
           surface_defaults.element_type == "tria" &&
           meshing_options.size_controls.empty() &&
           meshing_options.parameters.size() == 0U &&
           configured_options.parameters.contains("spacing_growth_rate") &&
           configured_options.parameters.contains("spacing_feature_angle") &&
           configured_options.parameters.contains("spacing_proximity") &&
           configured_options.parameters.contains("growth_rate") &&
           configured_options.parameters.contains("distortion_angle") &&
           configured_options.parameters.contains("proximity")
           ? EXIT_SUCCESS
           : EXIT_FAILURE;
}
