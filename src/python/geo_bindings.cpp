// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "bindings.hpp"
#include "exceptions.hpp"

#include "sqmesh/geo/api.hpp"

#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <cstdint>
#include <string>

namespace {

std::size_t hash_topology_entity_id(sqmesh::geo::TopologyEntityId entity) noexcept
{
  return (static_cast<std::size_t>(entity.index) << 8U) ^
         static_cast<std::size_t>(entity.dimension);
}

py::object cast_primary_outer_boundary_loop(
  const sqmesh::geo::FaceBoundaryLoops &boundary
)
{
  const auto *loop = sqmesh::geo::primary_outer_boundary_loop(boundary);
  if(loop == nullptr) {
    return py::none();
  }

  return py::cast(*loop);
}

} // namespace

void bind_geo(py::module_ &module)
{
  namespace geo = sqmesh::geo;
  namespace base = sqmesh::base;

  auto geo_module = module.def_submodule("geo", "SQMesh geometry bindings.");
  geo_module.attr("invalid_topology_index") = py::int_(geo::invalid_topology_index);
  geo_module.attr("invalid_boundary_loop_index") =
    py::int_(geo::invalid_boundary_loop_index);

  py::enum_<geo::IgesWriteMode>(geo_module, "IgesWriteMode")
    .value("faces", geo::IgesWriteMode::faces)
    .value("brep", geo::IgesWriteMode::brep);

  py::enum_<geo::TopologyDimension>(geo_module, "TopologyDimension")
    .value("vertex", geo::TopologyDimension::vertex)
    .value("edge", geo::TopologyDimension::edge)
    .value("face", geo::TopologyDimension::face)
    .value("region", geo::TopologyDimension::region);

  py::enum_<geo::FaceBoundaryLoopKind>(geo_module, "FaceBoundaryLoopKind")
    .value("unknown", geo::FaceBoundaryLoopKind::unknown)
    .value("outer", geo::FaceBoundaryLoopKind::outer)
    .value("inner", geo::FaceBoundaryLoopKind::inner);

  py::class_<geo::TopologyEntityId>(geo_module, "TopologyEntityId")
    .def(py::init<>())
    .def(
      py::init<geo::TopologyDimension, std::uint32_t>(),
      py::arg("dimension"),
      py::arg("index")
    )
    .def_readwrite("dimension", &geo::TopologyEntityId::dimension)
    .def_readwrite("index", &geo::TopologyEntityId::index)
    .def("is_valid", &geo::is_valid)
    .def("__bool__", &geo::is_valid)
    .def(
      "__repr__",
      [](geo::TopologyEntityId self) {
        return "<sqmesh.geo.TopologyEntityId dimension=" +
               std::to_string(static_cast<int>(self.dimension)) +
               " index=" + std::to_string(self.index) + ">";
      }
    )
    .def("__hash__", &hash_topology_entity_id)
    .def(py::self == py::self)
    .def(py::self != py::self);

  py::class_<geo::TopologyEntityInfo>(geo_module, "TopologyEntityInfo")
    .def(py::init<>())
    .def_readwrite("entity", &geo::TopologyEntityInfo::entity)
    .def_readwrite("parent_count", &geo::TopologyEntityInfo::parent_count)
    .def_readwrite("child_count", &geo::TopologyEntityInfo::child_count);

  py::class_<geo::TopologySnapshot>(geo_module, "TopologySnapshot")
    .def(py::init<>())
    .def_readwrite("topology_revision", &geo::TopologySnapshot::topology_revision)
    .def_readwrite("regions", &geo::TopologySnapshot::regions)
    .def_readwrite("faces", &geo::TopologySnapshot::faces)
    .def_readwrite("edges", &geo::TopologySnapshot::edges)
    .def_readwrite("vertices", &geo::TopologySnapshot::vertices)
    .def("entities", [](const geo::TopologySnapshot &self, geo::TopologyDimension dim) {
      return self.entities(dim);
    })
    .def("entity_count", &geo::TopologySnapshot::entity_count);

  py::class_<geo::FaceUvBounds>(geo_module, "FaceUvBounds")
    .def(py::init<>())
    .def_readwrite("face", &geo::FaceUvBounds::face)
    .def_readwrite("u_min", &geo::FaceUvBounds::u_min)
    .def_readwrite("u_max", &geo::FaceUvBounds::u_max)
    .def_readwrite("v_min", &geo::FaceUvBounds::v_min)
    .def_readwrite("v_max", &geo::FaceUvBounds::v_max);

  py::class_<geo::FaceSample>(geo_module, "FaceSample")
    .def(py::init<>())
    .def_readwrite("face", &geo::FaceSample::face)
    .def_readwrite("u", &geo::FaceSample::u)
    .def_readwrite("v", &geo::FaceSample::v)
    .def_readwrite("position", &geo::FaceSample::position)
    .def_readwrite("normal", &geo::FaceSample::normal)
    .def_readwrite("normal_defined", &geo::FaceSample::normal_defined);

  py::class_<geo::FaceCurvatureSample>(geo_module, "FaceCurvatureSample")
    .def(py::init<>())
    .def_readwrite("face", &geo::FaceCurvatureSample::face)
    .def_readwrite("u", &geo::FaceCurvatureSample::u)
    .def_readwrite("v", &geo::FaceCurvatureSample::v)
    .def_readwrite("min_curvature", &geo::FaceCurvatureSample::min_curvature)
    .def_readwrite("max_curvature", &geo::FaceCurvatureSample::max_curvature)
    .def_readwrite("mean_curvature", &geo::FaceCurvatureSample::mean_curvature)
    .def_readwrite(
      "gaussian_curvature",
      &geo::FaceCurvatureSample::gaussian_curvature
    )
    .def_readwrite(
      "curvature_defined",
      &geo::FaceCurvatureSample::curvature_defined
    );

  py::class_<geo::FaceDerivatives>(geo_module, "FaceDerivatives")
    .def(py::init<>())
    .def_readwrite("face", &geo::FaceDerivatives::face)
    .def_readwrite("u", &geo::FaceDerivatives::u)
    .def_readwrite("v", &geo::FaceDerivatives::v)
    .def_readwrite("position", &geo::FaceDerivatives::position)
    .def_readwrite("du", &geo::FaceDerivatives::du)
    .def_readwrite("dv", &geo::FaceDerivatives::dv)
    .def_readwrite("duu", &geo::FaceDerivatives::duu)
    .def_readwrite("duv", &geo::FaceDerivatives::duv)
    .def_readwrite("dvv", &geo::FaceDerivatives::dvv)
    .def_readwrite("normal", &geo::FaceDerivatives::normal)
    .def_readwrite(
      "first_derivatives_defined",
      &geo::FaceDerivatives::first_derivatives_defined
    )
    .def_readwrite(
      "second_derivatives_defined",
      &geo::FaceDerivatives::second_derivatives_defined
    )
    .def_readwrite("normal_defined", &geo::FaceDerivatives::normal_defined);

  py::class_<geo::FaceProjection>(geo_module, "FaceProjection")
    .def(py::init<>())
    .def_readwrite("face", &geo::FaceProjection::face)
    .def_readwrite("input_point", &geo::FaceProjection::input_point)
    .def_readwrite("projected_point", &geo::FaceProjection::projected_point)
    .def_readwrite("u", &geo::FaceProjection::u)
    .def_readwrite("v", &geo::FaceProjection::v)
    .def_readwrite("distance", &geo::FaceProjection::distance)
    .def_readwrite("normal", &geo::FaceProjection::normal)
    .def_readwrite("normal_defined", &geo::FaceProjection::normal_defined);

  py::class_<geo::FaceUvMapping>(geo_module, "FaceUvMapping")
    .def(py::init<>())
    .def_readwrite("face", &geo::FaceUvMapping::face)
    .def_readwrite("input_point", &geo::FaceUvMapping::input_point)
    .def_readwrite("mapped_point", &geo::FaceUvMapping::mapped_point)
    .def_readwrite("u", &geo::FaceUvMapping::u)
    .def_readwrite("v", &geo::FaceUvMapping::v)
    .def_readwrite("distance", &geo::FaceUvMapping::distance);

  py::class_<geo::EdgeCurveInfo>(geo_module, "EdgeCurveInfo")
    .def(py::init<>())
    .def_readwrite("edge", &geo::EdgeCurveInfo::edge)
    .def_readwrite("parameter_min", &geo::EdgeCurveInfo::parameter_min)
    .def_readwrite("parameter_max", &geo::EdgeCurveInfo::parameter_max)
    .def_readwrite("start_point", &geo::EdgeCurveInfo::start_point)
    .def_readwrite("end_point", &geo::EdgeCurveInfo::end_point)
    .def_readwrite(
      "approximate_length",
      &geo::EdgeCurveInfo::approximate_length
    );

  py::class_<geo::EdgeTangentSample>(geo_module, "EdgeTangentSample")
    .def(py::init<>())
    .def_readwrite("edge", &geo::EdgeTangentSample::edge)
    .def_readwrite("parameter", &geo::EdgeTangentSample::parameter)
    .def_readwrite("position", &geo::EdgeTangentSample::position)
    .def_readwrite("derivative", &geo::EdgeTangentSample::derivative)
    .def_readwrite("tangent", &geo::EdgeTangentSample::tangent)
    .def_readwrite("speed", &geo::EdgeTangentSample::speed)
    .def_readwrite("tangent_defined", &geo::EdgeTangentSample::tangent_defined);

  py::class_<geo::EdgeCurveSamplingOptions>(
    geo_module,
    "EdgeCurveSamplingOptions"
  )
    .def(py::init<>())
    .def_readwrite(
      "target_segment_length",
      &geo::EdgeCurveSamplingOptions::target_segment_length
    )
    .def_readwrite(
      "min_segment_count",
      &geo::EdgeCurveSamplingOptions::min_segment_count
    );

  py::class_<geo::EdgeCurveSamples>(geo_module, "EdgeCurveSamples")
    .def(py::init<>())
    .def_readwrite("curve", &geo::EdgeCurveSamples::curve)
    .def_readwrite("samples", &geo::EdgeCurveSamples::samples);

  py::class_<geo::FaceBoundaryEdgeUse>(geo_module, "FaceBoundaryEdgeUse")
    .def(py::init<>())
    .def_readwrite("edge", &geo::FaceBoundaryEdgeUse::edge)
    .def_readwrite("start_vertex", &geo::FaceBoundaryEdgeUse::start_vertex)
    .def_readwrite("end_vertex", &geo::FaceBoundaryEdgeUse::end_vertex)
    .def_readwrite(
      "same_orientation_as_edge",
      &geo::FaceBoundaryEdgeUse::same_orientation_as_edge
    )
    .def_readwrite("is_seam", &geo::FaceBoundaryEdgeUse::is_seam)
    .def_readwrite("is_degenerate", &geo::FaceBoundaryEdgeUse::is_degenerate);

  py::class_<geo::FaceBoundaryLoop>(geo_module, "FaceBoundaryLoop")
    .def(py::init<>())
    .def_readwrite("kind", &geo::FaceBoundaryLoop::kind)
    .def_readwrite("closed", &geo::FaceBoundaryLoop::closed)
    .def_readwrite("continuous", &geo::FaceBoundaryLoop::continuous)
    .def_readwrite(
      "seam_edge_use_count",
      &geo::FaceBoundaryLoop::seam_edge_use_count
    )
    .def_readwrite(
      "degenerate_edge_use_count",
      &geo::FaceBoundaryLoop::degenerate_edge_use_count
    )
    .def_readwrite(
      "repeated_edge_use_count",
      &geo::FaceBoundaryLoop::repeated_edge_use_count
    )
    .def_readwrite("edge_uses", &geo::FaceBoundaryLoop::edge_uses)
    .def_readwrite("vertex_ids", &geo::FaceBoundaryLoop::vertex_ids);

  py::class_<geo::FaceBoundaryLoops>(geo_module, "FaceBoundaryLoops")
    .def(py::init<>())
    .def_readwrite("face", &geo::FaceBoundaryLoops::face)
    .def_readwrite("loops", &geo::FaceBoundaryLoops::loops)
    .def_readwrite(
      "primary_outer_loop_index",
      &geo::FaceBoundaryLoops::primary_outer_loop_index
    )
    .def_readwrite("outer_loop_count", &geo::FaceBoundaryLoops::outer_loop_count)
    .def_readwrite("inner_loop_count", &geo::FaceBoundaryLoops::inner_loop_count)
    .def_readwrite(
      "unknown_loop_count",
      &geo::FaceBoundaryLoops::unknown_loop_count
    )
    .def_readwrite(
      "closed_loop_count",
      &geo::FaceBoundaryLoops::closed_loop_count
    )
    .def_readwrite("open_loop_count", &geo::FaceBoundaryLoops::open_loop_count)
    .def_readwrite(
      "non_continuous_loop_count",
      &geo::FaceBoundaryLoops::non_continuous_loop_count
    )
    .def_readwrite("seam_loop_count", &geo::FaceBoundaryLoops::seam_loop_count)
    .def_readwrite(
      "seam_edge_use_count",
      &geo::FaceBoundaryLoops::seam_edge_use_count
    )
    .def_readwrite(
      "degenerate_loop_count",
      &geo::FaceBoundaryLoops::degenerate_loop_count
    )
    .def_readwrite(
      "degenerate_edge_use_count",
      &geo::FaceBoundaryLoops::degenerate_edge_use_count
    )
    .def_readwrite(
      "repeated_edge_loop_count",
      &geo::FaceBoundaryLoops::repeated_edge_loop_count
    )
    .def_readwrite(
      "repeated_edge_use_count",
      &geo::FaceBoundaryLoops::repeated_edge_use_count
    )
    .def_readwrite("has_holes", &geo::FaceBoundaryLoops::has_holes)
    .def_readwrite("has_seams", &geo::FaceBoundaryLoops::has_seams)
    .def("primary_outer_loop", &cast_primary_outer_boundary_loop);

  py::class_<geo::VertexView>(geo_module, "VertexView")
    .def(py::init<>())
    .def_readwrite("model_handle", &geo::VertexView::model_handle)
    .def_readwrite("context_handle", &geo::VertexView::context_handle)
    .def_readwrite("entity", &geo::VertexView::entity)
    .def_readwrite("edge_ids", &geo::VertexView::edge_ids);

  py::class_<geo::EdgeView>(geo_module, "EdgeView")
    .def(py::init<>())
    .def_readwrite("model_handle", &geo::EdgeView::model_handle)
    .def_readwrite("context_handle", &geo::EdgeView::context_handle)
    .def_readwrite("entity", &geo::EdgeView::entity)
    .def_readwrite("face_ids", &geo::EdgeView::face_ids)
    .def_readwrite("vertex_ids", &geo::EdgeView::vertex_ids)
    .def("curve_info", [](const geo::EdgeView &self) {
      geo::EdgeCurveInfo info;
      throw_on_status(geo::edge_curve_info(self, info));
      return info;
    })
    .def(
      "sample_tangent",
      [](const geo::EdgeView &self, double parameter) {
        geo::EdgeTangentSample sample;
        throw_on_status(geo::sample_edge_tangent(self, parameter, sample));
        return sample;
      },
      py::arg("parameter")
    )
    .def(
      "sample_curve",
      [](const geo::EdgeView &self, const geo::EdgeCurveSamplingOptions &options) {
        geo::EdgeCurveSamples samples;
        throw_on_status(geo::sample_edge_curve(self, options, samples));
        return samples;
      },
      py::arg("options") = geo::EdgeCurveSamplingOptions {}
    );

  py::class_<geo::FaceView>(geo_module, "FaceView")
    .def(py::init<>())
    .def_readwrite("model_handle", &geo::FaceView::model_handle)
    .def_readwrite("context_handle", &geo::FaceView::context_handle)
    .def_readwrite("entity", &geo::FaceView::entity)
    .def_readwrite("region_ids", &geo::FaceView::region_ids)
    .def_readwrite("edge_ids", &geo::FaceView::edge_ids)
    .def_readwrite("ordered_boundary", &geo::FaceView::ordered_boundary)
    .def("boundary_loops", [](const geo::FaceView &self) {
      geo::FaceBoundaryLoops boundary;
      throw_on_status(geo::face_boundary_loops(self, boundary));
      return boundary;
    })
    .def("uv_bounds", [](const geo::FaceView &self) {
      geo::FaceUvBounds bounds;
      throw_on_status(geo::face_uv_bounds(self, bounds));
      return bounds;
    })
    .def(
      "sample",
      [](const geo::FaceView &self, double u, double v) {
        geo::FaceSample sample;
        throw_on_status(geo::sample_face(self, u, v, sample));
        return sample;
      },
      py::arg("u"),
      py::arg("v")
    )
    .def(
      "sample_curvature",
      [](const geo::FaceView &self, double u, double v) {
        geo::FaceCurvatureSample sample;
        throw_on_status(geo::sample_face_curvature(self, u, v, sample));
        return sample;
      },
      py::arg("u"),
      py::arg("v")
    )
    .def(
      "sample_derivatives",
      [](const geo::FaceView &self, double u, double v) {
        geo::FaceDerivatives derivatives;
        throw_on_status(geo::sample_face_derivatives(self, u, v, derivatives));
        return derivatives;
      },
      py::arg("u"),
      py::arg("v")
    )
    .def(
      "project_point",
      [](const geo::FaceView &self, const geo::Point3 &point) {
        geo::FaceProjection projection;
        throw_on_status(geo::project_point_to_face(self, point, projection));
        return projection;
      },
      py::arg("point")
    )
    .def(
      "recover_uv",
      [](const geo::FaceView &self, const geo::Point3 &point) {
        geo::FaceUvMapping mapping;
        throw_on_status(geo::recover_face_uv(self, point, mapping));
        return mapping;
      },
      py::arg("point")
    );

  py::class_<geo::RegionView>(geo_module, "RegionView")
    .def(py::init<>())
    .def_readwrite("model_handle", &geo::RegionView::model_handle)
    .def_readwrite("context_handle", &geo::RegionView::context_handle)
    .def_readwrite("entity", &geo::RegionView::entity)
    .def_readwrite("face_ids", &geo::RegionView::face_ids);

  py::class_<geo::ModelView>(geo_module, "ModelView")
    .def(py::init<>())
    .def_readwrite("model_handle", &geo::ModelView::model_handle)
    .def_readwrite("context_handle", &geo::ModelView::context_handle)
    .def_readwrite("snapshot", &geo::ModelView::snapshot)
    .def_readwrite("regions", &geo::ModelView::regions)
    .def_readwrite("faces", &geo::ModelView::faces)
    .def_readwrite("edges", &geo::ModelView::edges)
    .def_readwrite("vertices", &geo::ModelView::vertices)
    .def_readwrite("root_face_ids", &geo::ModelView::root_face_ids)
    .def_readwrite("root_edge_ids", &geo::ModelView::root_edge_ids)
    .def_readwrite("root_vertex_ids", &geo::ModelView::root_vertex_ids)
    .def("entity_count", &geo::ModelView::entity_count)
    .def(
      "find_region",
      [](const geo::ModelView &self, geo::TopologyEntityId entity) -> py::object {
        const auto *result = self.find_region(entity);
        return result == nullptr ? py::none() : py::cast(*result);
      }
    )
    .def(
      "find_face",
      [](const geo::ModelView &self, geo::TopologyEntityId entity) -> py::object {
        const auto *result = self.find_face(entity);
        return result == nullptr ? py::none() : py::cast(*result);
      }
    )
    .def(
      "find_edge",
      [](const geo::ModelView &self, geo::TopologyEntityId entity) -> py::object {
        const auto *result = self.find_edge(entity);
        return result == nullptr ? py::none() : py::cast(*result);
      }
    )
    .def(
      "find_vertex",
      [](const geo::ModelView &self, geo::TopologyEntityId entity) -> py::object {
        const auto *result = self.find_vertex(entity);
        return result == nullptr ? py::none() : py::cast(*result);
      }
    );

  py::class_<geo::FeatureEdgeOptions>(geo_module, "FeatureEdgeOptions")
    .def(py::init<>())
    .def_readwrite(
      "feature_angle_degrees",
      &geo::FeatureEdgeOptions::feature_angle_degrees
    )
    .def_readwrite(
      "include_boundary_edges",
      &geo::FeatureEdgeOptions::include_boundary_edges
    )
    .def_readwrite(
      "include_non_manifold_edges",
      &geo::FeatureEdgeOptions::include_non_manifold_edges
    );

  py::class_<geo::FeatureEdgeReport>(geo_module, "FeatureEdgeReport")
    .def(py::init<>())
    .def_readwrite("edges", &geo::FeatureEdgeReport::edges)
    .def_readwrite("sharp_edge_count", &geo::FeatureEdgeReport::sharp_edge_count)
    .def_readwrite(
      "boundary_edge_count",
      &geo::FeatureEdgeReport::boundary_edge_count
    )
    .def_readwrite(
      "non_manifold_edge_count",
      &geo::FeatureEdgeReport::non_manifold_edge_count
    );

  py::class_<geo::ModelSummary>(geo_module, "ModelSummary")
    .def(py::init<>())
    .def_readwrite("entity_count", &geo::ModelSummary::entity_count)
    .def_readwrite("compound_count", &geo::ModelSummary::compound_count)
    .def_readwrite("compsolid_count", &geo::ModelSummary::compsolid_count)
    .def_readwrite("solid_count", &geo::ModelSummary::solid_count)
    .def_readwrite("shell_count", &geo::ModelSummary::shell_count)
    .def_readwrite("face_count", &geo::ModelSummary::face_count)
    .def_readwrite("wire_count", &geo::ModelSummary::wire_count)
    .def_readwrite("edge_count", &geo::ModelSummary::edge_count)
    .def_readwrite("vertex_count", &geo::ModelSummary::vertex_count);

  py::class_<geo::StepImportOptions>(geo_module, "StepImportOptions")
    .def(py::init<>())
    .def_readwrite(
      "system_length_unit",
      &geo::StepImportOptions::system_length_unit
    );

  py::class_<geo::IgesImportOptions>(geo_module, "IgesImportOptions")
    .def(py::init<>())
    .def_readwrite(
      "read_visible_only",
      &geo::IgesImportOptions::read_visible_only
    );

  py::class_<geo::StlImportOptions>(geo_module, "StlImportOptions")
    .def(py::init<>())
    .def_readwrite(
      "relative_merge_tolerance",
      &geo::StlImportOptions::relative_merge_tolerance
    );

  py::class_<geo::TopologyCheckReport>(geo_module, "TopologyCheckReport")
    .def(py::init<>())
    .def_readwrite("is_valid", &geo::TopologyCheckReport::is_valid)
    .def_readwrite(
      "free_edge_count",
      &geo::TopologyCheckReport::free_edge_count
    )
    .def_readwrite(
      "contiguous_edge_count",
      &geo::TopologyCheckReport::contiguous_edge_count
    )
    .def_readwrite(
      "multiple_edge_count",
      &geo::TopologyCheckReport::multiple_edge_count
    )
    .def_readwrite("model_summary", &geo::TopologyCheckReport::model_summary);

  py::class_<geo::TopoOptions>(geo_module, "TopoOptions")
    .def(py::init<>())
    .def_readwrite("tolerance", &geo::TopoOptions::tolerance)
    .def_readwrite("min_tolerance", &geo::TopoOptions::min_tolerance)
    .def_readwrite("max_tolerance", &geo::TopoOptions::max_tolerance)
    .def_readwrite("fix_degenerated", &geo::TopoOptions::fix_degenerated)
    .def_readwrite("fix_small_edges", &geo::TopoOptions::fix_small_edges)
    .def_readwrite("fix_small_faces", &geo::TopoOptions::fix_small_faces)
    .def_readwrite("sew_faces", &geo::TopoOptions::sew_faces)
    .def_readwrite("make_solids", &geo::TopoOptions::make_solids);

  py::class_<geo::TopoReport>(geo_module, "TopoReport")
    .def(py::init<>())
    .def_readwrite("before", &geo::TopoReport::before)
    .def_readwrite("after", &geo::TopoReport::after)
    .def_readwrite(
      "topology_revision_before",
      &geo::TopoReport::topology_revision_before
    )
    .def_readwrite(
      "topology_revision_after",
      &geo::TopoReport::topology_revision_after
    )
    .def_readwrite("modified", &geo::TopoReport::modified)
    .def_readwrite(
      "topology_identity_changed",
      &geo::TopoReport::topology_identity_changed
    )
    .def_readwrite(
      "free_edges_reduced",
      &geo::TopoReport::free_edges_reduced
    )
    .def_readwrite("sewing_performed", &geo::TopoReport::sewing_performed)
    .def_readwrite("sewing_modified", &geo::TopoReport::sewing_modified)
    .def_readwrite("shape_fix_performed", &geo::TopoReport::shape_fix_performed)
    .def_readwrite("shape_fix_modified", &geo::TopoReport::shape_fix_modified);

  geo_module.def("module_name", []() {
    return std::string(geo::module_name());
  });
  geo_module.def("cad_io_available", &geo::cad_io_available);
  geo_module.def(
    "create_placeholder_model",
    [](base::ContextHandle context_handle) {
      geo::ModelHandle model_handle = sqmesh::invalid_handle;
      throw_on_status(geo::create_placeholder_model(model_handle, context_handle));
      return model_handle;
    },
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  geo_module.def(
    "import_step",
    [](const std::string &path,
       const geo::StepImportOptions &options,
       base::ContextHandle context_handle) {
      geo::ModelHandle model_handle = sqmesh::invalid_handle;
      throw_on_status(geo::import_step(path, model_handle, options, context_handle));
      return model_handle;
    },
    py::arg("path"),
    py::arg("options") = geo::StepImportOptions {},
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  geo_module.def(
    "import_iges",
    [](const std::string &path,
       const geo::IgesImportOptions &options,
       base::ContextHandle context_handle) {
      geo::ModelHandle model_handle = sqmesh::invalid_handle;
      throw_on_status(geo::import_iges(path, model_handle, options, context_handle));
      return model_handle;
    },
    py::arg("path"),
    py::arg("options") = geo::IgesImportOptions {},
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  geo_module.def(
    "import_stl",
    [](const std::string &path,
       const geo::StlImportOptions &options,
       base::ContextHandle context_handle) {
      geo::ModelHandle model_handle = sqmesh::invalid_handle;
      throw_on_status(geo::import_stl(path, model_handle, options, context_handle));
      return model_handle;
    },
    py::arg("path"),
    py::arg("options") = geo::StlImportOptions {},
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  geo_module.def(
    "export_step",
    [](geo::ModelHandle model_handle,
       const std::string &path,
       double linear_tolerance,
       const std::string &unit_name,
       const std::string &schema,
       base::ContextHandle context_handle) {
      geo::StepExportOptions options;
      options.linear_tolerance = linear_tolerance;
      options.unit_name = unit_name;
      options.schema = schema;
      throw_on_status(geo::export_step(model_handle, path, options, context_handle));
    },
    py::arg("model_handle"),
    py::arg("path"),
    py::arg("linear_tolerance") = 0.0,
    py::arg("unit_name") = "",
    py::arg("schema") = "",
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  geo_module.def(
    "export_iges",
    [](geo::ModelHandle model_handle,
       const std::string &path,
       const std::string &unit_name,
       geo::IgesWriteMode write_mode,
       base::ContextHandle context_handle) {
      geo::IgesExportOptions options;
      options.unit_name = unit_name;
      options.write_mode = write_mode;
      throw_on_status(geo::export_iges(model_handle, path, options, context_handle));
    },
    py::arg("model_handle"),
    py::arg("path"),
    py::arg("unit_name") = "",
    py::arg("write_mode") = geo::IgesWriteMode::brep,
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  geo_module.def(
    "model_summary",
    [](geo::ModelHandle model_handle, base::ContextHandle context_handle) {
      geo::ModelSummary summary;
      throw_on_status(geo::model_summary(model_handle, summary, context_handle));
      return summary;
    },
    py::arg("model_handle"),
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  geo_module.def(
    "model_proxy_mesh",
    [](geo::ModelHandle model_handle, base::ContextHandle context_handle) {
      sqmesh::Handle mesh_handle = sqmesh::invalid_handle;
      throw_on_status(geo::model_proxy_mesh(model_handle, mesh_handle, context_handle));
      return mesh_handle;
    },
    py::arg("model_handle"),
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  geo_module.def(
    "topology_snapshot",
    [](geo::ModelHandle model_handle, base::ContextHandle context_handle) {
      geo::TopologySnapshot snapshot;
      throw_on_status(geo::topology_snapshot(model_handle, snapshot, context_handle));
      return snapshot;
    },
    py::arg("model_handle"),
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  geo_module.def(
    "model_view",
    [](geo::ModelHandle model_handle, base::ContextHandle context_handle) {
      geo::ModelView view;
      throw_on_status(geo::model_view(model_handle, view, context_handle));
      return view;
    },
    py::arg("model_handle"),
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  geo_module.def(
    "topology_children",
    [](geo::ModelHandle model_handle,
       geo::TopologyEntityId entity,
       base::ContextHandle context_handle) {
      std::vector<geo::TopologyEntityId> children;
      throw_on_status(geo::topology_children(model_handle, entity, children, context_handle));
      return children;
    },
    py::arg("model_handle"),
    py::arg("entity"),
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  geo_module.def(
    "topology_parents",
    [](geo::ModelHandle model_handle,
       geo::TopologyEntityId entity,
       base::ContextHandle context_handle) {
      std::vector<geo::TopologyEntityId> parents;
      throw_on_status(geo::topology_parents(model_handle, entity, parents, context_handle));
      return parents;
    },
    py::arg("model_handle"),
    py::arg("entity"),
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  geo_module.def(
    "face_uv_bounds",
    [](geo::ModelHandle model_handle,
       geo::TopologyEntityId face_entity,
       base::ContextHandle context_handle) {
      geo::FaceUvBounds bounds;
      throw_on_status(geo::face_uv_bounds(model_handle, face_entity, bounds, context_handle));
      return bounds;
    },
    py::arg("model_handle"),
    py::arg("face_entity"),
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  geo_module.def(
    "sample_face",
    [](geo::ModelHandle model_handle,
       geo::TopologyEntityId face_entity,
       double u,
       double v,
       base::ContextHandle context_handle) {
      geo::FaceSample sample;
      throw_on_status(geo::sample_face(model_handle, face_entity, u, v, sample, context_handle));
      return sample;
    },
    py::arg("model_handle"),
    py::arg("face_entity"),
    py::arg("u"),
    py::arg("v"),
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  geo_module.def(
    "sample_face_curvature",
    [](geo::ModelHandle model_handle,
       geo::TopologyEntityId face_entity,
       double u,
       double v,
       base::ContextHandle context_handle) {
      geo::FaceCurvatureSample sample;
      throw_on_status(
        geo::sample_face_curvature(model_handle, face_entity, u, v, sample, context_handle)
      );
      return sample;
    },
    py::arg("model_handle"),
    py::arg("face_entity"),
    py::arg("u"),
    py::arg("v"),
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  geo_module.def(
    "sample_face_derivatives",
    [](geo::ModelHandle model_handle,
       geo::TopologyEntityId face_entity,
       double u,
       double v,
       base::ContextHandle context_handle) {
      geo::FaceDerivatives derivatives;
      throw_on_status(
        geo::sample_face_derivatives(
          model_handle,
          face_entity,
          u,
          v,
          derivatives,
          context_handle
        )
      );
      return derivatives;
    },
    py::arg("model_handle"),
    py::arg("face_entity"),
    py::arg("u"),
    py::arg("v"),
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  geo_module.def(
    "project_point_to_face",
    [](geo::ModelHandle model_handle,
       geo::TopologyEntityId face_entity,
       const geo::Point3 &point,
       base::ContextHandle context_handle) {
      geo::FaceProjection projection;
      throw_on_status(
        geo::project_point_to_face(
          model_handle,
          face_entity,
          point,
          projection,
          context_handle
        )
      );
      return projection;
    },
    py::arg("model_handle"),
    py::arg("face_entity"),
    py::arg("point"),
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  geo_module.def(
    "recover_face_uv",
    [](geo::ModelHandle model_handle,
       geo::TopologyEntityId face_entity,
       const geo::Point3 &point,
       base::ContextHandle context_handle) {
      geo::FaceUvMapping mapping;
      throw_on_status(
        geo::recover_face_uv(model_handle, face_entity, point, mapping, context_handle)
      );
      return mapping;
    },
    py::arg("model_handle"),
    py::arg("face_entity"),
    py::arg("point"),
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  geo_module.def(
    "edge_curve_info",
    [](geo::ModelHandle model_handle,
       geo::TopologyEntityId edge_entity,
       base::ContextHandle context_handle) {
      geo::EdgeCurveInfo info;
      throw_on_status(geo::edge_curve_info(model_handle, edge_entity, info, context_handle));
      return info;
    },
    py::arg("model_handle"),
    py::arg("edge_entity"),
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  geo_module.def(
    "sample_edge_tangent",
    [](geo::ModelHandle model_handle,
       geo::TopologyEntityId edge_entity,
       double parameter,
       base::ContextHandle context_handle) {
      geo::EdgeTangentSample sample;
      throw_on_status(
        geo::sample_edge_tangent(
          model_handle,
          edge_entity,
          parameter,
          sample,
          context_handle
        )
      );
      return sample;
    },
    py::arg("model_handle"),
    py::arg("edge_entity"),
    py::arg("parameter"),
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  geo_module.def(
    "sample_edge_curve",
    [](geo::ModelHandle model_handle,
       geo::TopologyEntityId edge_entity,
       const geo::EdgeCurveSamplingOptions &options,
       base::ContextHandle context_handle) {
      geo::EdgeCurveSamples samples;
      throw_on_status(
        geo::sample_edge_curve(model_handle, edge_entity, options, samples, context_handle)
      );
      return samples;
    },
    py::arg("model_handle"),
    py::arg("edge_entity"),
    py::arg("options") = geo::EdgeCurveSamplingOptions {},
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  geo_module.def(
    "face_boundary_loops",
    [](geo::ModelHandle model_handle,
       geo::TopologyEntityId face_entity,
       base::ContextHandle context_handle) {
      geo::FaceBoundaryLoops boundary;
      throw_on_status(
        geo::face_boundary_loops(model_handle, face_entity, boundary, context_handle)
      );
      return boundary;
    },
    py::arg("model_handle"),
    py::arg("face_entity"),
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  geo_module.def(
    "feature_edges",
    [](geo::ModelHandle model_handle,
       const geo::FeatureEdgeOptions &options,
       base::ContextHandle context_handle) {
      geo::FeatureEdgeReport report;
      throw_on_status(geo::feature_edges(model_handle, report, options, context_handle));
      return report;
    },
    py::arg("model_handle"),
    py::arg("options") = geo::FeatureEdgeOptions {},
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  geo_module.def(
    "check_topology",
    [](geo::ModelHandle model_handle, base::ContextHandle context_handle) {
      geo::TopologyCheckReport report;
      throw_on_status(geo::check_topology(model_handle, report, context_handle));
      return report;
    },
    py::arg("model_handle"),
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  geo_module.def(
    "free_edge_count",
    [](geo::ModelHandle model_handle, base::ContextHandle context_handle) {
      std::size_t count = 0U;
      throw_on_status(geo::free_edge_count(model_handle, count, context_handle));
      return count;
    },
    py::arg("model_handle"),
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  geo_module.def(
    "topo",
    [](geo::ModelHandle model_handle,
       const geo::TopoOptions &options,
       base::ContextHandle context_handle) {
      geo::TopoReport report;
      throw_on_status(geo::topo(model_handle, report, options, context_handle));
      return report;
    },
    py::arg("model_handle"),
    py::arg("options") = geo::TopoOptions {},
    py::arg("context_handle") = sqmesh::invalid_handle
  );
  geo_module.def(
    "placeholder_model_summary",
    &geo::placeholder_model_summary
  );
  geo_module.def(
    "primary_outer_boundary_loop",
    &cast_primary_outer_boundary_loop,
    py::arg("boundary")
  );
}
