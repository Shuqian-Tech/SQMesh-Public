// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "sqmesh/sqmesh.hpp"

#include "../src/mesh/auto_cfd/auto_cfd_surface_pipeline.hpp"
#include "../src/mesh/sizing/mesh_size_controls.hpp"
#include "../src/mesh/framework/meshing_framework.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>

namespace {

bool expect(bool condition, const char *message)
{
  if(condition) {
    return true;
  }

  std::fprintf(stderr, "bottom_up_pipeline_smoke: %s\n", message);
  return false;
}

bool write_patch_stl(const std::filesystem::path &path)
{
  std::ofstream stream(path);
  if(!stream) {
    return false;
  }

  stream <<
    "solid square_patch\n"
    "  facet normal 0 0 1\n"
    "    outer loop\n"
    "      vertex 0 0 0\n"
    "      vertex 1 0 0\n"
    "      vertex 1 1 0\n"
    "    endloop\n"
    "  endfacet\n"
    "  facet normal 0 0 1\n"
    "    outer loop\n"
    "      vertex 0 0 0\n"
    "      vertex 1 1 0\n"
    "      vertex 0 1 0\n"
    "    endloop\n"
    "  endfacet\n"
    "endsolid square_patch\n";

  return static_cast<bool>(stream);
}

class PipelineProbe final : public sqmesh::mesh::detail::BottomUpMeshingPipeline
{
public:
  [[nodiscard]] std::string_view name() const noexcept override
  {
    return "BottomUp Pipeline Probe";
  }

  [[nodiscard]] const std::string &stage_order() const noexcept
  {
    return stage_order_;
  }

  [[nodiscard]] bool reused_shared_polyline() const noexcept
  {
    return reused_shared_polyline_;
  }

  [[nodiscard]] bool saw_local_edge_size() const noexcept
  {
    return saw_local_edge_size_;
  }

protected:
  void reset() noexcept override
  {
    BottomUpMeshingPipeline::reset();
    target_size_ = 0.0;
    stage_order_.clear();
    reused_shared_polyline_ = false;
    saw_local_edge_size_ = false;
  }

  [[nodiscard]] sqmesh::base::StatusCode on_configure(
    const sqmesh::mesh::ParameterDictionary &parameters
  ) override
  {
    double minimum_length = 0.0;
    double target_size = 0.0;
    const bool has_minimum_length =
      parameters.try_get_number("minimum_length", minimum_length);
    const bool has_target_size = parameters.try_get_number("target_size", target_size);

    if(has_target_size) {
      if(target_size <= 0.0) {
        return sqmesh::base::StatusCode::invalid_argument;
      }
      target_size_ = target_size;
      return sqmesh::base::StatusCode::ok;
    }

    if(has_minimum_length) {
      if(minimum_length <= 0.0) {
        return sqmesh::base::StatusCode::invalid_argument;
      }
      target_size_ = minimum_length;
      return sqmesh::base::StatusCode::ok;
    }

    return sqmesh::base::StatusCode::invalid_argument;
  }

  [[nodiscard]] sqmesh::mesh::detail::MeshingDimension pipeline_dimension(
  ) const noexcept override
  {
    return sqmesh::mesh::detail::MeshingDimension::surface;
  }

  [[nodiscard]] double pipeline_global_target_size() const noexcept override
  {
    return target_size_;
  }

  [[nodiscard]] sqmesh::base::StatusCode on_prepare_pipeline() override
  {
    if(model_view().faces.empty() || model_view().edges.empty()) {
      return sqmesh::base::StatusCode::unsupported;
    }

    const auto &first_edge = model_view().edges.front();

    const sqmesh::mesh::detail::SampledEdgePolyline *first_polyline = nullptr;
    auto status = ensure_edge_polyline(first_edge, first_polyline);
    if(status != sqmesh::base::StatusCode::ok) {
      return status;
    }

    const sqmesh::mesh::detail::SampledEdgePolyline *second_polyline = nullptr;
    status = ensure_edge_polyline(first_edge.entity, second_polyline);
    if(status != sqmesh::base::StatusCode::ok) {
      return status;
    }

    reused_shared_polyline_ = first_polyline != nullptr && first_polyline == second_polyline;
    saw_local_edge_size_ =
      std::abs(
        sqmesh::mesh::detail::effective_edge_target_size(size_controls(), first_edge) -
        0.2
      ) <= 1.0e-12;

    if(!reused_shared_polyline_ || first_polyline == nullptr ||
       first_polyline->positions.size() < 2U ||
       first_polyline->curve.approximate_length <= 0.0) {
      return sqmesh::base::StatusCode::internal_error;
    }

    return sqmesh::base::StatusCode::ok;
  }

  [[nodiscard]] sqmesh::base::StatusCode on_stage_vertices(
    sqmesh::mesh::Domain &output
  ) override
  {
    static_cast<void>(output);
    stage_order_.push_back('0');
    return sqmesh::base::StatusCode::ok;
  }

  [[nodiscard]] sqmesh::base::StatusCode on_stage_edges(
    sqmesh::mesh::Domain &output
  ) override
  {
    static_cast<void>(output);
    stage_order_.push_back('1');
    return sqmesh::base::StatusCode::ok;
  }

  [[nodiscard]] sqmesh::base::StatusCode on_stage_faces(
    sqmesh::mesh::Domain &output
  ) override
  {
    const auto node_entity_group = output.create_entity_group(
      {
        sqmesh::mesh::EntityOrder::node,
        "pipeline_nodes",
        sqmesh::mesh::invalid_index,
        false,
        sqmesh::mesh::EntityKind::node_point,
      }
    );
    static_cast<void>(output.add_node(node_entity_group, {0.0, 0.0, 0.0}));
    stage_order_.push_back('2');
    return sqmesh::base::StatusCode::ok;
  }

  [[nodiscard]] sqmesh::base::StatusCode on_stage_regions(
    sqmesh::mesh::Domain &output
  ) override
  {
    static_cast<void>(output);
    stage_order_.push_back('3');
    return sqmesh::base::StatusCode::ok;
  }

private:
  double target_size_ = 0.0;
  std::string stage_order_ {};
  bool reused_shared_polyline_ = false;
  bool saw_local_edge_size_ = false;
};

} // namespace

int main()
{
  sqmesh::mesh::detail::register_builtin_meshing_algorithms();

  auto auto_cfd_algorithm =
    sqmesh::mesh::detail::create_meshing_algorithm("Auto CFD Surface Mesher");
  if(!expect(
       auto_cfd_algorithm != nullptr,
       "the algorithm registry should create the Auto CFD Surface Mesher implementation"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       dynamic_cast<sqmesh::mesh::detail::BottomUpMeshingPipeline *>(
         auto_cfd_algorithm.get()
       ) != nullptr,
       "Auto CFD Surface Mesher should build on the optional bottom-up pipeline base"
     )) {
    return EXIT_FAILURE;
  }
  if(!expect(
       sqmesh::mesh::detail::create_meshing_algorithm("auto_cfd_surface_mesher") != nullptr,
       "the Auto CFD surface mesher alias should resolve through the registry"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::base::ContextHandle context = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::base::initialize(context) == sqmesh::base::StatusCode::ok,
       "initialize should create a runtime context"
     )) {
    return EXIT_FAILURE;
  }

  const auto stl_path = std::filesystem::path("sqmesh_bottom_up_pipeline_ascii.stl");
  std::filesystem::remove(stl_path);
  if(!expect(
       write_patch_stl(stl_path),
       "the bottom-up pipeline smoke test should be able to write its STL input"
     )) {
    return EXIT_FAILURE;
  }

  sqmesh::geo::ModelHandle model_handle = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::geo::import_stl(stl_path.string(), model_handle, {}, context) ==
         sqmesh::base::StatusCode::ok,
       "import_stl should create a geometry model for the bottom-up pipeline smoke test"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }

  sqmesh::geo::ModelView geometry_view;
  if(!expect(
       sqmesh::geo::model_view(model_handle, geometry_view, context) ==
         sqmesh::base::StatusCode::ok,
       "model_view should expose a usable geometry view for the bottom-up pipeline smoke test"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }

  if(!expect(
       !geometry_view.edges.empty(),
       "the STL patch should expose geometry edges for shared polyline sampling"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }

  sqmesh::mesh::ParameterDictionary parameters;
  parameters.set_number("target_size", 1.0);

  sqmesh::mesh::detail::MeshingRequest request;
  request.context_handle = context;
  request.session_handle = sqmesh::base::current_session(context);
  request.model_handle = model_handle;
  request.target_dimension = sqmesh::mesh::detail::MeshingDimension::surface;
  request.size_controls.topology_revision = geometry_view.snapshot.topology_revision;
  request.size_controls.add_local_size(geometry_view.edges.front().entity, 0.2);

  sqmesh::mesh::ParameterDictionary auto_cfd_parameters;
  auto_cfd_parameters.set_number("minimum_length", 0.1);
  auto_cfd_parameters.set_number("maximum_length", 0.5);
  auto_cfd_parameters.set_number("distortion_angle", 18.0);
  auto_cfd_parameters.set_number("growth_rate", 1.2);
  auto_cfd_parameters.set_boolean("proximity", true);

  sqmesh::mesh::detail::AutoCfdSurfaceParameters resolved_auto_cfd_parameters;
  if(!expect(
       sqmesh::mesh::detail::resolve_auto_cfd_surface_parameters(
         auto_cfd_parameters,
         resolved_auto_cfd_parameters
       ) == sqmesh::base::StatusCode::ok,
       "the Auto CFD foundation should accept the Part 1 sizing parameter contract"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }

  sqmesh::mesh::detail::ResolvedMeshSizeControls auto_cfd_size_controls;
  if(!expect(
       sqmesh::mesh::detail::resolve_mesh_size_controls(
         geometry_view,
         {},
         resolved_auto_cfd_parameters.maximum_length,
         auto_cfd_size_controls
       ) == sqmesh::base::StatusCode::ok,
       "the Auto CFD foundation should resolve mesh size controls against the shared geometry view"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }

  sqmesh::mesh::detail::AutoCfdSurfacePipelineState auto_cfd_state;
  if(!expect(
       sqmesh::mesh::detail::build_auto_cfd_surface_pipeline_state(
         geometry_view,
         auto_cfd_size_controls,
         resolved_auto_cfd_parameters,
         auto_cfd_state
       ) == sqmesh::base::StatusCode::ok,
       "the Auto CFD foundation should build durable anchor, curve, and face pipeline state"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }
  if(!expect(
       auto_cfd_state.topology_revision == geometry_view.snapshot.topology_revision,
       "the Auto CFD pipeline state should retain the source topology revision"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }
  if(!expect(
       auto_cfd_state.anchor_candidates.size() == geometry_view.vertices.size(),
       "the Auto CFD foundation should snapshot every topology vertex as an anchor candidate"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }
  if(!expect(
       auto_cfd_state.curve_work_items.size() == geometry_view.edges.size(),
       "the Auto CFD foundation should create curve work items for every geometry edge"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }
  if(!expect(
       auto_cfd_state.face_work_items.size() == geometry_view.faces.size(),
       "the Auto CFD foundation should create face work items for every geometry face"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }
  if(!expect(
       auto_cfd_state.sizing_field.proximity_enabled &&
         auto_cfd_state.sizing_field.curvature_sources.empty() &&
         auto_cfd_state.sizing_field.proximity_sources.empty(),
       "the Auto CFD foundation should store the sizing-field contract without inventing curvature or proximity sources yet"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }
  if(!expect(
       !auto_cfd_state.curve_work_items.empty() &&
         auto_cfd_state.curve_work_items.front().curve.approximate_length > 0.0,
       "the Auto CFD foundation should retain edge curve info for later edge-stage work"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }
  if(!expect(
       !auto_cfd_state.face_work_items.empty() &&
         auto_cfd_state.face_work_items.front().boundary.face ==
           auto_cfd_state.face_work_items.front().face,
       "the Auto CFD face-stage scaffold should preserve the ordered boundary snapshot for each face"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }
  if(!expect(
       auto_cfd_state.boundary.vertex_node_indices.size() == geometry_view.vertices.size() &&
         auto_cfd_state.boundary.edge_discretizations.size() == geometry_view.edges.size() &&
         !auto_cfd_state.boundary.nodes.empty(),
       "the Auto CFD pipeline should now retain a durable 1D boundary-node dictionary for later edge and face stages"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }

  PipelineProbe probe;
  sqmesh::mesh::Domain output;
  if(!expect(
       probe.execute(request, parameters, output) == sqmesh::base::StatusCode::ok,
       "the optional bottom-up pipeline should execute its shared model-view, size-control, and stage-dispatch flow"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }

  if(!expect(
       probe.stage_order() == "012",
       "surface bottom-up execution should dispatch only the 0D, 1D, and 2D stages"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }

  if(!expect(
       probe.reused_shared_polyline(),
       "shared edge polyline sampling should reuse the same cached sample set on repeated requests"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }

  if(!expect(
       probe.saw_local_edge_size(),
       "bottom-up size-control resolution should preserve local edge target sizes before stage dispatch"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }

  sqmesh::mesh::MeshHandle auto_cfd_mesh = sqmesh::invalid_handle;
  if(!expect(
       sqmesh::mesh::create_surface_mesh(
         model_handle,
         "Auto CFD Surface Mesher",
         auto_cfd_parameters,
         auto_cfd_mesh,
         context
       ) == sqmesh::base::StatusCode::unsupported,
       "the Auto CFD path should stay explicit about unsupported discrete fixtures even after the bounded 2D face stage is implemented on the tested OCC subset"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }
  if(!expect(
       auto_cfd_mesh == sqmesh::invalid_handle,
       "unsupported Auto CFD fixtures should not allocate a mesh handle"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }

  request.target_dimension = sqmesh::mesh::detail::MeshingDimension::volume;
  if(!expect(
       probe.execute(request, parameters, output) == sqmesh::base::StatusCode::unsupported,
       "the optional bottom-up pipeline should keep dimension support explicit instead of acting as a universal base"
     )) {
    std::filesystem::remove(stl_path);
    return EXIT_FAILURE;
  }

  std::filesystem::remove(stl_path);
  if(!expect(
       sqmesh::base::shutdown_all() == sqmesh::base::StatusCode::ok,
       "shutdown_all should release the runtime after the bottom-up pipeline smoke test"
     )) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
