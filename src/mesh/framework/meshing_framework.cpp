// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "../auto_cfd/auto_cfd_surface_mesher.hpp"
#include "../boundary_layer/boundary_layer_mesher.hpp"
#include "../tet/tet_volume_mesher.hpp"
#include "meshing_framework.hpp"

#include "../sizing/mesh_size_controls.hpp"

#include "core/runtime_registry.hpp"

#include <cmath>
#include <exception>
#include <mutex>
#include <string>
#include <unordered_map>
#include <utility>

namespace sqmesh::mesh::detail {
namespace {

constexpr std::string_view kDummyMesherName = "Dummy Mesher";
constexpr std::string_view kDummyMesherAlias = "dummy_mesher";
constexpr std::string_view kAutoCfdSurfaceMesherName = "Auto CFD Surface Mesher";
constexpr std::string_view kAutoCfdSurfaceMesherAlias = "auto_cfd_surface_mesher";
constexpr std::string_view kTetVolumeMesherName = "Tetrahedral Volume Mesher";
constexpr std::string_view kTetVolumeMesherAlias = "tetrahedral_volume_mesher";
constexpr std::string_view kBLMesherName = "Boundary Layer Mesher";
constexpr std::string_view kBLMesherAlias = "boundary_layer_mesher";
constexpr double kMinimumSharedEdgeTargetSize = 1.0e-9;

[[nodiscard]] std::string_view meshing_dimension_name(MeshingDimension dimension) noexcept
{
  switch(dimension) {
  case MeshingDimension::surface:
    return "surface";
  case MeshingDimension::volume:
    return "volume";
  default:
    return "unknown";
  }
}

[[nodiscard]] base::StatusCode validate_meshing_request(
  const MeshingRequest &request,
  MeshingDimension expected_dimension,
  std::string_view mesher_name
)
{
  if(request.target_dimension != expected_dimension) {
    return core::detail::publish_error(
      base::StatusCode::unsupported,
      std::string(mesher_name) + " only supports " +
        std::string(meshing_dimension_name(expected_dimension)) + " meshing requests."
    );
  }

  if(request.context_handle == sqmesh::invalid_handle) {
    return core::detail::publish_error(
      base::StatusCode::invalid_handle,
      std::string(mesher_name) + " requires a valid context handle."
    );
  }

  if(request.session_handle == sqmesh::invalid_handle) {
    return core::detail::publish_error(
      base::StatusCode::invalid_handle,
      std::string(mesher_name) + " requires a valid session handle."
    );
  }

  if(request.model_handle == sqmesh::invalid_handle) {
    return core::detail::publish_error(
      base::StatusCode::invalid_handle,
      std::string(mesher_name) + " requires a valid model handle."
    );
  }

  return core::detail::clear_error_state();
}

class MeshingAlgorithmRegistry final
{
public:
  [[nodiscard]] bool register_factory(
    std::string_view algorithm_name,
    MeshingAlgorithmFactory factory
  )
  {
    return factories_.emplace(std::string(algorithm_name), factory).second;
  }

  [[nodiscard]] MeshingAlgorithmPtr create(std::string_view algorithm_name) const
  {
    const auto factory_it = factories_.find(std::string(algorithm_name));
    if(factory_it == factories_.end()) {
      return nullptr;
    }
    return factory_it->second();
  }

  [[nodiscard]] bool contains(std::string_view algorithm_name) const
  {
    return factories_.find(std::string(algorithm_name)) != factories_.end();
  }

private:
  std::unordered_map<std::string, MeshingAlgorithmFactory> factories_ {};
};

MeshingAlgorithmRegistry &registry() noexcept
{
  static MeshingAlgorithmRegistry instance;
  return instance;
}

std::mutex &registry_mutex() noexcept
{
  static std::mutex instance;
  return instance;
}

class DummyMesher final : public MeshingAlgorithm
{
public:
  [[nodiscard]] std::string_view name() const noexcept override
  {
    return kDummyMesherName;
  }

protected:
  void reset() noexcept override
  {
    request_ = {};
  }

  [[nodiscard]] base::StatusCode on_initialize(
    const MeshingRequest &request
  ) override
  {
    if(request.target_dimension != MeshingDimension::volume) {
      return core::detail::publish_error(
        base::StatusCode::unsupported,
        "Dummy Mesher only supports volume meshing requests."
      );
    }
    if(request.context_handle == sqmesh::invalid_handle) {
      return core::detail::publish_error(
        base::StatusCode::invalid_handle,
        "Meshing requires a valid context handle."
      );
    }
    if(request.session_handle == sqmesh::invalid_handle) {
      return core::detail::publish_error(
        base::StatusCode::invalid_handle,
        "Meshing requires a valid session handle."
      );
    }
    if(request.model_handle == sqmesh::invalid_handle) {
      return core::detail::publish_error(
        base::StatusCode::invalid_handle,
        "Meshing requires a valid model handle."
      );
    }

    request_ = request;
    return core::detail::clear_error_state();
  }

  [[nodiscard]] base::StatusCode on_configure(
    const ParameterDictionary &parameters
  ) override
  {
    double minimum_length = 0.0;
    if(parameters.try_get_number("minimum_length", minimum_length) &&
       minimum_length <= 0.0) {
      return core::detail::publish_error(
        base::StatusCode::invalid_argument,
        "Parameter 'minimum_length' must be positive."
      );
    }

    double growth_rate = 0.0;
    if(parameters.try_get_number("growth_rate", growth_rate) &&
       growth_rate <= 0.0) {
      return core::detail::publish_error(
        base::StatusCode::invalid_argument,
        "Parameter 'growth_rate' must be positive."
      );
    }

    std::string_view element_type;
    if(parameters.try_get_text("element_type", element_type) &&
       element_type != "tet" &&
       element_type != "tetra" &&
       element_type != "tetrahedron") {
      return core::detail::publish_error(
        base::StatusCode::unsupported,
        "Dummy Mesher only supports tetrahedral elements."
      );
    }

    return core::detail::clear_error_state();
  }

  [[nodiscard]] base::StatusCode on_generate(Domain &output) override
  {
    static_cast<void>(request_);
    output = make_dummy_tetra_domain();
    return core::detail::clear_error_state();
  }

private:
  MeshingRequest request_ {};
};

MeshingAlgorithmPtr create_dummy_mesher()
{
  return std::make_unique<DummyMesher>();
}

} // namespace

base::StatusCode MeshingAlgorithm::execute(
  const MeshingRequest &request,
  const ParameterDictionary &parameters,
  Domain &output
) noexcept
{
  try {
    reset();

    auto status = on_initialize(request);
    if(status != base::StatusCode::ok) {
      return status;
    }

    status = on_configure(parameters);
    if(status != base::StatusCode::ok) {
      return status;
    }

    status = on_generate(output);
    if(status != base::StatusCode::ok) {
      return status;
    }

    return core::detail::clear_error_state();
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
      "An unknown meshing algorithm error occurred."
    );
  }
}

void MeshingAlgorithm::reset() noexcept
{
}

void SharedEdgePolylineCache::clear() noexcept
{
  polylines_.clear();
}

void SharedEdgePolylineCache::reserve(std::size_t edge_count)
{
  polylines_.reserve(edge_count);
}

const SampledEdgePolyline *SharedEdgePolylineCache::find(
  geo::TopologyEntityId edge
) const noexcept
{
  if(edge.dimension != geo::TopologyDimension::edge) {
    return nullptr;
  }

  const auto it = polylines_.find(edge.index);
  if(it == polylines_.end()) {
    return nullptr;
  }

  return &it->second;
}

base::StatusCode SharedEdgePolylineCache::ensure_polyline(
  const geo::ModelView &model_view,
  geo::TopologyEntityId edge,
  const ResolvedMeshSizeControls &size_controls,
  const SampledEdgePolyline *&polyline
)
{
  polyline = nullptr;
  if(edge.dimension != geo::TopologyDimension::edge) {
    return core::detail::publish_error(
      base::StatusCode::internal_error,
      "Shared edge polyline sampling requires a valid geometry edge id."
    );
  }

  const auto *edge_view = model_view.find_edge(edge);
  if(edge_view == nullptr) {
    return core::detail::publish_error(
      base::StatusCode::internal_error,
      "Shared edge polyline sampling could not resolve the requested geometry edge from the model view."
    );
  }

  return ensure_polyline(*edge_view, size_controls, polyline);
}

base::StatusCode SharedEdgePolylineCache::ensure_polyline(
  const geo::EdgeView &edge_view,
  const ResolvedMeshSizeControls &size_controls,
  const SampledEdgePolyline *&polyline
)
{
  polyline = nullptr;
  if(edge_view.entity.dimension != geo::TopologyDimension::edge) {
    return core::detail::publish_error(
      base::StatusCode::internal_error,
      "Shared edge polyline sampling requires an edge view."
    );
  }

  if(const auto *cached = find(edge_view.entity)) {
    polyline = cached;
    return core::detail::clear_error_state();
  }

  geo::EdgeCurveSamples sampled_curve;
  auto status = geo::sample_edge_curve(
    edge_view,
    {
      std::max(effective_edge_target_size(size_controls, edge_view), kMinimumSharedEdgeTargetSize),
      1U,
    },
    sampled_curve
  );
  if(status != base::StatusCode::ok) {
    return status;
  }

  if(sampled_curve.curve.parameter_max <= sampled_curve.curve.parameter_min ||
     sampled_curve.curve.approximate_length <= 0.0 ||
     !std::isfinite(sampled_curve.curve.approximate_length) ||
     sampled_curve.samples.size() < 2U) {
    return core::detail::publish_error(
      base::StatusCode::unsupported,
      "Shared edge polyline sampling requires edges with a valid parameter range and positive length."
    );
  }

  auto [it, inserted] = polylines_.emplace(edge_view.entity.index, SampledEdgePolyline {});
  auto &stored = it->second;
  if(inserted) {
    stored.curve = sampled_curve.curve;
    stored.positions.reserve(sampled_curve.samples.size());
    for(const auto &sample : sampled_curve.samples) {
      stored.positions.push_back(sample.position);
    }
  }

  polyline = &stored;
  return core::detail::clear_error_state();
}

void BottomUpMeshingPipeline::reset() noexcept
{
  request_ = {};
  model_view_ = {};
  size_controls_ = {};
  edge_polyline_cache_.clear();
}

base::StatusCode BottomUpMeshingPipeline::on_initialize(
  const MeshingRequest &request
)
{
  auto status = validate_meshing_request(request, pipeline_dimension(), name());
  if(status != base::StatusCode::ok) {
    return status;
  }

  request_ = request;
  return on_pipeline_request(request_);
}

base::StatusCode BottomUpMeshingPipeline::on_generate(Domain &output)
{
  edge_polyline_cache_.clear();

  auto status = geo::model_view(
    request_.model_handle,
    model_view_,
    request_.context_handle
  );
  if(status != base::StatusCode::ok) {
    return status;
  }

  edge_polyline_cache_.reserve(model_view_.edges.size());

  status = resolve_mesh_size_controls(
    model_view_,
    request_.size_controls,
    pipeline_global_target_size(),
    size_controls_
  );
  if(status != base::StatusCode::ok) {
    return status;
  }

  status = on_prepare_pipeline();
  if(status != base::StatusCode::ok) {
    return status;
  }

  return dispatch_stages(output);
}

const MeshingRequest &BottomUpMeshingPipeline::request() const noexcept
{
  return request_;
}

const geo::ModelView &BottomUpMeshingPipeline::model_view() const noexcept
{
  return model_view_;
}

const ResolvedMeshSizeControls &BottomUpMeshingPipeline::size_controls() const noexcept
{
  return size_controls_;
}

base::StatusCode BottomUpMeshingPipeline::ensure_edge_polyline(
  geo::TopologyEntityId edge,
  const SampledEdgePolyline *&polyline
)
{
  return edge_polyline_cache_.ensure_polyline(model_view_, edge, size_controls_, polyline);
}

base::StatusCode BottomUpMeshingPipeline::ensure_edge_polyline(
  const geo::EdgeView &edge_view,
  const SampledEdgePolyline *&polyline
)
{
  return edge_polyline_cache_.ensure_polyline(edge_view, size_controls_, polyline);
}

base::StatusCode BottomUpMeshingPipeline::on_pipeline_request(
  const MeshingRequest &request
)
{
  static_cast<void>(request);
  return core::detail::clear_error_state();
}

base::StatusCode BottomUpMeshingPipeline::on_prepare_pipeline()
{
  return core::detail::clear_error_state();
}

base::StatusCode BottomUpMeshingPipeline::on_stage_vertices(Domain &output)
{
  static_cast<void>(output);
  return core::detail::clear_error_state();
}

base::StatusCode BottomUpMeshingPipeline::on_stage_edges(Domain &output)
{
  static_cast<void>(output);
  return core::detail::clear_error_state();
}

base::StatusCode BottomUpMeshingPipeline::on_stage_faces(Domain &output)
{
  static_cast<void>(output);
  return core::detail::clear_error_state();
}

base::StatusCode BottomUpMeshingPipeline::on_stage_regions(Domain &output)
{
  static_cast<void>(output);
  return core::detail::clear_error_state();
}

base::StatusCode BottomUpMeshingPipeline::dispatch_stages(Domain &output)
{
  auto status = on_stage_vertices(output);
  if(status != base::StatusCode::ok) {
    return status;
  }

  status = on_stage_edges(output);
  if(status != base::StatusCode::ok) {
    return status;
  }

  status = on_stage_faces(output);
  if(status != base::StatusCode::ok) {
    return status;
  }

  if(request_.target_dimension != MeshingDimension::volume) {
    return core::detail::clear_error_state();
  }

  return on_stage_regions(output);
}

void register_builtin_meshing_algorithms()
{
  static const bool registered = []() {
    std::lock_guard<std::mutex> lock(registry_mutex());
    static_cast<void>(registry().register_factory(kDummyMesherName, &create_dummy_mesher));
    static_cast<void>(registry().register_factory(kDummyMesherAlias, &create_dummy_mesher));
    static_cast<void>(
      registry().register_factory(
        kAutoCfdSurfaceMesherName,
        &create_auto_cfd_surface_mesher
      )
    );
    static_cast<void>(
      registry().register_factory(
        kAutoCfdSurfaceMesherAlias,
        &create_auto_cfd_surface_mesher
      )
    );
    static_cast<void>(
      registry().register_factory(
        kTetVolumeMesherName,
        &create_native_tet_volume_mesher
      )
    );
    static_cast<void>(
      registry().register_factory(
        kTetVolumeMesherAlias,
        &create_native_tet_volume_mesher
      )
    );
    static_cast<void>(
      registry().register_factory(
        kBLMesherName,
        &create_boundary_layer_mesher
      )
    );
    static_cast<void>(
      registry().register_factory(
        kBLMesherAlias,
        &create_boundary_layer_mesher
      )
    );
    return true;
  }();

  static_cast<void>(registered);
}

bool register_meshing_algorithm(
  std::string_view algorithm_name,
  MeshingAlgorithmFactory factory
)
{
  std::lock_guard<std::mutex> lock(registry_mutex());
  return registry().register_factory(algorithm_name, factory);
}

MeshingAlgorithmPtr create_meshing_algorithm(std::string_view algorithm_name)
{
  std::lock_guard<std::mutex> lock(registry_mutex());
  return registry().create(algorithm_name);
}

bool is_meshing_algorithm_registered(std::string_view algorithm_name)
{
  std::lock_guard<std::mutex> lock(registry_mutex());
  return registry().contains(algorithm_name);
}

} // namespace sqmesh::mesh::detail
