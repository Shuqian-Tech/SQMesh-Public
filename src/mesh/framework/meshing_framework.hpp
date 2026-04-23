// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#pragma once

#include "../sizing/mesh_size_controls.hpp"
#include "sqmesh/mesh/api.hpp"

#include <array>
#include <cstdint>
#include <memory>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace sqmesh::mesh::detail {

enum class MeshingDimension : std::uint8_t {
  surface = 2,
  volume = 3,
};

struct MeshingRequest final {
  base::ContextHandle context_handle = sqmesh::invalid_handle;
  base::SessionHandle session_handle = sqmesh::invalid_handle;
  geo::ModelHandle model_handle = sqmesh::invalid_handle;
  MeshingDimension target_dimension = MeshingDimension::volume;
  MeshSizeControls size_controls {};
};

class MeshingAlgorithm
{
public:
  virtual ~MeshingAlgorithm() = default;

  [[nodiscard]] base::StatusCode execute(
    const MeshingRequest &request,
    const ParameterDictionary &parameters,
    Domain &output
  ) noexcept;

  [[nodiscard]] virtual std::string_view name() const noexcept = 0;
  [[nodiscard]] virtual std::string_view entity_group_prefix() const noexcept { return {}; }

protected:
  virtual void reset() noexcept;
  [[nodiscard]] virtual base::StatusCode on_initialize(
    const MeshingRequest &request
  ) = 0;
  [[nodiscard]] virtual base::StatusCode on_configure(
    const ParameterDictionary &parameters
  ) = 0;
  [[nodiscard]] virtual base::StatusCode on_generate(Domain &output) = 0;
};

// Phase 1.1 keeps the shared bottom-up layer deliberately narrow: cached edge
// polylines expose only 3D sampled positions plus the source curve's endpoint
// and length data. UV recovery, loop assembly, seam handling, and face-local
// triangulation remain with the consuming algorithm.
struct SampledEdgePolyline final {
  geo::EdgeCurveInfo curve {};
  std::vector<geo::Point3> positions {};
};

class SharedEdgePolylineCache final
{
public:
  void clear() noexcept;
  void reserve(std::size_t edge_count);

  [[nodiscard]] const SampledEdgePolyline *find(
    geo::TopologyEntityId edge
  ) const noexcept;

  [[nodiscard]] base::StatusCode ensure_polyline(
    const geo::ModelView &model_view,
    geo::TopologyEntityId edge,
    const ResolvedMeshSizeControls &size_controls,
    const SampledEdgePolyline *&polyline
  );

  [[nodiscard]] base::StatusCode ensure_polyline(
    const geo::EdgeView &edge_view,
    const ResolvedMeshSizeControls &size_controls,
    const SampledEdgePolyline *&polyline
  );

private:
  std::unordered_map<std::uint32_t, SampledEdgePolyline> polylines_ {};
};

// `BottomUpMeshingPipeline` is an optional intermediate base for algorithms
// that genuinely advance from geometry boundary information toward interior
// stages. It is not the universal meshing superclass for SQMesh.
class BottomUpMeshingPipeline : public MeshingAlgorithm
{
protected:
  void reset() noexcept override;

  [[nodiscard]] base::StatusCode on_initialize(
    const MeshingRequest &request
  ) override final;

  [[nodiscard]] base::StatusCode on_generate(Domain &output) override final;

  [[nodiscard]] const MeshingRequest &request() const noexcept;
  [[nodiscard]] const geo::ModelView &model_view() const noexcept;
  [[nodiscard]] const ResolvedMeshSizeControls &size_controls() const noexcept;

  [[nodiscard]] base::StatusCode ensure_edge_polyline(
    geo::TopologyEntityId edge,
    const SampledEdgePolyline *&polyline
  );

  [[nodiscard]] base::StatusCode ensure_edge_polyline(
    const geo::EdgeView &edge_view,
    const SampledEdgePolyline *&polyline
  );

  [[nodiscard]] virtual MeshingDimension pipeline_dimension() const noexcept = 0;
  [[nodiscard]] virtual double pipeline_global_target_size() const noexcept = 0;

  [[nodiscard]] virtual base::StatusCode on_pipeline_request(
    const MeshingRequest &request
  );

  [[nodiscard]] virtual base::StatusCode on_prepare_pipeline();
  [[nodiscard]] virtual base::StatusCode on_stage_vertices(Domain &output);
  [[nodiscard]] virtual base::StatusCode on_stage_edges(Domain &output);
  [[nodiscard]] virtual base::StatusCode on_stage_faces(Domain &output);
  [[nodiscard]] virtual base::StatusCode on_stage_regions(Domain &output);

private:
  [[nodiscard]] base::StatusCode dispatch_stages(Domain &output);

  MeshingRequest request_ {};
  geo::ModelView model_view_ {};
  ResolvedMeshSizeControls size_controls_ {};
  SharedEdgePolylineCache edge_polyline_cache_ {};
};

using MeshingAlgorithmPtr = std::unique_ptr<MeshingAlgorithm>;
using MeshingAlgorithmFactory = MeshingAlgorithmPtr (*)();

void register_builtin_meshing_algorithms();
[[nodiscard]] bool register_meshing_algorithm(
  std::string_view algorithm_name,
  MeshingAlgorithmFactory factory
);
[[nodiscard]] MeshingAlgorithmPtr create_meshing_algorithm(
  std::string_view algorithm_name
);
[[nodiscard]] bool is_meshing_algorithm_registered(
  std::string_view algorithm_name
);

} // namespace sqmesh::mesh::detail
