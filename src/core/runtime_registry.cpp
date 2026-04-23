// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "runtime_registry.hpp"

#include "../model/geometry_model_storage.hpp"
#include "../model/geometry_proxy_mesh.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdint>
#include <exception>
#include <mutex>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace sqmesh::core::detail {
namespace {

constexpr std::size_t kMaxErrorMessageSize = 512U;
constexpr std::uint64_t kHandleKindShift = 56U;
constexpr std::uint64_t kHandleSerialMask = (std::uint64_t{1} << kHandleKindShift) - 1U;

struct ErrorState final {
  base::StatusCode code = base::StatusCode::ok;
  std::array<char, kMaxErrorMessageSize> message{};
  std::size_t message_size = 0U;
};

thread_local ErrorState g_last_error;
thread_local base::ContextHandle g_current_context = sqmesh::invalid_handle;

void set_last_error(base::StatusCode code, std::string_view message) noexcept
{
  g_last_error.code = code;

  const auto max_size = g_last_error.message.size() - 1U;
  const auto size = message.size() < max_size ? message.size() : max_size;
  for(std::size_t index = 0; index < size; ++index) {
    g_last_error.message[index] = message[index];
  }
  g_last_error.message[size] = '\0';
  g_last_error.message_size = size;
}

base::StatusCode clear_last_error() noexcept
{
  set_last_error(base::StatusCode::ok, {});
  return base::StatusCode::ok;
}

base::StatusCode fail(
  base::StatusCode code,
  std::string_view message
) noexcept
{
  set_last_error(code, message);
  return code;
}

sqmesh::Handle make_handle(base::HandleKind kind, std::uint64_t serial) noexcept
{
  return (static_cast<sqmesh::Handle>(kind) << kHandleKindShift) |
         (serial & kHandleSerialMask);
}

base::HandleKind handle_kind(sqmesh::Handle handle) noexcept
{
  if(handle == sqmesh::invalid_handle) {
    return base::HandleKind::invalid;
  }

  return static_cast<base::HandleKind>((handle >> kHandleKindShift) & 0xFFU);
}

struct SessionRecord final {
  base::SessionHandle handle = sqmesh::invalid_handle;
  base::ContextHandle owner_context = sqmesh::invalid_handle;
};

struct ModelRecord final {
  geo::ModelHandle handle = sqmesh::invalid_handle;
  base::ContextHandle owner_context = sqmesh::invalid_handle;
  geo::ModelSummary summary {};
  model::detail::GeometryModelStoragePtr storage {};
  mesh::MeshHandle proxy_mesh = sqmesh::invalid_handle;
};

struct MeshRecord final {
  mesh::MeshHandle handle = sqmesh::invalid_handle;
  base::ContextHandle owner_context = sqmesh::invalid_handle;
  geo::ModelHandle source_model = sqmesh::invalid_handle;
  std::string algorithm_name {};
  mesh::MeshSummary summary {};
  std::shared_ptr<mesh::Domain> domain {};
};

using BorrowedModelStoragePtr = std::shared_ptr<const model::detail::GeometryModelStorage>;
using BorrowedMeshDomainPtr = std::shared_ptr<const mesh::Domain>;
using BorrowedMeshDomainMutablePtr = std::shared_ptr<mesh::Domain>;

struct ContextRecord final {
  base::ContextHandle handle = sqmesh::invalid_handle;
  base::SessionHandle current_session = sqmesh::invalid_handle;
  std::vector<base::SessionHandle> sessions;
  std::vector<geo::ModelHandle> models;
  std::vector<mesh::MeshHandle> meshes;
};

[[nodiscard]] std::shared_ptr<mesh::Domain> build_proxy_domain(
  const model::detail::GeometryCoarseProxyMesh &proxy_data
)
{
  // One shared node group for all proxy vertices, plus one entity group per
  // CAD face / CAD edge containing only its own triangles / segments. The
  // per-face partition lets downstream algorithms (curvature, proximity) run
  // without cross-face contamination.
  auto proxy_domain = std::make_shared<mesh::Domain>("proxy_domain");

  mesh::EntityGroupDefinition node_def;
  node_def.order = mesh::EntityOrder::node;
  node_def.name = "proxy_nodes";
  node_def.role = mesh::EntityGroupRole::geometric_proxy;
  const auto node_entity_group = proxy_domain->create_entity_group(std::move(node_def));

  proxy_domain->reserve_entity_group_storage(
    node_entity_group, proxy_data.nodes.size()
  );
  proxy_domain->set_source_topology_revision(proxy_data.topology_revision);

  std::vector<mesh::EntityRef> node_refs;
  node_refs.reserve(proxy_data.nodes.size());
  for(const auto &point : proxy_data.nodes) {
    node_refs.push_back(proxy_domain->add_node(node_entity_group, point));
  }

  // Per-CAD-face EntityGroups: one group per entry of face_triangle_ranges.
  for(std::size_t cad_face_idx = 0U;
      cad_face_idx < proxy_data.face_triangle_ranges.size();
      ++cad_face_idx) {
    const auto &range = proxy_data.face_triangle_ranges[cad_face_idx];
    if(range.count == 0U) {
      continue;
    }

    mesh::EntityGroupDefinition face_def;
    face_def.order = mesh::EntityOrder::face;
    face_def.name = "proxy_face_" + std::to_string(cad_face_idx);
    face_def.boundary = true;
    face_def.default_kind = mesh::EntityKind::face_triangle;
    face_def.role = mesh::EntityGroupRole::geometric_proxy;
    const auto face_entity_group =
      proxy_domain->create_entity_group(std::move(face_def));

    proxy_domain->reserve_entity_group_storage(
      face_entity_group, range.count, range.count * 3U
    );

    for(std::size_t offset = 0U; offset < range.count; ++offset) {
      const auto tri_index = range.offset + offset;
      if(tri_index >= proxy_data.triangles.size()) {
        continue;
      }
      const auto &tri = proxy_data.triangles[tri_index];
      const std::array<mesh::EntityRef, 3> face_nodes = {
        node_refs[tri[0]],
        node_refs[tri[1]],
        node_refs[tri[2]]
      };
      const auto face_ref =
        proxy_domain->add_triangle_face(face_entity_group, face_nodes);
      if(tri_index < proxy_data.triangle_face_owner.size()) {
        proxy_domain->set_face_topology_owner(
          face_ref, proxy_data.triangle_face_owner[tri_index]
        );
      }
    }
  }

  // Per-CAD-edge EntityGroups: one group per entry of edge_nodes.
  for(std::size_t edge_index = 0U;
      edge_index < proxy_data.edge_nodes.size();
      ++edge_index) {
    const auto &edge_seq = proxy_data.edge_nodes[edge_index];
    if(edge_seq.size() < 2U) {
      continue;
    }

    mesh::EntityGroupDefinition edge_def;
    edge_def.order = mesh::EntityOrder::edge;
    edge_def.name = "proxy_edge_" + std::to_string(edge_index);
    edge_def.default_kind = mesh::EntityKind::edge_line;
    edge_def.role = mesh::EntityGroupRole::geometric_proxy;
    const auto edge_entity_group =
      proxy_domain->create_entity_group(std::move(edge_def));

    const std::size_t segment_count = edge_seq.size() - 1U;
    proxy_domain->reserve_entity_group_storage(
      edge_entity_group, segment_count, segment_count * 2U
    );

    std::vector<mesh::EntityRef> current_edge_nodes;
    current_edge_nodes.reserve(edge_seq.size());
    for(const auto node_index : edge_seq) {
      if(node_index < node_refs.size()) {
        current_edge_nodes.push_back(node_refs[node_index]);
      }
    }

    auto owner_entity =
      geo::TopologyEntityId {geo::TopologyDimension::edge, geo::invalid_topology_index};
    if(edge_index < proxy_data.edge_topology_owner.size()) {
      owner_entity = proxy_data.edge_topology_owner[edge_index];
    }

    for(std::size_t segment_index = 0U;
        segment_index + 1U < current_edge_nodes.size();
        ++segment_index) {
      const auto edge_ref = proxy_domain->add_edge(
        edge_entity_group,
        {current_edge_nodes[segment_index], current_edge_nodes[segment_index + 1U]}
      );
      if(geo::is_valid(owner_entity)) {
        proxy_domain->set_edge_topology_owner(edge_ref, owner_entity);
      }
    }
  }

  return proxy_domain;
}

class Registry final {
public:
  base::StatusCode initialize(base::ContextHandle &context_handle)
  {
    const auto context = make_handle(base::HandleKind::context, next_serial_++);
    const auto session = make_handle(base::HandleKind::session, next_serial_++);

    ContextRecord context_record;
    context_record.handle = context;
    context_record.current_session = session;
    context_record.sessions.push_back(session);

    SessionRecord session_record;
    session_record.handle = session;
    session_record.owner_context = context;

    contexts_.emplace(context, context_record);
    sessions_.emplace(session, session_record);
    context_handle = context;
    g_current_context = context;

    return clear_last_error();
  }

  base::StatusCode shutdown(base::ContextHandle requested_context)
  {
    if(contexts_.empty()) {
      return fail(
        base::StatusCode::not_initialized,
        "The runtime is not initialized."
      );
    }

    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return g_last_error.code;
    }

    const auto context_it = contexts_.find(context);

    for(const auto session : context_it->second.sessions) {
      sessions_.erase(session);
    }
    for(const auto model : context_it->second.models) {
      models_.erase(model);
    }
    for(const auto mesh_handle : context_it->second.meshes) {
      meshes_.erase(mesh_handle);
    }
    contexts_.erase(context_it);
    if(g_current_context == context) {
      g_current_context = sqmesh::invalid_handle;
    }

    return clear_last_error();
  }

  base::StatusCode shutdown_all()
  {
    contexts_.clear();
    sessions_.clear();
    models_.clear();
    meshes_.clear();
    g_current_context = sqmesh::invalid_handle;
    return clear_last_error();
  }

  [[nodiscard]] bool is_initialized() const noexcept
  {
    return !contexts_.empty();
  }

  [[nodiscard]] bool has_context(base::ContextHandle context_handle) const noexcept
  {
    return contexts_.find(context_handle) != contexts_.end();
  }

  [[nodiscard]] base::ContextHandle current_context()
  {
    if(g_current_context == sqmesh::invalid_handle) {
      clear_last_error();
      return sqmesh::invalid_handle;
    }
    if(!has_context(g_current_context)) {
      g_current_context = sqmesh::invalid_handle;
      clear_last_error();
      return sqmesh::invalid_handle;
    }

    clear_last_error();
    return g_current_context;
  }

  base::SessionHandle current_session(base::ContextHandle requested_context)
  {
    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return sqmesh::invalid_handle;
    }

    clear_last_error();
    return contexts_.at(context).current_session;
  }

  base::StatusCode create_placeholder_model(
    base::ContextHandle requested_context,
    geo::ModelHandle &model_handle
  )
  {
    model_handle = sqmesh::invalid_handle;

    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return g_last_error.code;
    }

    auto &context_record = contexts_.at(context);
    const auto model = make_handle(base::HandleKind::model, next_serial_++);

    ModelRecord model_record;
    model_record.handle = model;
    model_record.owner_context = context;
    model_record.summary = geo::ModelSummary {};

    models_.emplace(model, std::move(model_record));
    context_record.models.push_back(model);
    model_handle = model;
    return clear_last_error();
  }

  base::StatusCode store_model(
    base::ContextHandle requested_context,
    model::detail::GeometryModelStoragePtr storage,
    geo::ModelSummary summary,
    geo::ModelHandle &model_handle
  )
  {
    model_handle = sqmesh::invalid_handle;

    if(!storage) {
      return fail(
        base::StatusCode::invalid_argument,
        "A geometry model storage object is required."
      );
    }

    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return g_last_error.code;
    }

    auto &context_record = contexts_.at(context);
    const auto model = make_handle(base::HandleKind::model, next_serial_++);

    ModelRecord model_record;
    model_record.handle = model;
    model_record.owner_context = context;
    model_record.summary = summary;
    model_record.storage = storage;

    model::detail::GeometryCoarseProxyMesh proxy_data;
    if(storage->coarse_proxy_mesh(proxy_data) == base::StatusCode::ok &&
       !proxy_data.empty()) {
      auto proxy_domain = build_proxy_domain(proxy_data);

      // Store the proxy domain as a MeshRecord
      const auto proxy_mesh_handle = make_handle(base::HandleKind::mesh, next_serial_++);
      const auto proxy_summary = proxy_domain->summary();

      MeshRecord proxy_mesh_record;
      proxy_mesh_record.handle = proxy_mesh_handle;
      proxy_mesh_record.owner_context = context;
      proxy_mesh_record.source_model = model;
      proxy_mesh_record.algorithm_name = "proxy";
      proxy_mesh_record.summary = proxy_summary;
      proxy_mesh_record.domain = std::move(proxy_domain);

      meshes_.emplace(proxy_mesh_handle, std::move(proxy_mesh_record));
      context_record.meshes.push_back(proxy_mesh_handle);
      model_record.proxy_mesh = proxy_mesh_handle;
    }

    models_.emplace(model, std::move(model_record));
    context_record.models.push_back(model);
    model_handle = model;
    return clear_last_error();
  }

  base::StatusCode update_model(
    base::ContextHandle requested_context,
    geo::ModelHandle model_handle,
    model::detail::GeometryModelStoragePtr storage,
    geo::ModelSummary summary
  )
  {
    if(!storage) {
      return fail(
        base::StatusCode::invalid_argument,
        "A geometry model storage object is required."
      );
    }

    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return g_last_error.code;
    }

    auto *model_record = require_model_record_mutable(model_handle, context);
    if(model_record == nullptr) {
      return g_last_error.code;
    }

    model_record->summary = summary;
    model_record->storage = std::move(storage);

    // If the proxy mesh was not created during initial store_model()
    // (e.g. geometry was broken before topo repair), retry now.
    if(model_record->proxy_mesh == sqmesh::invalid_handle) {
      model::detail::GeometryCoarseProxyMesh proxy_data;
      if(model_record->storage->coarse_proxy_mesh(proxy_data) == base::StatusCode::ok &&
         !proxy_data.empty()) {
        auto proxy_domain = build_proxy_domain(proxy_data);
        const auto proxy_mesh_handle = make_handle(base::HandleKind::mesh, next_serial_++);
        const auto proxy_summary = proxy_domain->summary();

        MeshRecord proxy_mesh_record;
        proxy_mesh_record.handle = proxy_mesh_handle;
        proxy_mesh_record.owner_context = context;
        proxy_mesh_record.source_model = model_handle;
        proxy_mesh_record.algorithm_name = "proxy";
        proxy_mesh_record.summary = proxy_summary;
        proxy_mesh_record.domain = std::move(proxy_domain);

        auto &context_record = contexts_.at(context);
        meshes_.emplace(proxy_mesh_handle, std::move(proxy_mesh_record));
        context_record.meshes.push_back(proxy_mesh_handle);
        model_record->proxy_mesh = proxy_mesh_handle;
      }
    }

    return clear_last_error();
  }

  BorrowedModelStoragePtr borrow_model_storage(
    base::ContextHandle requested_context,
    geo::ModelHandle model_handle
  )
  {
    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return {};
    }

    const auto *model_record = require_model_record(model_handle, context);
    if(model_record == nullptr) {
      return {};
    }
    if(!model_record->storage) {
      fail(
        base::StatusCode::unsupported,
        "The model handle is not backed by geometry storage."
      );
      return {};
    }

    clear_last_error();
    return model_record->storage;
  }

  model::detail::GeometryModelStoragePtr lookup_model_storage(
    base::ContextHandle requested_context,
    geo::ModelHandle model_handle
  )
  {
    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return {};
    }

    const auto *model_record = require_model_record(model_handle, context);
    if(model_record == nullptr) {
      return {};
    }

    clear_last_error();
    return model_record->storage;
  }

  base::StatusCode model_summary(
    base::ContextHandle requested_context,
    geo::ModelHandle model_handle,
    geo::ModelSummary &summary
  )
  {
    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return g_last_error.code;
    }

    const auto *model_record = require_model_record(model_handle, context);
    if(model_record == nullptr) {
      return g_last_error.code;
    }

    summary = model_record->summary;
    return clear_last_error();
  }

  base::StatusCode topology_snapshot(
    base::ContextHandle requested_context,
    geo::ModelHandle model_handle,
    geo::TopologySnapshot &snapshot
  )
  {
    snapshot = {};

    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return g_last_error.code;
    }

    const auto *model_record = require_model_record(model_handle, context);
    if(model_record == nullptr) {
      return g_last_error.code;
    }

    if(!model_record->storage) {
      return clear_last_error();
    }

    return model_record->storage->topology_snapshot(snapshot);
  }

  base::StatusCode topology_children(
    base::ContextHandle requested_context,
    geo::ModelHandle model_handle,
    geo::TopologyEntityId entity,
    std::vector<geo::TopologyEntityId> &children
  )
  {
    children.clear();

    if(!geo::is_valid(entity)) {
      return fail(
        base::StatusCode::invalid_argument,
        "Topology child lookup requires a valid topology entity id."
      );
    }

    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return g_last_error.code;
    }

    const auto *model_record = require_model_record(model_handle, context);
    if(model_record == nullptr) {
      return g_last_error.code;
    }

    if(!model_record->storage) {
      return fail(
        base::StatusCode::invalid_argument,
        "The geometry model does not expose topology entities."
      );
    }

    return model_record->storage->topology_children(entity, children);
  }

  base::StatusCode topology_parents(
    base::ContextHandle requested_context,
    geo::ModelHandle model_handle,
    geo::TopologyEntityId entity,
    std::vector<geo::TopologyEntityId> &parents
  )
  {
    parents.clear();

    if(!geo::is_valid(entity)) {
      return fail(
        base::StatusCode::invalid_argument,
        "Topology parent lookup requires a valid topology entity id."
      );
    }

    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return g_last_error.code;
    }

    const auto *model_record = require_model_record(model_handle, context);
    if(model_record == nullptr) {
      return g_last_error.code;
    }

    if(!model_record->storage) {
      return fail(
        base::StatusCode::invalid_argument,
        "The geometry model does not expose topology entities."
      );
    }

    return model_record->storage->topology_parents(entity, parents);
  }

  base::StatusCode face_uv_bounds(
    base::ContextHandle requested_context,
    geo::ModelHandle model_handle,
    geo::TopologyEntityId face_entity,
    geo::FaceUvBounds &bounds
  )
  {
    bounds = {};

    if(!geo::is_valid(face_entity)) {
      return fail(
        base::StatusCode::invalid_argument,
        "Face UV bounds lookup requires a valid topology entity id."
      );
    }

    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return g_last_error.code;
    }

    const auto *model_record = require_model_record(model_handle, context);
    if(model_record == nullptr) {
      return g_last_error.code;
    }

    if(!model_record->storage) {
      return fail(
        base::StatusCode::invalid_argument,
        "The geometry model does not expose face geometry."
      );
    }

    return model_record->storage->face_uv_bounds(face_entity, bounds);
  }

  base::StatusCode sample_face(
    base::ContextHandle requested_context,
    geo::ModelHandle model_handle,
    geo::TopologyEntityId face_entity,
    double u,
    double v,
    geo::FaceSample &sample
  )
  {
    sample = {};

    if(!geo::is_valid(face_entity)) {
      return fail(
        base::StatusCode::invalid_argument,
        "Face sampling requires a valid topology entity id."
      );
    }

    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return g_last_error.code;
    }

    const auto *model_record = require_model_record(model_handle, context);
    if(model_record == nullptr) {
      return g_last_error.code;
    }

    if(!model_record->storage) {
      return fail(
        base::StatusCode::invalid_argument,
        "The geometry model does not expose face geometry."
      );
    }

    return model_record->storage->sample_face(face_entity, u, v, sample);
  }

  base::StatusCode sample_face_curvature(
    base::ContextHandle requested_context,
    geo::ModelHandle model_handle,
    geo::TopologyEntityId face_entity,
    double u,
    double v,
    geo::FaceCurvatureSample &sample
  )
  {
    sample = {};

    if(!geo::is_valid(face_entity)) {
      return fail(
        base::StatusCode::invalid_argument,
        "Face curvature sampling requires a valid topology entity id."
      );
    }

    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return g_last_error.code;
    }

    const auto *model_record = require_model_record(model_handle, context);
    if(model_record == nullptr) {
      return g_last_error.code;
    }

    if(!model_record->storage) {
      return fail(
        base::StatusCode::invalid_argument,
        "The geometry model does not expose face curvature."
      );
    }

    return model_record->storage->sample_face_curvature(face_entity, u, v, sample);
  }

  base::StatusCode sample_face_derivatives(
    base::ContextHandle requested_context,
    geo::ModelHandle model_handle,
    geo::TopologyEntityId face_entity,
    double u,
    double v,
    geo::FaceDerivatives &sample
  )
  {
    sample = {};

    if(!geo::is_valid(face_entity)) {
      return fail(
        base::StatusCode::invalid_argument,
        "Face derivative sampling requires a valid topology entity id."
      );
    }

    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return g_last_error.code;
    }

    const auto *model_record = require_model_record(model_handle, context);
    if(model_record == nullptr) {
      return g_last_error.code;
    }

    if(!model_record->storage) {
      return fail(
        base::StatusCode::invalid_argument,
        "The geometry model does not expose face derivatives."
      );
    }

    return model_record->storage->sample_face_derivatives(face_entity, u, v, sample);
  }

  base::StatusCode project_point_to_face(
    base::ContextHandle requested_context,
    geo::ModelHandle model_handle,
    geo::TopologyEntityId face_entity,
    const geo::Point3 &point,
    geo::FaceProjection &projection
  )
  {
    projection = {};

    if(!geo::is_valid(face_entity)) {
      return fail(
        base::StatusCode::invalid_argument,
        "Face projection requires a valid topology entity id."
      );
    }

    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return g_last_error.code;
    }

    const auto *model_record = require_model_record(model_handle, context);
    if(model_record == nullptr) {
      return g_last_error.code;
    }

    if(!model_record->storage) {
      return fail(
        base::StatusCode::invalid_argument,
        "The geometry model does not expose face projection."
      );
    }

    return model_record->storage->project_point_to_face(face_entity, point, projection);
  }

  base::StatusCode recover_face_uv(
    base::ContextHandle requested_context,
    geo::ModelHandle model_handle,
    geo::TopologyEntityId face_entity,
    const geo::Point3 &point,
    geo::FaceUvMapping &mapping
  )
  {
    mapping = {};

    if(!geo::is_valid(face_entity)) {
      return fail(
        base::StatusCode::invalid_argument,
        "Face UV recovery requires a valid topology entity id."
      );
    }

    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return g_last_error.code;
    }

    const auto *model_record = require_model_record(model_handle, context);
    if(model_record == nullptr) {
      return g_last_error.code;
    }

    if(!model_record->storage) {
      return fail(
        base::StatusCode::invalid_argument,
        "The geometry model does not expose face inverse mapping."
      );
    }

    return model_record->storage->recover_face_uv(face_entity, point, mapping);
  }

  base::StatusCode recover_face_uv_from_edge(
    base::ContextHandle requested_context,
    geo::ModelHandle model_handle,
    geo::TopologyEntityId face_entity,
    geo::TopologyEntityId edge_entity,
    double edge_parameter,
    geo::FaceUvMapping &mapping
  )
  {
    mapping = {};

    if(!geo::is_valid(face_entity)) {
      return fail(
        base::StatusCode::invalid_argument,
        "Face UV recovery from edge requires a valid face entity id."
      );
    }
    if(!geo::is_valid(edge_entity)) {
      return fail(
        base::StatusCode::invalid_argument,
        "Face UV recovery from edge requires a valid edge entity id."
      );
    }

    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return g_last_error.code;
    }

    const auto *model_record = require_model_record(model_handle, context);
    if(model_record == nullptr) {
      return g_last_error.code;
    }

    if(!model_record->storage) {
      return fail(
        base::StatusCode::invalid_argument,
        "The geometry model does not expose edge pcurve UV mapping."
      );
    }

    return model_record->storage->recover_face_uv_from_edge(
      face_entity, edge_entity, edge_parameter, mapping
    );
  }

  base::StatusCode recover_face_uv_from_edge_use(
    base::ContextHandle requested_context,
    geo::ModelHandle model_handle,
    geo::TopologyEntityId face_entity,
    const geo::FaceBoundaryEdgeUse &edge_use,
    double edge_parameter,
    geo::FaceUvMapping &mapping
  )
  {
    mapping = {};

    if(!geo::is_valid(face_entity)) {
      return fail(
        base::StatusCode::invalid_argument,
        "Face UV recovery from edge use requires a valid face entity id."
      );
    }
    if(!geo::is_valid(edge_use.edge)) {
      return fail(
        base::StatusCode::invalid_argument,
        "Face UV recovery from edge use requires a valid edge entity id."
      );
    }

    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return g_last_error.code;
    }

    const auto *model_record = require_model_record(model_handle, context);
    if(model_record == nullptr) {
      return g_last_error.code;
    }

    if(!model_record->storage) {
      return fail(
        base::StatusCode::invalid_argument,
        "The geometry model does not expose face-local edge-use pcurve UV mapping."
      );
    }

    return model_record->storage->recover_face_uv_from_edge_use(
      face_entity, edge_use, edge_parameter, mapping
    );
  }

  base::StatusCode edge_curve_info(
    base::ContextHandle requested_context,
    geo::ModelHandle model_handle,
    geo::TopologyEntityId edge_entity,
    geo::EdgeCurveInfo &info
  )
  {
    info = {};

    if(!geo::is_valid(edge_entity)) {
      return fail(
        base::StatusCode::invalid_argument,
        "Edge geometry lookup requires a valid topology entity id."
      );
    }

    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return g_last_error.code;
    }

    const auto *model_record = require_model_record(model_handle, context);
    if(model_record == nullptr) {
      return g_last_error.code;
    }

    if(!model_record->storage) {
      return fail(
        base::StatusCode::invalid_argument,
        "The geometry model does not expose edge geometry."
      );
    }

    return model_record->storage->edge_curve_info(edge_entity, info);
  }

  base::StatusCode sample_edge_tangent(
    base::ContextHandle requested_context,
    geo::ModelHandle model_handle,
    geo::TopologyEntityId edge_entity,
    double parameter,
    geo::EdgeTangentSample &sample
  )
  {
    sample = {};

    if(!geo::is_valid(edge_entity)) {
      return fail(
        base::StatusCode::invalid_argument,
        "Edge tangent sampling requires a valid topology entity id."
      );
    }

    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return g_last_error.code;
    }

    const auto *model_record = require_model_record(model_handle, context);
    if(model_record == nullptr) {
      return g_last_error.code;
    }

    if(!model_record->storage) {
      return fail(
        base::StatusCode::invalid_argument,
        "The geometry model does not expose edge tangents."
      );
    }

    return model_record->storage->sample_edge_tangent(edge_entity, parameter, sample);
  }

  base::StatusCode face_boundary_loops(
    base::ContextHandle requested_context,
    geo::ModelHandle model_handle,
    geo::TopologyEntityId face_entity,
    geo::FaceBoundaryLoops &boundary
  )
  {
    boundary = {};

    if(!geo::is_valid(face_entity)) {
      return fail(
        base::StatusCode::invalid_argument,
        "Face boundary loop lookup requires a valid topology entity id."
      );
    }

    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return g_last_error.code;
    }

    const auto *model_record = require_model_record(model_handle, context);
    if(model_record == nullptr) {
      return g_last_error.code;
    }

    if(!model_record->storage) {
      return fail(
        base::StatusCode::invalid_argument,
        "The geometry model does not expose face boundary loops."
      );
    }

    return model_record->storage->face_boundary_loops(face_entity, boundary);
  }

  base::StatusCode feature_edges(
    base::ContextHandle requested_context,
    geo::ModelHandle model_handle,
    geo::FeatureEdgeReport &report,
    const geo::FeatureEdgeOptions &options
  )
  {
    report = {};

    if(options.feature_angle_degrees < 0.0 || options.feature_angle_degrees > 180.0) {
      return fail(
        base::StatusCode::invalid_argument,
        "Feature edge extraction requires feature_angle_degrees in [0, 180]."
      );
    }

    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return g_last_error.code;
    }

    const auto *model_record = require_model_record(model_handle, context);
    if(model_record == nullptr) {
      return g_last_error.code;
    }

    if(!model_record->storage) {
      return fail(
        base::StatusCode::invalid_argument,
        "The geometry model does not expose feature edge analysis."
      );
    }

    return model_record->storage->feature_edges(report, options);
  }

  base::StatusCode validate_model_handle(
    base::ContextHandle requested_context,
    geo::ModelHandle model_handle
  )
  {
    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return g_last_error.code;
    }

    if(require_model_record(model_handle, context) == nullptr) {
      return g_last_error.code;
    }

    return clear_last_error();
  }

  base::StatusCode store_mesh(
    base::ContextHandle requested_context,
    geo::ModelHandle model_handle,
    std::string_view algorithm_name,
    mesh::MeshSummary summary,
    std::shared_ptr<mesh::Domain> domain,
    mesh::MeshHandle &mesh_handle
  )
  {
    mesh_handle = sqmesh::invalid_handle;

    if(!domain) {
      return fail(
        base::StatusCode::invalid_argument,
        "A mesh domain storage object is required."
      );
    }

    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return g_last_error.code;
    }

    const auto *model_record = require_model_record(model_handle, context);
    if(model_record == nullptr) {
      return g_last_error.code;
    }
    if(model_record->storage) {
      geo::TopologySnapshot snapshot;
      const auto status = model_record->storage->topology_snapshot(snapshot);
      if(status != base::StatusCode::ok) {
        return status;
      }
      summary.source_topology_revision = snapshot.topology_revision;
      domain->set_source_topology_revision(snapshot.topology_revision);
    }
    else {
      summary.source_topology_revision = 0U;
      domain->set_source_topology_revision(0U);
    }

    auto &context_record = contexts_.at(context);
    const auto mesh_resource = make_handle(base::HandleKind::mesh, next_serial_++);

    MeshRecord mesh_record;
    mesh_record.handle = mesh_resource;
    mesh_record.owner_context = context;
    mesh_record.source_model = model_handle;
    mesh_record.algorithm_name = std::string(algorithm_name);
    mesh_record.summary = summary;
    mesh_record.domain = std::move(domain);

    meshes_.emplace(mesh_resource, std::move(mesh_record));
    context_record.meshes.push_back(mesh_resource);
    mesh_handle = mesh_resource;
    return clear_last_error();
  }

  base::StatusCode mesh_summary(
    base::ContextHandle requested_context,
    mesh::MeshHandle mesh_handle,
    mesh::MeshSummary &summary
  )
  {
    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return g_last_error.code;
    }

    auto *mesh_record = require_mesh_record_mutable(mesh_handle, context);
    if(mesh_record == nullptr) {
      return g_last_error.code;
    }
    const auto refresh_status = ensure_proxy_mesh_current(context, *mesh_record);
    if(refresh_status != base::StatusCode::ok) {
      return refresh_status;
    }

    summary = mesh_record->summary;
    return clear_last_error();
  }

  base::StatusCode mesh_domain_snapshot(
    base::ContextHandle requested_context,
    mesh::MeshHandle mesh_handle,
    mesh::Domain &domain
  )
  {
    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return g_last_error.code;
    }

    auto *mesh_record = require_mesh_record_mutable(mesh_handle, context);
    if(mesh_record == nullptr) {
      return g_last_error.code;
    }
    const auto refresh_status = ensure_proxy_mesh_current(context, *mesh_record);
    if(refresh_status != base::StatusCode::ok) {
      return refresh_status;
    }
    if(!mesh_record->domain) {
      return fail(
        base::StatusCode::internal_error,
        "The mesh handle does not reference a runtime mesh domain."
      );
    }

    domain = *mesh_record->domain;
    return clear_last_error();
  }

  BorrowedMeshDomainPtr borrow_mesh_domain(
    base::ContextHandle requested_context,
    mesh::MeshHandle mesh_handle
  )
  {
    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return {};
    }

    auto *mesh_record = require_mesh_record_mutable(mesh_handle, context);
    if(mesh_record == nullptr) {
      return {};
    }
    const auto refresh_status = ensure_proxy_mesh_current(context, *mesh_record);
    if(refresh_status != base::StatusCode::ok) {
      return {};
    }
    if(!mesh_record->domain) {
      fail(
        base::StatusCode::internal_error,
        "The mesh handle does not reference a runtime mesh domain."
      );
      return {};
    }

    clear_last_error();
    return mesh_record->domain;
  }

  BorrowedMeshDomainMutablePtr borrow_mesh_domain_mutable(
    base::ContextHandle requested_context,
    mesh::MeshHandle mesh_handle
  )
  {
    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return {};
    }

    auto *mesh_record = require_mesh_record_mutable(mesh_handle, context);
    if(mesh_record == nullptr) {
      return {};
    }
    const auto refresh_status = ensure_proxy_mesh_current(context, *mesh_record);
    if(refresh_status != base::StatusCode::ok) {
      return {};
    }
    if(!mesh_record->domain) {
      fail(
        base::StatusCode::internal_error,
        "The mesh handle does not reference a runtime mesh domain."
      );
      return {};
    }

    clear_last_error();
    return mesh_record->domain;
  }

  base::StatusCode model_proxy_mesh(
    base::ContextHandle requested_context,
    geo::ModelHandle model_handle,
    mesh::MeshHandle &mesh_handle
  )
  {
    mesh_handle = sqmesh::invalid_handle;

    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return g_last_error.code;
    }

    const auto *model_record = require_model_record(model_handle, context);
    if(model_record == nullptr) {
      return g_last_error.code;
    }

    if(model_record->proxy_mesh == sqmesh::invalid_handle) {
      return fail(
        base::StatusCode::unsupported,
        "The model does not have a proxy mesh."
      );
    }

    mesh_handle = model_record->proxy_mesh;
    auto *mesh_record = require_mesh_record_mutable(mesh_handle, context);
    if(mesh_record == nullptr) {
      return g_last_error.code;
    }
    const auto refresh_status = ensure_proxy_mesh_current(context, *mesh_record);
    if(refresh_status != base::StatusCode::ok) {
      return refresh_status;
    }
    return clear_last_error();
  }

  base::StatusCode refresh_mesh_summary(
    base::ContextHandle requested_context,
    mesh::MeshHandle mesh_handle
  )
  {
    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return g_last_error.code;
    }

    auto mesh_it = meshes_.find(mesh_handle);
    if(mesh_it == meshes_.end()) {
      return fail(
        base::StatusCode::invalid_handle,
        "Unknown mesh handle."
      );
    }
    if(mesh_it->second.owner_context != context) {
      return fail(
        base::StatusCode::owner_mismatch,
        "The mesh handle belongs to a different context."
      );
    }

    const auto refresh_status = ensure_proxy_mesh_current(context, mesh_it->second);
    if(refresh_status != base::StatusCode::ok) {
      return refresh_status;
    }
    if(!mesh_it->second.domain) {
      return fail(
        base::StatusCode::internal_error,
        "The mesh handle does not reference a runtime mesh domain."
      );
    }

    mesh_it->second.summary = mesh_it->second.domain->summary();
    return clear_last_error();
  }

  base::StatusCode mesh_node_count(
    base::ContextHandle requested_context,
    mesh::MeshHandle mesh_handle,
    std::size_t &count
  )
  {
    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return g_last_error.code;
    }

    auto *mesh_record = require_mesh_record_mutable(mesh_handle, context);
    if(mesh_record == nullptr) {
      return g_last_error.code;
    }
    const auto refresh_status = ensure_proxy_mesh_current(context, *mesh_record);
    if(refresh_status != base::StatusCode::ok) {
      return refresh_status;
    }

    count = mesh_record->summary.node_count;
    return clear_last_error();
  }

  base::StatusCode mesh_cell_count(
    base::ContextHandle requested_context,
    mesh::MeshHandle mesh_handle,
    std::size_t &count
  )
  {
    const auto context = require_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return g_last_error.code;
    }

    auto *mesh_record = require_mesh_record_mutable(mesh_handle, context);
    if(mesh_record == nullptr) {
      return g_last_error.code;
    }
    const auto refresh_status = ensure_proxy_mesh_current(context, *mesh_record);
    if(refresh_status != base::StatusCode::ok) {
      return refresh_status;
    }

    count = mesh_record->summary.cell_count;
    return clear_last_error();
  }

private:
  base::StatusCode ensure_proxy_mesh_current(
    base::ContextHandle owner_context,
    MeshRecord &mesh_record
  )
  {
    if(mesh_record.source_model == sqmesh::invalid_handle) {
      return clear_last_error();
    }

    auto *model_record = require_model_record_mutable(mesh_record.source_model, owner_context);
    if(model_record == nullptr) {
      return g_last_error.code;
    }
    if(model_record->proxy_mesh != mesh_record.handle || !model_record->storage) {
      return clear_last_error();
    }

    model::detail::GeometryCoarseProxyMesh proxy_data;
    const auto status = model_record->storage->coarse_proxy_mesh(proxy_data);
    if(status != base::StatusCode::ok) {
      return status;
    }
    if(proxy_data.empty()) {
      return clear_last_error();
    }

    if(mesh_record.domain &&
       mesh_record.domain->source_topology_revision() == proxy_data.topology_revision) {
      mesh_record.summary = mesh_record.domain->summary();
      return clear_last_error();
    }

    auto proxy_domain = build_proxy_domain(proxy_data);
    mesh_record.summary = proxy_domain->summary();
    mesh_record.domain = std::move(proxy_domain);
    return clear_last_error();
  }

  base::ContextHandle resolve_context_handle(base::ContextHandle requested_context)
  {
    if(requested_context != sqmesh::invalid_handle) {
      return requested_context;
    }
    if(g_current_context != sqmesh::invalid_handle && has_context(g_current_context)) {
      return g_current_context;
    }
    if(contexts_.empty()) {
      fail(
        base::StatusCode::not_initialized,
        "The runtime is not initialized."
      );
      return sqmesh::invalid_handle;
    }

    g_current_context = sqmesh::invalid_handle;
    fail(
      base::StatusCode::invalid_handle,
      "No current context is bound to the calling thread."
    );
    return sqmesh::invalid_handle;
  }

  base::ContextHandle require_context_handle(base::ContextHandle requested_context)
  {
    const auto context = resolve_context_handle(requested_context);
    if(context == sqmesh::invalid_handle) {
      return sqmesh::invalid_handle;
    }
    if(handle_kind(context) != base::HandleKind::context) {
      fail(
        base::StatusCode::invalid_handle,
        "Expected a context handle."
      );
      return sqmesh::invalid_handle;
    }

    const auto context_it = contexts_.find(context);
    if(context_it == contexts_.end()) {
      fail(
        base::StatusCode::invalid_handle,
        "Unknown context handle."
      );
      return sqmesh::invalid_handle;
    }

    return context;
  }

  const ModelRecord *require_model_record(
    geo::ModelHandle model_handle,
    base::ContextHandle owner_context
  )
  {
    if(handle_kind(model_handle) != base::HandleKind::model) {
      fail(
        base::StatusCode::invalid_handle,
        "Expected a model handle."
      );
      return nullptr;
    }

    const auto model_it = models_.find(model_handle);
    if(model_it == models_.end()) {
      fail(
        base::StatusCode::invalid_handle,
        "Unknown model handle."
      );
      return nullptr;
    }
    if(model_it->second.owner_context != owner_context) {
      fail(
        base::StatusCode::owner_mismatch,
        "The model handle belongs to a different context."
      );
      return nullptr;
    }

    return &model_it->second;
  }

  ModelRecord *require_model_record_mutable(
    geo::ModelHandle model_handle,
    base::ContextHandle owner_context
  )
  {
    if(handle_kind(model_handle) != base::HandleKind::model) {
      fail(
        base::StatusCode::invalid_handle,
        "Expected a model handle."
      );
      return nullptr;
    }

    const auto model_it = models_.find(model_handle);
    if(model_it == models_.end()) {
      fail(
        base::StatusCode::invalid_handle,
        "Unknown model handle."
      );
      return nullptr;
    }
    if(model_it->second.owner_context != owner_context) {
      fail(
        base::StatusCode::owner_mismatch,
        "The model handle belongs to a different context."
      );
      return nullptr;
    }

    return &model_it->second;
  }

  const MeshRecord *require_mesh_record(
    mesh::MeshHandle mesh_handle,
    base::ContextHandle owner_context
  )
  {
    if(handle_kind(mesh_handle) != base::HandleKind::mesh) {
      fail(
        base::StatusCode::invalid_handle,
        "Expected a mesh handle."
      );
      return nullptr;
    }

    const auto mesh_it = meshes_.find(mesh_handle);
    if(mesh_it == meshes_.end()) {
      fail(
        base::StatusCode::invalid_handle,
        "Unknown mesh handle."
      );
      return nullptr;
    }
    if(mesh_it->second.owner_context != owner_context) {
      fail(
        base::StatusCode::owner_mismatch,
        "The mesh handle belongs to a different context."
      );
      return nullptr;
    }

    return &mesh_it->second;
  }

  MeshRecord *require_mesh_record_mutable(
    mesh::MeshHandle mesh_handle,
    base::ContextHandle owner_context
  )
  {
    if(handle_kind(mesh_handle) != base::HandleKind::mesh) {
      fail(
        base::StatusCode::invalid_handle,
        "Expected a mesh handle."
      );
      return nullptr;
    }

    const auto mesh_it = meshes_.find(mesh_handle);
    if(mesh_it == meshes_.end()) {
      fail(
        base::StatusCode::invalid_handle,
        "Unknown mesh handle."
      );
      return nullptr;
    }
    if(mesh_it->second.owner_context != owner_context) {
      fail(
        base::StatusCode::owner_mismatch,
        "The mesh handle belongs to a different context."
      );
      return nullptr;
    }

    return &mesh_it->second;
  }

  std::uint64_t next_serial_ = 1U;
  std::unordered_map<base::ContextHandle, ContextRecord> contexts_;
  std::unordered_map<base::SessionHandle, SessionRecord> sessions_;
  std::unordered_map<geo::ModelHandle, ModelRecord> models_;
  std::unordered_map<mesh::MeshHandle, MeshRecord> meshes_;
};

Registry &registry() noexcept
{
  static Registry instance;
  return instance;
}

std::mutex &registry_mutex() noexcept
{
  static std::mutex instance;
  return instance;
}

template <typename Callback>
auto with_registry(Callback &&callback) noexcept -> decltype(callback(registry()))
{
  try {
    std::lock_guard<std::mutex> lock(registry_mutex());
    return callback(registry());
  }
  catch(const std::exception &exception) {
    set_last_error(base::StatusCode::internal_error, exception.what());
  }
  catch(...) {
    set_last_error(
      base::StatusCode::internal_error,
      "An unknown runtime error occurred."
    );
  }

  using ResultType = decltype(callback(registry()));
  if constexpr(std::is_same_v<ResultType, base::StatusCode>) {
    return base::StatusCode::internal_error;
  }
  else if constexpr(std::is_same_v<ResultType, bool>) {
    return false;
  }
  else if constexpr(std::is_same_v<ResultType, sqmesh::Handle>) {
    return sqmesh::invalid_handle;
  }
  else {
    return ResultType{};
  }
}

} // namespace

base::StatusCode initialize(base::ContextHandle &context_handle) noexcept
{
  context_handle = sqmesh::invalid_handle;
  return with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.initialize(context_handle);
    }
  );
}

base::StatusCode shutdown(base::ContextHandle context_handle) noexcept
{
  return with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.shutdown(context_handle);
    }
  );
}

base::StatusCode shutdown_all() noexcept
{
  return with_registry(
    [](Registry &registry_instance) {
      return registry_instance.shutdown_all();
    }
  );
}

bool is_initialized() noexcept
{
  return with_registry(
    [](Registry &registry_instance) {
      return registry_instance.is_initialized();
    }
  );
}

base::ContextHandle current_context() noexcept
{
  return with_registry(
    [](Registry &registry_instance) {
      return registry_instance.current_context();
    }
  );
}

base::SessionHandle current_session(base::ContextHandle context_handle) noexcept
{
  return with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.current_session(context_handle);
    }
  );
}

base::StatusCode clear_error_state() noexcept
{
  return clear_last_error();
}

base::StatusCode publish_error(
  base::StatusCode code,
  std::string_view message
) noexcept
{
  return fail(code, message);
}

base::StatusCode create_placeholder_model(
  geo::ModelHandle &model_handle,
  base::ContextHandle context_handle
) noexcept
{
  model_handle = sqmesh::invalid_handle;
  return with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.create_placeholder_model(context_handle, model_handle);
    }
  );
}

base::StatusCode store_model(
  model::detail::GeometryModelStoragePtr storage,
  geo::ModelSummary summary,
  geo::ModelHandle &model_handle,
  base::ContextHandle context_handle
) noexcept
{
  model_handle = sqmesh::invalid_handle;
  return with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.store_model(
        context_handle,
        std::move(storage),
        summary,
        model_handle
      );
    }
  );
}

base::StatusCode update_model(
  geo::ModelHandle model_handle,
  model::detail::GeometryModelStoragePtr storage,
  geo::ModelSummary summary,
  base::ContextHandle context_handle
) noexcept
{
  return with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.update_model(
        context_handle,
        model_handle,
        std::move(storage),
        summary
      );
    }
  );
}

model::detail::GeometryModelStoragePtr lookup_model_storage(
  geo::ModelHandle model_handle,
  base::ContextHandle context_handle
) noexcept
{
  return with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.lookup_model_storage(context_handle, model_handle);
    }
  );
}

base::StatusCode with_model_storage(
  geo::ModelHandle model_handle,
  ModelStorageBorrowCallback callback,
  base::ContextHandle context_handle
) noexcept
{
  if(!callback) {
    return fail(
      base::StatusCode::invalid_argument,
      "Runtime model borrowing requires a callback."
    );
  }

  const auto storage = with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.borrow_model_storage(context_handle, model_handle);
    }
  );
  if(!storage) {
    return g_last_error.code;
  }

  try {
    clear_last_error();
    const auto status = callback(*storage);
    if(status == base::StatusCode::ok) {
      return clear_last_error();
    }
    if(g_last_error.code != status) {
      set_last_error(status, {});
    }
    return status;
  }
  catch(const std::exception &exception) {
    return fail(base::StatusCode::internal_error, exception.what());
  }
  catch(...) {
    return fail(
      base::StatusCode::internal_error,
      "An unknown runtime error occurred while borrowing model storage."
    );
  }
}

base::StatusCode model_summary(
  geo::ModelHandle model_handle,
  geo::ModelSummary &summary,
  base::ContextHandle context_handle
) noexcept
{
  summary = {};
  return with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.model_summary(context_handle, model_handle, summary);
    }
  );
}

base::StatusCode topology_snapshot(
  geo::ModelHandle model_handle,
  geo::TopologySnapshot &snapshot,
  base::ContextHandle context_handle
) noexcept
{
  snapshot = {};
  return with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.topology_snapshot(context_handle, model_handle, snapshot);
    }
  );
}

base::StatusCode topology_children(
  geo::ModelHandle model_handle,
  geo::TopologyEntityId entity,
  std::vector<geo::TopologyEntityId> &children,
  base::ContextHandle context_handle
) noexcept
{
  children.clear();
  return with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.topology_children(
        context_handle,
        model_handle,
        entity,
        children
      );
    }
  );
}

base::StatusCode topology_parents(
  geo::ModelHandle model_handle,
  geo::TopologyEntityId entity,
  std::vector<geo::TopologyEntityId> &parents,
  base::ContextHandle context_handle
) noexcept
{
  parents.clear();
  return with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.topology_parents(
        context_handle,
        model_handle,
        entity,
        parents
      );
    }
  );
}

base::StatusCode face_uv_bounds(
  geo::ModelHandle model_handle,
  geo::TopologyEntityId face_entity,
  geo::FaceUvBounds &bounds,
  base::ContextHandle context_handle
) noexcept
{
  bounds = {};
  return with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.face_uv_bounds(
        context_handle,
        model_handle,
        face_entity,
        bounds
      );
    }
  );
}

base::StatusCode sample_face(
  geo::ModelHandle model_handle,
  geo::TopologyEntityId face_entity,
  double u,
  double v,
  geo::FaceSample &sample,
  base::ContextHandle context_handle
) noexcept
{
  sample = {};
  return with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.sample_face(
        context_handle,
        model_handle,
        face_entity,
        u,
        v,
        sample
      );
    }
  );
}

base::StatusCode sample_face_curvature(
  geo::ModelHandle model_handle,
  geo::TopologyEntityId face_entity,
  double u,
  double v,
  geo::FaceCurvatureSample &sample,
  base::ContextHandle context_handle
) noexcept
{
  sample = {};
  return with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.sample_face_curvature(
        context_handle,
        model_handle,
        face_entity,
        u,
        v,
        sample
      );
    }
  );
}

base::StatusCode sample_face_derivatives(
  geo::ModelHandle model_handle,
  geo::TopologyEntityId face_entity,
  double u,
  double v,
  geo::FaceDerivatives &sample,
  base::ContextHandle context_handle
) noexcept
{
  sample = {};
  return with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.sample_face_derivatives(
        context_handle,
        model_handle,
        face_entity,
        u,
        v,
        sample
      );
    }
  );
}

base::StatusCode project_point_to_face(
  geo::ModelHandle model_handle,
  geo::TopologyEntityId face_entity,
  const geo::Point3 &point,
  geo::FaceProjection &projection,
  base::ContextHandle context_handle
) noexcept
{
  projection = {};
  return with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.project_point_to_face(
        context_handle,
        model_handle,
        face_entity,
        point,
        projection
      );
    }
  );
}

base::StatusCode recover_face_uv(
  geo::ModelHandle model_handle,
  geo::TopologyEntityId face_entity,
  const geo::Point3 &point,
  geo::FaceUvMapping &mapping,
  base::ContextHandle context_handle
) noexcept
{
  mapping = {};
  return with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.recover_face_uv(
        context_handle,
        model_handle,
        face_entity,
        point,
        mapping
      );
    }
  );
}

base::StatusCode recover_face_uv_from_edge(
  geo::ModelHandle model_handle,
  geo::TopologyEntityId face_entity,
  geo::TopologyEntityId edge_entity,
  double edge_parameter,
  geo::FaceUvMapping &mapping,
  base::ContextHandle context_handle
) noexcept
{
  mapping = {};
  return with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.recover_face_uv_from_edge(
        context_handle,
        model_handle,
        face_entity,
        edge_entity,
        edge_parameter,
        mapping
      );
    }
  );
}

base::StatusCode recover_face_uv_from_edge_use(
  geo::ModelHandle model_handle,
  geo::TopologyEntityId face_entity,
  const geo::FaceBoundaryEdgeUse &edge_use,
  double edge_parameter,
  geo::FaceUvMapping &mapping,
  base::ContextHandle context_handle
) noexcept
{
  mapping = {};
  return with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.recover_face_uv_from_edge_use(
        context_handle,
        model_handle,
        face_entity,
        edge_use,
        edge_parameter,
        mapping
      );
    }
  );
}

base::StatusCode edge_curve_info(
  geo::ModelHandle model_handle,
  geo::TopologyEntityId edge_entity,
  geo::EdgeCurveInfo &info,
  base::ContextHandle context_handle
) noexcept
{
  info = {};
  return with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.edge_curve_info(
        context_handle,
        model_handle,
        edge_entity,
        info
      );
    }
  );
}

base::StatusCode sample_edge_tangent(
  geo::ModelHandle model_handle,
  geo::TopologyEntityId edge_entity,
  double parameter,
  geo::EdgeTangentSample &sample,
  base::ContextHandle context_handle
) noexcept
{
  sample = {};
  return with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.sample_edge_tangent(
        context_handle,
        model_handle,
        edge_entity,
        parameter,
        sample
      );
    }
  );
}

base::StatusCode face_boundary_loops(
  geo::ModelHandle model_handle,
  geo::TopologyEntityId face_entity,
  geo::FaceBoundaryLoops &boundary,
  base::ContextHandle context_handle
) noexcept
{
  boundary = {};
  return with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.face_boundary_loops(
        context_handle,
        model_handle,
        face_entity,
        boundary
      );
    }
  );
}

base::StatusCode feature_edges(
  geo::ModelHandle model_handle,
  geo::FeatureEdgeReport &report,
  const geo::FeatureEdgeOptions &options,
  base::ContextHandle context_handle
) noexcept
{
  report = {};
  return with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.feature_edges(
        context_handle,
        model_handle,
        report,
        options
      );
    }
  );
}

base::StatusCode validate_model_handle(
  geo::ModelHandle model_handle,
  base::ContextHandle context_handle
) noexcept
{
  return with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.validate_model_handle(context_handle, model_handle);
    }
  );
}

base::StatusCode store_mesh(
  geo::ModelHandle model_handle,
  std::string_view algorithm_name,
  mesh::MeshSummary summary,
  std::shared_ptr<mesh::Domain> domain,
  mesh::MeshHandle &mesh_handle,
  base::ContextHandle context_handle
) noexcept
{
  mesh_handle = sqmesh::invalid_handle;
  return with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.store_mesh(
        context_handle,
        model_handle,
        algorithm_name,
        summary,
        std::move(domain),
        mesh_handle
      );
    }
  );
}

base::StatusCode mesh_summary(
  mesh::MeshHandle mesh_handle,
  mesh::MeshSummary &summary,
  base::ContextHandle context_handle
) noexcept
{
  summary = {};
  return with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.mesh_summary(context_handle, mesh_handle, summary);
    }
  );
}

base::StatusCode mesh_domain_snapshot(
  mesh::MeshHandle mesh_handle,
  mesh::Domain &domain,
  base::ContextHandle context_handle
) noexcept
{
  return with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.mesh_domain_snapshot(context_handle, mesh_handle, domain);
    }
  );
}

base::StatusCode with_mesh_domain(
  mesh::MeshHandle mesh_handle,
  MeshDomainBorrowCallback callback,
  base::ContextHandle context_handle
) noexcept
{
  if(!callback) {
    return fail(
      base::StatusCode::invalid_argument,
      "Runtime mesh borrowing requires a callback."
    );
  }

  const auto domain = with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.borrow_mesh_domain(context_handle, mesh_handle);
    }
  );
  if(!domain) {
    return g_last_error.code;
  }

  try {
    clear_last_error();
    const auto status = callback(*domain);
    if(status == base::StatusCode::ok) {
      return clear_last_error();
    }
    if(g_last_error.code != status) {
      set_last_error(status, {});
    }
    return status;
  }
  catch(const std::exception &exception) {
    return fail(base::StatusCode::internal_error, exception.what());
  }
  catch(...) {
    return fail(
      base::StatusCode::internal_error,
      "An unknown runtime error occurred while borrowing mesh storage."
    );
  }
}

base::StatusCode with_mesh_domain_mutable(
  mesh::MeshHandle mesh_handle,
  MeshDomainMutableBorrowCallback callback,
  base::ContextHandle context_handle
) noexcept
{
  if(!callback) {
    return fail(
      base::StatusCode::invalid_argument,
      "Runtime mutable mesh borrowing requires a callback."
    );
  }

  const auto domain = with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.borrow_mesh_domain_mutable(context_handle, mesh_handle);
    }
  );
  if(!domain) {
    return g_last_error.code;
  }

  try {
    clear_last_error();
    const auto status = callback(*domain);
    if(status != base::StatusCode::ok) {
      if(g_last_error.code != status) {
        set_last_error(status, {});
      }
      return status;
    }

    // Auto-refresh the cached MeshSummary after successful mutable borrow
    return with_registry(
      [&](Registry &registry_instance) {
        return registry_instance.refresh_mesh_summary(context_handle, mesh_handle);
      }
    );
  }
  catch(const std::exception &exception) {
    return fail(base::StatusCode::internal_error, exception.what());
  }
  catch(...) {
    return fail(
      base::StatusCode::internal_error,
      "An unknown runtime error occurred while mutably borrowing mesh storage."
    );
  }
}

base::StatusCode refresh_mesh_summary(
  mesh::MeshHandle mesh_handle,
  base::ContextHandle context_handle
) noexcept
{
  return with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.refresh_mesh_summary(context_handle, mesh_handle);
    }
  );
}

base::StatusCode model_proxy_mesh(
  geo::ModelHandle model_handle,
  mesh::MeshHandle &mesh_handle,
  base::ContextHandle context_handle
) noexcept
{
  mesh_handle = sqmesh::invalid_handle;
  return with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.model_proxy_mesh(context_handle, model_handle, mesh_handle);
    }
  );
}
base::StatusCode mesh_node_count(
  mesh::MeshHandle mesh_handle,
  std::size_t &count,
  base::ContextHandle context_handle
) noexcept
{
  count = 0U;
  return with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.mesh_node_count(context_handle, mesh_handle, count);
    }
  );
}

base::StatusCode mesh_cell_count(
  mesh::MeshHandle mesh_handle,
  std::size_t &count,
  base::ContextHandle context_handle
) noexcept
{
  count = 0U;
  return with_registry(
    [&](Registry &registry_instance) {
      return registry_instance.mesh_cell_count(context_handle, mesh_handle, count);
    }
  );
}

base::StatusCode last_error_code() noexcept
{
  return g_last_error.code;
}

std::string_view last_error_message() noexcept
{
  return {g_last_error.message.data(), g_last_error.message_size};
}

const char *status_code_name(base::StatusCode code) noexcept
{
  switch(code) {
  case base::StatusCode::ok:
    return "ok";
  case base::StatusCode::invalid_argument:
    return "invalid_argument";
  case base::StatusCode::not_initialized:
    return "not_initialized";
  case base::StatusCode::invalid_handle:
    return "invalid_handle";
  case base::StatusCode::owner_mismatch:
    return "owner_mismatch";
  case base::StatusCode::io_error:
    return "io_error";
  case base::StatusCode::unsupported:
    return "unsupported";
  case base::StatusCode::internal_error:
    return "internal_error";
  default:
    return "unknown_status";
  }
}

} // namespace sqmesh::core::detail
