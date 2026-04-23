// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "../framework/meshing_framework.hpp"

#include "core/runtime_registry.hpp"

#include <exception>
#include <memory>
#include <utility>
#include <variant>

namespace sqmesh::mesh {

namespace {

void set_auto_cfd_proximity_parameters(
  ParameterDictionary &parameters,
  std::string_view prefix,
  bool self_proximity,
  bool pid_proximity,
  double maximum_normals_angle,
  double length_to_gap_ratio,
  double proximity_minimum_length
) noexcept
{
  const auto make_key = [&](std::string_view suffix) {
    std::string key(prefix);
    key += suffix;
    return key;
  };

  const bool proximity_enabled = self_proximity || pid_proximity;
  parameters.set_boolean(make_key("proximity"), proximity_enabled);
  parameters.set_boolean(make_key("self_proximity"), self_proximity);
  parameters.set_boolean(make_key("pid_proximity"), pid_proximity);
  parameters.set_number(make_key("maximum_normals_angle"), maximum_normals_angle);
  parameters.set_number(make_key("length_to_gap_ratio"), length_to_gap_ratio);
  parameters.set_number(
    make_key("proximity_minimum_length"),
    proximity_minimum_length
  );
}

base::StatusCode create_mesh_impl(
  geo::ModelHandle model_handle,
  std::string_view algorithm_name,
  const MeshingOptions &options,
  detail::MeshingDimension target_dimension,
  MeshHandle &mesh_handle,
  base::ContextHandle context_handle
) noexcept
{
  mesh_handle = sqmesh::invalid_handle;

  try {
    detail::register_builtin_meshing_algorithms();

    if(algorithm_name.empty()) {
      return core::detail::publish_error(
        base::StatusCode::invalid_argument,
        "Meshing requires a non-empty algorithm name."
      );
    }

    const auto status = core::detail::validate_model_handle(model_handle, context_handle);
    if(status != base::StatusCode::ok) {
      return status;
    }

    auto resolved_context = context_handle;
    if(resolved_context == sqmesh::invalid_handle) {
      resolved_context = base::current_context();
      if(resolved_context == sqmesh::invalid_handle) {
        return base::last_error_code();
      }
    }

    const auto session_handle = base::current_session(resolved_context);
    if(session_handle == sqmesh::invalid_handle) {
      return base::last_error_code();
    }

    auto algorithm = detail::create_meshing_algorithm(algorithm_name);
    if(!algorithm) {
      return core::detail::publish_error(
        base::StatusCode::unsupported,
        "Requested meshing algorithm is not registered."
      );
    }

    detail::MeshingRequest request;
    request.context_handle = resolved_context;
    request.session_handle = session_handle;
    request.model_handle = model_handle;
    request.target_dimension = target_dimension;
    request.size_controls = options.size_controls;

    // Try to get the model's Proxy Domain (the unified global Domain)
    MeshHandle proxy_handle = sqmesh::invalid_handle;
    auto proxy_status = core::detail::model_proxy_mesh(
      model_handle, proxy_handle, resolved_context
    );

    if(proxy_status == base::StatusCode::ok && proxy_handle != sqmesh::invalid_handle) {
      // Use the unified Proxy Domain: clear old entity_groups, then generate into it
      const auto prefix = algorithm->entity_group_prefix();

      auto mutable_status = core::detail::with_mesh_domain_mutable(
        proxy_handle,
        [&](Domain &domain) {
          // Remove previously generated entity_groups for this algorithm
          if(!prefix.empty()) {
            // Volume meshes are downstream of the generated surface mesh inside
            // the shared proxy-backed Domain, so a surface rerun must first
            // invalidate dependent volume entity_groups before replacing surface_*.
            if(prefix == "surface_") {
              domain.remove_entity_groups_with_prefix("volume_");
            }
            domain.remove_entity_groups_with_prefix(prefix);
          }
          // Execute the algorithm directly on the shared Domain
          return algorithm->execute(request, options.parameters, domain);
        },
        resolved_context
      );

      if(mutable_status != base::StatusCode::ok) {
        return mutable_status;
      }

      // Refresh the cached summary to reflect the new domain state
      (void)core::detail::refresh_mesh_summary(proxy_handle, resolved_context);

      mesh_handle = proxy_handle;
      return core::detail::clear_error_state();
    }

    // Fallback: no proxy mesh, create a standalone Domain (legacy path)
    Domain output_domain;
    const auto execute_status = algorithm->execute(
      request,
      options.parameters,
      output_domain
    );
    if(execute_status != base::StatusCode::ok) {
      return execute_status;
    }

    const auto output_summary = output_domain.summary();
    auto output_storage = std::make_shared<Domain>(std::move(output_domain));
    return core::detail::store_mesh(
      model_handle,
      algorithm->name(),
      output_summary,
      std::move(output_storage),
      mesh_handle,
      resolved_context
    );
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
      "An unknown mesh SDK error occurred."
    );
  }
}

} // namespace

ParameterValue::ParameterValue(std::int64_t value) : storage_(value)
{
}

ParameterValue::ParameterValue(double value) : storage_(value)
{
}

ParameterValue::ParameterValue(bool value) : storage_(value)
{
}

ParameterValue::ParameterValue(std::string value) : storage_(std::move(value))
{
}

ParameterValue::ParameterValue(const char *value)
  : storage_(value == nullptr ? std::string {} : std::string(value))
{
}

ParameterType ParameterValue::type() const noexcept
{
  switch(storage_.index()) {
  case 1U:
    return ParameterType::integer;
  case 2U:
    return ParameterType::number;
  case 3U:
    return ParameterType::boolean;
  case 4U:
    return ParameterType::text;
  default:
    return ParameterType::empty;
  }
}

bool ParameterValue::empty() const noexcept
{
  return std::holds_alternative<std::monostate>(storage_);
}

bool ParameterValue::is_integer() const noexcept
{
  return std::holds_alternative<std::int64_t>(storage_);
}

bool ParameterValue::is_number() const noexcept
{
  return std::holds_alternative<double>(storage_) ||
         std::holds_alternative<std::int64_t>(storage_);
}

bool ParameterValue::is_boolean() const noexcept
{
  return std::holds_alternative<bool>(storage_);
}

bool ParameterValue::is_text() const noexcept
{
  return std::holds_alternative<std::string>(storage_);
}

std::int64_t ParameterValue::integer(std::int64_t fallback) const noexcept
{
  if(const auto *value = std::get_if<std::int64_t>(&storage_)) {
    return *value;
  }
  return fallback;
}

double ParameterValue::number(double fallback) const noexcept
{
  if(const auto *value = std::get_if<double>(&storage_)) {
    return *value;
  }
  if(const auto *integer_value = std::get_if<std::int64_t>(&storage_)) {
    return static_cast<double>(*integer_value);
  }
  return fallback;
}

bool ParameterValue::boolean(bool fallback) const noexcept
{
  if(const auto *value = std::get_if<bool>(&storage_)) {
    return *value;
  }
  return fallback;
}

std::string_view ParameterValue::text() const noexcept
{
  if(const auto *value = std::get_if<std::string>(&storage_)) {
    return *value;
  }
  return {};
}

void ParameterDictionary::set(std::string key, ParameterValue value)
{
  entries_.insert_or_assign(std::move(key), std::move(value));
}

void ParameterDictionary::set_integer(std::string key, std::int64_t value)
{
  set(std::move(key), ParameterValue(value));
}

void ParameterDictionary::set_number(std::string key, double value)
{
  set(std::move(key), ParameterValue(value));
}

void ParameterDictionary::set_boolean(std::string key, bool value)
{
  set(std::move(key), ParameterValue(value));
}

void ParameterDictionary::set_text(std::string key, std::string value)
{
  set(std::move(key), ParameterValue(std::move(value)));
}

bool ParameterDictionary::contains(std::string_view key) const
{
  return entries_.find(key) != entries_.end();
}

const ParameterValue *ParameterDictionary::find(std::string_view key) const noexcept
{
  const auto entry_it = entries_.find(key);
  if(entry_it == entries_.end()) {
    return nullptr;
  }
  return &entry_it->second;
}

bool ParameterDictionary::try_get_integer(
  std::string_view key,
  std::int64_t &value
) const noexcept
{
  const auto *entry = find(key);
  if(entry == nullptr || !entry->is_integer()) {
    return false;
  }
  value = entry->integer();
  return true;
}

bool ParameterDictionary::try_get_number(
  std::string_view key,
  double &value
) const noexcept
{
  const auto *entry = find(key);
  if(entry == nullptr || !entry->is_number()) {
    return false;
  }
  value = entry->number();
  return true;
}

bool ParameterDictionary::try_get_boolean(
  std::string_view key,
  bool &value
) const noexcept
{
  const auto *entry = find(key);
  if(entry == nullptr || !entry->is_boolean()) {
    return false;
  }
  value = entry->boolean();
  return true;
}

bool ParameterDictionary::try_get_text(
  std::string_view key,
  std::string_view &value
) const noexcept
{
  const auto *entry = find(key);
  if(entry == nullptr || !entry->is_text()) {
    return false;
  }
  value = entry->text();
  return true;
}

std::size_t ParameterDictionary::size() const noexcept
{
  return entries_.size();
}

bool is_algorithm_registered(std::string_view algorithm_name) noexcept
{
  detail::register_builtin_meshing_algorithms();
  static_cast<void>(core::detail::clear_error_state());
  return detail::is_meshing_algorithm_registered(algorithm_name);
}

void apply_auto_cfd_spacing(
  MeshingOptions &options,
  const AutoCfdSpacingOptions &spacing
) noexcept
{
  options.parameters.set_number("spacing_growth_rate", spacing.growth_rate);
  options.parameters.set_number("spacing_feature_angle", spacing.feature_angle);
  options.parameters.set_number("spacing_minimum_length", spacing.minimum_length);
  options.parameters.set_number("spacing_maximum_length", spacing.maximum_length);
  options.parameters.set_number(
    "spacing_sharp_angle_limit",
    spacing.sharp_angle_limit
  );
  options.parameters.set_number(
    "spacing_sharp_angle_length",
    spacing.sharp_angle_length
  );
  set_auto_cfd_proximity_parameters(
    options.parameters,
    "spacing_",
    spacing.self_proximity,
    spacing.pid_proximity,
    spacing.maximum_normals_angle,
    spacing.length_to_gap_ratio,
    spacing.proximity_minimum_length
  );
}

void apply_auto_cfd_surface_defaults(
  MeshingOptions &options,
  const AutoCfdSurfaceDefaults &defaults
) noexcept
{
  options.parameters.set_text("element_type", defaults.element_type);
  options.parameters.set_number("growth_rate", defaults.growth_rate);
  options.parameters.set_number("distortion_angle", defaults.distortion_angle);
  options.parameters.set_number("minimum_length", defaults.minimum_length);
  options.parameters.set_number("maximum_length", defaults.maximum_length);
  set_auto_cfd_proximity_parameters(
    options.parameters,
    "",
    defaults.self_proximity,
    defaults.pid_proximity,
    defaults.maximum_normals_angle,
    defaults.length_to_gap_ratio,
    defaults.proximity_minimum_length
  );
}

base::StatusCode create_surface_mesh(
  geo::ModelHandle model_handle,
  std::string_view algorithm_name,
  const MeshingOptions &options,
  MeshHandle &mesh_handle,
  base::ContextHandle context_handle
) noexcept
{
  return create_mesh_impl(
    model_handle,
    algorithm_name,
    options,
    detail::MeshingDimension::surface,
    mesh_handle,
    context_handle
  );
}

base::StatusCode create_surface_mesh(
  geo::ModelHandle model_handle,
  std::string_view algorithm_name,
  const ParameterDictionary &parameters,
  MeshHandle &mesh_handle,
  base::ContextHandle context_handle
) noexcept
{
  MeshingOptions options;
  options.parameters = parameters;
  return create_surface_mesh(
    model_handle,
    algorithm_name,
    options,
    mesh_handle,
    context_handle
  );
}

base::StatusCode create_volume_mesh(
  geo::ModelHandle model_handle,
  std::string_view algorithm_name,
  const MeshingOptions &options,
  MeshHandle &mesh_handle,
  base::ContextHandle context_handle
) noexcept
{
  return create_mesh_impl(
    model_handle,
    algorithm_name,
    options,
    detail::MeshingDimension::volume,
    mesh_handle,
    context_handle
  );
}

base::StatusCode create_volume_mesh(
  geo::ModelHandle model_handle,
  std::string_view algorithm_name,
  const ParameterDictionary &parameters,
  MeshHandle &mesh_handle,
  base::ContextHandle context_handle
) noexcept
{
  MeshingOptions options;
  options.parameters = parameters;
  return create_volume_mesh(
    model_handle,
    algorithm_name,
    options,
    mesh_handle,
    context_handle
  );
}

base::StatusCode mesh_summary(
  MeshHandle mesh_handle,
  MeshSummary &summary,
  base::ContextHandle context_handle
) noexcept
{
  return core::detail::mesh_summary(mesh_handle, summary, context_handle);
}

base::StatusCode mesh_quality_report(
  MeshHandle mesh_handle,
  MeshQualityReport &report,
  base::ContextHandle context_handle
) noexcept
{
  try {
    report = {};
    return core::detail::with_mesh_domain(
      mesh_handle,
      [&](const Domain &domain) {
        report = domain.quality_report();
        return base::StatusCode::ok;
      },
      context_handle
    );
  }
  catch(const std::exception &exception) {
    report = {};
    return core::detail::publish_error(
      base::StatusCode::internal_error,
      exception.what()
    );
  }
  catch(...) {
    report = {};
    return core::detail::publish_error(
      base::StatusCode::internal_error,
      "An unknown mesh quality reporting error occurred."
    );
  }
}

base::StatusCode domain_snapshot(
  MeshHandle mesh_handle,
  Domain &domain,
  base::ContextHandle context_handle
) noexcept
{
  domain = Domain("");
  return core::detail::mesh_domain_snapshot(mesh_handle, domain, context_handle);
}

base::StatusCode nodes_count(
  MeshHandle mesh_handle,
  std::size_t &count,
  base::ContextHandle context_handle
) noexcept
{
  return core::detail::mesh_node_count(mesh_handle, count, context_handle);
}

base::StatusCode cells_count(
  MeshHandle mesh_handle,
  std::size_t &count,
  base::ContextHandle context_handle
) noexcept
{
  return core::detail::mesh_cell_count(mesh_handle, count, context_handle);
}

} // namespace sqmesh::mesh
