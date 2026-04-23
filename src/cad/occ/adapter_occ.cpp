// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#include "adapter.hpp"

#include "../../core/log.hpp"
#include "../../core/runtime_registry.hpp"
#include "../../model/geometry_model_storage.hpp"

#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepLProp_SLProps.hxx>
#include <BRep_Tool.hxx>
#include <BRepTools.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <Bnd_Box.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom_Surface.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <IGESControl_Reader.hxx>
#include <IGESControl_Writer.hxx>
#include <Interface_Static.hxx>
#include <STEPControl_Reader.hxx>
#include <STEPControl_StepModelType.hxx>
#include <STEPControl_Writer.hxx>
#include <GProp_GProps.hxx>
#include <ShapeExtend_Status.hxx>
#include <ShapeFix_Shape.hxx>
#include <Standard_Failure.hxx>
#include <Standard_Version.hxx>
#include <Poly_PolygonOnTriangulation.hxx>
#include <TopAbs_Orientation.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopAbs_State.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Pnt2d.hxx>

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace sqmesh::cad::occ {
namespace {

using sqmesh::model::detail::GeometryKernel;
using sqmesh::model::detail::GeometryCoarseProxyMesh;
using sqmesh::model::detail::GeometryModelStorage;
using sqmesh::model::detail::GeometryModelStoragePtr;

using TopologyAdjacencyLists = std::vector<std::vector<geo::TopologyEntityId>>;
using TopologyShapeLists = std::array<std::vector<TopoDS_Shape>, 4>;
constexpr double kPi = 3.14159265358979323846;

struct TopologyCache final {
  geo::TopologySnapshot snapshot {};
  std::array<TopologyAdjacencyLists, 4> children {};
  std::array<TopologyAdjacencyLists, 4> parents {};
  TopologyShapeLists shapes {};
  std::vector<geo::FaceBoundaryLoops> face_boundaries {};
  std::vector<std::uint8_t> face_boundary_ready {};
};

struct ShapeIndexData final {
  TopTools_IndexedMapOfShape normalized_map {};
  std::vector<TopoDS_Shape> representative_shapes {};
};

class OccModelStorage final : public GeometryModelStorage {
public:
  explicit OccModelStorage(TopoDS_Shape shape) : shape_(std::move(shape))
  {
  }

  [[nodiscard]] GeometryKernel kernel() const noexcept override
  {
    return GeometryKernel::occ;
  }

  [[nodiscard]] base::StatusCode coarse_proxy_mesh(
    GeometryCoarseProxyMesh &proxy_mesh
  ) const noexcept override;
  [[nodiscard]] base::StatusCode copy_shape(TopoDS_Shape &shape) const noexcept;

  [[nodiscard]] std::uint64_t topology_revision() const noexcept
  {
    return topology_revision_;
  }

  void set_shape(TopoDS_Shape shape) noexcept
  {
    std::lock_guard<std::mutex> lock(topology_mutex_);
    shape_ = std::move(shape);
    ++topology_revision_;
    topology_cache_ = {};
    proxy_cache_ = {};
    topology_cache_valid_ = false;
    proxy_cache_valid_ = false;
    surface_cache_valid_ = false;
  }

  [[nodiscard]] base::StatusCode topology_snapshot(
    geo::TopologySnapshot &snapshot
  ) const noexcept override;
  [[nodiscard]] base::StatusCode topology_children(
    geo::TopologyEntityId entity,
    std::vector<geo::TopologyEntityId> &children
  ) const noexcept override;
  [[nodiscard]] base::StatusCode topology_parents(
    geo::TopologyEntityId entity,
    std::vector<geo::TopologyEntityId> &parents
  ) const noexcept override;
  [[nodiscard]] base::StatusCode face_uv_bounds(
    geo::TopologyEntityId face_entity,
    geo::FaceUvBounds &bounds
  ) const noexcept override;
  [[nodiscard]] base::StatusCode sample_face(
    geo::TopologyEntityId face_entity,
    double u,
    double v,
    geo::FaceSample &sample
  ) const noexcept override;
  [[nodiscard]] base::StatusCode sample_face_curvature(
    geo::TopologyEntityId face_entity,
    double u,
    double v,
    geo::FaceCurvatureSample &sample
  ) const noexcept override;
  [[nodiscard]] base::StatusCode sample_face_derivatives(
    geo::TopologyEntityId face_entity,
    double u,
    double v,
    geo::FaceDerivatives &sample
  ) const noexcept override;
  [[nodiscard]] base::StatusCode project_point_to_face(
    geo::TopologyEntityId face_entity,
    const geo::Point3 &point,
    geo::FaceProjection &projection
  ) const noexcept override;
  [[nodiscard]] base::StatusCode recover_face_uv(
    geo::TopologyEntityId face_entity,
    const geo::Point3 &point,
    geo::FaceUvMapping &mapping
  ) const noexcept override;
  [[nodiscard]] base::StatusCode recover_face_uv_from_edge(
    geo::TopologyEntityId face_entity,
    geo::TopologyEntityId edge_entity,
    double edge_parameter,
    geo::FaceUvMapping &mapping
  ) const noexcept override;
  [[nodiscard]] base::StatusCode recover_face_uv_from_edge_use(
    geo::TopologyEntityId face_entity,
    const geo::FaceBoundaryEdgeUse &edge_use,
    double edge_parameter,
    geo::FaceUvMapping &mapping
  ) const noexcept override;
  [[nodiscard]] base::StatusCode edge_curve_info(
    geo::TopologyEntityId edge_entity,
    geo::EdgeCurveInfo &info
  ) const noexcept override;
  [[nodiscard]] base::StatusCode sample_edge_tangent(
    geo::TopologyEntityId edge_entity,
    double parameter,
    geo::EdgeTangentSample &sample
  ) const noexcept override;
  [[nodiscard]] base::StatusCode face_boundary_loops(
    geo::TopologyEntityId face_entity,
    geo::FaceBoundaryLoops &boundary
  ) const noexcept override;
  [[nodiscard]] base::StatusCode feature_edges(
    geo::FeatureEdgeReport &report,
    const geo::FeatureEdgeOptions &options
  ) const noexcept override;

private:
  void ensure_topology_cache_locked() const;
  [[nodiscard]] base::StatusCode ensure_proxy_cache_locked() const;

  // Cached BRepAdaptor_Surface to avoid per-call reconstruction.
  // Must be called while topology_mutex_ is held and topology_cache_ is valid.
  [[nodiscard]] const BRepAdaptor_Surface &cached_surface_locked(
    const TopoDS_Face &face,
    std::uint32_t face_index
  ) const {
    if(surface_cache_valid_ && surface_cache_face_ == face_index) {
      return surface_cache_;
    }
    surface_cache_ = BRepAdaptor_Surface(face, Standard_True);
    surface_cache_face_ = face_index;
    surface_cache_valid_ = true;
    return surface_cache_;
  }

  TopoDS_Shape shape_ {};
  mutable std::mutex topology_mutex_ {};
  mutable TopologyCache topology_cache_ {};
  mutable GeometryCoarseProxyMesh proxy_cache_ {};
  mutable BRepAdaptor_Surface surface_cache_ {};
  mutable std::uint32_t surface_cache_face_ = std::numeric_limits<std::uint32_t>::max();
  mutable bool topology_cache_valid_ = false;
  mutable bool proxy_cache_valid_ = false;
  mutable bool surface_cache_valid_ = false;
  std::uint64_t topology_revision_ = 1U;
};

class ScopedInterfaceCString final {
public:
  ScopedInterfaceCString(const char *name, std::string_view value) : name_(name)
  {
    if(value.empty()) {
      return;
    }

    const auto *current = Interface_Static::CVal(name_);
    previous_value_ = current != nullptr ? current : "";
    const std::string next_value(value);
    if(Interface_Static::SetCVal(name_, next_value.c_str())) {
      restore_previous_ = true;
    }
  }

  ~ScopedInterfaceCString()
  {
    if(restore_previous_) {
      static_cast<void>(Interface_Static::SetCVal(name_, previous_value_.c_str()));
    }
  }

private:
  const char *name_ = "";
  std::string previous_value_ {};
  bool restore_previous_ = false;
};

class ScopedInterfaceInteger final {
public:
  ScopedInterfaceInteger(const char *name, Standard_Integer value) : name_(name)
  {
    previous_value_ = Interface_Static::IVal(name_);
    if(Interface_Static::SetIVal(name_, value)) {
      restore_previous_ = true;
    }
  }

  ~ScopedInterfaceInteger()
  {
    if(restore_previous_) {
      static_cast<void>(Interface_Static::SetIVal(name_, previous_value_));
    }
  }

private:
  const char *name_ = "";
  Standard_Integer previous_value_ = 0;
  bool restore_previous_ = false;
};

[[nodiscard]] std::size_t count_shapes(
  const TopoDS_Shape &shape,
  TopAbs_ShapeEnum shape_kind
)
{
  TopTools_IndexedMapOfShape indexed_shapes;
  TopExp::MapShapes(shape, shape_kind, indexed_shapes);
  return static_cast<std::size_t>(indexed_shapes.Extent());
}

[[nodiscard]] constexpr std::size_t topology_dimension_index(
  geo::TopologyDimension dimension
) noexcept
{
  return static_cast<std::size_t>(dimension);
}

[[nodiscard]] TopoDS_Shape normalized_shape(const TopoDS_Shape &shape)
{
  return shape.Oriented(TopAbs_FORWARD);
}

[[nodiscard]] bool same_normalized_shape(
  const TopoDS_Shape &lhs,
  const TopoDS_Shape &rhs
) noexcept
{
  if(lhs.IsNull() || rhs.IsNull()) {
    return false;
  }

  return normalized_shape(lhs).IsSame(normalized_shape(rhs));
}

[[nodiscard]] ShapeIndexData collect_shape_index_data(
  const TopoDS_Shape &shape,
  TopAbs_ShapeEnum shape_kind
)
{
  ShapeIndexData data;
  for(TopExp_Explorer explorer(shape, shape_kind); explorer.More(); explorer.Next()) {
    const TopoDS_Shape original_shape = explorer.Current();
    const TopoDS_Shape normalized = normalized_shape(original_shape);
    if(data.normalized_map.FindIndex(normalized) > 0) {
      continue;
    }

    data.normalized_map.Add(normalized);
    data.representative_shapes.push_back(original_shape);
  }
  return data;
}

void initialize_topology_adjacency(
  TopologyCache &cache,
  const ShapeIndexData &regions,
  const ShapeIndexData &faces,
  const ShapeIndexData &edges,
  const ShapeIndexData &vertices
)
{
  cache.children[topology_dimension_index(geo::TopologyDimension::region)]
    .resize(regions.representative_shapes.size());
  cache.children[topology_dimension_index(geo::TopologyDimension::face)]
    .resize(faces.representative_shapes.size());
  cache.children[topology_dimension_index(geo::TopologyDimension::edge)]
    .resize(edges.representative_shapes.size());
  cache.children[topology_dimension_index(geo::TopologyDimension::vertex)]
    .resize(vertices.representative_shapes.size());

  cache.parents[topology_dimension_index(geo::TopologyDimension::region)]
    .resize(regions.representative_shapes.size());
  cache.parents[topology_dimension_index(geo::TopologyDimension::face)]
    .resize(faces.representative_shapes.size());
  cache.parents[topology_dimension_index(geo::TopologyDimension::edge)]
    .resize(edges.representative_shapes.size());
  cache.parents[topology_dimension_index(geo::TopologyDimension::vertex)]
    .resize(vertices.representative_shapes.size());

  cache.shapes[topology_dimension_index(geo::TopologyDimension::region)] =
    regions.representative_shapes;
  cache.shapes[topology_dimension_index(geo::TopologyDimension::face)] =
    faces.representative_shapes;
  cache.shapes[topology_dimension_index(geo::TopologyDimension::edge)] =
    edges.representative_shapes;
  cache.shapes[topology_dimension_index(geo::TopologyDimension::vertex)] =
    vertices.representative_shapes;
}

void add_topology_relation(
  TopologyCache &cache,
  geo::TopologyDimension parent_dimension,
  std::uint32_t parent_index,
  geo::TopologyDimension child_dimension,
  std::uint32_t child_index
)
{
  cache.children[topology_dimension_index(parent_dimension)][parent_index].push_back(
    {child_dimension, child_index}
  );
  cache.parents[topology_dimension_index(child_dimension)][child_index].push_back(
    {parent_dimension, parent_index}
  );
}

void build_child_relations(
  const std::vector<TopoDS_Shape> &parents,
  TopAbs_ShapeEnum child_shape_kind,
  const TopTools_IndexedMapOfShape &global_children,
  geo::TopologyDimension parent_dimension,
  geo::TopologyDimension child_dimension,
  TopologyCache &cache
)
{
  for(std::size_t parent_index = 0; parent_index < parents.size(); ++parent_index) {
    TopTools_IndexedMapOfShape child_shapes;
    TopExp::MapShapes(parents[parent_index], child_shape_kind, child_shapes);
    for(Standard_Integer child_offset = 1; child_offset <= child_shapes.Extent(); ++child_offset) {
      const auto global_child_index =
        global_children.FindIndex(normalized_shape(child_shapes.FindKey(child_offset)));
      if(global_child_index <= 0) {
        continue;
      }

      add_topology_relation(
        cache,
        parent_dimension,
        static_cast<std::uint32_t>(parent_index),
        child_dimension,
        static_cast<std::uint32_t>(global_child_index - 1)
      );
    }
  }
}

[[nodiscard]] std::vector<geo::TopologyEntityInfo> build_entity_infos(
  geo::TopologyDimension dimension,
  const TopologyAdjacencyLists &parents,
  const TopologyAdjacencyLists &children
)
{
  std::vector<geo::TopologyEntityInfo> entities;
  entities.reserve(children.size());

  for(std::size_t index = 0; index < children.size(); ++index) {
    entities.push_back(
      {
        {dimension, static_cast<std::uint32_t>(index)},
        parents[index].size(),
        children[index].size(),
      }
    );
  }

  return entities;
}

[[nodiscard]] TopologyCache build_topology_cache(
  const TopoDS_Shape &shape,
  std::uint64_t topology_revision
)
{
  TopologyCache cache;
  cache.snapshot.topology_revision = topology_revision;
  if(shape.IsNull()) {
    return cache;
  }

  const auto regions = collect_shape_index_data(shape, TopAbs_SOLID);
  const auto faces = collect_shape_index_data(shape, TopAbs_FACE);
  const auto edges = collect_shape_index_data(shape, TopAbs_EDGE);
  const auto vertices = collect_shape_index_data(shape, TopAbs_VERTEX);

  initialize_topology_adjacency(cache, regions, faces, edges, vertices);

  build_child_relations(
    regions.representative_shapes,
    TopAbs_FACE,
    faces.normalized_map,
    geo::TopologyDimension::region,
    geo::TopologyDimension::face,
    cache
  );
  build_child_relations(
    faces.representative_shapes,
    TopAbs_EDGE,
    edges.normalized_map,
    geo::TopologyDimension::face,
    geo::TopologyDimension::edge,
    cache
  );
  build_child_relations(
    edges.representative_shapes,
    TopAbs_VERTEX,
    vertices.normalized_map,
    geo::TopologyDimension::edge,
    geo::TopologyDimension::vertex,
    cache
  );

  cache.snapshot.regions = build_entity_infos(
    geo::TopologyDimension::region,
    cache.parents[topology_dimension_index(geo::TopologyDimension::region)],
    cache.children[topology_dimension_index(geo::TopologyDimension::region)]
  );
  cache.snapshot.faces = build_entity_infos(
    geo::TopologyDimension::face,
    cache.parents[topology_dimension_index(geo::TopologyDimension::face)],
    cache.children[topology_dimension_index(geo::TopologyDimension::face)]
  );
  cache.snapshot.edges = build_entity_infos(
    geo::TopologyDimension::edge,
    cache.parents[topology_dimension_index(geo::TopologyDimension::edge)],
    cache.children[topology_dimension_index(geo::TopologyDimension::edge)]
  );
  cache.snapshot.vertices = build_entity_infos(
    geo::TopologyDimension::vertex,
    cache.parents[topology_dimension_index(geo::TopologyDimension::vertex)],
    cache.children[topology_dimension_index(geo::TopologyDimension::vertex)]
  );
  cache.face_boundaries.resize(cache.snapshot.faces.size());
  cache.face_boundary_ready.assign(cache.snapshot.faces.size(), 0U);
  return cache;
}

[[nodiscard]] base::StatusCode internal_error(std::string_view message) noexcept;
[[nodiscard]] base::StatusCode publish_occ_failure(
  std::string_view fallback_message,
  const Standard_Failure &failure
) noexcept;
[[nodiscard]] base::StatusCode publish_std_failure(
  std::string_view fallback_message,
  const std::exception &exception
) noexcept;

[[nodiscard]] base::StatusCode copy_topology_neighbors(
  const std::array<TopologyAdjacencyLists, 4> &neighbors,
  geo::TopologyEntityId entity,
  std::vector<geo::TopologyEntityId> &output,
  std::string_view relation_name
) noexcept
{
  const auto &entities = neighbors[topology_dimension_index(entity.dimension)];
  if(entity.index >= entities.size()) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      relation_name == "child"
        ? "Requested topology child lookup references a non-existent entity."
        : "Requested topology parent lookup references a non-existent entity."
    );
  }

  output = entities[entity.index];
  return core::detail::clear_error_state();
}

[[nodiscard]] const TopoDS_Shape *lookup_topology_shape(
  const TopologyCache &cache,
  geo::TopologyEntityId entity
) noexcept
{
  const auto &shapes = cache.shapes[topology_dimension_index(entity.dimension)];
  if(entity.index >= shapes.size()) {
    return nullptr;
  }
  return &shapes[entity.index];
}

[[nodiscard]] geo::TopologyEntityId lookup_topology_entity_id(
  const TopologyCache &cache,
  geo::TopologyDimension dimension,
  const TopoDS_Shape &shape
) noexcept
{
  if(shape.IsNull()) {
    return {dimension, geo::invalid_topology_index};
  }

  const auto &shapes = cache.shapes[topology_dimension_index(dimension)];
  for(std::size_t index = 0; index < shapes.size(); ++index) {
    if(same_normalized_shape(shapes[index], shape)) {
      return {dimension, static_cast<std::uint32_t>(index)};
    }
  }

  return {dimension, geo::invalid_topology_index};
}

[[nodiscard]] geo::Point3 to_point3(const gp_Pnt &point) noexcept
{
  return {point.X(), point.Y(), point.Z()};
}

[[nodiscard]] geo::Vector3 to_vector3(const gp_Vec &vector) noexcept
{
  return {vector.X(), vector.Y(), vector.Z()};
}

[[nodiscard]] gp_Pnt to_gp_point(const geo::Point3 &point) noexcept
{
  return {point[0], point[1], point[2]};
}

[[nodiscard]] geo::Point3 transformed_point3(
  const gp_Pnt &point,
  const TopLoc_Location &location
) noexcept
{
  if(location.IsIdentity()) {
    return to_point3(point);
  }
  return to_point3(point.Transformed(location.Transformation()));
}

[[nodiscard]] geo::Point3 triangulation_node_point(
  const Handle(Poly_Triangulation) &triangulation,
  const TopLoc_Location &location,
  Standard_Integer node_index
) noexcept
{
  return transformed_point3(triangulation->Node(node_index), location);
}

[[nodiscard]] double point_distance_squared(
  const geo::Point3 &lhs,
  const geo::Point3 &rhs
) noexcept
{
  const double dx = lhs[0] - rhs[0];
  const double dy = lhs[1] - rhs[1];
  const double dz = lhs[2] - rhs[2];
  return dx * dx + dy * dy + dz * dz;
}

[[nodiscard]] double shape_diagonal_length(const TopoDS_Shape &shape)
{
  if(shape.IsNull()) {
    return 1.0;
  }

  Bnd_Box bounds;
  BRepBndLib::Add(shape, bounds, Standard_False);
  if(bounds.IsVoid()) {
    return 1.0;
  }

  Standard_Real x_min = 0.0;
  Standard_Real y_min = 0.0;
  Standard_Real z_min = 0.0;
  Standard_Real x_max = 0.0;
  Standard_Real y_max = 0.0;
  Standard_Real z_max = 0.0;
  bounds.Get(x_min, y_min, z_min, x_max, y_max, z_max);
  const double dx = x_max - x_min;
  const double dy = y_max - y_min;
  const double dz = z_max - z_min;
  const double diagonal = std::sqrt(dx * dx + dy * dy + dz * dz);
  return diagonal > 0.0 ? diagonal : 1.0;
}

void ensure_shape_face_triangulations(const TopoDS_Shape &shape)
{
  bool requires_triangulation = false;
  for(TopExp_Explorer explorer(shape, TopAbs_FACE); explorer.More(); explorer.Next()) {
    TopLoc_Location location;
    const Handle(Poly_Triangulation) triangulation =
      BRep_Tool::Triangulation(TopoDS::Face(explorer.Current()), location);
    if(triangulation.IsNull() || triangulation->NbTriangles() <= 0) {
      requires_triangulation = true;
      break;
    }
  }

  if(!requires_triangulation) {
    return;
  }

  constexpr double kStlLinearDeflection = 0.001;
  constexpr double kStlAngularDeflection = 0.3; // ~17 degrees
  BRepMesh_IncrementalMesh mesher(
    shape,
    kStlLinearDeflection,
    Standard_True,       // Relative
    kStlAngularDeflection,
    Standard_True        // InParallel
  );
  static_cast<void>(mesher);
}


[[nodiscard]] geo::TopologyEntityId more_specific_owner(
  geo::TopologyEntityId current,
  geo::TopologyEntityId candidate
) noexcept
{
  if(!geo::is_valid(candidate)) {
    return current;
  }
  if(!geo::is_valid(current)) {
    return candidate;
  }

  return topology_dimension_index(candidate.dimension) <
           topology_dimension_index(current.dimension)
         ? candidate
         : current;
}

[[nodiscard]] std::array<geo::TopologyEntityId, 2> edge_endpoint_entities(
  const TopologyCache &cache,
  geo::TopologyEntityId edge_entity
) noexcept
{
  std::array<geo::TopologyEntityId, 2> endpoints {};
  endpoints[0].index = geo::invalid_topology_index;
  endpoints[1].index = geo::invalid_topology_index;

  if(edge_entity.dimension != geo::TopologyDimension::edge) {
    return endpoints;
  }

  const auto *edge_shape = lookup_topology_shape(cache, edge_entity);
  if(edge_shape == nullptr) {
    return endpoints;
  }

  TopoDS_Vertex start_vertex_shape;
  TopoDS_Vertex end_vertex_shape;
  TopExp::Vertices(TopoDS::Edge(*edge_shape), start_vertex_shape, end_vertex_shape, true);
  endpoints[0] = lookup_topology_entity_id(
    cache,
    geo::TopologyDimension::vertex,
    start_vertex_shape
  );
  endpoints[1] = lookup_topology_entity_id(
    cache,
    geo::TopologyDimension::vertex,
    end_vertex_shape
  );
  return endpoints;
}

[[nodiscard]] base::StatusCode build_proxy_cache(
  const TopoDS_Shape &shape,
  const TopologyCache &topology_cache,
  std::uint64_t topology_revision,
  GeometryCoarseProxyMesh &proxy_cache
)
{
  proxy_cache = {};
  proxy_cache.topology_revision = topology_revision;
  proxy_cache.face_triangle_ranges.resize(topology_cache.snapshot.faces.size());
  if(shape.IsNull() || topology_cache.snapshot.faces.empty()) {
    return core::detail::clear_error_state();
  }

 
  ensure_shape_face_triangulations(shape);

  const double kMergeTol = std::max(shape_diagonal_length(shape) * 1.0e-9, 1.0e-12);
  const double kMergeTol2 = kMergeTol * kMergeTol;

  using BucketKey = std::tuple<long long, long long, long long>;
  std::map<BucketKey, std::vector<std::uint32_t>> xyz_bucket;

  auto bucket_key = [kMergeTol](double x, double y, double z) -> BucketKey {
    return {
      static_cast<long long>(std::floor(x / kMergeTol)),
      static_cast<long long>(std::floor(y / kMergeTol)),
      static_cast<long long>(std::floor(z / kMergeTol)),
    };
  };

  auto get_or_add_vertex = [&](
    double x, double y, double z,
    const std::array<double, 2> &uv,
    geo::TopologyEntityId uv_face,
    geo::TopologyEntityId owner
  ) -> std::uint32_t {
    auto [kx, ky, kz] = bucket_key(x, y, z);
    for(long long dx = -1; dx <= 1; ++dx) {
      for(long long dy = -1; dy <= 1; ++dy) {
        for(long long dz = -1; dz <= 1; ++dz) {
          auto it = xyz_bucket.find({kx + dx, ky + dy, kz + dz});
          if(it == xyz_bucket.end()) {
            continue;
          }
          for(const std::uint32_t id : it->second) {
            const double ddx = proxy_cache.nodes[id][0] - x;
            const double ddy = proxy_cache.nodes[id][1] - y;
            const double ddz = proxy_cache.nodes[id][2] - z;
            if(ddx * ddx + ddy * ddy + ddz * ddz < kMergeTol2) {
              proxy_cache.node_topology_owner[id] =
                more_specific_owner(proxy_cache.node_topology_owner[id], owner);
              return id;
            }
          }
        }
      }
    }

    const auto id = static_cast<std::uint32_t>(proxy_cache.nodes.size());
    proxy_cache.nodes.push_back({x, y, z});
    proxy_cache.node_uv.push_back(uv);
    proxy_cache.node_uv_face.push_back(uv_face);
    proxy_cache.node_topology_owner.push_back(owner);
    xyz_bucket[{kx, ky, kz}].push_back(id);
    return id;
  };

  // ---- Phase 1: Face triangulations -> global vertex pool ----
  struct FaceLocalMapping {
    geo::TopologyEntityId face_entity {};
    Handle(Poly_Triangulation) triangulation {};
    TopLoc_Location location {};
    bool reverse_winding = false;
    std::vector<std::uint32_t> local_to_global {};
  };

  std::vector<FaceLocalMapping> face_mappings;
  face_mappings.reserve(topology_cache.snapshot.faces.size());

  for(std::size_t fi = 0U; fi < topology_cache.snapshot.faces.size(); ++fi) {
    const geo::TopologyEntityId face_entity {
      geo::TopologyDimension::face, static_cast<std::uint32_t>(fi),
    };
    const auto *shape_entry = lookup_topology_shape(topology_cache, face_entity);
    if(shape_entry == nullptr) {
      return core::detail::publish_error(
        base::StatusCode::internal_error,
        "Proxy triangulation cache build failed because a face representative is missing from the topology cache."
      );
    }

    const TopoDS_Face face = TopoDS::Face(*shape_entry);
    TopLoc_Location location;
    const Handle(Poly_Triangulation) triangulation =
      BRep_Tool::Triangulation(face, location);
    if(triangulation.IsNull() || triangulation->NbTriangles() <= 0) {
      continue;
    }

    FaceLocalMapping fm;
    fm.face_entity = face_entity;
    fm.triangulation = triangulation;
    fm.location = location;
    fm.reverse_winding = face.Orientation() == TopAbs_REVERSED;
    fm.local_to_global.resize(
      static_cast<std::size_t>(triangulation->NbNodes()) + 1U
    );

    const bool has_uv = triangulation->HasUVNodes();
    proxy_cache.source_local_node_count +=
      static_cast<std::size_t>(triangulation->NbNodes());

    for(Standard_Integer ni = 1; ni <= triangulation->NbNodes(); ++ni) {
      const gp_Pnt pp = triangulation->Node(ni);
      double x = pp.X();
      double y = pp.Y();
      double z = pp.Z();
      if(!location.IsIdentity()) {
        location.Transformation().Transforms(x, y, z);
      }

      std::array<double, 2> uv {0.0, 0.0};
      if(has_uv) {
        const gp_Pnt2d uv_pt = triangulation->UVNode(ni);
        uv = {uv_pt.X(), uv_pt.Y()};
      }

      fm.local_to_global[static_cast<std::size_t>(ni)] =
        get_or_add_vertex(x, y, z, uv, face_entity, face_entity);
    }

    face_mappings.push_back(std::move(fm));
  }

  // ---- Phase 2: Upgrade topology ownership along edges ----
  for(const auto &fm : face_mappings) {
    const auto *shape_entry = lookup_topology_shape(topology_cache, fm.face_entity);
    const TopoDS_Face face = TopoDS::Face(*shape_entry);

    for(TopExp_Explorer wire_it(face, TopAbs_WIRE); wire_it.More(); wire_it.Next()) {
      BRepTools_WireExplorer explorer(TopoDS::Wire(wire_it.Current()), face);
      for(; explorer.More(); explorer.Next()) {
        const TopoDS_Edge current_edge = TopoDS::Edge(explorer.Current());
        if(BRep_Tool::Degenerated(current_edge)) {
          continue;
        }

        const auto edge_entity = lookup_topology_entity_id(
          topology_cache, geo::TopologyDimension::edge, current_edge
        );
        if(!geo::is_valid(edge_entity)) {
          continue;
        }

        const TopoDS_Edge oriented_edge =
          TopoDS::Edge(current_edge.Oriented(explorer.Orientation()));
        const Handle(Poly_PolygonOnTriangulation) polygon =
          BRep_Tool::PolygonOnTriangulation(
            oriented_edge, fm.triangulation, fm.location
          );
        if(polygon.IsNull() || polygon->NbNodes() < 2) {
          continue;
        }

        const auto endpoints = edge_endpoint_entities(topology_cache, edge_entity);

        for(Standard_Integer pi = 1; pi <= polygon->NbNodes(); ++pi) {
          const Standard_Integer local_idx = polygon->Nodes()(pi);
          if(local_idx < 1 ||
             static_cast<std::size_t>(local_idx) >= fm.local_to_global.size()) {
            continue;
          }

          const std::uint32_t gid =
            fm.local_to_global[static_cast<std::size_t>(local_idx)];

          geo::TopologyEntityId owner = edge_entity;
          if(pi == 1) {
            owner = more_specific_owner(owner, endpoints[0]);
          }
          else if(pi == polygon->NbNodes()) {
            owner = more_specific_owner(owner, endpoints[1]);
          }

          proxy_cache.node_topology_owner[gid] =
            more_specific_owner(proxy_cache.node_topology_owner[gid], owner);
        }
      }
    }
  }

  // ---- Phase 3: Build face triangles ----
  for(const auto &fm : face_mappings) {
    auto &range = proxy_cache.face_triangle_ranges[fm.face_entity.index];
    range.offset = proxy_cache.triangles.size();
    range.count = 0U;

    for(Standard_Integer ti = 1; ti <= fm.triangulation->NbTriangles(); ++ti) {
      Standard_Integer n0 = 0;
      Standard_Integer n1 = 0;
      Standard_Integer n2 = 0;
      fm.triangulation->Triangle(ti).Get(n0, n1, n2);
      if(fm.reverse_winding) {
        std::swap(n0, n1);
      }
      if(n0 < 1 || n1 < 1 || n2 < 1 ||
         static_cast<std::size_t>(n0) >= fm.local_to_global.size() ||
         static_cast<std::size_t>(n1) >= fm.local_to_global.size() ||
         static_cast<std::size_t>(n2) >= fm.local_to_global.size()) {
        return core::detail::publish_error(
          base::StatusCode::internal_error,
          "Proxy triangulation cache build failed because an OCC triangle references a non-existent local node."
        );
      }

      const std::array<std::uint32_t, 3> triangle {
        fm.local_to_global[static_cast<std::size_t>(n0)],
        fm.local_to_global[static_cast<std::size_t>(n1)],
        fm.local_to_global[static_cast<std::size_t>(n2)],
      };
      if(triangle[0] == triangle[1] || triangle[1] == triangle[2] ||
         triangle[2] == triangle[0]) {
        continue;
      }

      proxy_cache.triangles.push_back(triangle);
      proxy_cache.triangle_face_owner.push_back(fm.face_entity);
      ++range.count;
    }
  }

  // ---- Phase 4: Build edge_nodes ----
  // For each topological edge, find its polygon on the first incident face's
  // triangulation and map local indices to global IDs.
  std::vector<bool> edge_visited(topology_cache.snapshot.edges.size(), false);

  for(const auto &fm : face_mappings) {
    const auto *shape_entry = lookup_topology_shape(topology_cache, fm.face_entity);
    const TopoDS_Face face = TopoDS::Face(*shape_entry);

    for(TopExp_Explorer wire_it(face, TopAbs_WIRE); wire_it.More(); wire_it.Next()) {
      BRepTools_WireExplorer explorer(TopoDS::Wire(wire_it.Current()), face);
      for(; explorer.More(); explorer.Next()) {
        const TopoDS_Edge current_edge = TopoDS::Edge(explorer.Current());
        if(BRep_Tool::Degenerated(current_edge) ||
           BRep_Tool::IsClosed(current_edge, face)) {
          continue;
        }

        const auto edge_entity = lookup_topology_entity_id(
          topology_cache, geo::TopologyDimension::edge, current_edge
        );
        if(!geo::is_valid(edge_entity) ||
           edge_entity.index >= edge_visited.size()) {
          continue;
        }
        if(edge_visited[edge_entity.index]) {
          continue;
        }

        const TopoDS_Edge oriented_edge =
          TopoDS::Edge(current_edge.Oriented(explorer.Orientation()));
        const Handle(Poly_PolygonOnTriangulation) polygon =
          BRep_Tool::PolygonOnTriangulation(
            oriented_edge, fm.triangulation, fm.location
          );
        if(polygon.IsNull() || polygon->NbNodes() < 2) {
          continue;
        }

        // Build global node sequence for this edge.
        std::vector<std::uint32_t> global_nodes;
        global_nodes.reserve(static_cast<std::size_t>(polygon->NbNodes()));
        bool valid = true;
        for(Standard_Integer pi = 1; pi <= polygon->NbNodes(); ++pi) {
          const Standard_Integer local_idx = polygon->Nodes()(pi);
          if(local_idx < 1 ||
             static_cast<std::size_t>(local_idx) >= fm.local_to_global.size()) {
            valid = false;
            break;
          }
          global_nodes.push_back(
            fm.local_to_global[static_cast<std::size_t>(local_idx)]
          );
        }
        if(!valid || global_nodes.size() < 2U) {
          continue;
        }

        // Ensure canonical ordering: polygon start should match edge start vertex.
        const auto *rep_edge = lookup_topology_shape(topology_cache, edge_entity);
        if(rep_edge != nullptr && global_nodes.size() >= 2U) {
          TopoDS_Vertex v_start, v_end;
          TopExp::Vertices(
            TopoDS::Edge(*rep_edge), v_start, v_end, Standard_True
          );
          if(!v_start.IsNull() && !v_end.IsNull()) {
            const auto &first_pos = proxy_cache.nodes[global_nodes.front()];
            const auto &last_pos = proxy_cache.nodes[global_nodes.back()];
            const auto start_pt = to_point3(BRep_Tool::Pnt(v_start));
            const auto end_pt = to_point3(BRep_Tool::Pnt(v_end));
            const double direct_err =
              point_distance_squared(first_pos, start_pt) +
              point_distance_squared(last_pos, end_pt);
            const double reverse_err =
              point_distance_squared(first_pos, end_pt) +
              point_distance_squared(last_pos, start_pt);
            if(reverse_err < direct_err) {
              std::reverse(global_nodes.begin(), global_nodes.end());
            }
          }
        }

        proxy_cache.edge_nodes.push_back(std::move(global_nodes));
        proxy_cache.edge_topology_owner.push_back(edge_entity);
        edge_visited[edge_entity.index] = true;
      }
    }
  }

  return core::detail::clear_error_state();
}

struct FaceProjectionData final {
  gp_Pnt projected_point {};
  double u = 0.0;
  double v = 0.0;
  double distance = 0.0;
};

[[nodiscard]] bool oriented_face_normal(
  const TopoDS_Face &face,
  const BRepAdaptor_Surface &surface,
  double u,
  double v,
  gp_Vec &normal
)
{
  BRepLProp_SLProps props(surface, u, v, 1, surface.Tolerance());
  if(!props.IsNormalDefined()) {
    return false;
  }

  const auto direction = props.Normal();
  normal = gp_Vec(direction.X(), direction.Y(), direction.Z());
  if(face.Orientation() == TopAbs_REVERSED) {
    normal.Reverse();
  }
  return normal.SquareMagnitude() > 0.0;
}

[[nodiscard]] bool is_inside_face_material(TopAbs_State state) noexcept
{
  return state == TopAbs_IN || state == TopAbs_ON;
}

[[nodiscard]] Handle(Geom_Surface) world_surface_handle(const TopoDS_Face &face)
{
  TopLoc_Location location;
  Handle(Geom_Surface) surface = BRep_Tool::Surface(face, location);
  if(surface.IsNull()) {
    return surface;
  }

  if(!location.IsIdentity()) {
    surface = Handle(Geom_Surface)::DownCast(
      surface->Transformed(location.Transformation())
    );
  }

  return surface;
}

[[nodiscard]] base::StatusCode project_point_onto_face(
  const TopoDS_Face &face,
  const BRepAdaptor_Surface &surface,
  const geo::Point3 &point,
  FaceProjectionData &projection
)
{
  projection = {};

  const Handle(Geom_Surface) world_surface = world_surface_handle(face);
  if(world_surface.IsNull()) {
    return core::detail::publish_error(
      base::StatusCode::unsupported,
      "Point-to-face projection requires an OCC face with an accessible surface."
    );
  }

  const double tolerance = std::max(surface.Tolerance(), 1.0e-7);
  GeomAPI_ProjectPointOnSurf projector(
    to_gp_point(point),
    world_surface,
    surface.FirstUParameter(),
    surface.LastUParameter(),
    surface.FirstVParameter(),
    surface.LastVParameter(),
    tolerance
  );
  if(!projector.IsDone() || projector.NbPoints() <= 0) {
    return core::detail::publish_error(
      base::StatusCode::unsupported,
      "Point-to-face projection is only supported when OpenCASCADE can compute an orthogonal projection."
    );
  }

  BRepClass_FaceClassifier classifier;
  bool found = false;
  Standard_Integer best_index = 0;
  double best_distance = 0.0;
  double best_u = 0.0;
  double best_v = 0.0;

  for(Standard_Integer index = 1; index <= projector.NbPoints(); ++index) {
    double u = 0.0;
    double v = 0.0;
    projector.Parameters(index, u, v);
    classifier.Perform(face, gp_Pnt2d(u, v), tolerance);
    if(!is_inside_face_material(classifier.State())) {
      continue;
    }

    const double distance = projector.Distance(index);
    if(!found || distance < best_distance) {
      found = true;
      best_index = index;
      best_distance = distance;
      best_u = u;
      best_v = v;
    }
  }

  if(!found) {
    return core::detail::publish_error(
      base::StatusCode::unsupported,
      "Point-to-face projection currently requires the orthogonal solution to lie on the trimmed face material."
    );
  }

  projection.projected_point = projector.Point(best_index);
  projection.u = best_u;
  projection.v = best_v;
  projection.distance = best_distance;
  return core::detail::clear_error_state();
}

[[nodiscard]] bool sample_face_normal_near_edge(
  const TopoDS_Face &face,
  const TopoDS_Edge &edge,
  gp_Vec &normal
)
{
  BRepAdaptor_Surface surface(face, Standard_True);

  Standard_Real u_min = 0.0;
  Standard_Real u_max = 0.0;
  Standard_Real v_min = 0.0;
  Standard_Real v_max = 0.0;
  BRepTools::UVBounds(face, edge, u_min, u_max, v_min, v_max);

  const double u_mid = 0.5 * (u_min + u_max);
  const double v_mid = 0.5 * (v_min + v_max);
  if(oriented_face_normal(face, surface, u_mid, v_mid, normal)) {
    return true;
  }

  BRepTools::UVBounds(face, u_min, u_max, v_min, v_max);
  return oriented_face_normal(
    face,
    surface,
    0.5 * (u_min + u_max),
    0.5 * (v_min + v_max),
    normal
  );
}

[[nodiscard]] std::size_t count_wire_edges(const TopoDS_Wire &wire)
{
  std::size_t count = 0;
  for(TopExp_Explorer explorer(wire, TopAbs_EDGE); explorer.More(); explorer.Next()) {
    ++count;
  }
  return count;
}

[[nodiscard]] geo::FaceBoundaryLoopKind boundary_loop_kind(
  const TopoDS_Wire &wire,
  const TopoDS_Wire &outer_wire
) noexcept
{
  if(outer_wire.IsNull()) {
    return geo::FaceBoundaryLoopKind::unknown;
  }
  return wire.IsSame(outer_wire) ? geo::FaceBoundaryLoopKind::outer
                                 : geo::FaceBoundaryLoopKind::inner;
}

[[nodiscard]] base::StatusCode append_face_boundary_loop(
  const TopologyCache &cache,
  const TopoDS_Face &face,
  const TopoDS_Wire &wire,
  geo::FaceBoundaryLoopKind kind,
  geo::FaceBoundaryLoops &boundary
)
{
  geo::FaceBoundaryLoop loop;
  loop.kind = kind;

  const auto expected_edge_count = count_wire_edges(wire);
  loop.edge_uses.reserve(expected_edge_count);

  BRepTools_WireExplorer explorer(wire, face);
  for(; explorer.More(); explorer.Next()) {
    const TopoDS_Edge current_edge = TopoDS::Edge(explorer.Current());
    const auto edge_entity =
      lookup_topology_entity_id(cache, geo::TopologyDimension::edge, current_edge);
    if(!geo::is_valid(edge_entity)) {
      return core::detail::publish_error(
        base::StatusCode::internal_error,
        "Face boundary loop extraction failed because a wire edge is missing from the topology cache."
      );
    }

    TopoDS_Vertex start_vertex_shape;
    TopoDS_Vertex end_vertex_shape;
    TopExp::Vertices(
      TopoDS::Edge(current_edge.Oriented(explorer.Orientation())),
      start_vertex_shape,
      end_vertex_shape,
      true
    );

    const auto start_vertex =
      lookup_topology_entity_id(cache, geo::TopologyDimension::vertex, start_vertex_shape);
    const auto end_vertex =
      lookup_topology_entity_id(cache, geo::TopologyDimension::vertex, end_vertex_shape);
    if(!geo::is_valid(start_vertex) || !geo::is_valid(end_vertex)) {
      return core::detail::publish_error(
        base::StatusCode::internal_error,
        "Face boundary loop extraction failed because a wire vertex is missing from the topology cache."
      );
    }

    const auto *representative_edge = lookup_topology_shape(cache, edge_entity);
    if(representative_edge == nullptr) {
      return core::detail::publish_error(
        base::StatusCode::internal_error,
        "Face boundary loop extraction failed because an edge representative is missing from the topology cache."
      );
    }

    loop.edge_uses.push_back(
      {
        edge_entity,
        start_vertex,
        end_vertex,
        explorer.Orientation() == TopAbs_FORWARD,
        BRep_Tool::IsClosed(current_edge, face),
        BRep_Tool::Degenerated(current_edge),
      }
    );
  }

  if(loop.edge_uses.empty()) {
    return core::detail::publish_error(
      base::StatusCode::unsupported,
      "Face boundary loop extraction produced an empty edge sequence."
    );
  }
  if(loop.edge_uses.size() != expected_edge_count) {
    SQMESH_LOG_WARN(
      "Face boundary wire: BRepTools_WireExplorer walked {} edges but wire "
      "contains {} edge uses (likely seam or slightly malformed wire); "
      "continuing with the ordered sequence.",
      loop.edge_uses.size(), expected_edge_count);
  }

  if(!loop.edge_uses.empty()) {
    loop.closed = loop.edge_uses.front().start_vertex == loop.edge_uses.back().end_vertex;
  }

  boundary.loops.push_back(std::move(loop));
  return core::detail::clear_error_state();
}

[[nodiscard]] base::StatusCode build_face_boundary_cache_entry(
  TopologyCache &cache,
  geo::TopologyEntityId face_entity
)
{
  if(face_entity.dimension != geo::TopologyDimension::face) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Face boundary loop lookup requires a face topology entity id."
    );
  }
  if(face_entity.index >= cache.face_boundaries.size() ||
     face_entity.index >= cache.face_boundary_ready.size()) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Requested face entity does not exist in the geometry topology cache."
    );
  }
  if(cache.face_boundary_ready[face_entity.index] != 0U) {
    return core::detail::clear_error_state();
  }

  const auto *shape = lookup_topology_shape(cache, face_entity);
  if(shape == nullptr) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Requested face entity does not exist in the geometry topology cache."
    );
  }

  const TopoDS_Face face = TopoDS::Face(*shape);
  geo::FaceBoundaryLoops boundary;
  boundary.face = face_entity;

  const TopoDS_Wire outer_wire = BRepTools::OuterWire(face);
  if(!outer_wire.IsNull()) {
    const auto status = append_face_boundary_loop(
      cache,
      face,
      outer_wire,
      geo::FaceBoundaryLoopKind::outer,
      boundary
    );
    if(status != base::StatusCode::ok) {
      return status;
    }
  }

  for(TopExp_Explorer explorer(face, TopAbs_WIRE); explorer.More(); explorer.Next()) {
    const TopoDS_Wire wire = TopoDS::Wire(explorer.Current());
    if(!outer_wire.IsNull() && wire.IsSame(outer_wire)) {
      continue;
    }

    const auto status = append_face_boundary_loop(
      cache,
      face,
      wire,
      boundary_loop_kind(wire, outer_wire),
      boundary
    );
    if(status != base::StatusCode::ok) {
      return status;
    }
  }

  cache.face_boundaries[face_entity.index] = std::move(boundary);
  cache.face_boundary_ready[face_entity.index] = 1U;
  return core::detail::clear_error_state();
}

void OccModelStorage::ensure_topology_cache_locked() const
{
  if(topology_cache_valid_) {
    return;
  }

  topology_cache_ = build_topology_cache(shape_, topology_revision_);
  topology_cache_valid_ = true;
}

base::StatusCode OccModelStorage::ensure_proxy_cache_locked() const
{
  if(proxy_cache_valid_ && proxy_cache_.topology_revision == topology_revision_) {
    return core::detail::clear_error_state();
  }

  GeometryCoarseProxyMesh proxy_cache;
  const auto status =
    build_proxy_cache(shape_, topology_cache_, topology_revision_, proxy_cache);
  if(status != base::StatusCode::ok) {
    return status;
  }

  proxy_cache_ = std::move(proxy_cache);
  proxy_cache_valid_ = true;
  return core::detail::clear_error_state();
}

base::StatusCode OccModelStorage::coarse_proxy_mesh(
  GeometryCoarseProxyMesh &proxy_mesh
) const noexcept
{
  proxy_mesh = {};

  try {
    std::lock_guard<std::mutex> lock(topology_mutex_);
    ensure_topology_cache_locked();
    const auto status = ensure_proxy_cache_locked();
    if(status != base::StatusCode::ok) {
      return status;
    }

    proxy_mesh = proxy_cache_;
    return core::detail::clear_error_state();
  }
  catch(const Standard_Failure &failure) {
    return publish_occ_failure(
      "Proxy triangulation cache build failed inside OpenCASCADE.",
      failure
    );
  }
  catch(const std::exception &exception) {
    return publish_std_failure(
      "Proxy triangulation cache build failed unexpectedly.",
      exception
    );
  }
  catch(...) {
    return internal_error("Proxy triangulation cache build failed with an unknown error.");
  }
}

base::StatusCode OccModelStorage::copy_shape(TopoDS_Shape &shape) const noexcept
{
  shape.Nullify();

  try {
    std::lock_guard<std::mutex> lock(topology_mutex_);
    if(shape_.IsNull()) {
      return core::detail::publish_error(
        base::StatusCode::internal_error,
        "The OCC CAD adapter resolved a null model payload."
      );
    }

    shape = shape_;
    return core::detail::clear_error_state();
  }
  catch(const Standard_Failure &failure) {
    return publish_occ_failure(
      "OCC shape copy failed inside OpenCASCADE.",
      failure
    );
  }
}

base::StatusCode OccModelStorage::topology_snapshot(
  geo::TopologySnapshot &snapshot
) const noexcept
{
  snapshot = {};

  try {
    std::lock_guard<std::mutex> lock(topology_mutex_);
    ensure_topology_cache_locked();
    snapshot = topology_cache_.snapshot;
    return core::detail::clear_error_state();
  }
  catch(const Standard_Failure &failure) {
    return publish_occ_failure("Topology snapshot build failed inside OpenCASCADE.", failure);
  }
  catch(const std::exception &exception) {
    return publish_std_failure("Topology snapshot build failed unexpectedly.", exception);
  }
  catch(...) {
    return internal_error("Topology snapshot build failed with an unknown error.");
  }
}

base::StatusCode OccModelStorage::topology_children(
  geo::TopologyEntityId entity,
  std::vector<geo::TopologyEntityId> &children
) const noexcept
{
  children.clear();

  try {
    std::lock_guard<std::mutex> lock(topology_mutex_);
    ensure_topology_cache_locked();
    return copy_topology_neighbors(topology_cache_.children, entity, children, "child");
  }
  catch(const Standard_Failure &failure) {
    return publish_occ_failure("Topology child lookup failed inside OpenCASCADE.", failure);
  }
  catch(const std::exception &exception) {
    return publish_std_failure("Topology child lookup failed unexpectedly.", exception);
  }
  catch(...) {
    return internal_error("Topology child lookup failed with an unknown error.");
  }
}

base::StatusCode OccModelStorage::topology_parents(
  geo::TopologyEntityId entity,
  std::vector<geo::TopologyEntityId> &parents
) const noexcept
{
  parents.clear();

  try {
    std::lock_guard<std::mutex> lock(topology_mutex_);
    ensure_topology_cache_locked();
    return copy_topology_neighbors(topology_cache_.parents, entity, parents, "parent");
  }
  catch(const Standard_Failure &failure) {
    return publish_occ_failure("Topology parent lookup failed inside OpenCASCADE.", failure);
  }
  catch(const std::exception &exception) {
    return publish_std_failure("Topology parent lookup failed unexpectedly.", exception);
  }
  catch(...) {
    return internal_error("Topology parent lookup failed with an unknown error.");
  }
}

base::StatusCode OccModelStorage::face_uv_bounds(
  geo::TopologyEntityId face_entity,
  geo::FaceUvBounds &bounds
) const noexcept
{
  bounds = {};

  if(face_entity.dimension != geo::TopologyDimension::face) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Face UV bounds lookup requires a face topology entity id."
    );
  }

  try {
    std::lock_guard<std::mutex> lock(topology_mutex_);
    ensure_topology_cache_locked();

    const auto *shape = lookup_topology_shape(topology_cache_, face_entity);
    if(shape == nullptr) {
      return core::detail::publish_error(
        base::StatusCode::invalid_argument,
        "Requested face entity does not exist in the geometry topology cache."
      );
    }

    const TopoDS_Face face = TopoDS::Face(*shape);
    Standard_Real u_min = 0.0;
    Standard_Real u_max = 0.0;
    Standard_Real v_min = 0.0;
    Standard_Real v_max = 0.0;
    BRepTools::UVBounds(face, u_min, u_max, v_min, v_max);

    bounds.face = face_entity;
    bounds.u_min = u_min;
    bounds.u_max = u_max;
    bounds.v_min = v_min;
    bounds.v_max = v_max;
    return core::detail::clear_error_state();
  }
  catch(const Standard_Failure &failure) {
    return publish_occ_failure("Face UV bounds query failed inside OpenCASCADE.", failure);
  }
  catch(const std::exception &exception) {
    return publish_std_failure("Face UV bounds query failed unexpectedly.", exception);
  }
  catch(...) {
    return internal_error("Face UV bounds query failed with an unknown error.");
  }
}

base::StatusCode OccModelStorage::sample_face(
  geo::TopologyEntityId face_entity,
  double u,
  double v,
  geo::FaceSample &sample
) const noexcept
{
  sample = {};

  if(face_entity.dimension != geo::TopologyDimension::face) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Face sampling requires a face topology entity id."
    );
  }

  try {
    std::lock_guard<std::mutex> lock(topology_mutex_);
    ensure_topology_cache_locked();

    const auto *shape = lookup_topology_shape(topology_cache_, face_entity);
    if(shape == nullptr) {
      return core::detail::publish_error(
        base::StatusCode::invalid_argument,
        "Requested face entity does not exist in the geometry topology cache."
      );
    }
    const TopoDS_Face face = TopoDS::Face(*shape);
    const auto &surface = cached_surface_locked(face, face_entity.index);
    gp_Pnt point;
    surface.D0(u, v, point);

    sample.face = face_entity;
    sample.u = u;
    sample.v = v;
    sample.position = to_point3(point);

    gp_Vec normal;
    if(oriented_face_normal(face, surface, u, v, normal)) {
      sample.normal = to_vector3(normal);
      sample.normal_defined = true;
    }

    return core::detail::clear_error_state();
  }
  catch(const Standard_Failure &failure) {
    return publish_occ_failure("Face sampling failed inside OpenCASCADE.", failure);
  }
  catch(const std::exception &exception) {
    return publish_std_failure("Face sampling failed unexpectedly.", exception);
  }
  catch(...) {
    return internal_error("Face sampling failed with an unknown error.");
  }
}

base::StatusCode OccModelStorage::sample_face_curvature(
  geo::TopologyEntityId face_entity,
  double u,
  double v,
  geo::FaceCurvatureSample &sample
) const noexcept
{
  sample = {};

  if(face_entity.dimension != geo::TopologyDimension::face) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Face curvature sampling requires a face topology entity id."
    );
  }

  try {
    std::lock_guard<std::mutex> lock(topology_mutex_);
    ensure_topology_cache_locked();

    const auto *shape = lookup_topology_shape(topology_cache_, face_entity);
    if(shape == nullptr) {
      return core::detail::publish_error(
        base::StatusCode::invalid_argument,
        "Requested face entity does not exist in the geometry topology cache."
      );
    }

    const TopoDS_Face face = TopoDS::Face(*shape);
    Standard_Real u_min = 0.0;
    Standard_Real u_max = 0.0;
    Standard_Real v_min = 0.0;
    Standard_Real v_max = 0.0;
    BRepTools::UVBounds(face, u_min, u_max, v_min, v_max);
    if(u < u_min || u > u_max || v < v_min || v > v_max) {
      return core::detail::publish_error(
        base::StatusCode::invalid_argument,
        "Face curvature sampling requires coordinates inside the face UV bounds."
      );
    }

    BRepAdaptor_Surface surface(face, Standard_True);
    BRepLProp_SLProps properties(surface, u, v, 2, surface.Tolerance());

    sample.face = face_entity;
    sample.u = u;
    sample.v = v;

    if(properties.IsCurvatureDefined()) {
      sample.min_curvature = properties.MinCurvature();
      sample.max_curvature = properties.MaxCurvature();
      sample.mean_curvature = properties.MeanCurvature();
      sample.gaussian_curvature = properties.GaussianCurvature();
      sample.curvature_defined = true;
    }

    return core::detail::clear_error_state();
  }
  catch(const Standard_Failure &failure) {
    return publish_occ_failure("Face curvature sampling failed inside OpenCASCADE.", failure);
  }
  catch(const std::exception &exception) {
    return publish_std_failure("Face curvature sampling failed unexpectedly.", exception);
  }
  catch(...) {
    return internal_error("Face curvature sampling failed with an unknown error.");
  }
}

base::StatusCode OccModelStorage::sample_face_derivatives(
  geo::TopologyEntityId face_entity,
  double u,
  double v,
  geo::FaceDerivatives &sample
) const noexcept
{
  sample = {};

  if(face_entity.dimension != geo::TopologyDimension::face) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Face derivative sampling requires a face topology entity id."
    );
  }

  try {
    std::lock_guard<std::mutex> lock(topology_mutex_);
    ensure_topology_cache_locked();

    const auto *shape = lookup_topology_shape(topology_cache_, face_entity);
    if(shape == nullptr) {
      return core::detail::publish_error(
        base::StatusCode::invalid_argument,
        "Requested face entity does not exist in the geometry topology cache."
      );
    }
    const TopoDS_Face face = TopoDS::Face(*shape);
    const auto &surface = cached_surface_locked(face, face_entity.index);

    sample.face = face_entity;
    sample.u = u;
    sample.v = v;

    try {
      gp_Pnt d1_position;
      gp_Vec du;
      gp_Vec dv;
      surface.D1(u, v, d1_position, du, dv);
      sample.position = to_point3(d1_position);
      sample.du = to_vector3(du);
      sample.dv = to_vector3(dv);
      sample.first_derivatives_defined = true;
    }
    catch(const Standard_Failure &) {
    }

    gp_Vec normal;
    if(oriented_face_normal(face, surface, u, v, normal)) {
      sample.normal = to_vector3(normal);
      sample.normal_defined = true;
    }

    return core::detail::clear_error_state();
  }
  catch(const Standard_Failure &failure) {
    return publish_occ_failure("Face derivative sampling failed inside OpenCASCADE.", failure);
  }
  catch(const std::exception &exception) {
    return publish_std_failure("Face derivative sampling failed unexpectedly.", exception);
  }
  catch(...) {
    return internal_error("Face derivative sampling failed with an unknown error.");
  }
}

base::StatusCode OccModelStorage::project_point_to_face(
  geo::TopologyEntityId face_entity,
  const geo::Point3 &point,
  geo::FaceProjection &projection
) const noexcept
{
  projection = {};

  if(face_entity.dimension != geo::TopologyDimension::face) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Face projection requires a face topology entity id."
    );
  }

  try {
    std::lock_guard<std::mutex> lock(topology_mutex_);
    ensure_topology_cache_locked();

    const auto *shape = lookup_topology_shape(topology_cache_, face_entity);
    if(shape == nullptr) {
      return core::detail::publish_error(
        base::StatusCode::invalid_argument,
        "Requested face entity does not exist in the geometry topology cache."
      );
    }

    const TopoDS_Face face = TopoDS::Face(*shape);
    const auto &surface = cached_surface_locked(face, face_entity.index);
    FaceProjectionData data;
    const auto status = project_point_onto_face(face, surface, point, data);
    if(status != base::StatusCode::ok) {
      return status;
    }

    projection.face = face_entity;
    projection.input_point = point;
    projection.projected_point = to_point3(data.projected_point);
    projection.u = data.u;
    projection.v = data.v;
    projection.distance = data.distance;

    gp_Vec normal;
    if(oriented_face_normal(face, surface, data.u, data.v, normal)) {
      projection.normal = to_vector3(normal);
      projection.normal_defined = true;
    }

    return core::detail::clear_error_state();
  }
  catch(const Standard_Failure &failure) {
    return publish_occ_failure("Face projection failed inside OpenCASCADE.", failure);
  }
  catch(const std::exception &exception) {
    return publish_std_failure("Face projection failed unexpectedly.", exception);
  }
  catch(...) {
    return internal_error("Face projection failed with an unknown error.");
  }
}

base::StatusCode OccModelStorage::recover_face_uv(
  geo::TopologyEntityId face_entity,
  const geo::Point3 &point,
  geo::FaceUvMapping &mapping
) const noexcept
{
  mapping = {};

  if(face_entity.dimension != geo::TopologyDimension::face) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Face UV recovery requires a face topology entity id."
    );
  }

  try {
    std::lock_guard<std::mutex> lock(topology_mutex_);
    ensure_topology_cache_locked();

    const auto *shape = lookup_topology_shape(topology_cache_, face_entity);
    if(shape == nullptr) {
      return core::detail::publish_error(
        base::StatusCode::invalid_argument,
        "Requested face entity does not exist in the geometry topology cache."
      );
    }

    const TopoDS_Face face = TopoDS::Face(*shape);
    const auto &surface = cached_surface_locked(face, face_entity.index);
    FaceProjectionData data;
    const auto status = project_point_onto_face(face, surface, point, data);
    if(status != base::StatusCode::ok) {
      return status;
    }

    mapping.face = face_entity;
    mapping.input_point = point;
    mapping.mapped_point = to_point3(data.projected_point);
    mapping.u = data.u;
    mapping.v = data.v;
    mapping.distance = data.distance;
    return core::detail::clear_error_state();
  }
  catch(const Standard_Failure &failure) {
    return publish_occ_failure("Face UV recovery failed inside OpenCASCADE.", failure);
  }
  catch(const std::exception &exception) {
    return publish_std_failure("Face UV recovery failed unexpectedly.", exception);
  }
  catch(...) {
    return internal_error("Face UV recovery failed with an unknown error.");
  }
}

base::StatusCode OccModelStorage::recover_face_uv_from_edge(
  geo::TopologyEntityId face_entity,
  geo::TopologyEntityId edge_entity,
  double edge_parameter,
  geo::FaceUvMapping &mapping
) const noexcept
{
  mapping = {};

  if(face_entity.dimension != geo::TopologyDimension::face) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Face UV recovery from edge requires a face topology entity id."
    );
  }
  if(edge_entity.dimension != geo::TopologyDimension::edge) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Face UV recovery from edge requires an edge topology entity id."
    );
  }

  try {
    std::lock_guard<std::mutex> lock(topology_mutex_);
    ensure_topology_cache_locked();

    const auto *face_shape = lookup_topology_shape(topology_cache_, face_entity);
    if(face_shape == nullptr) {
      return core::detail::publish_error(
        base::StatusCode::invalid_argument,
        "Requested face entity does not exist in the geometry topology cache."
      );
    }
    const auto *edge_shape = lookup_topology_shape(topology_cache_, edge_entity);
    if(edge_shape == nullptr) {
      return core::detail::publish_error(
        base::StatusCode::invalid_argument,
        "Requested edge entity does not exist in the geometry topology cache."
      );
    }

    const TopoDS_Face face = TopoDS::Face(*face_shape);
    const TopoDS_Edge edge = TopoDS::Edge(*edge_shape);
    const geo::FaceBoundaryEdgeUse edge_use {
      edge_entity,
      {},
      {},
      true,
      BRep_Tool::IsClosed(edge, face),
      BRep_Tool::Degenerated(edge),
    };
    return recover_face_uv_from_edge_use(
      face_entity,
      edge_use,
      edge_parameter,
      mapping
    );
  }
  catch(const Standard_Failure &failure) {
    return publish_occ_failure(
      "Face UV recovery from edge pcurve failed inside OpenCASCADE.", failure
    );
  }
  catch(const std::exception &exception) {
    return publish_std_failure(
      "Face UV recovery from edge pcurve failed unexpectedly.", exception
    );
  }
  catch(...) {
    return internal_error(
      "Face UV recovery from edge pcurve failed with an unknown error."
    );
  }
}

base::StatusCode OccModelStorage::recover_face_uv_from_edge_use(
  geo::TopologyEntityId face_entity,
  const geo::FaceBoundaryEdgeUse &edge_use,
  double edge_parameter,
  geo::FaceUvMapping &mapping
) const noexcept
{
  mapping = {};

  if(face_entity.dimension != geo::TopologyDimension::face) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Face UV recovery from edge use requires a face topology entity id."
    );
  }
  if(edge_use.edge.dimension != geo::TopologyDimension::edge) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Face UV recovery from edge use requires an edge topology entity id."
    );
  }

  try {
    std::lock_guard<std::mutex> lock(topology_mutex_);
    ensure_topology_cache_locked();

    const auto *face_shape = lookup_topology_shape(topology_cache_, face_entity);
    if(face_shape == nullptr) {
      return core::detail::publish_error(
        base::StatusCode::invalid_argument,
        "Requested face entity does not exist in the geometry topology cache."
      );
    }
    const auto *edge_shape = lookup_topology_shape(topology_cache_, edge_use.edge);
    if(edge_shape == nullptr) {
      return core::detail::publish_error(
        base::StatusCode::invalid_argument,
        "Requested edge entity does not exist in the geometry topology cache."
      );
    }

    const TopoDS_Face face = TopoDS::Face(*face_shape);
    TopoDS_Edge edge = TopoDS::Edge(*edge_shape);
    edge.Orientation(
      edge_use.same_orientation_as_edge ? TopAbs_FORWARD : TopAbs_REVERSED
    );

    Standard_Real pcurve_first = 0.0;
    Standard_Real pcurve_last = 0.0;
    const Handle(Geom2d_Curve) pcurve =
      BRep_Tool::CurveOnSurface(edge, face, pcurve_first, pcurve_last);
    if(pcurve.IsNull()) {
      return core::detail::publish_error(
        base::StatusCode::unsupported,
        "Edge use does not have a parametric curve (pcurve) on the requested face."
      );
    }

    const double clamped_parameter =
      std::max(pcurve_first, std::min(pcurve_last, edge_parameter));
    const gp_Pnt2d uv_point = pcurve->Value(clamped_parameter);

    const Handle(Geom_Surface) surface = BRep_Tool::Surface(face);
    gp_Pnt surface_point;
    if(!surface.IsNull()) {
      surface_point = surface->Value(uv_point.X(), uv_point.Y());
    }

    mapping.face = face_entity;
    mapping.u = uv_point.X();
    mapping.v = uv_point.Y();
    mapping.mapped_point = to_point3(surface_point);
    mapping.distance = 0.0;
    return core::detail::clear_error_state();
  }
  catch(const Standard_Failure &failure) {
    return publish_occ_failure(
      "Face UV recovery from face-local edge use pcurve failed inside OpenCASCADE.",
      failure
    );
  }
  catch(const std::exception &exception) {
    return publish_std_failure(
      "Face UV recovery from face-local edge use pcurve failed unexpectedly.",
      exception
    );
  }
  catch(...) {
    return internal_error(
      "Face UV recovery from face-local edge use pcurve failed with an unknown error."
    );
  }
}

base::StatusCode OccModelStorage::edge_curve_info(
  geo::TopologyEntityId edge_entity,
  geo::EdgeCurveInfo &info
) const noexcept
{
  info = {};

  if(edge_entity.dimension != geo::TopologyDimension::edge) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Edge geometry lookup requires an edge topology entity id."
    );
  }

  try {
    std::lock_guard<std::mutex> lock(topology_mutex_);
    ensure_topology_cache_locked();

    const auto *shape = lookup_topology_shape(topology_cache_, edge_entity);
    if(shape == nullptr) {
      return core::detail::publish_error(
        base::StatusCode::invalid_argument,
        "Requested edge entity does not exist in the geometry topology cache."
      );
    }

    const TopoDS_Edge edge = TopoDS::Edge(shape->Oriented(TopAbs_FORWARD));
    BRepAdaptor_Curve curve(edge);
    const Standard_Real parameter_min = curve.FirstParameter();
    const Standard_Real parameter_max = curve.LastParameter();

    gp_Pnt start_point;
    gp_Pnt end_point;
    curve.D0(parameter_min, start_point);
    curve.D0(parameter_max, end_point);

    GProp_GProps properties;
    BRepGProp::LinearProperties(edge, properties);

    info.edge = edge_entity;
    info.parameter_min = parameter_min;
    info.parameter_max = parameter_max;
    info.start_point = to_point3(start_point);
    info.end_point = to_point3(end_point);
    info.approximate_length = properties.Mass();
    return core::detail::clear_error_state();
  }
  catch(const Standard_Failure &failure) {
    return publish_occ_failure("Edge geometry lookup failed inside OpenCASCADE.", failure);
  }
  catch(const std::exception &exception) {
    return publish_std_failure("Edge geometry lookup failed unexpectedly.", exception);
  }
  catch(...) {
    return internal_error("Edge geometry lookup failed with an unknown error.");
  }
}

base::StatusCode OccModelStorage::sample_edge_tangent(
  geo::TopologyEntityId edge_entity,
  double parameter,
  geo::EdgeTangentSample &sample
) const noexcept
{
  sample = {};

  if(edge_entity.dimension != geo::TopologyDimension::edge) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Edge tangent sampling requires an edge topology entity id."
    );
  }

  try {
    std::lock_guard<std::mutex> lock(topology_mutex_);
    ensure_topology_cache_locked();

    const auto *shape = lookup_topology_shape(topology_cache_, edge_entity);
    if(shape == nullptr) {
      return core::detail::publish_error(
        base::StatusCode::invalid_argument,
        "Requested edge entity does not exist in the geometry topology cache."
      );
    }

    const TopoDS_Edge edge = TopoDS::Edge(shape->Oriented(TopAbs_FORWARD));
    BRepAdaptor_Curve curve(edge);
    const Standard_Real parameter_min = curve.FirstParameter();
    const Standard_Real parameter_max = curve.LastParameter();
    if(parameter < parameter_min || parameter > parameter_max) {
      return core::detail::publish_error(
        base::StatusCode::invalid_argument,
        "Edge tangent sampling requires a parameter inside the edge parameter range."
      );
    }

    gp_Pnt position;
    curve.D0(parameter, position);

    sample.edge = edge_entity;
    sample.parameter = parameter;
    sample.position = to_point3(position);

    try {
      gp_Pnt d1_position;
      gp_Vec derivative;
      curve.D1(parameter, d1_position, derivative);
      sample.position = to_point3(d1_position);
      sample.derivative = to_vector3(derivative);
      sample.speed = derivative.Magnitude();
      if(sample.speed > curve.Tolerance()) {
        derivative.Divide(sample.speed);
        sample.tangent = to_vector3(derivative);
        sample.tangent_defined = true;
      }
    }
    catch(const Standard_Failure &) {
    }

    return core::detail::clear_error_state();
  }
  catch(const Standard_Failure &failure) {
    return publish_occ_failure("Edge tangent sampling failed inside OpenCASCADE.", failure);
  }
  catch(const std::exception &exception) {
    return publish_std_failure("Edge tangent sampling failed unexpectedly.", exception);
  }
  catch(...) {
    return internal_error("Edge tangent sampling failed with an unknown error.");
  }
}

base::StatusCode OccModelStorage::face_boundary_loops(
  geo::TopologyEntityId face_entity,
  geo::FaceBoundaryLoops &boundary
) const noexcept
{
  boundary = {};

  if(face_entity.dimension != geo::TopologyDimension::face) {
    return core::detail::publish_error(
      base::StatusCode::invalid_argument,
      "Face boundary loop lookup requires a face topology entity id."
    );
  }

  try {
    std::lock_guard<std::mutex> lock(topology_mutex_);
    ensure_topology_cache_locked();

    const auto status = build_face_boundary_cache_entry(topology_cache_, face_entity);
    if(status != base::StatusCode::ok) {
      return status;
    }

    boundary = topology_cache_.face_boundaries[face_entity.index];
    return core::detail::clear_error_state();
  }
  catch(const Standard_Failure &failure) {
    return publish_occ_failure(
      "Face boundary loop extraction failed inside OpenCASCADE.",
      failure
    );
  }
  catch(const std::exception &exception) {
    return publish_std_failure(
      "Face boundary loop extraction failed unexpectedly.",
      exception
    );
  }
  catch(...) {
    return internal_error("Face boundary loop extraction failed with an unknown error.");
  }
}

base::StatusCode OccModelStorage::feature_edges(
  geo::FeatureEdgeReport &report,
  const geo::FeatureEdgeOptions &options
) const noexcept
{
  report = {};

  try {
    std::lock_guard<std::mutex> lock(topology_mutex_);
    ensure_topology_cache_locked();

    const auto &edge_shapes = topology_cache_.shapes[topology_dimension_index(geo::TopologyDimension::edge)];
    const auto &edge_parents =
      topology_cache_.parents[topology_dimension_index(geo::TopologyDimension::edge)];

    report.edges.reserve(edge_shapes.size());

    for(std::size_t edge_index = 0; edge_index < edge_shapes.size(); ++edge_index) {
      const auto edge_entity =
        geo::TopologyEntityId {geo::TopologyDimension::edge, static_cast<std::uint32_t>(edge_index)};
      const auto &parents = edge_parents[edge_index];

      if(parents.size() <= 1U) {
        if(options.include_boundary_edges) {
          report.edges.push_back(edge_entity);
          ++report.boundary_edge_count;
        }
        continue;
      }

      if(parents.size() > 2U) {
        if(options.include_non_manifold_edges) {
          report.edges.push_back(edge_entity);
          ++report.non_manifold_edge_count;
        }
        continue;
      }

      const auto *edge_shape = lookup_topology_shape(topology_cache_, edge_entity);
      const auto *face0_shape = lookup_topology_shape(topology_cache_, parents[0]);
      const auto *face1_shape = lookup_topology_shape(topology_cache_, parents[1]);
      if(edge_shape == nullptr || face0_shape == nullptr || face1_shape == nullptr) {
        continue;
      }

      gp_Vec normal0;
      gp_Vec normal1;
      const bool normal0_defined =
        sample_face_normal_near_edge(TopoDS::Face(*face0_shape), TopoDS::Edge(*edge_shape), normal0);
      const bool normal1_defined =
        sample_face_normal_near_edge(TopoDS::Face(*face1_shape), TopoDS::Edge(*edge_shape), normal1);
      if(!normal0_defined || !normal1_defined) {
        continue;
      }

      const double magnitude_product = normal0.Magnitude() * normal1.Magnitude();
      if(magnitude_product <= 0.0) {
        continue;
      }

      double cosine = normal0.Dot(normal1) / magnitude_product;
      cosine = std::max(-1.0, std::min(1.0, cosine));
      const double angle_degrees = std::acos(cosine) * (180.0 / kPi);
      if(angle_degrees >= options.feature_angle_degrees) {
        report.edges.push_back(edge_entity);
        ++report.sharp_edge_count;
      }
    }

    return core::detail::clear_error_state();
  }
  catch(const Standard_Failure &failure) {
    return publish_occ_failure("Feature edge extraction failed inside OpenCASCADE.", failure);
  }
  catch(const std::exception &exception) {
    return publish_std_failure("Feature edge extraction failed unexpectedly.", exception);
  }
  catch(...) {
    return internal_error("Feature edge extraction failed with an unknown error.");
  }
}

[[nodiscard]] bool same_summary(
  const geo::ModelSummary &lhs,
  const geo::ModelSummary &rhs
) noexcept
{
  return lhs.entity_count == rhs.entity_count &&
         lhs.compound_count == rhs.compound_count &&
         lhs.compsolid_count == rhs.compsolid_count &&
         lhs.solid_count == rhs.solid_count &&
         lhs.shell_count == rhs.shell_count &&
         lhs.face_count == rhs.face_count &&
         lhs.wire_count == rhs.wire_count &&
         lhs.edge_count == rhs.edge_count &&
         lhs.vertex_count == rhs.vertex_count;
}

std::mutex &step_export_mutex() noexcept
{
  static std::mutex instance;
  return instance;
}

[[nodiscard]] geo::ModelSummary summarize_shape(const TopoDS_Shape &shape)
{
  geo::ModelSummary summary;
  if(shape.IsNull()) {
    return summary;
  }

  summary.compound_count = count_shapes(shape, TopAbs_COMPOUND);
  summary.compsolid_count = count_shapes(shape, TopAbs_COMPSOLID);
  summary.solid_count = count_shapes(shape, TopAbs_SOLID);
  summary.shell_count = count_shapes(shape, TopAbs_SHELL);
  summary.face_count = count_shapes(shape, TopAbs_FACE);
  summary.wire_count = count_shapes(shape, TopAbs_WIRE);
  summary.edge_count = count_shapes(shape, TopAbs_EDGE);
  summary.vertex_count = count_shapes(shape, TopAbs_VERTEX);
  summary.entity_count = summary.compound_count + summary.compsolid_count +
                         summary.solid_count + summary.shell_count +
                         summary.face_count + summary.wire_count +
                         summary.edge_count + summary.vertex_count;
  return summary;
}

[[nodiscard]] base::StatusCode invalid_argument(std::string_view message) noexcept
{
  return core::detail::publish_error(base::StatusCode::invalid_argument, message);
}

[[nodiscard]] base::StatusCode io_error(std::string_view message) noexcept
{
  return core::detail::publish_error(base::StatusCode::io_error, message);
}

[[nodiscard]] base::StatusCode unsupported(std::string_view message) noexcept
{
  return core::detail::publish_error(base::StatusCode::unsupported, message);
}

[[nodiscard]] base::StatusCode internal_error(std::string_view message) noexcept
{
  return core::detail::publish_error(base::StatusCode::internal_error, message);
}

[[nodiscard]] base::StatusCode validate_import_path(
  std::string_view path,
  std::string_view label
) noexcept
{
  if(path.empty()) {
    return invalid_argument(label);
  }
  return base::StatusCode::ok;
}

[[nodiscard]] base::StatusCode validate_step_import_options(
  const geo::StepImportOptions &options
) noexcept
{
  if(options.system_length_unit < 0.0) {
    return invalid_argument("STEP import requires a non-negative system length unit.");
  }
  if(options.tolerance < 0.0) {
    return invalid_argument("STEP import requires a non-negative tolerance.");
  }
  if(options.min_tolerance < 0.0) {
    return invalid_argument("STEP import requires a non-negative minimum tolerance.");
  }
  if(options.max_tolerance < 0.0) {
    return invalid_argument("STEP import requires a non-negative maximum tolerance.");
  }
  if(options.max_tolerance > 0.0 && options.max_tolerance < options.min_tolerance) {
    return invalid_argument(
      "STEP import requires max_tolerance to be greater than or equal to min_tolerance."
    );
  }
  return base::StatusCode::ok;
}

[[nodiscard]] base::StatusCode validate_step_export_options(
  const geo::StepExportOptions &options
) noexcept
{
  if(options.linear_tolerance < 0.0) {
    return invalid_argument("STEP export requires a non-negative linear tolerance.");
  }
  return base::StatusCode::ok;
}

[[nodiscard]] base::StatusCode validate_topo_options(
  const geo::TopoOptions &options
) noexcept
{
  if(options.tolerance < 0.0) {
    return invalid_argument("Topo repair requires a non-negative tolerance.");
  }
  if(options.min_tolerance < 0.0) {
    return invalid_argument("Topo repair requires a non-negative minimum tolerance.");
  }
  if(options.max_tolerance < 0.0) {
    return invalid_argument("Topo repair requires a non-negative maximum tolerance.");
  }
  if(options.max_tolerance > 0.0 && options.max_tolerance < options.min_tolerance) {
    return invalid_argument(
      "Topo repair requires max_tolerance to be greater than or equal to min_tolerance."
    );
  }
  return base::StatusCode::ok;
}

[[nodiscard]] base::StatusCode publish_occ_failure(
  std::string_view fallback_message,
  const Standard_Failure &failure
) noexcept
{
  const auto *message = failure.GetMessageString();
  if(message != nullptr && message[0] != '\0') {
    return internal_error(message);
  }
  return internal_error(fallback_message);
}

[[nodiscard]] base::StatusCode publish_std_failure(
  std::string_view fallback_message,
  const std::exception &exception
) noexcept
{
  const auto *message = exception.what();
  if(message != nullptr && message[0] != '\0') {
    return internal_error(message);
  }
  return internal_error(fallback_message);
}

[[nodiscard]] base::StatusCode store_occ_model(
  TopoDS_Shape shape,
  geo::ModelHandle &model_handle,
  base::ContextHandle context_handle
) noexcept
{
  model_handle = sqmesh::invalid_handle;
  if(shape.IsNull()) {
    return invalid_argument("A non-null OpenCASCADE shape is required.");
  }

  const auto summary = summarize_shape(shape);
  auto storage = std::make_shared<OccModelStorage>(std::move(shape));
  return core::detail::store_model(std::move(storage), summary, model_handle, context_handle);
}

[[nodiscard]] base::StatusCode resolve_occ_storage(
  geo::ModelHandle model_handle,
  base::ContextHandle context_handle,
  std::shared_ptr<OccModelStorage> &storage
) noexcept
{
  storage.reset();

  const GeometryModelStoragePtr generic_storage =
    core::detail::lookup_model_storage(model_handle, context_handle);
  if(!generic_storage) {
    const auto status = base::last_error_code();
    if(status != base::StatusCode::ok) {
      return status;
    }
    return unsupported("The model handle is not backed by CAD geometry.");
  }
  if(generic_storage->kernel() != GeometryKernel::occ) {
    return unsupported("The model handle is not backed by the OCC CAD adapter.");
  }

  storage = std::dynamic_pointer_cast<OccModelStorage>(generic_storage);
  if(!storage) {
    return internal_error("The OCC CAD adapter encountered an unexpected model payload.");
  }

  TopoDS_Shape shape;
  const auto copy_status = storage->copy_shape(shape);
  if(copy_status != base::StatusCode::ok) {
    return copy_status;
  }
  if(shape.IsNull()) {
    return internal_error("The OCC CAD adapter resolved a null model payload.");
  }

  return base::StatusCode::ok;
}

[[nodiscard]] base::StatusCode resolve_occ_shape(
  geo::ModelHandle model_handle,
  base::ContextHandle context_handle,
  TopoDS_Shape &shape
) noexcept
{
  shape.Nullify();
  return core::detail::with_model_storage(
    model_handle,
    [&](const GeometryModelStorage &storage) {
      if(storage.kernel() != GeometryKernel::occ) {
        return unsupported("The model handle is not backed by the OCC CAD adapter.");
      }

      const auto *occ_storage = dynamic_cast<const OccModelStorage *>(&storage);
      if(occ_storage == nullptr) {
        return internal_error(
          "The OCC CAD adapter encountered an unexpected model payload."
        );
      }
      return occ_storage->copy_shape(shape);
    },
    context_handle
  );
}

[[nodiscard]] geo::TopologyCheckReport inspect_topology(const TopoDS_Shape &shape)
{
  geo::TopologyCheckReport report;
  report.model_summary = summarize_shape(shape);
  if(shape.IsNull()) {
    return report;
  }

  BRepCheck_Analyzer analyzer(shape, Standard_False);
  report.is_valid = analyzer.IsValid();

  TopTools_IndexedDataMapOfShapeListOfShape edge_faces;
  TopExp::MapShapesAndUniqueAncestors(shape, TopAbs_EDGE, TopAbs_FACE, edge_faces);
  for(Standard_Integer index = 1; index <= edge_faces.Extent(); ++index) {
    const auto face_count =
      static_cast<std::size_t>(edge_faces.FindFromIndex(index).Extent());
    if(face_count <= 1U) {
      ++report.free_edge_count;
      continue;
    }
    if(face_count == 2U) {
      ++report.contiguous_edge_count;
      continue;
    }
    ++report.multiple_edge_count;
  }

  return report;
}

[[nodiscard]] geo::TopoOptions step_import_topo_options(
  const geo::StepImportOptions &options
) noexcept
{
  geo::TopoOptions topo_options;
  topo_options.tolerance = options.tolerance;
  topo_options.min_tolerance = options.min_tolerance;
  topo_options.max_tolerance = options.max_tolerance;
  topo_options.fix_degenerated = options.fix_degenerated;
  topo_options.fix_small_edges = options.fix_small_edges;
  topo_options.fix_small_faces = options.fix_small_faces;
  topo_options.sew_faces = options.sew_faces;
  topo_options.make_solids = options.make_solids;
  return topo_options;
}

[[nodiscard]] geo::TopoReport repair_occ_shape_in_place(
  TopoDS_Shape &repaired_shape,
  const geo::TopoOptions &options
)
{
  geo::TopoReport report;
  report.before = inspect_topology(repaired_shape);

  if(options.sew_faces) {
    BRepBuilderAPI_Sewing sewing(
      options.tolerance,
      Standard_True,
      Standard_True,
      Standard_True,
      Standard_False
    );
    sewing.SetTolerance(options.tolerance);
    if(options.min_tolerance > 0.0) {
      sewing.SetMinTolerance(options.min_tolerance);
    }
    if(options.max_tolerance > 0.0) {
      sewing.SetMaxTolerance(options.max_tolerance);
    }
    sewing.SetLocalTolerancesMode(Standard_True);
    sewing.Add(repaired_shape);
    sewing.Perform();

    report.sewing_performed = true;
    const TopoDS_Shape sewed_shape = sewing.SewedShape();
    if(!sewed_shape.IsNull()) {
      report.sewing_modified =
        !sewed_shape.IsSame(repaired_shape) ||
        sewing.NbFreeEdges() !=
          static_cast<Standard_Integer>(report.before.free_edge_count);
      repaired_shape = sewed_shape;
    }
  }

  const bool run_shape_fix = options.fix_degenerated || options.fix_small_edges ||
                             options.fix_small_faces || options.make_solids;
  if(run_shape_fix) {
    ShapeFix_Shape shape_fix(repaired_shape);
    shape_fix.SetPrecision(options.tolerance);
    if(options.min_tolerance > 0.0) {
      shape_fix.SetMinTolerance(options.min_tolerance);
    }
    if(options.max_tolerance > 0.0) {
      shape_fix.SetMaxTolerance(options.max_tolerance);
    }
    shape_fix.FixSolidMode() = options.make_solids ? 1 : 0;
    shape_fix.FixFreeShellMode() =
      (options.fix_degenerated || options.make_solids) ? 1 : 0;
    shape_fix.FixFreeFaceMode() = options.fix_small_faces ? 1 : 0;
    shape_fix.FixFreeWireMode() = options.fix_small_edges ? 1 : 0;
    shape_fix.FixSameParameterMode() = 1;
    shape_fix.FixVertexPositionMode() = options.fix_degenerated ? 1 : 0;
    shape_fix.FixVertexTolMode() = options.fix_degenerated ? 1 : 0;
    static_cast<void>(shape_fix.Perform());

    report.shape_fix_performed = true;
    report.shape_fix_modified = shape_fix.Status(ShapeExtend_DONE);
    repaired_shape = shape_fix.Shape();
  }

  report.after = inspect_topology(repaired_shape);
  report.free_edges_reduced =
    report.after.free_edge_count < report.before.free_edge_count;
  report.modified = report.sewing_modified || report.shape_fix_modified ||
                    report.free_edges_reduced ||
                    report.after.is_valid != report.before.is_valid ||
                    report.after.contiguous_edge_count != report.before.contiguous_edge_count ||
                    report.after.multiple_edge_count != report.before.multiple_edge_count ||
                    !same_summary(
                      report.before.model_summary,
                      report.after.model_summary
                    );
  return report;
}

[[nodiscard]] base::StatusCode persist_repaired_shape(
  geo::ModelHandle model_handle,
  base::ContextHandle context_handle,
  const std::shared_ptr<OccModelStorage> &storage,
  TopoDS_Shape repaired_shape
) noexcept
{
  if(!storage) {
    return internal_error("The OCC CAD adapter lost the model storage during topo repair.");
  }
  if(repaired_shape.IsNull()) {
    return internal_error("Topo repair produced a null OpenCASCADE shape.");
  }

  const auto summary = summarize_shape(repaired_shape);
  storage->set_shape(std::move(repaired_shape));
  return core::detail::update_model(model_handle, storage, summary, context_handle);
}

} // namespace

bool cad_io_available() noexcept
{
  return true;
}

base::StatusCode store_shape(
  const TopoDS_Shape &shape,
  geo::ModelHandle &model_handle,
  base::ContextHandle context_handle
) noexcept
{
  try {
    return store_occ_model(shape, model_handle, context_handle);
  }
  catch(const Standard_Failure &failure) {
    model_handle = sqmesh::invalid_handle;
    return publish_occ_failure("OpenCASCADE shape storage failed.", failure);
  }
  catch(const std::exception &exception) {
    model_handle = sqmesh::invalid_handle;
    return publish_std_failure("OpenCASCADE shape storage failed unexpectedly.", exception);
  }
  catch(...) {
    model_handle = sqmesh::invalid_handle;
    return internal_error("OpenCASCADE shape storage failed with an unknown error.");
  }
}

base::StatusCode import_step(
  std::string_view path,
  geo::ModelHandle &model_handle,
  const geo::StepImportOptions &options,
  base::ContextHandle context_handle
) noexcept
{
  model_handle = sqmesh::invalid_handle;

  const auto path_status = validate_import_path(path, "STEP import requires a non-empty path.");
  if(path_status != base::StatusCode::ok) {
    return path_status;
  }

  const auto options_status = validate_step_import_options(options);
  if(options_status != base::StatusCode::ok) {
    return options_status;
  }

  try {
    const std::string file_path(path);
    std::lock_guard<std::mutex> lock(step_export_mutex());
    STEPControl_Reader reader;
    ScopedInterfaceInteger read_ideas_guard("read.step.ideas", 1);
    ScopedInterfaceInteger read_nonmanifold_guard("read.step.nonmanifold", 1);

    const auto read_status = reader.ReadFile(file_path.c_str());
    if(read_status != IFSelect_RetDone) {
      return io_error("STEP import failed while reading the file.");
    }

    if(options.system_length_unit > 0.0) {
#if OCC_VERSION_HEX >= 0x070600
      reader.SetSystemLengthUnit(options.system_length_unit);
#else
      return unsupported(
        "STEP import system_length_unit requires OpenCASCADE 7.6 or newer."
      );
#endif
    }

    if(reader.TransferRoots() <= 0) {
      return io_error("STEP import transferred no shapes.");
    }

    TopoDS_Shape shape = reader.OneShape();
    BRepTools::Clean(shape);

    const auto topo_options = step_import_topo_options(options);
    static_cast<void>(repair_occ_shape_in_place(shape, topo_options));

    return store_occ_model(shape, model_handle, context_handle);
  }
  catch(const Standard_Failure &failure) {
    return publish_occ_failure("STEP import failed inside OpenCASCADE.", failure);
  }
  catch(const std::exception &exception) {
    return publish_std_failure("STEP import failed unexpectedly.", exception);
  }
  catch(...) {
    return internal_error("STEP import failed with an unknown error.");
  }
}

base::StatusCode import_iges(
  std::string_view path,
  geo::ModelHandle &model_handle,
  const geo::IgesImportOptions &options,
  base::ContextHandle context_handle
) noexcept
{
  model_handle = sqmesh::invalid_handle;

  const auto path_status = validate_import_path(path, "IGES import requires a non-empty path.");
  if(path_status != base::StatusCode::ok) {
    return path_status;
  }

  try {
    const std::string file_path(path);
    IGESControl_Reader reader;
    reader.SetReadVisible(options.read_visible_only);

    const auto read_status = reader.ReadFile(file_path.c_str());
    if(read_status != IFSelect_RetDone) {
      return io_error("IGES import failed while reading the file.");
    }

    if(reader.TransferRoots() <= 0) {
      return io_error("IGES import transferred no shapes.");
    }

    return store_occ_model(reader.OneShape(), model_handle, context_handle);
  }
  catch(const Standard_Failure &failure) {
    return publish_occ_failure("IGES import failed inside OpenCASCADE.", failure);
  }
  catch(const std::exception &exception) {
    return publish_std_failure("IGES import failed unexpectedly.", exception);
  }
  catch(...) {
    return internal_error("IGES import failed with an unknown error.");
  }
}

base::StatusCode export_step(
  geo::ModelHandle model_handle,
  std::string_view path,
  const geo::StepExportOptions &options,
  base::ContextHandle context_handle
) noexcept
{
  const auto path_status = validate_import_path(path, "STEP export requires a non-empty path.");
  if(path_status != base::StatusCode::ok) {
    return path_status;
  }

  const auto options_status = validate_step_export_options(options);
  if(options_status != base::StatusCode::ok) {
    return options_status;
  }

  try {
    TopoDS_Shape shape;
    const auto shape_status = resolve_occ_shape(model_handle, context_handle, shape);
    if(shape_status != base::StatusCode::ok) {
      return shape_status;
    }

    // OCC STEP export options flow through Interface_Static global state,
    // so keep the whole translation/write sequence serialized.
    std::lock_guard<std::mutex> lock(step_export_mutex());

    const std::string file_path(path);
    STEPControl_Writer writer;
    const ScopedInterfaceCString schema_guard("write.step.schema", options.schema);
    const ScopedInterfaceCString unit_guard("write.step.unit", options.unit_name);

    if(options.linear_tolerance > 0.0) {
      writer.SetTolerance(options.linear_tolerance);
    }

    const auto transfer_status = writer.Transfer(shape, STEPControl_AsIs);
    if(transfer_status != IFSelect_RetDone) {
      return io_error("STEP export failed while translating the OCC shape.");
    }

    const auto write_status = writer.Write(file_path.c_str());
    if(write_status != IFSelect_RetDone) {
      return io_error("STEP export failed while writing the output file.");
    }

    return base::StatusCode::ok;
  }
  catch(const Standard_Failure &failure) {
    return publish_occ_failure("STEP export failed inside OpenCASCADE.", failure);
  }
  catch(const std::exception &exception) {
    return publish_std_failure("STEP export failed unexpectedly.", exception);
  }
  catch(...) {
    return internal_error("STEP export failed with an unknown error.");
  }
}

base::StatusCode export_iges(
  geo::ModelHandle model_handle,
  std::string_view path,
  const geo::IgesExportOptions &options,
  base::ContextHandle context_handle
) noexcept
{
  const auto path_status = validate_import_path(path, "IGES export requires a non-empty path.");
  if(path_status != base::StatusCode::ok) {
    return path_status;
  }

  try {
    TopoDS_Shape shape;
    const auto shape_status = resolve_occ_shape(model_handle, context_handle, shape);
    if(shape_status != base::StatusCode::ok) {
      return shape_status;
    }

    const std::string unit_name = options.unit_name.empty() ? "MM" : std::string(options.unit_name);
    const std::string file_path(path);
    const Standard_Integer write_mode =
      options.write_mode == geo::IgesWriteMode::brep ? 1 : 0;

    IGESControl_Writer writer(unit_name.c_str(), write_mode);
    if(!writer.AddShape(shape)) {
      return io_error("IGES export could not translate the OCC shape.");
    }
    if(!writer.Write(file_path.c_str())) {
      return io_error("IGES export could not write the output file.");
    }

    return base::StatusCode::ok;
  }
  catch(const Standard_Failure &failure) {
    return publish_occ_failure("IGES export failed inside OpenCASCADE.", failure);
  }
  catch(const std::exception &exception) {
    return publish_std_failure("IGES export failed unexpectedly.", exception);
  }
  catch(...) {
    return internal_error("IGES export failed with an unknown error.");
  }
}

base::StatusCode check_topology(
  geo::ModelHandle model_handle,
  geo::TopologyCheckReport &report,
  base::ContextHandle context_handle
) noexcept
{
  report = {};

  try {
    TopoDS_Shape shape;
    const auto shape_status = resolve_occ_shape(model_handle, context_handle, shape);
    if(shape_status != base::StatusCode::ok) {
      return shape_status;
    }

    report = inspect_topology(shape);
    return base::StatusCode::ok;
  }
  catch(const Standard_Failure &failure) {
    return publish_occ_failure("Topology checking failed inside OpenCASCADE.", failure);
  }
  catch(const std::exception &exception) {
    return publish_std_failure("Topology checking failed unexpectedly.", exception);
  }
  catch(...) {
    return internal_error("Topology checking failed with an unknown error.");
  }
}

base::StatusCode free_edge_count(
  geo::ModelHandle model_handle,
  std::size_t &count,
  base::ContextHandle context_handle
) noexcept
{
  count = 0U;

  geo::TopologyCheckReport report;
  const auto status = sqmesh::cad::occ::check_topology(
    model_handle,
    report,
    context_handle
  );
  if(status != base::StatusCode::ok) {
    return status;
  }

  count = report.free_edge_count;
  return base::StatusCode::ok;
}

base::StatusCode topo(
  geo::ModelHandle model_handle,
  geo::TopoReport &report,
  const geo::TopoOptions &options,
  base::ContextHandle context_handle
) noexcept
{
  report = {};

  const auto options_status = validate_topo_options(options);
  if(options_status != base::StatusCode::ok) {
    return options_status;
  }

  try {
    std::shared_ptr<OccModelStorage> storage;
    const auto storage_status = resolve_occ_storage(model_handle, context_handle, storage);
    if(storage_status != base::StatusCode::ok) {
      return storage_status;
    }

    report.topology_revision_before = storage->topology_revision();
    TopoDS_Shape repaired_shape;
    const auto shape_status = storage->copy_shape(repaired_shape);
    if(shape_status != base::StatusCode::ok) {
      return shape_status;
    }
    report = repair_occ_shape_in_place(repaired_shape, options);
    report.topology_revision_after = report.topology_revision_before;
    report.topology_identity_changed = false;

    if(report.modified) {
      const auto persist_status = persist_repaired_shape(
        model_handle,
        context_handle,
        storage,
        std::move(repaired_shape)
      );
      if(persist_status != base::StatusCode::ok) {
        report.topology_revision_after = 0U;
        return persist_status;
      }
      report.topology_revision_after = storage->topology_revision();
      report.topology_identity_changed =
        report.topology_revision_after != report.topology_revision_before;
      return base::StatusCode::ok;
    }

    return base::StatusCode::ok;
  }
  catch(const Standard_Failure &failure) {
    return publish_occ_failure("Topo repair failed inside OpenCASCADE.", failure);
  }
  catch(const std::exception &exception) {
    return publish_std_failure("Topo repair failed unexpectedly.", exception);
  }
  catch(...) {
    return internal_error("Topo repair failed with an unknown error.");
  }
}

} // namespace sqmesh::cad::occ
