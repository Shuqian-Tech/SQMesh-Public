from __future__ import annotations

from . import _sqmesh
from . import base

_raw = _sqmesh.geo

IgesWriteMode = _raw.IgesWriteMode
TopologyDimension = _raw.TopologyDimension
FaceBoundaryLoopKind = _raw.FaceBoundaryLoopKind

TopologyEntityId = _raw.TopologyEntityId
TopologyEntityInfo = _raw.TopologyEntityInfo
TopologySnapshot = _raw.TopologySnapshot
FaceUvBounds = _raw.FaceUvBounds
FaceSample = _raw.FaceSample
FaceCurvatureSample = _raw.FaceCurvatureSample
FaceDerivatives = _raw.FaceDerivatives
FaceProjection = _raw.FaceProjection
FaceUvMapping = _raw.FaceUvMapping
EdgeCurveInfo = _raw.EdgeCurveInfo
EdgeTangentSample = _raw.EdgeTangentSample
EdgeCurveSamplingOptions = _raw.EdgeCurveSamplingOptions
EdgeCurveSamples = _raw.EdgeCurveSamples
FaceBoundaryEdgeUse = _raw.FaceBoundaryEdgeUse
FaceBoundaryLoop = _raw.FaceBoundaryLoop
FaceBoundaryLoops = _raw.FaceBoundaryLoops
VertexView = _raw.VertexView
EdgeView = _raw.EdgeView
FaceView = _raw.FaceView
RegionView = _raw.RegionView
ModelView = _raw.ModelView
FeatureEdgeOptions = _raw.FeatureEdgeOptions
FeatureEdgeReport = _raw.FeatureEdgeReport
ModelSummary = _raw.ModelSummary
StepImportOptions = _raw.StepImportOptions
IgesImportOptions = _raw.IgesImportOptions
StlImportOptions = _raw.StlImportOptions
TopologyCheckReport = _raw.TopologyCheckReport
TopoOptions = _raw.TopoOptions
TopoReport = _raw.TopoReport

INVALID_HANDLE = _sqmesh.invalid_handle
INVALID_TOPOLOGY_INDEX = _raw.invalid_topology_index
INVALID_BOUNDARY_LOOP_INDEX = _raw.invalid_boundary_loop_index


def _resolve_context_handle(
    context: base.Context | None = None,
    owner: "Model | None" = None,
) -> int:
    if context is not None:
        return base._coerce_context_handle(context)
    if owner is not None and owner.context is not None:
        return base._coerce_context_handle(owner.context)
    return INVALID_HANDLE


def _resolve_context_object(
    context: base.Context | None = None,
    owner: "Model | None" = None,
) -> base.Context | None:
    if context is not None:
        return context
    if owner is not None:
        return owner.context
    return base.current_context()


def _coerce_model_handle(model: "Model") -> int:
    if not isinstance(model, Model):
        raise TypeError("model must be a sqmesh.geo.Model")
    model._require_live_context()
    return model.handle


def _build_options(option_type, options=None, **kwargs):
    if options is None:
        options = option_type()
    elif not isinstance(options, option_type):
        raise TypeError(f"options must be a {option_type.__name__}")
    for key, value in kwargs.items():
        setattr(options, key, value)
    return options


class Model:
    def __init__(self, handle: int, context: base.Context | None = None) -> None:
        self._handle = handle
        self._context = context

    @property
    def handle(self) -> int:
        return self._handle

    @property
    def context(self) -> base.Context | None:
        return self._context

    def _require_live_context(self) -> None:
        if self._context is not None:
            self._context._require_open()

    def summary(self) -> ModelSummary:
        return model_summary(self)

    def topology_snapshot(self) -> TopologySnapshot:
        return topology_snapshot(self)

    def view(self) -> ModelView:
        return model_view(self)

    def topo(self, options: TopoOptions | None = None, **kwargs) -> TopoReport:
        return topo(self, options=options, **kwargs)

    def check_topology(self) -> TopologyCheckReport:
        return check_topology(self)

    def free_edge_count(self) -> int:
        return free_edge_count(self)

    def proxy_mesh(self):
        return model_proxy_mesh(self)

    def export_step(
        self,
        path: str,
        *,
        linear_tolerance: float = 0.0,
        unit_name: str = "",
        schema: str = "",
        context: base.Context | None = None,
    ) -> None:
        export_step(
            self,
            path,
            linear_tolerance=linear_tolerance,
            unit_name=unit_name,
            schema=schema,
            context=context,
        )

    def export_iges(
        self,
        path: str,
        *,
        unit_name: str = "",
        write_mode: IgesWriteMode = IgesWriteMode.brep,
        context: base.Context | None = None,
    ) -> None:
        export_iges(
            self,
            path,
            unit_name=unit_name,
            write_mode=write_mode,
            context=context,
        )

    def __repr__(self) -> str:
        return f"<sqmesh.geo.Model handle={self.handle}>"


def module_name() -> str:
    return _raw.module_name()


def cad_io_available() -> bool:
    return _raw.cad_io_available()


def create_placeholder_model(*, context: base.Context | None = None) -> Model:
    resolved = _resolve_context_object(context)
    handle = _raw.create_placeholder_model(_resolve_context_handle(context))
    return Model(handle, resolved)


def import_step(
    path: str,
    *,
    context: base.Context | None = None,
    options: StepImportOptions | None = None,
    **kwargs,
) -> Model:
    resolved = _resolve_context_object(context)
    handle = _raw.import_step(
        path,
        _build_options(StepImportOptions, options, **kwargs),
        _resolve_context_handle(context),
    )
    return Model(handle, resolved)


def import_iges(
    path: str,
    *,
    context: base.Context | None = None,
    options: IgesImportOptions | None = None,
    **kwargs,
) -> Model:
    resolved = _resolve_context_object(context)
    handle = _raw.import_iges(
        path,
        _build_options(IgesImportOptions, options, **kwargs),
        _resolve_context_handle(context),
    )
    return Model(handle, resolved)


def import_stl(
    path: str,
    *,
    context: base.Context | None = None,
    options: StlImportOptions | None = None,
    **kwargs,
) -> Model:
    resolved = _resolve_context_object(context)
    handle = _raw.import_stl(
        path,
        _build_options(StlImportOptions, options, **kwargs),
        _resolve_context_handle(context),
    )
    return Model(handle, resolved)


def export_step(
    model: Model,
    path: str,
    *,
    linear_tolerance: float = 0.0,
    unit_name: str = "",
    schema: str = "",
    context: base.Context | None = None,
) -> None:
    _raw.export_step(
        _coerce_model_handle(model),
        path,
        linear_tolerance,
        unit_name,
        schema,
        _resolve_context_handle(context, model),
    )


def export_iges(
    model: Model,
    path: str,
    *,
    unit_name: str = "",
    write_mode: IgesWriteMode = IgesWriteMode.brep,
    context: base.Context | None = None,
) -> None:
    _raw.export_iges(
        _coerce_model_handle(model),
        path,
        unit_name,
        write_mode,
        _resolve_context_handle(context, model),
    )


def model_summary(
    model: Model,
    *,
    context: base.Context | None = None,
) -> ModelSummary:
    return _raw.model_summary(
        _coerce_model_handle(model),
        _resolve_context_handle(context, model),
    )


def model_proxy_mesh(model: Model, *, context: base.Context | None = None):
    from . import mesh

    resolved = _resolve_context_object(context, model)
    handle = _raw.model_proxy_mesh(
        _coerce_model_handle(model),
        _resolve_context_handle(context, model),
    )
    return mesh.Mesh(handle, resolved)


def topology_snapshot(
    model: Model,
    *,
    context: base.Context | None = None,
) -> TopologySnapshot:
    return _raw.topology_snapshot(
        _coerce_model_handle(model),
        _resolve_context_handle(context, model),
    )


def model_view(model: Model, *, context: base.Context | None = None) -> ModelView:
    return _raw.model_view(
        _coerce_model_handle(model),
        _resolve_context_handle(context, model),
    )


def topology_children(
    model: Model,
    entity: TopologyEntityId,
    *,
    context: base.Context | None = None,
):
    return _raw.topology_children(
        _coerce_model_handle(model),
        entity,
        _resolve_context_handle(context, model),
    )


def topology_parents(
    model: Model,
    entity: TopologyEntityId,
    *,
    context: base.Context | None = None,
):
    return _raw.topology_parents(
        _coerce_model_handle(model),
        entity,
        _resolve_context_handle(context, model),
    )


def face_uv_bounds(
    model: Model,
    face_entity: TopologyEntityId,
    *,
    context: base.Context | None = None,
) -> FaceUvBounds:
    return _raw.face_uv_bounds(
        _coerce_model_handle(model),
        face_entity,
        _resolve_context_handle(context, model),
    )


def sample_face(
    model: Model,
    face_entity: TopologyEntityId,
    u: float,
    v: float,
    *,
    context: base.Context | None = None,
) -> FaceSample:
    return _raw.sample_face(
        _coerce_model_handle(model),
        face_entity,
        u,
        v,
        _resolve_context_handle(context, model),
    )


def sample_face_curvature(
    model: Model,
    face_entity: TopologyEntityId,
    u: float,
    v: float,
    *,
    context: base.Context | None = None,
) -> FaceCurvatureSample:
    return _raw.sample_face_curvature(
        _coerce_model_handle(model),
        face_entity,
        u,
        v,
        _resolve_context_handle(context, model),
    )


def sample_face_derivatives(
    model: Model,
    face_entity: TopologyEntityId,
    u: float,
    v: float,
    *,
    context: base.Context | None = None,
) -> FaceDerivatives:
    return _raw.sample_face_derivatives(
        _coerce_model_handle(model),
        face_entity,
        u,
        v,
        _resolve_context_handle(context, model),
    )


def project_point_to_face(
    model: Model,
    face_entity: TopologyEntityId,
    point,
    *,
    context: base.Context | None = None,
) -> FaceProjection:
    return _raw.project_point_to_face(
        _coerce_model_handle(model),
        face_entity,
        point,
        _resolve_context_handle(context, model),
    )


def recover_face_uv(
    model: Model,
    face_entity: TopologyEntityId,
    point,
    *,
    context: base.Context | None = None,
) -> FaceUvMapping:
    return _raw.recover_face_uv(
        _coerce_model_handle(model),
        face_entity,
        point,
        _resolve_context_handle(context, model),
    )


def edge_curve_info(
    model: Model,
    edge_entity: TopologyEntityId,
    *,
    context: base.Context | None = None,
) -> EdgeCurveInfo:
    return _raw.edge_curve_info(
        _coerce_model_handle(model),
        edge_entity,
        _resolve_context_handle(context, model),
    )


def sample_edge_tangent(
    model: Model,
    edge_entity: TopologyEntityId,
    parameter: float,
    *,
    context: base.Context | None = None,
) -> EdgeTangentSample:
    return _raw.sample_edge_tangent(
        _coerce_model_handle(model),
        edge_entity,
        parameter,
        _resolve_context_handle(context, model),
    )


def sample_edge_curve(
    model: Model,
    edge_entity: TopologyEntityId,
    *,
    context: base.Context | None = None,
    options: EdgeCurveSamplingOptions | None = None,
    **kwargs,
) -> EdgeCurveSamples:
    return _raw.sample_edge_curve(
        _coerce_model_handle(model),
        edge_entity,
        _build_options(EdgeCurveSamplingOptions, options, **kwargs),
        _resolve_context_handle(context, model),
    )


def face_boundary_loops(
    model: Model,
    face_entity: TopologyEntityId,
    *,
    context: base.Context | None = None,
) -> FaceBoundaryLoops:
    return _raw.face_boundary_loops(
        _coerce_model_handle(model),
        face_entity,
        _resolve_context_handle(context, model),
    )


def feature_edges(
    model: Model,
    *,
    context: base.Context | None = None,
    options: FeatureEdgeOptions | None = None,
    **kwargs,
) -> FeatureEdgeReport:
    return _raw.feature_edges(
        _coerce_model_handle(model),
        _build_options(FeatureEdgeOptions, options, **kwargs),
        _resolve_context_handle(context, model),
    )


def check_topology(
    model: Model,
    *,
    context: base.Context | None = None,
) -> TopologyCheckReport:
    return _raw.check_topology(
        _coerce_model_handle(model),
        _resolve_context_handle(context, model),
    )


def free_edge_count(
    model: Model,
    *,
    context: base.Context | None = None,
) -> int:
    return _raw.free_edge_count(
        _coerce_model_handle(model),
        _resolve_context_handle(context, model),
    )


def topo(
    model: Model,
    *,
    context: base.Context | None = None,
    options: TopoOptions | None = None,
    **kwargs,
) -> TopoReport:
    return _raw.topo(
        _coerce_model_handle(model),
        _build_options(TopoOptions, options, **kwargs),
        _resolve_context_handle(context, model),
    )


def placeholder_model_summary() -> ModelSummary:
    return _raw.placeholder_model_summary()


def primary_outer_boundary_loop(boundary: FaceBoundaryLoops):
    return _raw.primary_outer_boundary_loop(boundary)
