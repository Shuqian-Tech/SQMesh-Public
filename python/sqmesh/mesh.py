from __future__ import annotations

from . import _sqmesh
from . import base

_raw = _sqmesh.mesh

EntityOrder = _raw.EntityOrder
FaceSide = _raw.FaceSide
EntityKind = _raw.EntityKind
ParameterType = _raw.ParameterType
QualityStatus = _raw.QualityStatus
EntityGroupSemantic = _raw.EntityGroupSemantic
EntityGroupRole = _raw.EntityGroupRole
EntityGroupImportFormat = _raw.EntityGroupImportFormat
NastranEntityGroupSourceCard = _raw.NastranEntityGroupSourceCard
NastranPropertyCard = _raw.NastranPropertyCard
NastranMaterialCard = _raw.NastranMaterialCard
MshFormatVersion = _raw.MshFormatVersion
NastranFieldFormat = _raw.NastranFieldFormat

EntityRef = _raw.EntityRef
ConnectivitySpan = _raw.ConnectivitySpan
EntityHeader = _raw.EntityHeader
Node = _raw.Node
Edge = _raw.Edge
Face = _raw.Face
Cell = _raw.Cell
MeshCoreLayout = _raw.MeshCoreLayout
MeshSummary = _raw.MeshSummary
ParameterValue = _raw.ParameterValue
ParameterDictionary = _raw.ParameterDictionary
MeshSizeControl = _raw.MeshSizeControl
MeshSizeControls = _raw.MeshSizeControls
MeshingOptions = _raw.MeshingOptions
DomainStatistics = _raw.DomainStatistics
QualityMetricSummary = _raw.QualityMetricSummary
ElementQuality = _raw.ElementQuality
KindQualitySummary = _raw.KindQualitySummary
MeshQualityReport = _raw.MeshQualityReport
CgnsEntityGroupImportInfo = _raw.CgnsEntityGroupImportInfo
NastranEntityGroupImportInfo = _raw.NastranEntityGroupImportInfo
EntityGroupImportInfo = _raw.EntityGroupImportInfo
EntityGroupInfo = _raw.EntityGroupInfo
EntityGroup = _raw.EntityGroup
Domain = _raw.Domain
MshImportOptions = _raw.MshImportOptions
MshExportOptions = _raw.MshExportOptions
ObjImportOptions = _raw.ObjImportOptions
ObjExportOptions = _raw.ObjExportOptions
CgnsImportOptions = _raw.CgnsImportOptions
CgnsExportOptions = _raw.CgnsExportOptions
NastranImportOptions = _raw.NastranImportOptions
NastranExportOptions = _raw.NastranExportOptions

INVALID_HANDLE = _sqmesh.invalid_handle
INVALID_INDEX = _raw.invalid_index


def _resolve_context_handle(context=None, owner=None) -> int:
    if context is not None:
        return base._coerce_context_handle(context)
    if owner is not None and owner.context is not None:
        return base._coerce_context_handle(owner.context)
    return INVALID_HANDLE


def _resolve_context_object(context=None, owner=None):
    if context is not None:
        return context
    if owner is not None:
        return owner.context
    return base.current_context()


def _coerce_mesh_handle(mesh: "Mesh") -> int:
    if not isinstance(mesh, Mesh):
        raise TypeError("mesh must be a sqmesh.mesh.Mesh")
    mesh._require_live_context()
    return mesh.handle


def _coerce_parameter_value(value):
    if isinstance(value, bool):
        return ParameterValue(value)
    if isinstance(value, int):
        return ParameterValue(value)
    if isinstance(value, float):
        return ParameterValue(value)
    if isinstance(value, str):
        return ParameterValue(value)
    if isinstance(value, ParameterValue):
        return value
    raise TypeError("parameter values must be bool, int, float, str, or ParameterValue")


def _coerce_parameters(parameters=None) -> ParameterDictionary:
    if parameters is None:
        return ParameterDictionary()
    if isinstance(parameters, ParameterDictionary):
        return parameters
    if not hasattr(parameters, "items"):
        raise TypeError("parameters must be a dict-like object or ParameterDictionary")
    dictionary = ParameterDictionary()
    for key, value in parameters.items():
        if not isinstance(key, str):
            raise TypeError("parameter keys must be strings")
        dictionary.set(key, _coerce_parameter_value(value))
    return dictionary


def _coerce_size_controls(size_controls=None) -> MeshSizeControls:
    if size_controls is None:
        return MeshSizeControls()
    if isinstance(size_controls, MeshSizeControls):
        return size_controls
    controls = MeshSizeControls()
    if hasattr(size_controls, "items"):
        for entity, target_size in size_controls.items():
            controls.add_local_size(entity, float(target_size))
        return controls
    for entry in size_controls:
        if isinstance(entry, MeshSizeControl):
            controls.add_local_size(entry.entity, entry.target_size)
            continue
        entity, target_size = entry
        controls.add_local_size(entity, float(target_size))
    return controls


def _build_options(option_type, options=None, **kwargs):
    if options is None:
        options = option_type()
    elif not isinstance(options, option_type):
        raise TypeError(f"options must be a {option_type.__name__}")
    for key, value in kwargs.items():
        setattr(options, key, value)
    return options


def _coerce_meshing_input(
    *,
    parameters=None,
    size_controls=None,
    options=None,
):
    if options is not None:
        if not isinstance(options, MeshingOptions):
            raise TypeError("options must be a MeshingOptions")
        if parameters is not None or size_controls is not None:
            raise TypeError("options cannot be combined with parameters or size_controls")
        return options
    if size_controls is None:
        return _coerce_parameters(parameters)
    opts = MeshingOptions()
    opts.parameters = _coerce_parameters(parameters)
    opts.size_controls = _coerce_size_controls(size_controls)
    return opts


class Mesh:
    def __init__(self, handle: int, context=None) -> None:
        self._handle = handle
        self._context = context

    @property
    def handle(self) -> int:
        return self._handle

    @property
    def context(self):
        return self._context

    def _require_live_context(self) -> None:
        if self._context is not None:
            self._context._require_open()

    def summary(self) -> MeshSummary:
        return mesh_summary(self)

    def quality_report(self) -> MeshQualityReport:
        return mesh_quality_report(self)

    def domain_snapshot(self) -> Domain:
        return domain_snapshot(self)

    def nodes_count(self) -> int:
        return nodes_count(self)

    def cells_count(self) -> int:
        return cells_count(self)

    def export_msh(
        self,
        path: str,
        *,
        context=None,
        options: MshExportOptions | None = None,
        **kwargs,
    ) -> None:
        export_msh(self, path, context=context, options=options, **kwargs)

    def export_obj(
        self,
        path: str,
        *,
        context=None,
        options: ObjExportOptions | None = None,
        **kwargs,
    ) -> None:
        export_obj(self, path, context=context, options=options, **kwargs)

    def export_cgns(
        self,
        path: str,
        *,
        context=None,
        options: CgnsExportOptions | None = None,
        **kwargs,
    ) -> None:
        export_cgns(self, path, context=context, options=options, **kwargs)

    def export_nastran(
        self,
        path: str,
        *,
        context=None,
        options: NastranExportOptions | None = None,
        **kwargs,
    ) -> None:
        export_nastran(self, path, context=context, options=options, **kwargs)

    def __repr__(self) -> str:
        return f"<sqmesh.mesh.Mesh handle={self.handle}>"


def module_name() -> str:
    return _raw.module_name()


def make_dummy_tetra_domain() -> Domain:
    return _raw.make_dummy_tetra_domain()


def encode_entity_kind(order: EntityOrder, arity: int) -> int:
    return _raw.encode_entity_kind(order, arity)


def entity_order_from_kind(kind: EntityKind) -> EntityOrder:
    return _raw.entity_order_from_kind(kind)


def entity_arity(kind: EntityKind) -> int:
    return _raw.entity_arity(kind)


def is_boundary_entity(flags: int) -> bool:
    return _raw.is_boundary_entity(flags)


def is_interface_entity(flags: int) -> bool:
    return _raw.is_interface_entity(flags)


def entity_order_from_flags(flags: int) -> EntityOrder:
    return _raw.entity_order_from_flags(flags)


def mesh_core_layout() -> MeshCoreLayout:
    return _raw.mesh_core_layout()


def expected_mesh_core_layout() -> MeshCoreLayout:
    return _raw.expected_mesh_core_layout()


def same_mesh_core_layout(lhs: MeshCoreLayout, rhs: MeshCoreLayout) -> bool:
    return _raw.same_mesh_core_layout(lhs, rhs)


def matches_mesh_core_layout_baseline(layout: MeshCoreLayout) -> bool:
    return _raw.matches_mesh_core_layout_baseline(layout)


def layout_benchmark_bytes_per_cell(layout: MeshCoreLayout) -> int:
    return _raw.layout_benchmark_bytes_per_cell(layout)


def expected_layout_benchmark_bytes_per_cell() -> int:
    return _raw.expected_layout_benchmark_bytes_per_cell()


def is_algorithm_registered(algorithm_name: str) -> bool:
    return _raw.is_algorithm_registered(algorithm_name)


def cgns_io_available() -> bool:
    return _raw.cgns_io_available()


def import_msh(
    path: str,
    *,
    context=None,
    options: MshImportOptions | None = None,
    **kwargs,
) -> Mesh:
    resolved = _resolve_context_object(context)
    handle = _raw.import_msh(
        path,
        _build_options(MshImportOptions, options, **kwargs),
        _resolve_context_handle(context),
    )
    return Mesh(handle, resolved)


def export_msh(
    mesh: Mesh,
    path: str,
    *,
    context=None,
    options: MshExportOptions | None = None,
    **kwargs,
) -> None:
    _raw.export_msh(
        _coerce_mesh_handle(mesh),
        path,
        _build_options(MshExportOptions, options, **kwargs),
        _resolve_context_handle(context, mesh),
    )


def import_obj(
    path: str,
    *,
    context=None,
    options: ObjImportOptions | None = None,
    **kwargs,
) -> Mesh:
    resolved = _resolve_context_object(context)
    handle = _raw.import_obj(
        path,
        _build_options(ObjImportOptions, options, **kwargs),
        _resolve_context_handle(context),
    )
    return Mesh(handle, resolved)


def export_obj(
    mesh: Mesh,
    path: str,
    *,
    context=None,
    options: ObjExportOptions | None = None,
    **kwargs,
) -> None:
    _raw.export_obj(
        _coerce_mesh_handle(mesh),
        path,
        _build_options(ObjExportOptions, options, **kwargs),
        _resolve_context_handle(context, mesh),
    )


def import_cgns(
    path: str,
    *,
    context=None,
    options: CgnsImportOptions | None = None,
    **kwargs,
) -> Mesh:
    resolved = _resolve_context_object(context)
    handle = _raw.import_cgns(
        path,
        _build_options(CgnsImportOptions, options, **kwargs),
        _resolve_context_handle(context),
    )
    return Mesh(handle, resolved)


def export_cgns(
    mesh: Mesh,
    path: str,
    *,
    context=None,
    options: CgnsExportOptions | None = None,
    **kwargs,
) -> None:
    _raw.export_cgns(
        _coerce_mesh_handle(mesh),
        path,
        _build_options(CgnsExportOptions, options, **kwargs),
        _resolve_context_handle(context, mesh),
    )


def import_nastran(
    path: str,
    *,
    context=None,
    options: NastranImportOptions | None = None,
    **kwargs,
) -> Mesh:
    resolved = _resolve_context_object(context)
    handle = _raw.import_nastran(
        path,
        _build_options(NastranImportOptions, options, **kwargs),
        _resolve_context_handle(context),
    )
    return Mesh(handle, resolved)


def export_nastran(
    mesh: Mesh,
    path: str,
    *,
    context=None,
    options: NastranExportOptions | None = None,
    **kwargs,
) -> None:
    _raw.export_nastran(
        _coerce_mesh_handle(mesh),
        path,
        _build_options(NastranExportOptions, options, **kwargs),
        _resolve_context_handle(context, mesh),
    )


def create_surface_mesh(
    model,
    algorithm_name: str,
    *,
    context=None,
    parameters=None,
    size_controls=None,
    options: MeshingOptions | None = None,
) -> Mesh:
    from . import geo

    if not isinstance(model, geo.Model):
        raise TypeError("model must be a sqmesh.geo.Model")
    resolved = _resolve_context_object(context, model)
    meshing_input = _coerce_meshing_input(
        parameters=parameters,
        size_controls=size_controls,
        options=options,
    )
    handle = _raw.create_surface_mesh(
        geo._coerce_model_handle(model),
        algorithm_name,
        meshing_input,
        _resolve_context_handle(context, model),
    )
    return Mesh(handle, resolved)


def create_volume_mesh(
    model,
    algorithm_name: str,
    *,
    context=None,
    parameters=None,
    size_controls=None,
    options: MeshingOptions | None = None,
) -> Mesh:
    from . import geo

    if not isinstance(model, geo.Model):
        raise TypeError("model must be a sqmesh.geo.Model")
    resolved = _resolve_context_object(context, model)
    meshing_input = _coerce_meshing_input(
        parameters=parameters,
        size_controls=size_controls,
        options=options,
    )
    handle = _raw.create_volume_mesh(
        geo._coerce_model_handle(model),
        algorithm_name,
        meshing_input,
        _resolve_context_handle(context, model),
    )
    return Mesh(handle, resolved)


def mesh_summary(mesh: Mesh, *, context=None) -> MeshSummary:
    return _raw.mesh_summary(
        _coerce_mesh_handle(mesh),
        _resolve_context_handle(context, mesh),
    )


def mesh_quality_report(mesh: Mesh, *, context=None) -> MeshQualityReport:
    return _raw.mesh_quality_report(
        _coerce_mesh_handle(mesh),
        _resolve_context_handle(context, mesh),
    )


def domain_snapshot(mesh: Mesh, *, context=None) -> Domain:
    return _raw.domain_snapshot(
        _coerce_mesh_handle(mesh),
        _resolve_context_handle(context, mesh),
    )


def nodes_count(mesh: Mesh, *, context=None) -> int:
    return _raw.nodes_count(
        _coerce_mesh_handle(mesh),
        _resolve_context_handle(context, mesh),
    )


def cells_count(mesh: Mesh, *, context=None) -> int:
    return _raw.cells_count(
        _coerce_mesh_handle(mesh),
        _resolve_context_handle(context, mesh),
    )
