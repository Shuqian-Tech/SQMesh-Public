from __future__ import annotations

import weakref

from . import _sqmesh

_raw = _sqmesh.base

HandleKind = _raw.HandleKind
StatusCode = _raw.StatusCode

SQMeshError = _sqmesh.SQMeshError
InvalidArgumentError = _sqmesh.InvalidArgumentError
NotInitializedError = _sqmesh.NotInitializedError
InvalidHandleError = _sqmesh.InvalidHandleError
OwnerMismatchError = _sqmesh.OwnerMismatchError
IoError = _sqmesh.IoError
UnsupportedError = _sqmesh.UnsupportedError
InternalError = _sqmesh.InternalError

INVALID_HANDLE = _sqmesh.invalid_handle

_contexts: "weakref.WeakValueDictionary[int, Context]" = weakref.WeakValueDictionary()


def _remember_context(context: "Context") -> "Context":
    if context.handle != INVALID_HANDLE:
        _contexts[context.handle] = context
    return context


def _forget_context(handle: int) -> None:
    context = _contexts.pop(handle, None)
    if context is not None:
        context._invalidate()


def _bind_context_handle(handle: int) -> "Context | None":
    if handle == INVALID_HANDLE:
        return None
    context = _contexts.get(handle)
    if context is not None:
        return context
    return _remember_context(Context._from_handle(handle))


def _coerce_context_handle(context: "Context | None") -> int:
    if context is None:
        return INVALID_HANDLE
    if not isinstance(context, Context):
        raise TypeError("context must be a sqmesh.base.Context or None")
    context._require_open()
    return context.handle


class Context:
    def __init__(self) -> None:
        self._handle = _raw.initialize()
        _remember_context(self)

    @classmethod
    def _from_handle(cls, handle: int) -> "Context":
        context = cls.__new__(cls)
        context._handle = handle
        return context

    @property
    def handle(self) -> int:
        return self._handle

    @property
    def closed(self) -> bool:
        return self._handle == INVALID_HANDLE

    def _require_open(self) -> None:
        if self.closed:
            raise InvalidHandleError("Context handle is closed or no longer valid.")

    def _invalidate(self) -> None:
        self._handle = INVALID_HANDLE

    def close(self) -> None:
        if self._handle == INVALID_HANDLE:
            return
        handle = self._handle
        _raw.shutdown(handle)
        self._invalidate()
        _contexts.pop(handle, None)

    def __enter__(self) -> "Context":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    def __repr__(self) -> str:
        return f"<sqmesh.base.Context handle={self.handle}>"


def module_name() -> str:
    return _raw.module_name()


def status_code_name(code: StatusCode) -> str:
    return _raw.status_code_name(code)


def initialize() -> Context:
    return Context()


def shutdown(context: Context | None = None) -> None:
    if context is not None:
        context.close()
        return
    handle = _raw.current_context()
    _raw.shutdown()
    if handle != INVALID_HANDLE:
        _forget_context(handle)


def shutdown_all() -> None:
    _raw.shutdown_all()
    for context in list(_contexts.values()):
        context._invalidate()
    _contexts.clear()


def is_initialized() -> bool:
    return _raw.is_initialized()


def current_context() -> Context | None:
    return _bind_context_handle(_raw.current_context())


def current_session(context: Context | None = None) -> int:
    return _raw.current_session(_coerce_context_handle(context))


def last_error_code() -> StatusCode:
    return _raw.last_error_code()


def last_error_message() -> str:
    return _raw.last_error_message()
