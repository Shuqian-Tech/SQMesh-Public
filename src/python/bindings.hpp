// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
#pragma once

#include <pybind11/pybind11.h>

namespace py = pybind11;

void register_sqmesh_exceptions(py::module_ &module);
void bind_base(py::module_ &module);
void bind_geo(py::module_ &module);
void bind_mesh(py::module_ &module);
