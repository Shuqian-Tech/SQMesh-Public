# SQMesh Third-Party Licenses

SQMesh is released under **GNU Affero General Public License v3.0 or later**
(see `LICENSE`).  This document summarizes the licenses of the third-party
software that SQMesh incorporates, derives from, or interoperates with.

A high-level acknowledgements section is also maintained in `NOTICE`.

---

## 1. TetGen 1.6

- **Component**: Tetrahedral meshing core (derivative work)
- **Location in SQMesh**: `src/mesh/tet/core/` (entire directory; Shewchuk's predicates have been promoted to `src/core/predicates.cxx` and are shared with auto_cfd)
- **Upstream author**: Hang Si, Weierstrass Institute for Applied Analysis and Stochastics (WIAS), Berlin
- **Upstream homepage**: https://wias-berlin.de/software/tetgen/
- **License**: GNU Affero General Public License v3.0
- **Upstream license text**: https://wias-berlin.de/software/tetgen/1.5/doc/manual/manual002.html

SQMesh's tetrahedral meshing core is a derivative work of TetGen 1.6.  SQMesh
has refactored the original single-file source into topic-sorted translation
units, renamed identifiers to SQMesh naming conventions, and integrated the
core with SQMesh's Domain / EntityGroup / Entity data model.  All such
modifications remain governed by the AGPL-3.0 terms inherited from TetGen 1.6.

Per the AGPL-3.0 license, any redistribution of SQMesh (source or binary,
including network-service distribution under §13) must preserve the license
and make corresponding source code available.

## 2. Shewchuk's Robust Geometric Predicates

- **Component**: Exact geometric predicates (`orient2d`, `orient3d`, `incircle`, `insphere`, etc.)
- **Location in SQMesh**: `src/core/predicates.cxx` (declared in `src/core/predicates.hpp`)
- **Upstream author**: Jonathan Richard Shewchuk, Carnegie Mellon University
- **Upstream homepage**: https://www.cs.cmu.edu/~quake/robust.html
- **License**: Public domain

This file is included in SQMesh verbatim and carries no copyright encumbrance
of its own.  The surrounding SQMesh project is AGPL-3.0, but this file alone
may be used without any restriction.  It is shared between the tetrahedral
mesher (`src/mesh/tet/core/`) and the auto-CFD surface mesher
(`src/mesh/auto_cfd/`).

## 3. Gmsh 4.x

- **Component**: MSH file format reference only — no source code incorporated
- **Location in SQMesh**: Format knowledge encoded in `src/mesh/io/mesh_io.cpp`
- **Upstream authors**: Christophe Geuzaine, Jean-François Remacle, and Gmsh contributors
- **Upstream homepage**: https://gmsh.info/
- **License**: GNU General Public License v2.0 or later, with linking exception
- **Upstream license text**: https://gitlab.onelab.info/gmsh/gmsh/-/blob/master/LICENSE.txt

SQMesh implements MSH 2 / MSH 4 file I/O against the publicly documented
format specification.  No Gmsh source code has been copied into SQMesh.

## 4. OpenCASCADE Technology (OCCT)

- **Component**: Optional CAD kernel link dependency
- **Location in SQMesh**: Linked via adapter in `src/cad/occ/` when `SQMESH_ENABLE_OCC=ON`
- **Upstream vendor**: Open CASCADE S.A.S.
- **Upstream homepage**: https://www.opencascade.com/
- **License**: GNU Lesser General Public License v2.1 (with exception)

OCCT is an optional runtime dependency.  Users who build SQMesh with OCCT
support and redistribute the resulting artifacts must comply with OCCT's
LGPL-2.1 terms independently of SQMesh's AGPL-3.0 terms.

## 5. CGNS

- **Component**: Optional mesh I/O dependency
- **Location in SQMesh**: Linked via `src/mesh/io/` when `SQMESH_ENABLE_CGNS=ON`
- **Upstream homepage**: https://cgns.github.io/
- **License**: zlib-style (very permissive)

## 6. spdlog

- **Component**: Logging backend
- **Location in SQMesh**: Consumed as public headers in `src/core/log.hpp`
- **Upstream author**: Gabi Melman
- **Upstream homepage**: https://github.com/gabime/spdlog
- **License**: MIT License

spdlog is used via its public headers with no modifications.

## 7. fmt (transitive through spdlog)

- **Component**: Formatting library used by spdlog
- **Upstream author**: Victor Zverovich
- **License**: MIT License (with optional exception)

## 8. Python / pybind11 (optional)

- **Component**: Python binding infrastructure, only when `SQMESH_BUILD_PYTHON_BINDINGS=ON`
- **Upstream author**: Wenzel Jakob and contributors
- **License**: BSD 3-Clause

---

## Compatibility Summary

| Component                | License                | Compatible with SQMesh AGPL-3.0? |
|--------------------------|------------------------|----------------------------------|
| TetGen 1.6               | AGPL-3.0               | Yes — same license               |
| Shewchuk's predicates    | Public domain          | Yes                              |
| Gmsh (not incorporated)  | GPL-2.0+               | Not linked; reference only       |
| OCCT                     | LGPL-2.1 + exception   | Yes (dynamic link permitted)     |
| CGNS                     | zlib-style             | Yes                              |
| spdlog                   | MIT                    | Yes                              |
| fmt                      | MIT                    | Yes                              |
| pybind11                 | BSD 3-Clause           | Yes                              |

---

## Network Service Use (AGPL-3.0 §13)

Because SQMesh inherits AGPL-3.0 from TetGen 1.6, operators of network
services that let users interact with SQMesh (directly or through a
derivative product) must offer the corresponding source code of the
running version.  See AGPL-3.0 §13 for the exact obligation.

## Commercial Licensing

If AGPL-3.0 terms are unsuitable for a particular deployment, contact
Hang Si / WIAS Berlin to obtain a commercial license for TetGen, and
separately contact the SQMesh maintainers (Suzhou AI Lab & Shuqian Tech)
to arrange a commercial license for SQMesh's own code.
