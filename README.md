# SQMesh

**SQMesh** is an open-source meshing toolkit for CFD and FEA workflows,
maintained by Suzhou AI Lab and Shuqian Tech.

## Features

- **Geometry** — STEP / IGES / STL import, OCC-free access layer
  (topology traversal, UV queries, projection, curvature sampling)
- **Surface meshing** — Auto CFD Surface Mesher (sizing-field driven,
  honors `minimum_length` / `maximum_length` / `distortion_angle` /
  `growth_rate` / `proximity`)
- **Volume meshing** — Tetrahedral Volume Mesher
- **Boundary-layer meshing** — Boundary Layer Mesher with prismatic
  inflation for CFD workflows
- **Mesh I/O** — MSH 2.2 / 4.1, OBJ, NASTRAN, optional CGNS
- **C++17 SDK** with **Python bindings** (`sqmesh.base` / `sqmesh.geo` /
  `sqmesh.mesh`)

## Build

### Prerequisites

- CMake ≥ 3.20
- A C++17 compiler (MSVC 2019+ on Windows, GCC 9+ / Clang 10+ on Linux/macOS)
- Python ≥ 3.9 (only if `-DSQMESH_BUILD_PYTHON_BINDINGS=ON`)

### Third-party libraries

| Library | Required? | Used for |
| --- | --- | --- |
| [**OpenCASCADE**](https://dev.opencascade.org/) ≥ 7.5 | optional (recommended) | STEP / IGES import, full CAD-driven meshing pipeline |
| [**spdlog**](https://github.com/gabime/spdlog) | optional | Structured logging (falls back to stderr if absent) |
| [**CGNS**](https://cgns.github.io/) | optional | `.cgns` I/O (enable with `-DSQMESH_ENABLE_CGNS=ON`) |

On Windows the easiest route is [vcpkg](https://github.com/microsoft/vcpkg):

```bash
# one-time vcpkg bootstrap
git clone https://github.com/microsoft/vcpkg.git C:\vcpkg
C:\vcpkg\bootstrap-vcpkg.bat

# install SQMesh's dependencies for the x64 Windows triplet
C:\vcpkg\vcpkg install opencascade:x64-windows spdlog:x64-windows
# (optional) CGNS I/O
C:\vcpkg\vcpkg install cgns:x64-windows
```

On Linux / macOS, use the system package manager or vcpkg classic mode:

```bash
# Ubuntu / Debian
sudo apt install libocct-foundation-dev libocct-modeling-algorithms-dev \
                 libocct-modeling-data-dev libocct-data-exchange-dev \
                 libspdlog-dev libcgns-dev

# macOS (Homebrew)
brew install opencascade spdlog cgns
```

### Configure + build

Pass the vcpkg toolchain so CMake's `find_package(OpenCASCADE)` / `find_package(spdlog)`
resolve against the installed packages:

```bash
cmake -S . -B build_occ ^
  -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake ^
  -DVCPKG_TARGET_TRIPLET=x64-windows ^
  -DSQMESH_ENABLE_OCC=ON ^
  -DSQMESH_BUILD_TESTS=ON ^
  -DCMAKE_BUILD_TYPE=Release
cmake --build build_occ --config Release -j
ctest --test-dir build_occ -C Release --output-on-failure
```

On Linux / macOS drop the `VCPKG_TARGET_TRIPLET` and use forward slashes:

```bash
cmake -S . -B build_occ \
  -DCMAKE_TOOLCHAIN_FILE=$HOME/vcpkg/scripts/buildsystems/vcpkg.cmake \
  -DSQMESH_ENABLE_OCC=ON -DSQMESH_BUILD_TESTS=ON \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build_occ -j
```

If you installed OCCT / spdlog without vcpkg (system packages, manual install),
drop the toolchain line and instead point CMake at the install roots via
`CMAKE_PREFIX_PATH` / `OpenCASCADE_DIR` / `CASROOT`.

### Build options

| Option | Default | Effect |
| --- | --- | --- |
| `SQMESH_ENABLE_OCC` | `OFF` | Enable OpenCASCADE-backed STEP / IGES import. Without it the CAD adapter builds as a stub that returns `unsupported`. |
| `SQMESH_ENABLE_CGNS` | `OFF` | Build the CGNS I/O path. |
| `SQMESH_BUILD_TESTS` | `ON` | Build the `ctest` suite. |
| `SQMESH_BUILD_BENCHMARKS` | `OFF` | Build benchmark placeholders. |
| `SQMESH_BUILD_PYTHON_BINDINGS` | `OFF` | Build `sqmesh.base` / `sqmesh.geo` / `sqmesh.mesh`. |

## Usage

Mesh a STEP / IGES file end to end with the bundled example:

```bash
./build_occ/examples/Release/surface_mesh_example \
    path/to/model.step <min_length> <max_length> [distortion_angle] [growth_rate] [proximity]
```

The output is a surface mesh (`surface_mesh.obj`) in the current directory.
To generate a volume mesh with boundary layers, use `volume_mesh_example`:

```bash
./build_occ/examples/Release/volume_mesh_example \
    path/to/model.step <min_length> <max_length> <distortion_angle> \
    <growth_rate> <proximity> <auto_topo> \
    <bl_first_height> <bl_growth_rate> <bl_num_layers> \
    <tet_max_length> <tet_growth_rate> <material_point> [bl_mode] [bl_aspect]
```

Link against the C++ SDK in your own project by including `<sqmesh/sqmesh.hpp>`
and the `sqmesh` CMake target. See
[`examples/surface_mesh_example.cpp`](examples/surface_mesh_example.cpp) and
[`examples/volume_mesh_example.cpp`](examples/volume_mesh_example.cpp) for
complete references covering context setup, geometry import, surface meshing,
volume meshing, and mesh export.

## License

SQMesh is licensed under **GNU Affero General Public License v3.0 or later**
— see [`LICENSE`](LICENSE). Third-party components and attribution are
documented in [`NOTICE`](NOTICE) and
[`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md).

## Original Contributors

The original contributor list before this public release is preserved below:

![Original SQMesh contributors](docs/assets/original-contributors.png)
