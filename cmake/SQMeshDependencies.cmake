include_guard(GLOBAL)

set(SQMESH_OCC_FOUND OFF)
set(SQMESH_OCC_PACKAGE "")
set(SQMESH_OCC_INCLUDE_DIRS "")
set(SQMESH_OCC_LIBRARIES "")
set(SQMESH_OCC_VERSION "")
set(SQMESH_CGNS_FOUND OFF)
set(SQMESH_CGNS_PACKAGE "")
set(SQMESH_CGNS_INCLUDE_DIRS "")
set(SQMESH_CGNS_LIBRARIES "")
set(SQMESH_CGNS_VERSION "")

# ── spdlog ──────────────────────────────────────────────────────────────
find_package(spdlog CONFIG QUIET)
if(spdlog_FOUND)
  set(SQMESH_SPDLOG_FOUND ON)
  message(STATUS "Found spdlog: ${spdlog_VERSION}")
else()
  set(SQMESH_SPDLOG_FOUND OFF)
  message(STATUS "spdlog not found — logging will fall back to stderr")
endif()

if(SQMESH_ENABLE_OCC)
  find_package(OpenCASCADE QUIET)
  if(OpenCASCADE_FOUND)
    set(SQMESH_OCC_FOUND ON)
    set(SQMESH_OCC_PACKAGE "OpenCASCADE")
    set(SQMESH_OCC_INCLUDE_DIRS "${OpenCASCADE_INCLUDE_DIR}")
    set(
      SQMESH_OCC_LIBRARIES
      ${OpenCASCADE_ModelingData_LIBRARIES}
      ${OpenCASCADE_ModelingAlgorithms_LIBRARIES}
      ${OpenCASCADE_DataExchange_LIBRARIES}
    )
    set(
      SQMESH_OCC_VERSION
      "${OpenCASCADE_MAJOR_VERSION}.${OpenCASCADE_MINOR_VERSION}.${OpenCASCADE_MAINTENANCE_VERSION}"
    )
  else()
    find_package(OCC QUIET)
    if(OCC_FOUND)
      set(SQMESH_OCC_FOUND ON)
      set(SQMESH_OCC_PACKAGE "OCC")
      set(SQMESH_OCC_INCLUDE_DIRS ${OCC_INCLUDE_DIRS} ${OCC_INCLUDE_DIR})
      set(SQMESH_OCC_LIBRARIES ${OCC_LIBRARIES})
    else()
      find_path(
        SQMESH_OCC_MANUAL_INCLUDE_DIR
        NAMES TopoDS_Shape.hxx
        HINTS
          ENV CASROOT
        PATH_SUFFIXES include include/opencascade opencascade
      )

      set(_sqmesh_occ_manual_components
        TKernel
        TKMath
        TKG2d
        TKG3d
        TKGeomBase
        TKGeomAlgo
        TKTopAlgo
        TKBRep
        TKPrim
        TKBO
        TKBool
        TKShHealing
        TKMesh
        TKHLR
        TKOffset
        TKFeat
        TKFillet
        TKService
        TKXSBase
        TKSTEP
        TKSTEPAttr
        TKSTEP209
        TKSTEPBase
        TKIGES
      )
      set(_sqmesh_occ_manual_libraries "")
      set(_sqmesh_occ_manual_missing_component "")
      foreach(_sqmesh_occ_component IN LISTS _sqmesh_occ_manual_components)
        unset(_sqmesh_occ_component_library CACHE)
        unset(_sqmesh_occ_component_library)
        find_library(
          _sqmesh_occ_component_library
          NAMES ${_sqmesh_occ_component}
          HINTS
            ENV CASROOT
          PATH_SUFFIXES lib lib64
        )
        if(NOT _sqmesh_occ_component_library)
          set(_sqmesh_occ_manual_missing_component "${_sqmesh_occ_component}")
          break()
        endif()
        list(APPEND _sqmesh_occ_manual_libraries "${_sqmesh_occ_component_library}")
      endforeach()

      if(SQMESH_OCC_MANUAL_INCLUDE_DIR AND NOT _sqmesh_occ_manual_missing_component)
        set(SQMESH_OCC_FOUND ON)
        set(SQMESH_OCC_PACKAGE "manual-search")
        set(SQMESH_OCC_INCLUDE_DIRS "${SQMESH_OCC_MANUAL_INCLUDE_DIR}")
        set(SQMESH_OCC_LIBRARIES ${_sqmesh_occ_manual_libraries})

        if(EXISTS "${SQMESH_OCC_MANUAL_INCLUDE_DIR}/Standard_Version.hxx")
          file(
            STRINGS
            "${SQMESH_OCC_MANUAL_INCLUDE_DIR}/Standard_Version.hxx"
            _sqmesh_occ_complete_version
            REGEX "^#define[ \t]+OCC_VERSION_COMPLETE[ \t]+\"[0-9.]+\""
          )
          if(_sqmesh_occ_complete_version)
            string(
              REGEX REPLACE
              ".*OCC_VERSION_COMPLETE[ \t]+\"([0-9.]+)\".*"
              "\\1"
              SQMESH_OCC_VERSION
              "${_sqmesh_occ_complete_version}"
            )
          endif()
        endif()
      else()
        message(WARNING
          "SQMESH_ENABLE_OCC=ON but no OpenCASCADE package was found. "
          "The cad/occ target will build as a stub boundary. "
          "A linkable OCCT build must be discoverable via OpenCASCADE_DIR, "
          "CMAKE_PREFIX_PATH, CASROOT, or an equivalent vcpkg toolchain. "
          "Manual OCCT probing also failed"
          "${SQMESH_OCC_MANUAL_INCLUDE_DIR};"
          " missing component: ${_sqmesh_occ_manual_missing_component}."
        )
      endif()
    endif()
  endif()
endif()

if(SQMESH_OCC_LIBRARIES)
  list(REMOVE_DUPLICATES SQMESH_OCC_LIBRARIES)
endif()

if(SQMESH_OCC_INCLUDE_DIRS)
  list(REMOVE_DUPLICATES SQMESH_OCC_INCLUDE_DIRS)
endif()

if(SQMESH_ENABLE_CGNS)
  find_path(
    SQMESH_CGNS_INCLUDE_DIR
    NAMES cgnslib.h
    HINTS
      ENV CGNS_ROOT
      ENV CGNS_DIR
    PATH_SUFFIXES include
  )
  find_library(
    SQMESH_CGNS_LIBRARY
    NAMES cgns libcgns
    HINTS
      ENV CGNS_ROOT
      ENV CGNS_DIR
    PATH_SUFFIXES lib lib64
  )

  if(SQMESH_CGNS_INCLUDE_DIR AND SQMESH_CGNS_LIBRARY)
    set(SQMESH_CGNS_FOUND ON)
    set(SQMESH_CGNS_PACKAGE "manual-search")
    set(SQMESH_CGNS_INCLUDE_DIRS "${SQMESH_CGNS_INCLUDE_DIR}")
    set(SQMESH_CGNS_LIBRARIES "${SQMESH_CGNS_LIBRARY}")

    file(
      STRINGS
      "${SQMESH_CGNS_INCLUDE_DIR}/cgnslib.h"
      _sqmesh_cgns_dotvers
      REGEX "^#define[ \t]+CGNS_DOTVERS[ \t]+[0-9.]+"
    )
    if(_sqmesh_cgns_dotvers)
      string(
        REGEX REPLACE
        ".*CGNS_DOTVERS[ \t]+([0-9.]+).*"
        "\\1"
        SQMESH_CGNS_VERSION
        "${_sqmesh_cgns_dotvers}"
      )
    endif()
  else()
    message(WARNING
      "SQMESH_ENABLE_CGNS=ON but no CGNS development install was found. "
      "The mesh CGNS API will build as an unsupported stub boundary. "
      "Provide cgnslib.h and a linkable cgns library through CMAKE_PREFIX_PATH, "
      "CGNS_ROOT, CGNS_DIR, or equivalent discovery hints to enable the real path."
    )
  endif()
endif()

if(SQMESH_CGNS_LIBRARIES)
  list(REMOVE_DUPLICATES SQMESH_CGNS_LIBRARIES)
endif()

if(SQMESH_CGNS_INCLUDE_DIRS)
  list(REMOVE_DUPLICATES SQMESH_CGNS_INCLUDE_DIRS)
endif()
