// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
// This file is a derivative work of TetGen 1.6 by Hang Si
// (Weierstrass Institute for Applied Analysis and Stochastics, Berlin).
// TetGen 1.6 is distributed under the GNU Affero General Public License v3.0;
// SQMesh inherits and carries forward those terms.
// Upstream: https://wias-berlin.de/software/tetgen/
//
// SQMesh tet I/O carrier — extracted from upstream `tet_core.hpp` so that
// SQMesh-side refactoring of the tet pipeline can iterate independently
// of the legacy header. This file deliberately preserves the upstream
// class name `TetMeshData`, member layout and lifetime semantics so that
// tet core topic files and `tet_native_algorithm.cpp` compile unchanged
// against this declaration.
//
// `class TetMeshData` is the I/O contract between SQMesh and the
// `TetMeshCore` algorithm core: SQMesh marshals surface mesh data into
// the input pointlist/facetlist arrays before invoking
// `run_tet_mesh_core(...)`, then reads tetrahedronlist/pointlist back out
// into a SQMesh `Domain`. The standalone-CLI file readers/writers that
// historically lived on this class were dropped in earlier refactor
// steps; what remains is purely the data carrier plus init/cleanup.

#pragma once

#include <stdint.h>

namespace sqmesh::mesh::tet::detail {

class TetMeshData {

public:

  // A "polygon" describes a simple polygon (no holes). It is not necessarily
  //   convex. Each polygon contains a number of corners (points) and the same
  //   number of sides (edges).  The points of the polygon must be given in
  //   either counterclockwise or clockwise order and they form a ring, so
  //   every two consecutive points forms an edge of the polygon.
  typedef struct {
    int *vertexlist;
    int numberofvertices;
  } polygon;

  // A "facet" describes a polygonal region possibly with holes, edges, and
  //   points floating in it.  Each facet consists of a list of polygons and
  //   a list of hole points (which lie strictly inside holes).
  typedef struct {
    polygon *polygonlist;
    int numberofpolygons;
    double *holelist;
    int numberofholes;
  } facet;

  // A "voroedge" is an edge of the Voronoi diagram. It corresponds to a
  //   Delaunay face.  Each voroedge is either a line segment connecting
  //   two Voronoi vertices or a ray starting from a Voronoi vertex to an
  //   "infinite vertex".  'v1' and 'v2' are two indices pointing to the
  //   list of Voronoi vertices. 'v1' must be non-negative, while 'v2' may
  //   be -1 if it is a ray, in this case, the unit normal of this ray is
  //   given in 'vnormal'.
  typedef struct {
    int v1, v2;
    double vnormal[3];
  } voroedge;

  // A "vorofacet" is an facet of the Voronoi diagram. It corresponds to a
  //   Delaunay edge.  Each Voronoi facet is a convex polygon formed by a
  //   list of Voronoi edges, it may not be closed.  'c1' and 'c2' are two
  //   indices pointing into the list of Voronoi cells, i.e., the two cells
  //   share this facet.  'elist' is an array of indices pointing into the
  //   list of Voronoi edges, 'elist[0]' saves the number of Voronoi edges
  //   (including rays) of this facet.
  typedef struct {
    int c1, c2;
    int *elist;
  } vorofacet;


  // Additional parameters associated with an input (or mesh) vertex.
  //   These informations are provided by CAD libraries.
  typedef struct {
    double uv[2];
    int tag;
    int type; // 0, 1, or 2.
  } pointparam;

  // Callback functions for meshing PSCs.
  typedef double (* GetVertexParamOnEdge)(void*, int, int);
  typedef void (* GetSteinerOnEdge)(void*, int, double, double*);
  typedef void (* GetVertexParamOnFace)(void*, int, int, double*);
  typedef void (* GetEdgeSteinerParamOnFace)(void*, int, double, int, double*);
  typedef void (* GetSteinerOnFace)(void*, int, double*, double*);

  // A callback function for mesh refinement.
  typedef bool (* TetSizeFunc)(double*, double*, double*, double*, double*, double);

  // Items are numbered starting from 'firstnumber' (0 or 1), default is 0.
  int firstnumber;

  // Dimension of the mesh (2 or 3), default is 3.
  int mesh_dim;

  // Does the lines in .node file contain index or not, default is 1.
  int useindex;

  // 'pointlist':  An array of point coordinates.  The first point's x
  //   coordinate is at index [0] and its y coordinate at index [1], its
  //   z coordinate is at index [2], followed by the coordinates of the
  //   remaining points.  Each point occupies three REALs.
  // 'pointattributelist':  An array of point attributes.  Each point's
  //   attributes occupy 'numberofpointattributes' REALs.
  // 'pointmtrlist': An array of metric tensors at points. Each point's
  //   tensor occupies 'numberofpointmtr' REALs.
  // 'pointmarkerlist':  An array of point markers; one integer per point.
  // 'point2tetlist': An array of tetrahedra indices; one integer per point.
  double *pointlist;
  double *pointattributelist;
  double *pointmtrlist;
  int  *pointmarkerlist;
  int  *point2tetlist;
  pointparam *pointparamlist;
  int numberofpoints;
  int numberofpointattributes;
  int numberofpointmtrs;

  // 'tetrahedronlist':  An array of tetrahedron corners.  The first
  //   tetrahedron's first corner is at index [0], followed by its other
  //   corners, followed by six nodes on the edges of the tetrahedron if the
  //   second order option (-o2) is applied. Each tetrahedron occupies
  //   'numberofcorners' ints.  The second order nodes are ouput only.
  // 'tetrahedronattributelist':  An array of tetrahedron attributes.  Each
  //   tetrahedron's attributes occupy 'numberoftetrahedronattributes' REALs.
  // 'tetrahedronvolumelist':  An array of constraints, i.e. tetrahedron's
  //   volume; one double per element.  Input only.
  // 'neighborlist':  An array of tetrahedron neighbors; 4 ints per element.
  // 'tet2facelist':  An array of tetrahedron face indices; 4 ints per element.
  // 'tet2edgelist':  An array of tetrahedron edge indices; 6 ints per element.
  int  *tetrahedronlist;
  double *tetrahedronattributelist;
  double *tetrahedronvolumelist;
  int  *neighborlist;
  int  *tet2facelist;
  int  *tet2edgelist;
  int numberoftetrahedra;
  int numberofcorners;
  int numberoftetrahedronattributes;

  // 'facetlist':  An array of facets.  Each entry is a structure of facet.
  // 'facetmarkerlist':  An array of facet markers; one int per facet.
  facet *facetlist;
  int *facetmarkerlist;
  int numberoffacets;

  // 'holelist':  An array of holes (in volume).  Each hole is given by a
  //   seed (point) which lies strictly inside it. The first seed's x, y and z
  //   coordinates are at indices [0], [1] and [2], followed by the
  //   remaining seeds.  Three REALs per hole.
  double *holelist;
  int numberofholes;

  // 'regionlist': An array of regions (subdomains).  Each region is given by
  //   a seed (point) which lies strictly inside it. The first seed's x, y and
  //   z coordinates are at indices [0], [1] and [2], followed by the regional
  //   attribute at index [3], followed by the maximum volume at index [4].
  //   Five REALs per region.
  // Note that each regional attribute is used only if you select the 'A'
  //   switch, and each volume constraint is used only if you select the
  //   'a' switch (with no number following).
  double *regionlist;
  int numberofregions;

  // 'refine_elem_list': An array of tetrahedra to be refined.  The first
  //   tetrahedron's first corner is at index [0], followed by its other
  //   corners. Four integers per element.
  // 'refine_elem_vol_list':  An array of constraints, i.e. tetrahedron's
  //   volume; one double per element.
  int  *refine_elem_list;
  double *refine_elem_vol_list;
  int  numberofrefineelems;

  // 'facetconstraintlist':  An array of facet constraints.  Each constraint
  //   specifies a maximum area bound on the subfaces of that facet.  The
  //   first facet constraint is given by a facet marker at index [0] and its
  //   maximum area bound at index [1], followed by the remaining facet con-
  //   straints. Two REALs per facet constraint.  Note: the facet marker is
  //   actually an integer.
  double *facetconstraintlist;
  int numberoffacetconstraints;

  // 'segmentconstraintlist': An array of segment constraints. Each constraint
  //   specifies a maximum length bound on the subsegments of that segment.
  //   The first constraint is given by the two endpoints of the segment at
  //   index [0] and [1], and the maximum length bound at index [2], followed
  //   by the remaining segment constraints.  Three REALs per constraint.
  //   Note the segment endpoints are actually integers.
  double *segmentconstraintlist;
  int numberofsegmentconstraints;


  // 'trifacelist':  An array of face (triangle) corners.  The first face's
  //   three corners are at indices [0], [1] and [2], followed by the remaining
  //   faces.  Three ints per face.
  // 'trifacemarkerlist':  An array of face markers; one int per face.
  // 'o2facelist':  An array of second order nodes (on the edges) of the face.
  //   It is output only if the second order option (-o2) is applied. The
  //   first face's three second order nodes are at [0], [1], and [2],
  //   followed by the remaining faces.  Three ints per face.
  // 'face2tetlist':  An array of tetrahedra indices; 2 ints per face.
  // 'face2edgelist':  An array of edge indices; 3 ints per face.
  int *trifacelist;
  int *trifacemarkerlist;
  int *o2facelist;
  int *face2tetlist;
  int *face2edgelist;
  int numberoftrifaces;

  // 'edgelist':  An array of edge endpoints.  The first edge's endpoints
  //   are at indices [0] and [1], followed by the remaining edges.
  //   Two ints per edge.
  // 'edgemarkerlist':  An array of edge markers; one int per edge.
  // 'o2edgelist':  An array of midpoints of edges. It is output only if the
  //   second order option (-o2) is applied. One int per edge.
  // 'edge2tetlist':  An array of tetrahedra indices.  One int per edge.
  int *edgelist;
  int *edgemarkerlist;
  int *o2edgelist;
  int *edge2tetlist;
  int numberofedges;

  // 'vpointlist':  An array of Voronoi vertex coordinates (like pointlist).
  // 'vedgelist':  An array of Voronoi edges.  Each entry is a 'voroedge'.
  // 'vfacetlist':  An array of Voronoi facets. Each entry is a 'vorofacet'.
  // 'vcelllist':  An array of Voronoi cells.  Each entry is an array of
  //   indices pointing into 'vfacetlist'. The 0th entry is used to store
  //   the length of this array.
  double *vpointlist;
  voroedge *vedgelist;
  vorofacet *vfacetlist;
  int **vcelllist;
  int numberofvpoints;
  int numberofvedges;
  int numberofvfacets;
  int numberofvcells;


  // Variable (and callback functions) for meshing PSCs.
  void *geomhandle;
  GetVertexParamOnEdge getvertexparamonedge;
  GetSteinerOnEdge getsteineronedge;
  GetVertexParamOnFace getvertexparamonface;
  GetEdgeSteinerParamOnFace getedgesteinerparamonface;
  GetSteinerOnFace getsteineronface;

  // A callback function.
  TetSizeFunc tetunsuitable;

  // File I/O routines (load_*/save_* + readline helpers) have been removed.
  // SQMesh marshals data directly into `pointlist`/`facetlist` via
  // `tet_native_algorithm.cpp`; the standalone CLI file readers/writers
  // are dead code in library mode.

  static void init(polygon* p) {
    p->vertexlist = (int *) NULL;
    p->numberofvertices = 0;
  }

  static void init(facet* f) {
    f->polygonlist = (polygon *) NULL;
    f->numberofpolygons = 0;
    f->holelist = (double *) NULL;
    f->numberofholes = 0;
  }

  // Initialize routine.
  void initialize()
  {
    firstnumber = 0;
    mesh_dim = 3;
    useindex = 1;

    pointlist = (double *) NULL;
    pointattributelist = (double *) NULL;
    pointmtrlist = (double *) NULL;
    pointmarkerlist = (int *) NULL;
	point2tetlist = (int *) NULL;
    pointparamlist = (pointparam *) NULL;
    numberofpoints = 0;
    numberofpointattributes = 0;
    numberofpointmtrs = 0;

    tetrahedronlist = (int *) NULL;
    tetrahedronattributelist = (double *) NULL;
    tetrahedronvolumelist = (double *) NULL;
    neighborlist = (int *) NULL;
	tet2facelist = (int *) NULL;
	tet2edgelist = (int *) NULL;
    numberoftetrahedra = 0;
    numberofcorners = 4;
    numberoftetrahedronattributes = 0;

    trifacelist = (int *) NULL;
    trifacemarkerlist = (int *) NULL;
    o2facelist = (int *) NULL;
    face2tetlist = (int *) NULL;
	face2edgelist = (int *) NULL;
    numberoftrifaces = 0;

    edgelist = (int *) NULL;
    edgemarkerlist = (int *) NULL;
    o2edgelist = (int *) NULL;
    edge2tetlist = (int *) NULL;
    numberofedges = 0;

    facetlist = (facet *) NULL;
    facetmarkerlist = (int *) NULL;
    numberoffacets = 0;

    holelist = (double *) NULL;
    numberofholes = 0;

    regionlist = (double *) NULL;
    numberofregions = 0;

    refine_elem_list = (int *) NULL;
    refine_elem_vol_list = (double *) NULL;
    numberofrefineelems = 0;

    facetconstraintlist = (double *) NULL;
    numberoffacetconstraints = 0;
    segmentconstraintlist = (double *) NULL;
    numberofsegmentconstraints = 0;


    vpointlist = (double *) NULL;
    vedgelist = (voroedge *) NULL;
    vfacetlist = (vorofacet *) NULL;
    vcelllist = (int **) NULL;
    numberofvpoints = 0;
    numberofvedges = 0;
    numberofvfacets = 0;
    numberofvcells = 0;


    tetunsuitable = NULL;

    geomhandle = NULL;
    getvertexparamonedge = NULL;
    getsteineronedge = NULL;
    getvertexparamonface = NULL;
    getedgesteinerparamonface = NULL;
    getsteineronface = NULL;
  }

  // Free the memory allocated in 'TetMeshData'.  Note that it assumes that the
  //   memory was allocated by the "new" operator (C++).
  void clean_memory()
  {
    int i, j;

    if (pointlist != (double *) NULL) {
      delete [] pointlist;
    }
    if (pointattributelist != (double *) NULL) {
      delete [] pointattributelist;
    }
    if (pointmtrlist != (double *) NULL) {
      delete [] pointmtrlist;
    }
    if (pointmarkerlist != (int *) NULL) {
      delete [] pointmarkerlist;
    }
	if (point2tetlist != (int *) NULL) {
      delete [] point2tetlist;
    }
    if (pointparamlist != (pointparam *) NULL) {
      delete [] pointparamlist;
    }

    if (tetrahedronlist != (int *) NULL) {
      delete [] tetrahedronlist;
    }
    if (tetrahedronattributelist != (double *) NULL) {
      delete [] tetrahedronattributelist;
    }
    if (tetrahedronvolumelist != (double *) NULL) {
      delete [] tetrahedronvolumelist;
    }
    if (neighborlist != (int *) NULL) {
      delete [] neighborlist;
    }
    if (tet2facelist != (int *) NULL) {
	  delete [] tet2facelist;
	}
	if (tet2edgelist != (int *) NULL) {
	  delete [] tet2edgelist;
	}

    if (trifacelist != (int *) NULL) {
      delete [] trifacelist;
    }
    if (trifacemarkerlist != (int *) NULL) {
      delete [] trifacemarkerlist;
    }
    if (o2facelist != (int *) NULL) {
      delete [] o2facelist;
    }
    if (face2tetlist != (int *) NULL) {
      delete [] face2tetlist;
    }
	if (face2edgelist != (int *) NULL) {
      delete [] face2edgelist;
    }

    if (edgelist != (int *) NULL) {
      delete [] edgelist;
    }
    if (edgemarkerlist != (int *) NULL) {
      delete [] edgemarkerlist;
    }
    if (o2edgelist != (int *) NULL) {
      delete [] o2edgelist;
    }
    if (edge2tetlist != (int *) NULL) {
      delete [] edge2tetlist;
    }

    if (facetlist != (facet *) NULL) {
      facet *f;
      polygon *p;
      for (i = 0; i < numberoffacets; i++) {
        f = &facetlist[i];
        for (j = 0; j < f->numberofpolygons; j++) {
          p = &f->polygonlist[j];
          delete [] p->vertexlist;
        }
        delete [] f->polygonlist;
        if (f->holelist != (double *) NULL) {
          delete [] f->holelist;
        }
      }
      delete [] facetlist;
    }
    if (facetmarkerlist != (int *) NULL) {
      delete [] facetmarkerlist;
    }

    if (holelist != (double *) NULL) {
      delete [] holelist;
    }
    if (regionlist != (double *) NULL) {
      delete [] regionlist;
    }
    if (refine_elem_list != (int *) NULL) {
      delete [] refine_elem_list;
      if (refine_elem_vol_list != (double *) NULL) {
        delete [] refine_elem_vol_list;
      }
    }
    if (facetconstraintlist != (double *) NULL) {
      delete [] facetconstraintlist;
    }
    if (segmentconstraintlist != (double *) NULL) {
      delete [] segmentconstraintlist;
    }
    if (vpointlist != (double *) NULL) {
      delete [] vpointlist;
    }
    if (vedgelist != (voroedge *) NULL) {
      delete [] vedgelist;
    }
    if (vfacetlist != (vorofacet *) NULL) {
      for (i = 0; i < numberofvfacets; i++) {
        delete [] vfacetlist[i].elist;
      }
      delete [] vfacetlist;
    }
    if (vcelllist != (int **) NULL) {
      for (i = 0; i < numberofvcells; i++) {
        delete [] vcelllist[i];
      }
      delete [] vcelllist;
    }
  }

  // Constructor & destructor.
  TetMeshData() {initialize();}
  ~TetMeshData() {clean_memory();}

}; // class TetMeshData

} // namespace sqmesh::mesh::tet::detail
