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
#ifndef SQMESH_MESH_TET_CORE_HPP
#define SQMESH_MESH_TET_CORE_HPP

// The maximum number of characters in a file name (including the null).

#define FILENAMESIZE 1024

// The maximum number of chars in a line read from a file (including the null).

#define INPUTLINESIZE 2048

// C standard libraries to perform Input/output operations, general utililities,
//   manipulate strings and arrays, compute common mathematical operations,
//   get date and time information.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// The types 'intptr_t' and 'uintptr_t' are signed and unsigned integer types,
//   respectively. They are guaranteed to be the same width as a pointer.
//   They are defined in <stdint.h> by the C99 Standard.

#include <stdint.h>

// `class TetMeshData` has been extracted to `src/mesh/tet/tet_io.hpp` as
// part of the SQMesh-side tet refactor. The declaration is now the
// single source of truth there; this header only re-exposes it so that
// existing translation units (tet core topic files,
// `tet_native_algorithm.cpp`) compile without changes.
#include "../tet_io.hpp"


//============================================================================//
//                                                                            //
// TetMeshBehavior                                                             //
//                                                                            //
// A structure for maintaining the switches and parameters used by the tet core's   //
// internal data structure and algorithms.                                    //
//                                                                            //
// All switches and parameters are initialized with default values. They are  //
// set by the command line arguments (argc, argv).                            //
//                                                                            //
// NOTE: Some switches are incompatible with others. While some may depend    //
// on other switches. The routine parse_commandline() sets the switches from  //
// the command line (a list of strings) and checks the consistency of the     //
// applied switches.                                                          //
//                                                                            //
//============================================================================//

// `class TetMeshBehavior` has been extracted to `src/mesh/tet/tet_behavior.hpp`
// as the first step of the SQMesh-side tet refactor. The declaration is now
// the single source of truth there; this header only re-exposes it so that
// existing translation units (including the upstream tet core topic files) compile
// without changes.
#include "../tet_behavior.hpp"

//============================================================================//
//                                                                            //
// Robust Geometric predicates                                                //
//                                                                            //
// The following routines are the robust geometric predicates for orientation //
// test and point-in-sphere test implemented by Jonathan Shewchuk.            //
// He generously provided the source code in the public domain,               //
// http://www.cs.cmu.edu/~quake/robust.html.                                  //
// predicates.cxx is a C++ version of the original C code, now shared with    //
// auto_cfd from src/core/. See src/core/predicates.hpp.                      //
//                                                                            //
//============================================================================//

#include "../../../core/predicates.hpp"


namespace sqmesh::mesh::tet::detail {

//============================================================================//
//                                                                            //
// TetMeshCore    the tet core's internal mesh data structure.                       //
//                                                                            //
// It uses a tetrahedron-based mesh data structure. It implements elementary  //
// flip operations to locally modify the mesh. It implements basic meshing    //
// algorithms to create Delaunay tetrahedraliations, to perform boundary      //
// recovery, to place Steiner points in the mesh domain, and to optimize the  //
// quality of the mesh.                                                       //
//                                                                            //
//============================================================================//

class TetMeshCore {

public:

//============================================================================//
//                                                                            //
// Mesh data structure                                                        //
//                                                                            //
// A tetrahedral mesh T of a 3D piecewise linear complex (PLC) X is a 3D      //
// simplicial complex whose underlying space is equal to the space of X.  T   //
// contains a 2D subcomplex S which is a triangular mesh of the boundary of   //
// X. S contains a 1D subcomplex L which is a linear mesh of the boundary of  //
// S. Faces and edges in S and L are respectively called subfaces and segme-  //
// nts to distinguish them from others in T.                                  //
//                                                                            //
// the tet core uses a tetrahedron-based data structure.  It stores tetrahedra and  //
// vertices.  This data structure is pointer-based. Each tetrahedron contains //
// pointers to its vertices and adjacent tetrahedra. Each vertex holds its x-,//
// y-, z-coordinates, and a pointer to one of the tetrahedra having it.  Both //
// tetrahedra and vertices may contain user data.                             //
//                                                                            //
// Let T be a tetrahedralization. Each triangular face of T belongs to either //
// two or one tetrahedron.  In the latter case, it is an exterior boundary    //
// face of T. the tet core attaches tetrahedra (one-to-one) to such faces. All such //
// tetrahedra contain an "infinite vertex" (which has no geometric coordinates//
// ).  One can imagine such a vertex lies in 4D space and is visible by all   //
// exterior boundary faces simultaneously.  This extended set of tetrahedra   //
// (including the infinite vertex) becomes a tetrahedralization of a 3-sphere //
// that has no boundary in 3d. It has the nice property that every triangular //
// face is shared by exactly two tetrahedra.                                  //
//                                                                            //
// The current version of the tet core stores explicitly the subfaces and segments  //
// (which are in surface mesh S and the linear mesh L), respectively.  Extra  //
// pointers are allocated in tetrahedra and subfaces to point each other.     //
//                                                                            //
//============================================================================//

  // The tetrahedron data structure.  It includes the following fields:
  //   - a list of four adjoining tetrahedra;
  //   - a list of four vertices;
  //   - a pointer to a list of four subfaces (optional, for -p switch);
  //   - a pointer to a list of six segments  (optional, for -p switch);
  //   - a list of user-defined floating-point attributes (optional);
  //   - a volume constraint (optional, for -a switch);
  //   - an integer of element marker (and flags);
  // The structure of a tetrahedron is an array of pointers.  Its actual size
  //   (the length of the array) is determined at runtime.

  typedef double **tetrahedron;

  // The subface data structure.  It includes the following fields:
  //   - a list of three adjoining subfaces;
  //   - a list of three vertices;
  //   - a list of three adjoining segments;
  //   - two adjoining tetrahedra;
  //   - an area constraint (optional, for -q switch);
  //   - an integer for boundary marker;
  //   - an integer for type, flags, etc.

  typedef double **shellface;

  // The point data structure.  It includes the following fields:
  //   - x, y and z coordinates;
  //   - a list of user-defined point attributes (optional);
  //   - u, v coordinates (optional, for -s switch);
  //   - a metric tensor (optional, for -q or -m switch);
  //   - a pointer to an adjacent tetrahedron;
  //   - a pointer to a parent (or a duplicate) point;
  //   - a pointer to an adjacent subface or segment (optional, -p switch);
  //   - a pointer to a tet in background mesh (optional, for -m switch);
  //   - an integer for boundary marker (point index);
  //   - an integer for point type (and flags).
  //   - an integer for geometry tag (optional, for -s switch).
  // The structure of a point is an array of REALs.  Its acutal size is 
  //   determined at the runtime.

  typedef double *point;

//============================================================================//
//                                                                            //
// Handles                                                                    //
//                                                                            //
// Navigation and manipulation in a tetrahedralization are accomplished by    //
// operating on structures referred as ``handles". A handle is a pair (t,v),  //
// where t is a pointer to a tetrahedron, and v is a 4-bit integer, in the    //
// range from 0 to 11. v is called the ``version'' of a tetrahedron, it rep-  //
// resents a directed edge of a specific face of the tetrahedron.             //
//                                                                            //
// There are 12 even permutations of the four vertices, each of them corres-  //
// ponds to a directed edge (a version) of the tetrahedron.  The 12 versions  //
// can be grouped into 4 distinct ``edge rings'' in 4 ``oriented faces'' of   //
// this tetrahedron.  One can encode each version (a directed edge) into a    //
// 4-bit integer such that the two upper bits encode the index (from 0 to 2)  //
// of this edge in the edge ring, and the two lower bits encode the index (   //
// from 0 to 3) of the oriented face which contains this edge.                //
//                                                                            //
// The four vertices of a tetrahedron are indexed from 0 to 3 (according to   //
// their storage in the data structure).  Give each face the same index as    //
// the node opposite it in the tetrahedron.  Denote the edge connecting face  //
// i to face j as i/j. We number the twelve versions as follows:              //
//                                                                            //
//           |   edge 0     edge 1     edge 2                                 //
//   --------|--------------------------------                                //
//    face 0 |   0 (0/1)    4 (0/3)    8 (0/2)                                //
//    face 1 |   1 (1/2)    5 (1/3)    9 (1/0)                                //
//    face 2 |   2 (2/3)    6 (2/1)   10 (2/0)                                //
//    face 3 |   3 (3/0)    7 (3/1)   11 (3/2)                                //
//                                                                            //
// Similarly, navigation and manipulation in a (boundary) triangulation are   //
// done by using handles of triangles. Each handle is a pair (s, v), where s  //
// is a pointer to a triangle, and v is a version in the range from 0 to 5.   //
// Each version corresponds to a directed edge of this triangle.              //
//                                                                            //
// Number the three vertices of a triangle from 0 to 2 (according to their    //
// storage in the data structure). Give each edge the same index as the node  //
// opposite it in the triangle. The six versions of a triangle are:           //
//                                                                            //
//                 | edge 0   edge 1   edge 2                                 //
//  ---------------|--------------------------                                //
//   ccw orieation |   0        2        4                                    //
//    cw orieation |   1        3        5                                    //
//                                                                            //
// In the following, a 'triface' is a handle of tetrahedron, and a 'face' is  //
// a handle of a triangle.                                                    //
//                                                                            //
//============================================================================//

  class triface {
  public:
    tetrahedron *tet;
    int ver; // Range from 0 to 11.
    triface() : tet(0), ver(0) {}
    triface& operator=(const triface& t) {
      tet = t.tet; ver = t.ver;
      return *this;
    }
  };

  class face {
  public:
    shellface *sh;
    int shver; // Range from 0 to 5.
    face() : sh(0), shver(0) {}
    face& operator=(const face& s) {
      sh = s.sh; shver = s.shver;
      return *this;
    }
  };

//============================================================================//
//                                                                            //
// Arraypool                                                                  //
//                                                                            //
// A dynamic linear array. (It is written by J. Shewchuk)                     //
//                                                                            //
// Each arraypool contains an array of pointers to a number of blocks.  Each  //
// block contains the same fixed number of objects.  Each index of the array  //
// addresses a particular object in the pool. The most significant bits add-  //
// ress the index of the block containing the object. The less significant    //
// bits address this object within the block.                                 //
//                                                                            //
// 'objectbytes' is the size of one object in blocks; 'log2objectsperblock'   //
// is the base-2 logarithm of 'objectsperblock'; 'objects' counts the number  //
// of allocated objects; 'totalmemory' is the total memory in bytes.          //
//                                                                            //
//============================================================================//

  class arraypool {

  public:

    int objectbytes;
    int objectsperblock;
    int log2objectsperblock;
    int objectsperblockmark;
    int toparraylen;
    char **toparray;
    long objects;
    unsigned long totalmemory;

    void restart();
    void poolinit(int sizeofobject, int log2objperblk);
    char* getblock(int objectindex);
    void* lookup(int objectindex);
    int newindex(void **newptr);

    arraypool(int sizeofobject, int log2objperblk);
    ~arraypool();
  };

// fastlookup() -- A fast, unsafe operation. Return the pointer to the object
//   with a given index.  Note: The object's block must have been allocated,
//   i.e., by the function newindex().

#define fastlookup(pool, index) \
  (void *) ((pool)->toparray[(index) >> (pool)->log2objectsperblock] + \
            ((index) & (pool)->objectsperblockmark) * (pool)->objectbytes)

//============================================================================//
//                                                                            //
// Memorypool                                                                 //
//                                                                            //
// A structure for memory allocation. (It is written by J. Shewchuk)          //
//                                                                            //
// firstblock is the first block of items. nowblock is the block from which   //
//   items are currently being allocated. nextitem points to the next slab    //
//   of free memory for an item. deaditemstack is the head of a linked list   //
//   (stack) of deallocated items that can be recycled.  unallocateditems is  //
//   the number of items that remain to be allocated from nowblock.           //
//                                                                            //
// Traversal is the process of walking through the entire list of items, and  //
//   is separate from allocation.  Note that a traversal will visit items on  //
//   the "deaditemstack" stack as well as live items.  pathblock points to    //
//   the block currently being traversed.  pathitem points to the next item   //
//   to be traversed.  pathitemsleft is the number of items that remain to    //
//   be traversed in pathblock.                                               //
//                                                                            //
//============================================================================//

  class memorypool {

  public:

    void **firstblock, **nowblock;
    void *nextitem;
    void *deaditemstack;
    void **pathblock;
    void *pathitem;
    int  alignbytes;
    int  itembytes, itemwords;
    int  itemsperblock;
    long items, maxitems;
    int  unallocateditems;
    int  pathitemsleft;

    memorypool();
    memorypool(int, int, int, int);
    ~memorypool();
    
    void poolinit(int, int, int, int);
    void restart();
    void *alloc();
    void dealloc(void*);
    void traversalinit();
    void *traverse();
  };  

//============================================================================//
//                                                                            //
// badface                                                                    //
//                                                                            //
// Despite of its name, a 'badface' can be used to represent one of the       //
// following objects:                                                         //
//   - a face of a tetrahedron which is (possibly) non-Delaunay;              //
//   - an encroached subsegment or subface;                                   //
//   - a bad-quality tetrahedron, i.e, has too large radius-edge ratio;       //
//   - a sliver, i.e., has good radius-edge ratio but nearly zero volume;     //
//   - a recently flipped face (saved for undoing the flip later).            //
//                                                                            //
//============================================================================//

  class badface {
  public:
    triface tt; 
    face ss;
    double key, cent[6];  // circumcenter or cos(dihedral angles) at 6 edges.
    point forg, fdest, fapex, foppo, noppo;
    badface *nextitem; 
    badface() : key(0), forg(0), fdest(0), fapex(0), foppo(0), noppo(0),
      nextitem(0) {}
    void init() {
      key = 0.;
      for (int k = 0; k < 6; k++) cent[k] = 0.;
      tt.tet = NULL; tt.ver = 0;
      ss.sh = NULL; ss.shver = 0;
      forg = fdest = fapex = foppo = noppo = NULL;
      nextitem = NULL;
    }
  };

//============================================================================//
//                                                                            //
// insertvertexflags                                                          //
//                                                                            //
// A collection of flags that pass to the routine insertvertex().             //
//                                                                            //
//============================================================================//

  class insertvertexflags {

  public:

    int iloc;  // input/output.
    int bowywat, lawson;
    int splitbdflag, validflag, respectbdflag;
    int rejflag, chkencflag, cdtflag;
    int assignmeshsize;
    int sloc, sbowywat;
    
    // Used by Delaunay refinement.
    int collect_inial_cavity_flag;
    int ignore_near_vertex;
    int check_insert_radius;
    int refineflag; // 0, 1, 2, 3
    triface refinetet;
    face refinesh;
    int smlenflag; // for useinsertradius.
    double smlen; // for useinsertradius.
    point parentpt;

    void init() {
      iloc = bowywat = lawson = 0;
      splitbdflag = validflag = respectbdflag = 0;
      rejflag = chkencflag = cdtflag = 0;
      assignmeshsize = 0;
      sloc = sbowywat = 0;

      collect_inial_cavity_flag = 0;
      ignore_near_vertex = 0;
      check_insert_radius = 0;
      refineflag = 0;
      refinetet.tet = NULL;
      refinesh.sh = NULL;
      smlenflag = 0;
      smlen = 0.0;
      parentpt = NULL;
    }

    insertvertexflags() {
      init();
    }
  };

//============================================================================//
//                                                                            //
// flipconstraints                                                            //
//                                                                            //
// A structure of a collection of data (options and parameters) which pass    //
// to the edge flip function flipnm().                                        //
//                                                                            //
//============================================================================//

  class flipconstraints {

  public:

    // Elementary flip flags.
    int enqflag; // (= flipflag)
    int chkencflag;

    // Control flags
    int unflip;  // Undo the performed flips.
    int collectnewtets; // Collect the new tets created by flips.
    int collectencsegflag;

    // Optimization flags.
    int noflip_in_surface; // do not flip edges (not segment) in surface.
    int remove_ndelaunay_edge; // Remove a non-Delaunay edge.
    double bak_tetprism_vol; // The value to be minimized.
    double tetprism_vol_sum;
    int remove_large_angle; // Remove a large dihedral angle at edge.
    double cosdihed_in; // The input cosine of the dihedral angle (> 0).
    double cosdihed_out; // The improved cosine of the dihedral angle.
    double max_asp_out; // Max asp ratio after the improvement of dihedral angle.

    // Boundary recovery flags.
    int checkflipeligibility;
    point seg[2];  // A constraining edge to be recovered.
    point fac[3];  // A constraining face to be recovered.
    point remvert; // A vertex to be removed.


    flipconstraints() {
      enqflag = 0; 
      chkencflag = 0;

      unflip = 0;
      collectnewtets = 0;
      collectencsegflag = 0;

      noflip_in_surface = 0;
      remove_ndelaunay_edge = 0;
      bak_tetprism_vol = 0.0;
      tetprism_vol_sum = 0.0;
      remove_large_angle = 0;
      cosdihed_in = 0.0;
      cosdihed_out = 0.0;
      max_asp_out = 0.0;

      checkflipeligibility = 0;
      seg[0] = NULL;
      fac[0] = NULL;
      remvert = NULL;
    }
  };

//============================================================================//
//                                                                            //
// optparameters                                                              //
//                                                                            //
// Optimization options and parameters.                                       //
//                                                                            //
//============================================================================//

  class optparameters {

  public:

    // The one of goals of optimization.
    int max_min_volume;      // Maximize the minimum volume.
	int min_max_aspectratio; // Minimize the maximum aspect ratio. 
    int min_max_dihedangle;  // Minimize the maximum dihedral angle.

    // The initial and improved value.
    double initval, imprval;

    int numofsearchdirs;
    double searchstep;
    int maxiter;  // Maximum smoothing iterations (disabled by -1).
    int smthiter; // Performed iterations.


    optparameters() {
      max_min_volume = 0;
      min_max_aspectratio = 0;
      min_max_dihedangle = 0;

      initval = imprval = 0.0;

      numofsearchdirs = 10;
      searchstep = 0.01;
      maxiter = -1;   // Unlimited smoothing iterations.
      smthiter = 0;

    }
  };


//============================================================================//
//                                                                            //
// Labels (enumeration declarations) used by the tet core.                          //
//                                                                            //
//============================================================================//

  // Labels that signify the type of a vertex. 
  enum verttype {UNUSEDVERTEX, DUPLICATEDVERTEX, RIDGEVERTEX, /*ACUTEVERTEX,*/
                 FACETVERTEX, VOLVERTEX, FREESEGVERTEX, FREEFACETVERTEX, 
                 FREEVOLVERTEX, NREGULARVERTEX, DEADVERTEX};
 
  // Labels that signify the result of triangle-triangle intersection test.
  enum interresult {DISJOINT, INTERSECT, SHAREVERT, SHAREEDGE, SHAREFACE,
                    TOUCHEDGE, TOUCHFACE, ACROSSVERT, ACROSSEDGE, ACROSSFACE,
                    SELF_INTERSECT};

  // Labels that signify the result of point location.
  enum locateresult {UNKNOWN, OUTSIDE, INTETRAHEDRON, ONFACE, ONEDGE, ONVERTEX,
                     ENCVERTEX, ENCSEGMENT, ENCSUBFACE, NEARVERTEX, NONREGULAR,
                     INSTAR, BADELEMENT, NULLCAVITY, SHARPCORNER, FENSEDIN,
                     NONCOPLANAR, SELF_ENCROACH};

//============================================================================//
//                                                                            //
// Variables of the tet core                                                        //
//                                                                            //
//============================================================================//

  // Pointer to the input data (a set of nodes, a PLC, or a mesh).
  TetMeshData *in, *addin;

  // Pointer to the switches and parameters.
  TetMeshBehavior *b;

  // Pointer to a background mesh (contains size specification map).
  TetMeshCore *bgm;

  // Memorypools to store mesh elements (points, tetrahedra, subfaces, and
  //   segments) and extra pointers between tetrahedra, subfaces, and segments.
  memorypool *tetrahedrons, *subfaces, *subsegs, *points;
  memorypool *tet2subpool, *tet2segpool;

  // Memorypools to store bad-quality (or encroached) elements.
  memorypool *badtetrahedrons, *badsubfacs, *badsubsegs;
  memorypool *split_subfaces_pool, *split_segments_pool;
  arraypool  *unsplit_badtets, *unsplit_subfaces, *unsplit_segments;
  arraypool  *check_tets_list;

  badface *stack_enc_segments, *stack_enc_subfaces;
  
  // Bad quality subfaces are ordered by priority queues.
  badface *queuefront[64];
  badface *queuetail[64];
  int nextnonemptyq[64];
  int firstnonemptyq, recentq;

  // Bad quality tetrahedra are ordered by priority queues.
  memorypool *badqual_tets_pool;
  badface *bt_queuefront[64];
  badface *bt_queuetail[64];
  int bt_nextnonemptyq[64];
  int bt_firstnonemptyq, bt_recentq;
  
  // A memorypool to store faces to be flipped.
  memorypool *flippool;
  arraypool *later_unflip_queue, *unflipqueue;
  badface *flipstack, *unflip_queue_front, *unflip_queue_tail;

  // Arrays used for point insertion (the Bowyer-Watson algorithm).
  arraypool *cavetetlist, *cavebdrylist, *caveoldtetlist;
  arraypool *cave_oldtet_list; // only tetrahedron's
  arraypool *cavetetshlist, *cavetetseglist, *cavetetvertlist;
  arraypool *caveencshlist, *caveencseglist;
  arraypool *caveshlist, *caveshbdlist, *cavesegshlist;
  triface _bw_faces[4096]; // _bw_faces[64][64];

  // Stacks used for CDT construction and boundary recovery.
  arraypool *subsegstack, *subfacstack, *subvertstack;
  arraypool *skipped_segment_list, *skipped_facet_list;

  // Arrays of encroached segments and subfaces (for mesh refinement).
  arraypool *encseglist, *encshlist;

  // The map between facets to their vertices (for mesh refinement).
  int    number_of_facets;
  int   *idx2facetlist;
  point *facetverticeslist;
  int    *idx_segment_facet_list; // segment-to-facet map.
  int    *segment_facet_list;
  int    *idx_ridge_vertex_facet_list; // vertex-to-facet map.
  int    *ridge_vertex_facet_list;

  // The map between segments to their endpoints (for mesh refinement).
  int    segmentendpointslist_length;
  point  *segmentendpointslist;
  double *segment_info_list;
  int    *idx_segment_ridge_vertex_list; // are two ridge vertices form a segment?
  point  *segment_ridge_vertex_list;

  // The infinite vertex.
  point dummypoint;
  // The recently visited tetrahedron, subface.
  triface recenttet;
  face recentsh;

  // PI is the ratio of a circle's circumference to its diameter.
  static double PI;

  // The list of subdomains. (-A option).
  int subdomains;                                    // Number of subdomains.
  int *subdomain_markers;

  // Various variables.
  int numpointattrib;                          // Number of point attributes.
  int numelemattrib;                     // Number of tetrahedron attributes.
  int sizeoftensor;                     // Number of REALs per metric tensor.
  int pointmtrindex;           // Index to find the metric tensor of a point.
  int pointparamindex;       // Index to find the u,v coordinates of a point.
  int point2simindex;         // Index to find a simplex adjacent to a point.
  int pointmarkindex;            // Index to find boundary marker of a point.
  int pointinsradiusindex;  // Index to find the insertion radius of a point.
  int elemattribindex;          // Index to find attributes of a tetrahedron.
  int polarindex;                // Index to find the polar plane parameters.
  int volumeboundindex;       // Index to find volume bound of a tetrahedron.
  int elemmarkerindex;              // Index to find marker of a tetrahedron.
  int shmarkindex;             // Index to find boundary marker of a subface.
  int areaboundindex;               // Index to find area bound of a subface.
  int checksubsegflag;   // Are there segments in the tetrahedralization yet?
  int checksubfaceflag;  // Are there subfaces in the tetrahedralization yet?
  int boundary_recovery_flag;
  int checkconstraints;  // Are there variant (node, seg, facet) constraints?
  int nonconvex;                               // Is current mesh non-convex?
  int autofliplinklevel;        // The increase of link levels, default is 1.
  int useinsertradius;       // Save the insertion radius for Steiner points.
  long samples;               // Number of random samples for point location.
  unsigned long randomseed;                    // Current random number seed.
  double cosmaxdihed, cosmindihed;    // The cosine values of max/min dihedral.
  double cossmtdihed;     // The cosine value of a bad dihedral to be smoothed.
  double cosslidihed;      // The cosine value of the max dihedral of a sliver.
  double cos_large_dihed;   // The cosine value of large dihedral (135 degree).
  double opt_max_sliver_asp_ratio;              // = 10 x b->opt_max_asp_ratio.
  double minfaceang, minfacetdihed;     // The minimum input (dihedral) angles.
  double cos_facet_separate_ang_tol;
  double cos_collinear_ang_tol;
  double tetprism_vol_sum;   // The total volume of tetrahedral-prisms (in 4D).
  double longest;                          // The longest possible edge length.
  double minedgelength;                               // = longest * b->epsion.
  double xmax, xmin, ymax, ymin, zmax, zmin;         // Bounding box of points.

  // Options for mesh refinement.
  double big_radius_edge_ratio;           // calculated by qualitystatistics().
  double smallest_insradius;             // Save the smallest insertion radius.
  long elem_limit;
  long insert_point_count;                 // number of attempted insertions.
  long report_refine_progress;                      // the next report event.
  long last_point_count;         // number of points after last report event.
  long last_insertion_count; // number of insertions after last report event.

  // Counters.  
  long insegments;                               // Number of input segments.  
  long hullsize;                        // Number of exterior boundary faces.
  long meshedges;                                    // Number of mesh edges.
  long meshhulledges;                       // Number of boundary mesh edges.
  long steinerleft;                 // Number of Steiner points not yet used.
  long dupverts;                            // Are there duplicated vertices?
  long unuverts;                                // Are there unused vertices?
  long duplicated_facets_count;              // Are there duplicated facets.?
  long nonregularcount;                    // Are there non-regular vertices?
  long st_segref_count, st_facref_count, st_volref_count;  // Steiner points.
  long fillregioncount, cavitycount, cavityexpcount;
  long flip14count, flip26count, flipn2ncount;
  long flip23count, flip32count, flip44count, flip41count;
  long flip31count, flip22count;
  long opt_flips_count, opt_collapse_count, opt_smooth_count;
  long recover_delaunay_count;
  unsigned long totalworkmemory;      // Total memory used by working arrays.


//============================================================================//
//                                                                            //
// Mesh manipulation primitives                                               //
//                                                                            //
//============================================================================//

  // Fast lookup tables for mesh manipulation primitives.
  static int bondtbl[12][12], fsymtbl[12][12];
  static int esymtbl[12], enexttbl[12], eprevtbl[12];
  static int enextesymtbl[12], eprevesymtbl[12]; 
  static int eorgoppotbl[12], edestoppotbl[12];
  static int facepivot1[12], facepivot2[12][12];
  static int orgpivot[12], destpivot[12], apexpivot[12], oppopivot[12];
  static int tsbondtbl[12][6], stbondtbl[12][6];
  static int tspivottbl[12][6], stpivottbl[12][6];
  static int ver2edge[12], edge2ver[6], epivot[12];
  static int sorgpivot [6], sdestpivot[6], sapexpivot[6];
  static int snextpivot[6];

  void inittables();

  // Primitives for tetrahedra.
  inline tetrahedron encode(triface& t);
  inline tetrahedron encode2(tetrahedron* ptr, int ver);
  inline void decode(tetrahedron ptr, triface& t);
  inline tetrahedron* decode_tet_only(tetrahedron ptr);
  inline int  decode_ver_only(tetrahedron ptr);
  inline void bond(triface& t1, triface& t2);
  inline void dissolve(triface& t);
  inline void esym(triface& t1, triface& t2);
  inline void esymself(triface& t);
  inline void enext(triface& t1, triface& t2);
  inline void enextself(triface& t);
  inline void eprev(triface& t1, triface& t2);
  inline void eprevself(triface& t);
  inline void enextesym(triface& t1, triface& t2);
  inline void enextesymself(triface& t);
  inline void eprevesym(triface& t1, triface& t2);
  inline void eprevesymself(triface& t);
  inline void eorgoppo(triface& t1, triface& t2);
  inline void eorgoppoself(triface& t);
  inline void edestoppo(triface& t1, triface& t2);
  inline void edestoppoself(triface& t);
  inline void fsym(triface& t1, triface& t2);
  inline void fsymself(triface& t);
  inline void fnext(triface& t1, triface& t2);
  inline void fnextself(triface& t);
  inline point org (triface& t);
  inline point dest(triface& t);
  inline point apex(triface& t);
  inline point oppo(triface& t);
  inline void setorg (triface& t, point p);
  inline void setdest(triface& t, point p);
  inline void setapex(triface& t, point p);
  inline void setoppo(triface& t, point p);
  inline double elemattribute(tetrahedron* ptr, int attnum);
  inline void setelemattribute(tetrahedron* ptr, int attnum, double value);
  inline double* get_polar(tetrahedron* ptr);
  inline double get_volume(tetrahedron* ptr);
  inline double volumebound(tetrahedron* ptr);
  inline void setvolumebound(tetrahedron* ptr, double value);
  inline int  elemindex(tetrahedron* ptr);
  inline void setelemindex(tetrahedron* ptr, int value);
  inline int  elemmarker(tetrahedron* ptr);
  inline void setelemmarker(tetrahedron* ptr, int value);
  inline void infect(triface& t);
  inline void uninfect(triface& t);
  inline bool infected(triface& t);
  inline void marktest(triface& t);
  inline void unmarktest(triface& t);
  inline bool marktested(triface& t);
  inline void markface(triface& t);
  inline void unmarkface(triface& t);
  inline bool facemarked(triface& t);
  inline void markedge(triface& t);
  inline void unmarkedge(triface& t);
  inline bool edgemarked(triface& t);
  inline void marktest2(triface& t);
  inline void unmarktest2(triface& t);
  inline bool marktest2ed(triface& t);
  inline int  elemcounter(triface& t);
  inline void setelemcounter(triface& t, int value);
  inline void increaseelemcounter(triface& t);
  inline void decreaseelemcounter(triface& t);
  inline bool ishulltet(triface& t);
  inline bool isdeadtet(triface& t);
 
  // Primitives for subfaces and subsegments.
  inline void sdecode(shellface sptr, face& s);
  inline shellface sencode(face& s);
  inline shellface sencode2(shellface *sh, int shver);
  inline void spivot(face& s1, face& s2);
  inline void spivotself(face& s);
  inline void sbond(face& s1, face& s2);
  inline void sbond1(face& s1, face& s2);
  inline void sdissolve(face& s);
  inline point sorg(face& s);
  inline point sdest(face& s);
  inline point sapex(face& s);
  inline void setsorg(face& s, point pointptr);
  inline void setsdest(face& s, point pointptr);
  inline void setsapex(face& s, point pointptr);
  inline void sesym(face& s1, face& s2);
  inline void sesymself(face& s);
  inline void senext(face& s1, face& s2);
  inline void senextself(face& s);
  inline void senext2(face& s1, face& s2);
  inline void senext2self(face& s);
  inline double areabound(face& s);
  inline void setareabound(face& s, double value);
  inline int shellmark(face& s);
  inline void setshellmark(face& s, int value);
  inline void sinfect(face& s);
  inline void suninfect(face& s);
  inline bool sinfected(face& s);
  inline void smarktest(face& s);
  inline void sunmarktest(face& s);
  inline bool smarktested(face& s);
  inline void smarktest2(face& s);
  inline void sunmarktest2(face& s);
  inline bool smarktest2ed(face& s);
  inline void smarktest3(face& s);
  inline void sunmarktest3(face& s);
  inline bool smarktest3ed(face& s);
  inline void setfacetindex(face& f, int value);
  inline int  getfacetindex(face& f);
  inline bool isdeadsh(face& s);

  // Primitives for interacting tetrahedra and subfaces.
  inline void tsbond(triface& t, face& s);
  inline void tsdissolve(triface& t);
  inline void stdissolve(face& s);
  inline void tspivot(triface& t, face& s);
  inline void stpivot(face& s, triface& t);

  // Primitives for interacting tetrahedra and segments.
  inline void tssbond1(triface& t, face& seg);
  inline void sstbond1(face& s, triface& t);
  inline void tssdissolve1(triface& t);
  inline void sstdissolve1(face& s);
  inline void tsspivot1(triface& t, face& s);
  inline void sstpivot1(face& s, triface& t);

  // Primitives for interacting subfaces and segments.
  inline void ssbond(face& s, face& edge);
  inline void ssbond1(face& s, face& edge);
  inline void ssdissolve(face& s);
  inline void sspivot(face& s, face& edge);

  // Primitives for points.
  inline int  pointmark(point pt);
  inline void setpointmark(point pt, int value);
  inline enum verttype pointtype(point pt);
  inline void setpointtype(point pt, enum verttype value);
  inline int  pointgeomtag(point pt);
  inline void setpointgeomtag(point pt, int value);
  inline double pointgeomuv(point pt, int i);
  inline void setpointgeomuv(point pt, int i, double value);
  inline void pinfect(point pt);
  inline void puninfect(point pt);
  inline bool pinfected(point pt);
  inline void pmarktest(point pt);
  inline void punmarktest(point pt);
  inline bool pmarktested(point pt);
  inline void pmarktest2(point pt);
  inline void punmarktest2(point pt);
  inline bool pmarktest2ed(point pt);
  inline void pmarktest3(point pt);
  inline void punmarktest3(point pt);
  inline bool pmarktest3ed(point pt);
  inline tetrahedron point2tet(point pt);
  inline void setpoint2tet(point pt, tetrahedron value);
  inline shellface point2sh(point pt);
  inline void setpoint2sh(point pt, shellface value);
  inline point point2ppt(point pt);
  inline void setpoint2ppt(point pt, point value);
  inline tetrahedron point2bgmtet(point pt);
  inline void setpoint2bgmtet(point pt, tetrahedron value);
  inline void setpointinsradius(point pt, double value);
  inline double getpointinsradius(point pt);
  inline bool issteinerpoint(point pt);

  // Advanced primitives.
  inline void point2tetorg(point pt, triface& t);
  inline void point2shorg(point pa, face& s);
  inline point farsorg(face& seg);
  inline point farsdest(face& seg);

//============================================================================//
//                                                                            //
//  Memory managment                                                          //
//                                                                            //
//============================================================================//

  void tetrahedrondealloc(tetrahedron*);
  tetrahedron *tetrahedrontraverse();
  tetrahedron *alltetrahedrontraverse();
  void shellfacedealloc(memorypool*, shellface*);
  shellface *shellfacetraverse(memorypool*);
  void pointdealloc(point);
  point pointtraverse();

  void makeindex2pointmap(point*&);
  void makepoint2submap(memorypool*, int*&, face*&);
  void maketetrahedron(triface*);
  void maketetrahedron2(triface*, point, point, point, point);
  void makeshellface(memorypool*, face*);
  void makepoint(point*, enum verttype);

  void initializepools();

//============================================================================//
//                                                                            //
// Advanced geometric predicates and calculations                             //
//                                                                            //
// the routine insphere_s() implements a simplified symbolic perturbation     //
// scheme from Edelsbrunner, et al [*].  Hence the point-in-sphere test never //
// returns a zero. The idea is to perturb the weights of vertices in 4D.      //
//                                                                            //
// The routine tri_edge_test() determines whether or not a triangle and an    //
// edge intersect in 3D. If they do cross, their intersection type is also    //
// reported. This test is a combination of n 3D orientation tests (3 < n < 9).//
// It uses the robust orient3d() test to make the branch decisions.           //
//                                                                            //
// There are several routines to calculate geometrical quantities, e.g.,      //
// circumcenters, angles, dihedral angles, face normals, face areas, etc.     //
// They are implemented using floating-point arithmetics.                     //
//                                                                            //
//============================================================================//

  // Symbolic perturbations (robust)
  double insphere_s(double*, double*, double*, double*, double*);
  double orient4d_s(double*, double*, double*, double*, double*, 
                  double, double, double, double, double);

  // An embedded 2-dimensional geometric predicate (non-robust)
  double incircle3d(point pa, point pb, point pc, point pd);

  // Triangle-edge intersection test (robust)
  int tri_edge_2d(point, point, point, point, point, point, int, int*, int*);
  int tri_edge_tail(point,point,point,point,point,point,double,double,int,int*,int*);
  int tri_edge_test(point, point, point, point, point, point, int, int*, int*);

    // Triangle-triangle intersection test (robust)
  int tri_edge_inter_tail(point, point, point, point, point, double, double);
  int tri_tri_inter(point, point, point, point, point, point);

  // Linear algebra functions
  inline double dot(double* v1, double* v2);
  inline void cross(double* v1, double* v2, double* n);
  bool lu_decmp(double lu[4][4], int n, int* ps, double* d, int N);
  void lu_solve(double lu[4][4], int n, int* ps, double* b, int N);

  // Geometric calculations (non-robust)
  double orient3dfast(double *pa, double *pb, double *pc, double *pd);
  inline double norm2(double x, double y, double z);
  inline double distance(double* p1, double* p2);
  inline double distance2(double* p1, double* p2);
  void facenormal(point pa, point pb, point pc, double *n, int pivot, double *lav);
  double facedihedral(double* pa, double* pb, double* pc1, double* pc2);
  double triarea(double* pa, double* pb, double* pc);
  double interiorangle(double* o, double* p1, double* p2, double* n);
  double cos_interiorangle(double* o, double* p1, double* p2);
  void projpt2edge(double* p, double* e1, double* e2, double* prj);
  void projpt2face(double* p, double* f1, double* f2, double* f3, double* prj);
  bool circumsphere(double*, double*, double*, double*, double* cent, double* radius);
  bool orthosphere(double*,double*,double*,double*,double,double,double,double,double*,double*);
  void planelineint(double*, double*, double*, double*, double*, double*, double*);
  int  linelineint(double*, double*, double*, double*, double*, double*, double*, double*);
  double tetprismvol(double* pa, double* pb, double* pc, double* pd);
  bool calculateabovepoint(arraypool*, point*, point*, point*);
  void calculateabovepoint4(point, point, point, point);

//============================================================================//
//                                                                            //
// Local mesh transformations                                                 //
//                                                                            //
// A local transformation replaces a set of tetrahedra with another set that  //
// partitions the same space and boundaries.                                  //
//                                                                            //
// In 3D, the most straightforward local transformations are the elementary   //
// flips performed within the convex hull of five vertices: 2-to-3, 3-to-2,   //
// 1-to-4, and 4-to-1 flips. The numbers indicate the number of tetrahedra    //
// before and after each flip.  The 1-to-4 and 4-to-1 flip involve inserting  //
// or deleting a vertex, respectively.                                        //
//                                                                            //
// There are complex local transformations that are a combination of element- //
// ary flips. For example, a 4-to-4 flip, which replaces two coplanar edges,  //
// combines a 2-to-3 flip and a 3-to-2 flip. Note that the first 2-to-3 flip  //
// will temporarily create a degenerate tetrahedron removed immediately by    //
// the followed 3-to-2 flip.  More generally, an n-to-m flip, where n > 3,    //
// m = (n - 2) * 2, which removes an edge, can be done by first performing a  //
// sequence of (n - 3) 2-to-3 flips followed by a 3-to-2 flip.                //
//                                                                            //
// The routines flip23(), flip32(), and flip41() perform the three elementray //
// flips. The flip14() is available inside the routine insertpoint().         //
//                                                                            //
// The routines flipnm() and flipnm_post() implement a generalized edge flip  //
// algorithm that uses elementary flips.                                      //
//                                                                            //
// The routine insertpoint() implements the Bowyer-Watson's cavity algorithm  //
// to insert a vertex.  It works for arbitrary tetrahedralization,  either    //
// Delaunay, or constrained Delaunay, or non-Delaunay.                        //
//                                                                            //
//============================================================================//

  void flippush(badface*&, triface*);

  // The elementary flips.
  void flip23(triface*, int, flipconstraints* fc);
  void flip32(triface*, int, flipconstraints* fc);
  void flip41(triface*, int, flipconstraints* fc);

  // A generalized edge flip.
  int flipnm(triface*, int n, int level, int, flipconstraints* fc);
  int flipnm_post(triface*, int n, int nn, int, flipconstraints* fc);

  // Point insertion.
  int  insertpoint(point, triface*, face*, face*, insertvertexflags*);
  void insertpoint_abort(face*, insertvertexflags*);

//============================================================================//
//                                                                            //
// Delaunay tetrahedralization                                                //
//                                                                            //
// The routine incrementaldelaunay() implemented two incremental algorithms   //
// for constructing Delaunay tetrahedralizations (DTs):  the Bowyer-Watson    //
// (B-W) algorithm and the incremental flip algorithm of Edelsbrunner and     //
// Shah, "Incremental topological flipping works for regular triangulation,"  //
// Algorithmica, 15:233-241, 1996.                                            //
//                                                                            //
// The routine incrementalflip() implements the flip algorithm of [Edelsbrun- //
// ner and Shah, 1996].  It flips a queue of locally non-Delaunay faces (in   //
// arbitrary order).  The success is guaranteed when the Delaunay tetrahedra- //
// lization is constructed incrementally by adding one vertex at a time.      //
//                                                                            //
// The routine locate() finds a tetrahedron contains a new point in current   //
// DT. It uses a simple stochastic walk algorithm: starting from an arbitrary //
// tetrahedron in DT, it finds the destination by visit one tetrahedron at a  //
// time, randomly chooses a tetrahedron if there are more than one choices.   //
// This algorithm terminates due to Edelsbrunner's acyclic theorem.           //
//                                                                            //
// Choose a good starting tetrahedron is crucial to the speed of the walk.    //
// the tet core initially uses the "jump-and-walk" algorithm of Muecke, E.P., Saias,//
// I., and Zhu, B. "Fast Randomized Point Location Without Preprocessing." In //
// Proceedings of the 12th ACM Symposium on Computational Geometry, 274-283,  //
// 1996.  It first randomly samples several tetrahedra in the DT and then     //
// choosing the closet one to start walking.                                  //
//                                                                            //
// The above algorithm slows download dramatically as the number of points    //
// grows -- reported in Amenta, N., Choi, S. and Rote, G., "Incremental       //
// construction con {BRIO}," In Proceedings of 19th ACM Symposium on Computa- //
// tional Geometry, 211-219, 2003. On the other hand, Liu and Snoeyink showed //
// that the point location could be made in constant time if the points are   //
// pre-sorted so that the nearby points in space have nearby indices, then    //
// adding the points in this order. They sorted the points along the 3D       //
// Hilbert curve.                                                             //
//                                                                            //
// The routine hilbert_sort3() sorts a set of 3D points along the 3D Hilbert  //
// curve. It recursively splits a point set according to the Hilbert indices  //
// mapped to the subboxes of the bounding box of the point set. The Hilbert   //
// indices is calculated by Butz's algorithm in 1971. An excellent exposition //
// of this algorithm can be found in the paper of Hamilton, C., "Compact      //
// Hilbert Indices", Technical Report CS-2006-07, Computer Science, Dalhousie //
// University, 2006 (the Section 2). My implementation also referenced Steven //
// Witham's performance of "Hilbert walk" (hopefully, it is still available   //
// at http://www.tiac.net/~sw/2008/10/Hilbert/).                              //
//                                                                            //
// the tet core sorts the points using the method in the paper of Boissonnat,J.-D., //
// Devillers, O. and Hornus, S. "Incremental Construction of the Delaunay     //
// Triangulation and the Delaunay Graph in Medium Dimension," In Proceedings  //
// of the 25th ACM Symposium on Computational Geometry, 2009.  It first       //
// randomly sorts the points into subgroups using the Biased Randomized       //
// Insertion Ordering (BRIO) of Amenta et al 2003, then sorts the points in   //
// each subgroup along the 3D Hilbert curve.  Inserting points in this order  //
// ensure a randomized "sprinkling" of the points over the domain, while      //
// sorting of each subset provides locality.                                  //
//                                                                            //
//============================================================================//

  void transfernodes();

  // Point sorting.
  int  transgc[8][3][8], tsb1mod3[8];
  void hilbert_init(int n);
  int  hilbert_split(point* vertexarray, int arraysize, int gc0, int gc1,
                     double, double, double, double, double, double);
  void hilbert_sort3(point* vertexarray, int arraysize, int e, int d,
                     double, double, double, double, double, double, int depth);
  void brio_multiscale_sort(point*,int,int threshold,double ratio,int* depth);

  // Point location.
  unsigned long randomnation(unsigned int choices);
  void randomsample(point searchpt, triface *searchtet);
  enum locateresult locate(point searchpt, triface *searchtet, int chkencflag = 0);

  // Incremental Delaunay construction.
  enum locateresult locate_dt(point searchpt, triface *searchtet);
  int  insert_vertex_bw(point, triface*, insertvertexflags*);
  void initialdelaunay(point pa, point pb, point pc, point pd);
  void incrementaldelaunay(clock_t&);

//============================================================================//
//                                                                            //
// Surface triangulation                                                      //
//                                                                            //
//============================================================================//

  void flipshpush(face*);
  void flip22(face*, int, int);
  void flip31(face*, int);
  long lawsonflip();
  int sinsertvertex(point newpt, face*, face*, int iloc, int bowywat, int);
  int sremovevertex(point delpt, face*, face*, int lawson);

  enum locateresult slocate(point, face*, int, int, int);
  enum interresult sscoutsegment(face*, point, int, int, int);
  void scarveholes(int, double*);
  int triangulate(int, arraypool*, arraypool*, int, double*);

  void unifysegments();
  void identifyinputedges(point*);
  void mergefacets();
  void meshsurface();


//============================================================================//
//                                                                            //
// Constrained Delaunay tetrahedralization                                    //
//                                                                            //
// A constrained Delaunay tetrahedralization (CDT) is a variation of a Delau- //
// nay tetrahedralization (DT) that respects the boundary of a 3D PLC (mesh   //
// domain).  A crucial difference between a CDT and a DT is that triangles in //
// the PLC's polygons are not required to be locally Delaunay, which frees    //
// the CDT to respect the PLC's polygons better. CDTs have optimal properties //
// similar to those of DTs.                                                   //
//                                                                            //
// Steiner Points and Steiner CDTs. It is well-known that even a simple 3D    //
// polyhedron may not have a tetrahedralization which only uses its vertices. //
// Some extra points, so-called "Steiner points" are needed to form a tetrah- //
// edralization of such polyhedron.  A Steiner CDT of a 3D PLC is a CDT       //
// containing Steiner points. the tet core generates Steiner CDTs.                  //
//                                                                            //
// The routine constraineddelaunay() creates a (Steiner) CDT of the PLC       //
// (including Steiner points). It has two steps, (1) segment recovery and (2) //
// facet (polygon) recovery.                                                  //
//                                                                            //
// The routine delaunizesegments() implements the segment recovery algorithm  //
// of Si, H., and Gaertner, K. "Meshing Piecewise Linear Complexes by         //
// Constrained Delaunay Tetrahedralizations," In Proceedings of the 14th      //
// International Meshing Roundtable, 147--163, 2005.  It adds Steiner points  //
// into non-Delaunay segments until all subsegments appear together in a DT.  //
// The running time of this algorithm is proportional to the number of        //
// Steiner points.                                                            //
//                                                                            //
// There are two incremental facet recovery algorithms: the cavity re-        //
// triangulation algorithm of Si, H., and Gaertner, K. "3D Boundary Recovery  //
// by Constrained Delaunay Tetrahedralization," International Journal for     //
// Numerical Methods in Engineering, 85:1341-1364, 2011, and the flip         //
// algorithm of Shewchuk, J. "Updating and Constructing Constrained Delaunay  //
// and Constrained Regular Triangulations by Flips." In Proceedings of the    //
// 19th ACM Symposium on Computational Geometry, 86-95, 2003.                 //
//                                                                            //
// Although no Steiner point is needed in step (2), a facet with non-coplanar //
// vertices might need Steiner points. It is discussed in the paper of Si, H.,//
// and  Shewchuk, J., "Incrementally Constructing and Updating Constrained    //
// Delaunay Tetrahedralizations with Finite Precision Coordinates." In        //
// Proceedings of the 21th International Meshing Roundtable, 2012.            //
//                                                                            //
// Our implementation of the facet recovery algorithms recovers a "missing    //
// region" at a time. Each missing region is a subset of connected interiors  //
// of a polygon. The routine formcavity() creates the cavity of crossing      //
// tetrahedra of the missing region. The cavity re-triangulation algorithm is //
// implemented by three subroutines, delaunizecavity(), fillcavity(), and     //
// carvecavity(). Since it may fail due to non-coplanar vertices, the         //
// subroutine restorecavity() is used to restore the original cavity.         //
//                                                                            //
// The routine flipinsertfacet() implements the flip algorithm. The sub-      //
// routine flipcertify() is used to maintain the priority queue of flips.     //
// The routine refineregion() is called when the facet recovery algorithm     //
// fails to recover a missing region. It inserts Steiner points to refine the //
// missing region. To avoid inserting Steiner points very close to existing   //
// segments.  The classical encroachment rules of the Delaunay refinement     //
// algorithm are used to choose the Steiner points.  The routine              //
// constrainedfacets() does the facet recovery by using either the cavity re- //
// triangulation algorithm (default) or the flip algorithm. It results in a   //
// CDT of the (modified) PLC (including Steiner points).                      //
//                                                                            //
//============================================================================//

  enum interresult finddirection(triface* searchtet, point endpt);
  enum interresult scoutsegment(point, point, face*, triface*, point*, 
                                arraypool*);
  int  getsteinerptonsegment(face* seg, point refpt, point steinpt);
  void delaunizesegments();

  int  scoutsubface(face* searchsh,triface* searchtet,int shflag);
  void formregion(face*, arraypool*, arraypool*, arraypool*);
  int  scoutcrossedge(triface& crosstet, arraypool*, arraypool*);
  bool formcavity(triface*, arraypool*, arraypool*, arraypool*, arraypool*, 
                  arraypool*, arraypool*);
  // Facet recovery by cavity re-triangulation [Si and Gaertner 2011].
  void delaunizecavity(arraypool*, arraypool*, arraypool*, arraypool*, 
                       arraypool*, arraypool*);
  bool fillcavity(arraypool*, arraypool*, arraypool*, arraypool*,
                  arraypool*, arraypool*, triface* crossedge);
  void carvecavity(arraypool*, arraypool*, arraypool*);
  void restorecavity(arraypool*, arraypool*, arraypool*, arraypool*);
  // Facet recovery by flips [Shewchuk 2003].
  void flipcertify(triface *chkface, badface **pqueue, point, point, point);
  void flipinsertfacet(arraypool*, arraypool*, arraypool*, arraypool*);

  int  insertpoint_cdt(point, triface*, face*, face*, insertvertexflags*,
                       arraypool*, arraypool*, arraypool*, arraypool*,
                       arraypool*, arraypool*);
  void refineregion(face&, arraypool*, arraypool*, arraypool*, arraypool*,
                    arraypool*, arraypool*);
  void constrainedfacets();  

  void constraineddelaunay(clock_t&);

//============================================================================//
//                                                                            //
// Constrained tetrahedralizations.                                           //
//                                                                            //
//============================================================================//

  void sort_2pts(point p1, point p2, point ppt[2]);
  void sort_3pts(point p1, point p2, point p3, point ppt[3]);

  bool is_collinear_at(point mid, point left, point right);
  bool is_segment(point p1, point p2);
  bool valid_constrained_f23(triface&, point pd, point pe);
  bool valid_constrained_f32(triface*, point pa, point pb);
  
  int checkflipeligibility(int fliptype, point, point, point, point, point,
                           int level, int edgepivot, flipconstraints* fc);

  int removeedgebyflips(triface*, flipconstraints*);
  int removefacebyflips(triface*, flipconstraints*);

  int recoveredgebyflips(point, point, face*, triface*, int fullsearch, int& idir);
  int add_steinerpt_in_schoenhardtpoly(triface*, int, int, int chkencflag);
  int add_steinerpt_in_segment(face*, int searchlevel, int& idir); 
  int add_steinerpt_to_recover_edge(point, point, face*, int, int, int& idir);
  int recoversegments(arraypool*, int fullsearch, int steinerflag);

  int recoverfacebyflips(point,point,point,face*,triface*,int&,point*,point*);
  int recoversubfaces(arraypool*, int steinerflag);

  int getvertexstar(int, point searchpt, arraypool*, arraypool*, arraypool*);
  int getedge(point, point, triface*);
  int reduceedgesatvertex(point startpt, arraypool* endptlist);
  int removevertexbyflips(point steinerpt);

  int smoothpoint(point smtpt, arraypool*, int ccw, optparameters *opm);
  int suppressbdrysteinerpoint(point steinerpt);
  int suppresssteinerpoints();

  void recoverboundary(clock_t&);

//============================================================================//
//                                                                            //
// Mesh reconstruction                                                        //
//                                                                            //
//============================================================================//

  void carveholes();

  void reconstructmesh();

  int  search_face(point p0, point p1, point p2, triface &tetloop);
  int  search_edge(point p0, point p1, triface &tetloop);
  int  scout_point(point, triface*, int randflag);
  double getpointmeshsize(point, triface*, int iloc);
  void interpolatemeshsize();

  void insertconstrainedpoints(point *insertarray, int arylen, int rejflag);
  void insertconstrainedpoints(TetMeshData *addio);

  void collectremovepoints(arraypool *remptlist);
  void meshcoarsening();

//============================================================================//
//                                                                            //
// Mesh refinement                                                            //
//                                                                            //
// The purpose of mesh refinement is to obtain a tetrahedral mesh with well-  //
// -shaped tetrahedra and appropriate mesh size.  It is necessary to insert   //
// new Steiner points to achieve this property. The questions are (1) how to  //
// choose the Steiner points? and (2) how to insert them?                     //
//                                                                            //
// Delaunay refinement is a technique first developed by Chew [1989] and      //
// Ruppert [1993, 1995] to generate quality triangular meshes in the plane.   //
// It provides guarantee on the smallest angle of the triangles.  Rupper's    //
// algorithm guarantees that the mesh is size-optimal (to within a constant   //
// factor) among all meshes with the same quality.                            //
//   Shewchuk generalized Ruppert's algorithm into 3D in his PhD thesis       //
// [Shewchuk 1997]. A short version of his algorithm appears in "Tetrahedral  //
// Mesh Generation by Delaunay Refinement," In Proceedings of the 14th ACM    //
// Symposium on Computational Geometry, 86-95, 1998.  It guarantees that all  //
// tetrahedra of the output mesh have a "radius-edge ratio" (equivalent to    //
// the minimal face angle) bounded. However, it does not remove slivers, a    //
// type of very flat tetrahedra which can have no small face angles but have  //
// very small (and large) dihedral angles. Moreover, it may not terminate if  //
// the input PLC contains "sharp features", e.g., two edges (or two facets)   //
// meet at an acute angle (or dihedral angle).                                //
//                                                                            //
// the tet core uses the basic Delaunay refinement scheme to insert Steiner points. //
// While it always maintains a constrained Delaunay mesh.  The algorithm is   //
// described in Si, H., "Adaptive Constrained Delaunay Mesh Generation,"      //
// International Journal for Numerical Methods in Engineering, 75:856-880.    //
// This algorithm always terminates and sharp features are easily preserved.  //
// The mesh has good quality (same as Shewchuk's Delaunay refinement algori-  //
// thm) in the bulk of the mesh domain. Moreover, it supports the generation  //
// of adaptive mesh according to a (isotropic) mesh sizing function.          //
//                                                                            //
//============================================================================//

  void makesegmentendpointsmap();
  double set_ridge_vertex_protecting_ball(point);
  double get_min_angle_at_ridge_vertex(face* seg);
  double get_min_diahedral_angle(face* seg);
  void create_segment_info_list();

  void makefacetverticesmap();
  void create_segment_facet_map();

  int  ridge_vertices_adjacent(point, point);
  int  facet_ridge_vertex_adjacent(face *, point);
  int  segsegadjacent(face *, face *);
  int  segfacetadjacent(face *checkseg, face *checksh);
  int  facetfacetadjacent(face *, face *);
  bool is_sharp_segment(face* seg);
  bool does_seg_contain_acute_vertex(face* seg);
  bool create_a_shorter_edge(point steinerpt, point nearpt);

  void enqueuesubface(memorypool*, face*);
  void enqueuetetrahedron(triface*);

  bool check_encroachment(point pa, point pb, point checkpt);
  bool check_enc_segment(face *chkseg, point *pencpt);  
  bool get_steiner_on_segment(face* seg, point encpt, point newpt);
  bool split_segment(face *splitseg, point encpt, double *param, int qflag, int, int*);
  void repairencsegs(double *param, int qflag, int chkencflag);

  bool get_subface_ccent(face *chkfac, double *ccent);
  bool check_enc_subface(face *chkfac, point *pencpt, double *ccent, double *radius);
  bool check_subface(face *chkfac, double *ccent, double radius, double *param);
  void enqueue_subface(face *bface, point encpt, double *ccent, double *param);
  badface* top_subface();
  void dequeue_subface();
  void parallel_shift(point pa, point pb, point pc, point pt, double* ppt);
  enum locateresult locate_on_surface(point searchpt, face* searchsh);
  bool split_subface(face *splitfac, point encpt, double *ccent, double*, int, int, int*);
  void repairencfacs(double *param, int qflag, int chkencflag);

  bool check_tetrahedron(triface *chktet, double* param, int& qflag);
  bool checktet4split(triface *chktet, double* param, int& qflag);
  enum locateresult locate_point_walk(point searchpt, triface*, int chkencflag);
  bool split_tetrahedron(triface*, double*, int, int, insertvertexflags &ivf);
  void repairbadtets(double queratio, int chkencflag);

  void delaunayrefinement();

//============================================================================//
//                                                                            //
// Mesh optimization                                                          //
//                                                                            //
//============================================================================//

  long lawsonflip3d(flipconstraints *fc);
  void recoverdelaunay();

  int  get_seg_laplacian_center(point mesh_vert, double target[3]);
  int  get_surf_laplacian_center(point mesh_vert, double target[3]);
  int  get_laplacian_center(point mesh_vert, double target[3]);
  bool move_vertex(point mesh_vert, double target[3]);
  void smooth_vertices();

  bool get_tet(point, point, point, point, triface *);
  bool get_tetqual(triface *chktet, point oppo_pt, badface *bf);
  bool get_tetqual(point, point, point, point, badface *bf);
  void enqueue_badtet(badface *bf);
  badface* top_badtet();
  void dequeue_badtet();

  bool add_steinerpt_to_repair(badface *bf, bool bSmooth);
  bool flip_edge_to_improve(triface *sliver_edge, double& improved_cosmaxd);
  bool repair_tet(badface *bf, bool bFlips, bool bSmooth, bool bSteiners);
  long repair_badqual_tets(bool bFlips, bool bSmooth, bool bSteiners);
  void improve_mesh();

//============================================================================//
//                                                                            //
// Mesh check and statistics                                                  //
//                                                                            //
//============================================================================//

  // Mesh validations.
  int check_mesh(int topoflag);
  int check_shells();
  int check_segments();
  int check_delaunay(int perturb = 1);
  int check_regular(int);
  int check_conforming(int);

  //  Mesh statistics.
  void printfcomma(unsigned long n);
  void qualitystatistics();
  void memorystatistics();
  void statistics();

//============================================================================//
//                                                                            //
// Mesh output                                                                //
//                                                                            //
//============================================================================//

  void jettisonnodes();
  void highorder();
  void indexelements();
  void numberedges();
  void outnodes(TetMeshData*);
  void outmetrics(TetMeshData*);
  void outelements(TetMeshData*);
  void outfaces(TetMeshData*);
  void outhullfaces(TetMeshData*);
  void outsubfaces(TetMeshData*);
  void outedges(TetMeshData*);
  void outsubsegments(TetMeshData*);
  void outneighbors(TetMeshData*);
  void outvoronoi(TetMeshData*);




//============================================================================//
//                                                                            //
// Constructor & destructor                                                   //
//                                                                            //
//============================================================================//

  void initialize_tet_mesh_core()
  {
    in  = addin = NULL;
    b   = NULL;
    bgm = NULL;

    tetrahedrons = subfaces = subsegs = points = NULL;
    tet2segpool = tet2subpool = NULL;
    dummypoint = NULL;

    badtetrahedrons = badsubfacs = badsubsegs = NULL;
    split_segments_pool = split_subfaces_pool = NULL;
    unsplit_badtets = unsplit_subfaces = unsplit_segments = NULL;
    check_tets_list = NULL;
    badqual_tets_pool = NULL;

    stack_enc_segments = stack_enc_subfaces = NULL;
  
    flippool = NULL;
    flipstack = unflip_queue_front = unflip_queue_tail = NULL;
    later_unflip_queue = unflipqueue = NULL;

    cavetetlist = cavebdrylist = caveoldtetlist = NULL;
    cave_oldtet_list = NULL;
    cavetetshlist = cavetetseglist = cavetetvertlist = NULL;
    caveencshlist = caveencseglist = NULL;
    caveshlist = caveshbdlist = cavesegshlist = NULL;

    subsegstack = subfacstack = subvertstack = NULL;
    skipped_segment_list = skipped_facet_list = NULL;
    
    encseglist = encshlist = NULL;
    
    number_of_facets = 0;
    idx2facetlist = NULL;
    facetverticeslist = NULL;
    idx_segment_facet_list = NULL;
    segment_facet_list = NULL;
    idx_ridge_vertex_facet_list = NULL;
    ridge_vertex_facet_list = NULL;

    segmentendpointslist_length = 0;
    segmentendpointslist = NULL;
    segment_info_list = NULL;
    idx_segment_ridge_vertex_list = NULL;
    segment_ridge_vertex_list = NULL;

    subdomains = 0;
    subdomain_markers = NULL;

    numpointattrib = numelemattrib = 0;
    sizeoftensor = 0;
    pointmtrindex = 0;
    pointparamindex = 0;
    pointmarkindex = 0;
    point2simindex = 0;
    pointinsradiusindex = 0;
    elemattribindex = 0;
    polarindex = 0;
    volumeboundindex = 0;
    shmarkindex = 0;
    areaboundindex = 0;
    checksubsegflag = 0;
    checksubfaceflag = 0;
    boundary_recovery_flag = 0;
    checkconstraints = 0;
    nonconvex = 0;
    autofliplinklevel = 1;
    useinsertradius = 0;
    samples = 0l;
    randomseed = 1l;
    minfaceang = minfacetdihed = PI;
    cos_facet_separate_ang_tol = cos(179.9/180.*PI);
    cos_collinear_ang_tol = cos(179.9/180.*PI);
    tetprism_vol_sum = 0.0;
    longest = minedgelength = 0.0;
    xmax = xmin = ymax = ymin = zmax = zmin = 0.0;

    smallest_insradius = 1.e+30;
    big_radius_edge_ratio = 100.0;
    elem_limit = 0;
    insert_point_count = 0l;
    report_refine_progress = 0l;
    last_point_count = 0l;
    last_insertion_count = 0l;

    insegments = 0l;
    hullsize = 0l;
    meshedges = meshhulledges = 0l;
    steinerleft = -1;
    dupverts = 0l;
    unuverts = 0l;
    duplicated_facets_count = 0l;
    nonregularcount = 0l;
    st_segref_count = st_facref_count = st_volref_count = 0l;
    fillregioncount = cavitycount = cavityexpcount = 0l;
    flip14count = flip26count = flipn2ncount = 0l;
    flip23count = flip32count = flip44count = flip41count = 0l;
    flip22count = flip31count = 0l;
    recover_delaunay_count = 0l;
    opt_flips_count = opt_collapse_count = opt_smooth_count = 0l;
    totalworkmemory = 0l;

  } // TetMeshCore()

  void freememory()
  {
    if (bgm != NULL) {
      delete bgm;
    }

    if (points != (memorypool *) NULL) {
      delete points;
      delete [] dummypoint;
    }
    if (tetrahedrons != (memorypool *) NULL) {
      delete tetrahedrons;
    }
    if (subfaces != (memorypool *) NULL) {
      delete subfaces;
      delete subsegs;
    }
    if (tet2segpool != NULL) {
      delete tet2segpool;
      delete tet2subpool;
    }

    if (badtetrahedrons) {
      delete badtetrahedrons;
    }
    if (badsubfacs) {
      delete badsubfacs;
    }
    if (badsubsegs) {
      delete badsubsegs;
    }
    if (unsplit_badtets) {
      delete unsplit_badtets;
    }
    if (check_tets_list) {
      delete check_tets_list;
    }

    if (flippool != NULL) {
      delete flippool;
      delete later_unflip_queue;
      delete unflipqueue;
    }

    if (cavetetlist != NULL) {
      delete cavetetlist;
      delete cavebdrylist;
      delete caveoldtetlist;
      delete cavetetvertlist;
      delete cave_oldtet_list;
    }

    if (caveshlist != NULL) {
      delete caveshlist;
      delete caveshbdlist;
      delete cavesegshlist;
      delete cavetetshlist;
      delete cavetetseglist;
      delete caveencshlist;
      delete caveencseglist;
    }

    if (subsegstack != NULL) {
      delete subsegstack;
      delete subfacstack;
      delete subvertstack;
    }

    if (idx2facetlist != NULL) {
      delete [] idx2facetlist;
      delete [] facetverticeslist;
      delete [] idx_segment_facet_list;
      delete [] segment_facet_list;
      delete [] idx_ridge_vertex_facet_list;
      delete [] ridge_vertex_facet_list;
    }

    if (segmentendpointslist != NULL) {
      delete [] segmentendpointslist;
      delete [] idx_segment_ridge_vertex_list;
      delete [] segment_ridge_vertex_list;
    }

    if (segment_info_list != NULL) {
      delete [] segment_info_list;
    }

    if (subdomain_markers != NULL) {
      delete [] subdomain_markers;
    }

    initialize_tet_mesh_core();
  }

  TetMeshCore()
  {
    initialize_tet_mesh_core();
  }

  ~TetMeshCore()
  {
    freememory();
  } // ~TetMeshCore()

};                                               // End of class TetMeshCore.

//============================================================================//
//                                                                            //
// run_tet_mesh_core()    Interface for using the tet core's library to generate       //
//                     Delaunay tetrahedralizations, constrained Delaunay     //
//                     tetrahedralizations, quality tetrahedral meshes.       //
//                                                                            //
// 'in' is an object of 'TetMeshData' containing a PLC or a previously generated //
// tetrahedral mesh you want to refine.  'out' is another object of 'TetMeshData'//
// for returing the generated tetrahedral mesh. If it is a NULL pointer, the  //
// output mesh is saved to file(s). If 'bgmin' != NULL, it contains a back-   //
// ground mesh defining a mesh size function.                                 //
//                                                                            //
//============================================================================//

void run_tet_mesh_core(TetMeshBehavior *b, TetMeshData *in, TetMeshData *out,
                    TetMeshData *addin = NULL, TetMeshData *bgmin = NULL);

void run_tet_mesh_core(char *switches, TetMeshData *in, TetMeshData *out,
                    TetMeshData *addin = NULL, TetMeshData *bgmin = NULL);

//============================================================================//
//                                                                            //
// terminate_tet_core()    Terminate the tet core with a given exit code.              //
//                                                                            //
//============================================================================//


inline void terminate_tet_core(TetMeshCore * /*m*/, int x)
{
  throw x;
}

//============================================================================//
//                                                                            //
// Primitives for tetrahedra                                                  //
//                                                                            //
//============================================================================//

// encode()  compress a handle into a single pointer.  It relies on the 
//   assumption that all addresses of tetrahedra are aligned to sixteen-
//   byte boundaries, so that the last four significant bits are zero.

inline TetMeshCore::tetrahedron TetMeshCore::encode(triface& t) {
  return (tetrahedron) ((uintptr_t) (t).tet | (uintptr_t) (t).ver);
}

inline TetMeshCore::tetrahedron TetMeshCore::encode2(tetrahedron* ptr, int ver) {
  return (tetrahedron) ((uintptr_t) (ptr) | (uintptr_t) (ver));
}

// decode()  converts a pointer to a handle. The version is extracted from
//   the four least significant bits of the pointer.

inline void TetMeshCore::decode(tetrahedron ptr, triface& t) {
  (t).ver = (int) ((uintptr_t) (ptr) & (uintptr_t) 15);
  (t).tet = (tetrahedron *) ((uintptr_t) (ptr) ^ (uintptr_t) (t).ver);
}

inline TetMeshCore::tetrahedron* TetMeshCore::decode_tet_only(tetrahedron ptr)
{
  return (tetrahedron *) ((((uintptr_t) ptr) >> 4) << 4);
}

inline int TetMeshCore::decode_ver_only(tetrahedron ptr)
{
  return (int) ((uintptr_t) (ptr) & (uintptr_t) 15);
}

// bond()  connects two tetrahedra together. (t1,v1) and (t2,v2) must 
//   refer to the same face and the same edge. 

inline void TetMeshCore::bond(triface& t1, triface& t2) {
  t1.tet[t1.ver & 3] = encode2(t2.tet, bondtbl[t1.ver][t2.ver]);
  t2.tet[t2.ver & 3] = encode2(t1.tet, bondtbl[t2.ver][t1.ver]);
}


// dissolve()  a bond (from one side).

inline void TetMeshCore::dissolve(triface& t) {
  t.tet[t.ver & 3] = NULL;
}

// enext()  finds the next edge (counterclockwise) in the same face.

inline void TetMeshCore::enext(triface& t1, triface& t2) {
  t2.tet = t1.tet;
  t2.ver = enexttbl[t1.ver]; // (t1.ver + 4) % 12;
}

inline void TetMeshCore::enextself(triface& t) {
  t.ver = enexttbl[t.ver]; // (t.ver + 4) % 12;
}

// eprev()   finds the next edge (clockwise) in the same face.

inline void TetMeshCore::eprev(triface& t1, triface& t2) {
  t2.tet = t1.tet;
  t2.ver = eprevtbl[t1.ver]; // (t1.ver + 8) % 12;
}

inline void TetMeshCore::eprevself(triface& t) {
  t.ver = eprevtbl[t.ver]; // (t.ver + 8) % 12;
}

// esym()  finds the reversed edge.  It is in the other face of the
//   same tetrahedron.

inline void TetMeshCore::esym(triface& t1, triface& t2) {
  (t2).tet = (t1).tet;
  (t2).ver = esymtbl[(t1).ver];
}

inline void TetMeshCore::esymself(triface& t) {
  (t).ver = esymtbl[(t).ver];
}

// enextesym()  finds the reversed edge of the next edge. It is in the other
//   face of the same tetrahedron. It is the combination esym() * enext(). 

inline void TetMeshCore::enextesym(triface& t1, triface& t2) {
  t2.tet = t1.tet;
  t2.ver = enextesymtbl[t1.ver];
}

inline void TetMeshCore::enextesymself(triface& t) {
  t.ver = enextesymtbl[t.ver];
}

// eprevesym()  finds the reversed edge of the previous edge.

inline void TetMeshCore::eprevesym(triface& t1, triface& t2) {
  t2.tet = t1.tet;
  t2.ver = eprevesymtbl[t1.ver];
}

inline void TetMeshCore::eprevesymself(triface& t) {
  t.ver = eprevesymtbl[t.ver];
}

// eorgoppo()    Finds the opposite face of the origin of the current edge.
//               Return the opposite edge of the current edge.

inline void TetMeshCore::eorgoppo(triface& t1, triface& t2) {
  t2.tet = t1.tet;
  t2.ver = eorgoppotbl[t1.ver];
}

inline void TetMeshCore::eorgoppoself(triface& t) {
  t.ver = eorgoppotbl[t.ver];
}

// edestoppo()    Finds the opposite face of the destination of the current 
//                edge. Return the opposite edge of the current edge.

inline void TetMeshCore::edestoppo(triface& t1, triface& t2) {
  t2.tet = t1.tet;
  t2.ver = edestoppotbl[t1.ver];
}

inline void TetMeshCore::edestoppoself(triface& t) {
  t.ver = edestoppotbl[t.ver];
}

// fsym()  finds the adjacent tetrahedron at the same face and the same edge.

inline void TetMeshCore::fsym(triface& t1, triface& t2) {
  decode((t1).tet[(t1).ver & 3], t2);
  t2.ver = fsymtbl[t1.ver][t2.ver];
}


#define fsymself(t) \
  t1ver = (t).ver; \
  decode((t).tet[(t).ver & 3], (t));\
  (t).ver = fsymtbl[t1ver][(t).ver]

// fnext()  finds the next face while rotating about an edge according to
//   a right-hand rule. The face is in the adjacent tetrahedron.  It is
//   the combination: fsym() * esym().

inline void TetMeshCore::fnext(triface& t1, triface& t2) {
  decode(t1.tet[facepivot1[t1.ver]], t2);
  t2.ver = facepivot2[t1.ver][t2.ver];
}


#define fnextself(t) \
  t1ver = (t).ver; \
  decode((t).tet[facepivot1[(t).ver]], (t)); \
  (t).ver = facepivot2[t1ver][(t).ver]


// The following primtives get or set the origin, destination, face apex,
//   or face opposite of an ordered tetrahedron.

inline TetMeshCore::point TetMeshCore::org(triface& t) {
  return (point) (t).tet[orgpivot[(t).ver]];
}

inline TetMeshCore::point TetMeshCore:: dest(triface& t) {
  return (point) (t).tet[destpivot[(t).ver]];
}

inline TetMeshCore::point TetMeshCore:: apex(triface& t) {
  return (point) (t).tet[apexpivot[(t).ver]];
}

inline TetMeshCore::point TetMeshCore:: oppo(triface& t) {
  return (point) (t).tet[oppopivot[(t).ver]];
}

inline void TetMeshCore:: setorg(triface& t, point p) {
  (t).tet[orgpivot[(t).ver]] = (tetrahedron) (p);
}

inline void TetMeshCore:: setdest(triface& t, point p) {
  (t).tet[destpivot[(t).ver]] = (tetrahedron) (p);
}

inline void TetMeshCore:: setapex(triface& t, point p) {
  (t).tet[apexpivot[(t).ver]] = (tetrahedron) (p);
}

inline void TetMeshCore:: setoppo(triface& t, point p) {
  (t).tet[oppopivot[(t).ver]] = (tetrahedron) (p);
}

#define setvertices(t, torg, tdest, tapex, toppo) \
  (t).tet[orgpivot[(t).ver]] = (tetrahedron) (torg);\
  (t).tet[destpivot[(t).ver]] = (tetrahedron) (tdest); \
  (t).tet[apexpivot[(t).ver]] = (tetrahedron) (tapex); \
  (t).tet[oppopivot[(t).ver]] = (tetrahedron) (toppo)


inline double* TetMeshCore::get_polar(tetrahedron* ptr)
{
  return &(((double *) (ptr))[polarindex]);
}
inline double TetMeshCore::get_volume(tetrahedron* ptr)
{
  return ((double *) (ptr))[polarindex + 4];
}

// Check or set a tetrahedron's attributes.

inline double TetMeshCore::elemattribute(tetrahedron* ptr, int attnum) {
  return ((double *) (ptr))[elemattribindex + attnum];
}

inline void TetMeshCore::setelemattribute(tetrahedron* ptr, int attnum, 
  double value) {
  ((double *) (ptr))[elemattribindex + attnum] = value;
}

// Check or set a tetrahedron's maximum volume bound.

inline double TetMeshCore::volumebound(tetrahedron* ptr) {
  return ((double *) (ptr))[volumeboundindex];
}

inline void TetMeshCore::setvolumebound(tetrahedron* ptr, double value) {
  ((double *) (ptr))[volumeboundindex] = value;
}

// Get or set a tetrahedron's index (only used for output).
//    These two routines use the reserved slot ptr[10].

inline int TetMeshCore::elemindex(tetrahedron* ptr) {
  int *iptr = (int *) &(ptr[10]);
  return iptr[0];
}

inline void TetMeshCore::setelemindex(tetrahedron* ptr, int value) {
  int *iptr = (int *) &(ptr[10]);
  iptr[0] = value;
}

// Get or set a tetrahedron's marker. 
//   Set 'value = 0' cleans all the face/edge flags.

inline int TetMeshCore::elemmarker(tetrahedron* ptr) {
  return ((int *) (ptr))[elemmarkerindex];
}

inline void TetMeshCore::setelemmarker(tetrahedron* ptr, int value) {
  ((int *) (ptr))[elemmarkerindex] = value;
}

// infect(), infected(), uninfect() -- primitives to flag or unflag a
//   tetrahedron. The last bit of the element marker is flagged (1)
//   or unflagged (0).

inline void TetMeshCore::infect(triface& t) {
  ((int *) (t.tet))[elemmarkerindex] |= 1;
}

inline void TetMeshCore::uninfect(triface& t) {
  ((int *) (t.tet))[elemmarkerindex] &= ~1;
}

inline bool TetMeshCore::infected(triface& t) {
  return (((int *) (t.tet))[elemmarkerindex] & 1) != 0;
}

// marktest(), marktested(), unmarktest() -- primitives to flag or unflag a
//   tetrahedron.  Use the second lowerest bit of the element marker.

inline void TetMeshCore::marktest(triface& t) {
  ((int *) (t.tet))[elemmarkerindex] |= 2;
}

inline void TetMeshCore::unmarktest(triface& t) {
  ((int *) (t.tet))[elemmarkerindex] &= ~2;
}
    
inline bool TetMeshCore::marktested(triface& t) {
  return (((int *) (t.tet))[elemmarkerindex] & 2) != 0;
}

// markface(), unmarkface(), facemarked() -- primitives to flag or unflag a
//   face of a tetrahedron.  From the last 3rd to 6th bits are used for
//   face markers, e.g., the last third bit corresponds to loc = 0. 

inline void TetMeshCore::markface(triface& t) {
  ((int *) (t.tet))[elemmarkerindex] |= (4 << (t.ver & 3));
}

inline void TetMeshCore::unmarkface(triface& t) {
  ((int *) (t.tet))[elemmarkerindex] &= ~(4 << (t.ver & 3));
}

inline bool TetMeshCore::facemarked(triface& t) {
  return (((int *) (t.tet))[elemmarkerindex] & (4 << (t.ver & 3))) != 0;
}

// markedge(), unmarkedge(), edgemarked() -- primitives to flag or unflag an
//   edge of a tetrahedron.  From the last 7th to 12th bits are used for
//   edge markers, e.g., the last 7th bit corresponds to the 0th edge, etc. 
//   Remark: The last 7th bit is marked by 2^6 = 64.

inline void TetMeshCore::markedge(triface& t) {
  ((int *) (t.tet))[elemmarkerindex] |= (int) (64 << ver2edge[(t).ver]);
}

inline void TetMeshCore::unmarkedge(triface& t) {
  ((int *) (t.tet))[elemmarkerindex] &= ~(int) (64 << ver2edge[(t).ver]);
}

inline bool TetMeshCore::edgemarked(triface& t) {
  return (((int *) (t.tet))[elemmarkerindex] & 
           (int) (64 << ver2edge[(t).ver])) != 0;
}

// marktest2(), unmarktest2(), marktest2ed() -- primitives to flag and unflag
//   a tetrahedron. The 13th bit (2^12 = 4096) is used for this flag.

inline void TetMeshCore::marktest2(triface& t) {
  ((int *) (t.tet))[elemmarkerindex] |= (int) (4096);
}

inline void TetMeshCore::unmarktest2(triface& t) {
  ((int *) (t.tet))[elemmarkerindex] &= ~(int) (4096);
}

inline bool TetMeshCore::marktest2ed(triface& t) {
  return (((int *) (t.tet))[elemmarkerindex] & (int) (4096)) != 0;
}

// elemcounter(), setelemcounter() -- primitives to read or ser a (small)
//   integer counter in this tet. It is saved from the 16th bit. On 32 bit
//   system, the range of the counter is [0, 2^15 = 32768]. 

inline int TetMeshCore::elemcounter(triface& t) {
  return (((int *) (t.tet))[elemmarkerindex]) >> 16;
}

inline void TetMeshCore::setelemcounter(triface& t, int value) {
  int c = ((int *) (t.tet))[elemmarkerindex];
  // Clear the old counter while keep the other flags.
  c &= 65535; // sum_{i=0^15} 2^i
  c |= (value << 16);
  ((int *) (t.tet))[elemmarkerindex] = c;
}

inline void TetMeshCore::increaseelemcounter(triface& t) {
  int c = elemcounter(t);
  setelemcounter(t, c + 1);
}

inline void TetMeshCore::decreaseelemcounter(triface& t) {
  int c = elemcounter(t);
  setelemcounter(t, c - 1);
}

// ishulltet()  tests if t is a hull tetrahedron.

inline bool TetMeshCore::ishulltet(triface& t) {
  return (point) (t).tet[7] == dummypoint;
}

// isdeadtet()  tests if t is a tetrahedron is dead.

inline bool TetMeshCore::isdeadtet(triface& t) {
  return ((t.tet == NULL) || (t.tet[4] == NULL));
}

//============================================================================//
//                                                                            //
// Primitives for subfaces and subsegments                                    //
//                                                                            //
//============================================================================//

// Each subface contains three pointers to its neighboring subfaces, with
//   edge versions.  To save memory, both information are kept in a single
//   pointer. To make this possible, all subfaces are aligned to eight-byte
//   boundaries, so that the last three bits of each pointer are zeros. An
//   edge version (in the range 0 to 5) is compressed into the last three
//   bits of each pointer by 'sencode()'.  'sdecode()' decodes a pointer,
//   extracting an edge version and a pointer to the beginning of a subface.

inline void TetMeshCore::sdecode(shellface sptr, face& s) {
  s.shver = (int) ((uintptr_t) (sptr) & (uintptr_t) 7);
  s.sh = (shellface *) ((uintptr_t) (sptr) ^ (uintptr_t) (s.shver));
}

inline TetMeshCore::shellface TetMeshCore::sencode(face& s) {
  return (shellface) ((uintptr_t) s.sh | (uintptr_t) s.shver);
}

inline TetMeshCore::shellface TetMeshCore::sencode2(shellface *sh, int shver) {
  return (shellface) ((uintptr_t) sh | (uintptr_t) shver);
}

// sbond() bonds two subfaces (s1) and (s2) together. s1 and s2 must refer
//   to the same edge. No requirement is needed on their orientations.

inline void TetMeshCore::sbond(face& s1, face& s2) 
{
  s1.sh[s1.shver >> 1] = sencode(s2);
  s2.sh[s2.shver >> 1] = sencode(s1);
}

// sbond1() bonds s1 <== s2, i.e., after bonding, s1 is pointing to s2,
//   but s2 is not pointing to s1.  s1 and s2 must refer to the same edge.
//   No requirement is needed on their orientations.

inline void TetMeshCore::sbond1(face& s1, face& s2) 
{
  s1.sh[s1.shver >> 1] = sencode(s2);
}

// Dissolve a subface bond (from one side).  Note that the other subface
//   will still think it's connected to this subface.

inline void TetMeshCore::sdissolve(face& s)
{
  s.sh[s.shver >> 1] = NULL;
}

// spivot() finds the adjacent subface (s2) for a given subface (s1).
//   s1 and s2 share at the same edge.

inline void TetMeshCore::spivot(face& s1, face& s2) 
{
  shellface sptr = s1.sh[s1.shver >> 1];
  sdecode(sptr, s2);
}

inline void TetMeshCore::spivotself(face& s) 
{
  shellface sptr = s.sh[s.shver >> 1];
  sdecode(sptr, s);
}

// These primitives determine or set the origin, destination, or apex
//   of a subface with respect to the edge version.

inline TetMeshCore::point TetMeshCore::sorg(face& s) 
{
  return (point) s.sh[sorgpivot[s.shver]];
}

inline TetMeshCore::point TetMeshCore::sdest(face& s) 
{
  return (point) s.sh[sdestpivot[s.shver]];
}

inline TetMeshCore::point TetMeshCore::sapex(face& s) 
{
  return (point) s.sh[sapexpivot[s.shver]];
}

inline void TetMeshCore::setsorg(face& s, point pointptr) 
{
  s.sh[sorgpivot[s.shver]] = (shellface) pointptr;
}

inline void TetMeshCore::setsdest(face& s, point pointptr) 
{
  s.sh[sdestpivot[s.shver]] = (shellface) pointptr;
}

inline void TetMeshCore::setsapex(face& s, point pointptr) 
{
  s.sh[sapexpivot[s.shver]] = (shellface) pointptr;
}

#define setshvertices(s, pa, pb, pc)\
  setsorg(s, pa);\
  setsdest(s, pb);\
  setsapex(s, pc)

// sesym()  reserves the direction of the lead edge.

inline void TetMeshCore::sesym(face& s1, face& s2) 
{
  s2.sh = s1.sh;
  s2.shver = (s1.shver ^ 1);  // Inverse the last bit.
}

inline void TetMeshCore::sesymself(face& s) 
{
  s.shver ^= 1;
}

// senext()  finds the next edge (counterclockwise) in the same orientation
//   of this face.

inline void TetMeshCore::senext(face& s1, face& s2) 
{
  s2.sh = s1.sh;
  s2.shver = snextpivot[s1.shver];
}

inline void TetMeshCore::senextself(face& s) 
{
  s.shver = snextpivot[s.shver];
}

inline void TetMeshCore::senext2(face& s1, face& s2) 
{
  s2.sh = s1.sh;
  s2.shver = snextpivot[snextpivot[s1.shver]];
}

inline void TetMeshCore::senext2self(face& s) 
{
  s.shver = snextpivot[snextpivot[s.shver]];
}


// Check or set a subface's maximum area bound.

inline double TetMeshCore::areabound(face& s) 
{
  return ((double *) (s.sh))[areaboundindex];
}

inline void TetMeshCore::setareabound(face& s, double value) 
{
  ((double *) (s.sh))[areaboundindex] = value;
}

// These two primitives read or set a shell marker.  Shell markers are used
//   to hold user boundary information.

inline int TetMeshCore::shellmark(face& s) 
{
  return ((int *) (s.sh))[shmarkindex];
}

inline void TetMeshCore::setshellmark(face& s, int value) 
{
  ((int *) (s.sh))[shmarkindex] = value;
}



// sinfect(), sinfected(), suninfect() -- primitives to flag or unflag a
//   subface. The last bit of ((int *) ((s).sh))[shmarkindex+1] is flagged.

inline void TetMeshCore::sinfect(face& s) 
{
  ((int *) ((s).sh))[shmarkindex+1] = 
    (((int *) ((s).sh))[shmarkindex+1] | (int) 1);
}

inline void TetMeshCore::suninfect(face& s) 
{
  ((int *) ((s).sh))[shmarkindex+1] = 
    (((int *) ((s).sh))[shmarkindex+1] & ~(int) 1);
}

// Test a subface for viral infection.

inline bool TetMeshCore::sinfected(face& s) 
{
  return (((int *) ((s).sh))[shmarkindex+1] & (int) 1) != 0;
}

// smarktest(), smarktested(), sunmarktest() -- primitives to flag or unflag
//   a subface. The last 2nd bit of the integer is flagged.

inline void TetMeshCore::smarktest(face& s) 
{
  ((int *) ((s).sh))[shmarkindex+1] = 
    (((int *)((s).sh))[shmarkindex+1] | (int) 2);
}

inline void TetMeshCore::sunmarktest(face& s) 
{
  ((int *) ((s).sh))[shmarkindex+1] = 
    (((int *)((s).sh))[shmarkindex+1] & ~(int)2);
}

inline bool TetMeshCore::smarktested(face& s) 
{
  return ((((int *) ((s).sh))[shmarkindex+1] & (int) 2) != 0);
}

// smarktest2(), smarktest2ed(), sunmarktest2() -- primitives to flag or 
//   unflag a subface. The last 3rd bit of the integer is flagged.

inline void TetMeshCore::smarktest2(face& s) 
{
  ((int *) ((s).sh))[shmarkindex+1] = 
    (((int *)((s).sh))[shmarkindex+1] | (int) 4);
}

inline void TetMeshCore::sunmarktest2(face& s) 
{
  ((int *) ((s).sh))[shmarkindex+1] = 
    (((int *)((s).sh))[shmarkindex+1] & ~(int)4);
}

inline bool TetMeshCore::smarktest2ed(face& s) 
{
  return ((((int *) ((s).sh))[shmarkindex+1] & (int) 4) != 0);
}

// The last 4th bit of ((int *) ((s).sh))[shmarkindex+1] is flagged.

inline void TetMeshCore::smarktest3(face& s) 
{
  ((int *) ((s).sh))[shmarkindex+1] = 
    (((int *)((s).sh))[shmarkindex+1] | (int) 8);
}

inline void TetMeshCore::sunmarktest3(face& s) 
{
  ((int *) ((s).sh))[shmarkindex+1] = 
    (((int *)((s).sh))[shmarkindex+1] & ~(int)8);
}

inline bool TetMeshCore::smarktest3ed(face& s) 
{
  return ((((int *) ((s).sh))[shmarkindex+1] & (int) 8) != 0);
}


// Each facet has a unique index (automatically indexed). Starting from '0'.
// We save this index in the same field of the shell type. 

inline void TetMeshCore::setfacetindex(face& s, int value)
{
  ((int *) (s.sh))[shmarkindex + 2] = value;
}

inline int TetMeshCore::getfacetindex(face& s)
{
  return ((int *) (s.sh))[shmarkindex + 2];
}

// Tests if the subface (subsegment) s is dead.

inline bool TetMeshCore::isdeadsh(face& s) {
  return ((s.sh == NULL) || (s.sh[3] == NULL));
}

//============================================================================//
//                                                                            //
// Primitives for interacting between tetrahedra and subfaces                 //
//                                                                            //
//============================================================================//

// tsbond() bond a tetrahedron (t) and a subface (s) together.
// Note that t and s must be the same face and the same edge. Moreover,
//   t and s have the same orientation. 
// Since the edge number in t and in s can be any number in {0,1,2}. We bond
//   the edge in s which corresponds to t's 0th edge, and vice versa.

inline void TetMeshCore::tsbond(triface& t, face& s)
{
  if ((t).tet[9] == NULL) {
    // Allocate space for this tet.
    (t).tet[9] = (tetrahedron) tet2subpool->alloc();
    // Initialize.
    for (int i = 0; i < 4; i++) {
      ((shellface *) (t).tet[9])[i] = NULL;
    }
  }
  // Bond t <== s.
  ((shellface *) (t).tet[9])[(t).ver & 3] = 
    sencode2((s).sh, tsbondtbl[t.ver][s.shver]);
  // Bond s <== t.
  s.sh[9 + ((s).shver & 1)] = 
    (shellface) encode2((t).tet, stbondtbl[t.ver][s.shver]);
}

// tspivot() finds a subface (s) abutting on the given tetrahdera (t).
//   Return s.sh = NULL if there is no subface at t. Otherwise, return
//   the subface s, and s and t must be at the same edge wth the same
//   orientation.

inline void TetMeshCore::tspivot(triface& t, face& s) 
{
  if ((t).tet[9] == NULL) {
    (s).sh = NULL;
    return;
  }
  // Get the attached subface s.
  sdecode(((shellface *) (t).tet[9])[(t).ver & 3], (s));
  (s).shver = tspivottbl[t.ver][s.shver];
}

// Quickly check if the handle (t, v) is a subface.
#define issubface(t) \
  ((t).tet[9] && ((t).tet[9])[(t).ver & 3])

// stpivot() finds a tetrahedron (t) abutting a given subface (s).
//   Return the t (if it exists) with the same edge and the same
//   orientation of s.

inline void TetMeshCore::stpivot(face& s, triface& t) 
{
  decode((tetrahedron) s.sh[9 + (s.shver & 1)], t);
  if ((t).tet == NULL) {
    return;
  }
  (t).ver = stpivottbl[t.ver][s.shver];
}

// Quickly check if this subface is attached to a tetrahedron.

#define isshtet(s) \
  ((s).sh[9 + ((s).shver & 1)])

// tsdissolve() dissolve a bond (from the tetrahedron side).

inline void TetMeshCore::tsdissolve(triface& t) 
{
  if ((t).tet[9] != NULL) {
    ((shellface *) (t).tet[9])[(t).ver & 3] = NULL;
  }
}

// stdissolve() dissolve a bond (from the subface side).

inline void TetMeshCore::stdissolve(face& s) 
{
  (s).sh[9] = NULL;
  (s).sh[10] = NULL;
}

//============================================================================//
//                                                                            //
// Primitives for interacting between subfaces and segments                   //
//                                                                            //
//============================================================================//

// ssbond() bond a subface to a subsegment.

inline void TetMeshCore::ssbond(face& s, face& edge) 
{
  s.sh[6 + (s.shver >> 1)] = sencode(edge);
  edge.sh[0] = sencode(s);
}

inline void TetMeshCore::ssbond1(face& s, face& edge) 
{
  s.sh[6 + (s.shver >> 1)] = sencode(edge);
  //edge.sh[0] = sencode(s);
}

// ssdisolve() dissolve a bond (from the subface side)

inline void TetMeshCore::ssdissolve(face& s) 
{
  s.sh[6 + (s.shver >> 1)] = NULL;
}

// sspivot() finds a subsegment abutting a subface.

inline void TetMeshCore::sspivot(face& s, face& edge) 
{
  sdecode((shellface) s.sh[6 + (s.shver >> 1)], edge);
}

// Quickly check if the edge is a subsegment.

#define isshsubseg(s) \
  ((s).sh[6 + ((s).shver >> 1)])

//============================================================================//
//                                                                            //
// Primitives for interacting between tetrahedra and segments                 //
//                                                                            //
//============================================================================//

inline void TetMeshCore::tssbond1(triface& t, face& s)
{
  if ((t).tet[8] == NULL) {
    // Allocate space for this tet.
    (t).tet[8] = (tetrahedron) tet2segpool->alloc();
    // Initialization.
    for (int i = 0; i < 6; i++) {
      ((shellface *) (t).tet[8])[i] = NULL;
    }
  }
  ((shellface *) (t).tet[8])[ver2edge[(t).ver]] = sencode((s)); 
}

inline void TetMeshCore::sstbond1(face& s, triface& t) 
{
  ((tetrahedron *) (s).sh)[9] = encode(t);
}

inline void TetMeshCore::tssdissolve1(triface& t)
{
  if ((t).tet[8] != NULL) {
    ((shellface *) (t).tet[8])[ver2edge[(t).ver]] = NULL;
  }
}

inline void TetMeshCore::sstdissolve1(face& s) 
{
  ((tetrahedron *) (s).sh)[9] = NULL;
}

inline void TetMeshCore::tsspivot1(triface& t, face& s)
{
  if ((t).tet[8] != NULL) {
    sdecode(((shellface *) (t).tet[8])[ver2edge[(t).ver]], s);
  } else {
    (s).sh = NULL;
  }
}

// Quickly check whether 't' is a segment or not.

#define issubseg(t) \
  ((t).tet[8] && ((t).tet[8])[ver2edge[(t).ver]])

inline void TetMeshCore::sstpivot1(face& s, triface& t) 
{
  decode((tetrahedron) s.sh[9], t);
}

//============================================================================//
//                                                                            //
// Primitives for points                                                      //
//                                                                            //
//============================================================================//

inline int TetMeshCore::pointmark(point pt) { 
  return ((int *) (pt))[pointmarkindex]; 
}

inline void TetMeshCore::setpointmark(point pt, int value) {
  ((int *) (pt))[pointmarkindex] = value;
}


// These two primitives set and read the type of the point.

inline enum TetMeshCore::verttype TetMeshCore::pointtype(point pt) {
  return (enum verttype) (((int *) (pt))[pointmarkindex + 1] >> (int) 8);
}

inline void TetMeshCore::setpointtype(point pt, enum verttype value) {
  ((int *) (pt))[pointmarkindex + 1] = 
    ((int) value << 8) + (((int *) (pt))[pointmarkindex + 1] & (int) 255);
}

// pinfect(), puninfect(), pinfected() -- primitives to flag or unflag
//   a point. The last bit of the integer '[pointindex+1]' is flagged.

inline void TetMeshCore::pinfect(point pt) {
  ((int *) (pt))[pointmarkindex + 1] |= (int) 1;
}

inline void TetMeshCore::puninfect(point pt) {
  ((int *) (pt))[pointmarkindex + 1] &= ~(int) 1;
}

inline bool TetMeshCore::pinfected(point pt) {
  return (((int *) (pt))[pointmarkindex + 1] & (int) 1) != 0;
}

// pmarktest(), punmarktest(), pmarktested() -- more primitives to 
//   flag or unflag a point. 

inline void TetMeshCore::pmarktest(point pt) {
  ((int *) (pt))[pointmarkindex + 1] |= (int) 2;
}

inline void TetMeshCore::punmarktest(point pt) {
  ((int *) (pt))[pointmarkindex + 1] &= ~(int) 2;
}

inline bool TetMeshCore::pmarktested(point pt) {
  return (((int *) (pt))[pointmarkindex + 1] & (int) 2) != 0;
}

inline void TetMeshCore::pmarktest2(point pt) {
  ((int *) (pt))[pointmarkindex + 1] |= (int) 4;
}

inline void TetMeshCore::punmarktest2(point pt) {
  ((int *) (pt))[pointmarkindex + 1] &= ~(int) 4;
}

inline bool TetMeshCore::pmarktest2ed(point pt) {
  return (((int *) (pt))[pointmarkindex + 1] & (int) 4) != 0;
}

inline void TetMeshCore::pmarktest3(point pt) {
  ((int *) (pt))[pointmarkindex + 1] |= (int) 8;
}

inline void TetMeshCore::punmarktest3(point pt) {
  ((int *) (pt))[pointmarkindex + 1] &= ~(int) 8;
}

inline bool TetMeshCore::pmarktest3ed(point pt) {
  return (((int *) (pt))[pointmarkindex + 1] & (int) 8) != 0;
}

// Read and set the geometry tag of the point (used by -s option).

inline int TetMeshCore::pointgeomtag(point pt) { 
  return ((int *) (pt))[pointmarkindex + 2]; 
}

inline void TetMeshCore::setpointgeomtag(point pt, int value) {
  ((int *) (pt))[pointmarkindex + 2] = value;
}

// Read and set the u,v coordinates of the point (used by -s option).

inline double TetMeshCore::pointgeomuv(point pt, int i) {
  return pt[pointparamindex + i];
}

inline void TetMeshCore::setpointgeomuv(point pt, int i, double value) {
  pt[pointparamindex + i] = value;
}



// These following primitives set and read a pointer to a tetrahedron
//   a subface/subsegment, a point, or a tet of background mesh.

inline TetMeshCore::tetrahedron TetMeshCore::point2tet(point pt) {
  return ((tetrahedron *) (pt))[point2simindex];
}

inline void TetMeshCore::setpoint2tet(point pt, tetrahedron value) {
  ((tetrahedron *) (pt))[point2simindex] = value;
}

inline TetMeshCore::point TetMeshCore::point2ppt(point pt) {
  return (point) ((tetrahedron *) (pt))[point2simindex + 1];
}

inline void TetMeshCore::setpoint2ppt(point pt, point value) {
  ((tetrahedron *) (pt))[point2simindex + 1] = (tetrahedron) value;
}

inline TetMeshCore::shellface TetMeshCore::point2sh(point pt) {
  return (shellface) ((tetrahedron *) (pt))[point2simindex + 2];
}

inline void TetMeshCore::setpoint2sh(point pt, shellface value) {
  ((tetrahedron *) (pt))[point2simindex + 2] = (tetrahedron) value;
}


inline TetMeshCore::tetrahedron TetMeshCore::point2bgmtet(point pt) {
  return ((tetrahedron *) (pt))[point2simindex + 3];
}

inline void TetMeshCore::setpoint2bgmtet(point pt, tetrahedron value) {
  ((tetrahedron *) (pt))[point2simindex + 3] = value;
}


// The primitives for saving and getting the insertion radius.
inline void TetMeshCore::setpointinsradius(point pt, double value)
{
  pt[pointinsradiusindex] = value;
}

inline double TetMeshCore::getpointinsradius(point pt)
{
  return pt[pointinsradiusindex];
}

inline bool TetMeshCore::issteinerpoint(point pt) {
 return (pointtype(pt) == FREESEGVERTEX) || (pointtype(pt) == FREEFACETVERTEX)
        || (pointtype(pt) == FREEVOLVERTEX);
}

// point2tetorg()    Get the tetrahedron whose origin is the point.

inline void TetMeshCore::point2tetorg(point pa, triface& searchtet)
{
  decode(point2tet(pa), searchtet);
  if ((point) searchtet.tet[4] == pa) {
    searchtet.ver = 11;
  } else if ((point) searchtet.tet[5] == pa) {
    searchtet.ver = 3;
  } else if ((point) searchtet.tet[6] == pa) {
    searchtet.ver = 7;
  } else {
    searchtet.ver = 0;
  }
}

// point2shorg()    Get the subface/segment whose origin is the point.

inline void TetMeshCore::point2shorg(point pa, face& searchsh)
{
  sdecode(point2sh(pa), searchsh);
  if ((point) searchsh.sh[3] == pa) {
    searchsh.shver = 0;
  } else if ((point) searchsh.sh[4] == pa) {
    searchsh.shver = (searchsh.sh[5] != NULL ? 2 : 1); 
  } else {
    searchsh.shver = 4;
  }
}

// farsorg()    Return the origin of the subsegment.
// farsdest()   Return the destination of the subsegment.

inline TetMeshCore::point TetMeshCore::farsorg(face& s)
{
  face travesh, neighsh;

  travesh = s;
  while (1) {
    senext2(travesh, neighsh);
    spivotself(neighsh); 
    if (neighsh.sh == NULL) break;
    if (sorg(neighsh) != sorg(travesh)) sesymself(neighsh);
    senext2(neighsh, travesh); 
  }
  return sorg(travesh);
}

inline TetMeshCore::point TetMeshCore::farsdest(face& s) 
{
  face travesh, neighsh;

  travesh = s;
  while (1) {
    senext(travesh, neighsh);
    spivotself(neighsh); 
    if (neighsh.sh == NULL) break;
    if (sdest(neighsh) != sdest(travesh)) sesymself(neighsh);
    senext(neighsh, travesh); 
  }
  return sdest(travesh);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Linear algebra operators.                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// dot() returns the dot product: v1 dot v2.
inline double TetMeshCore::dot(double* v1, double* v2) 
{
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

// cross() computes the cross product: n = v1 cross v2.
inline void TetMeshCore::cross(double* v1, double* v2, double* n) 
{
  n[0] =   v1[1] * v2[2] - v2[1] * v1[2];
  n[1] = -(v1[0] * v2[2] - v2[0] * v1[2]);
  n[2] =   v1[0] * v2[1] - v2[0] * v1[1];
}

// distance() computes the Euclidean distance between two points.
inline double TetMeshCore::distance(double* p1, double* p2)
{
  return sqrt((p2[0] - p1[0]) * (p2[0] - p1[0]) +
              (p2[1] - p1[1]) * (p2[1] - p1[1]) +
              (p2[2] - p1[2]) * (p2[2] - p1[2]));
}

inline double TetMeshCore::distance2(double* p1, double* p2)
{
  return norm2(p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]);
}

inline double TetMeshCore::norm2(double x, double y, double z)
{
  return (x) * (x) + (y) * (y) + (z) * (z);
}



} // namespace sqmesh::mesh::tet::detail

#endif // #ifndef SQMESH_MESH_TET_CORE_HPP

