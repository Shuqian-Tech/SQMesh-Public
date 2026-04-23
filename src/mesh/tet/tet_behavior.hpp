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
// SQMesh tet behavior carrier — extracted from upstream `tet_core.hpp` so that
// SQMesh-side refactoring of the tet pipeline can iterate independently of
// the legacy header. This file deliberately preserves the upstream class
// name `TetMeshBehavior`, member layout and default values so that
// tet core topic files compiles unchanged against this declaration.
//
// As the refactor progresses, callers should migrate to the SQMesh-style
// wrapper that will live in this directory; until then this header is the
// single source of truth for the configuration carrier.

#pragma once

namespace sqmesh::mesh::tet::detail {

class TetMeshBehavior {

public:

  // Switches of the tet core.
  int plc;                                                         // '-p', 0.
  int psc;                                                         // '-s', 0.
  int refine;                                                      // '-r', 0.
  int quality;                                                     // '-q', 0.
  int nobisect;                                                    // '-Y', 0.
  int cdt;                                                         // '-D', 0.
  int cdtrefine;                                                  // '-D#', 7.
  int coarsen;                                                     // '-R', 0.
  int weighted;                                                    // '-w', 0.
  int brio_hilbert;                                                // '-b', 1.
  int flipinsert;                                                  // '-L', 0.
  int metric;                                                      // '-m', 0.
  int varvolume;                                                   // '-a', 0.
  int fixedvolume;                                                 // '-a', 0.
  int regionattrib;                                                // '-A', 0.
  int insertaddpoints;                                             // '-i', 0.
  int diagnose;                                                    // '-d', 0.
  int convex;                                                      // '-c', 0.
  int nomergefacet;                                                // '-M', 0.
  int nomergevertex;                                               // '-M', 0.
  int noexact;                                                     // '-X', 0.
  int nostaticfilter;                                              // '-X', 0.
  int zeroindex;                                                   // '-z', 0.
  int facesout;                                                    // '-f', 0.
  int edgesout;                                                    // '-e', 0.
  int neighout;                                                    // '-n', 0.
  int voroout;                                                     // '-v', 0.
  int meditview;                                                   // '-g', 0.
  int vtkview;                                                     // '-k', 0.
  int vtksurfview;                                                 // '-k', 0.
  int nobound;                                                     // '-B', 0.
  int nonodewritten;                                               // '-N', 0.
  int noelewritten;                                                // '-E', 0.
  int nofacewritten;                                               // '-F', 0.
  int noiterationnum;                                              // '-I', 0.
  int nojettison;                                                  // '-J', 0.
  int docheck;                                                     // '-C', 0.
  int quiet;                                                       // '-Q', 0.
  int nowarning;                                                   // '-W', 0.
  int verbose;                                                     // '-V', 0.

  // Parameters of the tet core.
  int vertexperblock;                                           // '-x', 4092.
  int tetrahedraperblock;                                       // '-x', 8188.
  int shellfaceperblock;                                        // '-x', 2044.
  int supsteiner_level;                                           // '-Y/', 2.
  int addsteiner_algo;                                           // '-Y//', 1.
  int coarsen_param;                                               // '-R', 0.
  int weighted_param;                                              // '-w', 0.
  int fliplinklevel;                                                    // -1.
  int flipstarsize;                                                     // -1.
  int fliplinklevelinc;                                                 //  1.
  int opt_max_flip_level;                                          // '-O', 3.
  int opt_scheme;                                                // '-O/#', 7.
  int opt_iterations;                                             // -O//#, 3.
  int smooth_cirterion;                                              // -s, 1.
  int smooth_maxiter;                                                // -s, 7.
  int delmaxfliplevel;                                                   // 1.
  int order;                                                       // '-o', 1.
  int reversetetori;                                              // '-o/', 0.
  int steinerleft;                                                 // '-S', 0.
  int unflip_queue_limit;                                      // '-U#', 1000.
  int no_sort;                                                           // 0.
  int hilbert_order;                                           // '-b///', 52.
  int hilbert_limit;                                             // '-b//'  8.
  int brio_threshold;                                              // '-b' 64.
  double brio_ratio;                                             // '-b/' 0.125.
  double epsilon;                                               // '-T', 1.0e-8.
  double facet_separate_ang_tol;                                 // '-p', 179.9.
  double collinear_ang_tol;                                     // '-p/', 179.9.
  double facet_small_ang_tol;                                   // '-p//', 15.0.
  double maxvolume;                                               // '-a', -1.0.
  double maxvolume_length;                                        // '-a', -1.0.
  double minratio;                                                 // '-q', 0.0.
  double opt_max_asp_ratio;                                           // 1000.0.
  double opt_max_edge_ratio;                                           // 100.0.
  double mindihedral;                                              // '-q', 5.0.
  double optmaxdihedral;                                          // -o/# 177.0.
  double metric_scale;                                              // -m#, 1.0.
  double smooth_alpha;                                             // '-s', 0.3.
  double coarsen_percent;                                         // -R1/#, 1.0.
  double elem_growth_ratio;             // Growth ratio of # elements, -r#, 0.0.
  double refine_progress_ratio;                                   // -r/#, 0.333.

  // Strings of command line arguments and input/output file names.
  char commandline[1024];
  char infilename[1024];
  char outfilename[1024];
  char addinfilename[1024];
  char bgmeshfilename[1024];

  // Read an additional tetrahedral mesh and treat it as holes [2018-07-30].
  int hole_mesh;                                                   // '-H', 0.
  char hole_mesh_filename[1024];

  // The input object of the tet core. They are recognized by either the input
  //   file extensions or by the specified options.
  // Currently the following objects are supported:
  //   - NODES, a list of nodes (.node);
  //   - POLY, a piecewise linear complex (.poly or .smesh);
  //   - OFF, a polyhedron (.off, Geomview's file format);
  //   - PLY, a polyhedron (.ply, file format from gatech, only ASCII);
  //   - STL, a surface mesh (.stl, stereolithography format);
  //   - MEDIT, a surface mesh (.mesh, Medit's file format);
  //   - MESH, a tetrahedral mesh (.ele).
  // If no extension is available, the imposed command line switch
  //   (-p or -r) implies the object.
  enum objecttype {NODES, POLY, OFF, PLY, STL, MEDIT, VTK, MESH, NEU_MESH} object;

  // Initialize all variables.
  TetMeshBehavior()
  {
    plc = 0;
    psc = 0;
    refine = 0;
    quality = 0;
    nobisect = 0;
    cdt = 0; // set by -D (without a number following it)
    cdtrefine = 7; // default, set by -D#
    coarsen = 0;
    metric = 0;
    weighted = 0;
    brio_hilbert = 1;
    flipinsert = 0;
    varvolume = 0;
    fixedvolume = 0;
    noexact = 0;
    nostaticfilter = 0;
    insertaddpoints = 0;
    regionattrib = 0;
    diagnose = 0;
    convex = 0;
    zeroindex = 0;
    facesout = 0;
    edgesout = 0;
    neighout = 0;
    voroout = 0;
    meditview = 0;
    vtkview = 0;
    vtksurfview = 0;
    nobound = 0;
    nonodewritten = 0;
    noelewritten = 0;
    nofacewritten = 0;
    noiterationnum = 0;
    nomergefacet = 0;
    nomergevertex = 0;
    nojettison = 0;
    docheck = 0;
    quiet = 0;
    nowarning = 0;
    verbose = 0;

    vertexperblock = 4092;
    tetrahedraperblock = 8188;
    shellfaceperblock = 4092;
    supsteiner_level = 2;
    addsteiner_algo = 1;
    coarsen_param = 0;
    weighted_param = 0;
    fliplinklevel = -1;
    flipstarsize = -1;
    fliplinklevelinc = 1;
    opt_scheme = 7;
    opt_max_flip_level = 3;
    opt_iterations = 3;
    delmaxfliplevel = 1;
    order = 1;
    reversetetori = 0;
    steinerleft = -1;
    unflip_queue_limit = 1000;
    no_sort = 0;
    hilbert_order = 52; //-1;
    hilbert_limit = 8;
    brio_threshold = 64;
    brio_ratio = 0.125;
    facet_separate_ang_tol = 179.9;
    collinear_ang_tol = 179.9;
    facet_small_ang_tol = 15.0;
    maxvolume = -1.0;
    maxvolume_length = -1.0;
    minratio = 2.0;
    opt_max_asp_ratio = 1000.;
    opt_max_edge_ratio = 100.;
    mindihedral = 3.5;
    optmaxdihedral = 177.00;
    epsilon = 1.0e-8;
    coarsen_percent = 1.0;
    metric_scale = 1.0; // -m#
    elem_growth_ratio = 0.0; // -r#
    refine_progress_ratio = 0.333; // -r/#
    object = NODES;

    smooth_cirterion = 3; // -s# default smooth surface and volume vertices.
    smooth_maxiter = 7;   // set by -s#/7
    smooth_alpha = 0.3;   // relax parameter, set by -s#/#/0.3

    commandline[0] = '\0';
    infilename[0] = '\0';
    outfilename[0] = '\0';
    addinfilename[0] = '\0';
    bgmeshfilename[0] = '\0';

    hole_mesh = 0;
    hole_mesh_filename[0] = '\0';

  }

}; // class TetMeshBehavior

} // namespace sqmesh::mesh::tet::detail
