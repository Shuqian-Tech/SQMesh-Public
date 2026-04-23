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
#include "tet_core.hpp"

namespace sqmesh::mesh::tet::detail {

//== main_cxx ================================================================//
//                                                                            //
//                                                                            //

//============================================================================//
//                                                                            //
// run_tet_mesh_core()    The interface for users using tet core library to        //
//                     generate tetrahedral meshes with all features.         //
//                                                                            //
// The sequence is roughly as follows.  Many of these steps can be skipped,   //
// depending on the command line switches.                                    //
//                                                                            //
// - Initialize constants and parse the command line.                         //
// - Read the vertices from a file and either                                 //
//   - run_tet_mesh_core them (no -r), or                                        //
//   - read an old mesh from files and reconstruct it (-r).                   //
// - Insert the boundary segments and facets (-p or -Y).                      //
// - Read the holes (-p), regional attributes (-pA), and regional volume      //
//   constraints (-pa).  Carve the holes and concavities, and spread the      //
//   regional attributes and volume constraints.                              //
// - Enforce the constraints on minimum quality bound (-q) and maximum        //
//   volume (-a), and a mesh size function (-m).                              //
// - Optimize the mesh wrt. specified quality measures (-O and -o).           //
// - Write the output files and print the statistics.                         //
// - Check the consistency of the mesh (-C).                                  //
//                                                                            //
//============================================================================//

void run_tet_mesh_core(TetMeshBehavior *b, TetMeshData *in, TetMeshData *out,
                    TetMeshData *addin, TetMeshData *bgmin)
{
  TetMeshCore m;
  clock_t tv[13], ts[6]; // Timing informations (defined in time.h)
  double cps = (double) CLOCKS_PER_SEC;

  tv[0] = clock();
 
  m.b = b;
  m.in = in;
  m.addin = addin;

  if (b->metric && bgmin && (bgmin->numberofpoints > 0)) {
    m.bgm = new TetMeshCore(); // Create an empty background mesh.
    m.bgm->b = b;
    m.bgm->in = bgmin;
  }

  m.initializepools();
  m.transfernodes();


  tv[1] = clock();

  if (b->refine) { // -r
    m.reconstructmesh();
  } else { // -p
    m.incrementaldelaunay(ts[0]);
  }

  tv[2] = clock();

  if (!b->quiet) {
    if (b->refine) {
      printf("Mesh reconstruction seconds:  %g\n", ((double)(tv[2]-tv[1])) / cps);
    } else {
      printf("Delaunay seconds:  %g\n", ((double)(tv[2]-tv[1])) / cps);
      if (b->verbose) {
        printf("  Point sorting seconds:  %g\n", ((double)(ts[0]-tv[1])) / cps);

      }
    }
  }

  if (b->plc && !b->refine) { // -p
    m.meshsurface();

    ts[0] = clock();

    if (!b->quiet) {
      printf("Surface mesh seconds:  %g\n", ((double)(ts[0]-tv[2])) / cps);
    }
  }


  tv[3] = clock();

  if ((b->metric) && (m.bgm != NULL)) { // -m
    m.bgm->initializepools();
    m.bgm->transfernodes();
    m.bgm->reconstructmesh();

    ts[0] = clock();

    if (!b->quiet) {
      printf("Background mesh reconstruct seconds:  %g\n",
             ((double)(ts[0] - tv[3])) / cps);
    }

    if (b->metric) { // -m
      m.interpolatemeshsize();

      ts[1] = clock();

      if (!b->quiet) {
        printf("Size interpolating seconds:  %g\n",((double)(ts[1]-ts[0])) / cps);
      }
    }
  }

  tv[4] = clock();

  if (b->plc && !b->refine) { // -p
    if (!b->cdt) { // no -D
      m.recoverboundary(ts[0]);
    } else {
      m.constraineddelaunay(ts[0]);
    }

    ts[1] = clock();

    if (!b->quiet) {
      if (!b->cdt) { // no -D
        printf("Boundary recovery ");
      } else {
        printf("Constrained Delaunay ");
      }
      printf("seconds:  %g\n", ((double)(ts[1] - tv[4])) / cps);
      if (b->verbose) {
        printf("  Segment recovery seconds:  %g\n",((double)(ts[0]-tv[4]))/ cps);
        printf("  Facet recovery seconds:  %g\n", ((double)(ts[1]-ts[0])) / cps);
      }
    }

    if (m.skipped_facet_list != NULL) {
      if (!b->quiet) {
        printf("\n!!! %ld input triangles are skipped due to self-intersections.\n",
               m.skipped_facet_list->objects);
      }

      delete m.skipped_facet_list;
      m.skipped_facet_list = NULL;

      if (!b->nonodewritten) m.outnodes(out);
      if (!b->noelewritten)  m.outelements(out);
      if (!b->nofacewritten) m.outsubfaces(out);
      if (!b->nofacewritten) m.outsubsegments(out);

      terminate_tet_core(NULL, 3); // This is not a normal exit.
    }

    if (b->diagnose) { // -d
      if (!b->quiet) {
        printf("\nThe input surface mesh is correct.\n");
      }
      return;
    }

    m.carveholes();

    ts[2] = clock();

    if (!b->quiet) {
      printf("Exterior tets removal seconds:  %g\n",((double)(ts[2]-ts[1]))/cps);
    }

    ts[3] = clock();

    if ((!b->cdt || b->nobisect) && (b->supsteiner_level > 0)) { // no -D, -Y/1
      if (m.subvertstack->objects > 0l) {
        m.suppresssteinerpoints();
        if (!b->quiet) {
          printf("Steiner suppression seconds:  %g\n", ((double)(ts[3]-ts[2]))/cps);
        }
      }
    }
    
    if ((b->nobisect > 1)) { // -YY
      if ((m.st_segref_count > 0) || (m.st_facref_count > 0)) {
        if (!b->nonodewritten) m.outnodes(out);
        if (!b->noelewritten)  m.outelements(out);
        if (!b->nofacewritten) m.outsubfaces(out);
        if (!b->nofacewritten) m.outsubsegments(out);
        printf("!! Boundary contains Steiner points (-YY option). Program stopped.\n");
        terminate_tet_core(&m, 200);
      }
    }
  }

  tv[5] = clock();

  if (b->metric || b->coarsen) { // -m or -R
    m.meshcoarsening();
  }

  tv[6] = clock();

  if (!b->quiet) {
    if (b->metric || b->coarsen) {
      printf("Mesh coarsening seconds:  %g\n", ((double)(tv[6] - tv[5])) / cps);
    }
  }

  if (b->plc || (b->refine && b->quality && (in->refine_elem_list == NULL))) {
    if (!b->quiet) {
      printf("Recovering Delaunayness...\n");
    }
    m.recoverdelaunay();
  }

  tv[7] = clock();

  if (b->plc || (b->refine && b->quality && (in->refine_elem_list == NULL))) {
    if (!b->quiet) {
      printf("Delaunay recovery seconds:  %g\n", ((double)(tv[7] - tv[6]))/cps);
    }
  }

  if ((b->plc || b->refine) && b->insertaddpoints) { // -i
    if ((addin != NULL) && (addin->numberofpoints > 0)) {
      m.insertconstrainedpoints(addin); 
    }
  }

  tv[8] = clock();

  if (!b->quiet) {
    if ((b->plc || b->refine) && b->insertaddpoints) { // -i
      if ((addin != NULL) && (addin->numberofpoints > 0)) {
        printf("Constrained points seconds:  %g\n", ((double)(tv[8]-tv[7]))/cps);
      }
    }
  }
  if (b->quality) { // -q
    m.delaunayrefinement();    
  }

  tv[9] = clock();

  if (!b->quiet) {
    if (b->quality) {
      printf("Refinement seconds:  %g\n", ((double)(tv[9] - tv[8])) / cps);
    }
  }

  if ((b->plc || b->quality) &&
      (b->smooth_maxiter > 0) &&
      ((m.st_volref_count > 0) || (m.st_facref_count > 0))) {
    m.smooth_vertices(); // m.optimizemesh(ts[0]);
  }

  tv[10] = clock();

  if (!b->quiet) {
    if ((b->plc || b->quality) &&
        (b->smooth_maxiter > 0) &&
        ((m.st_volref_count > 0) || (m.st_facref_count > 0))) {
      printf("Mesh smoothing seconds:  %g\n", ((double)(tv[10] - tv[9])) / cps);
    }
  }

  if (b->plc || b->quality) {
    m.improve_mesh();
  }

  tv[11] = clock();

  if (!b->quiet) {
    if (b->plc || b->quality) {
      printf("Mesh improvement seconds:  %g\n", ((double)(tv[11] - tv[10])) / cps);
    }
  }

  if (!b->nojettison && ((m.dupverts > 0) || (m.unuverts > 0)
      || (b->refine && (in->numberofcorners == 10)))) {
    m.jettisonnodes();
  }


  if ((b->order == 2) && !b->convex) {
    m.highorder();
  }

  if (!b->quiet) {
    printf("\n");
  }

  if (out != (TetMeshData *) NULL) {
    out->firstnumber = in->firstnumber;
    out->mesh_dim = in->mesh_dim;
  }

  if (b->nonodewritten || b->noiterationnum) {
    if (!b->quiet) {
      printf("NOT writing a .node file.\n");
    }
  } else {
    m.outnodes(out);
  }

  if (b->noelewritten) {
    if (!b->quiet) {
      printf("NOT writing an .ele file.\n");
    }
    m.indexelements();
  } else {
    if (m.tetrahedrons->items > 0l) {
      m.outelements(out);
    }
  }

  if (b->nofacewritten) {
    if (!b->quiet) {
      printf("NOT writing an .face file.\n");
    }
  } else {
    if (b->facesout) {
      if (m.tetrahedrons->items > 0l) {
        m.outfaces(out);  // Output all faces.
      }
    } else {
      if (b->plc || b->refine) {
        if (m.subfaces->items > 0l) {
          m.outsubfaces(out); // Output boundary faces.
        }
      } else {
        if (m.tetrahedrons->items > 0l) {
          m.outhullfaces(out); // Output convex hull faces.
        }
      }
    }
  }


  if (b->nofacewritten) {
    if (!b->quiet) {
      printf("NOT writing an .edge file.\n");
    }
  } else {
    if (b->edgesout) { // -e
      m.outedges(out); // output all mesh edges. 
    } else {
      if (b->plc || b->refine) {
        m.outsubsegments(out); // output subsegments.
      }
    }
  }

  if ((b->plc || b->refine) && b->metric) { // -m
    m.outmetrics(out);
  }

  if (b->neighout) {
    m.outneighbors(out);
  }

  if (b->voroout) {
    m.outvoronoi(out);
  }


  tv[12] = clock();

  if (!b->quiet) {
    printf("\nOutput seconds:  %g\n", ((double)(tv[12] - tv[11])) / cps);
    printf("Total running seconds:  %g\n", ((double)(tv[12] - tv[0])) / cps);
  }

  if (b->docheck) {
    m.check_mesh(0);
    if (b->plc || b->refine) {
      m.check_shells();
      m.check_segments();
    }
    if (b->docheck > 1) {
      m.check_delaunay();
    }
  }

  if (!b->quiet) {
    m.statistics();
  }
}


} // namespace sqmesh::mesh::tet::detail
