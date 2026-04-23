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

//== output_cxx ==============================================================//
//                                                                            //
//                                                                            //

//============================================================================//
//                                                                            //
// jettisonnodes()    Jettison unused or duplicated vertices.                 //
//                                                                            //
// Unused points are those input points which are outside the mesh domain or  //
// have no connection (isolated) to the mesh.  Duplicated points exist for    //
// example if the input PLC is read from a .stl mesh file (marked during the  //
// Delaunay tetrahedralization step. This routine remove these points from    //
// points list. All existing points are reindexed.                            //
//                                                                            //
//============================================================================//

void TetMeshCore::jettisonnodes()
{
  point pointloop;
  bool jetflag;
  int oldidx, newidx;
  int remcount;

  if (!b->quiet) {
    printf("Jettisoning redundant points.\n");
  }

  points->traversalinit();
  pointloop = pointtraverse();
  oldidx = newidx = 0; // in->firstnumber;
  remcount = 0;
  while (pointloop != (point) NULL) {
    jetflag = (pointtype(pointloop) == DUPLICATEDVERTEX) || 
      (pointtype(pointloop) == UNUSEDVERTEX);
    if (jetflag) {
      // It is a duplicated or unused point, delete it.
      pointdealloc(pointloop);
      remcount++;
    } else {
      // Re-index it.
      setpointmark(pointloop, newidx + in->firstnumber);
      if (in->pointmarkerlist != (int *) NULL) {
        if (oldidx < in->numberofpoints) {
          // Re-index the point marker as well.
          in->pointmarkerlist[newidx] = in->pointmarkerlist[oldidx];
        }
      }
      newidx++;
    }
    oldidx++;
    pointloop = pointtraverse();
  }
  if (b->verbose) {
    printf("  %ld duplicated vertices are removed.\n", dupverts);
    printf("  %ld unused vertices are removed.\n", unuverts);
  }
  dupverts = 0l;
  unuverts = 0l;

  // The following line ensures that dead items in the pool of nodes cannot
  //   be allocated for the new created nodes. This ensures that the input
  //   nodes will occur earlier in the output files, and have lower indices.
  points->deaditemstack = (void *) NULL;
}

//============================================================================//
//                                                                            //
// highorder()   Create extra nodes for quadratic subparametric elements.     //
//                                                                            //
// 'highordertable' is an array (size = numberoftetrahedra * 6) for storing   //
// high-order nodes of each tetrahedron.  This routine is used only when -o2  //
// switch is used.                                                            //
//                                                                            //
//============================================================================//

void TetMeshCore::highorder()
{
  triface tetloop, worktet, spintet;
  point *extralist, *adjextralist;
  point torg, tdest, newpoint;
  int highorderindex;
  int t1ver;
  int i, j;

  if (!b->quiet) {
    printf("Adding vertices for second-order tetrahedra.\n");
  }

  // Initialize the 'highordertable'.
  point *highordertable = new point[tetrahedrons->items * 6];
  if (highordertable == (point *) NULL) {
    terminate_tet_core(this, 1);
  }

  // This will overwrite the slot for element markers.
  highorderindex = 11;

  // The following line ensures that dead items in the pool of nodes cannot
  //   be allocated for the extra nodes associated with high order elements.
  //   This ensures that the primary nodes (at the corners of elements) will
  //   occur earlier in the output files, and have lower indices, than the
  //   extra nodes.
  points->deaditemstack = (void *) NULL;

  // Assign an entry for each tetrahedron to find its extra nodes. At the
  //   mean while, initialize all extra nodes be NULL.
  i = 0;
  tetrahedrons->traversalinit();
  tetloop.tet = tetrahedrontraverse();
  while (tetloop.tet != (tetrahedron *) NULL) {
    tetloop.tet[highorderindex] = (tetrahedron) &highordertable[i];
    for (j = 0; j < 6; j++) {
      highordertable[i + j] = (point) NULL;
    }
    i += 6;
    tetloop.tet = tetrahedrontraverse();
  }

  // To create a unique node on each edge. Loop over all tetrahedra, and
  //   look at the six edges of each tetrahedron.  If the extra node in
  //   the tetrahedron corresponding to this edge is NULL, create a node
  //   for this edge, at the same time, set the new node into the extra
  //   node lists of all other tetrahedra sharing this edge.  
  tetrahedrons->traversalinit();
  tetloop.tet = tetrahedrontraverse();
  while (tetloop.tet != (tetrahedron *) NULL) {
    // Get the list of extra nodes.
    extralist = (point *) tetloop.tet[highorderindex];
    worktet.tet = tetloop.tet;
    for (i = 0; i < 6; i++) {
      if (extralist[i] == (point) NULL) {
        // Go to the ith-edge.
        worktet.ver = edge2ver[i];
        // Create a new point in the middle of this edge.
        torg = org(worktet);
        tdest = dest(worktet);
        makepoint(&newpoint, FREEVOLVERTEX);
        for (j = 0; j < 3 + numpointattrib; j++) {
          newpoint[j] = 0.5 * (torg[j] + tdest[j]);
        }
        // Interpolate its metrics.
        for (j = 0; j < in->numberofpointmtrs; j++) {
          newpoint[pointmtrindex + j] = 
            0.5 * (torg[pointmtrindex + j] + tdest[pointmtrindex + j]);
        }
        // Set this point into all extra node lists at this edge.
        spintet = worktet;
        while (1) {
          if (!ishulltet(spintet)) {
            adjextralist = (point *) spintet.tet[highorderindex];
            adjextralist[ver2edge[spintet.ver]] = newpoint;
          }
          fnextself(spintet);
          if (spintet.tet == worktet.tet) break;
        }
      } // if (!extralist[i])
    } // i
    tetloop.tet = tetrahedrontraverse();
  }

  delete [] highordertable;
}

//============================================================================//
//                                                                            //
// indexelements()    Index all tetrahedra.                                   //
//                                                                            //
// Many output functions require that the tetrahedra are indexed.  This       //
// routine is called when -E option is used.                                  //
//                                                                            //
//============================================================================//

void TetMeshCore::indexelements()
{
  triface worktet;
  int eindex = b->zeroindex ? 0 : in->firstnumber; // firstindex;
  tetrahedrons->traversalinit();
  worktet.tet = tetrahedrontraverse();
  while (worktet.tet != NULL) {
    setelemindex(worktet.tet, eindex);
    eindex++;
	if (b->metric) { // -m option
	  // Update the point-to-tet map, so that every point is pointing
	  //   to a real tet, not a fictious one. Used by .p2t file.
	  tetrahedron tptr = encode(worktet);
	  for (int i = 0; i < 4; i++) {
	    setpoint2tet((point) (worktet.tet[4 + i]), tptr);
	  }
	}
    worktet.tet = tetrahedrontraverse();
  }
}

//============================================================================//
//                                                                            //
// numberedges()    Count the number of edges, save in "meshedges".           //
//                                                                            //
// This routine is called when '-p' or '-r', and '-E' options are used.  The  //
// total number of edges depends on the genus of the input surface mesh.      //
//                                                                            //
// NOTE:  This routine must be called after outelements().  So all elements   //
// have been indexed.                                                         //
//                                                                            //
//============================================================================//

void TetMeshCore::numberedges()
{
  triface worktet, spintet;
  int ishulledge;
  int t1ver;
  int i;

  meshedges = meshhulledges = 0l;

  tetrahedrons->traversalinit();
  worktet.tet = tetrahedrontraverse();
  while (worktet.tet != NULL) {
    for (i = 0; i < 6; i++) {
      worktet.ver = edge2ver[i];
      ishulledge = 0;
      fnext(worktet, spintet);
      do {
        if (!ishulltet(spintet)) {
          if (elemindex(spintet.tet) < elemindex(worktet.tet)) break;
        } else {
          ishulledge = 1;
        }
        fnextself(spintet);
      } while (spintet.tet != worktet.tet);
      if (spintet.tet == worktet.tet) {
        meshedges++;
        if (ishulledge) meshhulledges++;
      }
    }
    infect(worktet);
    worktet.tet = tetrahedrontraverse();
  }
}

//============================================================================//
//                                                                            //
// outnodes()    Output the points to a .node file or a TetMeshData structure.   //
//                                                                            //
// Note: each point has already been numbered on input (the first index is    //
// 'in->firstnumber').                                                        //
//                                                                            //
//============================================================================//

void TetMeshCore::outnodes(TetMeshData* out)
{
  FILE *outfile = NULL;
  char outnodefilename[FILENAMESIZE];
  face parentsh;
  point pointloop;
  int nextras, bmark, marker = 0, weightDT = 0; 
  int coordindex, attribindex;
  int pointnumber, firstindex;
  int index, i;

  if (out == (TetMeshData *) NULL) {
    strcpy(outnodefilename, b->outfilename);
    strcat(outnodefilename, ".node");
  }

  if (!b->quiet) {
    if (out == (TetMeshData *) NULL) {
      printf("Writing %s.\n", outnodefilename);
    } else {
      printf("Writing nodes.\n");
    }
  }

  nextras = numpointattrib;
  if (b->weighted) { // -w
    if (b->weighted_param == 0) weightDT = 1; // Weighted DT.
  }

  bmark = !b->nobound && in->pointmarkerlist;

  if (out == (TetMeshData *) NULL) {
    outfile = fopen(outnodefilename, "w");
    if (outfile == (FILE *) NULL) {
      printf("File I/O Error:  Cannot create file %s.\n", outnodefilename);
      terminate_tet_core(this, 1);
    }
    // Number of points, number of dimensions, number of point attributes,
    //   and number of boundary markers (zero or one).
    //fprintf(outfile, "%ld  %d  %d  %d\n", points->items, 3, nextras, bmark);
    // [2020-01-16] added save flag (for viewing Steiner points).
    //fprintf(outfile, "%ld  %d  %d  %d  1\n", points->items, 3, nextras, bmark);
    fprintf(outfile, "%ld  %d  %d  %d\n", points->items, 3, nextras, bmark);
  } else {
    // Allocate space for 'pointlist';
    out->pointlist = new double[points->items * 3];
    if (out->pointlist == (double *) NULL) {
      printf("Error:  Out of memory.\n");
      terminate_tet_core(this, 1);
    }
    // Allocate space for 'pointattributelist' if necessary;
    if (nextras > 0) {
      out->pointattributelist = new double[points->items * nextras];
      if (out->pointattributelist == (double *) NULL) {
        printf("Error:  Out of memory.\n");
        terminate_tet_core(this, 1);
      }
    }
    // Allocate space for 'pointmarkerlist' if necessary;
    if (bmark) {
      out->pointmarkerlist = new int[points->items];
      if (out->pointmarkerlist == (int *) NULL) {
        printf("Error:  Out of memory.\n");
        terminate_tet_core(this, 1);
      }
    }
    if (b->psc) {
      out->pointparamlist = new TetMeshData::pointparam[points->items];
      if (out->pointparamlist == NULL) {
        printf("Error:  Out of memory.\n");
        terminate_tet_core(this, 1);
      }
    }
    out->numberofpoints = points->items;
    out->numberofpointattributes = nextras;
    coordindex = 0;
    attribindex = 0;
  }
  
  // Determine the first index (0 or 1).
  firstindex = b->zeroindex ? 0 : in->firstnumber;

  points->traversalinit();
  pointloop = pointtraverse();
  pointnumber = firstindex; // in->firstnumber;
  index = 0;
  while (pointloop != (point) NULL) {
    if (bmark) {
      // Default the vertex has a zero marker.
      marker = 0;
      // Is it an input vertex?
      if (index < in->numberofpoints) {
        // Input point's marker is directly copied to output.
        marker = in->pointmarkerlist[index];       
      } else {
        if ((pointtype(pointloop) == FREESEGVERTEX) ||
            (pointtype(pointloop) == FREEFACETVERTEX)) {
          sdecode(point2sh(pointloop), parentsh);
          if (parentsh.sh != NULL) {
            marker = shellmark(parentsh);
          }
        } // if (pointtype(...))
      }
    }
    if (out == (TetMeshData *) NULL) {
      // Point number, x, y and z coordinates.
      fprintf(outfile, "%4d    %.17g  %.17g  %.17g", pointnumber,
              pointloop[0], pointloop[1], pointloop[2]);
      for (i = 0; i < nextras; i++) {
        // Write an attribute.
        if ((i == 0) && weightDT) {          
          fprintf(outfile, "  %.17g", pointloop[0] * pointloop[0] +
             pointloop[1] * pointloop[1] + pointloop[2] * pointloop[2] 
             - pointloop[3 + i]);
        } else { 
          fprintf(outfile, "  %.17g", pointloop[3 + i]);
        }
      }
      if (bmark) {
        // Write the boundary marker.
        fprintf(outfile, "    %d", marker);
      }
      if (b->psc) {
        fprintf(outfile, "  %.8g  %.8g  %d", pointgeomuv(pointloop, 0),
                pointgeomuv(pointloop, 1), pointgeomtag(pointloop));
        if (pointtype(pointloop) == RIDGEVERTEX) {
          fprintf(outfile, "  0");
        //} else if (pointtype(pointloop) == ACUTEVERTEX) {
        //  fprintf(outfile, "  0");
        } else if (pointtype(pointloop) == FREESEGVERTEX) {
          fprintf(outfile, "  1");
        } else if (pointtype(pointloop) == FREEFACETVERTEX) {
          fprintf(outfile, "  2");
        } else if (pointtype(pointloop) == FREEVOLVERTEX) {
          fprintf(outfile, "  3");
        } else {
          fprintf(outfile, "  -1"); // Unknown type.
        }
      }
      // // [2020-01-16] Write vertex flags
      // if (pointnumber > in->numberofpoints) {
      //   fprintf(outfile, " 16"); // A Steiner point.
      // } else {
      //   fprintf(outfile, " 0");
      // }
      fprintf(outfile, "\n");
    } else {
      // X, y, and z coordinates.
      out->pointlist[coordindex++] = pointloop[0];
      out->pointlist[coordindex++] = pointloop[1];
      out->pointlist[coordindex++] = pointloop[2];
      // Point attributes.
      for (i = 0; i < nextras; i++) {
        // Output an attribute.
        if ((i == 0) && weightDT) {
          out->pointattributelist[attribindex++] = 
            pointloop[0] * pointloop[0] + pointloop[1] * pointloop[1] + 
            pointloop[2] * pointloop[2] - pointloop[3 + i];
        } else {
          out->pointattributelist[attribindex++] = pointloop[3 + i];
        }
      }
      if (bmark) {
        // Output the boundary marker.  
        out->pointmarkerlist[index] = marker;
      }
      if (b->psc) {
        out->pointparamlist[index].uv[0] = pointgeomuv(pointloop, 0);
        out->pointparamlist[index].uv[1] = pointgeomuv(pointloop, 1);
        out->pointparamlist[index].tag = pointgeomtag(pointloop);
        if (pointtype(pointloop) == RIDGEVERTEX) {
          out->pointparamlist[index].type = 0;
        //} else if (pointtype(pointloop) == ACUTEVERTEX) {
        //  out->pointparamlist[index].type = 0;
        } else if (pointtype(pointloop) == FREESEGVERTEX) {
          out->pointparamlist[index].type = 1;
        } else if (pointtype(pointloop) == FREEFACETVERTEX) {
          out->pointparamlist[index].type = 2;
        } else if (pointtype(pointloop) == FREEVOLVERTEX) {
          out->pointparamlist[index].type = 3;
        } else {
          out->pointparamlist[index].type = -1; // Unknown type.
        }
      }
    }
    pointloop = pointtraverse();
    pointnumber++; 
    index++;
  }

  if (out == (TetMeshData *) NULL) {
    fprintf(outfile, "# Generated by %s\n", b->commandline);
    fclose(outfile);
  }
}

//============================================================================//
//                                                                            //
// outmetrics()    Output the metric to a file (*.mtr) or a TetMeshData obj.     //
//                                                                            //
//============================================================================//

void TetMeshCore::outmetrics(TetMeshData* out)
{
  FILE *outfile = NULL;
  char outmtrfilename[FILENAMESIZE];
  point ptloop;
  int mtrindex = 0;
  int i;
  int msize = (sizeoftensor - useinsertradius);
  if (msize == 0) {
    return;
  }

  if (out == (TetMeshData *) NULL) {
    strcpy(outmtrfilename, b->outfilename);
    strcat(outmtrfilename, ".mtr");
  }

  if (!b->quiet) {
    if (out == (TetMeshData *) NULL) {
      printf("Writing %s.\n", outmtrfilename);
    } else {
      printf("Writing metrics.\n");
    }
  }

  if (out == (TetMeshData *) NULL) {
    outfile = fopen(outmtrfilename, "w");
    if (outfile == (FILE *) NULL) {
      printf("File I/O Error:  Cannot create file %s.\n", outmtrfilename);
      terminate_tet_core(this, 3);
    }
    // Number of points, number of point metrices,
    fprintf(outfile, "%ld  %d\n", points->items, msize);
  } else {
    // Allocate space for 'pointmtrlist'.
    out->numberofpointmtrs = msize; 
    out->pointmtrlist = new double[points->items * msize];
    if (out->pointmtrlist == (double *) NULL) {
      terminate_tet_core(this, 1);
    }
  }

  points->traversalinit();
  ptloop = pointtraverse();
  while (ptloop != (point) NULL) {
    if (out == (TetMeshData *) NULL) {
      for (i = 0; i < msize; i++) {
        fprintf(outfile, " %-16.8e", ptloop[pointmtrindex + i]);
      }
      fprintf(outfile, "\n");
    } else {
      for (i = 0; i < msize; i++) {
        out->pointmtrlist[mtrindex++] = ptloop[pointmtrindex + i];
      }
    }
    ptloop = pointtraverse();
  }

  // Output the point-to-tet map.
  if (out == (TetMeshData *) NULL) {
    strcpy(outmtrfilename, b->outfilename);
    strcat(outmtrfilename, ".p2t");
  }

  if (!b->quiet) {
    if (out == (TetMeshData *) NULL) {
      printf("Writing %s.\n", outmtrfilename);
    } else {
      printf("Writing point-to-tet map.\n");
    }
  }

  if (out == (TetMeshData *) NULL) {
    outfile = fopen(outmtrfilename, "w");
    if (outfile == (FILE *) NULL) {
      printf("File I/O Error:  Cannot create file %s.\n", outmtrfilename);
      terminate_tet_core(this, 3);
    }
    // Number of points,
    //fprintf(outfile, "%ld\n", points->items);
  } else {
    // Allocate space for 'point2tetlist'.
    out->point2tetlist = new int[points->items];
    if (out->point2tetlist == (int *) NULL) {
      terminate_tet_core(this, 1);
    }
  }

  // The list of tetrahedra must be indexed.
  if (bgm != NULL) {
    bgm->indexelements();
  }
  // Determine the first index (0 or 1).
  int firstindex = b->zeroindex ? 0 : in->firstnumber;
  int pointindex = firstindex;
  i = 0;

  triface parenttet;
  points->traversalinit();
  ptloop = pointtraverse();
  while (ptloop != (point) NULL) {
    if (bgm != NULL) {
	  bgm->decode(point2bgmtet(ptloop), parenttet);
	} else {
	  decode(point2tet(ptloop), parenttet);
	}
    if (out == (TetMeshData *) NULL) {
      fprintf(outfile, "%d  %d\n", pointindex, elemindex(parenttet.tet));
    } else {
      out->point2tetlist[i] = elemindex(parenttet.tet);
    }
	pointindex++;
	i++;
    ptloop = pointtraverse();
  }

  if (out == (TetMeshData *) NULL) {
    fprintf(outfile, "# Generated by %s\n", b->commandline);
    fclose(outfile);
  }
}

//============================================================================//
//                                                                            //
// outelements()    Output the tetrahedra to an .ele file or a TetMeshData       //
//                  structure.                                                //
//                                                                            //
// This routine also indexes all tetrahedra (exclusing hull tets) (from in->  //
// firstnumber). The total number of mesh edges is counted in 'meshedges'.    //
//                                                                            //
//============================================================================//

void TetMeshCore::outelements(TetMeshData* out)
{
  FILE *outfile = NULL;
  char outelefilename[FILENAMESIZE];
  tetrahedron* tptr;
  point p1, p2, p3, p4;
  point *extralist;
  double *talist = NULL;
  int *tlist = NULL;
  long ntets;
  int firstindex, shift;
  int pointindex, attribindex;
  int highorderindex = 11; 
  int elementnumber;
  int eextras;
  int i;

  if (out == (TetMeshData *) NULL) {
    strcpy(outelefilename, b->outfilename);
    strcat(outelefilename, ".ele");
  }

  if (!b->quiet) {
    if (out == (TetMeshData *) NULL) {
      printf("Writing %s.\n", outelefilename);
    } else {
      printf("Writing elements.\n");
    }
  }

  // The number of tets excluding hull tets.
  ntets = tetrahedrons->items - hullsize;

  eextras = numelemattrib;
  if (out == (TetMeshData *) NULL) {
    outfile = fopen(outelefilename, "w");
    if (outfile == (FILE *) NULL) {
      printf("File I/O Error:  Cannot create file %s.\n", outelefilename);
      terminate_tet_core(this, 1);
    }
    // Number of tetras, points per tetra, attributes per tetra.
    fprintf(outfile, "%ld  %d  %d\n", ntets, b->order == 1 ? 4 : 10, eextras);
  } else {
    // Allocate memory for output tetrahedra.
    out->tetrahedronlist = new int[ntets * (b->order == 1 ? 4 : 10)];
    if (out->tetrahedronlist == (int *) NULL) {
      printf("Error:  Out of memory.\n");
      terminate_tet_core(this, 1);
    }
    // Allocate memory for output tetrahedron attributes if necessary.
    if (eextras > 0) {
      out->tetrahedronattributelist = new double[ntets * eextras];
      if (out->tetrahedronattributelist == (double *) NULL) {
        printf("Error:  Out of memory.\n");
        terminate_tet_core(this, 1);
      }
    }
    out->numberoftetrahedra = ntets;
    out->numberofcorners = b->order == 1 ? 4 : 10;
    out->numberoftetrahedronattributes = eextras;
    tlist = out->tetrahedronlist;
    talist = out->tetrahedronattributelist;
    pointindex = 0;
    attribindex = 0;
  }

  // Determine the first index (0 or 1).
  firstindex = b->zeroindex ? 0 : in->firstnumber;
  shift = 0; // Default no shift.
  if ((in->firstnumber == 1) && (firstindex == 0)) {
    shift = 1; // Shift the output indices by 1.
  }

  tetrahedrons->traversalinit();
  tptr = tetrahedrontraverse();
  elementnumber = firstindex; // in->firstnumber;
  while (tptr != (tetrahedron *) NULL) {
    if (!b->reversetetori) {
      p1 = (point) tptr[4];
      p2 = (point) tptr[5];
    } else {
      p1 = (point) tptr[5];
      p2 = (point) tptr[4];
    }
    p3 = (point) tptr[6];
    p4 = (point) tptr[7];
    if (out == (TetMeshData *) NULL) {
      // Tetrahedron number, indices for four points.
      fprintf(outfile, "%5d   %5d %5d %5d %5d", elementnumber,
              pointmark(p1) - shift, pointmark(p2) - shift,
              pointmark(p3) - shift, pointmark(p4) - shift);
      if (b->order == 2) {
        extralist = (point *) tptr[highorderindex];
        // indices for six extra points.
        fprintf(outfile, "  %5d %5d %5d %5d %5d %5d",
          pointmark(extralist[0]) - shift, pointmark(extralist[1]) - shift,
          pointmark(extralist[2]) - shift, pointmark(extralist[3]) - shift,
          pointmark(extralist[4]) - shift, pointmark(extralist[5]) - shift);
      }
      for (i = 0; i < eextras; i++) {
        fprintf(outfile, "    %.17g", elemattribute(tptr, i));
      }
      fprintf(outfile, "\n");
    } else {
      tlist[pointindex++] = pointmark(p1) - shift;
      tlist[pointindex++] = pointmark(p2) - shift;
      tlist[pointindex++] = pointmark(p3) - shift;
      tlist[pointindex++] = pointmark(p4) - shift;
      if (b->order == 2) {
        extralist = (point *) tptr[highorderindex];
        tlist[pointindex++] = pointmark(extralist[0]) - shift;
        tlist[pointindex++] = pointmark(extralist[1]) - shift;
        tlist[pointindex++] = pointmark(extralist[2]) - shift;
        tlist[pointindex++] = pointmark(extralist[3]) - shift;
        tlist[pointindex++] = pointmark(extralist[4]) - shift;
        tlist[pointindex++] = pointmark(extralist[5]) - shift;
      }
      for (i = 0; i < eextras; i++) {
        talist[attribindex++] = elemattribute(tptr, i);
      }
    }
    // Remember the index of this element (for counting edges).
    setelemindex(tptr, elementnumber);
	if (b->metric) { // -m option
	  // Update the point-to-tet map, so that every point is pointing
	  //   to a real tet, not a fictious one. Used by .p2t file.
	  for (int i = 0; i < 4; i++) {
	    setpoint2tet((point) (tptr[4 + i]), (tetrahedron) tptr);
	  }
	}
    tptr = tetrahedrontraverse();
    elementnumber++;
  }


  if (out == (TetMeshData *) NULL) {
    fprintf(outfile, "# Generated by %s\n", b->commandline);
    fclose(outfile);
  }
}

//============================================================================//
//                                                                            //
// outfaces()    Output all faces to a .face file or a TetMeshData object.       //
//                                                                            //
// The total number of faces f can be calculated as following:  Let t be the  //
// total number of tets. Since each tet has 4 faces, the number t * 4 counts  //
// each interior face twice and each hull face once. So f = (t * 4 + h) / 2,  //
// where h is the total number of hull faces (which is known).                //
//                                                                            //
//============================================================================//

void TetMeshCore::outfaces(TetMeshData* out)
{
  FILE *outfile = NULL;
  char facefilename[FILENAMESIZE];
  triface tface, tsymface;
  face checkmark;
  point torg, tdest, tapex;
  long ntets, faces;
  int *elist = NULL, *emlist = NULL;
  int neigh1 = 0, neigh2 = 0;
  int marker = 0;
  int firstindex, shift;
  int facenumber;
  int index = 0;

  // For -o2 option.
  triface workface;
  point *extralist, pp[3] = {0,0,0}; 
  int highorderindex = 11; 
  int o2index = 0, i;

  // For -nn option.
  int *tet2facelist = NULL;
  int tidx; 

  if (out == (TetMeshData *) NULL) {
    strcpy(facefilename, b->outfilename);
    strcat(facefilename, ".face");
  }

  if (!b->quiet) {
    if (out == (TetMeshData *) NULL) {
      printf("Writing %s.\n", facefilename);
    } else {
      printf("Writing faces.\n");
    }
  }

  ntets = tetrahedrons->items - hullsize;
  faces = (ntets * 4l + hullsize) / 2l;

  if (out == (TetMeshData *) NULL) {
    outfile = fopen(facefilename, "w");
    if (outfile == (FILE *) NULL) {
      printf("File I/O Error:  Cannot create file %s.\n", facefilename);
      terminate_tet_core(this, 1);
    }
    fprintf(outfile, "%ld  %d\n", faces, !b->nobound);
  } else {
    // Allocate memory for 'trifacelist'.
    out->trifacelist = new int[faces * 3];
    if (out->trifacelist == (int *) NULL) {
      printf("Error:  Out of memory.\n");
      terminate_tet_core(this, 1);
    }
    if (b->order == 2) {
      out->o2facelist = new int[faces * 3];
    }
    // Allocate memory for 'trifacemarkerlist' if necessary.
    if (!b->nobound) {
      out->trifacemarkerlist = new int[faces];
      if (out->trifacemarkerlist == (int *) NULL) {
        printf("Error:  Out of memory.\n");
        terminate_tet_core(this, 1);
      }
    }
    if (b->neighout > 1) {
      // '-nn' switch.
      out->face2tetlist = new int[faces * 2];
      if (out->face2tetlist == (int *) NULL) {
        printf("Error:  Out of memory.\n");
        terminate_tet_core(this, 1);
      }
    }
    out->numberoftrifaces = faces;
    elist = out->trifacelist;
    emlist = out->trifacemarkerlist;
  }

  if (b->neighout > 1) { // -nn option
    // Output the tetrahedron-to-face map.
    tet2facelist = new int[ntets * 4];
  }

  // Determine the first index (0 or 1).
  firstindex = b->zeroindex ? 0 : in->firstnumber;
  shift = 0; // Default no shiftment.
  if ((in->firstnumber == 1) && (firstindex == 0)) {
    shift = 1; // Shift the output indices by 1.
  }

  tetrahedrons->traversalinit();
  tface.tet = tetrahedrontraverse();
  facenumber = firstindex; // in->firstnumber;
  // To loop over the set of faces, loop over all tetrahedra, and look at
  //   the four faces of each one. If its adjacent tet is a hull tet,
  //   operate on the face, otherwise, operate on the face only if the
  //   current tet has a smaller index than its neighbor.
  while (tface.tet != (tetrahedron *) NULL) {
    for (tface.ver = 0; tface.ver < 4; tface.ver ++) {
      fsym(tface, tsymface);
      if (ishulltet(tsymface) || 
          (elemindex(tface.tet) < elemindex(tsymface.tet))) {
        torg = org(tface);
        tdest = dest(tface);
        tapex = apex(tface);
        if (b->order == 2) { // -o2
          // Get the three extra vertices on edges.
          extralist = (point *) (tface.tet[highorderindex]);
          // The extra vertices are on edges opposite the corners.
          enext(tface, workface);
          for (i = 0; i < 3; i++) {
            pp[i] = extralist[ver2edge[workface.ver]];
            enextself(workface);
          }
        }
        if (!b->nobound) {
          // Get the boundary marker of this face.
          if (b->plc || b->refine) { 
            // Shell face is used.
            tspivot(tface, checkmark);
            if (checkmark.sh == NULL) {
              marker = 0;  // It is an inner face. It's marker is 0.
            } else {
              marker = shellmark(checkmark);
            }
          } else {
            // Shell face is not used, only distinguish outer and inner face.
            marker = (int) ishulltet(tsymface);
          }
        }
        if (b->neighout > 1) {
          // '-nn' switch. Output adjacent tets indices.
          if (!ishulltet(tface)) {
            neigh1 = elemindex(tface.tet);
          } else {
            neigh1 = -1;
          }
          if (!ishulltet(tsymface)) {
            neigh2 = elemindex(tsymface.tet);
          } else {
            neigh2 = -1;  
          }
		  // Fill the tetrahedron-to-face map.
		  tidx = elemindex(tface.tet) - firstindex;
		  tet2facelist[tidx * 4 + tface.ver] = facenumber;
		  if (!ishulltet(tsymface)) {
		    tidx = elemindex(tsymface.tet) - firstindex;
			tet2facelist[tidx * 4 + (tsymface.ver & 3)] = facenumber;
		  }
        }
        if (out == (TetMeshData *) NULL) {
          // Face number, indices of three vertices.
          fprintf(outfile, "%5d   %4d  %4d  %4d", facenumber,
                  pointmark(torg) - shift, pointmark(tdest) - shift,
                  pointmark(tapex) - shift);
          if (b->order == 2) { // -o2
            fprintf(outfile, "  %4d  %4d  %4d", pointmark(pp[0]) - shift, 
                    pointmark(pp[1]) - shift, pointmark(pp[2]) - shift);
          }
          if (!b->nobound) {
            // Output a boundary marker.
            fprintf(outfile, "  %d", marker);
          }
          if (b->neighout > 1) {
            fprintf(outfile, "    %5d  %5d", neigh1, neigh2);
          }
          fprintf(outfile, "\n");
        } else {
          // Output indices of three vertices.
          elist[index++] = pointmark(torg) - shift;
          elist[index++] = pointmark(tdest) - shift;
          elist[index++] = pointmark(tapex) - shift;
          if (b->order == 2) { // -o2
            out->o2facelist[o2index++] = pointmark(pp[0]) - shift;
            out->o2facelist[o2index++] = pointmark(pp[1]) - shift;
            out->o2facelist[o2index++] = pointmark(pp[2]) - shift;
          }
          if (!b->nobound) {
            emlist[facenumber - in->firstnumber] = marker;
          }
          if (b->neighout > 1) {
            out->face2tetlist[(facenumber - in->firstnumber) * 2]     = neigh1;
            out->face2tetlist[(facenumber - in->firstnumber) * 2 + 1] = neigh2;
          }
        }
        facenumber++;
      }
    }
    tface.tet = tetrahedrontraverse();
  }

  if (out == (TetMeshData *) NULL) {
    fprintf(outfile, "# Generated by %s\n", b->commandline);
    fclose(outfile);
  }

  if (b->neighout > 1) { // -nn option
    // Output the tetrahedron-to-face map.
	if (out == (TetMeshData *) NULL) {
	  strcpy(facefilename, b->outfilename);
      strcat(facefilename, ".t2f");
    }
	if (!b->quiet) {
      if (out == (TetMeshData *) NULL) {
        printf("Writing %s.\n", facefilename);
      } else {
        printf("Writing tetrahedron-to-face map.\n");
      }
    }
	if (out == (TetMeshData *) NULL) {
      outfile = fopen(facefilename, "w");
      for (tidx = 0; tidx < ntets; tidx++) {
        index = tidx * 4;
        fprintf(outfile, "%4d  %d %d %d %d\n", tidx + in->firstnumber,
                tet2facelist[index], tet2facelist[index+1], 
                tet2facelist[index+2], tet2facelist[index+3]);
      }
      fclose(outfile);
      delete [] tet2facelist;
    } else {
	  // Simply copy the address of the list to the output.
      out->tet2facelist = tet2facelist;
    }
  }
}

//============================================================================//
//                                                                            //
// outhullfaces()    Output hull faces to a .face file or a TetMeshData object.  //
//                                                                            //
// The normal of each face is pointing to the outside of the domain.          //
//                                                                            //
//============================================================================//

void TetMeshCore::outhullfaces(TetMeshData* out)
{
  FILE *outfile = NULL;
  char facefilename[FILENAMESIZE];
  triface hulltet;
  point torg, tdest, tapex;
  int *elist = NULL;
  int firstindex, shift;
  int facenumber;
  int index;

  if (out == (TetMeshData *) NULL) {
    strcpy(facefilename, b->outfilename);
    strcat(facefilename, ".face");
  }

  if (!b->quiet) {
    if (out == (TetMeshData *) NULL) {
      printf("Writing %s.\n", facefilename);
    } else {
      printf("Writing faces.\n");
    }
  }

  if (out == (TetMeshData *) NULL) {
    outfile = fopen(facefilename, "w");
    if (outfile == (FILE *) NULL) {
      printf("File I/O Error:  Cannot create file %s.\n", facefilename);
      terminate_tet_core(this, 1);
    }
    fprintf(outfile, "%ld  0\n", hullsize);
  } else {
    // Allocate memory for 'trifacelist'.
    out->trifacelist = new int[hullsize * 3];
    if (out->trifacelist == (int *) NULL) {
      printf("Error:  Out of memory.\n");
      terminate_tet_core(this, 1);
    }
    out->numberoftrifaces = hullsize;
    elist = out->trifacelist;
    index = 0;
  }

  // Determine the first index (0 or 1).
  firstindex = b->zeroindex ? 0 : in->firstnumber;
  shift = 0; // Default no shiftment.
  if ((in->firstnumber == 1) && (firstindex == 0)) {
    shift = 1; // Shift the output indices by 1.
  }

  tetrahedrons->traversalinit();
  hulltet.tet = alltetrahedrontraverse();
  facenumber = firstindex;
  while (hulltet.tet != (tetrahedron *) NULL) {
    if (ishulltet(hulltet)) {
      torg = (point) hulltet.tet[4];
      tdest = (point) hulltet.tet[5];
      tapex = (point) hulltet.tet[6];
      if (out == (TetMeshData *) NULL) {
        // Face number, indices of three vertices.
        fprintf(outfile, "%5d   %4d  %4d  %4d", facenumber,
                pointmark(torg) - shift, pointmark(tdest) - shift,
                pointmark(tapex) - shift);
        fprintf(outfile, "\n");
      } else {
        // Output indices of three vertices.
        elist[index++] = pointmark(torg) - shift;
        elist[index++] = pointmark(tdest) - shift;
        elist[index++] = pointmark(tapex) - shift;
      }
      facenumber++;
    }
    hulltet.tet = alltetrahedrontraverse();
  }

  if (out == (TetMeshData *) NULL) {
    fprintf(outfile, "# Generated by %s\n", b->commandline);
    fclose(outfile);
  }
}

//============================================================================//
//                                                                            //
// outsubfaces()    Output subfaces (i.e. boundary faces) to a .face file or  //
//                  a TetMeshData structure.                                     //
//                                                                            //
// The boundary faces are found in 'subfaces'. For listing triangle vertices  //
// in the same sense for all triangles in the mesh, the direction determined  //
// by right-hand rule is pointer to the inside of the volume.                 //
//                                                                            //
//============================================================================//

void TetMeshCore::outsubfaces(TetMeshData* out)
{
  FILE *outfile = NULL;
  char facefilename[FILENAMESIZE];
  int *elist = NULL;
  int *emlist = NULL;
  int index = 0, index1 = 0, index2 = 0;
  triface abuttingtet;
  face faceloop;
  point torg, tdest, tapex;
  int marker = 0;
  int firstindex, shift;
  int neigh1 = 0, neigh2 = 0;
  int facenumber;

  // For -o2 option.
  triface workface;
  point *extralist, pp[3] = {0,0,0}; 
  int highorderindex = 11;
  int o2index = 0, i;

  int t1ver; // used by fsymself()

  if (out == (TetMeshData *) NULL) {
    strcpy(facefilename, b->outfilename);
    strcat(facefilename, ".face");
  }

  if (!b->quiet) {
    if (out == (TetMeshData *) NULL) {
      printf("Writing %s.\n", facefilename);
    } else {
      printf("Writing faces.\n");
    }
  }

  if (out == (TetMeshData *) NULL) {
    outfile = fopen(facefilename, "w");
    if (outfile == (FILE *) NULL) {
      printf("File I/O Error:  Cannot create file %s.\n", facefilename);
      terminate_tet_core(this, 3);
    }
    // Number of subfaces.
    fprintf(outfile, "%ld  %d\n", subfaces->items, !b->nobound);
  } else {
    // Allocate memory for 'trifacelist'.
    out->trifacelist = new int[subfaces->items * 3];
    if (out->trifacelist == (int *) NULL) {
      terminate_tet_core(this, 1);
    }
    if (b->order == 2) {
      out->o2facelist = new int[subfaces->items * 3];
    }
    if (!b->nobound) {
      // Allocate memory for 'trifacemarkerlist'.
      out->trifacemarkerlist = new int[subfaces->items];
      if (out->trifacemarkerlist == (int *) NULL) {
        terminate_tet_core(this, 1);
      }
    }
    if (b->neighout > 1) {
      // '-nn' switch.
      out->face2tetlist = new int[subfaces->items * 2];
      if (out->face2tetlist == (int *) NULL) {
        terminate_tet_core(this, 1);
      }
    }
    out->numberoftrifaces = subfaces->items;
    elist = out->trifacelist;
    emlist = out->trifacemarkerlist;
  }

  // Determine the first index (0 or 1).
  firstindex = b->zeroindex ? 0 : in->firstnumber;
  shift = 0; // Default no shiftment.
  if ((in->firstnumber == 1) && (firstindex == 0)) {
    shift = 1; // Shift the output indices by 1.
  }

  subfaces->traversalinit();
  faceloop.sh = shellfacetraverse(subfaces);
  facenumber = firstindex; // in->firstnumber;
  while (faceloop.sh != (shellface *) NULL) {
    stpivot(faceloop, abuttingtet);
    // If there is a tetrahedron containing this subface, orient it so
    //   that the normal of this face points to inside of the volume by
    //   right-hand rule.
    if (abuttingtet.tet != NULL) {
      if (ishulltet(abuttingtet)) {
        fsymself(abuttingtet);
      }
    }
    if (abuttingtet.tet != NULL) {
      torg = org(abuttingtet);
      tdest = dest(abuttingtet);
      tapex = apex(abuttingtet);
      if (b->order == 2) { // -o2
        // Get the three extra vertices on edges.
        extralist = (point *) (abuttingtet.tet[highorderindex]);
        workface = abuttingtet;
        for (i = 0; i < 3; i++) {
          pp[i] = extralist[ver2edge[workface.ver]];
          enextself(workface);
        }
      }
    } else {
      // This may happen when only a surface mesh be generated.
      torg = sorg(faceloop);
      tdest = sdest(faceloop);
      tapex = sapex(faceloop);
      if (b->order == 2) { // -o2
        // There is no extra node list available.
        pp[0] = torg;
        pp[1] = tdest;
        pp[2] = tapex;
      }
    }
    if (!b->nobound) {
      marker = shellmark(faceloop);
    }
    if (b->neighout > 1) {
      // '-nn' switch. Output adjacent tets indices.
      neigh1 = -1;
      neigh2 = -1;
      stpivot(faceloop, abuttingtet);
      if (abuttingtet.tet != NULL) {
        if (!ishulltet(abuttingtet)) {
          neigh1 = elemindex(abuttingtet.tet);
        }
        fsymself(abuttingtet);
        if (!ishulltet(abuttingtet)) {
          neigh2 = elemindex(abuttingtet.tet);
        }
      }
    }
    if (out == (TetMeshData *) NULL) {
      fprintf(outfile, "%5d   %4d  %4d  %4d", facenumber,
              pointmark(torg) - shift, pointmark(tdest) - shift,
              pointmark(tapex) - shift);
      if (b->order == 2) { // -o2
        fprintf(outfile, "  %4d  %4d  %4d", pointmark(pp[0]) - shift, 
                pointmark(pp[1]) - shift, pointmark(pp[2]) - shift);
      }
      if (!b->nobound) {
        fprintf(outfile, "    %d", marker);
      }
      if (b->neighout > 1) {
        fprintf(outfile, "    %5d  %5d", neigh1, neigh2);
      }
      fprintf(outfile, "\n");
    } else {
      // Output three vertices of this face;
      elist[index++] = pointmark(torg) - shift;
      elist[index++] = pointmark(tdest) - shift;
      elist[index++] = pointmark(tapex) - shift;
      if (b->order == 2) { // -o2
        out->o2facelist[o2index++] = pointmark(pp[0]) - shift;
        out->o2facelist[o2index++] = pointmark(pp[1]) - shift;
        out->o2facelist[o2index++] = pointmark(pp[2]) - shift;
      }
      if (!b->nobound) {
        emlist[index1++] = marker;
      }
      if (b->neighout > 1) {
        out->face2tetlist[index2++] = neigh1;
        out->face2tetlist[index2++] = neigh2;
      }
    }
    facenumber++;
    faceloop.sh = shellfacetraverse(subfaces);
  }

  if (out == (TetMeshData *) NULL) {
    fprintf(outfile, "# Generated by %s\n", b->commandline);
    fclose(outfile);
  }
}

//============================================================================//
//                                                                            //
// outedges()    Output all edges to a .edge file or a TetMeshData object.       //
//                                                                            //
// Note: This routine must be called after outelements(),  so that the total  //
// number of edges 'meshedges' has been counted.                              //
//                                                                            //
//============================================================================//

void TetMeshCore::outedges(TetMeshData* out)
{
  FILE *outfile = NULL;
  char edgefilename[FILENAMESIZE];
  triface tetloop, worktet, spintet;
  face checkseg;
  point torg, tdest;
  int ishulledge;
  int firstindex, shift;
  int edgenumber, marker;
  int index = 0, index1 = 0, index2 = 0;
  int t1ver;
  int i;

  // For -o2 option.
  point *extralist, pp = NULL; 
  int highorderindex = 11;
  int o2index = 0;

  // For -nn option.
  int *tet2edgelist = NULL;
  int tidx;

  if (out == (TetMeshData *) NULL) {
    strcpy(edgefilename, b->outfilename);
    strcat(edgefilename, ".edge");
  }

  if (!b->quiet) {
    if (out == (TetMeshData *) NULL) {
      printf("Writing %s.\n", edgefilename);
    } else {
      printf("Writing edges.\n");
    }
  }

  if (meshedges == 0l) {
    if (nonconvex) {
      numberedges();  // Count the edges.
    } else {
      // Use Euler's characteristic to get the numbe of edges.
      // It states V - E + F - C = 1, hence E = V + F - C - 1.
      long tsize = tetrahedrons->items - hullsize;
      long fsize = (tsize * 4l + hullsize) / 2l;
      long vsize = points->items - dupverts - unuverts;
      if (b->weighted) vsize -= nonregularcount;
      meshedges = vsize + fsize - tsize - 1;
    }
  }
  meshhulledges = 0l; // It will be counted.

  if (out == (TetMeshData *) NULL) {
    outfile = fopen(edgefilename, "w");
    if (outfile == (FILE *) NULL) {
      printf("File I/O Error:  Cannot create file %s.\n", edgefilename);
      terminate_tet_core(this, 1);
    }
    // Write the number of edges, boundary markers (0 or 1).
    fprintf(outfile, "%ld  %d\n", meshedges, !b->nobound);
  } else {
    // Allocate memory for 'edgelist'.
    out->numberofedges = meshedges;
    out->edgelist = new int[meshedges * 2];
    if (out->edgelist == (int *) NULL) {
      printf("Error:  Out of memory.\n");
      terminate_tet_core(this, 1);
    }
    if (b->order == 2) { // -o2 switch
      out->o2edgelist = new int[meshedges];
    }
    if (!b->nobound) {
      out->edgemarkerlist = new int[meshedges];
    }
    if (b->neighout > 1) { // '-nn' switch.
      out->edge2tetlist = new int[meshedges];
    }
  }

  if (b->neighout > 1) { // -nn option
    // Output the tetrahedron-to-edge map.
	long tsize = tetrahedrons->items - hullsize;
    tet2edgelist = new int[tsize * 6];
  }

  // Determine the first index (0 or 1).
  firstindex = b->zeroindex ? 0 : in->firstnumber;
  shift = 0; // Default no shiftment.
  if ((in->firstnumber == 1) && (firstindex == 0)) {
    shift = 1; // Shift (reduce) the output indices by 1.
  }

  tetrahedrons->traversalinit();
  tetloop.tet = tetrahedrontraverse();
  edgenumber = firstindex; // in->firstnumber;
  while (tetloop.tet != (tetrahedron *) NULL) {
    // Count the number of Voronoi faces. 
    worktet.tet = tetloop.tet;
    for (i = 0; i < 6; i++) {
      worktet.ver = edge2ver[i];
      ishulledge = 0;
      fnext(worktet, spintet);
      do {
        if (!ishulltet(spintet)) {
          if (elemindex(spintet.tet) < elemindex(worktet.tet)) break;
        } else {
          ishulledge = 1;
        }
        fnextself(spintet);
      } while (spintet.tet != worktet.tet);
      if (spintet.tet == worktet.tet) {
        // Found a new edge.
        if (ishulledge) meshhulledges++;
        torg = org(worktet);
        tdest = dest(worktet);
        if (b->order == 2) { // -o2
          // Get the extra vertex on this edge.
          extralist = (point *) worktet.tet[highorderindex];
          pp = extralist[ver2edge[worktet.ver]];
        }
        if (out == (TetMeshData *) NULL) {
          fprintf(outfile, "%5d   %4d  %4d", edgenumber,
                  pointmark(torg) - shift, pointmark(tdest) - shift);
          if (b->order == 2) { // -o2
            fprintf(outfile, "  %4d", pointmark(pp) - shift);
          }
        } else {
          // Output three vertices of this face;
          out->edgelist[index++] = pointmark(torg) - shift;
          out->edgelist[index++] = pointmark(tdest) - shift;
          if (b->order == 2) { // -o2
            out->o2edgelist[o2index++] = pointmark(pp) - shift;
          }
        }
        if (!b->nobound) {
          if (b->plc || b->refine) {
            // Check if the edge is a segment.
            tsspivot1(worktet, checkseg);
            if (checkseg.sh != NULL) {
              marker = shellmark(checkseg);
            } else {
              marker = 0;  // It's not a segment.
            }
          } else {
            // Mark it if it is a hull edge.
            marker = ishulledge ? 1 : 0;
          }
          if (out == (TetMeshData *) NULL) {
            fprintf(outfile, "  %d", marker);
          } else {
            out->edgemarkerlist[index1++] = marker;
          }
        }
        if (b->neighout > 1) { // '-nn' switch.
          if (out == (TetMeshData *) NULL) {
            fprintf(outfile, "  %d", elemindex(tetloop.tet));
          } else {
            out->edge2tetlist[index2++] = elemindex(tetloop.tet);
          }
		  // Fill the tetrahedron-to-edge map.
		  spintet = worktet;
		  while (1) {
		    if (!ishulltet(spintet)) {
			  tidx = elemindex(spintet.tet) - firstindex;
			  tet2edgelist[tidx * 6 + ver2edge[spintet.ver]] = edgenumber;
			}
		    fnextself(spintet);
			if (spintet.tet == worktet.tet) break;
		  }
        }
        if (out == (TetMeshData *) NULL) {
          fprintf(outfile, "\n");
        }
        edgenumber++;
      }
    }
    tetloop.tet = tetrahedrontraverse();
  }

  if (out == (TetMeshData *) NULL) {
    fprintf(outfile, "# Generated by %s\n", b->commandline);
    fclose(outfile);
  }

  if (b->neighout > 1) { // -nn option
    long tsize = tetrahedrons->items - hullsize;

    if (b->facesout) { // -f option
      // Build the face-to-edge map (use the tet-to-edge map).
	  long fsize = (tsize * 4l + hullsize) / 2l;
	  int *face2edgelist = new int[fsize * 3];

	  tetrahedrons->traversalinit();
      tetloop.tet = tetrahedrontraverse();
      int facenumber = 0; // firstindex; // in->firstnumber;
	  while (tetloop.tet != (tetrahedron *) NULL) {
	    for (tetloop.ver = 0; tetloop.ver < 4; tetloop.ver++) {
		  fsym(tetloop, spintet);
		  if (ishulltet(spintet) || 
              (elemindex(tetloop.tet) < elemindex(spintet.tet))) {
            // The three edges of this face are ordered such that the
			//    first edge is opposite to the first vertex of this face
			//    that appears in the .face file, and so on.
			tidx = elemindex(tetloop.tet) - firstindex; 
			worktet = tetloop;
            for (i = 0; i < 3; i++) {
			  enextself(worktet); // The edge opposite to vertex i.
			  int eidx = tet2edgelist[tidx * 6 + ver2edge[worktet.ver]];
			  face2edgelist[facenumber * 3 + i] = eidx;
			}
		    facenumber++;
		  }
		}
	    tetloop.tet = tetrahedrontraverse();
	  }

	  // Output the face-to-edge map.
	  if (out == (TetMeshData *) NULL) {
	    strcpy(edgefilename, b->outfilename);
        strcat(edgefilename, ".f2e");
      }
	  if (!b->quiet) {
        if (out == (TetMeshData *) NULL) {
          printf("Writing %s.\n", edgefilename);
        } else {
          printf("Writing face-to-edge map.\n");
        }
      }
	  if (out == (TetMeshData *) NULL) {
	    outfile = fopen(edgefilename, "w");
	    for (tidx = 0; tidx < fsize; tidx++) { // Re-use `tidx'
	      i = tidx * 3;
		  fprintf(outfile, "%4d  %d %d %d\n", tidx + in->firstnumber,
		          face2edgelist[i], face2edgelist[i+1], face2edgelist[i+2]);
	    }
	    fclose(outfile);
	    delete [] face2edgelist;
	  } else {
	    // Simply copy the address of the list to the output.
	    out->face2edgelist = face2edgelist;
      }
    } // if (b->facesout)

    // Output the tetrahedron-to-edge map.
	if (out == (TetMeshData *) NULL) {
	  strcpy(edgefilename, b->outfilename);
      strcat(edgefilename, ".t2e");
    }
	if (!b->quiet) {
      if (out == (TetMeshData *) NULL) {
        printf("Writing %s.\n", edgefilename);
      } else {
        printf("Writing tetrahedron-to-edge map.\n");
      }
    }
	if (out == (TetMeshData *) NULL) {
	  outfile = fopen(edgefilename, "w");
	  for (tidx = 0; tidx < tsize; tidx++) {
	    i = tidx * 6;
		fprintf(outfile, "%4d  %d %d %d %d %d %d\n", tidx + in->firstnumber,
		        tet2edgelist[i], tet2edgelist[i+1], tet2edgelist[i+2],  
				tet2edgelist[i+3], tet2edgelist[i+4], tet2edgelist[i+5]);
	  }
	  fclose(outfile);
	  delete [] tet2edgelist;
	} else {
	  // Simply copy the address of the list to the output.
	  out->tet2edgelist = tet2edgelist;
    }
  }
}

//============================================================================//
//                                                                            //
// outsubsegments()    Output segments to a .edge file or a structure.        //
//                                                                            //
//============================================================================//

void TetMeshCore::outsubsegments(TetMeshData* out)
{
  FILE *outfile = NULL;
  char edgefilename[FILENAMESIZE];
  int *elist = NULL;
  int index, i;
  face edgeloop;
  point torg, tdest;
  int firstindex, shift;
  int marker;
  int edgenumber;

  // For -o2 option.
  triface workface, spintet;
  point *extralist, pp = NULL; 
  int highorderindex = 11;
  int o2index = 0;

  // For -nn option.
  int neigh = -1;
  int index2 = 0;

  int t1ver; // used by fsymself()

  if (out == (TetMeshData *) NULL) {
    strcpy(edgefilename, b->outfilename);
    strcat(edgefilename, ".edge");
  }

  if (!b->quiet) {
    if (out == (TetMeshData *) NULL) {
      printf("Writing %s.\n", edgefilename);
    } else {
      printf("Writing edges.\n");
    }
  }

  if (out == (TetMeshData *) NULL) {
    outfile = fopen(edgefilename, "w");
    if (outfile == (FILE *) NULL) {
      printf("File I/O Error:  Cannot create file %s.\n", edgefilename);
      terminate_tet_core(this, 3);
    }
    // Number of subsegments.
    fprintf(outfile, "%ld  1\n", subsegs->items);
  } else {
    // Allocate memory for 'edgelist'.
    out->edgelist = new int[subsegs->items * (b->order == 1 ? 2 : 3)];
    if (out->edgelist == (int *) NULL) {
      terminate_tet_core(this, 1);
    }
    if (b->order == 2) {
      out->o2edgelist = new int[subsegs->items];
    }
    out->edgemarkerlist = new int[subsegs->items];
    if (out->edgemarkerlist == (int *) NULL) {
      terminate_tet_core(this, 1);
    }
    if (b->neighout > 1) {
      out->edge2tetlist = new int[subsegs->items];
    }
    out->numberofedges = subsegs->items;
    elist = out->edgelist;
  }

  // Determine the first index (0 or 1).
  firstindex = b->zeroindex ? 0 : in->firstnumber;
  shift = 0; // Default no shiftment.
  if ((in->firstnumber == 1) && (firstindex == 0)) {
    shift = 1; // Shift the output indices by 1.
  }
  index = 0;
  i = 0;

  subsegs->traversalinit();
  edgeloop.sh = shellfacetraverse(subsegs);
  edgenumber = firstindex; // in->firstnumber;
  while (edgeloop.sh != (shellface *) NULL) {
    torg = sorg(edgeloop);
    tdest = sdest(edgeloop);
    if ((b->order == 2) || (b->neighout > 1)) {
      sstpivot1(edgeloop, workface);
      if (workface.tet != NULL) {
        // We must find a non-hull tet.
        if (ishulltet(workface)) {
          spintet = workface;
          while (1) {
            fnextself(spintet);
            if (!ishulltet(spintet)) break;
            if (spintet.tet == workface.tet) break;
          }
          workface = spintet;
        }
      }
    }
    if (b->order == 2) { // -o2
      // Get the extra vertex on this edge.
      if (workface.tet != NULL) {
        extralist = (point *) workface.tet[highorderindex];
        pp = extralist[ver2edge[workface.ver]];
      } else {
        pp = torg; // There is no extra node available.
      }
    }
    if (b->neighout > 1) { // -nn
      if (workface.tet != NULL) {
        neigh = elemindex(workface.tet);
      } else {
        neigh = -1;
      }
    }
    marker = shellmark(edgeloop);
    if (marker == 0) {
      marker = 1; // Default marker of a boundary edge is 1. 
    }
    if (out == (TetMeshData *) NULL) {
      fprintf(outfile, "%5d   %4d  %4d", edgenumber,
              pointmark(torg) - shift, pointmark(tdest) - shift);
      if (b->order == 2) { // -o2
        fprintf(outfile, "  %4d", pointmark(pp) - shift);
      }
      fprintf(outfile, "  %d", marker);
      if (b->neighout > 1) { // -nn
        fprintf(outfile, "  %4d", neigh);
      }
      fprintf(outfile, "\n");
    } else {
      // Output three vertices of this face;
      elist[index++] = pointmark(torg) - shift;
      elist[index++] = pointmark(tdest) - shift;
      if (b->order == 2) { // -o2
        out->o2edgelist[o2index++] = pointmark(pp) - shift;
      }
      out->edgemarkerlist[i++] = marker;
      if (b->neighout > 1) { // -nn
        out->edge2tetlist[index2++] = neigh;
      }
    }
    edgenumber++;
    edgeloop.sh = shellfacetraverse(subsegs);
  }

  if (out == (TetMeshData *) NULL) {
    fprintf(outfile, "# Generated by %s\n", b->commandline);
    fclose(outfile);
  }
}

//============================================================================//
//                                                                            //
// outneighbors()    Output tet neighbors to a .neigh file or a structure.    //
//                                                                            //
//============================================================================//

void TetMeshCore::outneighbors(TetMeshData* out)
{
  FILE *outfile = NULL;
  char neighborfilename[FILENAMESIZE];
  int *nlist = NULL;
  int index = 0;
  triface tetloop, tetsym;
  int neighbori[4];
  int firstindex;
  int elementnumber;
  long ntets;

  if (out == (TetMeshData *) NULL) {
    strcpy(neighborfilename, b->outfilename);
    strcat(neighborfilename, ".neigh");
  }

  if (!b->quiet) {
    if (out == (TetMeshData *) NULL) {
      printf("Writing %s.\n", neighborfilename);
    } else {
      printf("Writing neighbors.\n");
    }
  }

  ntets = tetrahedrons->items - hullsize;

  if (out == (TetMeshData *) NULL) {
    outfile = fopen(neighborfilename, "w");
    if (outfile == (FILE *) NULL) {
      printf("File I/O Error:  Cannot create file %s.\n", neighborfilename);
      terminate_tet_core(this, 1);
    }
    // Number of tetrahedra, four faces per tetrahedron.
    fprintf(outfile, "%ld  %d\n", ntets, 4);
  } else {
    // Allocate memory for 'neighborlist'.
    out->neighborlist = new int[ntets * 4];
    if (out->neighborlist == (int *) NULL) {
      printf("Error:  Out of memory.\n");
      terminate_tet_core(this, 1);
    }
    nlist = out->neighborlist;
  }

  // Determine the first index (0 or 1).
  firstindex = b->zeroindex ? 0 : in->firstnumber;

  tetrahedrons->traversalinit();
  tetloop.tet = tetrahedrontraverse();
  elementnumber = firstindex; // in->firstnumber;
  while (tetloop.tet != (tetrahedron *) NULL) {
    for (tetloop.ver = 0; tetloop.ver < 4; tetloop.ver++) {
      fsym(tetloop, tetsym);
      if (!ishulltet(tetsym)) {
        neighbori[tetloop.ver] = elemindex(tetsym.tet);
      } else {
        neighbori[tetloop.ver] = -1;
      }
    }
    if (out == (TetMeshData *) NULL) {
      // Tetrahedra number, neighboring tetrahedron numbers.
      fprintf(outfile, "%4d    %4d  %4d  %4d  %4d\n", elementnumber,
              neighbori[0], neighbori[1], neighbori[2], neighbori[3]);
    } else {
      nlist[index++] = neighbori[0];
      nlist[index++] = neighbori[1];
      nlist[index++] = neighbori[2];
      nlist[index++] = neighbori[3];
    }
    tetloop.tet = tetrahedrontraverse();
    elementnumber++;
  }

  if (out == (TetMeshData *) NULL) {
    fprintf(outfile, "# Generated by %s\n", b->commandline);
    fclose(outfile);
  }
}

//============================================================================//
//                                                                            //
// outvoronoi()    Output the Voronoi diagram to .v.node, .v.edge, v.face,    //
//                 and .v.cell.                                               //
//                                                                            //
// The Voronoi diagram is the geometric dual of the Delaunay triangulation.   //
// The Voronoi vertices are the circumcenters of Delaunay tetrahedra.  Each   //
// Voronoi edge connects two Voronoi vertices at two sides of a common Dela-  //
// unay face. At a face of convex hull, it becomes a ray (goto the infinity). //
// A Voronoi face is the convex hull of all Voronoi vertices around a common  //
// Delaunay edge. It is a closed polygon for any internal Delaunay edge. At a //
// ridge, it is unbounded.  Each Voronoi cell is the convex hull of all Vor-  //
// onoi vertices around a common Delaunay vertex. It is a polytope for any    //
// internal Delaunay vertex. It is an unbounded polyhedron for a Delaunay     //
// vertex belonging to the convex hull.                                       //
//                                                                            //
// NOTE: This routine is only used when the input is only a set of point.     //
// Comment: Special thanks to Victor Liu for finding and fixing few bugs.     //
//                                                                            //
//============================================================================//

void TetMeshCore::outvoronoi(TetMeshData* out)
{
  FILE *outfile = NULL;
  char outfilename[FILENAMESIZE];
  TetMeshData::voroedge *vedge = NULL;
  TetMeshData::vorofacet *vfacet = NULL;
  arraypool *tetlist, *ptlist;
  triface tetloop, worktet, spintet, firsttet;
  point pt[4], ploop, neipt;
  double ccent[3], infvec[3], vec1[3], vec2[3], L;
  long ntets, faces, edges;
  int *indexarray, *fidxs, *eidxs;
  int arraysize, *vertarray = NULL;
  int vpointcount, vedgecount, vfacecount, tcount;
  int ishullvert, ishullface;
  int index, shift, end1, end2;
  int i, j;

  int t1ver; // used by fsymself()

  // Output Voronoi vertices to .v.node file.
  if (out == (TetMeshData *) NULL) {
    strcpy(outfilename, b->outfilename);
    strcat(outfilename, ".v.node");
  }

  if (!b->quiet) {
    if (out == (TetMeshData *) NULL) {
      printf("Writing %s.\n", outfilename);
    } else {
      printf("Writing Voronoi vertices.\n");
    }
  }

  // Determine the first index (0 or 1).
  shift = (b->zeroindex ? 0 : in->firstnumber);

  // Each face and edge of the tetrahedral mesh will be indexed for indexing
  //   the Voronoi edges and facets. Indices of faces and edges are saved in
  //   each tetrahedron (including hull tets).

  // Allocate the total space once.
  indexarray = new int[tetrahedrons->items * 10];

  // Allocate space (10 integers) into each tetrahedron. It re-uses the slot
  //   for element markers, flags.
  i = 0;
  tetrahedrons->traversalinit();
  tetloop.tet = alltetrahedrontraverse();
  while (tetloop.tet != NULL) {
    tetloop.tet[11] = (tetrahedron) &(indexarray[i * 10]);
    i++;
    tetloop.tet = alltetrahedrontraverse();
  }

  // The number of tetrahedra (excluding hull tets) (Voronoi vertices).
  ntets = tetrahedrons->items - hullsize;
  // The number of Delaunay faces (Voronoi edges).
  faces = (4l * ntets + hullsize) / 2l;
  // The number of Delaunay edges (Voronoi faces).
  long vsize = points->items - dupverts - unuverts;
  if (b->weighted) vsize -= nonregularcount;
  if (!nonconvex) {
    edges = vsize + faces - ntets - 1;
  } else {
    if (meshedges == 0l) {
      numberedges(); // Count edges.
    }
    edges = meshedges; 
  }

  if (out == (TetMeshData *) NULL) {
    outfile = fopen(outfilename, "w");
    if (outfile == (FILE *) NULL) {
      printf("File I/O Error:  Cannot create file %s.\n", outfilename);
      terminate_tet_core(this, 3);
    }
    // Number of voronoi points, 3 dim, no attributes, no marker.
    fprintf(outfile, "%ld  3  0  0\n", ntets);
  } else {
    // Allocate space for 'vpointlist'.
    out->numberofvpoints = (int) ntets;
    out->vpointlist = new double[out->numberofvpoints * 3];
    if (out->vpointlist == (double *) NULL) {
      terminate_tet_core(this, 1);
    }
  }

  // Output Voronoi vertices (the circumcenters of tetrahedra). 
  tetrahedrons->traversalinit();
  tetloop.tet = tetrahedrontraverse();
  vpointcount = 0; // The (internal) v-index always starts from 0. 
  index = 0;
  while (tetloop.tet != (tetrahedron *) NULL) {
    for (i = 0; i < 4; i++) {
      pt[i] = (point) tetloop.tet[4 + i];
      setpoint2tet(pt[i], encode(tetloop));
    }
    if (b->weighted) {
      orthosphere(pt[0], pt[1], pt[2], pt[3], pt[0][3], pt[1][3], pt[2][3], 
                  pt[3][3], ccent, NULL);
    } else {
      circumsphere(pt[0], pt[1], pt[2], pt[3], ccent, NULL);
    }
    if (out == (TetMeshData *) NULL) {
      fprintf(outfile, "%4d  %16.8e %16.8e %16.8e\n", vpointcount + shift,
              ccent[0], ccent[1], ccent[2]);
    } else {
      out->vpointlist[index++] = ccent[0];
      out->vpointlist[index++] = ccent[1];
      out->vpointlist[index++] = ccent[2];
    }
    setelemindex(tetloop.tet, vpointcount);
    vpointcount++;
    tetloop.tet = tetrahedrontraverse();
  }

  if (out == (TetMeshData *) NULL) {
    fprintf(outfile, "# Generated by %s\n", b->commandline);
    fclose(outfile);
  }

  // Output Voronoi edges to .v.edge file.
  if (out == (TetMeshData *) NULL) {
    strcpy(outfilename, b->outfilename);
    strcat(outfilename, ".v.edge");
  }
  
  if (!b->quiet) {
    if (out == (TetMeshData *) NULL) {
      printf("Writing %s.\n", outfilename);
    } else {
      printf("Writing Voronoi edges.\n");
    }
  }

  if (out == (TetMeshData *) NULL) {
    outfile = fopen(outfilename, "w");
    if (outfile == (FILE *) NULL) {
      printf("File I/O Error:  Cannot create file %s.\n", outfilename);
      terminate_tet_core(this, 3);
    }
    // Number of Voronoi edges, no marker.
    fprintf(outfile, "%ld  0\n", faces);
  } else {
    // Allocate space for 'vpointlist'.
    out->numberofvedges = (int) faces;
    out->vedgelist = new TetMeshData::voroedge[out->numberofvedges];
  }

  // Output the Voronoi edges. 
  tetrahedrons->traversalinit();
  tetloop.tet = tetrahedrontraverse();
  vedgecount = 0; // D-Face (V-edge) index (from zero).
  index = 0; // The Delaunay-face index.
  while (tetloop.tet != (tetrahedron *) NULL) {
    // Count the number of Voronoi edges. Look at the four faces of each
    //   tetrahedron. Count the face if the tetrahedron's index is
    //   smaller than its neighbor's or the neighbor is outside.
    end1 = elemindex(tetloop.tet);
    for (tetloop.ver = 0; tetloop.ver < 4; tetloop.ver++) {
      fsym(tetloop, worktet);
      if (ishulltet(worktet) || 
          (elemindex(tetloop.tet) < elemindex(worktet.tet))) {
        // Found a Voronoi edge. Operate on it.
        if (out == (TetMeshData *) NULL) {
          fprintf(outfile, "%4d  %4d", vedgecount + shift, end1 + shift);
        } else {
          vedge = &(out->vedgelist[index++]);
          vedge->v1 = end1 + shift;
        }
        if (!ishulltet(worktet)) {
          end2 = elemindex(worktet.tet);
        } else {
          end2 = -1;
        }
        // Note that end2 may be -1 (worktet.tet is outside).
        if (end2 == -1) {
          // Calculate the out normal of this hull face.
          pt[0] = dest(worktet);
          pt[1] = org(worktet);
          pt[2] = apex(worktet);
          for (j = 0; j < 3; j++) vec1[j] = pt[1][j] - pt[0][j];
          for (j = 0; j < 3; j++) vec2[j] = pt[2][j] - pt[0][j];
          cross(vec1, vec2, infvec);
          // Normalize it.
          L = sqrt(infvec[0] * infvec[0] + infvec[1] * infvec[1]
                   + infvec[2] * infvec[2]);
          if (L > 0) for (j = 0; j < 3; j++) infvec[j] /= L;
          if (out == (TetMeshData *) NULL) {
            fprintf(outfile, " -1");
            fprintf(outfile, " %g %g %g\n", infvec[0], infvec[1], infvec[2]);
          } else {
            vedge->v2 = -1;
            vedge->vnormal[0] = infvec[0];
            vedge->vnormal[1] = infvec[1];
            vedge->vnormal[2] = infvec[2];
          }
        } else {
          if (out == (TetMeshData *) NULL) {
            fprintf(outfile, " %4d\n", end2 + shift);
          } else {
            vedge->v2 = end2 + shift;
            vedge->vnormal[0] = 0.0;
            vedge->vnormal[1] = 0.0;
            vedge->vnormal[2] = 0.0;
          }
        }
        // Save the V-edge index in this tet and its neighbor.
        fidxs = (int *) (tetloop.tet[11]);
        fidxs[tetloop.ver] = vedgecount;
        fidxs = (int *) (worktet.tet[11]);
        fidxs[worktet.ver & 3] = vedgecount;
        vedgecount++;
      }
    } // tetloop.ver
    tetloop.tet = tetrahedrontraverse();
  }

  if (out == (TetMeshData *) NULL) {
    fprintf(outfile, "# Generated by %s\n", b->commandline);
    fclose(outfile);
  }

  // Output Voronoi faces to .v.face file.
  if (out == (TetMeshData *) NULL) {
    strcpy(outfilename, b->outfilename);
    strcat(outfilename, ".v.face");
  }
  
  if (!b->quiet) {
    if (out == (TetMeshData *) NULL) {
      printf("Writing %s.\n", outfilename);
    } else {
      printf("Writing Voronoi faces.\n");
    }
  }

  if (out == (TetMeshData *) NULL) {
    outfile = fopen(outfilename, "w");
    if (outfile == (FILE *) NULL) {
      printf("File I/O Error:  Cannot create file %s.\n", outfilename);
      terminate_tet_core(this, 3);
    }
    // Number of Voronoi faces.
    fprintf(outfile, "%ld  0\n", edges);
  } else {
    out->numberofvfacets = edges;
    out->vfacetlist = new TetMeshData::vorofacet[out->numberofvfacets];
    if (out->vfacetlist == (TetMeshData::vorofacet *) NULL) {
      terminate_tet_core(this, 1);
    }
  }

  // Output the Voronoi facets.
  tetrahedrons->traversalinit();
  tetloop.tet = tetrahedrontraverse();
  vfacecount = 0; // D-edge (V-facet) index (from zero).
  while (tetloop.tet != (tetrahedron *) NULL) {
    // Count the number of Voronoi faces. Look at the six edges of each
    //   tetrahedron. Count the edge only if the tetrahedron's index is
    //   smaller than those of all other tetrahedra that share the edge.
    worktet.tet = tetloop.tet;
    for (i = 0; i < 6; i++) {
      worktet.ver = edge2ver[i];
      // Count the number of faces at this edge. If the edge is a hull edge,
      //   the face containing dummypoint is also counted.
      //ishulledge = 0; // Is it a hull edge.
      tcount = 0;
      firsttet = worktet;
      spintet = worktet;
      while (1) {
        tcount++;
        fnextself(spintet);
        if (spintet.tet == worktet.tet) break;
        if (!ishulltet(spintet)) {
          if (elemindex(spintet.tet) < elemindex(worktet.tet)) break;
        } else {
          //ishulledge = 1;
          if (apex(spintet) == dummypoint) {
            // We make this V-edge appear in the end of the edge list.
            fnext(spintet, firsttet); 
          }
        }
      } // while (1)
      if (spintet.tet == worktet.tet) {
        // Found a Voronoi facet. Operate on it.
        pt[0] = org(worktet);
        pt[1] = dest(worktet);
        end1 = pointmark(pt[0]) - in->firstnumber; // V-cell index
        end2 = pointmark(pt[1]) - in->firstnumber;
        if (out == (TetMeshData *) NULL) {
          fprintf(outfile, "%4d  %4d %4d  %-2d ", vfacecount + shift, 
                  end1 + shift, end2 + shift, tcount);
        } else {
          vfacet = &(out->vfacetlist[vfacecount]);
          vfacet->c1 = end1 + shift;
          vfacet->c2 = end2 + shift;
          vfacet->elist = new int[tcount + 1];
          vfacet->elist[0] = tcount;
          index = 1;
        }
        // Output V-edges of this V-facet.
        spintet = firsttet; //worktet;
        while (1) {
          fidxs = (int *) (spintet.tet[11]);
          if (apex(spintet) != dummypoint) {
            vedgecount = fidxs[spintet.ver & 3];
            ishullface = 0;
          } else {
            ishullface = 1; // It's not a real face.
          }
          if (out == (TetMeshData *) NULL) {
            fprintf(outfile, " %d", !ishullface ? (vedgecount + shift) : -1); 
          } else {
            vfacet->elist[index++] = !ishullface ? (vedgecount + shift) : -1;
          }
          // Save the V-facet index in this tet at this edge.
          eidxs = &(fidxs[4]);
          eidxs[ver2edge[spintet.ver]] = vfacecount;
          // Go to the next face.
          fnextself(spintet);
          if (spintet.tet == firsttet.tet) break;
        } // while (1)
        if (out == (TetMeshData *) NULL) {
          fprintf(outfile, "\n");
        }
        vfacecount++;
      } // if (spintet.tet == worktet.tet)
    } // if (i = 0; i < 6; i++)
    tetloop.tet = tetrahedrontraverse();
  }

  if (out == (TetMeshData *) NULL) {
    fprintf(outfile, "# Generated by %s\n", b->commandline);
    fclose(outfile);
  }

  // Output Voronoi cells to .v.cell file.
  if (out == (TetMeshData *) NULL) {
    strcpy(outfilename, b->outfilename);
    strcat(outfilename, ".v.cell");
  }
  
  if (!b->quiet) {
    if (out == (TetMeshData *) NULL) {
      printf("Writing %s.\n", outfilename);
    } else {
      printf("Writing Voronoi cells.\n");
    }
  }

  if (out == (TetMeshData *) NULL) {
    outfile = fopen(outfilename, "w");
    if (outfile == (FILE *) NULL) {
      printf("File I/O Error:  Cannot create file %s.\n", outfilename);
      terminate_tet_core(this, 3);
    }
    // Number of Voronoi cells.
    fprintf(outfile, "%ld\n", points->items - unuverts - dupverts);
  } else {
    out->numberofvcells = points->items - unuverts - dupverts;
    out->vcelllist = new int*[out->numberofvcells];
    if (out->vcelllist == (int **) NULL) {
      terminate_tet_core(this, 1);
    }
  }

  // Output Voronoi cells.
  tetlist = cavetetlist;
  ptlist = cavetetvertlist;
  points->traversalinit();
  ploop = pointtraverse();
  vpointcount = 0;
  while (ploop != (point) NULL) {
    if ((pointtype(ploop) != UNUSEDVERTEX) &&
        (pointtype(ploop) != DUPLICATEDVERTEX) &&
        (pointtype(ploop) != NREGULARVERTEX)) { 
      getvertexstar(1, ploop, tetlist, ptlist, NULL);
      // Mark all vertices. Check if it is a hull vertex.
      ishullvert = 0;
      for (i = 0; i < ptlist->objects; i++) {
        neipt = * (point *) fastlookup(ptlist, i);
        if (neipt != dummypoint) {
          pinfect(neipt);
        } else {
          ishullvert = 1;
        }
      }
      tcount = (int) ptlist->objects;
      if (out == (TetMeshData *) NULL) {
        fprintf(outfile, "%4d  %-2d ", vpointcount + shift, tcount);
      } else {
        arraysize = tcount;
        vertarray = new int[arraysize + 1];
        out->vcelllist[vpointcount] = vertarray;
        vertarray[0] = tcount;
        index = 1;
      }
      // List Voronoi facets bounding this cell.
      for (i = 0; i < tetlist->objects; i++) {
        worktet = * (triface *) fastlookup(tetlist, i);
        // Let 'worktet' be [a,b,c,d] where d = ploop.
        for (j = 0; j < 3; j++) {
          neipt = org(worktet); // neipt is a, or b, or c
          // Skip the dummypoint.
          if (neipt != dummypoint) {
            if (pinfected(neipt)) {
              // It's not processed yet.
              puninfect(neipt);
              // Go to the DT edge [a,d], or [b,d], or [c,d]. 
              esym(worktet, spintet);
              enextself(spintet);
              // Get the V-face dual to this edge.
              eidxs = (int *) spintet.tet[11];
              vfacecount = eidxs[4 + ver2edge[spintet.ver]];
              if (out == (TetMeshData *) NULL) {
                fprintf(outfile, " %d", vfacecount + shift);
              } else {
                vertarray[index++] = vfacecount + shift;
              }
            }
          }
          enextself(worktet);
        } // j
      } // i
      if (ishullvert) {
        // Add a hull facet (-1) to the facet list.
        if (out == (TetMeshData *) NULL) {
          fprintf(outfile, " -1");
        } else {
          vertarray[index++] = -1;
        }
      }
      if (out == (TetMeshData *) NULL) {
        fprintf(outfile, "\n");
      }
      tetlist->restart();
      ptlist->restart();
      vpointcount++;
    }
    ploop = pointtraverse();
  }

  // Delete the space for face/edge indices.
  delete [] indexarray;

  if (out == (TetMeshData *) NULL) {
    fprintf(outfile, "# Generated by %s\n", b->commandline);
    fclose(outfile);
  }
}

//============================================================================//
//                                                                            //


//                                                                            //
//                                                                            //
//== output_cxx ==============================================================//

} // namespace sqmesh::mesh::tet::detail
