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

//== delaunay_cxx ============================================================//
//                                                                            //
//                                                                            //

//============================================================================//
//                                                                            //
// transfernodes()    Read the vertices from the input (TetMeshData).            //
//                                                                            //
// Transferring all points from input ('in->pointlist') to the tet core's 'points'. //
// All points are indexed (the first point index is 'in->firstnumber'). Each  //
// point's type is initialized as UNUSEDVERTEX. The bounding box (xmax, xmin, //
// ...) and the diameter (longest) of the point set are calculated.           //
//                                                                            //
//============================================================================//

void TetMeshCore::transfernodes()
{
  point pointloop;
  double x, y, z, w, mtr;
  int coordindex;
  int attribindex;
  int mtrindex;
  int i, j;

  // Read the points.
  coordindex = 0;
  attribindex = 0;
  mtrindex = 0;
  for (i = 0; i < in->numberofpoints; i++) {
    makepoint(&pointloop, UNUSEDVERTEX);
    // Read the point coordinates.
    x = pointloop[0] = in->pointlist[coordindex++];
    y = pointloop[1] = in->pointlist[coordindex++];
    z = pointloop[2] = in->pointlist[coordindex++];
    // Read the point attributes. (Including point weights.)
    for (j = 0; j < in->numberofpointattributes; j++) {
      pointloop[3 + j] = in->pointattributelist[attribindex++];
    }
    // Read the point metric tensor.
    for (j = 0; j < in->numberofpointmtrs; j++) {
      mtr = in->pointmtrlist[mtrindex++] * b->metric_scale;
      pointloop[pointmtrindex + j] = mtr; // in->pointmtrlist[mtrindex++];
    }
    if (b->weighted) { // -w option
      if (in->numberofpointattributes > 0) {
        // The first point attribute is its weight.
        //w = in->pointattributelist[in->numberofpointattributes * i];
        w = pointloop[3];
      } else {
        // No given weight available. Default choose the maximum
        //   absolute value among its coordinates.        
        w = fabs(x);
        if (w < fabs(y)) w = fabs(y);
        if (w < fabs(z)) w = fabs(z);
      }
      if (b->weighted_param == 0) {
        pointloop[3] = x * x + y * y + z * z - w; // Weighted DT.
      } else { // -w1 option
        pointloop[3] = w;  // Regular tetrahedralization.
      }
    }
    // Determine the smallest and largest x, y and z coordinates.
    if (i == 0) {
      xmin = xmax = x;
      ymin = ymax = y;
      zmin = zmax = z;
    } else {
      xmin = (x < xmin) ? x : xmin;
      xmax = (x > xmax) ? x : xmax;
      ymin = (y < ymin) ? y : ymin;
      ymax = (y > ymax) ? y : ymax;
      zmin = (z < zmin) ? z : zmin;
      zmax = (z > zmax) ? z : zmax;
    }
  }

  x = xmax - xmin;
  y = ymax - ymin;
  z = zmax - zmin;

  exactinit(b->verbose, b->noexact, b->nostaticfilter, x, y, z);

  // Use the number of points as the random seed.
  srand(in->numberofpoints);

  // 'longest' is the largest possible edge length formed by input vertices.
  longest = sqrt(x * x + y * y + z * z);
  if (longest == 0.0) {
    printf("Error:  The point set is trivial.\n");
    terminate_tet_core(this, 10);
  }
  // Two identical points are distinguished by 'minedgelength'.
  minedgelength = longest * b->epsilon;
}

//============================================================================//
//                                                                            //
// hilbert_init()    Initialize the Gray code permutation table.              //
//                                                                            //
// The table 'transgc' has 8 x 3 x 8 entries. It contains all possible Gray   //
// code sequences traveled by the 1st order Hilbert curve in 3 dimensions.    //
// The first column is the Gray code of the entry point of the curve, and     //
// the second column is the direction (0, 1, or 2, 0 means the x-axis) where  //
// the exit point of curve lies.                                              //
//                                                                            //
// The table 'tsb1mod3' contains the numbers of trailing set '1' bits of the  //
// indices from 0 to 7, modulo by '3'. The code for generating this table is  //
// from: http://graphics.stanford.edu/~seander/bithacks.html.                 //
//                                                                            //
//============================================================================//

void TetMeshCore::hilbert_init(int n)
{
  int gc[8], N, mask, travel_bit;
  int e, d, f, k, g;
  int v, c;
  int i;

  N = (n == 2) ? 4 : 8;
  mask = (n == 2) ? 3 : 7;

  // Generate the Gray code sequence.
  for (i = 0; i < N; i++) {
    gc[i] = i ^ (i >> 1);
  }

  for (e = 0; e < N; e++) {
    for (d = 0; d < n; d++) {
      // Calculate the end point (f).
      f = e ^ (1 << d);  // Toggle the d-th bit of 'e'.
      // travel_bit = 2**p, the bit we want to travel. 
      travel_bit = e ^ f;
      for (i = 0; i < N; i++) {
        // // Rotate gc[i] left by (p + 1) % n bits.
        k = gc[i] * (travel_bit * 2);
        g = ((k | (k / N)) & mask);
        // Calculate the permuted Gray code by xor with the start point (e).
        transgc[e][d][i] = (g ^ e);
      }
    } // d
  } // e

  // Count the consecutive '1' bits (trailing) on the right.
  tsb1mod3[0] = 0;
  for (i = 1; i < N; i++) {
    v = ~i; // Count the 0s.
    v = (v ^ (v - 1)) >> 1; // Set v's trailing 0s to 1s and zero rest
    for (c = 0; v; c++) {
      v >>= 1;
    }
    tsb1mod3[i] = c % n;
  }
}

//============================================================================//
//                                                                            //
// hilbert_sort3()    Sort points using the 3d Hilbert curve.                 //
//                                                                            //
//============================================================================//

int TetMeshCore::hilbert_split(point* vertexarray,int arraysize,int gc0,int gc1,
                              double bxmin, double bxmax, double bymin, double bymax, 
                              double bzmin, double bzmax)
{
  point swapvert;
  int axis, d;
  double split;
  int i, j;


  // Find the current splitting axis. 'axis' is a value 0, or 1, or 2, which 
  //   correspoding to x-, or y- or z-axis.
  axis = (gc0 ^ gc1) >> 1; 

  // Calulate the split position along the axis.
  if (axis == 0) {
    split = 0.5 * (bxmin + bxmax);
  } else if (axis == 1) {
    split = 0.5 * (bymin + bymax);
  } else { // == 2
    split = 0.5 * (bzmin + bzmax);
  }

  // Find the direction (+1 or -1) of the axis. If 'd' is +1, the direction
  //   of the axis is to the positive of the axis, otherwise, it is -1.
  d = ((gc0 & (1<<axis)) == 0) ? 1 : -1;


  // Partition the vertices into left- and right-arrays such that left points
  //   have Hilbert indices lower than the right points.
  i = 0;
  j = arraysize - 1;

  // Partition the vertices into left- and right-arrays.
  if (d > 0) {
    do {
      for (; i < arraysize; i++) {      
        if (vertexarray[i][axis] >= split) break;
      }
      for (; j >= 0; j--) {
        if (vertexarray[j][axis] < split) break;
      }
      // Is the partition finished?
      if (i == (j + 1)) break;
      // Swap i-th and j-th vertices.
      swapvert = vertexarray[i];
      vertexarray[i] = vertexarray[j];
      vertexarray[j] = swapvert;
      // Continue patitioning the array;
    } while (true);
  } else {
    do {
      for (; i < arraysize; i++) {      
        if (vertexarray[i][axis] <= split) break;
      }
      for (; j >= 0; j--) {
        if (vertexarray[j][axis] > split) break;
      }
      // Is the partition finished?
      if (i == (j + 1)) break;
      // Swap i-th and j-th vertices.
      swapvert = vertexarray[i];
      vertexarray[i] = vertexarray[j];
      vertexarray[j] = swapvert;
      // Continue patitioning the array;
    } while (true);
  }

  return i;
}

void TetMeshCore::hilbert_sort3(point* vertexarray, int arraysize, int e, int d, 
                               double bxmin, double bxmax, double bymin, double bymax, 
                               double bzmin, double bzmax, int depth)
{
  double x1, x2, y1, y2, z1, z2;
  int p[9], w, e_w, d_w, k, ei, di;
  int n = 3, mask = 7;

  p[0] = 0;
  p[8] = arraysize;

  // Sort the points according to the 1st order Hilbert curve in 3d.
  p[4] = hilbert_split(vertexarray, p[8], transgc[e][d][3], transgc[e][d][4], 
                       bxmin, bxmax, bymin, bymax, bzmin, bzmax);
  p[2] = hilbert_split(vertexarray, p[4], transgc[e][d][1], transgc[e][d][2], 
                       bxmin, bxmax, bymin, bymax, bzmin, bzmax);
  p[1] = hilbert_split(vertexarray, p[2], transgc[e][d][0], transgc[e][d][1], 
                       bxmin, bxmax, bymin, bymax, bzmin, bzmax);
  p[3] = hilbert_split(&(vertexarray[p[2]]), p[4] - p[2], 
                       transgc[e][d][2], transgc[e][d][3], 
                       bxmin, bxmax, bymin, bymax, bzmin, bzmax) + p[2];
  p[6] = hilbert_split(&(vertexarray[p[4]]), p[8] - p[4], 
                       transgc[e][d][5], transgc[e][d][6], 
                       bxmin, bxmax, bymin, bymax, bzmin, bzmax) + p[4];
  p[5] = hilbert_split(&(vertexarray[p[4]]), p[6] - p[4], 
                       transgc[e][d][4], transgc[e][d][5], 
                       bxmin, bxmax, bymin, bymax, bzmin, bzmax) + p[4];
  p[7] = hilbert_split(&(vertexarray[p[6]]), p[8] - p[6], 
                       transgc[e][d][6], transgc[e][d][7], 
                       bxmin, bxmax, bymin, bymax, bzmin, bzmax) + p[6];

  if (b->hilbert_order > 0) {
    // A maximum order is prescribed. 
    if ((depth + 1) == b->hilbert_order) {
      // The maximum prescribed order is reached.
      return;
    }
  }

  // Recursively sort the points in sub-boxes.
  for (w = 0; w < 8; w++) {
    // w is the local Hilbert index (NOT Gray code).
    // Sort into the sub-box either there are more than 2 points in it, or
    //   the prescribed order of the curve is not reached yet.
    //if ((p[w+1] - p[w] > b->hilbert_limit) || (b->hilbert_order > 0)) {
    if ((p[w+1] - p[w]) > b->hilbert_limit) {
      // Calculcate the start point (ei) of the curve in this sub-box.
      //   update e = e ^ (e(w) left_rotate (d+1)).
      if (w == 0) {
        e_w = 0;
      } else {
        //   calculate e(w) = gc(2 * floor((w - 1) / 2)).
        k = 2 * ((w - 1) / 2); 
        e_w = k ^ (k >> 1); // = gc(k).
      }
      k = e_w;
      e_w = ((k << (d+1)) & mask) | ((k >> (n-d-1)) & mask);
      ei = e ^ e_w;
      // Calulcate the direction (di) of the curve in this sub-box.
      //   update d = (d + d(w) + 1) % n
      if (w == 0) {
        d_w = 0;
      } else {
        d_w = ((w % 2) == 0) ? tsb1mod3[w - 1] : tsb1mod3[w];
      }
      di = (d + d_w + 1) % n;
      // Calculate the bounding box of the sub-box.
      if (transgc[e][d][w] & 1) { // x-axis
        x1 = 0.5 * (bxmin + bxmax);
        x2 = bxmax;
      } else {
        x1 = bxmin;
        x2 = 0.5 * (bxmin + bxmax);
      }
      if (transgc[e][d][w] & 2) { // y-axis
        y1 = 0.5 * (bymin + bymax);
        y2 = bymax;
      } else {
        y1 = bymin;
        y2 = 0.5 * (bymin + bymax);
      }
      if (transgc[e][d][w] & 4) { // z-axis
        z1 = 0.5 * (bzmin + bzmax);
        z2 = bzmax;
      } else {
        z1 = bzmin;
        z2 = 0.5 * (bzmin + bzmax);
      }
      hilbert_sort3(&(vertexarray[p[w]]), p[w+1] - p[w], ei, di, 
                    x1, x2, y1, y2, z1, z2, depth+1);
    } // if (p[w+1] - p[w] > 1)
  } // w
}

//============================================================================//
//                                                                            //
// brio_multiscale_sort()    Sort the points using BRIO and Hilbert curve.    //
//                                                                            //
//============================================================================//

void TetMeshCore::brio_multiscale_sort(point* vertexarray, int arraysize, 
                                      int threshold, double ratio, int *depth)
{
  int middle;

  middle = 0;
  if (arraysize >= threshold) {
    (*depth)++;
    middle = arraysize * ratio;
    brio_multiscale_sort(vertexarray, middle, threshold, ratio, depth);
  }
  // Sort the right-array (rnd-th round) using the Hilbert curve.
  hilbert_sort3(&(vertexarray[middle]), arraysize - middle, 0, 0, // e, d
                xmin, xmax, ymin, ymax, zmin, zmax, 0); // depth.
}

//============================================================================//
//                                                                            //
// randomnation()    Generate a random number between 0 and 'choices' - 1.    //
//                                                                            //
//============================================================================//

unsigned long TetMeshCore::randomnation(unsigned int choices)
{
  unsigned long newrandom;

  if (choices >= 714025l) {
    newrandom = (randomseed * 1366l + 150889l) % 714025l;
    randomseed = (newrandom * 1366l + 150889l) % 714025l;
    newrandom = newrandom * (choices / 714025l) + randomseed;
    if (newrandom >= choices) {
      return newrandom - choices;
    } else {
      return newrandom;
    }
  } else {
    randomseed = (randomseed * 1366l + 150889l) % 714025l;
    return randomseed % choices;
  }
}

//============================================================================//
//                                                                            //
// randomsample()    Randomly sample the tetrahedra for point loation.        //
//                                                                            //
// Searching begins from one of handles:  the input 'searchtet', a recently   //
// encountered tetrahedron 'recenttet',  or from one chosen from a random     //
// sample.  The choice is made by determining which one's origin is closest   //
// to the point we are searching for.                                         //
//                                                                            //
//============================================================================//

void TetMeshCore::randomsample(point searchpt,triface *searchtet)
{
  tetrahedron *firsttet, *tetptr;
  point torg;
  void **sampleblock;
  uintptr_t alignptr;
  long sampleblocks, samplesperblock, samplenum;
  long tetblocks, i, j;
  double searchdist, dist;

  if (b->verbose > 2) {
    printf("      Random sampling tetrahedra for searching point %d.\n",
           pointmark(searchpt));
  }

  if (!nonconvex) {
    if (searchtet->tet == NULL) {
      // A null tet. Choose the recenttet as the starting tet.
      *searchtet = recenttet;
    }

    // 'searchtet' should be a valid tetrahedron. Choose the base face
    //   whose vertices must not be 'dummypoint'.
    searchtet->ver = 3;
    // Record the distance from its origin to the searching point.
    torg = org(*searchtet);
    searchdist = (searchpt[0] - torg[0]) * (searchpt[0] - torg[0]) +
                 (searchpt[1] - torg[1]) * (searchpt[1] - torg[1]) +
                 (searchpt[2] - torg[2]) * (searchpt[2] - torg[2]);

    // If a recently encountered tetrahedron has been recorded and has not
    //   been deallocated, test it as a good starting point.
    if (recenttet.tet != searchtet->tet) {
      recenttet.ver = 3;
      torg = org(recenttet);
      dist = (searchpt[0] - torg[0]) * (searchpt[0] - torg[0]) +
             (searchpt[1] - torg[1]) * (searchpt[1] - torg[1]) +
             (searchpt[2] - torg[2]) * (searchpt[2] - torg[2]);
      if (dist < searchdist) {
        *searchtet = recenttet;
        searchdist = dist;
      }
    }
  } else {
    // The mesh is non-convex. Do not use 'recenttet'.
    searchdist = longest;
  }

  // Select "good" candidate using k random samples, taking the closest one.
  //   The number of random samples taken is proportional to the fourth root
  //   of the number of tetrahedra in the mesh. 
  while (samples * samples * samples * samples < tetrahedrons->items) {
    samples++;
  }
  // Find how much blocks in current tet pool.
  tetblocks = (tetrahedrons->maxitems + b->tetrahedraperblock - 1) 
            / b->tetrahedraperblock;
  // Find the average samples per block. Each block at least have 1 sample.
  samplesperblock = 1 + (samples / tetblocks);
  sampleblocks = samples / samplesperblock;
  if (sampleblocks == 0) {
    sampleblocks = 1; // at least one sample block is needed.
  }
  sampleblock = tetrahedrons->firstblock;
  for (i = 0; i < sampleblocks; i++) {
    alignptr = (uintptr_t) (sampleblock + 1);
    firsttet = (tetrahedron *)
               (alignptr + (uintptr_t) tetrahedrons->alignbytes
               - (alignptr % (uintptr_t) tetrahedrons->alignbytes));
    for (j = 0; j < samplesperblock; j++) {
      if (i == tetblocks - 1) {
        // This is the last block.
        samplenum = randomnation((int)
                      (tetrahedrons->maxitems - (i * b->tetrahedraperblock)));
      } else {
        samplenum = randomnation(b->tetrahedraperblock);
      }
      tetptr = (tetrahedron *)
               (firsttet + (samplenum * tetrahedrons->itemwords));
      torg = (point) tetptr[4];
      if (torg != (point) NULL) {
        dist = (searchpt[0] - torg[0]) * (searchpt[0] - torg[0]) +
               (searchpt[1] - torg[1]) * (searchpt[1] - torg[1]) +
               (searchpt[2] - torg[2]) * (searchpt[2] - torg[2]);
        if (dist < searchdist) {
          searchtet->tet = tetptr;
          searchtet->ver = 11; // torg = org(t);
          searchdist = dist;
        }
      } else {
        // A dead tet. Re-sample it.
        if (i != tetblocks - 1) j--;
      }
    }
    sampleblock = (void **) *sampleblock;
  }
}

//============================================================================//
//                                                                            //
// locate()    Find a tetrahedron containing a given point.                   //
//                                                                            //
// Begins its search from 'searchtet', assume there is a line segment L from  //
// a vertex of 'searchtet' to the query point 'searchpt', and simply walk     //
// towards 'searchpt' by traversing all faces intersected by L.               //
//                                                                            //
// On completion, 'searchtet' is a tetrahedron that contains 'searchpt'. The  //
// returned value indicates one of the following cases:                       //
//   - ONVERTEX, the search point lies on the origin of 'searchtet'.          //
//   - ONEDGE, the search point lies on an edge of 'searchtet'.               //
//   - ONFACE, the search point lies on a face of 'searchtet'.                //
//   - INTET, the search point lies in the interior of 'searchtet'.           //
//   - OUTSIDE, the search point lies outside the mesh. 'searchtet' is a      //
//              hull face which is visible by the search point.               //
//                                                                            //
// WARNING: This routine is designed for convex triangulations, and will not  //
// generally work after the holes and concavities have been carved.           //
//                                                                            //
//============================================================================//

enum TetMeshCore::locateresult
  TetMeshCore::locate_dt(point searchpt, triface* searchtet)
{
  //enum {ORGMOVE, DESTMOVE, APEXMOVE} nextmove;
  double ori, oriorg, oridest, oriapex;
  enum locateresult loc = OUTSIDE;
  point toppo;
  int s, i;

  if (searchtet->tet == NULL) {
    searchtet->tet = recenttet.tet;
  }

  if (ishulltet(*searchtet)) {
    // Get its adjacent tet (inside the hull).
    searchtet->tet = decode_tet_only(searchtet->tet[3]);
  }

  // Let searchtet be the face such that 'searchpt' lies above to it.
  for (searchtet->ver = 0; searchtet->ver < 4; searchtet->ver++) {
    ori = orient3d(org(*searchtet), dest(*searchtet), apex(*searchtet), searchpt);
    if (ori < 0.0) break;
  }

  if (searchtet->ver == 4) {
    terminate_tet_core(this, 2);
  }

  // Walk through tetrahedra to locate the point.
  do {

    toppo = oppo(*searchtet);
    
    // Check if the vertex is we seek.
    if (toppo == searchpt) {
      // Adjust the origin of searchtet to be searchpt.
      esymself(*searchtet);
      eprevself(*searchtet);
      loc = ONVERTEX; // return ONVERTEX;
      break;
    }
    
    // We enter from one of serarchtet's faces, which face do we exit?
    // Randomly choose one of three faces (containig  toppo) of this tet.
    s = rand() % 3; // s \in \{0,1,2\}
    for (i = 0; i < s; i++) enextself(*searchtet);

    oriorg = orient3d(dest(*searchtet), apex(*searchtet), toppo, searchpt);
    if (oriorg < 0) {
      //nextmove = ORGMOVE;
      enextesymself(*searchtet);
    } else {
      oridest = orient3d(apex(*searchtet), org(*searchtet), toppo, searchpt);
      if (oridest < 0) {
        //nextmove = DESTMOVE;
        eprevesymself(*searchtet);
      } else {
        oriapex = orient3d(org(*searchtet), dest(*searchtet), toppo, searchpt);
        if (oriapex < 0) {
          //nextmove = APEXMOVE;
          esymself(*searchtet);
        } else {
          // oriorg >= 0, oridest >= 0, oriapex >= 0 ==> found the point.
          // The point we seek must be on the boundary of or inside this
          //   tetrahedron. Check for boundary cases first.
          if (oriorg == 0) {
            // Go to the face opposite to origin.
            enextesymself(*searchtet);
            if (oridest == 0) {
              eprevself(*searchtet); // edge oppo->apex
              if (oriapex == 0) {
                // oppo is duplicated with p.
                loc = ONVERTEX; // return ONVERTEX;
                break;
              }
              loc = ONEDGE; // return ONEDGE;
              break;
            }
            if (oriapex == 0) {
              enextself(*searchtet); // edge dest->oppo
              loc = ONEDGE; // return ONEDGE;
              break;
            }
            loc = ONFACE; // return ONFACE;
            break;
          }
          if (oridest == 0) {
            // Go to the face opposite to destination.
            eprevesymself(*searchtet);
            if (oriapex == 0) {
              eprevself(*searchtet); // edge oppo->org
              loc = ONEDGE; // return ONEDGE;
              break;
            }
            loc = ONFACE; // return ONFACE;
            break;
          }
          if (oriapex == 0) {
            // Go to the face opposite to apex
            esymself(*searchtet);
            loc = ONFACE; // return ONFACE;
            break;
          }
          loc = INTETRAHEDRON;
          break;
        }
      }
    } // if (locateflag)

    // Move to the next tet adjacent to the selected face.
    decode(searchtet->tet[searchtet->ver & 3], *searchtet); // fsymself

    if (ishulltet(*searchtet)) {
      loc = OUTSIDE; // return OUTSIDE;
      break;
    }

  } while (true);

  return loc;
}

enum TetMeshCore::locateresult 
  TetMeshCore::locate(point searchpt, triface* searchtet, int chkencflag)
{
  point torg, tdest, tapex, toppo;
  enum {ORGMOVE, DESTMOVE, APEXMOVE} nextmove;
  double ori, oriorg, oridest, oriapex;
  enum locateresult loc = OUTSIDE;
  //int t1ver;
  int s;

  torg = tdest = tapex = toppo = NULL;

  if (searchtet->tet == NULL) {
    // A null tet. Choose the recenttet as the starting tet.
    searchtet->tet = recenttet.tet;
  }

  // Check if we are in the outside of the convex hull.
  if (ishulltet(*searchtet)) {
    // Get its adjacent tet (inside the hull).
    searchtet->tet = decode_tet_only(searchtet->tet[3]);
  }

  // Let searchtet be the face such that 'searchpt' lies above to it.
  for (searchtet->ver = 0; searchtet->ver < 4; searchtet->ver++) {
    torg = org(*searchtet);
    tdest = dest(*searchtet);
    tapex = apex(*searchtet);
    ori = orient3d(torg, tdest, tapex, searchpt); 
    if (ori < 0.0) break;
  }
  if (searchtet->ver == 4) {
    terminate_tet_core(this, 2);
  }

  // Walk through tetrahedra to locate the point.
  while (true) {
    toppo = oppo(*searchtet);
    
    // Check if the vertex is we seek.
    if (toppo == searchpt) {
      // Adjust the origin of searchtet to be searchpt.
      esymself(*searchtet);
      eprevself(*searchtet);
      loc = ONVERTEX; // return ONVERTEX;
      break;
    }

    // We enter from one of serarchtet's faces, which face do we exit?
    oriorg = orient3d(tdest, tapex, toppo, searchpt); 
    oridest = orient3d(tapex, torg, toppo, searchpt);
    oriapex = orient3d(torg, tdest, toppo, searchpt);

    // Now decide which face to move. It is possible there are more than one
    //   faces are viable moves. If so, randomly choose one.
    if (oriorg < 0) {
      if (oridest < 0) {
        if (oriapex < 0) {
          // All three faces are possible.
          s = randomnation(3); // 's' is in {0,1,2}.
          if (s == 0) {
            nextmove = ORGMOVE;
          } else if (s == 1) {
            nextmove = DESTMOVE;
          } else {
            nextmove = APEXMOVE;
          }
        } else {
          // Two faces, opposite to origin and destination, are viable.
          //s = randomnation(2); // 's' is in {0,1}.
          if (randomnation(2)) {
            nextmove = ORGMOVE;
          } else {
            nextmove = DESTMOVE;
          }
        }
      } else {
        if (oriapex < 0) {
          // Two faces, opposite to origin and apex, are viable.
          //s = randomnation(2); // 's' is in {0,1}.
          if (randomnation(2)) {
            nextmove = ORGMOVE;
          } else {
            nextmove = APEXMOVE;
          }
        } else {
          // Only the face opposite to origin is viable.
          nextmove = ORGMOVE;
        }
      }
    } else {
      if (oridest < 0) {
        if (oriapex < 0) {
          // Two faces, opposite to destination and apex, are viable.
          //s = randomnation(2); // 's' is in {0,1}.
          if (randomnation(2)) {
            nextmove = DESTMOVE;
          } else {
            nextmove = APEXMOVE;
          }
        } else {
          // Only the face opposite to destination is viable.
          nextmove = DESTMOVE;
        }
      } else {
        if (oriapex < 0) {
          // Only the face opposite to apex is viable.
          nextmove = APEXMOVE;
        } else {
          // The point we seek must be on the boundary of or inside this
          //   tetrahedron. Check for boundary cases.
          if (oriorg == 0) {
            // Go to the face opposite to origin.
            enextesymself(*searchtet);
            if (oridest == 0) {
              eprevself(*searchtet); // edge oppo->apex
              if (oriapex == 0) {
                // oppo is duplicated with p.
                loc = ONVERTEX; // return ONVERTEX;
                break;
              }
              loc = ONEDGE; // return ONEDGE;
              break;
            }
            if (oriapex == 0) {
              enextself(*searchtet); // edge dest->oppo
              loc = ONEDGE; // return ONEDGE;
              break;
            }
            loc = ONFACE; // return ONFACE;
            break;
          }
          if (oridest == 0) {
            // Go to the face opposite to destination.
            eprevesymself(*searchtet);
            if (oriapex == 0) {
              eprevself(*searchtet); // edge oppo->org
              loc = ONEDGE; // return ONEDGE;
              break;
            }
            loc = ONFACE; // return ONFACE;
            break;
          }
          if (oriapex == 0) {
            // Go to the face opposite to apex
            esymself(*searchtet);
            loc = ONFACE; // return ONFACE;
            break;
          }
          loc = INTETRAHEDRON; // return INTETRAHEDRON;
          break;
        }
      }
    }
    
    // Move to the selected face.
    if (nextmove == ORGMOVE) {
      enextesymself(*searchtet);
    } else if (nextmove == DESTMOVE) {
      eprevesymself(*searchtet);
    } else {
      esymself(*searchtet);
    }
    if (chkencflag) {
      // Check if we are walking across a subface.
      if (issubface(*searchtet)) {
        loc = ENCSUBFACE;
        break;
      }
    }
    // Move to the adjacent tetrahedron (maybe a hull tetrahedron).
    decode(searchtet->tet[searchtet->ver & 3], *searchtet); // fsymself
    if (ishulltet(*searchtet)) {
      loc = OUTSIDE; // return OUTSIDE;
      break;
    }

    // Retreat the three vertices of the base face.
    torg = org(*searchtet);
    tdest = dest(*searchtet);
    tapex = apex(*searchtet);

  } // while (true)

  return loc;
}

//============================================================================//
//                                                                            //
// insert_vertex_bw()    Insert a vertex using the Bowyer-Watson algorithm.   //
//                                                                            //
// This function is only used for initial Delaunay triangulation construction.//
// It improves the speed of incremental algorithm.                            //
//                                                                            //
//============================================================================//

int  TetMeshCore::insert_vertex_bw(point insertpt, triface *searchtet,
  insertvertexflags *ivf)
{
  tetrahedron **ptptr, *tptr;
  triface cavetet, spintet, neightet, neineitet, *parytet;
  triface oldtet, newtet; //, newneitet;
  point *pts; //, pa, pb, pc, *parypt;
  enum locateresult loc = OUTSIDE;
  double sign, ori;
  //double attrib, volume;
  bool enqflag;
  int t1ver;
  int i, j, k; //, s;

  if (b->verbose > 2) {
    printf("      Insert point %d\n", pointmark(insertpt));
  }

  // Locate the point.
  if (searchtet->tet != NULL) {
    loc = (enum locateresult) ivf->iloc;
  }

  if (loc == OUTSIDE) {
    if (searchtet->tet == NULL) {
      if (!b->weighted) {
        randomsample(insertpt, searchtet);
      } else {
        // Weighted DT. There may exist dangling vertex. 
        *searchtet = recenttet;
      }
    }
    loc = locate_dt(insertpt, searchtet);
  }

  ivf->iloc = (int) loc; // The return value.

  if (b->weighted) {
    if (loc != OUTSIDE) {
      // Check if this vertex is regular.
      pts = (point *) searchtet->tet;
      sign = orient4d_s(pts[4], pts[5], pts[6], pts[7], insertpt,
                        pts[4][3], pts[5][3], pts[6][3], pts[7][3],
                        insertpt[3]);
      if (sign > 0) {
        // This new vertex lies above the lower hull. Do not insert it.
        ivf->iloc = (int) NONREGULAR;
        return 0;
      }
    }
  }

  // Create the initial cavity C(p) which contains all tetrahedra that
  //   intersect p. It may include 1, 2, or n tetrahedra.

  if (loc == OUTSIDE) {
    infect(*searchtet);
    cave_oldtet_list->newindex((void **) &ptptr);
    *ptptr = searchtet->tet;
  } else if (loc == INTETRAHEDRON) {
    infect(*searchtet);
    cave_oldtet_list->newindex((void **) &ptptr);
    *ptptr = searchtet->tet;
  } else if (loc == ONFACE) {
    infect(*searchtet);
    cave_oldtet_list->newindex((void **) &ptptr);
    *ptptr = searchtet->tet;
    neightet.tet = decode_tet_only(searchtet->tet[searchtet->ver & 3]);
    infect(neightet);
    cave_oldtet_list->newindex((void **) &ptptr);
    *ptptr = neightet.tet;
  } else if (loc == ONEDGE) {

    // Add all adjacent boundary tets into list.
    spintet = *searchtet;
    while (1) {
      infect(spintet);
      cave_oldtet_list->newindex((void **) &ptptr);
      *ptptr = spintet.tet;
      fnextself(spintet);
      if (spintet.tet == searchtet->tet) break;
    } // while (1)
  } else if (loc == ONVERTEX) {
    // The point already exist. Do nothing and return.
    return 0;
  } 

  // Create the cavity C(p).

  for (i = 0; i < cave_oldtet_list->objects; i++) {
    ptptr = (tetrahedron **) fastlookup(cave_oldtet_list, i);
    cavetet.tet = *ptptr;
    for (cavetet.ver = 0; cavetet.ver < 4; cavetet.ver++) {
      neightet.tet = decode_tet_only(cavetet.tet[cavetet.ver]);
      if (!infected(neightet)) {
        // neightet.tet is current outside the cavity.
        enqflag = false;
        if (!marktested(neightet)) {
          if (!ishulltet(neightet)) {
            pts = (point *) neightet.tet;
            sign = insphere_s(pts[4], pts[5], pts[6], pts[7], insertpt);
            enqflag = (sign < 0.0);
          } else {
            pts = (point *) neightet.tet;
            ori = orient3d(pts[4], pts[5], pts[6], insertpt);
            if (ori < 0) {
              // A visible hull face.
              enqflag = true;
            } else if (ori == 0.) {
              // A coplanar hull face. We need to test if this hull face is
              //   Delaunay or not. We test if the adjacent tet (not faked)
              //   of this hull face is Delaunay or not.
              triface neineitet;
              neineitet.tet = decode_tet_only(neightet.tet[3]);
              pts = (point *) neineitet.tet;
              sign = insphere_s(pts[4],pts[5],pts[6],pts[7], insertpt);
              enqflag = (sign < 0.0);
            }
          }
          marktest(neightet);
        }
        if (enqflag) {
          infect(neightet);
          cave_oldtet_list->newindex((void **) &ptptr);
          *ptptr = neightet.tet;
        } else {
          // A boundary face.
          cavebdrylist->newindex((void **) &parytet);
          *parytet = cavetet;
        }
      } // if (!infected(neightet))
    }
  } // i

  // Create new tetrahedra to fill the cavity.
  int f_out = cavebdrylist->objects;
  int v_out = (f_out + 4) / 2;  


  triface *pcavetet;
  point V[3];
  int local_vcount = 0; // local index of vertex
  int sidx[3];

  static int row_v08_tbl[12] = {8,9,10,11,0,1,2,3,4,5,6,7};
  static int row_v11_tbl[12] = {8,9,10,11,0,1,2,3,4,5,6,7};
  static int col_v01_tbl[12] = {1,1,1,1,5,5,5,5,9,9,9,9};
  static int col_v02_tbl[12] = {2,2,2,2,6,6,6,6,10,10,10,10};
  static int col_v08_tbl[12] = {8,8,8,8,0,0,0,0,4,4,4,4};
  static int col_v11_tbl[12] = {11,11,11,11,3,3,3,3,7,7,7,7};

  triface *tmp_bw_faces = NULL;
  int shiftbits = 0;

  if (v_out < 64) {
    shiftbits = 6;
    tmp_bw_faces = _bw_faces;
  } else if (v_out < 1024) {
    // Dynamically allocate an array to store the adjacencies.
    int arysize = 1;
    int tmp = v_out;
    shiftbits = 1;
    while ((tmp >>= 1)) shiftbits++;
    arysize <<= shiftbits;
    tmp_bw_faces = new triface[arysize * arysize];
  } 

  if (v_out < 1024) {
    for (i = 0; i < f_out; i++) {
      pcavetet = (triface *) fastlookup(cavebdrylist, i);
      oldtet = *pcavetet;

      // Get the tet outside the cavity.
      decode(oldtet.tet[oldtet.ver], neightet);
      unmarktest(neightet);

      if (ishulltet(oldtet)) {
        // neightet.tet may be also a hull tet (=> oldtet is a hull edge).
        neightet.ver = epivot[neightet.ver];
        if ((apex(neightet) == dummypoint)) {
          hullsize++; // Create a new hull tet.
        }
      }
 
      // Create a new tet in the cavity.
      V[0] = dest(neightet);
      V[1] = org(neightet);
      V[2] = apex(neightet);
      maketetrahedron2(&newtet, V[1], V[0], insertpt, V[2]);
      //bond(newtet, neightet);
      newtet.tet[2] = encode2(neightet.tet, neightet.ver);
      neightet.tet[neightet.ver & 3] = encode2(newtet.tet, col_v02_tbl[neightet.ver]);

      // Fill the adjacency matrix, and count v_out.
      for (j = 0; j < 3; j++) {
        tptr = (tetrahedron *) point2tet(V[j]);
        if (((point *) tptr)[6] != insertpt) {
          // Found a unique vertex of the cavity.
          setpointgeomtag(V[j], local_vcount++);
          //local_vcount++;
          setpoint2tet(V[j], (tetrahedron) (newtet.tet));
        }
        sidx[j] = pointgeomtag(V[j]);
      } // j
      
      neightet.tet = newtet.tet;
      // Avoid using lookup tables.
      neightet.ver = 11;
      tmp_bw_faces[(sidx[1] << shiftbits) | sidx[0]] = neightet;
      neightet.ver = 1;
      tmp_bw_faces[(sidx[2] << shiftbits) | sidx[1]] = neightet;
      neightet.ver = 8;
      tmp_bw_faces[(sidx[0] << shiftbits) | sidx[2]] = neightet;
      
      *pcavetet = newtet;
    } // i // f_out

    // Set a handle for speeding point location.
    // Randomly pick a new tet.
    i = rand() % f_out;
    recenttet = * (triface *) fastlookup(cavebdrylist, i);    
    setpoint2tet(insertpt, (tetrahedron) (recenttet.tet));

    for (i = 0; i < f_out; i++) {
      neightet = * (triface *) fastlookup(cavebdrylist, i);
      if (neightet.tet[3] == NULL) {
        neightet.ver = 11;
        j = pointgeomtag(org(neightet));
        k = pointgeomtag(dest(neightet));
        neineitet = tmp_bw_faces[(k << shiftbits) | j];
        // bondtbl[i][j] = (j & 3) + (((i & 12) + (j & 12)) % 12);
        neightet.tet[3] = encode2(neineitet.tet, row_v11_tbl[neineitet.ver]);
        neineitet.tet[neineitet.ver & 3] = encode2(neightet.tet, col_v11_tbl[neineitet.ver]);
      }
      if (neightet.tet[1] == NULL) {
        neightet.ver = 1;
        j = pointgeomtag(org(neightet));
        k = pointgeomtag(dest(neightet));
        neineitet = tmp_bw_faces[(k << shiftbits) | j];
        neightet.tet[1] = encode2(neineitet.tet, neineitet.ver); // row_v01_tbl
        neineitet.tet[neineitet.ver & 3] = encode2(neightet.tet, col_v01_tbl[neineitet.ver]);
      }
      if (neightet.tet[0] == NULL) {
        neightet.ver = 8;
        j = pointgeomtag(org(neightet));
        k = pointgeomtag(dest(neightet));
        neineitet = tmp_bw_faces[(k << shiftbits) | j];
        // bondtbl[i][j] = (j & 3) + (((i & 12) + (j & 12)) % 12);
        neightet.tet[0] = encode2(neineitet.tet, row_v08_tbl[neineitet.ver]);
        neineitet.tet[neineitet.ver & 3] = encode2(neightet.tet, col_v08_tbl[neineitet.ver]);
      }    
    } // i
    
    if (v_out >= 64) {
      delete [] tmp_bw_faces;
    }
  } // v_out < 1024
  else {
    // Fill a very large cavity with original neighboring searching method.
    for (i = 0; i < f_out; i++) {
      pcavetet = (triface *) fastlookup(cavebdrylist, i);
      oldtet = *pcavetet;

      // Get the tet outside the cavity.
      decode(oldtet.tet[oldtet.ver], neightet);
      unmarktest(neightet);

      if (ishulltet(oldtet)) {
        // neightet.tet may be also a hull tet (=> oldtet is a hull edge).
        neightet.ver = epivot[neightet.ver];
        if ((apex(neightet) == dummypoint)) {
          hullsize++; // Create a new hull tet.
        }
      }
 
      // Create a new tet in the cavity.
      V[0] = dest(neightet);
      V[1] = org(neightet);
      V[2] = apex(neightet);
      maketetrahedron2(&newtet, V[1], V[0], insertpt, V[2]);
      //newtet.ver = 2; // esymself(newtet);
      //assert(oppo(newtet) == insertpt);
 
      //bond(newtet, neightet);
      newtet.tet[2] = encode2(neightet.tet, neightet.ver);
      neightet.tet[neightet.ver & 3] = encode2(newtet.tet, col_v02_tbl[neightet.ver]);

      // Fill the adjacency matrix, and count v_out.
      for (j = 0; j < 3; j++) {
        tptr = (tetrahedron *) point2tet(V[j]);
        if (((point *) tptr)[6] != insertpt) {
          // Found a unique vertex of the cavity.
          //setpointgeomtag(V[j], local_vcount);
          local_vcount++;
          setpoint2tet(V[j], (tetrahedron) (newtet.tet));
        }
        //sidx[j] = pointgeomtag(V[j]);
      } // j
    } // i, f_out
    
    // Set a handle for speeding point location.
    //recenttet = newtet;
    //setpoint2tet(insertpt, (tetrahedron) (newtet.tet));
    i = rand() % f_out;
    recenttet = * (triface *) fastlookup(cavebdrylist, i);
    // This is still an oldtet.
    fsymself(recenttet);
    fsymself(recenttet);
    setpoint2tet(insertpt, (tetrahedron) (recenttet.tet));

    for (i = 0; i < f_out; i++) {
      pcavetet = (triface *) fastlookup(cavebdrylist, i);
      oldtet = *pcavetet;

      fsym(oldtet, neightet);
      fsym(neightet, newtet);
      // Comment: oldtet and newtet must be at the same directed edge.
      // Connect the three other faces of this newtet.
      for (j = 0; j < 3; j++) {
        esym(newtet, neightet); // Go to the face.
        if (neightet.tet[neightet.ver & 3] == NULL) {
          // Find the adjacent face of this newtet.
          spintet = oldtet;
          while (1) {
            fnextself(spintet);
            if (!infected(spintet)) break;
          }
          fsym(spintet, neineitet);
          esymself(neineitet);
          bond(neightet, neineitet);
        }
        enextself(newtet);
        enextself(oldtet);
      } // j
    } // i
  } // fill cavity
  
  // C(p) is re-meshed successfully.

  // Delete the old tets in C(p).
  for (i = 0; i < cave_oldtet_list->objects; i++) {
    oldtet.tet = *(tetrahedron **) fastlookup(cave_oldtet_list, i);
    if (ishulltet(oldtet)) {
      hullsize--;
    }
    tetrahedrondealloc(oldtet.tet);
  }

  cave_oldtet_list->restart();
  cavebdrylist->restart();

  return 1;
}

//============================================================================//
//                                                                            //
// initialdelaunay()    Create an initial Delaunay tetrahedralization.        //
//                                                                            //
// The tetrahedralization contains only one tetrahedron abcd, and four hull   //
// tetrahedra. The points pa, pb, pc, and pd must be linearly independent.    //
//                                                                            //
//============================================================================//

void TetMeshCore::initialdelaunay(point pa, point pb, point pc, point pd)
{
  triface firsttet, tetopa, tetopb, tetopc, tetopd;
  triface worktet, worktet1;

  if (b->verbose > 2) {
    printf("      Create init tet (%d, %d, %d, %d)\n", pointmark(pa),
           pointmark(pb), pointmark(pc), pointmark(pd));
  }

  // Create the first tetrahedron.
  maketetrahedron2(&firsttet, pa, pb, pc, pd);
  //setvertices(firsttet, pa, pb, pc, pd);

  // Create four hull tetrahedra.
  maketetrahedron2(&tetopa, pb, pc, pd, dummypoint);
  //setvertices(tetopa, pb, pc, pd, dummypoint);
  maketetrahedron2(&tetopb, pc, pa, pd, dummypoint);
  //setvertices(tetopb, pc, pa, pd, dummypoint);
  maketetrahedron2(&tetopc, pa, pb, pd, dummypoint);
  //setvertices(tetopc, pa, pb, pd, dummypoint);
  maketetrahedron2(&tetopd, pb, pa, pc, dummypoint);
  //setvertices(tetopd, pb, pa, pc, dummypoint);

  hullsize += 4;

  // Connect hull tetrahedra to firsttet (at four faces of firsttet).
  bond(firsttet, tetopd);
  esym(firsttet, worktet);
  bond(worktet, tetopc); // ab
  enextesym(firsttet, worktet);
  bond(worktet, tetopa); // bc 
  eprevesym(firsttet, worktet);
  bond(worktet, tetopb); // ca

  // Connect hull tetrahedra together (at six edges of firsttet).
  esym(tetopc, worktet); 
  esym(tetopd, worktet1);
  bond(worktet, worktet1); // ab
  esym(tetopa, worktet);
  eprevesym(tetopd, worktet1);
  bond(worktet, worktet1); // bc
  esym(tetopb, worktet);
  enextesym(tetopd, worktet1);
  bond(worktet, worktet1); // ca
  eprevesym(tetopc, worktet);
  enextesym(tetopb, worktet1);
  bond(worktet, worktet1); // da
  eprevesym(tetopa, worktet);
  enextesym(tetopc, worktet1);
  bond(worktet, worktet1); // db
  eprevesym(tetopb, worktet);
  enextesym(tetopa, worktet1);
  bond(worktet, worktet1); // dc

  // Set the vertex type.
  if (pointtype(pa) == UNUSEDVERTEX) {
    setpointtype(pa, VOLVERTEX);
  }
  if (pointtype(pb) == UNUSEDVERTEX) {
    setpointtype(pb, VOLVERTEX);
  }
  if (pointtype(pc) == UNUSEDVERTEX) {
    setpointtype(pc, VOLVERTEX);
  }
  if (pointtype(pd) == UNUSEDVERTEX) {
    setpointtype(pd, VOLVERTEX);
  }

  setpoint2tet(pa, encode(firsttet));
  setpoint2tet(pb, encode(firsttet));
  setpoint2tet(pc, encode(firsttet));
  setpoint2tet(pd, encode(firsttet));

  setpoint2tet(dummypoint, encode(tetopa));

  // Remember the first tetrahedron.
  recenttet = firsttet;
}


//============================================================================//
//                                                                            //
// incrementaldelaunay()    Create a Delaunay tetrahedralization by           //
//                          the incremental approach.                         //
//                                                                            //
//============================================================================//


void TetMeshCore::incrementaldelaunay(clock_t& tv)
{
  triface searchtet;
  point *permutarray, swapvertex;
  double v1[3], v2[3], n[3];
  double bboxsize, bboxsize2, bboxsize3, ori;
  int randindex; 
  int ngroup = 0;
  int i, j;

  if (!b->quiet) {
    printf("Delaunizing vertices...\n");
  }
  // Form a random permuation (uniformly at random) of the set of vertices.
  permutarray = new point[in->numberofpoints];
  points->traversalinit();

  if (b->no_sort) {
    if (b->verbose) {
      printf("  Using the input order.\n"); 
    }
    for (i = 0; i < in->numberofpoints; i++) {
      permutarray[i] = (point) points->traverse();
    }
  } else {
    if (b->verbose) {
      printf("  Permuting vertices.\n"); 
    }
    srand(in->numberofpoints);
    for (i = 0; i < in->numberofpoints; i++) {
      randindex = rand() % (i + 1); // randomnation(i + 1);
      permutarray[i] = permutarray[randindex];
      permutarray[randindex] = (point) points->traverse();
    }
    if (b->brio_hilbert) { // -b option
      if (b->verbose) {
        printf("  Sorting vertices.\n"); 
      }
      hilbert_init(in->mesh_dim);
      brio_multiscale_sort(permutarray, in->numberofpoints, b->brio_threshold, 
                           b->brio_ratio, &ngroup);
    }
  }

  tv = clock(); // Remember the time for sorting points.

  // Calculate the diagonal size of its bounding box.
  bboxsize = sqrt(norm2(xmax - xmin, ymax - ymin, zmax - zmin));
  bboxsize2 = bboxsize * bboxsize;
  bboxsize3 = bboxsize2 * bboxsize;

  // Make sure the second vertex is not identical with the first one.
  i = 1;
  while ((distance(permutarray[0],permutarray[i])/bboxsize)<b->epsilon) {
    i++;
    if (i == in->numberofpoints - 1) {
      printf("Exception:  All vertices are (nearly) identical (Tol = %g).\n",
             b->epsilon);
      terminate_tet_core(this, 10);
    }
  }
  if (i > 1) {
    // Swap to move the non-identical vertex from index i to index 1.
    swapvertex = permutarray[i];
    permutarray[i] = permutarray[1];
    permutarray[1] = swapvertex;
  }

  // Make sure the third vertex is not collinear with the first two.
  i = 2;
  for (j = 0; j < 3; j++) {
    v1[j] = permutarray[1][j] - permutarray[0][j];
    v2[j] = permutarray[i][j] - permutarray[0][j];
  }
  cross(v1, v2, n);
  while ((sqrt(norm2(n[0], n[1], n[2])) / bboxsize2) < b->epsilon) {
    i++;
    if (i == in->numberofpoints - 1) {
      printf("Exception:  All vertices are (nearly) collinear (Tol = %g).\n",
             b->epsilon);
      terminate_tet_core(this, 10);
    }
    for (j = 0; j < 3; j++) {
      v2[j] = permutarray[i][j] - permutarray[0][j];
    }
    cross(v1, v2, n);
  }
  if (i > 2) {
    // Swap to move the non-identical vertex from index i to index 1.
    swapvertex = permutarray[i];
    permutarray[i] = permutarray[2];
    permutarray[2] = swapvertex;
  }

  // Make sure the fourth vertex is not coplanar with the first three.
  i = 3;
  ori = orient3dfast(permutarray[0], permutarray[1], permutarray[2], 
                     permutarray[i]);
  while ((fabs(ori) / bboxsize3) < b->epsilon) {
    i++;
    if (i == in->numberofpoints) {
      printf("Exception:  All vertices are coplanar (Tol = %g).\n",
             b->epsilon);
      terminate_tet_core(this, 10);
    }
    ori = orient3dfast(permutarray[0], permutarray[1], permutarray[2], 
                       permutarray[i]);
  }
  if (i > 3) {
    // Swap to move the non-identical vertex from index i to index 1.
    swapvertex = permutarray[i];
    permutarray[i] = permutarray[3];
    permutarray[3] = swapvertex;
  }

  // Orient the first four vertices in permutarray so that they follow the
  //   right-hand rule.
  if (ori > 0.0) {
    // Swap the first two vertices.
    swapvertex = permutarray[0];
    permutarray[0] = permutarray[1];
    permutarray[1] = swapvertex;
  }

  // Create the initial Delaunay tetrahedralization.
  initialdelaunay(permutarray[0], permutarray[1], permutarray[2],
                  permutarray[3]);

  if (b->verbose) {
    printf("  Incrementally inserting vertices.\n");
  }
  insertvertexflags ivf;
  flipconstraints fc;

  ivf.bowywat = 1; // Use Bowyer-Watson algorithm
  ivf.lawson = 0;


  for (i = 4; i < in->numberofpoints; i++) {
    if (pointtype(permutarray[i]) == UNUSEDVERTEX) {
      setpointtype(permutarray[i], VOLVERTEX);
    }
    if (b->brio_hilbert || b->no_sort) { // -b or -b/1
      // Start the last updated tet.
      searchtet.tet = recenttet.tet;
    } else { // -b0
      // Randomly choose the starting tet for point location.
      searchtet.tet = NULL;
    }
    ivf.iloc = (int) OUTSIDE;
    // Insert the vertex.
    if (!insert_vertex_bw(permutarray[i], &searchtet, &ivf)) {
      if (ivf.iloc == (int) ONVERTEX) {
        // The point already exists. Mark it and do nothing on it.
        swapvertex = org(searchtet);
        if (b->object != TetMeshBehavior::STL) {
          if (!b->quiet) {
            printf("Warning:  Point #%d is coincident with #%d. Ignored!\n",
                   pointmark(permutarray[i]), pointmark(swapvertex));
          }
        }
        setpoint2ppt(permutarray[i], swapvertex);
        setpointtype(permutarray[i], DUPLICATEDVERTEX);
        dupverts++;
      } else if (ivf.iloc == (int) NEARVERTEX) {
        // This should not happen by insert_point_bw().
        terminate_tet_core(this, 2); // report a bug.
      } else if (ivf.iloc == (int) NONREGULAR) {
        // The point is non-regular. Skipped.
        if (b->verbose) {
          printf("  Point #%d is non-regular, skipped.\n",
                 pointmark(permutarray[i]));
        }
        setpointtype(permutarray[i], NREGULARVERTEX);
        nonregularcount++;
      }
    }
  }


  
  delete [] permutarray;
}

//                                                                            //
//                                                                            //
//== delaunay_cxx ============================================================//

} // namespace sqmesh::mesh::tet::detail
