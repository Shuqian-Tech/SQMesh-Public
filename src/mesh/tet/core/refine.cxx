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

//== refine_cxx ==============================================================//
//                                                                            //
//                                                                            //

//============================================================================//
//                                                                            //
// makesegmentendpointsmap()    Create a map from a segment to its endpoints. //
//                                                                            //
// The map is saved in the array 'segmentendpointslist'. The length of this   //
// array is twice the number of segments.  Each segment is assigned a unique  //
// index (starting from 0).                                                   //
//                                                                            //
//============================================================================//

void TetMeshCore::makesegmentendpointsmap()
{
  arraypool *segptlist;
  face segloop, prevseg, nextseg;
  point eorg, edest, *parypt;
  int segindex = 0, idx = 0;
  int i;

  if (b->verbose > 0) {
    printf("  Creating the segment-endpoints map.\n");
  }
  segptlist = new arraypool(2 * sizeof(point), 10);

  // for creating ridge_vertex-to-segment map;
  // The index might start from 0 or 1.
  idx_segment_ridge_vertex_list = new int[points->items + 2];
  for (i = 0; i < points->items + 2; i++) {
    idx_segment_ridge_vertex_list[i] = 0;
  }

  // A segment s may have been split into many subsegments. Operate the one
  //   which contains the origin of s. Then mark the rest of subsegments.
  subsegs->traversalinit();
  segloop.sh = shellfacetraverse(subsegs);
  segloop.shver = 0;
  while (segloop.sh != NULL) {
    senext2(segloop, prevseg);
    spivotself(prevseg);
    if (prevseg.sh == NULL) {
      eorg = sorg(segloop);
      edest = sdest(segloop);
      setfacetindex(segloop, segindex);
      senext(segloop, nextseg);
      spivotself(nextseg);
      while (nextseg.sh != NULL) {
        setfacetindex(nextseg, segindex);
        nextseg.shver = 0;
        if (sorg(nextseg) != edest) sesymself(nextseg);
        edest = sdest(nextseg);
        // Go the next connected subsegment at edest.
        senextself(nextseg);
        spivotself(nextseg);
      }
      segptlist->newindex((void **) &parypt);
      parypt[0] = eorg;
      parypt[1] = edest;
      segindex++;
      // for creating adj_ridge_vertex_list;
      idx_segment_ridge_vertex_list[pointmark(eorg)]++;
      idx_segment_ridge_vertex_list[pointmark(edest)]++;
    }
    segloop.sh = shellfacetraverse(subsegs);
  }

  if (b->verbose) {
    printf("  Found %ld segments.\n", segptlist->objects);
  }

  segmentendpointslist_length = segptlist->objects;
  segmentendpointslist = new point[segptlist->objects * 2];

  totalworkmemory += (segptlist->objects * 2) * sizeof(point *);

  for (i = 0; i < segptlist->objects; i++) {
    parypt = (point *) fastlookup(segptlist, i);
    segmentendpointslist[idx++] = parypt[0];
    segmentendpointslist[idx++] = parypt[1];
  }

  // Create the adj_ridge_vertex_list.
  int j = idx_segment_ridge_vertex_list[0], k;
  idx_segment_ridge_vertex_list[0] = 0;
  for (i = 0; i < points->items + 1; i++) {
    k = idx_segment_ridge_vertex_list[i+1];
    idx_segment_ridge_vertex_list[i+1] = idx_segment_ridge_vertex_list[i] + j;
    j = k;
  }

  //assert(i == points->items+1);
  int total_count = idx_segment_ridge_vertex_list[i] + 1;
  segment_ridge_vertex_list = new point[total_count];
  for (i = 0; i < segptlist->objects; i++) {
    eorg  = segmentendpointslist[i*2];
    edest = segmentendpointslist[i*2+1];
    j = pointmark(eorg);
    k = pointmark(edest);
    segment_ridge_vertex_list[idx_segment_ridge_vertex_list[j]] = edest; //eorg;
    segment_ridge_vertex_list[idx_segment_ridge_vertex_list[k]] = eorg;  //edest;
    idx_segment_ridge_vertex_list[j]++;
    idx_segment_ridge_vertex_list[k]++;
  }

  // Counters in idx_adj_ridge_vertex_list[] are shifted by 1.
  for (i = points->items; i >= 0; i--) {
    idx_segment_ridge_vertex_list[i+1] = idx_segment_ridge_vertex_list[i];
  }
  idx_segment_ridge_vertex_list[0] = 0;


  delete segptlist;
}

//============================================================================//
//                                                                            //
// set_ridge_vertex_protecting_ball()    Calculate the protecting ball for a  //
//                                       given ridge vertex.                  //
//                                                                            //
//============================================================================//

double TetMeshCore::set_ridge_vertex_protecting_ball(point ridge_pt)
{
  double rv = getpointinsradius(ridge_pt);
  if (rv == 0.) {
    double mindist = 1.e+30, dist;
    int idx = pointmark(ridge_pt);
    for (int i = idx_segment_ridge_vertex_list[idx];
           i < idx_segment_ridge_vertex_list[idx+1]; i++) {
      dist = distance(ridge_pt, segment_ridge_vertex_list[i]);
      if (mindist > dist) mindist = dist;
    }
    rv = mindist * 0.95; // mindist / 3.0; // refer to J. Shewchuk
    setpointinsradius(ridge_pt, rv);
  }
  return rv;
}

//============================================================================//
//                                                                            //
// get_min_diahedral_angle()    Calculate the minimum (interior) dihedral     //
//                              angle a given segment.                        //
//                                                                            //
//============================================================================//

double TetMeshCore::get_min_diahedral_angle(face* seg)
{
  triface adjtet, spintet;
  face startsh, neighsh;
  point pa, pb, pc1, pc2;
  double n1[3], n2[3];
  double n1len, n2len;
  double costheta; //, ori;
  double theta, sum_theta, minang = 2.0 * PI;
  int t1ver;
  
  pa = sorg(*seg);
  pb = sdest(*seg);
  spivot(*seg, startsh);
  if (startsh.sh == NULL) {
    // This segment is not connected by any facet.
    sstpivot1(*seg, adjtet);
    if (adjtet.tet != NULL) {
      // This segment is completely inside the volume.
      return 360.; // 2*pi.
    }
  } else {
    if (sorg(startsh) != pa) sesymself(startsh);
    stpivot(startsh, adjtet);
  }
  if (adjtet.tet == NULL) {
    // This segment is not inserted (recovered) yet.
    return 0.;
  }


  sum_theta = 0.;
  spintet = adjtet;
  while (true) {
    if (!ishulltet(spintet)) {
      // Increase the interior dihedral angle (sum_theta).
      pc1 = apex(spintet);
      pc2 = oppo(spintet);
      facenormal(pa, pb, pc1, n1, 1, NULL);
      facenormal(pa, pb, pc2, n2, 1, NULL);
      n1len = sqrt(dot(n1, n1));
      n2len = sqrt(dot(n2, n2));
      costheta = dot(n1, n2) / (n1len * n2len);
      // Be careful rounding error!
      if (costheta > 1.0) {
      costheta = 1.0;
      } else if (costheta < -1.0) {
      costheta = -1.0;
      }
      theta = acos(costheta);
      sum_theta += theta;
    }
    // Go to the next adjacent tetrahedron at this segment.
    fnextself(spintet);
    // Check if we meet a subface.
    tspivot(spintet, neighsh);
    if ((neighsh.sh != NULL) && (sum_theta > 0.)) {
      // Update the smallest dihedral angle.
      if (sum_theta < minang) minang = sum_theta;
      sum_theta = 0.; // clear it
    }
    if (spintet.tet == adjtet.tet) break;
  }
      
  double mindihedang = minang / PI * 180.;
  return mindihedang;
}

//============================================================================//
//                                                                            //
// get_min_angle_at_ridge_vertex()    Calculate the minimum face angle at a   //
//                                    given ridge vertex.                     //
//                                                                            //
//============================================================================//

double TetMeshCore::get_min_angle_at_ridge_vertex(face* seg)
{
  face startsh, spinsh, neighsh;
  point pa, pb, pc;
  double theta, sum_theta, minang = 2.0 * PI;
  
  pa = sorg(*seg);
  spivot(*seg, startsh);
  if (startsh.sh == NULL) {
    // This segment does not belong to any facet.
    return 360.; // 2*pi.
  } else {
    if (sorg(startsh) != pa) sesymself(startsh);
  }

  spinsh = startsh;
  while (spinsh.sh != NULL) {
    sum_theta = 0.;
    neighsh = spinsh;
    while (true) {
      pb = sdest(neighsh);
      pc = sapex(neighsh);
      theta = interiorangle(pa, pb, pc, NULL);
      sum_theta += theta;
      senext2self(neighsh);
      if (isshsubseg(neighsh)) break;
      spivotself(neighsh);
      if (sorg(neighsh) != pa) sesymself(neighsh);
    }
    if (sum_theta < minang) {
      minang = sum_theta;
    }
    // Go to the next facet at this segment.
    spivotself(spinsh);
    if (spinsh.sh == startsh.sh) break;
    if (spinsh.sh == NULL) break; // A single facet may happen.
    if (sorg(spinsh) != pa) sesymself(spinsh);
  }

  return minang / PI * 180.;
}

//============================================================================//
//                                                                            //
// create_segment_info_list()    Calculate the minimum dihedral angle at a    //
//                               a given segment.                             //
//                                                                            //
// segment_info_list = new double[segmentendpointslist_length * 4];           //
//  - [0] min_dihedral_angle (degree) at this segment,                        //
//  - [1] min_protect_cylinder_radius at this segment (for bookkeeping only), //
//  - [2] min_seg_seg_angle (degree) at its endpoint [0],                     //
//  - [3] min_seg_seg_angle (degree) at its endpoint [1].                     //
//                                                                            //
// This function must be called after makesegmentendpointsmap(). The number   //
//   of unique segments (segmentendpointslist_length) is calculated.          //
//                                                                            //
//============================================================================//

void TetMeshCore::create_segment_info_list()
{
  face min_dihedral_ang_seg;
  point min_face_ang_vertex;
  double min_dihedral_ang = 360.;
  double min_face_ang = 360.;
  
  if (b->verbose > 0) {
    printf("  Creating the segment_info_list.\n");
  }
  if (segment_info_list != NULL) {
    delete [] segment_info_list;
  }

  if (subsegs->items == 0) {
    return; // There is no segments.
  }

  int count = (segmentendpointslist_length + 1) * 4;
  segment_info_list = new double[count];
  for (int i = 0; i < count; i++) {
    segment_info_list[i] = 0.;
  }

  // Loop through the list of segments.
  face segloop;
  subsegs->traversalinit();
  segloop.sh = shellfacetraverse(subsegs);
  while (segloop.sh != NULL) {
    int segidx = getfacetindex(segloop);
    // Check if this segment has been already calulcated.
    double *values = &(segment_info_list[segidx * 4]);

    // The min_diahedral_angle at this segment is in (0, 2pi].
    if (values[0] == 0.) {
      // Get the smallest dihedral angle at this segment.
      values[0] = get_min_diahedral_angle(&segloop);
      if (values[0] < min_dihedral_ang) {
        min_dihedral_ang = values[0];
        min_dihedral_ang_seg = segloop;
      }
    }

    point *endpts = &(segmentendpointslist[segidx * 2]);

    for (int k = 0; k < 2; k++) {
      segloop.shver = 0;
      if (values[2+k] == 0.) {
        if (sorg(segloop) != endpts[k]) {
          sesymself(segloop);
        }
        if (sorg(segloop) == endpts[k]) {
          // Get the min face angle at vertex endpts[0].
          values[2+k] = get_min_angle_at_ridge_vertex(&segloop);
          if (values[2+k] < min_face_ang) {
            min_face_ang = values[2+k];
            min_face_ang_vertex = endpts[k];
          }
        }
      }
    }

    segloop.sh = shellfacetraverse(subsegs);
  }

  if (b->verbose) {
    printf("  min_dihedral angle = %g degree, at segment [%d,%d]\n",
           min_dihedral_ang, pointmark(sorg(min_dihedral_ang_seg)),
           pointmark(sdest(min_dihedral_ang_seg)));
    printf("  min face angle = %g degree, at vertex %d\n",
           min_face_ang, pointmark(min_face_ang_vertex));
  }
  
}

//============================================================================//
//                                                                            //
// makefacetverticesmap()    Create a map from facet to its vertices.         //
//                                                                            //
// All facets will be indexed (starting from 0).  The map is saved in two     //
// global arrays: 'idx2facetlist' and 'facetverticeslist'.                    //
//                                                                            //
//============================================================================//

void TetMeshCore::makefacetverticesmap()
{
  arraypool *facetvertexlist, *vertlist, **paryvertlist;
  face subloop, neighsh, *parysh, *parysh1;
  point pa, *ppt, *parypt;
  verttype vt;
  int facetindex, totalvertices;
  unsigned long max_facet_size = 0l;
  int max_facet_idx = 0;
  int i, j, k;

  if (b->verbose) {
    printf("  Creating the facet vertices map.\n");
  }

  facetvertexlist = new arraypool(sizeof(arraypool *), 10);
  facetindex = totalvertices = 0;

  // The index might start from 0 or 1.
  idx_ridge_vertex_facet_list = new int[points->items + 2];
  for (i = 0; i < points->items + 2; i++) {
    idx_ridge_vertex_facet_list[i] = 0;
  }

  subfaces->traversalinit();
  subloop.sh = shellfacetraverse(subfaces);
  while (subloop.sh != NULL) {
    if (!sinfected(subloop)) {
      // A new facet. Create its vertices list.
      vertlist = new arraypool(sizeof(point *), 8);
      ppt = (point *) &(subloop.sh[3]);
      for (k = 0; k < 3; k++) {
        vt = pointtype(ppt[k]);
        //if ((vt != FREESEGVERTEX) && (vt != FREEFACETVERTEX)) {
        if (vt == RIDGEVERTEX) {
          pinfect(ppt[k]);
          vertlist->newindex((void **) &parypt);
          *parypt = ppt[k];
          // for creating ridge_vertex-to-facet map.
          idx_ridge_vertex_facet_list[pointmark(ppt[k])]++;
        }
      }
      sinfect(subloop);
      caveshlist->newindex((void **) &parysh);
      *parysh = subloop;
      for (i = 0; i < caveshlist->objects; i++) {
        parysh = (face *) fastlookup(caveshlist, i);
        setfacetindex(*parysh, facetindex);
        for (j = 0; j < 3; j++) {
          if (!isshsubseg(*parysh)) {
            spivot(*parysh, neighsh);
            if (!sinfected(neighsh)) {
              pa = sapex(neighsh);
              if (!pinfected(pa)) {
                vt = pointtype(pa);
                //if ((vt != FREESEGVERTEX) && (vt != FREEFACETVERTEX)) {
                if (vt == RIDGEVERTEX) {
                  pinfect(pa);
                  vertlist->newindex((void **) &parypt);
                  *parypt = pa;
                  // for creating ridge_vertex-to-facet map.
                  idx_ridge_vertex_facet_list[pointmark(pa)]++;
                }
              }
              sinfect(neighsh);
              caveshlist->newindex((void **) &parysh1);
              *parysh1 = neighsh;
            }
          }
          senextself(*parysh);
        }
      } // i
      totalvertices += (int) vertlist->objects;
      if (max_facet_size < vertlist->objects) {
        max_facet_size = vertlist->objects;
        max_facet_idx = facetindex;
      }
      // Uninfect facet vertices.
      for (k = 0; k < vertlist->objects; k++) {
        parypt = (point *) fastlookup(vertlist, k);
        puninfect(*parypt);
      }
      caveshlist->restart();
      // Save this vertex list.
      facetvertexlist->newindex((void **) &paryvertlist);
      *paryvertlist = vertlist;
      facetindex++;
    } 
    subloop.sh = shellfacetraverse(subfaces);
  }

  // All subfaces are infected. Uninfect them.
  subfaces->traversalinit();
  subloop.sh = shellfacetraverse(subfaces);
  while (subloop.sh != NULL) {
    suninfect(subloop);
    subloop.sh = shellfacetraverse(subfaces);
  }

  if (b->verbose) {
    printf("  Found %ld facets. Max facet idx(%d), size(%ld)\n",
           facetvertexlist->objects, max_facet_idx, max_facet_size);
  }

  number_of_facets = facetindex;
  idx2facetlist = new int[facetindex + 1];
  facetverticeslist = new point[totalvertices];

  // create ridge_vertex-to-facet map.
  j = idx_ridge_vertex_facet_list[0]; //k;
  idx_ridge_vertex_facet_list[0] = 0;
  for (i = 0; i < points->items + 1; i++) {
    k = idx_ridge_vertex_facet_list[i+1];
    idx_ridge_vertex_facet_list[i+1] = idx_ridge_vertex_facet_list[i] + j;
    j = k;
  }

  int total_count = idx_ridge_vertex_facet_list[i] + 1;
  ridge_vertex_facet_list = new int[total_count];

  // Bookkeeping
  totalworkmemory += ((facetindex + 1) * sizeof(int) + 
                      totalvertices * sizeof(point *));

  idx2facetlist[0] = 0;
  for (i = 0, k = 0; i < facetindex; i++) {
    paryvertlist = (arraypool **) fastlookup(facetvertexlist, i);
    vertlist = *paryvertlist;
    idx2facetlist[i + 1] = (idx2facetlist[i] + (int) vertlist->objects);
    for (j = 0; j < vertlist->objects; j++) {
      parypt = (point *) fastlookup(vertlist, j);
      facetverticeslist[k] = *parypt;
      k++;
      // create ridge_vertex-to-facet map.
      int ridge_idx = pointmark(*parypt); // index of this ridge vertex
      // 'i' is the current facet index.
      ridge_vertex_facet_list[idx_ridge_vertex_facet_list[ridge_idx]] = i;
      // for the next facet index of this ridge vertex.
      idx_ridge_vertex_facet_list[ridge_idx]++;
    }
  }

  // Counters in idx_ridge_vertex_facet_list[] are shifted by 1.
  for (i = points->items; i >= 0; i--) {
    idx_ridge_vertex_facet_list[i+1] = idx_ridge_vertex_facet_list[i];
  }
  idx_ridge_vertex_facet_list[0] = 0;


  // Free the lists.
  for (i = 0; i < facetvertexlist->objects; i++) {
    paryvertlist = (arraypool **) fastlookup(facetvertexlist, i);
    vertlist = *paryvertlist;
    delete vertlist;
  }
  delete facetvertexlist;
}

//============================================================================//
//                                                                            //
// create_segment_facet_map()    Create the map from segments to adjacent     //
//                               facets.                                      //
//                                                                            //
//============================================================================//

void TetMeshCore::create_segment_facet_map()
{
  if (b->verbose > 0) {
    printf("  Creating the segment-to-facets map.\n");
  }
  if (idx_segment_facet_list != NULL) {
    delete [] idx_segment_facet_list;
    delete [] segment_facet_list;
  }

  face startsh, spinsh;
  face segloop;
  int segindex, facetidx;
  int totalcount = 0;
  int i;

  // both segment-index and facet-index start from zero.
  idx_segment_facet_list = new int[segmentendpointslist_length + 1];
  for (i = 0; i < segmentendpointslist_length + 1; i++) {
    idx_segment_facet_list[i] = 0;
  }

  subsegs->traversalinit();
  segloop.sh = shellfacetraverse(subsegs);
  while (segloop.sh != NULL) {
    segindex = getfacetindex(segloop);
    if (idx_segment_facet_list[segindex] == 0) {
      // Count the number of facets at this segment.
      spivot(segloop, startsh);
      spinsh = startsh;
      while (spinsh.sh != NULL) {
        idx_segment_facet_list[segindex]++;
        spivotself(spinsh);
        if (spinsh.sh == startsh.sh) break;
      }
      totalcount += idx_segment_facet_list[segindex];
    }
    segloop.sh = shellfacetraverse(subsegs);
  }

  // A working list.
  bool *bflags = new bool[segmentendpointslist_length + 1];

  // Have got the totalcount, fill the starting indices into the list.
  int j = idx_segment_facet_list[0], k;
  idx_segment_facet_list[0] = 0;
  //for (i = 0; i < segmentendpointslist_length + 1; i++) {
  for (i = 0; i < segmentendpointslist_length; i++) {
    k = idx_segment_facet_list[i+1];
    idx_segment_facet_list[i+1] = idx_segment_facet_list[i] + j;
    j = k;
    bflags[i] = false;
  }

  segment_facet_list = new int[totalcount + 1];
  subsegs->traversalinit();
  segloop.sh = shellfacetraverse(subsegs);
  while (segloop.sh != NULL) {
    segindex = getfacetindex(segloop);
    if (!bflags[segindex]) {
      spivot(segloop, startsh);
      spinsh = startsh;
      while (spinsh.sh != NULL) {
        facetidx = getfacetindex(spinsh);
        segment_facet_list[idx_segment_facet_list[segindex]] = facetidx;
        idx_segment_facet_list[segindex]++; // for the next one
        spivotself(spinsh);
        if (spinsh.sh == startsh.sh) break;
      }
      bflags[segindex] = true;
    }
    segloop.sh = shellfacetraverse(subsegs);
  }

  // Counters in idx_segment_facet_list[] are shifted by 1.
  for (i = segmentendpointslist_length - 1; i >= 0; i--) {
    idx_segment_facet_list[i+1] = idx_segment_facet_list[i];
  }
  idx_segment_facet_list[0] = 0;


  delete [] bflags;
}

//============================================================================//
//                                                                            //
// ridge_vertices_adjacent()    Check if two ridge vertices are connected by  //
//                              an input segment.                             //
//                                                                            //
//============================================================================//

int TetMeshCore::ridge_vertices_adjacent(point e1, point e2)
{
  int idx = pointmark(e1);
  int acount = idx_segment_ridge_vertex_list[idx+1]-idx_segment_ridge_vertex_list[idx];
  for (int i = 0; i < acount; i++) {
    if (segment_ridge_vertex_list[idx_segment_ridge_vertex_list[idx]+i] == e2) {
      return 1; // adjacent.
    }
  }
  return 0; // not adjacent.
}

//============================================================================//
//                                                                            //
// facet_ridge_vertex_adjacent()    Check if a facet and a ridge vertex is    //
//                                  adjacent by an input segment.             //
//                                                                            //
//============================================================================//

int TetMeshCore::facet_ridge_vertex_adjacent(face *chkfac, point chkpt)
{
  int ridge_idx = pointmark(chkpt);
  int facet_idx = getfacetindex(*chkfac);
  for (int i = idx_ridge_vertex_facet_list[ridge_idx];
           i < idx_ridge_vertex_facet_list[ridge_idx+1]; i++) {
    if (ridge_vertex_facet_list[i] == facet_idx) {
      return 1; // They are adjacent.
    }
  }
  return 0;
}

//============================================================================//
//                                                                            //
// segsegadjacent()    Check whether two segments, or a segment and a facet,  //
//                     or two facets are adjacent to each other.              //
//                                                                            //
//============================================================================//

int TetMeshCore::segsegadjacent(face *seg1, face *seg2)
{
  int segidx1 = getfacetindex(*seg1);
  int segidx2 = getfacetindex(*seg2);

  if (segidx1 == segidx2) {
    return 2; // Adjacent. They are the same segment.
  }

  point pa1 = segmentendpointslist[segidx1 * 2];
  point pb1 = segmentendpointslist[segidx1 * 2 + 1];
  point pa2 = segmentendpointslist[segidx2 * 2];
  point pb2 = segmentendpointslist[segidx2 * 2 + 1];

  if ((pa1 == pa2) || (pa1 == pb2) || (pb1 == pa2) || (pb1 == pb2)) {
    return 1; // Adjacent.
  }
  return 0; // not adjacent
}

//============================================================================//
//                                                                            //
// segfacetadjacent()    Check whether a segment and a facet are adjacent or  //
//                       not.                                                 //
//                                                                            //
//============================================================================//

int TetMeshCore::segfacetadjacent(face *subseg, face *subsh)
{
  int seg_idx   = getfacetindex(*subseg);
  int facet_idx = getfacetindex(*subsh);
  for (int i = idx_segment_facet_list[seg_idx];
           i < idx_segment_facet_list[seg_idx+1]; i++) {
    if (segment_facet_list[i] == facet_idx) {
      return 1; // They are adjacent.
    }
  }
  return 0;
}

//============================================================================//
//                                                                            //
// facetfacetadjacent()    Check whether two facets are adjacent or not.      //
//                                                                            //
//============================================================================//

int TetMeshCore::facetfacetadjacent(face *subsh1, face *subsh2)
{
  int count = 0, i;

  int fidx1 = getfacetindex(*subsh1);
  int fidx2 = getfacetindex(*subsh2);

  if (fidx1 == fidx2) {
    return 2; // Adjacent. They are the same facet.
  }

  for (i = idx2facetlist[fidx1]; i < idx2facetlist[fidx1+1]; i++) {
    pinfect(facetverticeslist[i]);
  }

  for (i = idx2facetlist[fidx2]; i < idx2facetlist[fidx2+1]; i++) {
    if (pinfected(facetverticeslist[i])) count++;
  }

  // Uninfect the vertices.
  for (i = idx2facetlist[fidx1]; i < idx2facetlist[fidx1+1]; i++) {
    puninfect(facetverticeslist[i]);
  }

  if (count > 0) {
    return 1;
  } else {
    return 0;
  }
}

//============================================================================//
//                                                                            //
// is_sharp_segment()    Check whether a given segment is sharp or not.       //
//                                                                            //
//============================================================================//

bool TetMeshCore::is_sharp_segment(face *seg)
{
  int segidx = getfacetindex(*seg);
  double mindihedang = segment_info_list[segidx*4];
  return mindihedang < 72.; // (in theory) < 72 degree is sufficient.
}

//============================================================================//
//                                                                            //
// does_seg_contain_acute_vertex()    Check whether one of the endpoints of a //
//                                    given segment is a sharp corner.        //
//                                                                            //
//============================================================================//

bool TetMeshCore::does_seg_contain_acute_vertex(face* seg)
{
  int segidx = getfacetindex(*seg);
  point *ppt = &(segmentendpointslist[segidx * 2]);
  double ang = 180.;
  // Get the smallest angle at its endpoints.
  for (int i = 0; i < 2; i++) {
    if ((ppt[i] == sorg(*seg)) || (ppt[i] == sdest(*seg))) {
      if (segment_info_list[segidx * 4 + 2 + i] < ang) {
        ang = segment_info_list[segidx * 4 + 2 + i];
      }
    }
  }
  return ang < 60.;
}

//============================================================================//
//                                                                            //
// create_a_shorter_edge()    Can we create an edge (which is shorter than    //
//                            minedgelength) between the two given vertices?  //
//                                                                            //
//============================================================================//

bool TetMeshCore::create_a_shorter_edge(point steinerpt, point nearpt)
{
  bool createflag = false; // default, do not create a shorter edge.

  enum verttype nearpt_type = pointtype(nearpt);
  enum verttype steiner_type = pointtype(steinerpt);

  if (nearpt_type == RIDGEVERTEX) {
    if (steiner_type == FREESEGVERTEX) {
      // Create a shorter edge if the Steiner point does not on an adjacent
      //   segment of this ridge vertex.
      face parentseg;
      sdecode(point2sh(steinerpt), parentseg);
      int segidx = getfacetindex(parentseg);
      point pa = segmentendpointslist[segidx * 2];
      point pb = segmentendpointslist[segidx * 2 + 1];
      if ((pa != nearpt) && (pb != nearpt)) {
        createflag = true; // create a shorter edge.
      }
    } else if (steiner_type == FREEFACETVERTEX) {
      // Create a shorter edge if the Steiner point does not on an adjacent
      //   facet of this ridge vertex.
      face parentsh;
      sdecode(point2sh(steinerpt), parentsh);
      if (!facet_ridge_vertex_adjacent(&parentsh, nearpt)) {
        createflag = true; // create a shorter edge.
      }
    }
  } else if (nearpt_type == FREESEGVERTEX) {
    if (steiner_type == FREESEGVERTEX) {
      // Check if they are on the same segment.
      face seg1, seg2;
      sdecode(point2sh(steinerpt), seg1);
      sdecode(point2sh(nearpt), seg2);
      int sidx1 = getfacetindex(seg1);
      int sidx2 = getfacetindex(seg2);
      if (sidx1 != sidx2) {
        createflag = true; // create a shorter edge.
      }
    } else if (steiner_type == FREEFACETVERTEX) {
      face parentseg, paresntsh;
      sdecode(point2sh(steinerpt), paresntsh);
      sdecode(point2sh(nearpt), parentseg);
      if (!segfacetadjacent(&parentseg, &paresntsh)) {
        createflag = true; // create a shorter edge.
      }
    }
  } else if (nearpt_type == FREEFACETVERTEX) {
    if (steiner_type == FREESEGVERTEX) {
      //assert(0); // to debug...
      face parentseg, paresntsh;
      sdecode(point2sh(nearpt), paresntsh);
      sdecode(point2sh(steinerpt), parentseg);
      if (!segfacetadjacent(&parentseg, &paresntsh)) {
        createflag = true; // create a shorter edge.
      }
    } else if (steiner_type == FREEFACETVERTEX) {
      // Create a short edge if they are on two different facets.
      face paresntsh1, paresntsh2;
      sdecode(point2sh(nearpt), paresntsh1);
      sdecode(point2sh(steinerpt), paresntsh2);
      int sidx1 = getfacetindex(paresntsh1);
      int sidx2 = getfacetindex(paresntsh2);
      if (sidx1 != sidx2) {
        createflag = true; // create a shorter edge.
      }
    }
  }

  return createflag;
}

//============================================================================//
//                                                                            //
// enqueuesubface()    Queue a subface or a subsegment for encroachment check.//
//                                                                            //
//============================================================================//

void TetMeshCore::enqueuesubface(memorypool *pool, face *chkface)
{
  if (!smarktest2ed(*chkface)) {
    smarktest2(*chkface); // Only queue it once.
    face *queface = (face *) pool->alloc();
    *queface = *chkface;
  }
}

//============================================================================//
//                                                                            //
// enqueuetetrahedron()    Queue a tetrahedron for quality check.             //
//                                                                            //
//============================================================================//

void TetMeshCore::enqueuetetrahedron(triface *chktet)
{
  if (!marktest2ed(*chktet)) {
    marktest2(*chktet); // Only queue it once.
    triface *quetet = (triface *) badtetrahedrons->alloc();
    *quetet = *chktet;
  }
}

//============================================================================//
//                                                                            //
// check_encroachment()    Check whether a given point encroaches upon a line //
//                         segment or not.                                    //
//                                                                            //
// 'checkpt' should not be dummypoint.                                        //
//                                                                            //
//============================================================================//

bool TetMeshCore::check_encroachment(point pa, point pb, point checkpt)
{
  // dot = (pa->checkpt) * (pb->checkpt)
  double d = (pa[0] - checkpt[0]) * (pb[0] - checkpt[0])
         + (pa[1] - checkpt[1]) * (pb[1] - checkpt[1])
         + (pa[2] - checkpt[2]) * (pb[2] - checkpt[2]);
  return d < 0.; // cos\theta < 0. ==> 90 < theta <= 180 degree.
}

//============================================================================//
//                                                                            //
// check_enc_segment()    Is a given segment encroached?                      //
//                                                                            //
//============================================================================//

bool TetMeshCore::check_enc_segment(face *chkseg, point *pencpt)
{
  point *ppt = (point *) &(chkseg->sh[3]);

  if (*pencpt != NULL) {
    return check_encroachment(ppt[0], ppt[1], *pencpt);
  }

  triface searchtet, spintet;
  point encpt = NULL, tapex;
  double prjpt[3]; // The projection point from encpt to segment.
  double minprjdist = 0., prjdist;
  int t1ver;

  sstpivot1(*chkseg, searchtet);
  spintet = searchtet;
  while (1) {
    tapex = apex(spintet);
    if (tapex != dummypoint) {
      if (check_encroachment(ppt[0], ppt[1], tapex)) {
        // Find one encroaching vertex. Calculate its projection distance
        projpt2edge(tapex, ppt[0], ppt[1], prjpt);
        prjdist = distance(tapex, prjpt);
        if (encpt == NULL) {
          encpt = tapex;
          minprjdist = prjdist;
        } else {
          if (prjdist < minprjdist) {
            encpt = tapex;
            minprjdist = prjdist;
          }
        }
      }
    }
    fnextself(spintet);
    if (spintet.tet == searchtet.tet) break;
  }

  if (encpt != NULL) {
    *pencpt = encpt; // Return this enc point.
    return true;
  }

  return false; // do not split it.
}

//============================================================================//
//                                                                            //
// get_steiner_on_segment()    Get the Steiner point to split a given segment.//
//                                                                            //
//============================================================================//

bool TetMeshCore::get_steiner_on_segment(face* seg, point refpt, point steinpt)
{
  point ei = sorg(*seg);
  point ej = sdest(*seg);
  //if (*prefpt == NULL) {
  //  // Check if this segment is encroached by some existing vertices.
  //  assert(0); // to do ...
  //}
  // Is this segment contains an acute seg-seg angle?
  bool acute_flag = false;
  int i;

  if ((refpt) != NULL) {
    // This segment is encroched by an existing vertex.
    double L, L1, t;
  
    if (pointtype(refpt) == FREESEGVERTEX) {
      face parentseg;
      sdecode(point2sh(refpt), parentseg);
      int sidx1 = getfacetindex(parentseg);
      point far_pi = segmentendpointslist[sidx1 * 2];
      point far_pj = segmentendpointslist[sidx1 * 2 + 1];
      int sidx2 = getfacetindex(*seg);
      point far_ei = segmentendpointslist[sidx2 * 2];
      point far_ej = segmentendpointslist[sidx2 * 2 + 1];
      if ((far_pi == far_ei) || (far_pj == far_ei)) {
        // Two segments are adjacent at far_ei!
        // Create a Steiner point at the intersection of the segment
        //   [far_ei, far_ej] and the sphere centered at far_ei with 
        //   radius |far_ei - refpt|.
        L = distance(far_ei, far_ej);
        L1 = distance(far_ei, refpt);
        t = L1 / L;
        for (i = 0; i < 3; i++) {
          steinpt[i] = far_ei[i] + t * (far_ej[i] - far_ei[i]);
        }
        double lfs_at_steiner = distance(refpt, steinpt);
        //double dist_to_ei = distance(steinpt, ei);
        double dist_to_ej = distance(steinpt, ej);
        if (/*(dist_to_ei < lfs_at_steiner) ||*/
            (dist_to_ej < lfs_at_steiner)) {
          // Split the point at the middle.
          for (i = 0; i < 3; i++) {
            steinpt[i] = ei[i] + 0.5 * (ej[i] - ei[i]);
          }
        }
        set_ridge_vertex_protecting_ball(far_ei);
        acute_flag = true;
      } else if ((far_pi == far_ej) || (far_pj == far_ej)) {
        // Two segments are adjacent at far_ej!
        L = distance(far_ei, far_ej);
        L1 = distance(far_ej, refpt);
        t = L1 / L;
        for (i = 0; i < 3; i++) {
          steinpt[i] = far_ej[i] + t * (far_ei[i] - far_ej[i]);
        }
        double lfs_at_steiner = distance(refpt, steinpt);
        double dist_to_ei = distance(steinpt, ei);
        //double dist_to_ej = distance(steinpt, ej);
        if ((dist_to_ei < lfs_at_steiner) /*||
            (dist_to_ej < lfs_at_steiner)*/) {
          // Split the point at the middle.
          for (i = 0; i < 3; i++) {
            steinpt[i] = ei[i] + 0.5 * (ej[i] - ei[i]);
          }
        }
        set_ridge_vertex_protecting_ball(far_ej);
        acute_flag = true;
      } else {
        // Cut the segment by the projection point of refpt.
        projpt2edge(refpt, ei, ej, steinpt);
        double lfs_at_steiner = distance(refpt, steinpt);
        double dist_to_ei = distance(steinpt, ei);
        double dist_to_ej = distance(steinpt, ej);
        if ((dist_to_ei < lfs_at_steiner) ||
            (dist_to_ej < lfs_at_steiner)) {
          // Split the point at the middle.
          for (i = 0; i < 3; i++) {
            steinpt[i] = ei[i] + 0.5 * (ej[i] - ei[i]);
          }
        }
      }
    } else if (pointtype(refpt) == RIDGEVERTEX) {
      int sidx2 = getfacetindex(*seg);
      point far_ei = segmentendpointslist[sidx2 * 2];
      point far_ej = segmentendpointslist[sidx2 * 2 + 1];
      if (ridge_vertices_adjacent(far_ei, refpt)) {
        // Thjey are adjacent at far_ei.
        // Create a Steiner point at the intersection of the segment
        //   [far_ei, far_ej] and the sphere centered at far_ei with 
        //   radius |far_ei - refpt|.
        L = distance(far_ei, far_ej);
        L1 = distance(far_ei, refpt);
        t = L1 / L;
        for (i = 0; i < 3; i++) {
          steinpt[i] = far_ei[i] + t * (far_ej[i] - far_ei[i]);
        }
        double lfs_at_steiner = distance(refpt, steinpt);
        //double dist_to_ei = distance(steinpt, ei);
        double dist_to_ej = distance(steinpt, ej);
        if (/*(dist_to_ei < lfs_at_steiner) ||*/
            (dist_to_ej < lfs_at_steiner)) {
          // Split the point at the middle.
          for (i = 0; i < 3; i++) {
            steinpt[i] = ei[i] + 0.5 * (ej[i] - ei[i]);
          }
        }
        set_ridge_vertex_protecting_ball(far_ei);
        acute_flag = true;
      } else if (ridge_vertices_adjacent(far_ej, refpt)) {
        // Calulate a new point.
        L = distance(far_ei, far_ej);
        L1 = distance(far_ej, refpt);
        t = L1 / L;
        for (i = 0; i < 3; i++) {
          steinpt[i] = far_ej[i] + t * (far_ei[i] - far_ej[i]);
        }
        double lfs_at_steiner = distance(refpt, steinpt);
        double dist_to_ei = distance(steinpt, ei);
        //double dist_to_ej = distance(steinpt, ej);
        if ((dist_to_ei < lfs_at_steiner) /*||
            (dist_to_ej < lfs_at_steiner)*/) {
          // Split the point at the middle.
          for (i = 0; i < 3; i++) {
            steinpt[i] = ei[i] + 0.5 * (ej[i] - ei[i]);
          }
        }
        set_ridge_vertex_protecting_ball(far_ej);
        acute_flag = true;
      } else {
        // Cut the segment by the projection point of refpt.
        projpt2edge(refpt, ei, ej, steinpt);
        double lfs_at_steiner = distance(refpt, steinpt);
        double dist_to_ei = distance(steinpt, ei);
        double dist_to_ej = distance(steinpt, ej);
        if ((dist_to_ei < lfs_at_steiner) ||
            (dist_to_ej < lfs_at_steiner)) {
          // Split the point at the middle.
          for (i = 0; i < 3; i++) {
            steinpt[i] = ei[i] + 0.5 * (ej[i] - ei[i]);
          }
        }
       }
    } else if (pointtype(refpt) == FREEFACETVERTEX) {
      // Cut the segment by the projection point of refpt.
      projpt2edge(refpt, ei, ej, steinpt);
      double lfs_at_steiner = distance(refpt, steinpt);
      double dist_to_ei = distance(steinpt, ei);
      double dist_to_ej = distance(steinpt, ej);
      if ((dist_to_ei < lfs_at_steiner) ||
          (dist_to_ej < lfs_at_steiner)) {
        // Split the point at the middle.
        for (i = 0; i < 3; i++) {
          steinpt[i] = ei[i] + 0.5 * (ej[i] - ei[i]);
        }
      }
    } else {
      // Cut the segment by the projection point of refpt.
      projpt2edge(refpt, ei, ej, steinpt);
      // Make sure that steinpt is not too close to ei and ej.
      double lfs_at_steiner = distance(refpt, steinpt);
      double dist_to_ei = distance(steinpt, ei);
      double dist_to_ej = distance(steinpt, ej);
      if ((dist_to_ei < lfs_at_steiner) ||
          (dist_to_ej < lfs_at_steiner)) {
        // Split the point at the middle.
        for (i = 0; i < 3; i++) {
          steinpt[i] = ei[i] + 0.5 * (ej[i] - ei[i]);
        }
      }
    }

    // Make sure that steinpt is not too close to ei and ej.
  } else {
    // Split the point at the middle.
    for (i = 0; i < 3; i++) {
      steinpt[i] = ei[i] + 0.5 * (ej[i] - ei[i]);
    }
  }


  return acute_flag;
}

//============================================================================//
//                                                                            //
// split_segment()    Split a given segment.                                  //
//                                                                            //
// If param != NULL, it contains the circumcenter and its insertion radius    //
// of a bad quality tetrahedron. Tghis circumcenter is rejected since it      //
// encroaches upon this segment.                                              //
//                                                                            //
//============================================================================//

bool TetMeshCore::split_segment(face *splitseg, point encpt, double *param,
                               int qflag, int chkencflag, int *iloc)
{
  triface searchtet;
  face searchsh;
  point newpt;
  insertvertexflags ivf;


  insert_point_count++;
  if (!b->quiet && (b->refine_progress_ratio > 0)) {
    if (insert_point_count >= report_refine_progress) {
      printf("  %ld insertions, added %ld points",
             insert_point_count - last_insertion_count,
             points->items - last_point_count);
      last_point_count = points->items; // update it.
      last_insertion_count = insert_point_count;
      if (check_tets_list->objects > 0l) {
        printf(", %ld tetrahedra in queue.\n", check_tets_list->objects);
      } else if (split_subfaces_pool->items > 0l) {
        printf(", %ld subfaces in queue.\n", split_subfaces_pool->items);
      } else {
        printf(", %ld segments in queue.\n", split_segments_pool->items);
      }
      // The next report event
      report_refine_progress *= (1. + b->refine_progress_ratio);
    }
  }
  // Is this segment shared by two facets form an acute dihedral angle?
  int segidx = getfacetindex(*splitseg);
  bool is_sharp = is_sharp_segment(splitseg);

  if (!qflag && (encpt == NULL)) {
    // The split of this segment is due to a rejected ccent of a bad quality
    //   subface or a tetrahedron.
    if (is_sharp) {
      // Do not split a sharp segment.
      *iloc = (int) SHARPCORNER;
      return false;
    }
    // Do not split this segment if one of its endpoints is a sharp corner.
    if (does_seg_contain_acute_vertex(splitseg)) {
      *iloc = (int) SHARPCORNER;
      return false;
    }
  }
  
  // We need to know whether the segment of the new point is adjacent
  //   to another segment which contains the encroached point (encpt).
  makepoint(&newpt, FREESEGVERTEX);
  get_steiner_on_segment(splitseg, encpt, newpt);

  // For create_a_shorter_edge() called in insertpoint().
  setpoint2sh(newpt, sencode(*splitseg));
  
  // Split the segment by the Bowyer-Watson algorithm.
  sstpivot1(*splitseg, searchtet);
  ivf.iloc = (int) ONEDGE;
  ivf.bowywat = 3; // Use Bowyer-Watson, preserve subsegments and subfaces;
  ivf.validflag = 1; // Validate the B-W cavity.
  ivf.lawson = 2; // Do flips to recover Delaunayness.
  ivf.rejflag = 0;     // Do not check encroachment of new segments/facets.
  if (b->metric) {
    ivf.rejflag |= 4;  // Do check encroachment of protecting balls.
  }
  ivf.chkencflag = chkencflag;
  ivf.sloc = (int) INSTAR; // ivf.iloc;
  ivf.sbowywat = 3; // ivf.bowywat;  // Surface mesh options.
  ivf.splitbdflag = 1;
  ivf.respectbdflag = 1;
  ivf.assignmeshsize = b->metric;

  ivf.smlenflag = useinsertradius; // Return the distance to its nearest vertex.
  // Reject a near Steiner point on this segment when:
  //  - it is only encroached by a rejected circumcenter, or
  //  - the insertion of the reject ccent is not due to mesh size (qflag).
  if (!qflag) { //if (!is_adjacent || !qflag) {
    ivf.check_insert_radius = useinsertradius;
  }
  ivf.parentpt = NULL;
  
  if (insertpoint(newpt, &searchtet, &searchsh, splitseg, &ivf)) {
    st_segref_count++;
    if (steinerleft > 0) steinerleft--;
    if (useinsertradius) {
      double rv = 0.0; // param[3]; // emin, maybe zero.
      
      if (is_sharp) {
        // A Steiner point on a sharp segment needs insertion radius.
        // Default use the distance to its neartest vertex.
        double L = ivf.smlen * 0.95; // (ivf.smlen / 3.);
        // Choose the larger one between param[3] and L
        rv = (param[3] > L ? param[3] : L);
        // Record the minimum insertion radius for this segment.
        double minradius = segment_info_list[segidx*4+1];
        if (minradius == 0.) {
          minradius = rv;
        } else {
          if (rv < minradius) minradius = rv;
        }
        segment_info_list[segidx*4+1] = minradius;
      }

      setpointinsradius(newpt, rv); // ivf.smlen
      setpoint2ppt(newpt, ivf.parentpt);
      if (ivf.smlen < smallest_insradius) { // rv?
        smallest_insradius = ivf.smlen;
      }
    }
    if (flipstack != NULL) {
      flipconstraints fc;
      fc.chkencflag = chkencflag;
      fc.enqflag = 2;
      lawsonflip3d(&fc);
      //unflipqueue->restart();
    }
    
    if (later_unflip_queue->objects > b->unflip_queue_limit) {
      recoverdelaunay();
    }
    
    *iloc = ivf.iloc;
    return true;
  } else {
    // Point is not inserted.
    if (ivf.iloc == (int) NEARVERTEX) {
      terminate_tet_core(this, 2); // report a bug.
    }


    pointdealloc(newpt);
    
    *iloc = ivf.iloc;
    return false;
  }
}

//============================================================================//
//                                                                            //
// repairencsegs()    Repair encroached (sub) segments.                       //
//                                                                            //
//============================================================================//

void TetMeshCore::repairencsegs(double *param, int qflag, int chkencflag)
{
  int split_count = 0, rej_count = 0;
  bool ref_segment = ((b->cdtrefine & 1) > 0); // -D1, -D3, -D5, -D7

  while (ref_segment &&
         ((badsubsegs->items > 0) || (split_segments_pool->items > 0))) {
  
    if (badsubsegs->items > 0) {
      badsubsegs->traversalinit();
      face *bface = (face *) badsubsegs->traverse();
      while (bface != NULL) {
        // A queued segment may have been deleted (split).
        if ((bface->sh != NULL) && (bface->sh[3] != NULL)) {
          // A queued segment may have been processed.
          if (smarktest2ed(*bface)) {
            sunmarktest2(*bface);
            point encpt = NULL;
            if (check_enc_segment(bface, &encpt)) {
              badface *bf = (badface *) split_segments_pool->alloc();
              bf->init();
              bf->ss = *bface;
              bf->forg = sorg(*bface);
              bf->fdest = sdest(*bface);
              bf->noppo = encpt;
              // Push it onto stack.
              bf->nextitem = stack_enc_segments;
              stack_enc_segments = bf;
            }
          }
        }
        bface = (face *) badsubsegs->traverse();
      } // while (bface != NULL)
      badsubsegs->restart();
    } // if (badsubsegs->items > 0)
  
    if (split_segments_pool->items == 0) break;

    // Stop if we have used the desried number of Steiner points.
    if (steinerleft == 0) break;
    // Stop if the desried number of tetrahedra is reached.
    if ((elem_limit > 0) &&
        ((tetrahedrons->items - hullsize) > elem_limit)) break;

    // Pop up an encroached segment.
    badface *bf = stack_enc_segments;
    stack_enc_segments = bf->nextitem;
    if ((bf->ss.sh != NULL) &&
        (sorg(bf->ss) == bf->forg) &&
        (sdest(bf->ss) == bf->fdest)) {
     int iloc = (int) UNKNOWN;
      split_count++;
      if (!split_segment(&(bf->ss), bf->noppo, param, qflag, chkencflag, &iloc)) {
        rej_count++;
      }
    }
    // Return this badface to the pool.
    split_segments_pool->dealloc((void *) bf);
  }

  if (b->verbose > 2) {
    printf("    Trying to split %d segments, %d were rejected.\n",
           split_count, rej_count);
  }

  if (badsubsegs->items > 0) {
    // Clean this list (due to ref_segment).
    badsubsegs->traversalinit();
    face *bface = (face *) badsubsegs->traverse();
    while (bface != NULL) {
      // A queued segment may have been deleted (split).
      if ((bface->sh != NULL) && (bface->sh[3] != NULL)) {
        // A queued segment may have been processed.
        if (smarktest2ed(*bface)) {
          sunmarktest2(*bface);
        }
      }
      bface = (face *) badsubsegs->traverse();
    } // while (bface != NULL)
    badsubsegs->restart();
  } // if (badsubsegs->items > 0)

  if (split_segments_pool->items > 0) {
    if (steinerleft == 0) {
      if (b->verbose) {
        printf("The desired number of Steiner points is reached.\n");
      }
    } else if (elem_limit > 0) {
      if (b->verbose) {
        printf("The desired number %ld of elements is reached.\n", elem_limit);
      }
    }
    split_segments_pool->restart();
    stack_enc_segments = NULL;
  }
}

//============================================================================//
//                                                                            //
// get_subface_ccent()    Calculate the circumcenter of the diametrical circ- //
//                        umsphere of a given subface.                        //
//                                                                            //
//============================================================================//

bool TetMeshCore::get_subface_ccent(face *chkfac, double *pos)
{
  point P = (point) chkfac->sh[3];
  point Q = (point) chkfac->sh[4];
  point R = (point) chkfac->sh[5];

  if (circumsphere(P, Q, R, NULL, pos, NULL)) {
    return true;
  } else {
    terminate_tet_core(this, 2);
    return false;
  }

}

//============================================================================//
//                                                                            //
// check_enc_subface()    Check if a given subface is encroached or not.      //
//                                                                            //
//============================================================================//

bool TetMeshCore::check_enc_subface(face *chkfac, point *pencpt, double *ccent,
                                   double *radius)
{
  triface adjtet;
  point encpt = NULL, pa, pb, pc, toppo;
  double prjpt[3], minprjdist = 0., prjdist;
  double ori;
  int t1ver;
  
  //get_subface_ccent(chkfac, ccent);
  double rd = distance(ccent, sorg(*chkfac));
  *radius = rd;

  if (*pencpt != NULL) {
    // This is only used during the insertion of a Steiner point.
    double len = distance(ccent, *pencpt);
    if ((fabs(len - rd) / rd) < 1e-3) len = rd; // Rounding.
    if (len < rd) {
      return true;
    }
    return false;
  }

  stpivot(*chkfac, adjtet);
  if (adjtet.tet == NULL) {
    // This subface is not attached to any tet.
    return false;
  }
  for (int i = 0; i < 2; i++) {
    toppo = oppo(adjtet);
    if (toppo != dummypoint) {
      double len = distance(ccent, toppo);
      //if ((fabs(len - rd) / rd) < b->epsilon) len = rd; // Rounding.
      if ((fabs(len - rd) / rd) < 1e-3) len = rd; // Rounding.
      if (len < rd) {
        int adjacent = 0; // not adjacent
        if (pointtype(toppo) == RIDGEVERTEX) {
          adjacent = facet_ridge_vertex_adjacent(chkfac, toppo);
        } else if (pointtype(toppo) == FREESEGVERTEX) {
          face parentseg;
          sdecode(point2sh(toppo), parentseg);
          adjacent = segfacetadjacent(&parentseg, chkfac);
        } else if (pointtype(toppo) == FREEFACETVERTEX) {
          face parentsh;
          sdecode(point2sh(toppo), parentsh);
          int facidx1 = getfacetindex(parentsh);
          int facidx2 = getfacetindex(*chkfac);
          if (facidx1 == facidx2) {
            adjacent = 1; // They are on the same facet.
          }
        }
        if (adjacent) {
          // They are adjacent and they are on the same facet.
          flippush(flipstack, &adjtet);
          return false;
        }
        pa = org(adjtet);
        pb = dest(adjtet);
        pc = apex(adjtet);
        projpt2face(toppo, pa, pb, pc, prjpt);
        ori = orient3d(pa, pb, toppo, prjpt);
        if (ori >= 0) {
          ori = orient3d(pb, pc, toppo, prjpt);
          if (ori >= 0) {
            ori = orient3d(pc, pa, toppo, prjpt);
            if (ori >= 0) {
              prjdist = distance(toppo, prjpt);
              if (encpt == NULL) {
                encpt = toppo;
                minprjdist = prjdist;
              } else {
                if (prjdist < minprjdist) {
                  encpt = toppo;
                  minprjdist = prjdist;
                }
              }
            } // if (ori >= 0)
          } // if (ori >= 0)
        } // if (ori >= 0)
      } // if (len < rd)
    }
    fsymself(adjtet);
  }

  if (encpt != NULL) {
    *pencpt = encpt;
    return true;
  }

  return false; // this subface is not encroached.
}

//============================================================================//
//                                                                            //
// check_subface()    Is a given subface in a bad shape (radius-edge ratio)?  //
//                                                                            //
//============================================================================//

bool TetMeshCore::check_subface(face *chkfac, double *ccent, double radius, double *param)
{

  // Get the shortest edge length.
  double emin = 1.e+30, dist;
  int shver = 0;
  for (chkfac->shver = 0; chkfac->shver < 3; chkfac->shver++) {
    dist = distance(sorg(*chkfac), sdest(*chkfac));
    if (dist < emin) {
      emin = dist;
      shver = chkfac->shver;
    }
  }
  chkfac->shver = shver;

  double ratio = radius / emin;
  if (ratio > b->minratio) {
    // Set a small value to protect this vertex (refer to J. Shewchuk).
    // Enlarge the insertion radius (due to small angle)
    point pa = sorg(*chkfac);
    point pb = sdest(*chkfac);
    double ra = getpointinsradius(pa);
    double rb = getpointinsradius(pb);
    if (ra > 0.) {
      if (ra > emin) {
        emin = ra;
      }
    }
    if (rb > 0.) {
      if (rb > emin) {
        emin = rb;
      }
    }

    param[3] = emin; // emin / 3.; // (emin * b->minratio);
    param[4] = ratio;
    param[5] = 0.; // not used.
    return true; // need to split it.
  }

  return false;
}

//============================================================================//
//                                                                            //
// enqueue_subface()    Push a badly-shaped subface into the priority queue.  //
//                                                                            //
//============================================================================//

void TetMeshCore::enqueue_subface(face *bface, point encpt, double *ccent, double *param)
{
  badface *bf = (badface *) split_subfaces_pool->alloc();
  bf->init();
  bf->ss = *bface;
  bf->forg  =  sorg(*bface);
  bf->fdest = sdest(*bface);
  bf->fapex = sapex(*bface);
  bf->noppo = encpt;
  int i;
  for (i = 0; i < 3; i++) bf->cent[i] = ccent[i];
  for (i = 3; i < 6; i++) bf->cent[i] = param[i];

  if (encpt != NULL) {
    // Push it into the encroaching stack.
    bf->nextitem = stack_enc_subfaces;
    stack_enc_subfaces = bf;
  } else {
    // Push it into the priority queue.
    double qual = 1.0;
    if (param[4] > 1.) {
      qual = 1.0 / param[4]; // 1 / radius_edge_ratio.
    }
    // Determine the appropriate queue to put the bad subface into.
    int queuenumber = 0;
    if (qual < 1) {
      queuenumber = (int) (64.0 * (1 - qual));
      if (queuenumber > 63) {
        queuenumber = 63;
      }
    } else {
      // It's not a bad shape; put the subface in the lowest-priority queue.
      queuenumber = 0;
    }
    
    // Are we inserting into an empty queue?
    if (queuefront[queuenumber] == (badface *) NULL) {
      // Yes, we are inserting into an empty queue.
      //   Will this become the highest-priority queue?
      if (queuenumber > firstnonemptyq) {
        // Yes, this is the highest-priority queue.
        nextnonemptyq[queuenumber] = firstnonemptyq;
        firstnonemptyq = queuenumber;
      } else {
        // No, this is not the highest-priority queue.
        //   Find the queue with next higher priority.
        int i = queuenumber + 1;
        while (queuefront[i] == (badface *) NULL) {
          i++;
        }
        // Mark the newly nonempty queue as following a higher-priority queue.
        nextnonemptyq[queuenumber] = nextnonemptyq[i];
        nextnonemptyq[i] = queuenumber;
      }
      // Put the bad subface at the beginning of the (empty) queue.
      queuefront[queuenumber] = bf;
    } else {
      // Add the bad tetrahedron to the end of an already nonempty queue.
      queuetail[queuenumber]->nextitem = bf;
    }
    // Maintain a pointer to the last subface of the queue.
    queuetail[queuenumber] = bf;
  }
}

// Return the subface at the front of the queue.
TetMeshCore::badface* TetMeshCore::top_subface()
{
  if (stack_enc_subfaces != NULL) {
    return stack_enc_subfaces;
  } else {
    // Keep a record of which queue was accessed in case dequeuebadtetra()
    //   is called later.
    recentq = firstnonemptyq;
    // If no queues are nonempty, return NULL.
    if (firstnonemptyq < 0) {
      return (badface *) NULL;
    } else {
      // Return the first tetrahedron of the highest-priority queue.
      return queuefront[firstnonemptyq];
    }
  }
}

//============================================================================//
//                                                                            //
// dequeue_subface()    Popup a badly-shaped subface from the priority queue. //
//                                                                            //
//============================================================================//

void TetMeshCore::dequeue_subface()
{
  badface *bf;
  int i;

  if (stack_enc_subfaces != NULL) {
    bf = stack_enc_subfaces;
    stack_enc_subfaces = bf->nextitem;
    // Return the bad subface to the pool.
    split_subfaces_pool->dealloc((void *) bf);
  } else {
    // If queues were empty last time topbadtetra() was called, do nothing.
    if (recentq >= 0) {
      // Find the tetrahedron last returned by topbadtetra().
      bf = queuefront[recentq];
      // Remove the tetrahedron from the queue.
      queuefront[recentq] = bf->nextitem;
      // If this queue is now empty, update the list of nonempty queues.
      if (bf == queuetail[recentq]) {
        // Was this the highest-priority queue?
        if (firstnonemptyq == recentq) {
          // Yes; find the queue with next lower priority.
          firstnonemptyq = nextnonemptyq[firstnonemptyq];
        } else {
          // No; find the queue with next higher priority.
          i = recentq + 1;
          while (queuefront[i] == (badface *) NULL) {
            i++;
          }
          nextnonemptyq[i] = nextnonemptyq[recentq];
        }
      }
      // Return the bad subface to the pool.
      split_subfaces_pool->dealloc((void *) bf);
    }
  }
}

//============================================================================//
//                                                                            //
// parallel_shift()    Parallel shift a triangle along its normal.            //
//                                                                            //
// Given a triangle (a, b, c), create a parallel triangle (pa, pb, pc) at a   //
// distance above (a, b, c).                                                  //
//                                                                            //
//============================================================================//

void TetMeshCore::parallel_shift(point pa, point pb, point pc,
                                point pt, double* ppt)
{
  // Get the normal and the average edge length of this triangle.
  double N[3], Lav;
  facenormal(pa, pb, pc, N, 1, &Lav);

  // Normalize the normal.
  double L = sqrt(N[0]*N[0]+N[1]*N[1]+N[2]*N[2]);
  N[0] /= L;
  N[1] /= L;
  N[2] /= L;
  
  // Calculate the shifted vertices.
  for (int i = 0; i < 3; i++) {
    ppt[0] = pt[0] + Lav * N[0];
    ppt[1] = pt[1] + Lav * N[1];
    ppt[2] = pt[2] + Lav * N[2];
  }

}

//============================================================================//
//                                                                            //
// locate_on_surface()    Locate a vertex in a facet.                         //
//                                                                            //
//============================================================================//

enum TetMeshCore::locateresult
TetMeshCore::locate_on_surface(point searchpt, face* searchsh)
{
  enum locateresult loc = OUTSIDE;

  triface searchtet;
  stpivot(*searchsh, searchtet);
  if (ishulltet(searchtet)) {
    sesymself(*searchsh);
    stpivot(*searchsh, searchtet);
  }

  // Select an edge such that pt lies to CCW of it.
  point pa, pb, pc;
  double toppo[3]; // a parallel-shifted point
  double n1[3], n2[3], cosang;
  int t1ver; // used by fnextself()
  int i;
  
  for (i = 0; i < 3; i++) {
    pa =  org(searchtet);
    pb = dest(searchtet);
    pc = apex(searchtet);
    parallel_shift(pa, pb, pc, pa, toppo);
    if (orient3d(pa, pb, toppo, searchpt) > 0) {
      break;
    }
    enextself(searchtet);
  }
  if (i == 3) {
    terminate_tet_core(this, 2);
  }
  
  while (true) {

    // Let E = [a,b,c] and p lies to the CCW of [a->b].
    // Make sure that the searching vertex and the current subface (a,b,c) are
    //   (nearly) coplanar. We check the dihedral angle between (a,b,c) and
    //   (a,b,searchpt). If it is within the tolerance of co-planar facets,
    //   then we continue the search, otherwise, the search is stopped.
    facenormal(pa, pb, pc, n1, 1, NULL);
    facenormal(pb, pa, searchpt, n2, 1, NULL);
    cosang = dot(n1, n2) / (sqrt(dot(n1, n1)) * sqrt(dot(n2, n2)));
    if (cosang > cos_facet_separate_ang_tol) {
      // The searching vertex is not coplanar with this subface.
      loc = NONCOPLANAR;
      break;
    }

    parallel_shift(pa, pb, pc, pc, toppo);
    double ori1 = orient3d(pb, pc, toppo, searchpt);
    double ori2 = orient3d(pc, pa, toppo, searchpt);

    if (ori1 > 0) {
      if (ori2 > 0) {
        //break; // Found.
        loc = ONFACE; break;
      } else if (ori2 < 0) {
        //E.ver = _eprev_tbl[E.ver];
        eprevself(searchtet);
      } else { // ori2 == 0
        //E.ver = _eprev_tbl[E.ver];
        //return LOC_ON_EDGE; // ONEDGE  p lies on edge [c,a]
        eprevself(searchtet);
        loc = ONEDGE; break;
      }
    } else if (ori1 < 0) {
      if (ori2 > 0) {
        //E.ver = _enext_tbl[E.ver];
        enextself(searchtet);
      } else if (ori2 < 0) {
        // Randomly choose one.
        if (rand() % 2) { // flipping a coin.
          //E.ver = _enext_tbl[E.ver];
          enextself(searchtet);
        } else {
          //E.ver = _eprev_tbl[E.ver];
          eprevself(searchtet);
        }
      } else { // ori2 == 0
        //E.ver = _enext_tbl[E.ver];
        enextself(searchtet);
      }
    } else { // ori1 == 0
      if (ori2 > 0) {
        //E.ver = _enext_tbl[E.ver]; // p lies on edge [b,c].
        //return LOC_ON_EDGE; // ONEDGE
        enextself(searchtet);
        loc = ONEDGE; break;
      } else if (ori2 < 0) {
        //E.ver = _eprev_tbl[E.ver];
        eprevself(searchtet);
      } else { // ori2 == 0
        //E.ver = _eprev_tbl[E.ver]; // p is coincident with apex.
        //return LOC_ON_VERT; // ONVERTEX Org(E)
        eprevself(searchtet);
        loc = ONVERTEX; break;
      }
    }

    // Check if we want to cross a segment.
    if (issubseg(searchtet)) {
      loc = ENCSEGMENT; break;
    }

    // Goto the adjacent subface at this subedge.
    int fcount = 0;
    while (fcount < 100000) {
      esymself(searchtet);
      if (issubface(searchtet)) break;
      fsymself(searchtet);
      fcount++;
    }
    if (!issubface(searchtet)) {
      terminate_tet_core(this, 2); // report a bug
    }

    // Update the vertices.
    pa =  org(searchtet);
    pb = dest(searchtet);
    pc = apex(searchtet);
    //toppo = oppo(searchtet);
  } // while (true)

  tspivot(searchtet, *searchsh);

  return loc;
}

//============================================================================//
//                                                                            //
// split_subface()    Split a subface.                                        //
//                                                                            //
// param[6], it contains the following data:                                  //
//   [0],[1],[2] - the location of a rejected circumcent,                     //
//   [3] - the samllest edge length ( = insertion radius)                     //
//   [4] - ratio-edge ratio (of this subface).                                //
//         If it is zero, it is an encroached subface.                        //
//   [5] - no used.                                                           //
///                                                                           //
//============================================================================//

bool TetMeshCore::split_subface(face *splitfac, point encpt, double *ccent,
                               double *param, int qflag, int chkencflag, int *iloc)
{
  triface searchtet;
  face searchsh;
  insertvertexflags ivf;
  point newpt, bak_pts[3], *ppt;
  bool is_adjacent = false;
  bool splitflag = false;  // Indicate if any Steiner point is added.
  int i;

  insert_point_count++;
  if (!b->quiet && (b->refine_progress_ratio > 0.)) {
    if (insert_point_count >= report_refine_progress) {
      printf("  %ld insertions, added %ld points",
             insert_point_count - last_insertion_count,
             points->items - last_point_count);
      last_point_count = points->items; // update it.
      last_insertion_count = insert_point_count;
      if (check_tets_list->objects > 0l) {
        printf(", %ld tetrahedra in queue.\n", check_tets_list->objects);
      } else {
        printf(", %ld subfaces in queue.\n", split_subfaces_pool->items);
      }
      // The next report event
      report_refine_progress *= (1. + b->refine_progress_ratio);
    }
  }

  // Check if this subface is adjacent to a sharp segment, i.e., it is incident
  //   by two facets which form an acute dihedral angle.
  face checkface = *splitfac;
  face checkseg;
  for (i = 0; i < 3; i++) {
    sspivot(checkface, checkseg);
    if (checkseg.sh != NULL) {
      if (is_sharp_segment(&checkseg)) {
        is_adjacent = true;
        break;
      }
    }
    senext2self(checkface);
  }

  if (is_adjacent) {
    // Only split it either it is a bad quality triangle, or due to the
    //   qflag, i.e., mesh size requirement.
    if (!qflag) {
      if (encpt != NULL) {
        *iloc = (int) SHARPCORNER;
        return false; // reject splitting this subface.
      } else {
        if (param[4] == 0.0) {
          // It is not a bad quality subface.
          *iloc = (int) SHARPCORNER;
          return false; // reject splitting this subface.
        }
      }
    }
  } // if (is_adjacent)


  // Deciding the inserting point.
  if (encpt != NULL) {
    // Insert at the projection of the encpt on the facet.
    double pos[3];
    ppt = (point *) &(splitfac->sh[3]);
    projpt2face(encpt, ppt[0], ppt[1], ppt[2], pos);
    makepoint(&newpt, FREEFACETVERTEX);
    for (i = 0; i < 3; i++) newpt[i] = pos[i];

    //if (is_adjacent) {
    // Check whether this new position is too close to an existing vertex.
    double prjdist = distance(encpt, newpt);
    double dist, mindist = 1.e+30;
    for (i = 0; i < 3; i++) {
      dist = distance(ppt[i], newpt);
      if (dist < mindist) mindist = dist;
    }
    if (mindist < prjdist) {
      // Use the circumcenter of this triange instead of the proj of encpt.
      for (i = 0; i < 3; i++) newpt[i] = ccent[i];
    }
    //}
  } else {
    // Split the subface at its circumcenter.
    makepoint(&newpt, FREEFACETVERTEX);
    for (i = 0; i < 3; i++) newpt[i] = ccent[i];
  }

  // This info is needed by create_a_shorter_edge() (called in insertpoint()).
  setpoint2sh(newpt, sencode(*splitfac));


  searchsh = *splitfac;
  ivf.iloc = (int) locate_on_surface(newpt, &searchsh);

  if (ivf.iloc == (int) ENCSEGMENT) {
    // Point lies in the outside of the facet.
    pointdealloc(newpt);
    *iloc = FENSEDIN; // it is a fested in vertex.
    return splitflag;
  } else if (ivf.iloc == (int) ONVERTEX) {
    pointdealloc(newpt);
    *iloc = ONVERTEX;
    return splitflag;
  } else if (ivf.iloc == (int) NONCOPLANAR) {
    pointdealloc(newpt);
    *iloc = NONCOPLANAR;
    return splitflag;
  }

  if ((ivf.iloc != (int) ONFACE) && (ivf.iloc != (int) ONEDGE)) {
    terminate_tet_core(this, 2); // report a bug
  }

  // Insert the point.
  stpivot(searchsh, searchtet);
  ivf.bowywat = 3; // Use Bowyer-Watson. Preserve subsegments and subfaces;
  ivf.lawson = 2;
  ivf.rejflag = 1; // Do check the encroachment of segments.
  if (b->metric) {
    ivf.rejflag |= 4;  // Do check encroachment of protecting balls.
  }
  ivf.chkencflag = (chkencflag & (~1));
  ivf.sloc = (int) INSTAR; // ivf.iloc;
  ivf.sbowywat = 3; // ivf.bowywat;
  ivf.splitbdflag = 1;
  ivf.validflag = 1;
  ivf.respectbdflag = 1;
  ivf.assignmeshsize = b->metric;

  ivf.refineflag = 2;
  ivf.refinesh = *splitfac;

  ivf.smlenflag = useinsertradius; // Update the insertion radius.

  // Reject a near Steiner point on this subface when:
  //  - the insertion of the reject ccent is not due to mesh size (qflag).
  if (!qflag) {
    ivf.check_insert_radius = useinsertradius;
  }
  //if (is_adjacent) {
  //  ivf.parentpt = encpt; // This allows to insert a shorter edge.
  //} else {
    ivf.parentpt = NULL; // init
  //}

  if (insertpoint(newpt, &searchtet, &searchsh, NULL, &ivf)) {
    st_facref_count++;
    if (steinerleft > 0) steinerleft--;
    if (useinsertradius) {
      double rv = 0.0; // param[3]; // emin, maybe zero.
      
      if (is_adjacent) { // if (encpt != NULL) {
        // A sharp (dihedral) angle is involved.
        // Insertion radius must be > 0.
        double L = (ivf.smlen / 3.);
        // Choose the larger one between param[3] and L
        rv = (param[3] > L ? param[3] : L);
      }

      setpointinsradius(newpt, rv);
      setpoint2ppt(newpt, ivf.parentpt);
      if (smallest_insradius > ivf.smlen) {
        smallest_insradius = ivf.smlen;
      }
    }
    if (flipstack != NULL) {
      flipconstraints fc;
      fc.chkencflag = (chkencflag & (~1)); //chkencflag;
      fc.enqflag = 2;
      lawsonflip3d(&fc);
      //unflipqueue->restart();
    }

    if (later_unflip_queue->objects > b->unflip_queue_limit) {
      recoverdelaunay();
    }
    
    *iloc = ivf.iloc;
    return true;
  }

  // Point is not inserted.
  pointdealloc(newpt);

  if (ivf.iloc == (int) ENCSEGMENT) {
    // Bakup the split subface.
    ppt = (point *) &(splitfac->sh[3]);
    for (i = 0; i < 3; i++) bak_pts[i] = ppt[i];

    bool ref_segment = ((b->cdtrefine & 1) > 0); // -D1, -D3, -D5, or -D7

    if (ref_segment || qflag) {
      // Select an encroached segment and split it.
      for (i = 0; i < encseglist->objects; i++) {
        //face *paryseg = (face *) fastlookup(encseglist, i);
        badface *bf = (badface *) fastlookup(encseglist, i);
        if ((bf->ss.sh == NULL) ||
            (sorg(bf->ss) != bf->forg) ||
            (sdest(bf->ss) != bf->fdest)) continue; // Skip this segment.
        int tmp_iloc;
        if (split_segment(&(bf->ss), NULL, param, qflag, (chkencflag | 1), &tmp_iloc)) {
          // A Steiner point is inserted on an encroached segment.
          // Check if this subface is split as well.
          if ((splitfac->sh == NULL) || (splitfac->sh[3] == NULL)) {
            splitflag = true; break;
          } else {
            ppt = (point *) &(splitfac->sh[3]);
            if ((ppt[0] != bak_pts[0]) ||
                (ppt[1] != bak_pts[1]) ||
                (ppt[2] != bak_pts[2])) {
              splitflag = true; break;
            }
          }
        }
      }
    } // if (ref_segment)
    encseglist->restart();
    // Some segments may be encroached.
    if (badsubsegs->items > 0) {
      //repairencsegs(param, qflag, (chkencflag | 1));
      repairencsegs(param, 0, (chkencflag | 1)); // qflag = 0
    }
    // Check if this subface is split as well.
    if ((splitfac->sh == NULL) || (splitfac->sh[3] == NULL)) {
      splitflag = true;
    } else {
      ppt = (point *) &(splitfac->sh[3]);
      if ((ppt[0] != bak_pts[0]) ||
          (ppt[1] != bak_pts[1]) ||
          (ppt[2] != bak_pts[2])) {
        splitflag = true;
      }
    }
  } else if (ivf.iloc == (int) NEARVERTEX) {
    terminate_tet_core(this, 2); // report a bug
  }

  *iloc = ivf.iloc;
  return splitflag;
}

//============================================================================//
//                                                                            //
// repairencfacs()    Repair encroached subfaces.                             //
//                                                                            //
//============================================================================//

void TetMeshCore::repairencfacs(double *param, int qflag, int chkencflag)
{
  point encpt = NULL;
  double ccent[3], radius; //, param[6] = {0.,};
  int split_count = 0, rej_count = 0;
  //int qflag = 0;
  int i;

  bool ref_subface = ((b->cdtrefine & 2) > 0); // -D2, -D3, -D6, -D7

  // This function may be called from split_tetrahedron(). In this case, the
  //   insertion radius of the rejected circumcenter is stored in param[3].
  //   The check_subface() will return the insertion radius of the circumcenter
  //   of a bad quality subface also in param[3].
  double tet_emin = param[3];

  while (ref_subface &&
         ((badsubfacs->items > 0) || (split_subfaces_pool->items > 0))) {

    if (badsubfacs->items > 0) {
      badsubfacs->traversalinit();
      face *bface = (face *) badsubfacs->traverse();
      while (bface != NULL) {
        // A queued subface may have been deleted (split).
        if ((bface->sh != NULL) && (bface->sh[3] != NULL)) {
          // A queued subface may have been processed.
          if (smarktest2ed(*bface)) {
            sunmarktest2(*bface);
            for (i = 3; i < 6; i++) param[i] = 0.; // Clear previous values.
            if (get_subface_ccent(bface, ccent)) {
              encpt = NULL;
              if (check_enc_subface(bface, &encpt, ccent, &radius)) {
                param[3] = tet_emin; // maybe zero.
                enqueue_subface(bface, encpt, ccent, param);
              } else {
                if (check_subface(bface, ccent, radius, param)) {
                  if (tet_emin > 0) {
                    // Use the larger one.
                    param[3] = (param[3] > tet_emin ? param[3] : tet_emin);
                  }
                  enqueue_subface(bface, NULL, ccent, param);
                }
              }
            } else {
              // report a bug.
              terminate_tet_core(this, 2);
            }
          }
        }
        bface = (face *) badsubfacs->traverse();
      } // while (bface != NULL)

      badsubfacs->restart(); // clear this pool

      // check_enc_subface() may find some non-Delaunay subfaces.
      if (flippool->items > 0) {
        flipconstraints fc;
        fc.chkencflag = chkencflag;
        fc.enqflag = 2;
        lawsonflip3d(&fc);
      }
    } // if (badsubfacs->items > 0)

    if (split_subfaces_pool->items == 0) break;

    // Stop if we have used the desried number of Steiner points.
    if (steinerleft == 0) break;
    // Stop if the desried number of tetrahedra is reached.
    if ((elem_limit > 0) &&
        ((tetrahedrons->items - hullsize) > elem_limit)) break;


    badface *bf = top_subface();

    if ((bf->ss.sh != NULL) &&
        ( sorg(bf->ss) == bf->forg) &&
        (sdest(bf->ss) == bf->fdest) &&
        (sapex(bf->ss) == bf->fapex)) {
      // Try to split this subface.
      encpt = bf->noppo; // The encroaching vertex.
      for (i = 0; i < 3; i++) ccent[i] = bf->cent[i];
      for (i = 3; i < 6; i++) param[i] = bf->cent[i];
      split_count++;
      
      int iloc = (int) UNKNOWN;
      if (!split_subface(&bf->ss, encpt, ccent, param, qflag, chkencflag, &iloc)) {
        rej_count++;
        if (qflag || ((param[4] > (3. * b->minratio)) && (iloc != SHARPCORNER))) {
          // Queue a unsplit (bad quality) subface.
          badface *bt = NULL;
          unsplit_subfaces->newindex((void **) &bt);
          //bt->init();
          *bt = *bf;
        }
      }
    }
    dequeue_subface();
  } // while ((badsubfacs->items > 0) || (split_subfaces_pool->items > 0))

  if (b->verbose > 3) {
    printf("      Tried to split %d subfaces, %d were rejected.\n",
           split_count, rej_count);
  }
  param[3] = tet_emin; // Restore this value.

  if (badsubfacs->items > 0) {
    // Clean this list (due to the ref_subface flag)
    badsubfacs->traversalinit();
    face *bface = (face *) badsubfacs->traverse();
    while (bface != NULL) {
      // A queued subface may have been deleted (split).
      if ((bface->sh != NULL) && (bface->sh[3] != NULL)) {
        // A queued subface may have been processed.
        if (smarktest2ed(*bface)) {
          sunmarktest2(*bface);
        }
      }
      bface = (face *) badsubfacs->traverse();
    } // while (bface != NULL)
    badsubfacs->restart(); // clear this pool
  } // if (badsubfacs->items > 0)

  if (split_subfaces_pool->items > 0) {
    if (steinerleft == 0) {
      if (b->verbose) {
        printf("The desired number of Steiner points is reached.\n");
      }
    } else if (elem_limit > 0) {
      if (b->verbose) {
        printf("The desired number %ld of elements is reached.\n", elem_limit);
      }
    }
    split_subfaces_pool->restart(); // Clear this pool.
    unsplit_subfaces->restart();
    stack_enc_subfaces = NULL;
  }
}

//============================================================================//
//                                                                            //
// check_tetrahedron()    Check if the tet needs to be split.                 //
//                                                                            //
// "param[6]" returns the following data:                                     //
//   [0],[1],[2] - the location of the new point                              //
//   [3] - the samllest edge length ( = insertion radius)                     //
//   [4] - the radius-edge ratio                                              //
//   [5] - (optional) edge ratio                                              //
//                                                                            //
// "chktet" returns the shortest edge of this tet.                            //
//                                                                            //
//============================================================================//

bool TetMeshCore::check_tetrahedron(triface *chktet, double* param, int &qflag)
{
  point pd = (point) chktet->tet[7];
  if (pd == dummypoint) {
    return false; // Do not split a hull tet.
  }

  point pa = (point) chktet->tet[4];
  point pb = (point) chktet->tet[5];
  point pc = (point) chktet->tet[6];


  double D = orient3dexact(pa, pb, pc, pd); // =6*vol

  if (D >= 0.0) {
    // A degenerated tetrahedron.
    terminate_tet_core(this, 2);
  }

  qflag = 0; // default

  double elen[6];
  double vol = -D / 6.0;
  double emin = 0., ratio = 0.;

  // Calculate the circumcenter of this tet.
  point P = pa, Q = pb, R = pc, S = pd;

  double U[3], V[3], W[3], Z[3]; // variables.
       
  double hp = P[0]*P[0] + P[1]*P[1] + P[2]*P[2]; // - wp
  double hq = Q[0]*Q[0] + Q[1]*Q[1] + Q[2]*Q[2]; // - wq
  double hr = R[0]*R[0] + R[1]*R[1] + R[2]*R[2]; // - wr
  double hs = S[0]*S[0] + S[1]*S[1] + S[2]*S[2]; // - wr

  U[0] = hp; U[1] = P[1]; U[2] = P[2];
  V[0] = hq; V[1] = Q[1]; V[2] = Q[2];
  W[0] = hr; W[1] = R[1]; W[2] = R[2];
  Z[0] = hs; Z[1] = S[1]; Z[2] = S[2];

  double D1 = orient3d(U, V, W, Z);
        
  U[0] = P[0]; U[1] = hp; //U[2] = P[2];
  V[0] = Q[0]; V[1] = hq; //V[2] = Q[2];
  W[0] = R[0]; W[1] = hr; //W[2] = R[2];
  Z[0] = S[0]; Z[1] = hs; //Z[2] = S[2];

  double D2 = orient3d(U, V, W, Z);

  /*U[0] = P[0];*/ U[1] = P[1]; U[2] = hp;
  /*V[0] = Q[0];*/ V[1] = Q[1]; V[2] = hq;
  /*W[0] = R[0];*/ W[1] = R[1]; W[2] = hr;
  /*Z[0] = S[0];*/ Z[1] = S[1]; Z[2] = hs;

  double D3 = orient3d(U, V, W, Z);

  double DD = D * 2.;

  param[0] = D1 / DD;
  param[1] = D2 / DD;
  param[2] = D3 / DD;


  param[4] = 1.0; // default a good ratio.
  param[5] = vol;

  elen[0] = distance2(pc, pd);
  elen[1] = distance2(pd, pa);
  elen[2] = distance2(pa, pb);
  elen[3] = distance2(pb, pc);
  elen[4] = distance2(pb, pd);
  elen[5] = distance2(pa, pc);

  // Find the shortest edge.
  emin = elen[0];
  int eidx = 0;
  for (int i = 1; i < 6; i++) {
    if (emin > elen[i]) {
      emin = elen[i]; eidx = i;
    }
  }
  emin = sqrt(emin);
  // Let chktet be the shortest edge in this tet.
  chktet->ver = edge2ver[eidx];

  // check mesh size (qflag).
  if (b->varvolume || b->fixedvolume) { // -a#
    if (b->fixedvolume) {
      if (vol > b->maxvolume) {
        // set the insertion radius, use the smaller one between the
        //   smallest edge length of this tet and mesh size;
        emin = (emin < b->maxvolume_length ? emin : b->maxvolume_length);
        qflag = 1;
      }
    }
    if (!qflag && b->varvolume) {
      double volbnd = volumebound(chktet->tet);
      if ((volbnd > 0.0) && (vol > volbnd)) {
        // set the insertion radius;
        double msize = pow(volbnd, 1./3.) / 3.;
        emin = (emin < msize ? emin : msize);
        qflag = 1;
      }
    }
  } // -a#

  if (!qflag && b->metric) { // -m
    //int eidx = 0;
    for (int i = 0; i < 6; i++) {
      elen[i] = sqrt(elen[i]);
    }
    if (pa[pointmtrindex] > 0) {
      // Get the longest edge {pa, pd}, {pa, pb}, {pa, pc}
      double maxelen = elen[1]; //eidx = 1;
      if (maxelen < elen[2]) {maxelen = elen[2]; /*eidx = 2;*/}
      if (maxelen < elen[5]) {maxelen = elen[5]; /*eidx = 5;*/}
      maxelen /= 2.0;
      if (maxelen > pa[pointmtrindex]) {
        emin = (emin < pa[pointmtrindex] ? emin : pa[pointmtrindex]);
        //emax = maxelen;
        qflag = 1;
      }
    }
    if (!qflag && (pb[pointmtrindex] > 0)) {
      // Get the longest edge at pb.
      double maxelen = elen[2]; //eidx = 2;
      if (maxelen < elen[3]) {maxelen = elen[3]; /*eidx = 3;*/}
      if (maxelen < elen[4]) {maxelen = elen[4]; /*eidx = 4;*/}
      maxelen /= 2.0;
      if (maxelen > pb[pointmtrindex]) {
        emin = (emin < pb[pointmtrindex] ? emin : pb[pointmtrindex]);
        //emax = maxelen;
        qflag = 1;
      }
    }
    if (!qflag && (pc[pointmtrindex] > 0)) {
      // Get the longest edge at pc.
      double maxelen = elen[0]; //eidx = 0;
      if (maxelen < elen[3]) {maxelen = elen[3]; /*eidx = 3;*/}
      if (maxelen < elen[5]) {maxelen = elen[5]; /*eidx = 5;*/}
      maxelen /= 2.0;
      if (maxelen > pc[pointmtrindex]) {
        emin = (emin < pc[pointmtrindex] ? emin : pc[pointmtrindex]);
        //emax = maxelen;
        qflag = 1;
      }
    }
    if (!qflag && (pd[pointmtrindex] > 0)) {
      // Get the longest edge at pd.
      double maxelen = elen[0]; //eidx = 0;
      if (maxelen < elen[1]) {maxelen = elen[1]; /*eidx = 1;*/}
      if (maxelen < elen[4]) {maxelen = elen[4]; /*eidx = 4;*/}
      maxelen /= 2.0;
      if (maxelen > pd[pointmtrindex]) {
        emin = (emin < pd[pointmtrindex] ? emin : pd[pointmtrindex]);
        //emax = maxelen;
        qflag = 1;
      }
    }
  } // if (!qflag && b->metric) // -m

  if (qflag) {
    param[3] = emin; // The desired mesh size.
    //param[4] = 1.0; // ratio; // = 0.
    //param[5] = vol;
    return true;
  }

  if (b->minratio > 1.0) {
    double radius = distance(param, pa);

    ratio = radius / emin;


    if (ratio > b->minratio) {
      //qflag = 0;
      // The smallest insertion radius should be at least larger than
      //   the smallest edge length (==> graded mesh size).
      point pa = org(*chktet);
      point pb = dest(*chktet);
      double ra = getpointinsradius(pa);
      double rb = getpointinsradius(pb);
      if ((ra > 0.) && (ra > emin)) {
        emin = ra; // the relaxed (enlarged) insertion radius.
      }
      if ((rb > 0.) && (rb > emin)) {
        emin = rb; // the relaxed (enlarged) insertion radius.
      }

      param[3] = emin; // (emin * b->minratio);
      param[4] = ratio;
      //param[5] = vol;
      return true;
    }
  }

  return false; // no need to split this tetrahedron.
}

//============================================================================//
//                                                                            //
// checktet4split()    Check if a given tet has a bad shape.                  //
//                                                                            //
//============================================================================//

bool TetMeshCore::checktet4split(triface *chktet, double* param, int& qflag)
{
  point pa, pb, pc, pd, *ppt;
  double vda[3], vdb[3], vdc[3];
  double vab[3], vbc[3], vca[3];
  double N[4][3], L[4], cosd[6], elen[6];
  double maxcosd, vol, volbnd, rd, Lmax, Lmin;
  double A[4][4], rhs[4], D;
  int indx[4];
  int i, j;

  if (b->convex) { // -c
    // Skip this tet if it lies in the exterior.
    if (elemattribute(chktet->tet, numelemattrib - 1) == -1.0) {
    return 0;
    }
  }

  qflag = 0;
  for (i = 0; i < 6; i++) param[i] = 0.;

  pd = (point) chktet->tet[7];
  if (pd == dummypoint) {
    return 0; // Do not split a hull tet.
  }

  pa = (point) chktet->tet[4];
  pb = (point) chktet->tet[5];
  pc = (point) chktet->tet[6];


  // Get the edge vectors vda: d->a, vdb: d->b, vdc: d->c.
  // Set the matrix A = [vda, vdb, vdc]^T.
  for (i = 0; i < 3; i++) A[0][i] = vda[i] = pa[i] - pd[i];
  for (i = 0; i < 3; i++) A[1][i] = vdb[i] = pb[i] - pd[i];
  for (i = 0; i < 3; i++) A[2][i] = vdc[i] = pc[i] - pd[i];
  
  // Get the other edge vectors.
  for (i = 0; i < 3; i++) vab[i] = pb[i] - pa[i];
  for (i = 0; i < 3; i++) vbc[i] = pc[i] - pb[i];
  for (i = 0; i < 3; i++) vca[i] = pa[i] - pc[i];

  if (!lu_decmp(A, 3, indx, &D, 0)) {
    // Is it a degenerated tet (vol = 0).
    double D = orient3dexact(pa, pb, pc, pd); // =6*vol
    if (D >= 0.0) {
      // A degenerated tetrahedron.
      terminate_tet_core(this, 2);
    }
    // We temporarily leave this tet. It should be fixed by mesh improvement.
    return false;
  }

  // Calculate the circumcenter and radius of this tet.
  rhs[0] = 0.5 * dot(vda, vda);
  rhs[1] = 0.5 * dot(vdb, vdb);
  rhs[2] = 0.5 * dot(vdc, vdc);
  lu_solve(A, 3, indx, rhs, 0);

  for (i = 0; i < 3; i++) param[i] = pd[i] + rhs[i];
  rd = sqrt(dot(rhs, rhs));

  // Check volume if '-a#' and '-a' options are used.
  if (b->varvolume || b->fixedvolume) {
    vol = fabs(A[indx[0]][0] * A[indx[1]][1] * A[indx[2]][2]) / 6.0;
    if (b->fixedvolume) {
      if (vol > b->maxvolume) {
        qflag = 1;
      }
    }
    if (!qflag && b->varvolume) {
      volbnd = volumebound(chktet->tet);
      if ((volbnd > 0.0) && (vol > volbnd)) {
        qflag = 1;
      }
    }
    if (qflag == 1) {
      return true;
    }
  }

  if (b->metric) { // -m option. Check mesh size.
    // Check if the ccent lies outside one of the prot.balls at vertices.
    ppt = (point *) &(chktet->tet[4]);
    for (i = 0; i < 4; i++) {
      if (ppt[i][pointmtrindex] > 0) {
        if (rd > ppt[i][pointmtrindex]) {
          qflag = 1; // Enforce mesh size.
          return true;
        }
      }
    }
  }

  if (in->tetunsuitable != NULL) {
    // Execute the user-defined meshing sizing evaluation.
    if ((*(in->tetunsuitable))(pa, pb, pc, pd, NULL, 0)) {
      return true;
    }
  }


  // Check the radius-edge ratio. Set by -q#.
  if (b->minratio > 0) {
    elen[0] = dot(vdc, vdc);
    elen[1] = dot(vda, vda);
    elen[2] = dot(vab, vab);
    elen[3] = dot(vbc, vbc);
    elen[4] = dot(vdb, vdb);
    elen[5] = dot(vca, vca);

    Lmax = Lmin = elen[0];
    int eidx = 0;
    for (i = 1; i < 6; i++) {
      Lmax = (Lmax < elen[i] ? elen[i] : Lmax);
      //Lmin = (Lmin > elen[i] ? elen[i] : Lmin);
      if (Lmin > elen[i]) {
        Lmin = elen[i]; eidx = i;
      }
    }
    // Let chktet be the shortest edge in this tet.
    chktet->ver = edge2ver[eidx];

    //Lmax = sqrt(Lmax);
    Lmin = sqrt(Lmin);
    D = rd / Lmin;
    if (D > b->minratio) {
      // A bad radius-edge ratio.
      param[3] = Lmin;
      param[4] = D;
      param[5] = sqrt(Lmax) / Lmin; // edge ratio.
      return true;
    }
  } // if (b->minratio > 0)

  // Check the minimum dihedral angle. Set by -q/#.
  if (b->mindihedral > 0) {
    // Compute the 4 face normals (N[0], ..., N[3]).
    for (j = 0; j < 3; j++) {
      for (i = 0; i < 3; i++) N[j][i] = 0.0;
      N[j][j] = 1.0;  // Positive means the inside direction
      lu_solve(A, 3, indx, N[j], 0);
    }
    for (i = 0; i < 3; i++) N[3][i] = - N[0][i] - N[1][i] - N[2][i];
    // Normalize the normals.
    for (i = 0; i < 4; i++) {
      L[i] = sqrt(dot(N[i], N[i]));
      if (L[i] == 0) {
        terminate_tet_core(this, 2);
      }
      for (j = 0; j < 3; j++) N[i][j] /= L[i];
    }
    // Calculate the six dihedral angles.
    cosd[0] = -dot(N[0], N[1]); // Edge cd, bd, bc.
    cosd[1] = -dot(N[0], N[2]);
    cosd[2] = -dot(N[0], N[3]);
    cosd[3] = -dot(N[1], N[2]); // Edge ad, ac
    cosd[4] = -dot(N[1], N[3]);
    cosd[5] = -dot(N[2], N[3]); // Edge ab
    // Get the smallest dihedral angle.
    //maxcosd = mincosd = cosd[0];
    maxcosd = cosd[0];
    for (i = 1; i < 6; i++) {
      //if (cosd[i] > maxcosd) maxcosd = cosd[i];
      maxcosd = (cosd[i] > maxcosd ? cosd[i] : maxcosd);
      //mincosd = (cosd[i] < mincosd ? cosd[i] : maxcosd);
    }
    if (maxcosd > cosmindihed) {
      // A bad dihedral angle.
      return true;
    }
  } // if (b->mindihedral > 0)

  return 0;
}

//============================================================================//
//                                                                            //
// locate_point_walk()    Locate a point by line searching.                   //
//                                                                            //
//============================================================================//

enum TetMeshCore::locateresult
  TetMeshCore::locate_point_walk(point searchpt, triface* searchtet, int chkencflag)
{
  // Construct the starting point to be the barycenter of 'searchtet'.
  double startpt[3];
  point *ppt = (point *) &(searchtet->tet[4]);
  for (int i = 0; i < 3; i++) {
    startpt[i] = (ppt[0][i] + ppt[1][i] + ppt[2][i] + ppt[3][i]) / 4.;
  }

  point torg, tdest, tapex, toppo;
  double ori, oriorg, oridest, oriapex;
  enum locateresult loc = OUTSIDE;
  enum {ORGMOVE, DESTMOVE, APEXMOVE} nextmove;

  for (searchtet->ver = 0; searchtet->ver < 4; searchtet->ver++) {
    torg = org(*searchtet);
    tdest = dest(*searchtet);
    tapex = apex(*searchtet);
    ori = orient3d(torg, tdest, tapex, searchpt);
    if (ori < 0) break;
  }

  if (searchtet->ver == 4) {
    terminate_tet_core(this, 2);
  }
  int max_visited_tets = 10000; // tetrahedrons->items;

  // Walk through tetrahedra to locate the point.
  while (max_visited_tets > 0) {
    toppo = oppo(*searchtet);

    // Check if the vertex is we seek.
    if (toppo == searchpt) {
      // Adjust the origin of searchtet to be searchpt.
      esymself(*searchtet);
      eprevself(*searchtet);
      loc = ONVERTEX; // return ONVERTEX;
      break;
    }

    // We enter from the crruent face of `serarchtet', which face do we exit?
    // Find the next face which is intersect with the line (startpt->searchpt).
    oriorg  = orient3d(tdest, tapex, toppo, searchpt);
    oridest = orient3d(tapex,  torg, toppo, searchpt);
    oriapex = orient3d( torg, tdest, toppo, searchpt);

    if (oriorg < 0) {
      if (oridest < 0) {
        if (oriapex < 0) {
          // All three faces are possible.
          if (tri_edge_test(tdest,tapex,toppo,startpt,searchpt,NULL,0,NULL,NULL)) {
            nextmove = ORGMOVE;
          } else if (tri_edge_test(tapex,torg,toppo,startpt,searchpt,NULL,0,NULL,NULL)) {
            nextmove = DESTMOVE;
          } else if (tri_edge_test(torg,tdest,toppo,startpt,searchpt,NULL,0,NULL,NULL)) {
            nextmove = APEXMOVE;
          } else {
            int s = randomnation(3); // 's' is in {0,1,2}.
            if (s == 0) {
              nextmove = ORGMOVE;
            } else if (s == 1) {
              nextmove = DESTMOVE;
            } else {
              nextmove = APEXMOVE;
            }
          }
        } else {
          // Two faces, opposite to origin and destination, are viable.
          if (tri_edge_test(tdest,tapex,toppo,startpt,searchpt,NULL,0,NULL,NULL)) {
            nextmove = ORGMOVE;
          } else if (tri_edge_test(tapex,torg,toppo,startpt,searchpt,NULL,0,NULL,NULL)) {
            nextmove = DESTMOVE;
          } else {
            //s = randomnation(2); // 's' is in {0,1}.
            if (randomnation(2)) {
              nextmove = ORGMOVE;
            } else {
              nextmove = DESTMOVE;
            }
          }
        }
      } else {
        if (oriapex < 0) {
          // Two faces, opposite to origin and apex, are viable.
          if (tri_edge_test(tdest,tapex,toppo,startpt,searchpt,NULL,0,NULL,NULL)) {
            nextmove = ORGMOVE;
          } else if (tri_edge_test(torg,tdest,toppo,startpt,searchpt,NULL,0,NULL,NULL)) {
            nextmove = APEXMOVE;
          } else {
            //s = randomnation(2); // 's' is in {0,1}.
            if (randomnation(2)) {
              nextmove = ORGMOVE;
            } else {
              nextmove = APEXMOVE;
            }
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
          if (tri_edge_test(tapex,torg,toppo,startpt,searchpt,NULL,0,NULL,NULL)) {
            nextmove = DESTMOVE;
          } else if (tri_edge_test(torg,tdest,toppo,startpt,searchpt,NULL,0,NULL,NULL)) {
            nextmove = APEXMOVE;
          } else {
            //s = randomnation(2); // 's' is in {0,1}.
            if (randomnation(2)) {
              nextmove = DESTMOVE;
            } else {
              nextmove = APEXMOVE;
            }
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
    //fsymself(*searchtet);
    //if (oppo(*searchtet) == dummypoint) {
    //  loc = OUTSIDE; // return OUTSIDE;
    //  break;
    //}
    decode(searchtet->tet[searchtet->ver & 3], *searchtet); // fsymself
    if (ishulltet(*searchtet)) {
      loc = OUTSIDE; // return OUTSIDE;
      break;
    }
    max_visited_tets--;

    // Retreat the three vertices of the base face.
    torg = org(*searchtet);
    tdest = dest(*searchtet);
    tapex = apex(*searchtet);
  } // while (true)

  return loc;
}

//============================================================================//
//                                                                            //
// splittetrahedron()    Split a tetrahedron.                                 //
//                                                                            //
//============================================================================//

bool TetMeshCore::split_tetrahedron(triface* splittet, // the tet to be split.
  double *param, // param[6], it contains the following data
               // [0],[1],[2] - the location of the new point
               // [3] - the samllest edge length ( = insertion radius)
               // [4] - radius-edge ratio
               // [5] - its volume
  int qflag,   // split due to mesh size enforcement.
  int chkencflag,
  insertvertexflags &ivf)
{
  triface searchtet;
  point newpt, bak_pts[4], *ppt;
  bool splitflag = false;
  int i;


  insert_point_count++;
  if (!b->quiet && (b->refine_progress_ratio > 0.)) {
    if (insert_point_count >= report_refine_progress) {
      printf("  %ld insertions, added %ld points, %ld tetrahedra in queue.\n",
             insert_point_count - last_insertion_count,
             points->items - last_point_count,
             check_tets_list->objects);
      last_point_count = points->items; // update it.
      last_insertion_count = insert_point_count;
      // The next report event
      report_refine_progress *= (1. + b->refine_progress_ratio);
    }
  }

  makepoint(&newpt, FREEVOLVERTEX);
  for (i = 0; i < 3; i++) newpt[i] = param[i];

  // Locate the new point. Starting from an interior point 'q' of the
  //   splittet. We perform a walk from q to the 'newpt', stop walking
  //   either we hit a subface or enter OUTSIDE.
  searchtet = *splittet;
  ivf.iloc = (int) OUTSIDE;
  //ivf.iloc = locate(newpt, &searchtet, 1); // 'chkencflag' = 1.
  ivf.iloc = locate_point_walk(newpt, &searchtet, 1); // 'chkencflag' = 1.


  if ((ivf.iloc == (int) ENCSUBFACE) || (ivf.iloc == (int) OUTSIDE)) {
    // The circumcenter 'c' is not visible from 'q' (the interior of the tet).
    pointdealloc(newpt);  // Do not insert this vertex.


    ivf.iloc = (int) FENSEDIN;
    return splitflag;
  } // if (ivf.iloc == (int) ENCSUBFACE)

  // Use Bowyer-Watson algorithm. Preserve subsegments and subfaces;
  ivf.bowywat = 3;
  ivf.lawson = 2;
  ivf.rejflag = 3;  // Do check for encroached segments and subfaces.
  if (b->metric) {
    ivf.rejflag |= 4; // Reject it if it lies in some protecting balls.
  }
  ivf.chkencflag = (chkencflag & (~3)); // chkencflag;
  ivf.sloc = ivf.sbowywat = 0; // No use.
  ivf.splitbdflag = 0; // No use (its an interior vertex).
  ivf.validflag = 1;
  ivf.respectbdflag = 1;
  ivf.assignmeshsize = b->metric;

  // Mesh refinement options.
  ivf.refineflag = 1;
  ivf.refinetet = *splittet;
  // get the shortest edge length to the new point.
  ivf.smlenflag = useinsertradius;
  if (!qflag) {
    // Avoid creating an unnecessarily short edge.
    ivf.check_insert_radius = useinsertradius;
  } else {
    ivf.check_insert_radius = 0;
  }
  ivf.parentpt = NULL; // init.

  if (insertpoint(newpt, &searchtet, NULL, NULL, &ivf)) {
    // Vertex is inserted.
    st_volref_count++;
    if (steinerleft > 0) steinerleft--;
    if (useinsertradius) {
      // Save the shortest edge between: emin and ivf.smlen
      double rv = 0.0; // ivf.smlen;
      if (param[3] > 0.0) { // The smallest edge length of this tet.
        rv = (param[3] < ivf.smlen ? param[3] : ivf.smlen);
      }
      setpointinsradius(newpt, rv); // ivf.smlen
      setpoint2ppt(newpt, ivf.parentpt);
      if (ivf.smlen < smallest_insradius) { // ivf.smlen
        smallest_insradius = ivf.smlen;
      }
    }
    if (flipstack != NULL) {                 
      flipconstraints fc;
      fc.chkencflag = (chkencflag & (~3)); //chkencflag;
      fc.enqflag = 2;
      lawsonflip3d(&fc);
      //unflipqueue->restart();
    }

    if (later_unflip_queue->objects > b->unflip_queue_limit) {
      recoverdelaunay();
    }

    return true;
  } 

  // Point is not inserted.
  pointdealloc(newpt);

  if (ivf.iloc == (int) ENCSEGMENT) {
    if (!b->nobisect) { //if (!b->nobisect && qflag) { // no -Y
      // bakup the vertices of this tet.
      ppt = (point *) &(splittet->tet[4]);
      for (i = 0; i < 4; i++) bak_pts[i] = ppt[i];
    
      bool ref_segment = ((b->cdtrefine & 1) > 0);
    
      if (ref_segment || qflag) {
        for (i = 0; i < encseglist->objects; i++) {
          //face *paryseg = (face *) fastlookup(encseglist, i);
          badface *bf = (badface *) fastlookup(encseglist, i);
          if ((bf->ss.sh == NULL) ||
              (sorg(bf->ss) != bf->forg) ||
              (sdest(bf->ss) != bf->fdest)) {
            continue; // Skip this segment.
          }
          int tmp_iloc;
          if (split_segment(&(bf->ss), NULL, param, qflag, (chkencflag | 3), &tmp_iloc)) {
            // A Steienr point is inserted on a segment.
            // Check if this tet is split as well.
            if ((splittet->tet == NULL) || (splittet->tet[4] == NULL)) {
              splitflag = true; // The tet is split as well.
            } else {
              ppt = (point *) &(splittet->tet[4]);
              if ((ppt[0] != bak_pts[0]) ||
                  (ppt[1] != bak_pts[1]) ||
                  (ppt[2] != bak_pts[2]) ||
                  (ppt[3] != bak_pts[3])) {
                splitflag = true; // The tet is split as well.
              }
            }
            if (splitflag) {
              break; // This tetrahedron is split.
            }
          }
        } // i
      } // if (ref_segment ||qflag)
      encseglist->restart();
      // Some segments may need to be repaired.
      if (badsubsegs->items > 0) {
        //repairencsegs(param, qflag, (chkencflag | 3)); // Queue new enroached subsegments and subfaces.
        repairencsegs(param, 0, (chkencflag | 3)); // qflag = 0
      }
      // Some subfaces may need to be repaired.
      if (badsubfacs->items > 0) {
        //repairencfacs(param, qflag, (chkencflag | 2)); // Queue new encroached subfaces.
        repairencfacs(param, 0, (chkencflag | 2)); // qflag = 0
        if (unsplit_subfaces->objects > 0) {
          unsplit_subfaces->restart(); // clear this list;
        }
      }
      // Check if this tet is split as well.
      if ((splittet->tet == NULL) || (splittet->tet[4] == NULL)) {
        splitflag = true; // The tet is split as well.
      } else {
        ppt = (point *) &(splittet->tet[4]);
        if ((ppt[0] != bak_pts[0]) ||
            (ppt[1] != bak_pts[1]) ||
            (ppt[2] != bak_pts[2]) ||
            (ppt[3] != bak_pts[3])) {
          splitflag = true; // The tet is split as well.
        }
      }
    } else { // if (!b->nobisect) { // no -Y
      encseglist->restart();
    }
  } else if (ivf.iloc == (int) ENCSUBFACE) {
    if (!b->nobisect) { //if (!b->nobisect && qflag) { // no -Y
      // bakup the vertices of this tet.
      ppt = (point *) &(splittet->tet[4]);
      for (i = 0; i < 4; i++) bak_pts[i] = ppt[i];
    
      bool ref_subface = ((b->cdtrefine & 2) > 0);
    
      if (ref_subface || qflag) {
        // This rejected Steiner point may encroach upon more than one subfaces.
        //   We split the one which contains the projection of this rejected
        //   Steiner point. Moreover, there may be many subfaces.
        triface adjtet;
        point pa, pb, pc, toppo;
        double prjpt[3], ori;
        int scount = 0;
        int t1ver;
        
        // Clean the bad radius-edge ratio, so split_subface() knows that
        //   the split of this subface is due to a rejected tet ccenter.
        param[4] = 0.0;
    
        for (i = 0; i < encshlist->objects; i++) {
          badface *bface = (badface *) fastlookup(encshlist, i);
          // This subface may be split.
          if ((bface->ss.sh == NULL) ||
              (sorg(bface->ss) != bface->forg) ||
              (sdest(bface->ss) != bface->fdest) ||
              (sapex(bface->ss) != bface->fapex)) {
            continue;
          }
          stpivot(bface->ss, adjtet);
          if (ishulltet(adjtet)) {
            fsymself(adjtet);
          }
          toppo = oppo(adjtet); // used by orient3d()
          //assert(toppo != dummypoint);
          pa = org(adjtet);
          pb = dest(adjtet);
          pc = apex(adjtet);
          projpt2face(param, pa, pb, pc, prjpt);
          ori = orient3d(pa, pb, toppo, prjpt);
          if (ori >= 0) {
            ori = orient3d(pb, pc, toppo, prjpt);
            if (ori >= 0) {
              ori = orient3d(pc, pa, toppo, prjpt);
              if (ori >= 0) {
                scount++;
                // Found such a subface, try to split it.
                int tmp_iloc;
                split_subface(&(bface->ss), NULL, bface->cent, param, qflag,
                              chkencflag | 2, &tmp_iloc);
                // This subface may not be split while some encroached subsegments
                //   might be split.
                // Check if this tet is split as well.
                if ((splittet->tet == NULL) || (splittet->tet[4] == NULL)) {
                  splitflag = true; // The tet is split as well.
                } else {
                  ppt = (point *) &(splittet->tet[4]);
                  if ((ppt[0] != bak_pts[0]) ||
                      (ppt[1] != bak_pts[1]) ||
                      (ppt[2] != bak_pts[2]) ||
                      (ppt[3] != bak_pts[3])) {
                    splitflag = true; // The tet is split as well.
                  }
                }
                if (splitflag) {
                  break;
                }
              } // if (ori >= 0)
            }
          }
        } // i
        if (scount == 0) {
          // Not such subface is found! This can happen due to the existence
          //   of small angles and non-Delaunay elements.
          // Select an encroached subface and split it.
          for (i = 0; i < encshlist->objects; i++) {
            badface *bface = (badface *) fastlookup(encshlist, i);
            if ((bface->ss.sh == NULL) ||
                (sorg(bface->ss) != bface->forg) ||
                (sdest(bface->ss) != bface->fdest) ||
                (sapex(bface->ss) != bface->fapex)) {
              continue;
            }
            //if (get_subface_ccent(&(bface->ss), ccent)) {
            int tmp_iloc;
            split_subface(&(bface->ss), NULL, bface->cent, param, qflag,
                          chkencflag | 2, &tmp_iloc);
            // Check if this tet is split as well.
            if ((splittet->tet == NULL) || (splittet->tet[4] == NULL)) {
              splitflag = true; // The tet is split as well.
            } else {
              ppt = (point *) &(splittet->tet[4]);
              if ((ppt[0] != bak_pts[0]) ||
                  (ppt[1] != bak_pts[1]) ||
                  (ppt[2] != bak_pts[2]) ||
                  (ppt[3] != bak_pts[3])) {
                splitflag = true; // The tet is split as well.
              }
            }
            if (splitflag) {
              break; // This tetrahedron is split.
            }
          }
        } // if (scount == 0)
      } // if (ref_subface)
      encshlist->restart(); // Clear the list.
      // Some subfaces may need to be repaired.
      if (badsubfacs->items > 0) {
        //repairencfacs(param, qflag, (chkencflag | 2)); // Queue new encroached subfaces.
        repairencfacs(param, 0, (chkencflag | 2)); // qflag = 0
        if (unsplit_subfaces->objects > 0) {
          unsplit_subfaces->restart(); // clear this list.
        }
      }
      // Check if this tet is split as well.
      if ((splittet->tet == NULL) || (splittet->tet[4] == NULL)) {
        splitflag = true; // The tet is split as well.
      } else {
        ppt = (point *) &(splittet->tet[4]);
        if ((ppt[0] != bak_pts[0]) ||
            (ppt[1] != bak_pts[1]) ||
            (ppt[2] != bak_pts[2]) ||
            (ppt[3] != bak_pts[3])) {
          splitflag = true; // The tet is split as well.
        }
      }
    } else { // if (!b->nobisect)
      encshlist->restart();
    }
  }

  return splitflag;
}

//============================================================================//
//                                                                            //
// repairbadtets()    Repair bad quality tetrahedra.                          //
//                                                                            //
//============================================================================//

void TetMeshCore::repairbadtets(double queratio, int chkencflag)
{
  triface *bface, *quetet, *last_quetet;
  triface checktet;
  double param[6] = {0.,};
  int qflag = 0;
  int i;

  while ((badtetrahedrons->items > 0) || (check_tets_list->objects > 0)) {

    if (badtetrahedrons->items > 0) {
      badtetrahedrons->traversalinit();
      bface = (triface *) badtetrahedrons->traverse();
      while (bface != NULL) {
        check_tets_list->newindex((void **) &quetet);
        *quetet = *bface;
        bface = (triface *) badtetrahedrons->traverse();
      }
      badtetrahedrons->restart();
    }
    
    // Stop if we have used the desried number of Steiner points.
    if (steinerleft == 0) break;
    // Stop if the desried number of tetrahedra is reached.
    if ((elem_limit > 0) &&
        ((tetrahedrons->items - hullsize) > elem_limit)) break;


    // Randomly select a tet to split.
    i = rand() % check_tets_list->objects;
    quetet = (triface *) fastlookup(check_tets_list, i);
    checktet = *quetet;
    
    // Fill the current position by the last tet in the list.
    i = check_tets_list->objects - 1;
    last_quetet = (triface *) fastlookup(check_tets_list, i);
    *quetet = *last_quetet;
    check_tets_list->objects--;

    if (!isdeadtet(checktet)) {
      if (marktest2ed(checktet)) {
        unmarktest2(checktet);
        //if (check_tetrahedron(&checktet, param, qflag)) {
        if (checktet4split(&checktet, param, qflag)) {
          bool splitflag = false;
          insertvertexflags ivf;
          splitflag = split_tetrahedron(&checktet, param, qflag, chkencflag, ivf);
          if (!splitflag) {
            if (qflag || (param[4] > queratio)) { // radius-edge ratio
              badface *bt = NULL;
              unsplit_badtets->newindex((void **) &bt);
              bt->init();
              bt->tt = checktet;
              bt->forg  =  org(checktet);
              bt->fdest = dest(checktet);
              bt->fapex = apex(checktet);
              bt->foppo = oppo(checktet);
              for (i = 0; i < 6; i++) bt->cent[i] = param[i];
              bt->key = (double) qflag;
            }
          }
        }
      }
    } // if (!isdeadtet(checktet)) {

  } // while ((badtetrahedrons->items > 0) || (check_tets_list->objects > 0))

  if (check_tets_list->objects > 0) {
    if (steinerleft == 0) {
      if (b->verbose) {
        printf("The desired number of Steiner points is reached.\n");
      }
    } else if (elem_limit > 0) {
      if (b->verbose) {
        printf("The desired number %ld of elements is reached.\n", elem_limit);
      }
    }
    //split_tets_pool->restart(); // Clear this pool.
    // Unmark all unchecked tetrahedra.
    for (i = 0; i < check_tets_list->objects; i++) {
      quetet = (triface *) fastlookup(check_tets_list, i);
      if (!isdeadtet(*quetet)) {
        unmarktest2(*quetet);
      }
    }
    check_tets_list->restart();
  }
}

//============================================================================//
//                                                                            //
// delaunayrefinement()    Refine the mesh by Delaunay refinement.            //
//                                                                            //
//============================================================================//

void TetMeshCore::delaunayrefinement()
{
  triface checktet;
  face checksh;
  face checkseg;
  long steinercount;
  double param[6] = {0., 0., 0., 0., 0., 0.};
  int qflag = 0;
  int chkencflag = 0;
  int i;

  long bak_segref_count, bak_facref_count, bak_volref_count;

  if (!b->quiet) {
    printf("Refining mesh...\n");
  }

  if (b->verbose) {
    printf("  Min radius-edge ratio = %g.\n", b->minratio);
    if (b->mindihedral > 0.) {
      printf("  Min dihedral   angle  = %g.\n", b->mindihedral);
    }
    if (b->fixedvolume) {
      printf("  Max tet volume  = %g.\n", b->maxvolume);
    }
    //printf("  Min Edge length = %g.\n", b->minedgelength);
  }
  // Used in locate_point_on_surface();
  cos_facet_separate_ang_tol = cos(b->facet_separate_ang_tol/180.*PI); // -p/#
  // Used in function is_collinear_at(mid, left, right);
  cos_collinear_ang_tol = cos(b->collinear_ang_tol/180.*PI); // -p///#

  // The cosine value of the min dihedral angle (-q/#) for tetrahedra.
  cosmindihed = cos(b->mindihedral / 180.0 * PI);

  steinerleft = b->steinerleft;  // Upperbound of # Steiner points (by -S#).
  if (steinerleft > 0) {
    // Check if we've already used up the given number of Steiner points.
    steinercount = st_segref_count + st_facref_count + st_volref_count;
    if (steinercount < steinerleft) {
      steinerleft -= steinercount;
    } else {
      if (!b->quiet) {
        printf("\nWarning:  ");
        printf("The desired number of Steiner points (%d) has reached.\n\n",
               b->steinerleft);
      }
      return; // No more Steiner points.
    }
  }

  if (b->refine && (b->elem_growth_ratio > 0.0)) { // -r#
    int ntet = in->numberoftetrahedra; // tetrahedrons->items - hullsize;
    elem_limit = ntet * (1.0 + b->elem_growth_ratio);
  }

  if (b->refine_progress_ratio > 0) { // -r/# default is 0.333
    insert_point_count = 0l;
    last_insertion_count = 0l;
    last_point_count = points->items;
    report_refine_progress = points->items * (1. + b->refine_progress_ratio);
  }

  if (!b->nobisect) { // no -Y.
    if (segmentendpointslist == NULL) {
      makesegmentendpointsmap(); // create ridge_vertex-to-segment map.
    }
    create_segment_info_list();
    makefacetverticesmap();      // create ridge_vertex-to-facet map.
    create_segment_facet_map();  // vreate segment-to-facet map.
  }


  // Begin of memory allocation ===============================================
  // Initialize the pools and priority queues.
  long bls = b->shellfaceperblock;
  long blt = b->tetrahedraperblock;

  badsubsegs = new memorypool(sizeof(face), 256, sizeof(void *), 0);
  badsubfacs = new memorypool(sizeof(face), 256, sizeof(void *), 0);
  badtetrahedrons = new memorypool(sizeof(triface), blt, sizeof(void *), 0);
  
  split_segments_pool = new memorypool(sizeof(badface), bls, sizeof(void *), 0);
  split_subfaces_pool = new memorypool(sizeof(badface), bls, sizeof(void *), 0);

  long est_size = blt;
  int log2objperblk = 0;
  while (est_size >>= 1) log2objperblk++;
  if (log2objperblk < 10) log2objperblk = 10; // At least 1024.

  check_tets_list = new arraypool(sizeof(triface), log2objperblk);
  
  unsplit_segments = new arraypool(sizeof(badface), 10);
  unsplit_subfaces = new arraypool(sizeof(badface), 10);
  unsplit_badtets = new arraypool(sizeof(badface), 10);
  
  stack_enc_segments = stack_enc_subfaces = NULL;

  for (i = 0; i < 64; i++) {
    queuefront[i] = NULL;
  }
  firstnonemptyq = -1;
  recentq = -1;

  encseglist = new arraypool(sizeof(badface), 8);
  encshlist  = new arraypool(sizeof(badface), 8);
  // End of memory allocation =================================================


  // with -r and an .elem file ================================================
  if (b->refine && (in->refine_elem_list != NULL)) {
    if (b->verbose) {
      printf("  Refining a list of given elements.\n");
    }
    //assert(b->varvolume > 0); // -a option must be used.
    chkencflag = 4; // Check bad tetrahedra.
    steinercount = points->items;

    double queratio = b->minratio > 2. ? b->minratio : 2.0;
    queratio *= 2.0; // queratio; // increase this value.

    // Create a map from indices to points.
    point *idx2verlist;
    makeindex2pointmap(idx2verlist);

    int *elelist = in->refine_elem_list;
    int elem;

    for (elem = 0; elem < in->numberofrefineelems; elem++) {
      point p1 = idx2verlist[elelist[elem*4]];
      point p2 = idx2verlist[elelist[elem*4+1]];
      point p3 = idx2verlist[elelist[elem*4+2]];
      point p4 = idx2verlist[elelist[elem*4+3]];
      
      if (!get_tet(p1, p2, p3, p4, &checktet)) {
        continue;
      }

      double volume_limit;
      if (in->refine_elem_vol_list != NULL) {
        volume_limit = in->refine_elem_vol_list[i];
      } else {
        point *ppt = (point *) &(checktet.tet[4]);
        double volume = orient3dfast(ppt[1], ppt[0], ppt[2], ppt[3]) / 6.;
        volume_limit = volume / 3.;
      }
      setvolumebound(checktet.tet, volume_limit);

      //assert(check_tets_list->objects == 0l);
      triface *quetet;
      marktest2(checktet);
      check_tets_list->newindex((void **) &quetet);
      *quetet = checktet;
      
      int maxiter = 2, iter;

      for (iter = 0; iter < maxiter; iter++) {
        repairbadtets(queratio, chkencflag);

        if (later_unflip_queue->objects > 0l) {
          recoverdelaunay();
        }

        // Split unsplit tetrahedra
        long badtetcount = 0, splitcount = 0;
        int j;

        for (i = 0; i < unsplit_badtets->objects; i++) {
          badface *bt = (badface *) fastlookup(unsplit_badtets, i);
          if ((bt->tt.tet != NULL) &&
              ( org(bt->tt) == bt->forg ) &&
              (dest(bt->tt) == bt->fdest) &&
              (apex(bt->tt) == bt->fapex) &&
              (oppo(bt->tt) == bt->foppo)) {

            if (steinerleft == 0) break;
            if (elem_limit > 0) {
              if ((tetrahedrons->items - hullsize) > elem_limit) {
                break;
              }
            }

            // Count a live tet.
            badtetcount++;
            insertvertexflags ivf;
            qflag = (int) bt->key;
            point *ppt = (point *) &(bt->tt.tet[4]);
            for (j = 0; j < 3; j++) {
              param[j] = (ppt[0][j]+ppt[1][j]+ppt[2][j]+ppt[3][j]) / 4.0;
            }
            for (; j < 6; j++) {
              param[j] = bt->cent[j];
            }
            if (split_tetrahedron(&bt->tt, param, qflag, chkencflag, ivf)) {
              splitcount++;
            }

            if (badtetrahedrons->items > 0) {
              // Push new bad quality tetrahedron into queue.
              badtetrahedrons->traversalinit();
              triface *bface = (triface *) badtetrahedrons->traverse();
              while (bface != NULL) {
                check_tets_list->newindex((void **) &quetet);
                *quetet = *bface;
                bface = (triface *) badtetrahedrons->traverse();
              }
              badtetrahedrons->restart();
            }
          }
        } // i

        unsplit_badtets->restart();

        if (splitcount == 0) break;
      } // iter

      if (check_tets_list->objects > 0) {
        // Clean the list.
        for (i = 0; i < check_tets_list->objects; i++) {
          quetet = (triface *) fastlookup(check_tets_list, i);
          if (!isdeadtet(*quetet)) {
            unmarktest2(*quetet);
          }
        }
        check_tets_list->restart();
      }

      if (steinerleft == 0) break;
      if (elem_limit > 0) {
        if ((tetrahedrons->items - hullsize) > elem_limit) {
          break;
        }
      }

    } // elem

    if (b->verbose) {
      printf("  Added %ld Steiner points.\n", points->items - steinercount);
    }
    delete [] idx2verlist;
  } // if (b->refine && (in->refine_elem_list != NULL))
  // with -r and an .elem file ================================================

  bool force_quit_refinement = false;

  if (steinerleft == 0) {
    force_quit_refinement = true;
  } else if (elem_limit > 0) {
    if ((tetrahedrons->items - hullsize) > elem_limit) {
      force_quit_refinement = true;
    }
  }

  if (!b->nobisect) { // no -Y
    bool ref_segment = ((b->cdtrefine & 1) > 0); // -D1, -D3, -D5, or -D7

    if (ref_segment && !force_quit_refinement) {
      if (b->verbose) {
        printf("  Splitting encroached subsegments.\n");
      }

      chkencflag = 1; // Only check encroaching subsegments.
      steinercount = points->items;
    
      // Add all segments into the pool.
      subsegs->traversalinit();
      checkseg.sh = shellfacetraverse(subsegs);
      while (checkseg.sh != (shellface *) NULL) {
        //enqueuesubface(badsubsegs, &checkseg);
        point encpt = NULL;
        if (check_enc_segment(&checkseg, &encpt)) {
          badface *bf = (badface *) split_segments_pool->alloc();
          bf->init();
          bf->ss = checkseg;
          bf->forg = sorg(checkseg);
          bf->fdest = sdest(checkseg);
          bf->noppo = encpt;
          // Push it into stack.
          bf->nextitem = stack_enc_segments;
          stack_enc_segments = bf;
        }
        checkseg.sh = shellfacetraverse(subsegs);
      }

      // Split all encroached segments.
      for (i = 0; i < 6; i++) param[i] = 0.0;
      qflag = 0;

      repairencsegs(param, qflag, chkencflag);

      if (b->verbose) {
        printf("  Added %ld Steiner points.\n", points->items - steinercount);
      }

      if (later_unflip_queue->objects > 0l) {
        recoverdelaunay();
      }
    } // if (ref_segment)

    bool ref_surface = ((b->cdtrefine & 2) > 0); // -D2, -D3, or -D7

    if (ref_surface && !force_quit_refinement) {
      if (b->verbose) {
        printf("  Splitting encroached and bad quality subfaces.\n");
      }

      chkencflag = 2; // only check encroaching subfaces.
      steinercount = points->items;
      bak_segref_count = st_segref_count;
      bak_facref_count = st_facref_count;

      // Add all subfaces into the pool.
      double ccent[3], radius;
      point encpt = NULL;

      subfaces->traversalinit();
      checksh.sh = shellfacetraverse(subfaces);
      while (checksh.sh != (shellface *) NULL) {
        //enqueuesubface(badsubfacs, &checksh);
        if (get_subface_ccent(&checksh, ccent)) {
          encpt = NULL;
          for (i = 3; i < 6; i++) param[i] = 0.0;
          if (check_enc_subface(&checksh, &encpt, ccent, &radius)) {
            enqueue_subface(&checksh, encpt, ccent, param);
          } else {
            if (check_subface(&checksh, ccent, radius, param)) {
              enqueue_subface(&checksh, NULL, ccent, param);
            }
          }
        } else {
          terminate_tet_core(this, 2); // report a bug.
        }
        checksh.sh = shellfacetraverse(subfaces);
      }

      // check_enc_subface() may find some non-Delaunay faces.
      if (flippool->items > 0) {
        flipconstraints fc;
        fc.chkencflag = chkencflag;
        fc.enqflag = 2;
        lawsonflip3d(&fc);
      }

      // Split all encroached subfaces.
      for (i = 0; i < 6; i++) param[i] = 0.0;
      qflag = 0;
    
      int maxiter = 3, iter;
    
      for (iter = 0; iter < maxiter; iter++) {

        if (b->verbose > 1) {
          printf("  iter = %d\n", iter+1);
        }
        long iter_steinercount = points->items;
        long iter_segref_count = st_segref_count;
        long iter_facref_count = st_facref_count;

        repairencfacs(param, qflag, chkencflag);

        if (b->verbose > 1) {
          printf("  Added %ld (%ld,%ld) Steiner points.\n",
                 points->items-iter_steinercount,
                 st_segref_count-iter_segref_count,
                 st_facref_count-iter_facref_count);
        }

        if (unsplit_subfaces->objects > 0) {
          if (b->verbose > 1) {
            printf("  splitting %ld unsplit subfaces\n", unsplit_subfaces->objects);
          }
          int scount = 0; // Count the split subfaces.
          
          for (i = 0; i < unsplit_subfaces->objects; i++) {
            badface *bf = (badface *) fastlookup(unsplit_subfaces, i);
            if ((bf->ss.sh != NULL) &&
                ( sorg(bf->ss) == bf->forg) &&
                (sdest(bf->ss) == bf->fdest) &&
                (sapex(bf->ss) == bf->fapex)) {
              // Try to split it in its barycenter.
              int iloc, j;
              for (j = 0; j < 3; j++) {
                ccent[j] = (bf->forg[j] + bf->fdest[j] + bf->fapex[j]) / 3.;
              }
              for (j = 3; j < 6; j++) param[j] = bf->cent[j];
              encpt = bf->noppo; // The encroaching vertex.
              if (split_subface(&bf->ss, encpt, ccent, param, qflag, chkencflag, &iloc)) {
                scount++;
              }
            }
          } // i

          unsplit_subfaces->restart();

          if (b->verbose > 1) {
            printf("  Split %d subfaces.\n", scount);
          }
        } else {
          break; // no unsplit subfaces.
        } // if (unsplit_subfaces->objects > 0)
      } // iter

      if (b->verbose) {
        printf("  Added %ld (%ld,%ld) Steiner points.\n",
               points->items-steinercount, st_segref_count-bak_segref_count,
               st_facref_count-bak_facref_count);
      }

      if (badsubfacs->items > 0) {
        // Clean this pool.
        badsubfacs->traversalinit();
        face *bface = (face *) badsubfacs->traverse();
        while (bface != NULL) {
          if ((bface->sh != NULL) && (bface->sh[3] != NULL)) {
            if (smarktest2ed(*bface)) {
              sunmarktest2(*bface);
            }
          }
          bface = (face *) badsubfacs->traverse();
        }
        badsubfacs->restart();
      }

      if (later_unflip_queue->objects > 0l) {
        recoverdelaunay();
      }
    } // if (ref_subface)
  } // if (!b->nobisect)

  if (((b->cdtrefine & 4) > 0) && !force_quit_refinement) { // -D4, -D5, or -D7
    // Begin of adding Steiner points in volume ===============================

    // A Steiner point can be added only if it does not encroach upon any
    //   boundary segment or subface.
    if (b->verbose) {
      printf("  Splitting bad quality tets.\n");
    }

    for (i = 0; i < 6; i++) param[i] = 0.0;
    qflag = 0;

    // Add all tetrahedra (no hull tets) into the pool.
    triface *quetet;
    tetrahedrons->traversalinit();
    checktet.tet = tetrahedrontraverse();
    while (checktet.tet != NULL) {
      marktest2(checktet);
      check_tets_list->newindex((void **) &quetet);
      *quetet = checktet;
      checktet.tet = tetrahedrontraverse();
    }


    chkencflag = 4; // Check bad tetrahedra.

    double queratio = b->minratio > 2. ? b->minratio : 2.0;
    queratio *= 2.0; // queratio; // increase this value.

    int maxiter = 3, iter;

    for (iter = 0; iter < maxiter; iter++) {
      steinercount = points->items;
      bak_segref_count = st_segref_count;
      bak_facref_count = st_facref_count;
      bak_volref_count = st_volref_count;

      if (b->verbose > 1) {
        printf("    iter = %d: queratio = %g\n", iter + 1, queratio);
      }
    
      // Split all bad quality tetrahedra.
      repairbadtets(queratio, chkencflag);

      if (b->verbose) {
        printf("  Added %ld (%ld,%ld,%ld) Steiner points.\n",
               points->items - steinercount,
               st_segref_count - bak_segref_count,
               st_facref_count - bak_facref_count,
               st_volref_count - bak_volref_count);
      }


      if (later_unflip_queue->objects > 0l) {
        recoverdelaunay();
      }

      if (unsplit_badtets->objects == 0) break;

      //queratio *= 2.0; // queratio; // increase this value.

      // Split unsplit tetrahedra
      long badtetcount = 0, splitcount = 0;
      int j;

      for (i = 0; i < unsplit_badtets->objects; i++) {
        badface *bt = (badface *) fastlookup(unsplit_badtets, i);
        if ((bt->tt.tet != NULL) &&
            ( org(bt->tt) == bt->forg ) &&
            (dest(bt->tt) == bt->fdest) &&
            (apex(bt->tt) == bt->fapex) &&
            (oppo(bt->tt) == bt->foppo)) {
          if (steinerleft == 0) break;
          if (elem_limit > 0) {
            if ((tetrahedrons->items - hullsize) > elem_limit) {
              break;
            }
          }

          // Count a live tet.
          badtetcount++;
          insertvertexflags ivf;
          qflag = (int) bt->key;
          point *ppt = (point *) &(bt->tt.tet[4]);
          for (j = 0; j < 3; j++) {
            param[j] = (ppt[0][j]+ppt[1][j]+ppt[2][j]+ppt[3][j]) / 4.0;
          }
          for (; j < 6; j++) {
            param[j] = bt->cent[j];
          }
          if (split_tetrahedron(&bt->tt, param, qflag, chkencflag, ivf)) {
            splitcount++;
          }
          if (badtetrahedrons->items > 0) {
            // Push new bad quality tetrahedron into queue.
            badtetrahedrons->traversalinit();
            triface *bface = (triface *) badtetrahedrons->traverse();
            while (bface != NULL) {
              check_tets_list->newindex((void **) &quetet);
              *quetet = *bface;
              bface = (triface *) badtetrahedrons->traverse();
            }
            badtetrahedrons->restart();
          }
        }
      } // i
      unsplit_badtets->restart();


      if (splitcount == 0) break;
    } // iter

    if (check_tets_list->objects > 0) {
      // Unmark all unchecked tetrahedra.
      for (i = 0; i < check_tets_list->objects; i++) {
        quetet = (triface *) fastlookup(check_tets_list, i);
        if (!isdeadtet(*quetet)) {
          unmarktest2(*quetet);
        }
      }
      check_tets_list->restart();
    }

    // End of adding Steiner points in volume =================================
  } // if ((b->cdtrefine & 4) > 0) {


  delete encseglist;
  delete encshlist;
  encseglist = NULL;
  encshlist = NULL;

  totalworkmemory += (badsubsegs->maxitems * badsubsegs->itembytes);
  delete badsubsegs;
  badsubsegs = NULL;
  totalworkmemory += (badsubfacs->maxitems * badsubfacs->itembytes);
  delete badsubfacs;
  badsubfacs = NULL;
  totalworkmemory += (split_subfaces_pool->maxitems * split_subfaces_pool->itembytes);
  delete split_subfaces_pool;
  split_subfaces_pool = NULL;
  delete split_segments_pool;
  split_segments_pool = NULL;
  delete unsplit_segments;
  delete unsplit_subfaces;
  unsplit_segments = unsplit_subfaces = NULL;

  totalworkmemory += (badtetrahedrons->maxitems*badtetrahedrons->itembytes);
  delete badtetrahedrons;
  badtetrahedrons = NULL;
  totalworkmemory += (check_tets_list->totalmemory);
  delete check_tets_list;
  check_tets_list = NULL;
  delete unsplit_badtets;
  unsplit_badtets = NULL;
}

//                                                                            //
//                                                                            //
//== refine_cxx ==============================================================//

} // namespace sqmesh::mesh::tet::detail
