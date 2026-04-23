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

//== steiner_cxx =============================================================//
//                                                                            //
//                                                                            //

void TetMeshCore::sort_2pts(point p1, point p2, point ppt[2])
{
  if (pointmark(p1) < pointmark(p2)) {
    ppt[0] = p1;
    ppt[1] = p2;
  } else {
    ppt[0] = p2;
    ppt[1] = p1;
  }
}

void TetMeshCore::sort_3pts(point p1, point p2, point p3, point ppt[3])
{
  int i1 = pointmark(p1);
  int i2 = pointmark(p2);
  int i3 = pointmark(p3);

  if (i1 < i2) {
    if (i1 < i3) {
      ppt[0] = p1;
      if (i2 < i3) {
        ppt[1] = p2;
        ppt[2] = p3;
      } else {
        ppt[1] = p3;
        ppt[2] = p2;
      }
    } else {
      ppt[0] = p3;
      ppt[1] = p1;
      ppt[2] = p2;
    }
  } else { // i1 > i2
    if (i2 < i3) {
      ppt[0] = p2;
      if (i1 < i3) {
        ppt[1] = p1;
        ppt[2] = p3;
      } else {
        ppt[1] = p3;
        ppt[2] = p1;
      }
    } else {
      ppt[0] = p3;
      ppt[1] = p2;
      ppt[2] = p1;
    }
  }
}


//============================================================================//
//                                                                            //
// is_collinear_at()    Check if three vertices (from left to right): left,   //
//                      mid, and right are collinear.                         //
//                                                                            //
//============================================================================//

bool TetMeshCore::is_collinear_at(point mid, point left, point right)
{
  double v1[3], v2[3];
  
  v1[0] =  left[0] - mid[0];
  v1[1] =  left[1] - mid[1];
  v1[2] =  left[2] - mid[2];
  
  v2[0] = right[0] - mid[0];
  v2[1] = right[1] - mid[1];
  v2[2] = right[2] - mid[2];
  
  double L1 = sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);
  double L2 = sqrt(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]);
  double D = (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);

  double cos_ang = D / (L1 * L2);
  return cos_ang < cos_collinear_ang_tol;
}

//============================================================================//
//                                                                            //
// is_segment()    Check if the two vertices are endpoints of a segment.      //
//                                                                            //
//============================================================================//

bool TetMeshCore::is_segment(point p1, point p2)
{
  if (pointtype(p1) == RIDGEVERTEX) {
    if (pointtype(p2) == RIDGEVERTEX) {
      // Check if p2 is connect to p1.
      int idx = pointmark(p1);
      for (int i = idx_segment_ridge_vertex_list[idx];
             i < idx_segment_ridge_vertex_list[idx+1]; i++) {
        if (segment_ridge_vertex_list[i] == p2) {
          return true;
        }
      }
    } else if (pointtype(p2) == FREESEGVERTEX) {
      // Check if the segment contains p2 has one if its endpoints be p1.
      face parsentseg;
      sdecode(point2sh(p2), parsentseg);
      int segidx = getfacetindex(parsentseg);
      if ((segmentendpointslist[segidx*2] == p1) ||
          (segmentendpointslist[segidx*2+1] == p1)) {
        return true;
      }
    }
  } else {
    if (pointtype(p1) == FREESEGVERTEX) {
      if (pointtype(p2) == RIDGEVERTEX) {
        face parsentseg;
        sdecode(point2sh(p1), parsentseg);
        int segidx = getfacetindex(parsentseg);
        if ((segmentendpointslist[segidx*2] == p2) ||
            (segmentendpointslist[segidx*2+1] == p2)) {
          return true;
        }
      } else if (pointtype(p2) == FREESEGVERTEX) {
        face parsentseg1, parsentseg2;
        sdecode(point2sh(p1), parsentseg1);
        sdecode(point2sh(p2), parsentseg2);
        int segidx1 = getfacetindex(parsentseg1);
        int segidx2 = getfacetindex(parsentseg2);
        if (segidx1 == segidx2) {
          return true;
        }
      }
    }
  }
  
  return false;
}

//============================================================================//
//                                                                            //
// valid_constrained_f23()    Validate a 2-3 flip.                            //
//                                                                            //
// The purpose of the following check is to avoid creating a degenrated face  //
//   (and subface) whose three vertices are nearly on one segment or on two   //
//   nearly collinear segments.                                               //
//                                                                            //
// "checktet" is a face (a,b,c) which is 2-3 flippable, and (d,e) will be     //
//   the new edge after this flip.                                            //
//                                                                            //
// return true if this 2-3 flip is good, otherwise, return false.             //
//                                                                            //
//============================================================================//

bool TetMeshCore::valid_constrained_f23(triface& checktet, point pd, point pe)
{
  bool validflag = true;

  triface spintet;
  face checkseg1, checkseg2;
  point checkpt;

  for (int k = 0; k < 3; k++) {
    checkpt = org(checktet);
    esym(checktet, spintet);
    enextself(spintet); // [x, d], x = a,b,c
    tsspivot1(spintet, checkseg1);
    bool isseg = (checkseg1.sh != NULL);
    if (!isseg && boundary_recovery_flag) {
      isseg = is_segment(checkpt, pd);
    }
    if (isseg) {
      fsym(checktet, spintet);
      esymself(spintet);
      eprevself(spintet);
      tsspivot1(spintet, checkseg2);
      isseg = (checkseg2.sh != NULL);
      if (!isseg && boundary_recovery_flag) {
        isseg = is_segment(checkpt, pe);
      }
      if (isseg) {
        if (pointtype(checkpt) == FREESEGVERTEX) {
          // In this case, the two subsegments (checkseg1, checkseg2)
          //   must belong to the same segment, do not flip.
          validflag = false;
          break;
        } else {
          // Check if three vertices are nearly collinear. The middle
          //   vertex is checkpt.
          if ((checkpt != dummypoint) &&
              (pe != dummypoint) &&
              (pd != dummypoint)) {
            if (is_collinear_at(checkpt, pe, pd)) {
              validflag = false;
              break;
            }
          }
        }
      } // if (isseg)
    } // if (isseg)
    enextself(checktet);
  } // k

  return validflag;
}

//============================================================================//
//                                                                            //
// valid_constrained_f32()    Validate a 3-2 flip.                            //
//                                                                            //
// Avoid creating a degenerated tetrahedral face whose three vertices are on  //
//   one (sub)segment. abtets[0], abdtets[1], abtets[2] are three tets        //
//   at the flipping edge (a,b), the new face will be (c, d, e).              //
//  The only new face we will create is (c,d,e), make sure that it is not     //
//  a (nearly) degenerated face. If the vertex c is RIDGEVEETEX or            //
//  FREESEGVERTEX, then the edges (c, d) and (c, e) should not on one segment.//
//  The same for the vertex d and e.                                          //
//                                                                            //
// return true if this 3-2 flip is good, otherwise, return false.             //
//                                                                            //
//============================================================================//

bool TetMeshCore::valid_constrained_f32(triface* abtets, point pa, point pb)
{
  bool validflag = true; // default.

  triface spintet;
  face checksegs[3]; // edges: [c,d], [d,e], and [e,c]
  point chkpt, leftpt, rightpt;

  // Check edges [c,d], [d,e], and [e,c]
  for (int k = 0; k < 3; k++) { // [a,b,c], [a,b,d], [a,b,e]
    enext(abtets[k], spintet);
    esymself(spintet);
    eprevself(spintet); // [c,d], [d,e], and [e,c]
    tsspivot1(spintet, checksegs[k]);
    // Ignore a temporaray segment (used in recoversubfaces()).
    if (checksegs[k].sh != NULL) {
      if (smarktest2ed(checksegs[k])) {
        checksegs[k].sh = NULL;
      }
    }
  } // k

  for (int k = 0; k < 3; k++) {
    chkpt   = apex(abtets[k]);         // pc
    leftpt  = apex(abtets[(k+2)%3]);   // pe
    rightpt = apex(abtets[(k+1)%3]);   // pd
    bool isseg = (checksegs[k].sh != NULL); // [c,d]
    if (!isseg && boundary_recovery_flag) {
      isseg = is_segment(chkpt, rightpt);
    }
    if (isseg) {
      isseg = (checksegs[(k+2)%3].sh != NULL); // [e,c]
      if (!isseg && boundary_recovery_flag) {
        isseg = is_segment(chkpt, leftpt);
      }
      if (isseg) {
        if (pointtype(chkpt) == FREESEGVERTEX) {
          validflag = false;
          break;
        } else {
          if ((chkpt != dummypoint) &&
              (leftpt != dummypoint) &&
              (rightpt != dummypoint)) {
            if (is_collinear_at(chkpt, leftpt, rightpt)) {
              validflag = false;
              break;
            }
          }
        }
      }
    }
  } // k

  return validflag;
}

//============================================================================//
//                                                                            //
// checkflipeligibility()    A call back function for boundary recovery.      //
//                                                                            //
// 'fliptype' indicates which elementary flip will be performed: 1 : 2-to-3,  //
// and 2 : 3-to-2, respectively.                                              //
//                                                                            //
// 'pa, ..., pe' are the vertices involved in this flip, where [a,b,c] is     //
// the flip face, and [d,e] is the flip edge. NOTE: 'pc' may be 'dummypoint', //
// other points must not be 'dummypoint'.                                     //
//                                                                            //
//============================================================================//

int TetMeshCore::checkflipeligibility(int fliptype, point pa, point pb, 
                                     point pc, point pd, point pe,
                                     int level, int edgepivot,
                                     flipconstraints* fc)
{
  point tmppts[3];
  enum interresult dir;
  int types[2], poss[4];
  int intflag;
  int rejflag = 0;
  int i;

  if (fc->seg[0] != NULL) {
    // A constraining edge is given (e.g., for edge recovery).
    if (fliptype == 1) {
      // A 2-to-3 flip: [a,b,c] => [e,d,a], [e,d,b], [e,d,c].
      tmppts[0] = pa;
      tmppts[1] = pb;
      tmppts[2] = pc;
      for (i = 0; i < 3 && !rejflag; i++) {
        if (tmppts[i] != dummypoint) {
          // Test if the face [e,d,#] intersects the edge.
          intflag = tri_edge_test(pe, pd, tmppts[i], fc->seg[0], fc->seg[1], 
                                  NULL, 1, types, poss);
          if (intflag == 2) {
            // They intersect at a single point.
            dir = (enum interresult) types[0];
            if (dir == ACROSSFACE) {
              // The interior of [e,d,#] intersect the segment.
              rejflag = 1;
            } else if (dir == ACROSSEDGE) {
              if (poss[0] == 0) {
                // The interior of [e,d] intersect the segment.
                // Since [e,d] is the newly created edge. Reject this flip.
                rejflag = 1; 
              }
            }
            else { 
              if ((dir == ACROSSVERT) || (dir == TOUCHEDGE) || 
                  (dir == TOUCHFACE)) {
                // should be a self-intersection.
                rejflag = 1;
              }
            } // dir
          } else if (intflag == 4) {
            // They may intersect at either a point or a line segment.
            dir = (enum interresult) types[0];
            if (dir == ACROSSEDGE) {
              if (poss[0] == 0) {
                // The interior of [e,d] intersect the segment.
                // Since [e,d] is the newly created edge. Reject this flip.
                rejflag = 1;
              }
            }
            else if (dir == ACROSSFACE) {
              //assert(0); // This should be not possible.
              terminate_tet_core(this, 2);
            } 
            else {
              if ((dir == ACROSSVERT) || (dir == TOUCHEDGE) || 
                  (dir == TOUCHFACE)) {
                // This should be caused by a self-intersection.
                rejflag = 1; // Do not flip.
              }
            }
          }
        } // if (tmppts[0] != dummypoint)
      } // i
    } else if (fliptype == 2) {
      // A 3-to-2 flip: [e,d,a], [e,d,b], [e,d,c] => [a,b,c]
      if (pc != dummypoint) {
        // Check if the new face [a,b,c] intersect the edge in its interior.
        intflag = tri_edge_test(pa, pb, pc, fc->seg[0], fc->seg[1], NULL, 
                                1, types, poss);
        if (intflag == 2) {
          // They intersect at a single point.
          dir = (enum interresult) types[0];
          if (dir == ACROSSFACE) {
            // The interior of [a,b,c] intersect the segment.
            rejflag = 1; // Do not flip.
          }
        } else if (intflag == 4) {
          // [a,b,c] is coplanar with the edge. 
          dir = (enum interresult) types[0];
          if (dir == ACROSSEDGE) {
            // The boundary of [a,b,c] intersect the segment.            
            rejflag = 1; // Do not flip.
          }
        }
      } // if (pc != dummypoint)
    }
  } // if (fc->seg[0] != NULL)

  if ((fc->fac[0] != NULL) && !rejflag) {
    // A constraining face is given (e.g., for face recovery).
    if (fliptype == 1) {
      // A 2-to-3 flip.
      // Test if the new edge [e,d] intersects the face.
      intflag = tri_edge_test(fc->fac[0], fc->fac[1], fc->fac[2], pe, pd, 
                              NULL, 1, types, poss);
      if (intflag == 2) {
        // They intersect at a single point.
        dir = (enum interresult) types[0];
        if (dir == ACROSSFACE) {
          rejflag = 1;
        } else if (dir == ACROSSEDGE) {
          rejflag = 1;
        } 
      } else if (intflag == 4) {
        // The edge [e,d] is coplanar with the face.
        // There may be two intersections.
        for (i = 0; i < 2 && !rejflag; i++) {
          dir = (enum interresult) types[i];
          if (dir == ACROSSFACE) {
            rejflag = 1;
          } else if (dir == ACROSSEDGE) {
            rejflag = 1;
          }
        }
      }
    } // if (fliptype == 1)
  } // if (fc->fac[0] != NULL)

  if ((fc->remvert != NULL) && !rejflag) {
    // The vertex is going to be removed. Do not create a new edge which
    //   contains this vertex.
    if (fliptype == 1) {
      // A 2-to-3 flip.
      if ((pd == fc->remvert) || (pe == fc->remvert)) {
        rejflag = 1;
      }
    }
  }

  if (fc->remove_large_angle && !rejflag) {
    // Remove a large dihedral angle. Do not create a new small angle.
    badface bf; // used by get_tetqual(...)
    double cosmaxd = 0, diff;
    if (fliptype == 1) {
      // We assume that neither 'a' nor 'b' is dummypoint.
      // A 2-to-3 flip: [a,b,c] => [e,d,a], [e,d,b], [e,d,c].
      // The new tet [e,d,a,b] will be flipped later. Only two new tets:
      //   [e,d,b,c] and [e,d,c,a] need to be checked.
      if ((pc != dummypoint) && (pe != dummypoint) && (pd != dummypoint)) {
        double min_cosmaxd = 1.0, max_asp = 0; // record the worst quality.
        // Get the largest dihedral angle of [e,d,b,c].
        if (get_tetqual(pe, pd, pb, pc, &bf)) {
          cosmaxd = bf.cent[0];
          diff = cosmaxd - fc->cosdihed_in;
          if (fabs(diff/fc->cosdihed_in) < b->epsilon) diff = 0.0; // Rounding.
        } else {
          diff = 0.0; // no improve.
        }
        if (diff <= 0) { //if (cosmaxd <= fc->cosdihed_in) {
          rejflag = 1;
        } else {
          // Record the largest new angle.
          min_cosmaxd = (min_cosmaxd < cosmaxd ? min_cosmaxd : cosmaxd);
          max_asp = (max_asp > bf.key ? max_asp : bf.key);
          // Get the largest dihedral angle of [e,d,c,a].
          if (get_tetqual(pe, pd, pc, pa, &bf)) {
            cosmaxd = bf.cent[0];
            diff = cosmaxd - fc->cosdihed_in;
            if (fabs(diff/fc->cosdihed_in) < b->epsilon) diff = 0.0; // Rounding.
          } else {
            diff = 0.0; // no improve.
          }
          if (diff <= 0) { //if (cosmaxd <= fc->cosdihed_in) {
            rejflag = 1;
          } else {
            // Record the largest new angle.
            min_cosmaxd = (min_cosmaxd < cosmaxd ? min_cosmaxd : cosmaxd);
            max_asp = (max_asp > bf.key ? max_asp : bf.key);
            // save the worst quality.
            fc->cosdihed_out = (fc->cosdihed_out < min_cosmaxd ? fc->cosdihed_out : min_cosmaxd);
            fc->max_asp_out = (fc->max_asp_out > max_asp ? fc->max_asp_out : max_asp);
          }
        }
      } // if (pc != dummypoint && ...)
    } else if (fliptype == 2) {
      // A 3-to-2 flip: [e,d,a], [e,d,b], [e,d,c] => [a,b,c]
      // We assume that neither 'e' nor 'd' is dummypoint.
      if (level == 0) {
        // Both new tets [a,b,c,d] and [b,a,c,e] are new tets.
        if ((pa != dummypoint) && (pb != dummypoint) && (pc != dummypoint)) {
          double min_cosmaxd = 1.0, max_asp = 0; // record the worst quality.
          // Get the largest dihedral angle of [a,b,c,d].
          if (get_tetqual(pa, pb, pc, pd, &bf)) {
            cosmaxd = bf.cent[0];
            diff = cosmaxd - fc->cosdihed_in;
            if (fabs(diff/fc->cosdihed_in) < b->epsilon) diff = 0.0; // Rounding.
          } else {
            diff = 0.0; // no improve.
          }
          if (diff <= 0) { //if (cosmaxd <= fc->cosdihed_in) {
            rejflag = 1;
          } else {
            // Record the largest new angle.
            min_cosmaxd = (min_cosmaxd < cosmaxd ? min_cosmaxd : cosmaxd);
            max_asp = (max_asp > bf.key ? max_asp : bf.key);
            // Get the largest dihedral angle of [b,a,c,e].
            if (get_tetqual(pb, pa, pc, pe, &bf)) {
              cosmaxd = bf.cent[0];
              diff = cosmaxd - fc->cosdihed_in;
              if (fabs(diff/fc->cosdihed_in) < b->epsilon) diff = 0.0; // Rounding.
            } else {
              diff = 0.0; // no improve.
            }
            if (diff <= 0) { //if (cosmaxd <= fc->cosdihed_in) {
              rejflag = 1;
            } else {
              // Record the largest new angle.
              min_cosmaxd = (min_cosmaxd < cosmaxd ? min_cosmaxd : cosmaxd);
              max_asp = (max_asp > bf.key ? max_asp : bf.key);
              // save the worst quality.
              fc->cosdihed_out = (fc->cosdihed_out < min_cosmaxd ? fc->cosdihed_out : min_cosmaxd);
              fc->max_asp_out = (fc->max_asp_out > max_asp ? fc->max_asp_out : max_asp);
            }
          }
        }
      } else { // level > 0
        if (edgepivot == 1) {
          // The new tet [a,b,c,d] will be flipped. Only check [b,a,c,e].
          if ((pa != dummypoint) && (pb != dummypoint) && (pc != dummypoint)) {
            // Get the largest dihedral angle of [b,a,c,e].
            if (get_tetqual(pb, pa, pc, pe, &bf)) {
              cosmaxd = bf.cent[0];
              diff = cosmaxd - fc->cosdihed_in;
              if (fabs(diff/fc->cosdihed_in) < b->epsilon) diff = 0.0; // Rounding.
            } else {
              diff = 0.0; // no improve.
            }
            if (diff <= 0) { //if (cosmaxd <= fc->cosdihed_in) {
              rejflag = 1;
            } else {
              // Record the largest new angle.
              // save the worst quality.
              fc->cosdihed_out = (fc->cosdihed_out < cosmaxd ? fc->cosdihed_out : cosmaxd);
              fc->max_asp_out = (fc->max_asp_out > bf.key ? fc->max_asp_out : bf.key);
            }
          }
        } else {
          // The new tet [b,a,c,e] will be flipped. Only check [a,b,c,d].
          if ((pa != dummypoint) && (pb != dummypoint) && (pc != dummypoint)) {
            // Get the largest dihedral angle of [b,a,c,e].
            if (get_tetqual(pa, pb, pc, pd, &bf)) {
              cosmaxd = bf.cent[0];
              diff = cosmaxd - fc->cosdihed_in;
              if (fabs(diff/fc->cosdihed_in) < b->epsilon) diff = 0.0; // Rounding.
            } else {
              diff = 0.0; // no improve.
            }
            if (diff <= 0) { //if (cosmaxd <= fc->cosdihed_in) {
              rejflag = 1;
            } else {
              // Record the largest new angle.
              // save the worst quality.
              fc->cosdihed_out = (fc->cosdihed_out < cosmaxd ? fc->cosdihed_out : cosmaxd);
              fc->max_asp_out = (fc->max_asp_out > bf.key ? fc->max_asp_out : bf.key);
            }
          }
        } // edgepivot
      } // level
    }
  } // if (fc->remove_large_angle && !rejflag)

  return rejflag;
}

//============================================================================//
//                                                                            //
// removeedgebyflips()    Attempt to remove an edge by flips.                 //
//                                                                            //
// 'flipedge' is a non-convex or flat edge [a,b,#,#] to be removed.           //
//                                                                            //
// The return value is a positive integer, it indicates whether the edge is   //
// removed or not.  A value "2" means the edge is removed, otherwise, the     //
// edge is not removed and the value (must >= 3) is the current number of     //
// tets in the edge star.                                                     //
//                                                                            //
//============================================================================//

int TetMeshCore::removeedgebyflips(triface *flipedge, flipconstraints* fc)
{
  triface *abtets, spintet;
  int t1ver; 
  int n, nn, i;


  if (checksubsegflag) {
    // Do not flip a segment.
    if (issubseg(*flipedge)) {
      if (fc->collectencsegflag) {
        face checkseg, *paryseg;
        tsspivot1(*flipedge, checkseg);
        if (!sinfected(checkseg)) {
          // Queue this segment in list.
          sinfect(checkseg);                
          caveencseglist->newindex((void **) &paryseg);
          *paryseg = checkseg;
        }
      }
      return 0;
    }
  }

  // Count the number of tets at edge [a,b].
  int subface_count = 0; // count the # of subfaces at this edge.
  n = 0;
  spintet = *flipedge;
  while (1) {
    if (issubface(spintet)) subface_count++;
    n++;
    fnextself(spintet);
    if (spintet.tet == flipedge->tet) break;
  }
  if (n < 3) {
    // It is only possible when the mesh contains inverted tetrahedra.  
    terminate_tet_core(this, 2); // Report a bug
  }

  if (fc->noflip_in_surface) {
    if (subface_count > 0) {
      return 0;
    }
  }

  //if ((b->flipstarsize > 0) && (n > (b->flipstarsize+4))) {
  if ((b->flipstarsize > 0) && (n > b->flipstarsize)) {
    // The star size exceeds the limit.
    return 0; // Do not flip it.
  }

  // Allocate spaces.
  abtets = new triface[n];
  // Collect the tets at edge [a,b].
  spintet = *flipedge;
  for (i = 0; i < n; i++) {
    abtets[i] = spintet;
    setelemcounter(abtets[i], 1);
    fnextself(spintet);
  }


  // Try to flip the edge (level = 0, edgepivot = 0).
  nn = flipnm(abtets, n, 0, 0, fc);


  if (nn > 2) {
    // Edge is not flipped. Unmarktest the remaining tets in Star(ab).
    for (i = 0; i < nn; i++) {
      setelemcounter(abtets[i], 0);
    }
    // Restore the input edge (needed by Lawson's flip).
    *flipedge = abtets[0];
  }

  // Release the temporary allocated spaces.
  // NOTE: fc->unflip must be 0.
  int bakunflip = fc->unflip;
  fc->unflip = 0;
  flipnm_post(abtets, n, nn, 0, fc);
  fc->unflip = bakunflip;

  delete [] abtets;

  return nn; 
}

//============================================================================//
//                                                                            //
// removefacebyflips()    Remove a face by flips.                             //
//                                                                            //
// Return 1 if the face is removed. Otherwise, return 0.                      //
//                                                                            //
// ASSUMPTION: 'flipface' must not be a subface or a hull face.               //
//                                                                            //
//============================================================================//

int TetMeshCore::removefacebyflips(triface *flipface, flipconstraints* fc)
{
  triface fliptets[3], flipedge;
  point pa, pb, pc, pd, pe;
  double ori;
  int reducflag = 0;

  fliptets[0] = *flipface;
  fsym(*flipface, fliptets[1]);
  pa = org(fliptets[0]);
  pb = dest(fliptets[0]);
  pc = apex(fliptets[0]);
  pd = oppo(fliptets[0]);
  pe = oppo(fliptets[1]);

  ori = orient3d(pa, pb, pd, pe);
  if (ori > 0) {
    ori = orient3d(pb, pc, pd, pe);
    if (ori > 0) {
      ori = orient3d(pc, pa, pd, pe);
      if (ori > 0) {
        // Found a 2-to-3 flip.
        reducflag = 1;
      } else {
        eprev(*flipface, flipedge); // [c,a]
      }
    } else {
      enext(*flipface, flipedge); // [b,c]
    }
  } else {
    flipedge = *flipface; // [a,b]
  }

  if (reducflag) {
    triface checkface = fliptets[0];
    if (!valid_constrained_f23(checkface, pd, pe)) {
      return 0; //reducflag = 0;
    }
  }

  if (reducflag) {
    // A 2-to-3 flip is found.
    flip23(fliptets, 0, fc);
    return 1;
  } else {
    // Try to flip the selected edge of this face.
    if (removeedgebyflips(&flipedge, fc) == 2) {
      if (b->verbose > 3) {
        printf("      Face is removed by removing an edge.\n");
      }
      return 1;
    }
  }

  // Face is not removed.
  return 0;
}

//============================================================================//
//                                                                            //
// recoveredgebyflips()    Recover an edge in current tetrahedralization.     //
//                                                                            //
// If the edge is recovered, 'searchtet' returns a tet containing the edge.   //
//                                                                            //
// If the parameter 'fullsearch' is set, it tries to flip any face or edge    //
// that intersects the recovering edge.  Otherwise, only the face or edge     //
// which is visible by 'startpt' is tried.                                    //
//                                                                            //
// The parameter 'sedge' is used to report self-intersection. If it is not    //
// a NULL, it is EITHER a segment OR a subface that contains this edge.       //
//                                                                            //
// This routine assumes that the tetrahedralization is convex.                //
//                                                                            //
//============================================================================//

int TetMeshCore::recoveredgebyflips(point startpt, point endpt, face *sedge,
                                   triface* searchtet, int fullsearch, int& idir)
{
  flipconstraints fc;
  enum interresult dir;

  idir = (int) DISJOINT; // init.

  fc.seg[0] = startpt;
  fc.seg[1] = endpt;
  fc.checkflipeligibility = 1;

  // The mainloop of the edge reocvery.
  while (1) { // Loop I

    // Search the edge from 'startpt'.
    point2tetorg(startpt, *searchtet);
    dir = finddirection(searchtet, endpt);

    if (dir == ACROSSVERT) {
      if (dest(*searchtet) == endpt) {
        return 1; // Edge is recovered.
      } else {
        if (sedge != NULL) {
          // It is a segment or a subedge (an edge of a facet).
          // Check and report if there exists a self-intersection.
          insertvertexflags ivf;
          bool intersect_flag = false;
          point nearpt = dest(*searchtet);
          ivf.iloc = ONVERTEX;
          
          if (sedge->sh[5] == NULL) {
            // It is a segment.
            if (!issteinerpoint(nearpt)) {
              // It is an input point.
              if (!b->quiet && !b->nowarning) {
                int segidx = getfacetindex(*sedge);
                point p1 = segmentendpointslist[segidx*2];
                point p2 = segmentendpointslist[segidx*2+1];
                point tmppt = NULL;
                if (is_segment(p1, nearpt)) tmppt = p1;
                else if (is_segment(p2, nearpt)) tmppt = p2;
                if (tmppt != NULL) {
                  printf("Warning:  Two line segments are %s overlapping.\n",
                         ivf.iloc == NEARVERTEX ? "nearly" : "exactly");
                  printf("  1st: [%d,%d].\n", pointmark(p1), pointmark(p2));
                  printf("  2nd: [%d,%d].\n", pointmark(tmppt), pointmark(nearpt));
                } else {
                  printf("Warning:  A vertex lies %s on a line segment.\n",
                         ivf.iloc == NEARVERTEX ? "nearly" : "exactly");
                  printf("  1st: [%d,%d].\n", pointmark(p1), pointmark(p2));
                  printf("  2nd: [%d].\n", pointmark(nearpt));
                }
              }
              intersect_flag = true;
            } else {
              if (pointtype(nearpt) == FREESEGVERTEX) {
                // Check if two segments are (nearly) intersecting.
                int segidx = getfacetindex(*sedge);
                face parsentseg;
                sdecode(point2sh(nearpt), parsentseg);
                int segidx2 = getfacetindex(parsentseg);
                if (segidx2 != segidx) {
                  if (!b->quiet && !b->nowarning) { // -no -Q no -W
                    point p1 = segmentendpointslist[segidx*2];
                    point p2 = segmentendpointslist[segidx*2+1];
                    point p3 = segmentendpointslist[segidx2*2];
                    point p4 = segmentendpointslist[segidx2*2+1];
                    printf("Warning:  Two line segments are %s crossing.\n",
                           ivf.iloc == NEARVERTEX ? "nearly" : "exactly");
                    printf("  1st: [%d,%d].\n", pointmark(p1), pointmark(p2));
                    printf("  2nd: [%d,%d].\n", pointmark(p3), pointmark(p4));
                  }
                  intersect_flag = true;
                } else {
                  //if (ivf.iloc == ONVERTEX) {
                  terminate_tet_core(this, 2); // This should not be possible.
                  //}
                }
              } else if (pointtype(nearpt) == FREEFACETVERTEX) {
                // This case is very unlikely.
                terminate_tet_core(this, 2); // to debug...
                if (!b->quiet && !b->nowarning) { // -no -Q no -W
                  //face parsentsh;
                  //sdecode(point2sh(nearpt), parsentsh);
                  printf("Warning:  A segment and a facet intersect.\n");
                }
                intersect_flag = true;
              } else {
                // other cases...
                terminate_tet_core(this, 2); // to be checked.
              }
            }
          } else {
            // It is an edge of a facet.
            if (!issteinerpoint(nearpt)) {
              if (!b->quiet && !b->nowarning) { // no "-Q -W"
                point p1 =  sorg(*sedge);
                point p2 = sdest(*sedge);
                point p3 = sapex(*sedge);
                printf("Warning:  A vertex lies on a facet.\n");
                printf("  vertex: [%d]\n", pointmark(nearpt));
                printf("  facet triangle: [%d,%d,%d], tag(%d).\n",
                       pointmark(p1), pointmark(p2), pointmark(p3),
                       shellmark(*sedge));
              }
              intersect_flag = true;
            } else {
              // A Steiner point.
              if (pointtype(nearpt) == FREESEGVERTEX) {
                // A facet and a segment is intersecting.
                if (!b->quiet && !b->nowarning) {
                  printf("Warning:  A facet and a segment intersect.\n");
                  printf("  ...\n");
                }
                intersect_flag = true;
              } else if (pointtype(nearpt) == FREEFACETVERTEX) {
                // Check if two facets are intersecting.
                if (!b->quiet && !b->nowarning) {
                  printf("Warning:  Two facets intersect.\n");
                  printf("  ...\n");
                }
                intersect_flag = true;
              } else {
                // A FREEVOLVERTEX.
                // This is not a real self-intersection.
                terminate_tet_core(this, 2); // check this case.
              }
            }
          }
          
          if (intersect_flag) {
            idir = (int) SELF_INTERSECT;
          }
        } // if (sedge != NULL)
        return 0;
      }
    } // if (dir == ACROSSVERT)

    // The edge is missing. 

    // Try to remove the first intersecting face/edge.
    enextesymself(*searchtet); // Go to the opposite face.

    if (dir == ACROSSFACE) {
      if (checksubfaceflag) {
        if (issubface(*searchtet)) {
          if (sedge) {
            // A self-intersection is detected.
            if (!b->quiet && !b->nowarning) {
              bool is_seg = (sedge->sh[5] == NULL);
              if (is_seg) {
                face fac; tspivot(*searchtet, fac);
                int segidx = getfacetindex(*sedge);
                point p1 = segmentendpointslist[segidx*2];
                point p2 = segmentendpointslist[segidx*2+1];
                printf("Warning:  A segment and a facet exactly intersect.\n");
                printf("  seg  : [%d,%d].\n", pointmark(p1), pointmark(p2));
                printf("  facet triangle: [%d,%d,%d] tag(%d).\n",
                       pointmark(sorg(fac)), pointmark(sdest(fac)),
                       pointmark(sapex(fac)), shellmark(fac));
              } else {
                // It is a subedge of a facet.
                point *ppt = (point *) &(sedge->sh[3]);
                printf("Warning:  Two facets exactly intersect.\n");
                printf("  1st facet triangle: [%d,%d,%d] tag(%d).\n",
                       pointmark(ppt[0]), pointmark(ppt[1]),
                       pointmark(ppt[2]), shellmark(*sedge));
                face fac; tspivot(*searchtet, fac);
                ppt = (point *) &(fac.sh[3]);
                printf("  2nd facet triangle: [%d,%d,%d] tag(%d).\n",
                       pointmark(ppt[0]), pointmark(ppt[1]),
                       pointmark(ppt[2]), shellmark(fac));
              }
            }
            idir = (int) SELF_INTERSECT;
          }
          return 0;
        } // if (issubface(*searchtet))
      }
      // Try to flip a crossing face.
      if (removefacebyflips(searchtet, &fc)) {
        continue;
      }
    } else if (dir == ACROSSEDGE) {
      if (checksubsegflag) {
        if (issubseg(*searchtet)) {
          if (sedge) {
            // A self-intersection is detected.
            if (!b->quiet && !b->nowarning) { // no -Q, -W
              bool is_seg = (sedge->sh[5] == NULL);
              if (is_seg) {
                face seg; tsspivot1(*searchtet, seg);
                int segidx  = getfacetindex(*sedge);
                int segidx2 = getfacetindex(seg);
                if (segidx != segidx2) {
                  point p1 = segmentendpointslist[segidx*2];
                  point p2 = segmentendpointslist[segidx*2+1];
                  point p3 = segmentendpointslist[segidx2*2];
                  point p4 = segmentendpointslist[segidx2*2+1];
                  printf("Warning:  Two segments exactly intersect.\n");
                  printf("  1st seg  [%d,%d] tag(%d).\n",
                         pointmark(p1), pointmark(p2), shellmark(*sedge));
                  printf("  2nd seg: [%d,%d] tag(%d).\n",
                         pointmark(p3), pointmark(p4), shellmark(seg));
                } else {
                  terminate_tet_core(this, 2);
                }
              } else {
                // It is a subedge of a facet.
                point *ppt = (point *) &(sedge->sh[3]);
                printf("Warning:  A facet and a segment exactly intersect.\n");
                printf("  facet triangle: [%d,%d,%d] tag(%d).\n",
                       pointmark(ppt[0]), pointmark(ppt[1]),
                       pointmark(ppt[2]), shellmark(*sedge));
                face seg; tsspivot1(*searchtet, seg);
                ppt = (point *) &(seg.sh[3]);
                printf("  seg: [%d,%d] tag(%d).\n",
                       pointmark(ppt[0]), pointmark(ppt[1]), shellmark(seg));
              }
            }
            idir = (int) SELF_INTERSECT;
          }
          return 0;
        }
      }
      // Try to flip an intersecting edge.
      if (removeedgebyflips(searchtet, &fc) == 2) {
        continue;
      }
    } else {
      terminate_tet_core(this, 2); // report a bug
    }

    // The edge is missing.

    if (fullsearch) {
      // Try to flip one of the faces/edges which intersects the edge.
      triface neightet, spintet;
      point pa, pb, pc, pd;
      badface bakface;
      enum interresult dir1;
      int types[2], poss[4], pos = 0;
      int success = 0;
      int t1ver; 
      int i, j;

      // Loop through the sequence of intersecting faces/edges from
      //   'startpt' to 'endpt'.
      point2tetorg(startpt, *searchtet);
      dir = finddirection(searchtet, endpt);

      // Go to the face/edge intersecting the searching edge.
      enextesymself(*searchtet); // Go to the opposite face.
      // This face/edge has been tried in previous step.

      while (1) { // Loop I-I

        // Find the next intersecting face/edge.
        fsymself(*searchtet);
        if (dir == ACROSSFACE) {
          neightet = *searchtet;
          j = (neightet.ver & 3); // j is the current face number.
          for (i = j + 1; i < j + 4; i++) {
            neightet.ver = (i % 4);
            pa = org(neightet);
            pb = dest(neightet);
            pc = apex(neightet);
            pd = oppo(neightet); // The above point.
            if (tri_edge_test(pa,pb,pc,startpt,endpt, pd, 1, types, poss)) {
              dir = (enum interresult) types[0];
              pos = poss[0];
              break;
            } else {
              dir = DISJOINT;
              pos = 0;
            }
          } // i
          // There must be an intersection face/edge.
          if (dir == DISJOINT) {
            terminate_tet_core(this, 2);
          }
        } else if (dir == ACROSSEDGE) {
          while (1) { // Loop I-I-I
            // Check the two opposite faces (of the edge) in 'searchtet'.  
            for (i = 0; i < 2; i++) {
              if (i == 0) {
                enextesym(*searchtet, neightet);
              } else {
                eprevesym(*searchtet, neightet);
              }
              pa = org(neightet);
              pb = dest(neightet);
              pc = apex(neightet);
              pd = oppo(neightet); // The above point.
              if (tri_edge_test(pa,pb,pc,startpt,endpt,pd,1, types, poss)) {
                dir = (enum interresult) types[0];
                pos = poss[0];
                break; // for loop
              } else {
                dir = DISJOINT;
                pos = 0;
              }
            } // i
            if (dir != DISJOINT) {
              // Find an intersection face/edge.
              break;  // Loop I-I-I
            }
            // No intersection. Rotate to the next tet at the edge.
            fnextself(*searchtet);
          } // while (1) // Loop I-I-I
        } else {
          terminate_tet_core(this, 2); // Report a bug
        }

        // Adjust to the intersecting edge/vertex.
        for (i = 0; i < pos; i++) {
          enextself(neightet);
        }

        if (dir == SHAREVERT) {
          // Check if we have reached the 'endpt'.
          pd = org(neightet);
          if (pd == endpt) {
            // Failed to recover the edge.
            break; // Loop I-I
          } else {
            return 0;
          }
        }

        // The next to be flipped face/edge.
        *searchtet = neightet;

        // Bakup this face (tetrahedron).
        bakface.forg = org(*searchtet);
        bakface.fdest = dest(*searchtet);
        bakface.fapex = apex(*searchtet);
        bakface.foppo = oppo(*searchtet);

        // Try to flip this intersecting face/edge.
        if (dir == ACROSSFACE) {
          if (checksubfaceflag) {
            if (issubface(*searchtet)) {
              return 0;
            }
          }
          if (removefacebyflips(searchtet, &fc)) {
            success = 1;
            break; // Loop I-I 
          }
        } else if (dir == ACROSSEDGE) {
          if (checksubsegflag) {
            if (issubseg(*searchtet)) {
              return 0;
            }
          }
          if (removeedgebyflips(searchtet, &fc) == 2) {
            success = 1;
            break; // Loop I-I
          }
        } else if (dir == ACROSSVERT) {
          return 0;
        } else {
          terminate_tet_core(this, 2); 
        }

        // The face/edge is not flipped.
        if ((searchtet->tet == NULL) ||
            (org(*searchtet) != bakface.forg) ||
            (dest(*searchtet) != bakface.fdest) ||
            (apex(*searchtet) != bakface.fapex) ||
            (oppo(*searchtet) != bakface.foppo)) {
          // 'searchtet' was flipped. We must restore it.
          point2tetorg(bakface.forg, *searchtet);
          dir1 = finddirection(searchtet, bakface.fdest);
          if (dir1 == ACROSSVERT) {
            if (dest(*searchtet) == bakface.fdest) {
              spintet = *searchtet;
              while (1) {
                if (apex(spintet) == bakface.fapex) {
                  // Found the face.
                  *searchtet = spintet;
                  break;
                }
                fnextself(spintet);
                if (spintet.tet == searchtet->tet) {
                  searchtet->tet = NULL;
                  break; // Not find.
                }
	          } // while (1)
              if (searchtet->tet != NULL) {
                if (oppo(*searchtet) != bakface.foppo) {
                  fsymself(*searchtet);
                  if (oppo(*searchtet) != bakface.foppo) {
                    // The original (intersecting) tet has been flipped.
                    searchtet->tet = NULL;
                    break; // Not find.
                  }
                }
              }
            } else {
              searchtet->tet = NULL; // Not find.
            }
          } else {
            searchtet->tet = NULL; // Not find.
          }
          if (searchtet->tet == NULL) {
            success = 0; // This face/edge has been destroyed.
            break; // Loop I-I 
          }
        }
      } // while (1) // Loop I-I

      if (success) {
        // One of intersecting faces/edges is flipped.
        continue;
      }

    } // if (fullsearch)

    // The edge is missing.
    break; // Loop I

  } // while (1) // Loop I

  return 0;
}

//============================================================================//
//                                                                            //
// add_steinerpt_in_schoenhardtpoly()    Insert a Steiner point in a Schoen-  //
//                                       hardt polyhedron.                    //
//                                                                            //
// 'abtets' is an array of n tets which all share at the edge [a,b]. Let the  //
// tets are [a,b,p0,p1], [a,b,p1,p2], ..., [a,b,p_(n-2),p_(n-1)].  Moreover,  //
// the edge [p0,p_(n-1)] intersects all of the tets in 'abtets'.  A special   //
// case is that the edge [p0,p_(n-1)] is coplanar with the edge [a,b].        //
// Such set of tets arises when we want to recover an edge from 'p0' to 'p_   //
// (n-1)', and the number of tets at [a,b] can not be reduced by any flip.    //
//                                                                            //
//============================================================================//

int TetMeshCore::add_steinerpt_in_schoenhardtpoly(triface *abtets, int n,
  int splitsliverflag, int chkencflag)
{
  triface worktet, *parytet;
  triface faketet1, faketet2;
  point pc, pd, steinerpt;
  insertvertexflags ivf;
  optparameters opm;
  double vcd[3], sampt[3], smtpt[3];
  double maxminvol = 0.0, minvol = 0.0, ori;
  int success, maxidx = 0;
  int it, i;


  if (splitsliverflag) {
    // randomly pick a tet.
    int idx = rand() % n;

    // Calulcate the barycenter of this tet.
    point pa = org(abtets[idx]);
    point pb = dest(abtets[idx]);
    pc = apex(abtets[idx]);
    pd = oppo(abtets[idx]);

    makepoint(&steinerpt, FREEVOLVERTEX);
    for (i = 0; i < 3; i++) {
      steinerpt[i] = (pa[i] + pb[i] + pc[i] + pd[i]) / 4.;
    }


    worktet = abtets[idx];
    ivf.iloc = (int) OUTSIDE; // need point location.
    ivf.bowywat = 1;
    //ivf.lawson = 0;
    ivf.lawson = 2; // Do flips to recover Delaunayness.
    ivf.rejflag = 0;
    ivf.chkencflag = chkencflag;
    ivf.sloc = 0;
    ivf.sbowywat = 0;
    ivf.splitbdflag = 0;
    ivf.validflag = 1;
    ivf.respectbdflag = 1;
    ivf.assignmeshsize = b->metric;

    if (insertpoint(steinerpt, &worktet, NULL, NULL, &ivf)) {
      // The vertex has been inserted.
      if (flipstack != NULL) {
        recoverdelaunay();
      }
      st_volref_count++;
      if (steinerleft > 0) steinerleft--;
      return 1;
    } else {
      // Not inserted.
      pointdealloc(steinerpt);
      return 0;
    }
  } // if (splitsliverflag)

  pc = apex(abtets[0]);   // pc = p0
  pd = oppo(abtets[n-1]); // pd = p_(n-1)

  // Find an optimial point in edge [c,d]. It is visible by all outer faces
  //   of 'abtets', and it maxmizes the min volume.

  // initialize the list of 2n boundary faces.
  for (i = 0; i < n; i++) {    
    edestoppo(abtets[i], worktet); // [p_i,p_i+1,a]
    cavetetlist->newindex((void **) &parytet);
    *parytet = worktet;
    eorgoppo(abtets[i], worktet);  // [p_i+1,p_i,b]
    cavetetlist->newindex((void **) &parytet);
    *parytet = worktet;
  }

  int N = 100;
  double stepi = 0.01;

  // Search the point along the edge [c,d].
  for (i = 0; i < 3; i++) vcd[i] = pd[i] - pc[i];

  // Sample N points in edge [c,d].
  for (it = 1; it < N; it++) {
    for (i = 0; i < 3; i++) {
      sampt[i] = pc[i] + (stepi * (double) it) * vcd[i];
    }
    for (i = 0; i < cavetetlist->objects; i++) {
      parytet = (triface *) fastlookup(cavetetlist, i);
      ori = orient3d(dest(*parytet), org(*parytet), apex(*parytet), sampt);
      if (i == 0) {
        minvol = ori;
      } else {
        if (minvol > ori) minvol = ori;
      }
    } // i
    if (it == 1) {
      maxminvol = minvol;
      maxidx = it;
    } else {
      if (maxminvol < minvol) {
        maxminvol = minvol;
        maxidx = it;
      } 
    }
  } // it

  if (maxminvol <= 0) {
    cavetetlist->restart();
    return 0;
  }

  for (i = 0; i < 3; i++) {
    smtpt[i] = pc[i] + (stepi * (double) maxidx) * vcd[i];
  }

  // Create two faked tets to hold the two non-existing boundary faces:
  //   [d,c,a] and [c,d,b].
  maketetrahedron(&faketet1);
  setvertices(faketet1, pd, pc, org(abtets[0]), dummypoint);
  cavetetlist->newindex((void **) &parytet);
  *parytet = faketet1;
  maketetrahedron(&faketet2);
  setvertices(faketet2, pc, pd, dest(abtets[0]), dummypoint);
  cavetetlist->newindex((void **) &parytet);
  *parytet = faketet2;

  // Point smooth options.
  opm.max_min_volume = 1;
  opm.numofsearchdirs = 20;
  opm.searchstep = 0.001;  
  opm.maxiter = 100; // Limit the maximum iterations.
  opm.initval = 0.0; // Initial volume is zero.

  // Try to relocate the point into the inside of the polyhedron.
  success = smoothpoint(smtpt, cavetetlist, 1, &opm);

  if (success) {
    while (opm.smthiter == 100) {
      // It was relocated and the prescribed maximum iteration reached. 
      // Try to increase the search stepsize.
      opm.searchstep *= 10.0;
      //opm.maxiter = 100; // Limit the maximum iterations.
      opm.initval = opm.imprval;
      opm.smthiter = 0; // Init.
      smoothpoint(smtpt, cavetetlist, 1, &opm);  
    }
  } // if (success)

  // Delete the two faked tets.
  tetrahedrondealloc(faketet1.tet);
  tetrahedrondealloc(faketet2.tet);

  cavetetlist->restart();

  if (success) {
    // Insert this Steiner point.

    // Insert the Steiner point.
    makepoint(&steinerpt, FREEVOLVERTEX);
    for (i = 0; i < 3; i++) steinerpt[i] = smtpt[i];

    // Insert the created Steiner point.
    for (i = 0; i < n; i++) {
      infect(abtets[i]);
      caveoldtetlist->newindex((void **) &parytet);
      *parytet = abtets[i];
    }
    worktet = abtets[0]; // No need point location.
    ivf.iloc = (int) INSTAR;
    ivf.chkencflag = chkencflag;
    ivf.assignmeshsize = b->metric;
    if (ivf.assignmeshsize) {
      // Search the tet containing 'steinerpt' for size interpolation.
      locate(steinerpt, &(abtets[0]));
      worktet = abtets[0];
    }

    // Insert the new point into the tetrahedralization T.
    if (insertpoint(steinerpt, &worktet, NULL, NULL, &ivf)) {
      // The vertex has been inserted.
      st_volref_count++;
      if (steinerleft > 0) steinerleft--;
      return 1;
    } else {
      // Not inserted.
      pointdealloc(steinerpt);
      return 0;
    }
  }

  //if (!success) {
    return 0;
  //}
}

//============================================================================//
//                                                                            //
// add_steinerpt_in_segment()    Add a Steiner point inside a segment.        //
//                                                                            //
//============================================================================//

int TetMeshCore::add_steinerpt_in_segment(face* misseg, int searchlevel, int& idir)
{
  triface searchtet;
  face *paryseg, candseg;
  point startpt, endpt, pc, pd;
  flipconstraints fc;
  enum interresult dir;
  double P[3], Q[3], tp, tq;
  double len, smlen = 0, split = 0, split_q = 0;
  int success;
  int i;

  startpt = sorg(*misseg);
  endpt = sdest(*misseg);

  idir = DISJOINT; // init.
  
  // sort the vertices
  //if (pointmark(startpt) > pointmark(endpt)) {
  //  endpt = sorg(*misseg);
  //  startpt = sdest(*misseg);
  //}
  

  fc.seg[0] = startpt;
  fc.seg[1] = endpt;
  fc.checkflipeligibility = 1;
  fc.collectencsegflag = 1;

  point2tetorg(startpt, searchtet);
  dir = finddirection(&searchtet, endpt);
  if (dir == ACROSSVERT) {
    return 0;
  }

  // Try to flip the first intersecting face/edge.
  enextesymself(searchtet); // Go to the opposite face.

  int bak_fliplinklevel = b->fliplinklevel;
  b->fliplinklevel = searchlevel;

  if (dir == ACROSSFACE) {
    // A face is intersected with the segment. Try to flip it.
    success = removefacebyflips(&searchtet, &fc);
  } else if (dir == ACROSSEDGE) {
    // An edge is intersected with the segment. Try to flip it.
    success = removeedgebyflips(&searchtet, &fc);
  }

  split = 0;
  for (i = 0; i < caveencseglist->objects; i++) {
    paryseg = (face *) fastlookup(caveencseglist, i);
    suninfect(*paryseg);
    // Calculate the shortest edge between the two lines.
    pc = sorg(*paryseg);
    pd = sdest(*paryseg);
    
    // sort the vertices
    //if (pointmark(pc) > pointmark(pd)) {
    //  pd = sorg(*paryseg);
    //  pc = sdest(*paryseg);
    //}
    
    tp = tq = 0;
    if (linelineint(startpt, endpt, pc, pd, P, Q, &tp, &tq)) {
      // Does the shortest edge lie between the two segments? 
      // Round tp and tq.
      if ((tp > 0) && (tq < 1)) {
        if (tp < 0.5) {
          if (tp < (b->epsilon * 1e+3)) tp = 0.0;
        } else {
          if ((1.0 - tp) < (b->epsilon * 1e+3)) tp = 1.0;
        }
      }
      if ((tp <= 0) || (tp >= 1)) continue; 
      if ((tq > 0) && (tq < 1)) {
        if (tq < 0.5) {
          if (tq < (b->epsilon * 1e+3)) tq = 0.0;
        } else {
          if ((1.0 - tq) < (b->epsilon * 1e+3)) tq = 1.0;
        }
      }
      if ((tq <= 0) || (tq >= 1)) continue;
      // It is a valid shortest edge. Calculate its length.
      len = distance(P, Q);
      if (split == 0) {
        smlen = len;
        split = tp;
        split_q = tq;
        candseg = *paryseg;
      } else {
        if (len < smlen) {
          smlen = len;
          split = tp;
          split_q = tq;
          candseg = *paryseg;
        }
      }
    }
  }

  caveencseglist->restart();
  b->fliplinklevel = bak_fliplinklevel;

  if (split == 0) {
    // Found no crossing segment. 
    return 0;
  }

  face splitsh;
  face splitseg;
  point steinerpt, *parypt;
  insertvertexflags ivf;

  if (b->addsteiner_algo == 1) {
    // Split the segment at the closest point to a near segment.
    makepoint(&steinerpt, FREESEGVERTEX);
    for (i = 0; i < 3; i++) {
      steinerpt[i] = startpt[i] + split * (endpt[i] - startpt[i]);
    }
  } else { // b->addsteiner_algo == 2
    for (i = 0; i < 3; i++) {
      P[i] = startpt[i] + split * (endpt[i] - startpt[i]);
    }
    pc = sorg(candseg);
    pd = sdest(candseg);
    for (i = 0; i < 3; i++) {
      Q[i] = pc[i] + split_q * (pd[i] - pc[i]);
    }
    makepoint(&steinerpt, FREEVOLVERTEX);
    for (i = 0; i < 3; i++) {
      steinerpt[i] = 0.5 * (P[i] + Q[i]);
    }
  }

  // Check if the two segments are nearly crossing each other.
  pc = sorg(candseg);
  pd = sdest(candseg);
  if (is_collinear_at(steinerpt, pc, pd)) { // -p///#, default 179.9 degree
    if (!b->quiet && !b->nowarning) { // no -Q, -W
      int segidx = getfacetindex(*misseg);
      point p1 = segmentendpointslist[segidx*2];
      point p2 = segmentendpointslist[segidx*2+1];
      int segidx2 = getfacetindex(candseg);
      point p3 = segmentendpointslist[segidx2*2];
      point p4 = segmentendpointslist[segidx2*2+1];
      printf("Warning:  Two line segments are almost crossing.\n");
      printf("  1st: [%d,%d].\n", pointmark(p1), pointmark(p2));
      printf("  2nd: [%d,%d].\n", pointmark(p3), pointmark(p4));
    }

    // calculate a new angle tolerance.
    double collinear_ang = interiorangle(steinerpt, pc, pd, NULL) / PI * 180.;
    double ang_diff = collinear_ang - b->collinear_ang_tol;
    double new_ang_tol = collinear_ang + ang_diff / 180.;
 
    if (new_ang_tol < 180.0) { // no -Q, -W
      // Reduce the angle tolerance to detect collinear event.
      if (!b->quiet && !b->nowarning) {
        printf("  Reducing collinear tolerance from %g to %g degree.\n",
               b->collinear_ang_tol, new_ang_tol);
      }
      b->collinear_ang_tol = new_ang_tol;
      cos_collinear_ang_tol = cos(b->collinear_ang_tol / 180.0 * PI);
    } else {
      // Report a self-intersection event due to epsilon.
      if (!b->quiet && !b->nowarning) { // no -Q, -W
        printf("  Cannot reduce the current collinear tolerance (=%g degree).\n",
               b->collinear_ang_tol);
      }
      idir = SELF_INTERSECT;
      pointdealloc(steinerpt);
      return 0;
    }
  }

  // We need to locate the point. Start searching from 'searchtet'.
  if (split < 0.5) {
    point2tetorg(startpt, searchtet);
  } else {
    point2tetorg(endpt, searchtet);
  }
  if (b->addsteiner_algo == 1) {
    splitseg = *misseg;
    spivot(*misseg, splitsh);
    // for create_a_shorter_edge().
    setpoint2sh(steinerpt, sencode(*misseg));
  } else {
    splitsh.sh = NULL;
    splitseg.sh = NULL;
  }
  ivf.iloc = (int) OUTSIDE;
  ivf.bowywat = 1;
  //ivf.lawson = 0;
  ivf.lawson = 2; // Do flips to recover Delaunayness.
  ivf.rejflag = 0;
  ivf.chkencflag = 0;
  ivf.sloc = (int) ONEDGE;
  ivf.sbowywat = 1; // split surface mesh separately, new subsegments are
                    //   pushed into "subsegstack".
  ivf.splitbdflag = 0;
  ivf.validflag = 1;
  ivf.respectbdflag = 1;
  ivf.assignmeshsize = b->metric; 

  if (insertpoint(steinerpt, &searchtet, &splitsh, &splitseg, &ivf)) {
    if (flipstack != NULL) {
      recoverdelaunay();
    }
  } else {
    pointdealloc(steinerpt);
    return 0;
  }

  if (b->addsteiner_algo == 1) {
    // Save this Steiner point (for removal).
    //   Re-use the array 'subvertstack'.
    subvertstack->newindex((void **) &parypt);
    *parypt = steinerpt;
    st_segref_count++;
  } else { // b->addsteiner_algo == 2
    // Queue the segment for recovery.
    subsegstack->newindex((void **) &paryseg);
    *paryseg = *misseg; 
    st_volref_count++;
  }
  if (steinerleft > 0) steinerleft--;

  return 1;
}

//============================================================================//
//                                                                            //
// addsteiner4recoversegment()    Add a Steiner point for recovering a seg.   //
//                                                                            //
// Tries to add a Steiner point in the volume (near this segment) which will  //
// help to recover this segment.  This segment itself is not split.           //
//                                                                            //
// 'splitsliverflag' is a parameter used in the subroutine add_steiner_in_    //
// schonhardpoly().                                                           //
//                                                                            //
//============================================================================//

int TetMeshCore::add_steinerpt_to_recover_edge(point startpt, point endpt,
  face* misseg, int splitsegflag, int splitsliverflag, int& idir)
{
  triface *abtets, searchtet, spintet;
  face splitsh;
  face *paryseg;
  point pa, pb, pd, steinerpt, *parypt;
  enum interresult dir;
  insertvertexflags ivf;
  int types[2], poss[4];
  int n, endi, success;
  int t1ver;
  int i;

  idir = (int) DISJOINT;

  if (misseg != NULL) {
    startpt = sorg(*misseg);
    if (pointtype(startpt) == FREESEGVERTEX) {
      sesymself(*misseg);
      startpt = sorg(*misseg);
    }
    endpt = sdest(*misseg);
  }

  
  point2tetorg(startpt, searchtet);
  dir = finddirection(&searchtet, endpt);


  if (dir == ACROSSVERT) {
    if (dest(searchtet) == endpt) {
      // This edge exists.
      if ((misseg != NULL) && (subsegstack != NULL)) {
        // Add the missing segment back to the recovering list.
        subsegstack->newindex((void **) &paryseg);
        *paryseg = *misseg;
      }
      return 1;
    } else {
      // This edge crosses a vertex (not endpt).
      bool intersect_flag = false; // return
      if (misseg != NULL) {
        // Check whether there exists a self-intersection.
        point nearpt = dest(searchtet);
        ivf.iloc = ONVERTEX;
        // report_seg_vertex_intersect(misseg, dest(searchtet), ONVERTEX);
        int segidx = getfacetindex(*misseg);
        point p1 = segmentendpointslist[segidx*2];
        point p2 = segmentendpointslist[segidx*2+1];

        if (!issteinerpoint(nearpt)) {
          // It is an input point.
          if (!b->quiet && !b->nowarning) {
            point tmppt = NULL;
            if (is_segment(p1, nearpt)) tmppt = p1;
            else if (is_segment(p2, nearpt)) tmppt = p2;
            if (tmppt != NULL) {
              printf("Warning:  Two line segments are %s overlapping.\n",
                     ivf.iloc == NEARVERTEX ? "nearly" : "exactly");
              printf("  1st: [%d,%d].\n", pointmark(p1), pointmark(p2));
              printf("  2nd: [%d,%d].\n", pointmark(tmppt), pointmark(nearpt));
            } else {
              printf("Warning:  A vertex lies %s on a line segment.\n",
                     ivf.iloc == NEARVERTEX ? "nearly" : "exactly");
              printf("  segment: [%d,%d].\n", pointmark(p1), pointmark(p2));
              printf("  vertex : [%d].\n", pointmark(nearpt));
            }
          }
          intersect_flag = true;
        } else {
          if (pointtype(nearpt) == FREESEGVERTEX) {
            // Check if two segments are exactly intersecting.
            face parsentseg;
            sdecode(point2sh(nearpt), parsentseg);
            int segidx2 = getfacetindex(parsentseg);
            if (segidx2 != segidx) {
              if (!b->quiet && !b->nowarning) {
                point p3 = segmentendpointslist[segidx2*2];
                point p4 = segmentendpointslist[segidx2*2+1];
                printf("Warning:  Two line segments are %s crossing.\n",
                       ivf.iloc == NEARVERTEX ? "nearly" : "exactly");
                printf("  1st: [%d,%d].\n", pointmark(p1), pointmark(p2));
                printf("  2nd: [%d,%d].\n", pointmark(p3), pointmark(p4));
              }
              intersect_flag = true;
            } else {
              if (ivf.iloc == ONVERTEX) {
                terminate_tet_core(this, 2); // This should not be possible.
              }
            }
          } else {
            // other cases...
            terminate_tet_core(this, 2);
          }
        }
      } // if (misseg != NULL)
      if (intersect_flag) {
        idir = (int) SELF_INTERSECT;
      }
      return 0;
    }
  } // if (dir == ACROSSVERT) {

  enextself(searchtet);

  if (dir == ACROSSFACE) {
    // The segment is crossing at least 3 faces. Find the common edge of
    //   the first 3 crossing faces.
    esymself(searchtet);
    fsym(searchtet, spintet);
    pd = oppo(spintet);

    if (pd == endpt) {
      if (misseg != NULL) {
        // Calclate the smallest angle between (a,b,c) and (startpt, endpt).
        triface tmptet;
        double ang, collinear_ang = 0.;
        for (int k = 0; k < 3; k++) {
          ang = interiorangle(org(searchtet), startpt, endpt, NULL); // in [0, PI]
          if (ang > collinear_ang) {
            collinear_ang = ang;
            tmptet = searchtet; // org(tmptet)
          }
          enextself(searchtet);
        }
        collinear_ang = collinear_ang / PI * 180.; // in degree

        if (collinear_ang > b->collinear_ang_tol) { // -p///#, default 179.9 degree
          // Report a self-intersection event due to epsilon.
          if (!b->quiet && !b->nowarning) { // no -Q, -W
            point nearpt = org(tmptet);
            ivf.iloc = NEARVERTEX;
            // report_seg_vertex_intersect(misseg, dest(searchtet), ONVERTEX);
            int segidx = getfacetindex(*misseg);
            point p1 = segmentendpointslist[segidx*2];
            point p2 = segmentendpointslist[segidx*2+1];
              
            if (!issteinerpoint(nearpt)) {
              point tmppt = NULL;
              if (is_segment(p1, nearpt)) tmppt = p1;
              else if (is_segment(p2, nearpt)) tmppt = p2;
              if (tmppt != NULL) {
                printf("Warning:  Two line segments are %s overlapping.\n",
                       ivf.iloc == NEARVERTEX ? "nearly" : "exactly");
                printf("  1st: [%d,%d].\n", pointmark(p1), pointmark(p2));
                printf("  2nd: [%d,%d].\n", pointmark(tmppt), pointmark(nearpt));
              } else {
                printf("Warning:  A vertex lies %s on a line segment.\n",
                       ivf.iloc == NEARVERTEX ? "nearly" : "exactly");
                printf("  segment: [%d,%d].\n", pointmark(p1), pointmark(p2));
                printf("  vertex : [%d].\n", pointmark(nearpt));
              }
            } else {
              if (pointtype(nearpt) == FREESEGVERTEX) {
                // Check if two segments are nearly intersecting.
                face parsentseg;
                sdecode(point2sh(nearpt), parsentseg);
                int segidx2 = getfacetindex(parsentseg);
                if (segidx2 != segidx) {
                  //if (!b->quiet && !b->nowarning) {
                  point p3 = segmentendpointslist[segidx2*2];
                  point p4 = segmentendpointslist[segidx2*2+1];
                  printf("Warning:  Two line segments are %s crossing.\n",
                         ivf.iloc == NEARVERTEX ? "nearly" : "exactly");
                  printf("  1st: [%d,%d].\n", pointmark(p1), pointmark(p2));
                  printf("  2nd: [%d,%d].\n", pointmark(p3), pointmark(p4));
                  //}
                  //intersect_flag = true;
                } else {
                  //if (ivf.iloc == ONVERTEX) {
                  terminate_tet_core(this, 2); // This should not be possible.
                  //}
                }
              } else {
                // Other case to report.
                //  assert(0); // to do...
                terminate_tet_core(this, 2);
              }
            }
          }
          
          // calculate a new angle tolerance.
          double ang_diff = collinear_ang - b->collinear_ang_tol;
          double new_ang_tol = collinear_ang + ang_diff / 180.;
          
          if (new_ang_tol < 180.) {
            // Reduce the angle tolerance to detect collinear event.
            if (!b->quiet && !b->nowarning) {
              printf("  Reducing collinear tolerance from %g to %g degree.\n",
                     b->collinear_ang_tol, new_ang_tol);
            }
            b->collinear_ang_tol = new_ang_tol;
            cos_collinear_ang_tol = cos(b->collinear_ang_tol / 180. * PI);
          
            // This segment can be recovered by a 2-3 flip.
            if (subsegstack != NULL) {
              // Add the missing segment back to the recovering list.
              subsegstack->newindex((void **) &paryseg);
              *paryseg = *misseg;
            }
            return 1;
          } else {
            if (!b->quiet && !b->nowarning) {
              printf("  Cannot reduce the current collinear tolerance (=%g degree).\n",
                     b->collinear_ang_tol);
            }
            idir = (int) SELF_INTERSECT;
            return 0;
          }
        } else {
          // This segment can be recovered by a 2-3 flip.
          if (subsegstack != NULL) {
            // Add the missing segment back to the recovering list.
            subsegstack->newindex((void **) &paryseg);
            *paryseg = *misseg;
          }
          return 1;
        }
      } else {
        // This edge (not a segment) can be recovered by a 2-3 flip.
        return 1;
      }
    } // if (pd == endpt)

    if (issubface(searchtet)) {
      if (misseg != NULL) {
        terminate_tet_core(this, 2);
        // Report a segment and a facet intersect.
        if (!b->quiet && !b->nowarning) {
          face fac; tspivot(searchtet, fac);
          int segidx = getfacetindex(*misseg);
          point p1 = segmentendpointslist[segidx*2];
          point p2 = segmentendpointslist[segidx*2+1];
          printf("Warning:  A segment and a facet exactly intersect.\n");
          printf("  segment  : [%d,%d].\n", pointmark(p1), pointmark(p2));
          printf("  facet triangle: [%d,%d,%d] tag(%d).\n",
                 pointmark(org(searchtet)), pointmark(dest(searchtet)),
                 pointmark(apex(searchtet)), shellmark(fac));
        }
        idir = (int) SELF_INTERSECT;
      }
      return 0;
    } // if (issubface(searchtet))

    for (i = 0; i < 3; i++) {
      pa = org(spintet);
      pb = dest(spintet);
      if (tri_edge_test(pa, pb, pd, startpt, endpt, NULL, 1, types, poss)) {
        break; // Found the edge.
      }
      enextself(spintet);
      eprevself(searchtet);
    }
    esymself(searchtet);
  }
  else { // dir == ACROSSEDGE;
    if (issubseg(searchtet)) {
      terminate_tet_core(this, 2);
      if (misseg != NULL) {
        // Report a self_intersection.
        //bool intersect_flag = false;
        //point nearpt = dest(searchtet);
        ivf.iloc = ONVERTEX;
        // report_seg_vertex_intersect(misseg, dest(searchtet), ONVERTEX);
        int segidx = getfacetindex(*misseg);
        point p1 = segmentendpointslist[segidx*2];
        point p2 = segmentendpointslist[segidx*2+1];
        face parsentseg;
        //sdecode(point2sh(nearpt), parsentseg);
        tsspivot1(searchtet, parsentseg);
        int segidx2 = getfacetindex(parsentseg);
        if (segidx2 != segidx) {
          if (!b->quiet && !b->nowarning) {
            point p3 = segmentendpointslist[segidx2*2];
            point p4 = segmentendpointslist[segidx2*2+1];
            printf("Warning:  Two line segments are %s crossing.\n",
                   ivf.iloc == NEARVERTEX ? "nearly" : "exactly");
            printf("  1st: [%d,%d].\n", pointmark(p1), pointmark(p2));
            printf("  2nd: [%d,%d].\n", pointmark(p3), pointmark(p4));
          }
          //intersect_flag = true;
        } else {
          if (ivf.iloc == ONVERTEX) {
            terminate_tet_core(this, 2); // This should not be possible.
          }
        }
        idir = (int) SELF_INTERSECT;
      } // if (misseg != NULL)
      return 0;
    }
  }

  if (!splitsegflag) {
    // Try to recover this segment by adding Steiner points near it.

    spintet = searchtet;
    n = 0; endi = -1;
    while (1) {
      // Check if the endpt appears in the star.
      if (apex(spintet) == endpt) {
        endi = n; // Remember the position of endpt.
      }
      n++; // Count a tet in the star.
      fnextself(spintet);
      if (spintet.tet == searchtet.tet) break;
    }

    if (endi > 0) {
      // endpt is also in the edge star
      // Get all tets in the edge star.
      abtets = new triface[n];
      spintet = searchtet;
      for (i = 0; i < n; i++) {
        abtets[i] = spintet;
        fnextself(spintet);
      }

      success = 0;

      if (dir == ACROSSFACE) {
        // Find a Steiner points inside the polyhedron.
        if (add_steinerpt_in_schoenhardtpoly(abtets, endi, splitsliverflag, 0)) {
          success = 1;
        }
      } else if (dir == ACROSSEDGE) {
        // PLC check.
        if (issubseg(searchtet)) {
          terminate_tet_core(this, 2);
        }
        if (n > 4) {
          // In this case, 'abtets' is separated by the plane (containing the
          //   two intersecting edges) into two parts, P1 and P2, where P1
          //   consists of 'endi' tets: abtets[0], abtets[1], ...,
          //   abtets[endi-1], and P2 consists of 'n - endi' tets:
          //   abtets[endi], abtets[endi+1], abtets[n-1].
          if (endi > 2) { // P1
            // There are at least 3 tets in the first part.
            if (add_steinerpt_in_schoenhardtpoly(abtets, endi, splitsliverflag, 0)) {
              success++;
            }
          }
          if ((n - endi) > 2) { // P2
            // There are at least 3 tets in the first part.
            if (add_steinerpt_in_schoenhardtpoly(&(abtets[endi]), n - endi, splitsliverflag, 0)) {
              success++;
            }
          }
        } else {
          // In this case, a 4-to-4 flip should be re-cover the edge [c,d].
          //   However, there will be invalid tets (either zero or negtive
          //   volume). Otherwise, [c,d] should already be recovered by the
          //   recoveredge() function.
        }
      } else {
        terminate_tet_core(this, 2);
      }

      delete [] abtets;

      if (success && (misseg != NULL)) {
        // Add the missing segment back to the recovering list.
        subsegstack->newindex((void **) &paryseg);
        *paryseg = *misseg;
      }
      
      if (success) {
        return 1;
      }
    } // if (endi > 0)

    return 0;
  } // if (!splitsegflag)

  if (b->verbose > 3) {
    printf("      Recover segment (%d, %d) by splitting it.\n",
           pointmark(startpt), pointmark(endpt));
  }
  steinerpt = NULL;

  if (b->addsteiner_algo > 0) { // -Y/1 or -Y/2
    if (add_steinerpt_in_segment(misseg, 3, idir)) {
      return 1;
    }
    if (idir == SELF_INTERSECT) {
      return 0;
    }
    sesymself(*misseg);
    if (add_steinerpt_in_segment(misseg, 3, idir)) {
      return 1;
    }
    sesymself(*misseg);
    if (idir == SELF_INTERSECT) {
      return 0;
    }
  }


  // Let the face [a,b,d] be the first intersecting face of the segment
  //   [startpt, endpt]. We add the interseting point.
  double ip[3], u;
  point2tetorg(startpt, searchtet);
  dir = finddirection(&searchtet, endpt);
  if (dir == ACROSSVERT) {
    if (dest(searchtet) == endpt) {
      // This edge exists.
      if (misseg != NULL) {
        // Add the missing segment back to the recovering list.
        subsegstack->newindex((void **) &paryseg);
        *paryseg = *misseg;
      }
      return 1;
    } else {
      // This should be a self-intersection.
      if (misseg != NULL) {
        terminate_tet_core(this, 2);
        // report_seg_vertex_intersect(misseg, dest(searchtet), ONVERTEX);
        idir = (int) SELF_INTERSECT;
      }
      return 0;
    }
  }

  enextself(searchtet);
  pa = org(searchtet);
  pb = dest(searchtet);
  pd = oppo(searchtet);

  // Calculate the intersection of the face [a,b,d] and the segment.
  //planelineint(pa, pb, pd, startpt, endpt, ip, &u);

  point fpt[3], ept[2];
  sort_3pts(pa, pb, pd, fpt);
  sort_2pts(startpt, endpt, ept);
  planelineint(fpt[0], fpt[1], fpt[2], ept[0], ept[1], ip, &u);

  if ((u > 0) && (u < 1)) {
    // Create a Steiner point.
    makepoint(&steinerpt, FREESEGVERTEX);
    for (i = 0; i < 3; i++) steinerpt[i] = ip[i];
    
    // for create_a_shorter_edge().
    setpoint2sh(steinerpt, sencode(*misseg));
    
    esymself(searchtet); // The crossing face/edge.
    spivot(*misseg, splitsh);
    if (dir == ACROSSFACE) {
      //ivf.iloc = (int) ONFACE;
      ivf.refineflag = 4; // Check if the crossing face is removed.
    } else {
      //ivf.iloc = (int) ONEDGE;
      ivf.refineflag = 8; // Check if the crossing edge is removed.
    }
    ivf.iloc = (int) OUTSIDE; // do point location.
    ivf.refinetet = searchtet; // The crossing face/edge.
    ivf.bowywat = 1;
    // ivf.lawson = 0;
    ivf.lawson = 2; // Recover Delaunay after inserting this vertex.
    ivf.rejflag = 0;
    ivf.chkencflag = 0;
    ivf.sloc = (int) ONEDGE;
    ivf.sbowywat = 1; // split surface mesh separately, new subsegments are
                      //   pushed into "subsegstack".
    ivf.splitbdflag = 0;
    ivf.validflag = 1;
    ivf.respectbdflag = 1;
    ivf.assignmeshsize = b->metric;

    if (insertpoint(steinerpt, &searchtet, &splitsh, misseg, &ivf)) {
      if (flipstack != NULL) {
        recoverdelaunay();
      }
      
      // Save this Steiner point (for removal).
      //   Re-use the array 'subvertstack'.
      subvertstack->newindex((void **) &parypt);
      *parypt = steinerpt;

      st_segref_count++;
      if (steinerleft > 0) steinerleft--;

      return 1;
    } else {
      // Check if this failure is due to a self-intersection.
      if ((ivf.iloc == ONVERTEX) || (ivf.iloc == NEARVERTEX)) {
        if (misseg != NULL) {
          // report_seg_vertex_intersect(misseg, nearpt, ivf.iloc);
          int segidx = getfacetindex(*misseg);
          point p1 = segmentendpointslist[segidx*2];
          point p2 = segmentendpointslist[segidx*2+1];
          bool intersect_flag = false;
          point nearpt = org(searchtet);
          if (!issteinerpoint(nearpt)) {
            // 'nearpt' is an input vertex.
            if (!b->quiet && !b->nowarning) {
              point tmppt = NULL;
              if (is_segment(p1, nearpt)) tmppt = p1;
              else if (is_segment(p2, nearpt)) tmppt = p2;
              if (tmppt != NULL) {
                // Two input segments are nearly overlapping.
                printf("Warning:  Two line segments are %s overlapping.\n",
                       ivf.iloc == NEARVERTEX ? "nearly" : "exactly");
                printf("  1st: [%d,%d].\n", pointmark(p1), pointmark(p2));
                printf("  2nd: [%d,%d].\n", pointmark(tmppt), pointmark(nearpt));
              } else {
                // An input vertex is very close to a segment.
                printf("Warning:  A vertex lies %s on a line segment.\n",
                       ivf.iloc == NEARVERTEX ? "nearly" : "exactly");
                printf("  segment: [%d,%d].\n", pointmark(p1), pointmark(p2));
                printf("  vertex : [%d].\n", pointmark(nearpt));
              }
            } // if (!b->quiet && !b->nowarning)
            intersect_flag = true;
          } else {
            if (pointtype(nearpt) == FREESEGVERTEX) {
              // Check if two segments are nearly intersecting.
              face parsentseg;
              sdecode(point2sh(nearpt), parsentseg);
              int segidx2 = getfacetindex(parsentseg);
              if (segidx2 != segidx) {
                point p3 = segmentendpointslist[segidx2*2];
                point p4 = segmentendpointslist[segidx2*2+1];
                printf("Warning:  Two line segments are %s crossing.\n",
                       ivf.iloc == NEARVERTEX ? "nearly" : "exactly");
                printf("  1st: [%d,%d].\n", pointmark(p1), pointmark(p2));
                printf("  2nd: [%d,%d].\n", pointmark(p3), pointmark(p4));
                intersect_flag = true;
              }
            } else {
              // report other cases.
              // to do...
              terminate_tet_core(this, 2);
            }
          }
          if (intersect_flag) {
            if (!b->quiet && !b->nowarning) {
              if (ivf.iloc == NEARVERTEX) {
                double dd = distance(steinerpt, nearpt);
                double new_dd = minedgelength - dd / longest;
                double new_eps = new_dd / longest;
                printf("You can ignore this warning by using -T%e (default is %e) option.\n",
                       new_eps, b->epsilon);
                printf("  This will allow a short edge (len = %.17g) (default limit is %g).\n",
                       dd, minedgelength);
              }
            }
            // A self-intersection is detected.
            idir = (int) SELF_INTERSECT;
          } // if (intersect_flag)
        } // if (misseg != NULL)
      } // if ((ivf.iloc == ONVERTEX) || (ivf.iloc == NEARVERTEX))

      // The vertex is not inserted.
      pointdealloc(steinerpt);
      steinerpt = NULL;
    }
  } // if ((u > 0) && (u < 1))

  return 0; // Failed to reocver this segment.

  // [2020-05-02] The following code is skipped.
  if (steinerpt == NULL) {
    // Split the segment at its midpoint.
    makepoint(&steinerpt, FREESEGVERTEX);
    for (i = 0; i < 3; i++) {
      steinerpt[i] = 0.5 * (startpt[i] + endpt[i]);
    }

    // We need to locate the point.
    spivot(*misseg, splitsh);
    ivf.iloc = (int) OUTSIDE;
    ivf.bowywat = 1;
    //ivf.lawson = 0;
    ivf.lawson = 2; // do flip to recover locally Delaunay faces.
    ivf.rejflag = 0;
    ivf.chkencflag = 0;
    ivf.sloc = (int) ONEDGE;
    ivf.sbowywat = 1; // mesh surface separately
    ivf.splitbdflag = 0;
    ivf.validflag = 1;
    ivf.respectbdflag = 1;
    ivf.assignmeshsize = b->metric; 
    if (insertpoint(steinerpt, &searchtet, &splitsh, misseg, &ivf)) {
      if (flipstack != NULL) {
        recoverdelaunay();
      }
    } else {
      terminate_tet_core(this, 2);
    }
  } // if (endi > 0)

  // Save this Steiner point (for removal).
  //   Re-use the array 'subvertstack'.
  subvertstack->newindex((void **) &parypt);
  *parypt = steinerpt;

  st_segref_count++;
  if (steinerleft > 0) steinerleft--;

  return 1;
}

//============================================================================//
//                                                                            //
// recoversegments()    Recover all segments.                                 //
//                                                                            //
// All segments need to be recovered are in 'subsegstack'.                    //
//                                                                            //
// This routine first tries to recover each segment by only using flips. If   //
// no flip is possible, and the flag 'steinerflag' is set, it then tries to   //
// insert Steiner points near or in the segment.                              //
//                                                                            //
//============================================================================//

int TetMeshCore::recoversegments(arraypool *misseglist, int fullsearch,
                                int steinerflag)
{
  triface searchtet, spintet;
  face sseg, *paryseg;
  point startpt, endpt;
  int success, idir;
  int t1ver;

  long bak_inpoly_count = st_volref_count; 
  long bak_segref_count = st_segref_count;

  if (b->verbose > 1) {
    printf("    Recover segments [%s level = %2d] #:  %ld.\n",
           (b->fliplinklevel > 0) ? "fixed" : "auto",
           (b->fliplinklevel > 0) ? b->fliplinklevel : autofliplinklevel,
           subsegstack->objects);
  }

  // Loop until 'subsegstack' is empty.
  while (subsegstack->objects > 0l) {
    // seglist is used as a stack.
    subsegstack->objects--;
    paryseg = (face *) fastlookup(subsegstack, subsegstack->objects);
    sseg = *paryseg;

    // Check if this segment has been recovered.
    sstpivot1(sseg, searchtet);
    if (searchtet.tet != NULL) {
      continue; // Not a missing segment.
    }

    startpt = sorg(sseg);
    endpt = sdest(sseg);

    if (b->verbose > 2) {
      printf("      Recover segment (%d, %d).\n", pointmark(startpt), 
             pointmark(endpt));
    }

    success = 0;

    if (recoveredgebyflips(startpt, endpt, &sseg, &searchtet, 0, idir)) {
      success = 1;
    } else {
      // Try to recover it from the other direction.
      if ((idir != (int) SELF_INTERSECT) &&
          recoveredgebyflips(endpt, startpt, &sseg, &searchtet, 0, idir)) {
        success = 1;
      }
    }

    
    if (!success && fullsearch) {
      if ((idir != (int) SELF_INTERSECT) &&
          recoveredgebyflips(startpt, endpt, &sseg, &searchtet, fullsearch, idir)) {
        success = 1;
      }
    }
    
    if (success) {
      // Segment is recovered. Insert it.
      // Let the segment remember an adjacent tet.
      sstbond1(sseg, searchtet);
      // Bond the segment to all tets containing it.
      spintet = searchtet;
      do {
        tssbond1(spintet, sseg);
        fnextself(spintet);
      } while (spintet.tet != searchtet.tet);
    } else {
      if ((idir != (int) SELF_INTERSECT) && (steinerflag > 0)) {
        // Try to recover the segment but do not split it.
        if (add_steinerpt_to_recover_edge(startpt, endpt, &sseg, 0, 0, idir)) {
          success = 1;
        }
        if (!success && (idir != (int) SELF_INTERSECT) && (steinerflag > 1)) {
          // Split the segment.
          if (add_steinerpt_to_recover_edge(startpt, endpt, &sseg, 1, 0, idir)) {
            success = 1;
          }
        }
      }
    
      if (!success) {
        if (idir != (int) SELF_INTERSECT) {
          if (misseglist != NULL) {
            // Save this segment (recover it later).
            misseglist->newindex((void **) &paryseg);
            *paryseg = sseg;
          }
        } else {
          // Save this segment (do not recover it again).
          if (skipped_segment_list == NULL) {
            skipped_segment_list = new arraypool(sizeof(badface), 10);
          }
          badface *bf;
          skipped_segment_list->newindex((void **) &bf);
          bf->init();
          bf->ss = sseg;
          bf->forg  =  sorg(sseg);
          bf->fdest = sdest(sseg);
          bf->key = (double) shellmark(sseg);
          smarktest3(sseg);
          // Save all subfaces at this segment, do not recover them later.
          if (skipped_facet_list == NULL) {
            skipped_facet_list = new arraypool(sizeof(badface), 10);
          }
          face neighsh, spinsh;
          bf->ss.shver = 0;
          spivot(bf->ss, neighsh);
          spinsh = neighsh;
          while (spinsh.sh != NULL) {
            skipped_facet_list->newindex((void **) &bf);
            bf->init();
            bf->ss = spinsh;
            bf->forg  = (point) spinsh.sh[3];
            bf->fdest = (point) spinsh.sh[4];
            bf->fapex = (point) spinsh.sh[5];
            bf->key = (double) shellmark(spinsh);
            smarktest3(spinsh); // do not recover it.
            spivotself(spinsh);
            if (spinsh.sh == neighsh.sh) break;
          }
        }
      } // if (!success)
    }

  } // while (subsegstack->objects > 0l)

  if (steinerflag) {
    if (b->verbose > 1) {
      // Report the number of added Steiner points.
      if (st_volref_count > bak_inpoly_count) {
        printf("    Add %ld Steiner points in volume.\n", 
               st_volref_count - bak_inpoly_count);
      }
      if (st_segref_count > bak_segref_count) {
        printf("    Add %ld Steiner points in segments.\n", 
               st_segref_count - bak_segref_count);
      }
    }
  }

  return 0;
}

//============================================================================//
//                                                                            //
// recoverfacebyflips()    Recover a face by flips.                           //
//                                                                            //
// 'pa', 'pb', and 'pc' are the three vertices of this face.  This routine    //
// tries to recover it in the tetrahedral mesh. It is assumed that the three  //
// edges, i.e., pa->pb, pb->pc, and pc->pa all exist.                         //
//                                                                            //
// If the face is recovered, it is returned by 'searchtet'.                   //
//                                                                            //
// If 'searchsh' is not NULL, it is a subface to be recovered.  Its vertices  //
// must be pa, pb, and pc.  It is mainly used to check self-intersections.    //
// Another use of this subface is to split it when a Steiner point is found   //
// inside this subface.                                                       //
//                                                                            //
//============================================================================//

int TetMeshCore::recoverfacebyflips(point pa, point pb, point pc, 
                                   face *searchsh, triface* searchtet,
                                   int &dir, point *p1, point *p2)
{
  triface spintet, flipedge;
  point pd, pe;
  flipconstraints fc;
  int types[2], poss[4], intflag;
  int success;
  int t1ver; 
  int i, j;


  fc.fac[0] = pa;
  fc.fac[1] = pb;
  fc.fac[2] = pc;
  fc.checkflipeligibility = 1;
  
  dir = (int) DISJOINT;
  success = 0;

  for (i = 0; i < 3 && !success; i++) {
    while (1) {
      // Get a tet containing the edge [a,b].
      point2tetorg(fc.fac[i], *searchtet);
      finddirection(searchtet, fc.fac[(i+1)%3]);
      // Search the face [a,b,c]
      spintet = *searchtet;
      while (1) {
        if (apex(spintet) == fc.fac[(i+2)%3]) {
          // Found the face.
          *searchtet = spintet;
          // Return the face [a,b,c].
          for (j = i; j > 0; j--) {
            eprevself(*searchtet);
          }
          dir = (int) SHAREFACE;
          success = 1;
          break;
        }
        fnextself(spintet);
        if (spintet.tet == searchtet->tet) break;
      } // while (1)

      if (success) break;

      // The face is missing. Try to recover it.
      flipedge.tet = NULL;
      // Find a crossing edge of this face.
      spintet = *searchtet;
      while (1) {
        pd = apex(spintet);
        pe = oppo(spintet);
        if ((pd != dummypoint) && (pe != dummypoint)) {
          // Check if [d,e] intersects [a,b,c]
          intflag = tri_edge_test(pa, pb, pc, pd, pe, NULL, 1, types, poss);
          if (intflag > 0) {
            // By the assumption that all edges of the face exist, they can
            //   only intersect at a single point.
            if (intflag == 2) {
              // Go to the edge [d,e].
              edestoppo(spintet, flipedge); // [d,e,a,b]
              if (searchsh != NULL) {
                // Check the intersection type.
                dir = types[0]; // return this value.
                if ((types[0] == (int) ACROSSFACE) || 
                    (types[0] == (int) ACROSSEDGE)) {
                  // Check if [e,d] is a segment.
                  if (issubseg(flipedge)) {
                    // This subface intersects with a segment.
                    if (!b->quiet && !b->nowarning) {
                      if (!b->quiet && !b->nowarning) {
                        printf("Warning:  A segment and a facet intersect.\n");
                        face sseg; tsspivot1(flipedge, sseg);
                        int segidx = getfacetindex(sseg);
                        point p1 = segmentendpointslist[segidx*2];
                        point p2 = segmentendpointslist[segidx*2+1];
                        printf("  segment: [%d,%d] tag(%d).\n",
                               pointmark(p1), pointmark(p2), shellmark(sseg));
                        point *ppt = (point *) &(searchsh->sh[3]);
                        printf("  facet triangle: [%d,%d,%d] tag(%d)\n",
                               pointmark(ppt[0]), pointmark(ppt[1]),
                               pointmark(ppt[2]), shellmark(*searchsh));
                      }
                    }
                    dir = (int) SELF_INTERSECT;
                    return 0; // Found a self-intersection.
		          } else {
				    // Check if [e,d] is an edge of a subface.
					triface chkface = flipedge;
					while (1) {
					  if (issubface(chkface)) break;
					  fsymself(chkface);
					  if (chkface.tet == flipedge.tet) break;
					}
					if (issubface(chkface)) {
                      if (searchsh != NULL) {
                        // Two subfaces are intersecting.
                        if (!b->quiet && !b->nowarning) {
                          printf("Warning:  Found two facets intersect.\n");
                          point *ppt = (point *) &(searchsh->sh[3]);
                          printf("  1st facet triangle: [%d,%d,%d] tag(%d)\n",
                                 pointmark(ppt[0]), pointmark(ppt[1]),
                                 pointmark(ppt[2]), shellmark(*searchsh));
                          face fa; tspivot(chkface, fa);
                          ppt = (point *) &(fa.sh[3]);
                          printf("  2nd facet triangle: [%d,%d,%d] tag(%d)\n",
                                 pointmark(ppt[0]), pointmark(ppt[1]),
                                 pointmark(ppt[2]), shellmark(fa));
                        }
                        dir = (int) SELF_INTERSECT;
                      }
                      return 0; // Found a self-intersection.
					}
				  }
				} else if (types[0] == TOUCHFACE) {
                  // This is possible when a Steiner point was added on it.
				  point touchpt, *parypt;
                  if (poss[1] == 0) {
                    touchpt = pd; // pd is a coplanar vertex.
                  } else {
                    touchpt = pe; // pe is a coplanar vertex.
                  }
                  if (!issteinerpoint(touchpt)) {
                    if (!b->quiet && !b->nowarning) {
                      printf("Warning:  A vertex lies on a facet.\n");
                      printf("  vertex : [%d]\n", pointmark(touchpt));
                      point *ppt = (point *) &(searchsh->sh[3]);
                      printf("  facet triangle: [%d,%d,%d] tag(%d)\n",
                             pointmark(ppt[0]), pointmark(ppt[1]),
                             pointmark(ppt[2]), shellmark(*searchsh));
                    }
                    dir = (int) SELF_INTERSECT;
                    return 0;
                  } else if (pointtype(touchpt) == FREESEGVERTEX) {
                    if (!b->quiet && !b->nowarning) {
                      printf("Warning:  A segment and a facet intersect.\n");
                      face sseg;
                      sdecode(point2sh(touchpt), sseg);
                      int segidx = getfacetindex(sseg);
                      point p1 = segmentendpointslist[segidx*2];
                      point p2 = segmentendpointslist[segidx*2+1];
                      printf("  segment: [%d,%d] tag(%d).\n",
                             pointmark(p1), pointmark(p2), shellmark(sseg));
                      point *ppt = (point *) &(searchsh->sh[3]);
                      printf("  facet triangle: [%d,%d,%d] tag(%d)\n",
                             pointmark(ppt[0]), pointmark(ppt[1]),
                             pointmark(ppt[2]), shellmark(*searchsh));
                    }
                    dir = (int) SELF_INTERSECT;
                    return 0;
                  } else if (pointtype(touchpt) == FREEFACETVERTEX) {
                    if (!b->quiet && !b->nowarning) {
                      printf("Warning:  Found two facets intersect.\n");
                      point *ppt = (point *) &(searchsh->sh[3]);
                      printf("  1st facet triangle: [%d,%d,%d] tag(%d)\n",
                             pointmark(ppt[0]), pointmark(ppt[1]),
                             pointmark(ppt[2]), shellmark(*searchsh));
                      face fa;
                      sdecode(point2sh(touchpt), fa);
                      ppt = (point *) &(fa.sh[3]);
                      printf("  2nd facet triangle: [%d,%d,%d] tag(%d)\n",
                             pointmark(ppt[0]), pointmark(ppt[1]),
                             pointmark(ppt[2]), shellmark(fa));
                    }
                    dir = (int) SELF_INTERSECT;
                    return 0;
                  } else if (pointtype(touchpt) == FREEVOLVERTEX) {
                    // A volume Steiner point was added in this subface.
                    // Split this subface by this point.
                    face checksh, *parysh;
                    int siloc = (int) ONFACE;
                    int sbowat = 0; // Only split this subface. A 1-to-3 flip.
                    setpointtype(touchpt, FREEFACETVERTEX);
                    sinsertvertex(touchpt, searchsh, NULL, siloc, sbowat, 0);
                    st_volref_count--;
                    st_facref_count++;
                    // Queue this vertex for removal.
                    subvertstack->newindex((void **) &parypt);
                    *parypt = touchpt;
                    // Queue new subfaces for recovery.
                    // Put all new subfaces into stack for recovery.
                    for (i = 0; i < caveshbdlist->objects; i++) {
                      // Get an old subface at edge [a, b].
                      parysh = (face *) fastlookup(caveshbdlist, i);
                      spivot(*parysh, checksh); // The new subface [a, b, p].
                      // Do not recover a deleted new face (degenerated).
                      if (checksh.sh[3] != NULL) {
                        subfacstack->newindex((void **) &parysh);
                        *parysh = checksh;
                      }
                    }
                    // Delete the old subfaces in sC(p).
                    for (i = 0; i < caveshlist->objects; i++) {
                      parysh = (face *) fastlookup(caveshlist, i);
                      shellfacedealloc(subfaces, parysh->sh);
                    }
                    // Clear working lists.
                    caveshlist->restart();
                    caveshbdlist->restart();
                    cavesegshlist->restart();
                    // We can return this function.
                    searchsh->sh = NULL; // It has been split.
					return 1;
                  } else {
				    // Other cases may be due to a bug or a PLC error.
					//return report_selfint_face(pa, pb, pc, searchsh, &flipedge,
                    //                           intflag, types, poss);
                    terminate_tet_core(this, 2); // to debug...
                    dir = (int) SELF_INTERSECT;
                    return 0; // Found a self-intersection.
				  }
                } else {
                  // The other intersection types: ACROSSVERT, TOUCHEDGE, 
                  // SHAREVERTEX should not be possible or due to a PLC error.
                  //return report_selfint_face(pa, pb, pc, searchsh, &flipedge,
                  //                           intflag, types, poss);
                  terminate_tet_core(this, 2); // to report
                  dir = (int) SELF_INTERSECT;
                  return 0;
                }
              } // if (searchsh != NULL)
            } else { // intflag == 4. Coplanar case.
              // Found a mesh edge is coplanar with this subface.
              // It migh be caused by a self-intersection.
              terminate_tet_core(this, 2); // report this bug
            }
            break;
          } // if (intflag > 0)
        }
        fnextself(spintet);
        if (spintet.tet == searchtet->tet) {
          terminate_tet_core(this, 2);
        }
      } // while (1)
      // Try to flip the edge [d,e].
      // Remember a crossing edge.
      *p1 =  org(flipedge);
      *p2 = dest(flipedge);

      if (removeedgebyflips(&flipedge, &fc) == 2) {
        // A crossing edge is removed.
        continue; 
      }
      
      // Unable to remove a crossing edge of this face.
      break;
    } // while (1)
  } // i

  return success;
}

 
//============================================================================//
//                                                                            //
// recoversubfaces()    Recover all subfaces.                                 //
//                                                                            //
//============================================================================//

int TetMeshCore::recoversubfaces(arraypool *misshlist, int steinerflag)
{
  triface searchtet, neightet, spintet;
  face searchsh, neighsh, neineish, *parysh;
  face bdsegs[3];
  point startpt, endpt, apexpt, *parypt;
  point cross_e1 = NULL, cross_e2 = NULL; // endpoints of a crossing edge.
  point steinerpt;
  insertvertexflags ivf;
  int success, dir;
  int t1ver;
  int i, j;

  if (b->verbose > 1) {
    printf("    Recover subfaces [%s level = %2d] #:  %ld.\n",
           (b->fliplinklevel > 0) ? "fixed" : "auto",
           (b->fliplinklevel > 0) ? b->fliplinklevel : autofliplinklevel,
           subfacstack->objects);
  }

  // Loop until 'subfacstack' is empty.
  while (subfacstack->objects > 0l) {

    subfacstack->objects--;
    parysh = (face *) fastlookup(subfacstack, subfacstack->objects);
    searchsh = *parysh;

    if (searchsh.sh[3] == NULL) continue; // Skip a dead subface.
    if (smarktest3ed(searchsh)) continue; // Skip a self-intersected subface.

    stpivot(searchsh, neightet);
    if (neightet.tet != NULL) continue; // Skip a recovered subface.

    if (b->verbose > 2) {
      printf("      Recover subface (%d, %d, %d).\n",pointmark(sorg(searchsh)),
             pointmark(sdest(searchsh)), pointmark(sapex(searchsh)));
    }
    dir = (int) DISJOINT; // No self intersection is detected.

    // The three edges of the face need to be existed first.
    for (i = 0; i < 3; i++) {
      sspivot(searchsh, bdsegs[i]);
      if (bdsegs[i].sh != NULL) {
        // Check if this segment exist.
        sstpivot1(bdsegs[i], searchtet);
        if (searchtet.tet == NULL) {
          // This segment is not recovered yet. Try to recover it.
          success = 0;
          startpt = sorg(searchsh);
          endpt = sdest(searchsh);
          if (recoveredgebyflips(startpt, endpt, &bdsegs[i], &searchtet, 0, dir)) {
            success = 1;
          } else {
            if ((dir != (int) SELF_INTERSECT) &&
                recoveredgebyflips(endpt, startpt, &bdsegs[i], &searchtet, 0, dir)) {
              success = 1;
            }
          }
          if (success) {
            // Segment is recovered. Insert it.
            // Let the segment remember an adjacent tet.
            sstbond1(bdsegs[i], searchtet);
            // Bond the segment to all tets containing it.
            spintet = searchtet;
            do {
              tssbond1(spintet, bdsegs[i]);
              fnextself(spintet);
            } while (spintet.tet != searchtet.tet);
          } else {
            // An edge of this subface is missing. Can't recover this subface.
            // Delete any temporary segment that has been created.
            for (j = (i - 1); j >= 0; j--) {
              if (smarktest2ed(bdsegs[j])) {
                spivot(bdsegs[j], neineish);
                ssdissolve(neineish);
                spivot(neineish, neighsh);
                if (neighsh.sh != NULL) {
                  ssdissolve(neighsh);
                }
                sstpivot1(bdsegs[j], searchtet);
                spintet = searchtet;
                while (1) {
                  tssdissolve1(spintet);
                  fnextself(spintet);
                  if (spintet.tet == searchtet.tet) break;
                }
                shellfacedealloc(subsegs, bdsegs[j].sh);
              }
            } // j
            break; // i
          } // if (success) else
        } // if (searchtet.tet == NULL)
      } else {
        // This edge is not a segment.
        // Check whether it exists or not.
        success = 0;
        startpt = sorg(searchsh);
        endpt = sdest(searchsh);
        point2tetorg(startpt, searchtet);
        finddirection(&searchtet, endpt);
        if (dest(searchtet) == endpt) {
          success = 1; // Found this edge.
        } else {
          // The edge is missing. Try to recover it.
          if (recoveredgebyflips(startpt, endpt, &searchsh, &searchtet, 0, dir)) {
            success = 1;
          } else {
            if ((dir != (int) SELF_INTERSECT) &&
                recoveredgebyflips(endpt, startpt, &searchsh, &searchtet, 0, dir)) {
              success = 1;
            }
          }
        }

        if (success) {
          // This edge exists.
          if (issubseg(searchtet)) {
            // A segment already exists at this edge!
            //terminate_tet_core(this, 2); // to debug
            //dir = SELF_INTERSECT;
            // We contnue to recover this subface instead of reporting a
            //   SELF_INTERSECT event.
            // Eventually, we will find "a duplicated triangle" event.
          }
        }
        
        if (success && (dir != SELF_INTERSECT)) {
          // This edge exists.
          //if (!issubseg(searchtet)) {
            // Insert a temporary segment to protect this edge.
            makeshellface(subsegs, &(bdsegs[i]));
            setshvertices(bdsegs[i], startpt, endpt, NULL);
            smarktest2(bdsegs[i]); // It's a temporary segment.
            // Insert this segment into surface mesh.
            ssbond(searchsh, bdsegs[i]);
            spivot(searchsh, neighsh);
            if (neighsh.sh != NULL) {
              ssbond(neighsh, bdsegs[i]);
            }
            // Insert this segment into tetrahedralization.
            sstbond1(bdsegs[i], searchtet);
            // Bond the segment to all tets containing it.
            spintet = searchtet;
            do {
              tssbond1(spintet, bdsegs[i]);
              fnextself(spintet);
            } while (spintet.tet != searchtet.tet);
          //}
        } else {
          // An edge of this subface is missing. Can't recover this subface.
          // Delete any temporary segment that has been created.
          for (j = (i - 1); j >= 0; j--) {
            if (smarktest2ed(bdsegs[j])) { 
              spivot(bdsegs[j], neineish);
                ssdissolve(neineish);
                spivot(neineish, neighsh);
                if (neighsh.sh != NULL) {
                  ssdissolve(neighsh);
                }
              sstpivot1(bdsegs[j], searchtet);
                spintet = searchtet;
                while (1) {
                  tssdissolve1(spintet);
                  fnextself(spintet);
                  if (spintet.tet == searchtet.tet) break;
                }
              shellfacedealloc(subsegs, bdsegs[j].sh);
            }
          } // j

          break;
        }
      }
      senextself(searchsh);
    } // i

    if (i == 3) {
      // All edges of this subface exist (or have been recovered).
      // Recover the subface.
      startpt = sorg(searchsh);
      endpt   = sdest(searchsh);
      apexpt  = sapex(searchsh);

      success = recoverfacebyflips(startpt, endpt, apexpt,&searchsh, &searchtet,
                                   dir, &cross_e1, &cross_e2);

      // Delete any temporary segment that has been created.
      for (j = 0; j < 3; j++) {
        if (smarktest2ed(bdsegs[j])) { 
          spivot(bdsegs[j], neineish);
            ssdissolve(neineish);
            spivot(neineish, neighsh);
            if (neighsh.sh != NULL) {
              ssdissolve(neighsh);
            }
          sstpivot1(bdsegs[j], neightet);
            spintet = neightet;
            while (1) {
              tssdissolve1(spintet);
              fnextself(spintet);
              if (spintet.tet == neightet.tet) break;
            }
          shellfacedealloc(subsegs, bdsegs[j].sh);
        }
      } // j

      if (success) {
        if (searchsh.sh != NULL) {
          // Face is recovered. Insert it.
          face chkface;
          tspivot(searchtet, chkface);
          if (chkface.sh == NULL) {
            tsbond(searchtet, searchsh);
            fsymself(searchtet);
            sesymself(searchsh);
            tsbond(searchtet, searchsh);
          } else {
            // A duplicated facet is found.
            if (shellmark(chkface) == shellmark(searchsh)) {
              if (!b->quiet && !b->nowarning) {
                point *ppt = (point *) &(searchsh.sh[3]);
                printf("Warning:  A duplicated triangle (%d,%d,%d) tag(%d) is ignored.\n",
                       pointmark(ppt[0]), pointmark(ppt[1]), pointmark(ppt[2]),
                       shellmark(searchsh));
              }
              duplicated_facets_count++;
              smarktest3(searchsh); // do not recover it.
              sinfect(searchsh); // it is an igonred duplicated facet.
            } else {
              if (!b->quiet && !b->nowarning) {
                point *ppt = (point *) &(chkface.sh[3]);
                printf("Warning:  Two facets are overlapping at triangle (%d,%d,%d).\n",
                       pointmark(ppt[0]), pointmark(ppt[1]), pointmark(ppt[2]));
                printf("  1st facet tag(%d).\n", shellmark(chkface));
                printf("  2nd facet tag(%d).\n", shellmark(searchsh));
              }
              dir = SELF_INTERSECT;
              success = 0;
            }
          }
        }
      } else {
        if ((dir != (int) SELF_INTERSECT) && steinerflag) {
          // Add a Steiner point at the barycenter of this subface.
          double ip[3], u;

          //planelineint(startpt, endpt, apexpt, cross_e1, cross_e2, ip, &u);

          point fpt[3], ept[2];
          sort_3pts(startpt, endpt, apexpt, fpt);
          sort_2pts(cross_e1, cross_e2, ept);
          planelineint(fpt[0], fpt[1], fpt[2], ept[0], ept[1], ip, &u);

          makepoint(&steinerpt, FREEFACETVERTEX);
          if ((u > 0.) && (u < 1.)) {
            for (j = 0; j < 3; j++) steinerpt[j] = ip[j];
            // Make sure that this Steiner point is inside the subface.
            if (is_collinear_at(steinerpt, startpt, endpt) ||
                is_collinear_at(steinerpt, endpt, apexpt) ||
                is_collinear_at(steinerpt, apexpt, startpt)) {
              // Add the barycenter of this missing subface.
              for (j = 0; j < 3; j++) {
                steinerpt[j] = (startpt[j] + endpt[j] + apexpt[j]) / 3.0;
              }
              // Avoid creating a very skinny triangle
              if (is_collinear_at(steinerpt, startpt, endpt) ||
                  is_collinear_at(steinerpt, endpt, apexpt) ||
                  is_collinear_at(steinerpt, apexpt, startpt)) {
                terminate_tet_core(this, 2);
              }
            }
          } else {
            // Add the barycenter of this missing subface.
            for (j = 0; j < 3; j++) {
              steinerpt[j] = (startpt[j] + endpt[j] + apexpt[j]) / 3.0;
            }
            // Avoid creating a very skinny triangle
            if (is_collinear_at(steinerpt, startpt, endpt) ||
                is_collinear_at(steinerpt, endpt, apexpt) ||
                is_collinear_at(steinerpt, apexpt, startpt)) {
              //assert(0); // to debug...
              terminate_tet_core(this, 2);
            }
          }

          // for create_a_shorter_edge().
          setpoint2sh(steinerpt, sencode(searchsh));

          ivf.init();
          point2tetorg(startpt, searchtet); // Start from 'searchtet'.
          ivf.iloc = (int) OUTSIDE; // Need point location.
          ivf.bowywat = 1;
          ivf.lawson = 2; // do recover delaunay.
          ivf.rejflag = 0;
          ivf.chkencflag = 0;
          ivf.sloc = (int) ONFACE; // "searchsh" must be the subface.
          ivf.sbowywat = 1; // split subface mesh separately, new subfaces
                            // are pushed into "subfacestack".
          ivf.splitbdflag = 0;
          ivf.validflag = 1;
          ivf.respectbdflag = 1;
          ivf.assignmeshsize = b->metric;

          if (insertpoint(steinerpt, &searchtet, &searchsh, NULL, &ivf)) {
            if (flipstack != NULL) {
              recoverdelaunay();
            }
              
            // Save this Steiner point (for removal).
            //   Re-use the array 'subvertstack'.
            subvertstack->newindex((void **) &parypt);
            *parypt = steinerpt;

            st_facref_count++;
            if (steinerleft > 0) steinerleft--;
              
            success = 1; // This subface has been split.
          } else {
            // Failed to insert this point.
            if (ivf.iloc == NEARVERTEX) {
              // Check if this subface is nearly "touched" by an existing
              //   vertex. If so, report an event.
              point chkpt = org(searchtet);
              double dist = distance(steinerpt, chkpt);
              if (dist < minedgelength) {
                if (!issteinerpoint(chkpt)) {
                  if (!b->quiet && !b->nowarning) { // -no -Q -W
                    printf("Warning:  A facet (%d,%d,%d) and a vertex %d are very close.\n",
                           pointmark(sorg(searchsh)), pointmark(sdest(searchsh)),
                           pointmark(sapex(searchsh)), pointmark(chkpt));
                    double dd = dist; // distance(steinerpt, nearpt);
                    //assert(dd > 0.);
                    //minedgelength = longest * b->epsilon;
                    //assert(dd < minedgelength);
                    double new_dd = minedgelength - dd / longest;
                    double new_eps = new_dd / longest;
                    printf("You can ignore this warning by using -T%e (default is %e) option.\n",
                           new_eps, b->epsilon);
                    printf("  This will allow a short edge (len = %.17g) (default limit is %g).\n",
                           dd, minedgelength);
                  }
                  dir = SELF_INTERSECT;
                }
              } else {
                // Report other types of possible (nearly) self-intersection.
                terminate_tet_core(this, 2);
                dir = SELF_INTERSECT;
              }
            }

            if ((dir != SELF_INTERSECT) && (steinerflag >= 2)) {
              if (ivf.iloc == NULLCAVITY) {
                // Collect a list of bad quality tets which prevent the
                //   insertion of this Steiner point.
                terminate_tet_core(this, 2);
                point2tetorg(startpt, searchtet);
                ivf.iloc = (int) OUTSIDE; // re-do point location.
                ivf.collect_inial_cavity_flag = 1;
                insertpoint(steinerpt, &searchtet, &searchsh, NULL, &ivf);
              } else {
                terminate_tet_core(this, 2); // report a bug.
              }
            } // if (steinerflag >= 2)

            pointdealloc(steinerpt);
            steinerpt = NULL;
            success = 0; // queue this subface.
          }
        } // if (steinerflag)
      }
    } else { // when i < 3
      // An edge (startpt, endpt) of this subface is missing.
      if ((dir != (int) SELF_INTERSECT) && (steinerflag > 0)) {
        // Split this edge by adding a Steiner point.
        // Find the first face/edge crossed by the edge (startpt, endpt).
        point2tetorg(startpt, searchtet);
        dir = finddirection(&searchtet, endpt);


        if (dir != (int) SELF_INTERSECT) {
          // Insert a Steiner point.
          double ip[3], u;

          enextself(searchtet);
          point pa =  org(searchtet);
          point pb = dest(searchtet);
          point pd = oppo(searchtet);

          //planelineint(pa, pb, pd, startpt, endpt, ip, &u);
          
          point fpt[3], ept[2];
          sort_3pts(pa, pb, pd, fpt);
          sort_2pts(startpt, endpt, ept);
          planelineint(fpt[0], fpt[1], fpt[2], ept[0], ept[1], ip, &u);

          makepoint(&steinerpt, FREEFACETVERTEX);
          for (j = 0; j < 3; j++) steinerpt[j] = ip[j];

          ivf.init();

          ivf.refinetet = searchtet; // bakup the crossing face/edge.

          triface tmptet = searchtet;
          ivf.iloc = locate(steinerpt, &tmptet);

          if (ivf.iloc == ONVERTEX) {
            // the origin of tmptet is co-incident with this Steiner point.
            searchtet = tmptet;
          }
          //else if (ivf.iloc == ONFACE) {
          //  searchtet = tmptet;
          //} else if (ivf.iloc == ONEDGE) {
          //  searchtet = tmptet;
          //}
          else {
            //assert(0); // to debug...
            // Make sure that we can split the crossing edge/face (a,b,d).
            if (dir == ACROSSFACE) {
              ivf.iloc = (int) ONFACE;
              //ivf.refineflag = 4; // Check if the crossing face is removed.
            } else if (dir == ACROSSEDGE) {
              ivf.iloc = (int) ONEDGE;
              //ivf.refineflag = 8; // Check if the crossing edge is removed.
            } else {
              terminate_tet_core(this, 2);
            }
            //ivf.iloc = (int) OUTSIDE; // do point location.
            //ivf.refinetet = searchtet; // The crossing face/edge.
          }

          ivf.bowywat = 1;
          ivf.lawson = 2; // do recover delaunay.
          ivf.rejflag = 0;
          ivf.chkencflag = 0;
          ivf.sloc = (int) ONEDGE; // "searchsh" must be the subedge.
          ivf.sbowywat = 1; // split subface mesh separately, new subfaces
                            // are pushed into "subfacestack".
          ivf.splitbdflag = 0;
          ivf.validflag = 1;
          ivf.respectbdflag = 1;
          ivf.assignmeshsize = b->metric;

          //if (steinerflag >= 2) {
            // Skip NEARVERTEX. This may create a very short edge.
            //ivf.ignore_near_vertex = 1;
          //}

          // searchsh may contain a missing segment.
          // After splitting this subface, this segment must also be split.
          //   the two missing subsegments are stored in "subsegstack".
          face misseg, *splitseg = NULL;
          sspivot(searchsh, misseg);
          if (misseg.sh != NULL) {
            splitseg = &misseg;
            setpointtype(steinerpt, FREESEGVERTEX); // default is FREEFACETVERTEX.
            // for create_a_shorter_edge()
            setpoint2sh(steinerpt, sencode(misseg));
          } else {
            // for create_a_shorter_edge()
            setpoint2sh(steinerpt, sencode(searchsh));
          }

          bool splitseg_flag = (splitseg != NULL);
          int bak_iloc = ivf.iloc; // for collect_initial_cavity

          if (insertpoint(steinerpt, &searchtet, &searchsh, splitseg, &ivf)) {
            if (flipstack != NULL) {
              recoverdelaunay();
            }

            // Save this Steiner point (for removal).
            //   Re-use the array 'subvertstack'.
            subvertstack->newindex((void **) &parypt);
            *parypt = steinerpt;

            if (splitseg_flag) {
              st_segref_count++;
            } else {
              st_facref_count++;
            }
            if (steinerleft > 0) steinerleft--;

            success = 1; // This subface has been split.
          } else {
            // Failed to insert this point.
            if (ivf.iloc == NEARVERTEX) {
              // Check if this subface is nearly "touched" by an existing
              //   vertex. If so, report an event.
              point chkpt = org(searchtet);
              double dist = distance(steinerpt, chkpt); // for reporting.
              if (!issteinerpoint(chkpt)) {
                if (!b->quiet && !b->nowarning) {
                  if (splitseg_flag) {
                    printf("Warning:  A segment (%d,%d) and a vertex %d are very close.\n",
                           pointmark(sorg(searchsh)), pointmark(sdest(searchsh)),
                           pointmark(chkpt));
                  } else {
                    printf("Warning:  A facet (%d,%d,%d) and a vertex %d are very close.\n",
                           pointmark(sorg(searchsh)), pointmark(sdest(searchsh)),
                           pointmark(sapex(searchsh)), pointmark(chkpt));
                  }
                  //printf("  Will result a vert short edge (len=%.17g) (< %.17g)\n",
                  //       dist, minedgelength);
                  double dd = dist; // distance(steinerpt, nearpt);
                  //assert(dd > 0.);
                  //minedgelength = longest * b->epsilon;
                  //assert(dd < minedgelength);
                  double new_dd = minedgelength - dd / longest;
                  double new_eps = new_dd / longest;
                  printf("You can ignore this warning by using -T%e (default is %e) option.\n",
                         new_eps, b->epsilon);
                  printf("  This will allow a short edge (len = %.17g) (default limit is %g).\n",
                         dd, minedgelength);
                }
                dir = SELF_INTERSECT;
              }
            }

            if ((dir != SELF_INTERSECT) && (steinerflag >= 2)) {
              success = 0; // queue this subface.
              // Failed to split a crossing edge/face.
              if ((ivf.iloc == ONVERTEX) || (ivf.iloc == NEARVERTEX)) {
                // Get the existing vertex (must be a Steiner point).
                if (dir == ACROSSEDGE) {
                  int idir;
                  if (add_steinerpt_to_recover_edge(startpt, endpt, NULL, 0, 1, idir)) {
                    // A Steiner point is inserted.
                    // Push this subface back to stack, to recover it again.
                    subfacstack->newindex((void **) &parysh);
                    *parysh = searchsh;
                    success = 1;
                  }
                } else if (dir == ACROSSFACE)  {
                  // to do...
                  terminate_tet_core(this, 2);
                } else {
                  terminate_tet_core(this, 2); // not possible.
                }
              } else if (ivf.iloc == NULLCAVITY) {
                // Collect a list of bad quality tets which prevent the
                //   insertion of this Steiner point.
                terminate_tet_core(this, 2);
              } else {
                terminate_tet_core(this, 2); // report a bug.
              }
            } // if (steinerflag >= 2)
            
            pointdealloc(steinerpt);
            steinerpt = NULL;
            //success = 0; // queue this subface.
          }
        } // if (dir != SELF_INTERSECT)
      } // if ((dir != SELF_INTERSECT) && steinerflag > 0)
    } // if (i == 2) else

    if (success) continue; // recover the next subface.

    if (dir == (int) SELF_INTERSECT) {
      // Found a self-intersection. This subface cannot be recovered.
      // Save it in a separate list, and remove it from the subface pool.
      if (skipped_facet_list == NULL) {
        skipped_facet_list = new arraypool(sizeof(badface), 10);
      }
      badface *bf;
      skipped_facet_list->newindex((void **) &bf);
      bf->init();
      bf->ss = searchsh;
      bf->forg  = (point) searchsh.sh[3];
      bf->fdest = (point) searchsh.sh[4];
      bf->fapex = (point) searchsh.sh[5];
      bf->key = (double) shellmark(searchsh);
      smarktest3(searchsh); // do not recover it later.
      continue; // recover the next subface.
    }

    // This subface is missing.
    if (steinerflag >= 2) {
      terminate_tet_core(this, 2);
    } // if (steinerflag >= 2)

    // Save this subface to recover it later.
    misshlist->newindex((void **) &parysh);
    *parysh = searchsh;
  } // while (subfacstack->objects > 0l)

  return 0;
}

//============================================================================//
//                                                                            //
// getvertexstar()    Return the star of a vertex.                            //
//                                                                            //
// If the flag 'fullstar' is set, return the complete star of this vertex.    //
// Otherwise, only a part of the star which is bounded by facets is returned. //
//                                                                            //
// 'tetlist' returns the list of tets in the star of the vertex 'searchpt'.   //
// Every tet in 'tetlist' is at the face opposing to 'searchpt'.              //
//                                                                            //
// 'vertlist' returns the list of vertices in the star (exclude 'searchpt').  //
//                                                                            //
// 'shlist' returns the list of subfaces in the star. Each subface must face  //
// to the interior of this star.                                              //
//                                                                            //
//============================================================================//

int TetMeshCore::getvertexstar(int fullstar, point searchpt, arraypool* tetlist, 
                              arraypool* vertlist, arraypool* shlist)
{
  triface searchtet, neightet, *parytet;
  face checksh, *parysh;
  point pt, *parypt;
  int collectflag;
  int t1ver;
  int i, j;

  point2tetorg(searchpt, searchtet);

  // Go to the opposite face (the link face) of the vertex.
  enextesymself(searchtet);
  //assert(oppo(searchtet) == searchpt);
  infect(searchtet); // Collect this tet (link face).
  tetlist->newindex((void **) &parytet);
  *parytet = searchtet;
  if (vertlist != NULL) {
    // Collect three (link) vertices.
    j = (searchtet.ver & 3); // The current vertex index.
    for (i = 1; i < 4; i++) {
      pt = (point) searchtet.tet[4 + ((j + i) % 4)];
      pinfect(pt);
      vertlist->newindex((void **) &parypt);
      *parypt = pt;
    }
  }

  collectflag = 1;
  esym(searchtet, neightet);
  if (issubface(neightet)) {
    if (shlist != NULL) {
      tspivot(neightet, checksh);
      if (!sinfected(checksh)) {
        // Collect this subface (link edge).
        sinfect(checksh);
        shlist->newindex((void **) &parysh);
        *parysh = checksh;
      }
    } 
    if (!fullstar) {
      collectflag = 0;
    }
  }
  if (collectflag) {
    fsymself(neightet); // Goto the adj tet of this face.
    esymself(neightet); // Goto the oppo face of this vertex.
    // assert(oppo(neightet) == searchpt);
    infect(neightet); // Collect this tet (link face).
    tetlist->newindex((void **) &parytet);
    *parytet = neightet;
    if (vertlist != NULL) {
      // Collect its apex.
      pt = apex(neightet);
      pinfect(pt);
      vertlist->newindex((void **) &parypt);
      *parypt = pt;
    }
  } // if (collectflag)

  // Continue to collect all tets in the star.
  for (i = 0; i < tetlist->objects; i++) {
    searchtet = * (triface *) fastlookup(tetlist, i);
    // Note that 'searchtet' is a face opposite to 'searchpt', and the neighbor
    //   tet at the current edge is already collected.
    // Check the neighbors at the other two edges of this face.
    for (j = 0; j < 2; j++) {
      collectflag = 1;
      enextself(searchtet);
      esym(searchtet, neightet);
      if (issubface(neightet)) {
        if (shlist != NULL) {
          tspivot(neightet, checksh);
          if (!sinfected(checksh)) {
            // Collect this subface (link edge).
            sinfect(checksh);
            shlist->newindex((void **) &parysh);
            *parysh = checksh;
          }
        }
        if (!fullstar) {
          collectflag = 0;
        }
      }
      if (collectflag) {
        fsymself(neightet);
        if (!infected(neightet)) {
          esymself(neightet); // Go to the face opposite to 'searchpt'.
          infect(neightet);
          tetlist->newindex((void **) &parytet);
          *parytet = neightet;
          if (vertlist != NULL) {
            // Check if a vertex is collected.
            pt = apex(neightet);
            if (!pinfected(pt)) {
              pinfect(pt);
              vertlist->newindex((void **) &parypt);
              *parypt = pt;
            }
          }
        } // if (!infected(neightet))
      } // if (collectflag)
    } // j
  } // i


  // Uninfect the list of tets and vertices.
  for (i = 0; i < tetlist->objects; i++) {
    parytet = (triface *) fastlookup(tetlist, i);
    uninfect(*parytet);
  }

  if (vertlist != NULL) {
    for (i = 0; i < vertlist->objects; i++) {
      parypt = (point *) fastlookup(vertlist, i);
      puninfect(*parypt);
    }
  }

  if (shlist != NULL) {
    for (i = 0; i < shlist->objects; i++) {
      parysh = (face *) fastlookup(shlist, i);
      suninfect(*parysh);
    }
  }

  return (int) tetlist->objects;
}

//============================================================================//
//                                                                            //
// getedge()    Get a tetrahedron having the two endpoints.                   //
//                                                                            //
// The method here is to search the second vertex in the link faces of the    //
// first vertex. The global array 'cavetetlist' is re-used for searching.     //
//                                                                            //
// This function is used for the case when the mesh is non-convex. Otherwise, //
// the function finddirection() should be faster than this.                   //
//                                                                            //
//============================================================================//

int TetMeshCore::getedge(point e1, point e2, triface *tedge)
{
  triface searchtet, neightet, *parytet;
  point pt;
  int done;
  int i, j;

  if (e1 == NULL || e2 == NULL) {
    return 0;
  }
  if ((pointtype(e1) == UNUSEDVERTEX) ||
      (pointtype(e2) == UNUSEDVERTEX)) {
    return 0;
  }

  // Quickly check if 'tedge' is just this edge.
  if (!isdeadtet(*tedge)) {
    if (org(*tedge) == e1) {
      if (dest(*tedge) == e2) {
        return 1;
      }
    } else if (org(*tedge) == e2) {
      if (dest(*tedge) == e1) {
        esymself(*tedge);
        return 1;
      }
    }
  }

  // Search for the edge [e1, e2].
  point2tetorg(e1, *tedge);
  finddirection(tedge, e2);
  if (dest(*tedge) == e2) {
    return 1;
  } else {
    // Search for the edge [e2, e1].
    point2tetorg(e2, *tedge);
    finddirection(tedge, e1);
    if (dest(*tedge) == e1) {
      esymself(*tedge);
      return 1;
    }
  }


  // Go to the link face of e1.
  point2tetorg(e1, searchtet);
  enextesymself(searchtet);
  arraypool *tetlist = cavebdrylist;

  // Search e2.
  for (i = 0; i < 3; i++) {
    pt = apex(searchtet);
    if (pt == e2) {
      // Found. 'searchtet' is [#,#,e2,e1].
      eorgoppo(searchtet, *tedge); // [e1,e2,#,#].
      return 1;
    }
    enextself(searchtet);
  }

  // Get the adjacent link face at 'searchtet'.
  fnext(searchtet, neightet);
  esymself(neightet);
  // assert(oppo(neightet) == e1);
  pt = apex(neightet);
  if (pt == e2) {
    // Found. 'neightet' is [#,#,e2,e1].
    eorgoppo(neightet, *tedge); // [e1,e2,#,#].
    return 1;
  }

  // Continue searching in the link face of e1.
  infect(searchtet);
  tetlist->newindex((void **) &parytet);
  *parytet = searchtet;
  infect(neightet);
  tetlist->newindex((void **) &parytet);
  *parytet = neightet;

  done = 0;

  for (i = 0; (i < tetlist->objects) && !done; i++) {
    parytet = (triface *) fastlookup(tetlist, i);
    searchtet = *parytet;
    for (j = 0; (j < 2) && !done; j++) {
      enextself(searchtet);
      fnext(searchtet, neightet);
      if (!infected(neightet)) {        
        esymself(neightet);
        pt = apex(neightet);
        if (pt == e2) {
          // Found. 'neightet' is [#,#,e2,e1].
          eorgoppo(neightet, *tedge);
          done = 1;
        } else {
          infect(neightet);
          tetlist->newindex((void **) &parytet);
          *parytet = neightet;
        }
      }
    } // j
  } // i 

  // Uninfect the list of visited tets.
  for (i = 0; i < tetlist->objects; i++) {
    parytet = (triface *) fastlookup(tetlist, i);
    uninfect(*parytet);
  }
  tetlist->restart();

  return done;
}

//============================================================================//
//                                                                            //
// reduceedgesatvertex()    Reduce the number of edges at a given vertex.     //
//                                                                            //
// 'endptlist' contains the endpoints of edges connecting at the vertex.      //
//                                                                            //
//============================================================================//

int TetMeshCore::reduceedgesatvertex(point startpt, arraypool* endptlist)
{
  triface searchtet;
  point *pendpt, *parypt;
  enum interresult dir;
  flipconstraints fc;
  int reduceflag;
  int count;
  int n, i, j;


  fc.remvert = startpt;
  fc.checkflipeligibility = 1;

  while (1) {

    count = 0;

    for (i = 0; i < endptlist->objects; i++) {
      pendpt = (point *) fastlookup(endptlist, i);
      if (*pendpt == dummypoint) {
        continue; // Do not reduce a virtual edge.
      }
      reduceflag = 0;
      // Find the edge.
      if (nonconvex) {
        if (getedge(startpt, *pendpt, &searchtet)) {
          dir = ACROSSVERT;
        } else {
          // The edge does not exist (was flipped).
          dir = INTERSECT;
        }
      } else {
        point2tetorg(startpt, searchtet);
        dir = finddirection(&searchtet, *pendpt);
      }
      if (dir == ACROSSVERT) {
        if (dest(searchtet) == *pendpt) {
          // Do not flip a segment.
          if (!issubseg(searchtet)) {
            n = removeedgebyflips(&searchtet, &fc);
            if (n == 2) {
              reduceflag = 1;
            }
          }
        }
        else {
          terminate_tet_core(this, 2);
        }
      } else {
        // The edge has been flipped.
        reduceflag = 1;
      }
      if (reduceflag) {
        count++;
        // Move the last vertex into this slot.
        j = endptlist->objects - 1;
        parypt = (point *) fastlookup(endptlist, j);
        *pendpt = *parypt;
        endptlist->objects--;
        i--;
      }
    } // i

    if (count == 0) {
      // No edge is reduced.
      break;
    }

  } // while (1)

  return (int) endptlist->objects;
}

//============================================================================//
//                                                                            //
// removevertexbyflips()    Remove a vertex by flips.                         //
//                                                                            //
// This routine attempts to remove the given vertex 'rempt' (p) from the      //
// tetrahedralization (T) by a sequence of flips.                             //
//                                                                            //
// The algorithm used here is a simple edge reduce method. Suppose there are  //
// n edges connected at p. We try to reduce the number of edges by flipping   //
// any edge (not a segment) that is connecting at p.                          //
//                                                                            //
// Unless T is a Delaunay tetrahedralization, there is no guarantee that 'p'  //
// can be successfully removed.                                               //
//                                                                            //
//============================================================================//

int TetMeshCore::removevertexbyflips(point steinerpt)
{
  triface *fliptets = NULL, wrktets[4];
  triface searchtet, spintet, neightet;
  face parentsh, spinsh, checksh;
  face leftseg, rightseg, checkseg;
  point lpt = NULL, rpt = NULL, apexpt; //, *parypt;
  flipconstraints fc;
  enum verttype vt;
  enum locateresult loc;
  int valence, removeflag;
  int slawson;
  int t1ver;
  int n, i;

  vt = pointtype(steinerpt);


  if (vt == FREESEGVERTEX) {
    sdecode(point2sh(steinerpt), leftseg);
    leftseg.shver = 0;
    if (sdest(leftseg) == steinerpt) {
      senext(leftseg, rightseg);
      spivotself(rightseg);
      rightseg.shver = 0;
    } else {
      rightseg = leftseg;
      senext2(rightseg, leftseg);
      spivotself(leftseg);
      leftseg.shver = 0;
    }
    lpt = sorg(leftseg);
    rpt = sdest(rightseg);
    
    // Check if both leftseg and rightseg are recovered in tet mesh.
    sstpivot1(leftseg, neightet);
    if (neightet.tet == NULL) {
      return 0; // Do not remove this Steiner point.
    }
    sstpivot1(rightseg, neightet);
    if (neightet.tet == NULL) {
      return 0; // Do not remove this Steiner point.
    }
    
    if (b->verbose > 2) {
      printf("      Removing Steiner point %d in segment (%d, %d).\n",
             pointmark(steinerpt), pointmark(lpt), pointmark(rpt));

    }
  } else if (vt == FREEFACETVERTEX) {
    if (b->verbose > 2) {
      printf("      Removing Steiner point %d in facet.\n",
             pointmark(steinerpt));
    }
  } else if (vt == FREEVOLVERTEX) {
    if (b->verbose > 2) {
      printf("      Removing Steiner point %d in volume.\n",
             pointmark(steinerpt));
    }
  } else if (vt == VOLVERTEX) {
    if (b->verbose > 2) {
      printf("      Removing a point %d in volume.\n",
             pointmark(steinerpt));
    }
  } else {
    // It is not a Steiner point.
    return 0;
  }

  // Try to reduce the number of edges at 'p' by flips.
  getvertexstar(1, steinerpt, cavetetlist, cavetetvertlist, NULL);
  cavetetlist->restart(); // This list may be re-used.
  if (cavetetvertlist->objects > 3l) {
    valence = reduceedgesatvertex(steinerpt, cavetetvertlist);
  } else {
    valence = cavetetvertlist->objects;
  }
  cavetetvertlist->restart();

  removeflag = 0;

  if (valence == 4) {
    // Only 4 vertices (4 tets) left! 'p' is inside the convex hull of the 4
    //   vertices. This case is due to that 'p' is not exactly on the segment.
    point2tetorg(steinerpt, searchtet);
    loc = INTETRAHEDRON;
    removeflag = 1;
  } else if (valence == 5) {
    // There are 5 edges.
    if (vt == FREESEGVERTEX) {
      sstpivot1(leftseg, searchtet);
      if (org(searchtet) != steinerpt) {
        esymself(searchtet);
      }
      i = 0; // Count the numbe of tet at the edge [p,lpt].
      neightet.tet = NULL; // Init the face.
      spintet = searchtet;
      while (1) {
        i++;
        if (apex(spintet) == rpt) {
          // Remember the face containing the edge [lpt, rpt].
          neightet = spintet;
        }
        fnextself(spintet);
        if (spintet.tet == searchtet.tet) break;
      }
      if (i == 3) {
        // This case has been checked below.
      } else if (i == 4) {
        // There are 4 tets sharing at [p,lpt]. There must be 4 tets sharing
        //   at [p,rpt].  There must be a face [p, lpt, rpt].  
        if (apex(neightet) == rpt) {
          // The edge (segment) has been already recovered!  
          // Check if a 6-to-2 flip is possible (to remove 'p').
          // Let 'searchtet' be [p,d,a,b]
          esym(neightet, searchtet);
          enextself(searchtet);
          // Check if there are exactly three tets at edge [p,d].
          wrktets[0] = searchtet; // [p,d,a,b]
          for (i = 0; i < 2; i++) {
            fnext(wrktets[i], wrktets[i+1]); // [p,d,b,c], [p,d,c,a]
          }
          if (apex(wrktets[0]) == oppo(wrktets[2])) {
            loc = ONFACE;
            removeflag = 1;
          }
        }
      }
    } else if (vt == FREEFACETVERTEX) {
      // It is possible to do a 6-to-2 flip to remove the vertex.
      point2tetorg(steinerpt, searchtet);
      // Get the three faces of 'searchtet' which share at p.
      //    All faces has p as origin.
      wrktets[0] = searchtet;
      wrktets[1] = searchtet;
      esymself(wrktets[1]);
      enextself(wrktets[1]);
      wrktets[2] = searchtet;
      eprevself(wrktets[2]);
      esymself(wrktets[2]);
      // All internal edges of the six tets have valance either 3 or 4.
      // Get one edge which has valance 3.
      searchtet.tet = NULL;
      for (i = 0; i < 3; i++) {
        spintet = wrktets[i];
        valence = 0;
        while (1) {
          valence++;
          fnextself(spintet);
          if (spintet.tet == wrktets[i].tet) break;
        }
        if (valence == 3) {
          // Found the edge.
          searchtet = wrktets[i];
          break;
        }
      }
      // Note, we do not detach the three subfaces at p.
      // They will be removed within a 4-to-1 flip.
      loc = ONFACE;
      removeflag = 1;
    }
    //removeflag = 1;
  } 

  if (!removeflag) {
    if (vt == FREESEGVERTEX) { 
      // Check is it possible to recover the edge [lpt,rpt].
      // The condition to check is:  Whether each tet containing 'leftseg' is
      //   adjacent to a tet containing 'rightseg'.
      sstpivot1(leftseg, searchtet);
      if (org(searchtet) != steinerpt) {
        esymself(searchtet);
      }
      spintet = searchtet;
      while (1) {
        // Go to the bottom face of this tet.
        eprev(spintet, neightet);
        esymself(neightet);  // [steinerpt, p1, p2, lpt]
        // Get the adjacent tet.
        fsymself(neightet);  // [p1, steinerpt, p2, rpt]
        if (oppo(neightet) != rpt) {
          // Found a non-matching adjacent tet.
          break;
        }
		{
			// [2017-10-15] Check if the tet is inverted?
			point chkp1 = org(neightet);
			point chkp2 = apex(neightet);
			double chkori = orient3d(rpt, lpt, chkp1, chkp2);
			if (chkori >= 0.0) {
				// Either inverted or degenerated.
				break;
			}
		}
        fnextself(spintet);
        if (spintet.tet == searchtet.tet) {
          // 'searchtet' is [p,d,p1,p2].
          loc = ONEDGE;
          removeflag = 1;
          break;
        }
      }
    } // if (vt == FREESEGVERTEX)
  }

  if (!removeflag) {
    if (vt == FREESEGVERTEX) {
      // Check if the edge [lpt, rpt] exists.
      if (getedge(lpt, rpt, &searchtet)) {
        // We have recovered this edge. Shift the vertex into the volume.
        // We can recover this edge if the subfaces are not recovered yet.
        if (!checksubfaceflag) {
          // Remove the vertex from the surface mesh.
          //   This will re-create the segment [lpt, rpt] and re-triangulate
          //   all the facets at the segment.
          // Detach the subsegments from their surrounding tets.
          for (i = 0; i < 2; i++) {
            checkseg = (i == 0) ? leftseg : rightseg;
            sstpivot1(checkseg, neightet);
            spintet = neightet;
            while (1) {
              tssdissolve1(spintet);
              fnextself(spintet);
              if (spintet.tet == neightet.tet) break;
            }
            sstdissolve1(checkseg);
          } // i
          slawson = 1; // Do lawson flip after removal.
          spivot(rightseg, parentsh); // 'rightseg' has p as its origin.
          sremovevertex(steinerpt, &parentsh, &rightseg, slawson);
          // Clear the list for new subfaces.
          caveshbdlist->restart();
          // Insert the new segment.
          sstbond1(rightseg, searchtet);
          spintet = searchtet;
          while (1) {
            tssbond1(spintet, rightseg);
            fnextself(spintet);
            if (spintet.tet == searchtet.tet) break;
          }
          // The Steiner point has been shifted into the volume.
          setpointtype(steinerpt, FREEVOLVERTEX);          
          st_segref_count--;
          st_volref_count++;
          return 1;
        } // if (!checksubfaceflag)
      } // if (getedge(...))
    } // if (vt == FREESEGVERTEX)
  } // if (!removeflag)

  if (!removeflag) {
    return 0;
  }

  if (vt == FREESEGVERTEX) {
    // Detach the subsegments from their surronding tets.
    for (i = 0; i < 2; i++) {
      checkseg = (i == 0) ? leftseg : rightseg;
      sstpivot1(checkseg, neightet);
      spintet = neightet;
      while (1) {
        tssdissolve1(spintet);
        fnextself(spintet);
        if (spintet.tet == neightet.tet) break;
      }
      sstdissolve1(checkseg);
    } // i
    if (checksubfaceflag) {
      // Detach the subfaces at the subsegments from their attached tets.
      for (i = 0; i < 2; i++) {
        checkseg = (i == 0) ? leftseg : rightseg;
        spivot(checkseg, parentsh);
        if (parentsh.sh != NULL) {
          spinsh = parentsh;
          while (1) {
            stpivot(spinsh, neightet);
            if (neightet.tet != NULL) {
              tsdissolve(neightet);
            }
            sesymself(spinsh);
            stpivot(spinsh, neightet);
            if (neightet.tet != NULL) {
              tsdissolve(neightet);
            }
            stdissolve(spinsh);
            spivotself(spinsh); // Go to the next subface.
            if (spinsh.sh == parentsh.sh) break;
          }
        }
      } // i
    } // if (checksubfaceflag)
  }

  if (loc == INTETRAHEDRON) {
    // Collect the four tets containing 'p'.
    fliptets = new triface[4];
    fliptets[0] = searchtet; // [p,d,a,b]
    for (i = 0; i < 2; i++) {
      fnext(fliptets[i], fliptets[i+1]); // [p,d,b,c], [p,d,c,a]
    }
    eprev(fliptets[0], fliptets[3]);
    fnextself(fliptets[3]); // it is [a,p,b,c]
    eprevself(fliptets[3]);
    esymself(fliptets[3]); // [a,b,c,p].
    if (vt == FREEFACETVERTEX) {
      // [2018-03-08] Check if the last 4-to-1 flip is valid.
      // fliptets[0],[1],[2] are [p,d,a,b],[p,d,b,c],[p,d,c,a]
      triface checktet, chkface;
      for (i = 0; i < 3; i++) {
        enext(fliptets[i], checktet);
        esymself(checktet); // [a,d,b,p],[b,d,c,p],[c,d,a,p]
        int scount = 0; int k;
        for (k = 0; k < 3; k++) {
          esym(checktet, chkface);
          if (issubface(chkface)) scount++;
          enextself(checktet);
        }
        if (scount == 3) {
          break; // Found a tet which support a 3-to-1 flip.
        } else if (scount == 2) {
          // This is a strange configuration. Debug it.
          // Do not do this flip.
          delete [] fliptets;
          return 0;
        }
      }
      if (i == 3) {
        // No tet in [p,d,a,b],[p,d,b,c],[p,d,c,a] support it.
        int scount = 0;
        for (i = 0; i < 3; i++) {
          eprev(fliptets[i], checktet);
          esymself(checktet); // [p,a,b,d],[p,b,c,d],[p,c,a,d]
          if (issubface(chkface)) scount++;
        }
        if (scount != 3) {
          // Do not do this flip.
          delete [] fliptets;
          return 0;
        }
      }
    } // if (vt == FREEFACETVERTEX)
    flip41(fliptets, 1, &fc);
    //recenttet = fliptets[0];
  } else if (loc == ONFACE) {
    // Let the original two tets be [a,b,c,d] and [b,a,c,e]. And p is in
    //   face [a,b,c].  Let 'searchtet' be the tet [p,d,a,b].
    // Collect the six tets containing 'p'.
    fliptets = new triface[6];
    fliptets[0] = searchtet; // [p,d,a,b]
    for (i = 0; i < 2; i++) {
      fnext(fliptets[i], fliptets[i+1]); // [p,d,b,c], [p,d,c,a]
    }
    eprev(fliptets[0], fliptets[3]);
    fnextself(fliptets[3]); // [a,p,b,e]
    esymself(fliptets[3]);  // [p,a,e,b]
    eprevself(fliptets[3]); // [e,p,a,b]
    for (i = 3; i < 5; i++) {
      fnext(fliptets[i], fliptets[i+1]); // [e,p,b,c], [e,p,c,a]
    }
    if (vt == FREEFACETVERTEX) {
      // We need to determine the location of three subfaces at p.
      valence = 0; // Re-use it.
      for (i = 3; i < 6; i++) {
        if (issubface(fliptets[i])) valence++;
      }
      if (valence > 0) {
        // We must do 3-to-2 flip in the upper part. We simply re-arrange
        //   the six tets.
        for (i = 0; i < 3; i++) {
          esym(fliptets[i+3], wrktets[i]);
          esym(fliptets[i], fliptets[i+3]);
          fliptets[i] = wrktets[i];
        }
        // Swap the last two pairs, i.e., [1]<->[[2], and [4]<->[5]
        wrktets[1] = fliptets[1];
        fliptets[1] = fliptets[2];
        fliptets[2] = wrktets[1];
        wrktets[1] = fliptets[4];
        fliptets[4] = fliptets[5];
        fliptets[5] = wrktets[1];
      }
      // [2018-03-08] Check if the last 4-to-1 flip is valid.
      // fliptets[0],[1],[2] are [p,d,a,b],[p,d,b,c],[p,d,c,a]
      triface checktet, chkface;
      for (i = 0; i < 3; i++) {
        enext(fliptets[i], checktet);
        esymself(checktet); // [a,d,b,p],[b,d,c,p],[c,d,a,p]
        int scount = 0; int k;
        for (k = 0; k < 3; k++) {
          esym(checktet, chkface);
          if (issubface(chkface)) scount++;
          enextself(checktet);
        }
        if (scount == 3) {
          break; // Found a tet which support a 3-to-1 flip.
        } else if (scount == 2) {
          // This is a strange configuration. Debug it.
          // Do not do this flip.
          delete [] fliptets;
          return 0;
        }
      }
      if (i == 3) {
        // No tet in [p,d,a,b],[p,d,b,c],[p,d,c,a] support it.
        int scount = 0;
        for (i = 0; i < 3; i++) {
          eprev(fliptets[i], checktet);
          esymself(checktet); // [p,a,b,d],[p,b,c,d],[p,c,a,d]
          if (issubface(chkface)) scount++;
        }
        if (scount != 3) {
          // Do not do this flip.
          delete [] fliptets;
          return 0;
        }
      }
    } // vt == FREEFACETVERTEX
    // Remove p by a 6-to-2 flip, which is a combination of two flips:
    //   a 3-to-2 (deletes the edge [e,p]), and
    //   a 4-to-1 (deletes the vertex p).
    // First do a 3-to-2 flip on [e,p,a,b],[e,p,b,c],[e,p,c,a]. It creates
    //   two new tets: [a,b,c,p] and [b,a,c,e].  The new tet [a,b,c,p] is
    //   degenerate (has zero volume). It will be deleted in the followed
    //   4-to-1 flip.
    //flip32(&(fliptets[3]), 1, 0, 0);
    flip32(&(fliptets[3]), 1, &fc);
    // Second do a 4-to-1 flip on [p,d,a,b],[p,d,b,c],[p,d,c,a],[a,b,c,p].
    //   This creates a new tet [a,b,c,d].
    //flip41(fliptets, 1, 0, 0);
    flip41(fliptets, 1, &fc);
    //recenttet = fliptets[0];
  } else if (loc == ONEDGE) {
    // Let the original edge be [e,d] and p is in [e,d]. Assume there are n
    //   tets sharing at edge [e,d] originally.  We number the link vertices
    //   of [e,d]: p_0, p_1, ..., p_n-1. 'searchtet' is [p,d,p_0,p_1].
    // Count the number of tets at edge [e,p] and [p,d] (this is n).
    n = 0;
    spintet = searchtet;
    while (1) {
      n++;
      fnextself(spintet);
      if (spintet.tet == searchtet.tet) break;
    }
    // Collect the 2n tets containing 'p'.
    fliptets = new triface[2 * n];
    fliptets[0] = searchtet; // [p,b,p_0,p_1] 
    for (i = 0; i < (n - 1); i++) {
      fnext(fliptets[i], fliptets[i+1]); // [p,d,p_i,p_i+1].
    }
    eprev(fliptets[0], fliptets[n]);
    fnextself(fliptets[n]); // [p_0,p,p_1,e]
    esymself(fliptets[n]);  // [p,p_0,e,p_1]
    eprevself(fliptets[n]); // [e,p,p_0,p_1]
    for (i = n; i <  (2 * n - 1); i++) {
      fnext(fliptets[i], fliptets[i+1]); // [e,p,p_i,p_i+1].
    }
    // Remove p by a 2n-to-n flip, it is a sequence of n flips:
    // - Do a 2-to-3 flip on 
    //     [p_0,p_1,p,d] and 
    //     [p,p_1,p_0,e].
    //   This produces: 
    //     [e,d,p_0,p_1], 
    //     [e,d,p_1,p] (degenerated), and 
    //     [e,d,p,p_0] (degenerated).
    wrktets[0] = fliptets[0]; // [p,d,p_0,p_1]
    eprevself(wrktets[0]);    // [p_0,p,d,p_1]
    esymself(wrktets[0]);     // [p,p_0,p_1,d]
    enextself(wrktets[0]);    // [p_0,p_1,p,d] [0]
    wrktets[1] = fliptets[n]; // [e,p,p_0,p_1]
    enextself(wrktets[1]);    // [p,p_0,e,p_1]
    esymself(wrktets[1]);     // [p_0,p,p_1,e]
    eprevself(wrktets[1]);    // [p_1,p_0,p,e] [1]
    //flip23(wrktets, 1, 0, 0);
    flip23(wrktets, 1, &fc);
    // Save the new tet [e,d,p,p_0] (degenerated).
    fliptets[n] = wrktets[2];
    // Save the new tet [e,d,p_0,p_1].
    fliptets[0] = wrktets[0];
    // - Repeat from i = 1 to n-2: (n - 2) flips
    //   - Do a 3-to-2 flip on 
    //       [p,p_i,d,e], 
    //       [p,p_i,e,p_i+1], and 
    //       [p,p_i,p_i+1,d]. 
    //     This produces: 
    //       [d,e,p_i+1,p_i], and
    //       [e,d,p_i+1,p] (degenerated).
    for (i = 1; i < (n - 1); i++) {
      wrktets[0] = wrktets[1]; // [e,d,p_i,p] (degenerated).
      enextself(wrktets[0]);   // [d,p_i,e,p] (...)
      esymself(wrktets[0]);    // [p_i,d,p,e] (...) 
      eprevself(wrktets[0]);   // [p,p_i,d,e] (degenerated) [0].
      wrktets[1] = fliptets[n+i];  // [e,p,p_i,p_i+1]
      enextself(wrktets[1]);       // [p,p_i,e,p_i+1] [1]
      wrktets[2] = fliptets[i]; // [p,d,p_i,p_i+1]
      eprevself(wrktets[2]);    // [p_i,p,d,p_i+1]
      esymself(wrktets[2]);     // [p,p_i,p_i+1,d] [2]
      //flip32(wrktets, 1, 0, 0);
      flip32(wrktets, 1, &fc);
      // Save the new tet [e,d,p_i,p_i+1].         // FOR DEBUG ONLY
      fliptets[i] = wrktets[0]; // [d,e,p_i+1,p_i] // FOR DEBUG ONLY
      esymself(fliptets[i]);    // [e,d,p_i,p_i+1] // FOR DEBUG ONLY
    }
    // - Do a 4-to-1 flip on 
    //     [p,p_0,e,d],     [d,e,p_0,p],
    //     [p,p_0,d,p_n-1], [e,p_n-1,p_0,p], 
    //     [p,p_0,p_n-1,e], [p_0,p_n-1,d,p], and
    //     [e,d,p_n-1,p]. 
    //   This produces 
    //     [e,d,p_n-1,p_0] and 
    //     deletes p.
    wrktets[3] = wrktets[1];  // [e,d,p_n-1,p] (degenerated) [3]
    wrktets[0] = fliptets[n]; // [e,d,p,p_0] (degenerated)
    eprevself(wrktets[0]);    // [p,e,d,p_0] (...)
    esymself(wrktets[0]);     // [e,p,p_0,d] (...)
    enextself(wrktets[0]);    // [p,p_0,e,d] (degenerated) [0]
    wrktets[1] = fliptets[n-1];   // [p,d,p_n-1,p_0]
    esymself(wrktets[1]);         // [d,p,p_0,p_n-1]
    enextself(wrktets[1]);        // [p,p_0,d,p_n-1] [1]
    wrktets[2] = fliptets[2*n-1]; // [e,p,p_n-1,p_0]
    enextself(wrktets[2]);        // [p_p_n-1,e,p_0]
    esymself(wrktets[2]);         // [p_n-1,p,p_0,e]
    enextself(wrktets[2]);        // [p,p_0,p_n-1,e] [2]
    //flip41(wrktets, 1, 0, 0);
    flip41(wrktets, 1, &fc);
    // Save the new tet [e,d,p_n-1,p_0]             // FOR DEBUG ONLY
    fliptets[n-1] = wrktets[0];  // [e,d,p_n-1,p_0] // FOR DEBUG ONLY
    //recenttet = fliptets[0];
  }

  delete [] fliptets;

  if (vt == FREESEGVERTEX) {
    // Remove the vertex from the surface mesh.
    //   This will re-create the segment [lpt, rpt] and re-triangulate
    //   all the facets at the segment.
    // Only do lawson flip when subfaces are not recovery yet.
    slawson = (checksubfaceflag ? 0 : 1);
    spivot(rightseg, parentsh); // 'rightseg' has p as its origin.
    sremovevertex(steinerpt, &parentsh, &rightseg, slawson);

    // The original segment is returned in 'rightseg'. 
    rightseg.shver = 0;
    // Insert the new segment.
    point2tetorg(lpt, searchtet);
    finddirection(&searchtet, rpt);
    if (dest(searchtet) != rpt) {
      terminate_tet_core(this, 2);
    }
    sstbond1(rightseg, searchtet);
    spintet = searchtet;
    while (1) {
      tssbond1(spintet, rightseg);
      fnextself(spintet);
      if (spintet.tet == searchtet.tet) break;
    }

    if (checksubfaceflag) {
      // Insert subfaces at segment [lpt,rpt] into the tetrahedralization.
      spivot(rightseg, parentsh);
      if (parentsh.sh != NULL) {
        spinsh = parentsh;
        while (1) {
          if (sorg(spinsh) != lpt) {
            sesymself(spinsh);
          }
          apexpt = sapex(spinsh);
          // Find the adjacent tet of [lpt,rpt,apexpt];
          spintet = searchtet;
          while (1) {
            if (apex(spintet) == apexpt) {
              tsbond(spintet, spinsh);
              sesymself(spinsh); // Get to another side of this face.
              fsym(spintet, neightet);
              tsbond(neightet, spinsh);
              sesymself(spinsh); // Get back to the original side.
              break;
            }
            fnextself(spintet);
          }
          spivotself(spinsh);
          if (spinsh.sh == parentsh.sh) break;
        }
      }
    } // if (checksubfaceflag)

    // Clear the set of new subfaces.
    caveshbdlist->restart();
  } // if (vt == FREESEGVERTEX)

  // The point has been removed.
  if (pointtype(steinerpt) != UNUSEDVERTEX) {
    setpointtype(steinerpt, UNUSEDVERTEX);
    unuverts++;
  }
  if (vt != VOLVERTEX) {
    // Update the correspinding counters.
    if (vt == FREESEGVERTEX) {
      st_segref_count--;
    } else if (vt == FREEFACETVERTEX) {
      st_facref_count--;
    } else if (vt == FREEVOLVERTEX) {
      st_volref_count--;
    }
    if (steinerleft > 0) steinerleft++;
  }

  return 1;
}

//============================================================================//
//                                                                            //
// smoothpoint()    Moving a vertex to improve the mesh quality.              //
//                                                                            //
// 'smtpt' (p) is a point to be smoothed. Generally, it is a Steiner point.   //
// It may be not a vertex of the mesh.                                        //
//                                                                            //
// This routine tries to move 'p' inside its star until a selected objective  //
// function over all tetrahedra in the star is improved. The function may be  //
// the some quality measures, i.e., aspect ratio, maximum dihedral angel, or  //
// simply the volume of the tetrahedra.                                       //
//                                                                            //
// 'linkfacelist' contains the list of link faces of 'p'.  Since a link face  //
// has two orientations, ccw or cw, with respect to 'p'.  'ccw' indicates     //
// the orientation is ccw (1) or not (0).                                     //
//                                                                            //
// 'opm' is a structure contains the parameters of the objective function.    //
// It is needed by the evaluation of the function value.                      //
//                                                                            //
// The return value indicates weather the point is smoothed or not.           //
//                                                                            //
// ASSUMPTION: This routine assumes that all link faces are true faces, i.e,  //
// no face has 'dummypoint' as its vertex.                                    //
//                                                                            //
//============================================================================//

int TetMeshCore::smoothpoint(point smtpt, arraypool *linkfacelist, int ccw,
                            optparameters *opm)
{
  triface *parytet, *parytet1, swaptet;
  badface bf;
  point pa, pb, pc;
  double fcent[3], startpt[3], nextpt[3], bestpt[3];
  double oldval, minval = 0.0, val;
  double maxcosd; // oldang, newang;
  double ori, diff;
  int numdirs, iter;
  int i, j, k;

  // Decide the number of moving directions.
  numdirs = (int) linkfacelist->objects;
  if (numdirs > opm->numofsearchdirs) {
    numdirs = opm->numofsearchdirs; // Maximum search directions.
  }

  // Set the initial value.
  opm->imprval = opm->initval;
  iter = 0;

  for (i = 0; i < 3; i++) {
    bestpt[i] = startpt[i] = smtpt[i];
  }

  // Iterate until the obj function is not improved.
  while (1) {

    // Find the best next location.
    oldval = opm->imprval;

    for (i = 0; i < numdirs; i++) {
      // Randomly pick a link face (0 <= k <= objects - i - 1).
      k = (int) randomnation(linkfacelist->objects - i);
      parytet = (triface *) fastlookup(linkfacelist, k);
      // Calculate a new position from 'p' to the center of this face.
      pa = org(*parytet);
      pb = dest(*parytet);
      pc = apex(*parytet);
      for (j = 0; j < 3; j++) {
        fcent[j] = (pa[j] + pb[j] + pc[j]) / 3.0;
      }
      for (j = 0; j < 3; j++) {
        nextpt[j] = startpt[j] + opm->searchstep * (fcent[j] - startpt[j]);
      }
      // Calculate the largest minimum function value for the new location.
      for (j = 0; j < linkfacelist->objects; j++) {
        parytet = (triface *) fastlookup(linkfacelist, j);
        if (ccw) {
          pa = org(*parytet);
          pb = dest(*parytet);
        } else {
          pb = org(*parytet);
          pa = dest(*parytet);
        }
        pc = apex(*parytet);
        ori = orient3d(pa, pb, pc, nextpt);
        if (ori < 0.0) {
          // Calcuate the objective function value.
          if (opm->max_min_volume) {
            //val = -ori;
            val = - orient3dfast(pa, pb, pc, nextpt);
          } else if (opm->min_max_aspectratio) {
            get_tetqual(pa, pb, pc, nextpt, &bf);
            val = 1.0 / bf.key;
          } else if (opm->min_max_dihedangle) {
            get_tetqual(pa, pb, pc, nextpt, &bf);
            maxcosd = bf.cent[0];
            if (maxcosd < -1) maxcosd = -1.0; // Rounding.
            val = maxcosd + 1.0; // Make it be positive.
          } else {
            // Unknown objective function.
            val = 0.0;
          }
        } else { // ori >= 0.0;
          // An invalid new tet.
          // This may happen if the mesh contains inverted elements.
          if (opm->max_min_volume) {
            //val = -ori;
            val = - orient3dfast(pa, pb, pc, nextpt);
          } else {
            // Discard this point.
            break; // j
          }
        } // if (ori >= 0.0)
        // Stop looping when the object value is not improved.
        if (val <= opm->imprval) {
          break; // j
        } else {
          // Remember the smallest improved value.
          if (j == 0) {
            minval = val;
          } else {
            minval = (val < minval) ? val : minval;
          }
        }
      } // j
      if (j == linkfacelist->objects) {
        // The function value has been improved.
        opm->imprval = minval;
        // Save the new location of the point.
        for (j = 0; j < 3; j++) bestpt[j] = nextpt[j];
      }
      // Swap k-th and (object-i-1)-th entries.
      j = linkfacelist->objects - i - 1;
      parytet  = (triface *) fastlookup(linkfacelist, k);
      parytet1 = (triface *) fastlookup(linkfacelist, j);
      swaptet = *parytet1;
      *parytet1 = *parytet;
      *parytet = swaptet;
    } // i

    diff = opm->imprval - oldval;
    if (diff > 0.0) {
      // Is the function value improved effectively?
      if (opm->max_min_volume) {
        //if ((diff / oldval) < b->epsilon) diff = 0.0;
      } else if (opm->min_max_aspectratio) {
        if ((diff / oldval) < 1e-3) diff = 0.0;
      } else if (opm->min_max_dihedangle) {
        //oldang = acos(oldval - 1.0);
        //newang = acos(opm->imprval - 1.0);
        //if ((oldang - newang) < 0.00174) diff = 0.0; // about 0.1 degree.
      } else {
        // Unknown objective function.
        terminate_tet_core(this, 2);
      }
    }

    if (diff > 0.0) {
      // Yes, move p to the new location and continue.
      for (j = 0; j < 3; j++) startpt[j] = bestpt[j];
      iter++;
      if ((opm->maxiter > 0) && (iter >= opm->maxiter)) {
        // Maximum smoothing iterations reached.
        break;
      }
    } else {
      break;
    }

  } // while (1)

  if (iter > 0) {
    // The point has been smoothed.
    opm->smthiter = iter; // Remember the number of iterations.
    // The point has been smoothed. Update it to its new position.
    for (i = 0; i < 3; i++) smtpt[i] = startpt[i];
  }

  return iter;
}

//============================================================================//
//                                                                            //
// suppressbdrysteinerpoint()    Suppress a boundary Steiner point            //
//                                                                            //
//============================================================================//

int TetMeshCore::suppressbdrysteinerpoint(point steinerpt)
{
  face parentsh, spinsh, *parysh;
  face leftseg, rightseg;
  point lpt = NULL, rpt = NULL;
  int i;

  verttype vt = pointtype(steinerpt);

  if (vt == FREESEGVERTEX) {
    sdecode(point2sh(steinerpt), leftseg);
    leftseg.shver = 0;
    if (sdest(leftseg) == steinerpt) {
      senext(leftseg, rightseg);
      spivotself(rightseg);
      rightseg.shver = 0;
    } else {
      rightseg = leftseg;
      senext2(rightseg, leftseg);
      spivotself(leftseg);
      leftseg.shver = 0;
    }
    lpt = sorg(leftseg);
    rpt = sdest(rightseg);
    if (b->verbose > 2) {
      printf("      Suppressing Steiner point %d in segment (%d, %d).\n",
             pointmark(steinerpt), pointmark(lpt), pointmark(rpt));
    }
    // Get all subfaces at the left segment [lpt, steinerpt].
    spivot(leftseg, parentsh);
    if (parentsh.sh != NULL) {
      // It is not a dangling segment.
      spinsh = parentsh;
      while (1) {
        cavesegshlist->newindex((void **) &parysh);
        *parysh = spinsh;
        // Orient the face consistently. 
        if (sorg(*parysh)!= sorg(parentsh)) sesymself(*parysh);
        spivotself(spinsh);
        if (spinsh.sh == NULL) break;
        if (spinsh.sh == parentsh.sh) break;
      }
    }
    if (cavesegshlist->objects < 2) {
      // It is a single segment. Not handle it yet.
      cavesegshlist->restart();
      return 0;
    }
  } else if (vt == FREEFACETVERTEX) {
    if (b->verbose > 2) {
      printf("      Suppressing Steiner point %d from facet.\n",
             pointmark(steinerpt));
    }
    sdecode(point2sh(steinerpt), parentsh);
    // A facet Steiner point. There are exactly two sectors.
    for (i = 0; i < 2; i++) {
      cavesegshlist->newindex((void **) &parysh);
      *parysh = parentsh;
      sesymself(parentsh);
    }
  } else {
    return 0; // no need to suppress it.
  }

  triface searchtet, neightet, *parytet;
  point pa, pb, pc, pd;
  double v1[3], v2[3], len, u;

  double startpt[3] = {0,}, samplept[3] = {0,}, candpt[3] = {0,};
  double ori, minvol, smallvol;
  int samplesize;
  int it, j, k;

  int n = (int) cavesegshlist->objects;
  point *newsteiners = new point[n];
  for (i = 0; i < n; i++) newsteiners[i] = NULL;

  // Search for each sector an interior vertex. 
  for (i = 0; i < cavesegshlist->objects; i++) {
    parysh = (face *) fastlookup(cavesegshlist, i);
    stpivot(*parysh, searchtet);
    // Skip it if it is outside.
    if (ishulltet(searchtet)) continue;
    // Get the "half-ball". Tets in 'cavetetlist' all contain 'steinerpt' as
    //   opposite.  Subfaces in 'caveshlist' all contain 'steinerpt' as apex.
    //   Moreover, subfaces are oriented towards the interior of the ball.
    setpoint2tet(steinerpt, encode(searchtet));
    getvertexstar(0, steinerpt, cavetetlist, NULL, caveshlist);
    // Calculate the searching vector.
    pa = sorg(*parysh);
    pb = sdest(*parysh);
    pc = sapex(*parysh);
    facenormal(pa, pb, pc, v1, 1, NULL);
    len = sqrt(dot(v1, v1));
    v1[0] /= len;
    v1[1] /= len;
    v1[2] /= len;
    if (vt == FREESEGVERTEX) {
      parysh = (face *) fastlookup(cavesegshlist, (i + 1) % n);
      pd = sapex(*parysh);
      facenormal(pb, pa, pd, v2, 1, NULL);
      len = sqrt(dot(v2, v2));
      v2[0] /= len;
      v2[1] /= len;
      v2[2] /= len;
      // Average the two vectors.
      v1[0] = 0.5 * (v1[0] + v2[0]);
      v1[1] = 0.5 * (v1[1] + v2[1]);
      v1[2] = 0.5 * (v1[2] + v2[2]);
    }
    // Search the intersection of the ray starting from 'steinerpt' to
    //   the search direction 'v1' and the shell of the half-ball.
    // - Construct an endpoint.
    len = distance(pa, pb);
    v2[0] = steinerpt[0] + len * v1[0];
    v2[1] = steinerpt[1] + len * v1[1];
    v2[2] = steinerpt[2] + len * v1[2];
    for (j = 0; j < cavetetlist->objects; j++) {
      parytet = (triface *) fastlookup(cavetetlist, j);
      pa = org(*parytet);
      pb = dest(*parytet);
      pc = apex(*parytet);
      // Test if the ray startpt->v2 lies in the cone: where 'steinerpt'
      //   is the apex, and three sides are defined by the triangle 
      //   [pa, pb, pc].
      ori = orient3d(steinerpt, pa, pb, v2);
      if (ori >= 0) {
        ori = orient3d(steinerpt, pb, pc, v2);
        if (ori >= 0) {
          ori = orient3d(steinerpt, pc, pa, v2);
          if (ori >= 0) {
            // Found! Calculate the intersection.
            planelineint(pa, pb, pc, steinerpt, v2, startpt, &u);
            break;
          }
        }
      }
    } // j
    if (j == cavetetlist->objects) {
      break; // There is no intersection!! Debug is needed.
    }
    // Close the ball by adding the subfaces.
    for (j = 0; j < caveshlist->objects; j++) {
      parysh = (face *) fastlookup(caveshlist, j);
      stpivot(*parysh, neightet);
      cavetetlist->newindex((void **) &parytet);
      *parytet = neightet;
    }
    // Search a best point inside the segment [startpt, steinerpt].
    it = 0;
    samplesize = 100;
    v1[0] = steinerpt[0] - startpt[0];
    v1[1] = steinerpt[1] - startpt[1];
    v1[2] = steinerpt[2] - startpt[2];
    minvol = -1.0;
    while (it < 3) {
      for (j = 1; j < samplesize - 1; j++) {
        samplept[0] = startpt[0] + ((double) j / (double) samplesize) * v1[0];
        samplept[1] = startpt[1] + ((double) j / (double) samplesize) * v1[1];
        samplept[2] = startpt[2] + ((double) j / (double) samplesize) * v1[2];
        // Find the minimum volume for 'samplept'.
        smallvol = -1;
        for (k = 0; k < cavetetlist->objects; k++) {
          parytet = (triface *) fastlookup(cavetetlist, k);
          pa = org(*parytet);
          pb = dest(*parytet);
          pc = apex(*parytet);
          ori = orient3d(pb, pa, pc, samplept);
          {
            // [2017-10-15] Rounding
            double lab = distance(pa, pb);
            double lbc = distance(pb, pc);
            double lca = distance(pc, pa);
            double lv = (lab + lbc + lca) / 3.0;
            double l3 = lv*lv*lv;
            if (fabs(ori) / l3 < 1e-8) ori = 0.0;
          }
          if (ori <= 0) {
            break; // An invalid tet.
          }
          if (smallvol == -1) {
            smallvol = ori;
          } else {
            if (ori < smallvol) smallvol = ori;
          }
        } // k
        if (k == cavetetlist->objects) {
          // Found a valid point. Remember it.
          if (minvol == -1.0) {
            candpt[0] = samplept[0];
            candpt[1] = samplept[1];
            candpt[2] = samplept[2];
            minvol = smallvol;
          } else {
            if (minvol < smallvol) {
              // It is a better location. Remember it.
              candpt[0] = samplept[0];
              candpt[1] = samplept[1];
              candpt[2] = samplept[2];
              minvol = smallvol;
            } else {
              // No improvement of smallest volume. 
              // Since we are searching along the line [startpt, steinerpy],
              // The smallest volume can only be decreased later.
              break;
            }
          }
        }
      } // j
      if (minvol > 0) break; 
      samplesize *= 10;
      it++;
    } // while (it < 3)
    if (minvol == -1.0) {
      // Failed to find a valid point.
      cavetetlist->restart();
      caveshlist->restart();
      break;
    }
    // Create a new Steiner point inside this section.
    makepoint(&(newsteiners[i]), FREEVOLVERTEX);
    newsteiners[i][0] = candpt[0];
    newsteiners[i][1] = candpt[1];
    newsteiners[i][2] = candpt[2];
    cavetetlist->restart();
    caveshlist->restart();
  } // i

  if (i < cavesegshlist->objects) {
    // Failed to suppress the vertex.
    for (; i > 0; i--) {
      if (newsteiners[i - 1] != NULL) {
        pointdealloc(newsteiners[i - 1]);
      }
    }
    delete [] newsteiners;
    cavesegshlist->restart();
    return 0;
  }

  // First insert Steiner points into the mesh.    
  // 'cavesegshlist' will be used by insertpoint().
  //int nfaces = cavesegshlist->objects;
  face *segshlist = new face[n];
  for (i = 0; i < cavesegshlist->objects; i++) {
    segshlist[i] = * (face *) fastlookup(cavesegshlist, i);    
  }
  cavesegshlist->restart();

  for (i = 0; i < n; i++) {    
    //assert(caveoldtetlist->objects == 0); 
    //assert(cavetetlist->objects == 0);
    parysh = &(segshlist[i]);
    // 'parysh' is the face [lpt, steinerpt, #].
    stpivot(*parysh, searchtet);
    // Skip it if it is outside.
    if (ishulltet(searchtet)) continue;
    
    // Get the "half-ball". Tets in 'cavetetlist' all contain 'steinerpt' as
    //   opposite.  Subfaces in 'caveshlist' all contain 'steinerpt' as apex.
    //   Moreover, subfaces are oriented towards the interior of the ball.
    setpoint2tet(steinerpt, encode(searchtet));
    getvertexstar(0, steinerpt, cavetetlist, NULL, caveshlist);
    
    // Get all tets in this sector.
    for (int j = 0; j < cavetetlist->objects; j++) {
      neightet = * (triface *) fastlookup(cavetetlist, j);
      infect(neightet);
      caveoldtetlist->newindex((void **) &parytet);
      *parytet = neightet;
    }
    cavetetlist->restart();
    caveshlist->restart();

    insertvertexflags ivf;
    searchtet = neightet; // No need point location.
    ivf.iloc = (int) INSTAR;  // No need point location.
    // The following are default options.
    //ivf.bowywat = 0;
    //ivf.lawson = 0;
    //ivf.validflag = 0; // no need to validate cavity.
    //ivf.chkencflag = 0; //chkencflag;
    ivf.assignmeshsize = b->metric; 
    if (ivf.assignmeshsize) {
      // Search the tet containing 'steinerpt' for size interpolation.
      locate(newsteiners[i], &searchtet);
    }

    // Insert the new point into the tetrahedralization T.
    // Note that T is convex (nonconvex = 0).
    if (insertpoint(newsteiners[i], &searchtet, NULL, NULL, &ivf)) {
      // The vertex has been inserted.
      st_volref_count++; 
      if (steinerleft > 0) steinerleft--;
      //return 1;
    } else {
      // Not inserted. 
      //assert(0);
      pointdealloc(newsteiners[i]);
      newsteiners[i] = NULL;
      break; //return 0;
    }
  } // i

  delete [] segshlist;

  if (i < n) {
    //assert(0); 
    delete [] newsteiners;
    return 0;  
  }

  // Now remove the Steiner point from the segment.
  if (!removevertexbyflips(steinerpt)) {
    //assert(0);
    delete [] newsteiners;
    return 0;
  }

  // We've removed a Steiner points.
  setpointtype(steinerpt, UNUSEDVERTEX);
  unuverts++;
  
  int steinercount = 0;

  int bak_fliplinklevel = b->fliplinklevel;
  b->fliplinklevel = 100000; // Unlimited flip level.

  // Try to remove newly added Steiner points.
  for (i = 0; i < n; i++) {
    if (newsteiners[i] != NULL) {
      if (!removevertexbyflips(newsteiners[i])) {
        if (b->supsteiner_level > 0) { // Not -Y/0
          // Save it in subvertstack for removal.
          point *parypt;
          subvertstack->newindex((void **) &parypt);
          *parypt = newsteiners[i];
        }
        steinercount++;
      }
    }
  }

  b->fliplinklevel = bak_fliplinklevel;

  if (steinercount > 0) {
    if (b->verbose > 3) {
      printf("      Added %d interior Steiner points.\n", steinercount);
    }
  }

  delete [] newsteiners;

  return 1;
}


//============================================================================//
//                                                                            //
// suppresssteinerpoints()    Suppress Steiner points.                        //
//                                                                            //
// All Steiner points have been saved in 'subvertstack' in the routines       //
// carveholes() and suppresssteinerpoint().                                   //
// Each Steiner point is either removed or shifted into the interior.         //
//                                                                            //
//============================================================================//

int TetMeshCore::suppresssteinerpoints()
{

  if (!b->quiet) {
    printf("Suppressing Steiner points ...\n");
  }

  point rempt, *parypt;

  int bak_fliplinklevel = b->fliplinklevel;
  b->fliplinklevel = 100000; // Unlimited flip level.
  int suppcount = 0, remcount = 0;
  int i;

  // Try to suppress boundary Steiner points.
  for (i = 0; i < subvertstack->objects; i++) {
    parypt = (point *) fastlookup(subvertstack, i);
    rempt = *parypt;
    if (pointtype(rempt) != UNUSEDVERTEX) {
      if ((pointtype(rempt) == FREESEGVERTEX) || 
          (pointtype(rempt) == FREEFACETVERTEX)) {
        if (suppressbdrysteinerpoint(rempt)) {
          suppcount++;
        }
      }
    }
  } // i

  if (suppcount > 0) {
    if (b->verbose) {
      printf("  Suppressed %d boundary Steiner points.\n", suppcount);
    }
  }

  if (b->supsteiner_level > 0) { // -Y/1
    for (i = 0; i < subvertstack->objects; i++) {
      parypt = (point *) fastlookup(subvertstack, i);
      rempt = *parypt;
      if (pointtype(rempt) != UNUSEDVERTEX) {
        if (pointtype(rempt) == FREEVOLVERTEX) {
          if (removevertexbyflips(rempt)) {
            remcount++;
          }
        }
      }
    }
  }

  if (remcount > 0) {
    if (b->verbose) {
      printf("  Removed %d interior Steiner points.\n", remcount);
    }
  }

  b->fliplinklevel = bak_fliplinklevel;

  if (b->supsteiner_level > 1) { // -Y/2
    // Smooth interior Steiner points.
    optparameters opm;
    triface *parytet;
    point *ppt;
    double ori;
    int smtcount, count, ivcount;
    int nt, j;

    // Point smooth options.
    opm.max_min_volume = 1;
    opm.numofsearchdirs = 20;
    opm.searchstep = 0.001;
    opm.maxiter = 30; // Limit the maximum iterations.

    smtcount = 0;

    do {

      nt = 0;

      while (1) {
        count = 0;
        ivcount = 0; // Clear the inverted count.

        for (i = 0; i < subvertstack->objects; i++) {
          parypt = (point *) fastlookup(subvertstack, i);
          rempt = *parypt;
          if (pointtype(rempt) == FREEVOLVERTEX) {
            getvertexstar(1, rempt, cavetetlist, NULL, NULL);
            // Calculate the initial smallest volume (maybe zero or negative).
            for (j = 0; j < cavetetlist->objects; j++) {
              parytet = (triface *) fastlookup(cavetetlist, j);
              ppt = (point *) &(parytet->tet[4]);
              ori = orient3dfast(ppt[1], ppt[0], ppt[2], ppt[3]);
              if (j == 0) {
                opm.initval = ori;
              } else {
                if (opm.initval > ori) opm.initval = ori; 
              }
            }
            if (smoothpoint(rempt, cavetetlist, 1, &opm)) {
              count++;
            }
            if (opm.imprval <= 0.0) {
              ivcount++; // The mesh contains inverted elements.
            }
            cavetetlist->restart();
          }
        } // i

        smtcount += count;

        if (count == 0) {
          // No point has been smoothed.
          break;
        }

        nt++;
        if (nt > 2) {
          break; // Already three iterations.
        }
      } // while

      if (ivcount > 0) {
        // There are inverted elements!
        if (opm.maxiter > 0) {
          // Set unlimited smoothing steps. Try again.
          opm.numofsearchdirs = 30;
          opm.searchstep = 0.0001;
          opm.maxiter = -1;
          continue;
        }
      }

      break;
    } while (1); // Additional loop for (ivcount > 0)

    if (ivcount > 0) {
      printf("BUG Report!  The mesh contain inverted elements.\n");
    }

    if (b->verbose) {
      if (smtcount > 0) {
        printf("  Smoothed %d Steiner points.\n", smtcount); 
      }
    }
  } // -Y2

  subvertstack->restart();

  return 1;
}

//============================================================================//
//                                                                            //
// recoverboundary()    Recover segments and facets.                          //
//                                                                            //
//============================================================================//

void TetMeshCore::recoverboundary(clock_t& tv)
{
  arraypool *misseglist, *misshlist;
  arraypool *bdrysteinerptlist;
  face searchsh, *parysh;
  face searchseg, *paryseg;
  point rempt, *parypt;
  long ms; // The number of missing segments/subfaces.
  int nit; // The number of iterations.
  int s, i;

  // Counters.
  long bak_segref_count, bak_facref_count, bak_volref_count;

  if (!b->quiet) {
    printf("Recovering boundaries...\n");
  }

  boundary_recovery_flag = 1;
  cos_collinear_ang_tol = cos(b->collinear_ang_tol / 180. * PI);

  if (segmentendpointslist == NULL) {
    // We need segment adjacent information during flips.
    makesegmentendpointsmap();
  }


  if (b->verbose) {
    printf("  Recovering segments.\n");
  }

  // Segments will be introduced.
  checksubsegflag = 1;

  misseglist = new arraypool(sizeof(face), 8);
  bdrysteinerptlist = new arraypool(sizeof(point), 8);

  // In random order.
  subsegs->traversalinit();
  for (i = 0; i < subsegs->items; i++) {
    s = randomnation(i + 1);
    // Move the s-th seg to the i-th.
    subsegstack->newindex((void **) &paryseg);
    *paryseg = * (face *) fastlookup(subsegstack, s);
    // Put i-th seg to be the s-th.
    searchseg.sh = shellfacetraverse(subsegs);
    paryseg = (face *) fastlookup(subsegstack, s);
    *paryseg = searchseg;
  }

  // The init number of missing segments.
  ms = subsegs->items;
  nit = 0; 
  if (b->fliplinklevel < 0) {
    autofliplinklevel = 1; // Init value.
  }

  // First, trying to recover segments by only doing flips.
  while (1) {
    recoversegments(misseglist, 0, 0);

    if (misseglist->objects > 0) {
      if (b->fliplinklevel >= 0) {
        break;
      } else {
        if (misseglist->objects >= ms) {
          nit++;
          if (nit >= 3) {
            //break;
            // Do the last round with unbounded flip link level.
            b->fliplinklevel = 100000;
          }
        } else {
          ms = misseglist->objects;
          if (nit > 0) {
            nit--;
          }
        }
        for (i = 0; i < misseglist->objects; i++) {
          subsegstack->newindex((void **) &paryseg);
          *paryseg = * (face *) fastlookup(misseglist, i);
        }
        misseglist->restart();
        autofliplinklevel+=b->fliplinklevelinc;
      }
    } else {
      // All segments are recovered.
      break;
    }
  } // while (1)

  if (b->verbose) {
    printf("  %ld (%ld) segments are recovered (missing).\n", 
           subsegs->items - misseglist->objects, misseglist->objects);
  }

  if (misseglist->objects > 0) {
    // Second, trying to recover segments by doing more flips (fullsearch).
    while (misseglist->objects > 0) {
      ms = misseglist->objects;
      for (i = 0; i < misseglist->objects; i++) {
        subsegstack->newindex((void **) &paryseg);
        *paryseg = * (face *) fastlookup(misseglist, i);
      }
      misseglist->restart();

      recoversegments(misseglist, 1, 0);

      if (misseglist->objects < ms) {
        // The number of missing segments is reduced.
        continue;
      } else {
        break;
      }
    }
    if (b->verbose) {
      printf("  %ld (%ld) segments are recovered (missing).\n", 
             subsegs->items - misseglist->objects, misseglist->objects);
    }
  }

  //int bak_verbose = b->verbose;
  //if (b->verbose < 3) {
  //  b->verbose = 3; // debug...
  //}

  if (misseglist->objects > 0) {
    // Third, trying to recover segments by doing more flips (fullsearch)
    //   and adding Steiner points in the volume.

    if (b->verbose) {
      printf("  Recovering Delaunay.\n");
    }
    
    recoverdelaunay();
    
    if (b->verbose) {
      printf("  Recovering segments with Steiner points.\n");
    }

    while (misseglist->objects > 0) {
      ms = misseglist->objects;
      for (i = 0; i < misseglist->objects; i++) {
        subsegstack->newindex((void **) &paryseg);
        *paryseg = * (face *) fastlookup(misseglist, i);
      }
      misseglist->restart();

      //recoversegments(misseglist, 1, 1);
      recoversegments(misseglist, 0, 1); // no full search

      if (misseglist->objects < ms) {
        // The number of missing segments is reduced.
        continue;
      } else {
        break;
      }
    }
    if (b->verbose) {
      printf("  Added %ld Steiner points in volume.\n", st_volref_count);
    }
  }

  if (misseglist->objects > 0) {
    // Last, trying to recover segments by doing more flips (fullsearch),
    //   and adding Steiner points in the volume, and splitting segments.
    long bak_inpoly_count = st_volref_count; //st_inpoly_count;

    if (b->verbose) {
      printf("  Recovering Delaunay.\n");
    }
    
    recoverdelaunay();

    if (b->verbose) {
      printf("  Recovering segments with Steiner points.\n");
    }

    while (misseglist->objects > 0) {
      ms = misseglist->objects;
      for (i = 0; i < misseglist->objects; i++) {
        subsegstack->newindex((void **) &paryseg);
        *paryseg = * (face *) fastlookup(misseglist, i);
      }
      misseglist->restart();

      //recoversegments(misseglist, 1, 2);
      recoversegments(misseglist, 0, 2); // no full search

      if (misseglist->objects < ms) {
        // The number of missing segments is reduced.
        continue;
      } else {
        break;
      }
    } // while (misseglist->objects > 0)

    if (b->verbose) {
      printf("  Added %ld Steiner points in segments.\n", st_segref_count);
      if (st_volref_count > bak_inpoly_count) {
        printf("  Added another %ld Steiner points in volume.\n", 
               st_volref_count - bak_inpoly_count);
      }
    }

    // There may be un-recovered subsegments.
    if (misseglist->objects > 0l) {
      if (b->verbose) {
        printf("  !! %ld subsegments are missing.\n", misseglist->objects);
      }
    }
  }

  if (skipped_segment_list != NULL) {
    if (!b->quiet) {
      printf("  Skipped %ld segments due to intersections.\n",
             skipped_segment_list->objects);
    }
    delete skipped_segment_list;
  }


  //b->verbose = bak_verbose; // debug...

  if (st_segref_count > 0) {
    // Try to remove the Steiner points added in segments.
    if (b->verbose) {
      printf("  Suppressing %ld Steiner points in segments.\n", st_segref_count);
    }
    int bak_fliplinklevel = b->fliplinklevel;
    b->fliplinklevel = 20; // limit this value
    
    bak_segref_count = st_segref_count;
    bak_volref_count = st_volref_count;
    for (i = 0; i < subvertstack->objects; i++) {
      // Get the Steiner point.
      parypt = (point *) fastlookup(subvertstack, i);
      rempt = *parypt;
      if (!removevertexbyflips(rempt)) {
        // Save it in list.
        bdrysteinerptlist->newindex((void **) &parypt);
        *parypt = rempt;
      }
    }
    if (b->verbose) {
      if (st_segref_count < bak_segref_count) {
        if (bak_volref_count < st_volref_count) {
          printf("  Suppressed %ld Steiner points in segments.\n", 
                 st_volref_count - bak_volref_count);
        }
        if ((st_segref_count + (st_volref_count - bak_volref_count)) <
            bak_segref_count) {
          printf("  Removed %ld Steiner points in segments.\n", 
                 bak_segref_count - 
                   (st_segref_count + (st_volref_count - bak_volref_count)));
        }
      }
    }
    
    b->fliplinklevel = bak_fliplinklevel; // restore it.
    subvertstack->restart();
  }


  tv = clock();

  if (b->verbose) {
    printf("  Recovering facets.\n");
  }

  // Subfaces will be introduced.
  checksubfaceflag = 1;

  misshlist = new arraypool(sizeof(face), 8);

  // Randomly order the subfaces.
  subfaces->traversalinit();
  for (i = 0; i < subfaces->items; i++) {
    s = randomnation(i + 1);
    // Move the s-th subface to the i-th.
    subfacstack->newindex((void **) &parysh);
    *parysh = * (face *) fastlookup(subfacstack, s);
    // Put i-th subface to be the s-th.
    searchsh.sh = shellfacetraverse(subfaces);
    parysh = (face *) fastlookup(subfacstack, s);
    *parysh = searchsh;
  }

  ms = subfaces->items;
  nit = 0; 
  b->fliplinklevel = -1; // Init.
  if (b->fliplinklevel < 0) {
    autofliplinklevel = 1; // Init value.
  }

  while (1) {
    recoversubfaces(misshlist, 0);

    if (misshlist->objects > 0) {
      if (b->fliplinklevel >= 0) {
        break;
      } else {
        if (misshlist->objects >= ms) {
          nit++;
          if (nit >= 3) {
            //break;
            // Do the last round with unbounded flip link level.
            //b->fliplinklevel = 100000; // this can be very slow.
            if (autofliplinklevel < 30) {
              b->fliplinklevel = 30;
            } else {
              b->fliplinklevel = autofliplinklevel + 30;
            }
          }
        } else {
          ms = misshlist->objects;
          if (nit > 0) {
            nit--;
          }
        }
        for (i = 0; i < misshlist->objects; i++) {
          subfacstack->newindex((void **) &parysh);
          *parysh = * (face *) fastlookup(misshlist, i);
        }
        misshlist->restart();
        autofliplinklevel+=b->fliplinklevelinc;
      }
    } else {
      // All subfaces are recovered.
      break;
    }
  } // while (1)

  if (b->verbose) {
    printf("  %ld (%ld) subfaces are recovered (missing).\n",
           subfaces->items - misshlist->objects, misshlist->objects);
  }

  if (misshlist->objects > 0) {
    // There are missing subfaces. Add Steiner points.

    if (b->verbose) {
      printf("  Recovering Delaunay.\n");
    }
    
    recoverdelaunay();
    
    if (b->verbose) {
      printf("  Recovering facets with Steiner points.\n");
    }
    
    while (misshlist->objects > 0) {
      ms = misshlist->objects;
      for (i = 0; i < misshlist->objects; i++) {
        subfacstack->newindex((void **) &parysh);
        *parysh = * (face *) fastlookup(misshlist, i);
      }
      misshlist->restart();
      
      recoversubfaces(misshlist, 1);
      
      if (misshlist->objects < ms) {
        continue;
      } else {
        break;
      }
    }
    
    if (b->verbose) {
      printf("  %ld (%ld) subfaces are recovered (missing).\n",
             subfaces->items - misshlist->objects, misshlist->objects);
      printf("  Added %ld Steiner points in facets.\n", st_facref_count);
    }
  }

  if (misshlist->objects > 0) {
    long bak_steiner = st_facref_count;
  
    if (b->verbose) {
      printf("  Recovering Delaunay.\n");
    }
    
    recoverdelaunay();
  
    if (b->verbose) {
      printf("  Recovering facets with Steiner points.\n");
    }
  
    while (misshlist->objects > 0) {
      ms = misshlist->objects;
      for (i = 0; i < misshlist->objects; i++) {
        subfacstack->newindex((void **) &parysh);
        *parysh = * (face *) fastlookup(misshlist, i);
      }
      misshlist->restart();
      
      recoversubfaces(misshlist, 2); // steinerflag = 2;
      
      if (misshlist->objects < ms) {
        continue;
      } else {
        break;
      }
    }

    if (subsegstack->objects > 0) {
      // Save unrecovered subsegments.
      triface neightet;
      face checkseg;
      for (i = 0; i < subsegstack->objects; i++) {
        checkseg = * (face *) fastlookup(subsegstack, i);
        if ((checkseg.sh == NULL) ||
            (checkseg.sh[3] == NULL)) continue;
        // Check if this subsegment is missing.
        sstpivot1(checkseg, neightet);
        if (neightet.tet != NULL) continue;
        // Save a missing subsegment.
        misseglist->newindex((void **) &paryseg);
        *paryseg = checkseg;
      }
      subsegstack->restart();
    } // if (subsegstack->objects > 0)

    if (b->verbose) {
      printf("  %ld (%ld) subfaces are recovered (missing).\n",
             subfaces->items - misshlist->objects, misshlist->objects);
      printf("  Added %ld Steiner points in facets.\n",
             st_facref_count - bak_steiner);
    }
  }

  // There may be un-recovered subsegments.
  if (misshlist->objects > 0l) {
    if (b->verbose) {
      printf("  !! %ld subfaces are missing.\n", misshlist->objects);
    }
    terminate_tet_core(this, 2);
    // Save the list of missing subface.
    //missing_tri_list = new arraypool(sizeof(face), 8);
    //for (i = 0; i < misshlist->objects; i++) {
    //  missing_tri_list->newindex((void **) &parysh);
    //  *parysh = * (face *) fastlookup(misshlist, i);
    //}
    //misshlist->restart();
  }
  
  if (duplicated_facets_count > 0l) {
    if (b->verbose) {
      printf("  Deleting %ld duplicated facets.\n", duplicated_facets_count);
    }
    triface neightet, spintet;
    face faceloop, sfaces[256]; // *tmp_sfaces = NULL;
    face sseg;
    int snum, snum_limit = 256;
    int t1ver;
    subfaces->traversalinit();
    faceloop.sh = shellfacetraverse(subfaces);
    while (faceloop.sh != NULL) {
      if (sinfected(faceloop)) {
        // Delete an ignored duplicated subface.
        shellfacedealloc(subfaces, faceloop.sh);
      }
      if (!smarktest3ed(faceloop)) {
        faceloop.shver = 0;
        stpivot(faceloop, neightet);
        if (neightet.tet == NULL) {
          terminate_tet_core(this, 2);
        }
        // Update the subface connections at its three edges.
        for (int k= 0; k < 3; k++) {
          sspivot(faceloop, sseg);
          if (sseg.sh != NULL) {
            ssbond(faceloop, sseg); // Update segment connection.
          }
          // Get all subfaces at this edge.
          snum = 0;
          spintet = neightet;
          do {
            if (issubface(spintet)) {
              tspivot(spintet, sfaces[snum++]);
              if (snum > snum_limit) {
                // Unlikely to happen.
                terminate_tet_core(this, 2);
                //tmp_sfaces = new face[snum_limit * 2];
              }
            }
            fnextself(spintet);
          } while (spintet.tet != neightet.tet);
          // Re-create the face ring.
          for (int j = 0; j < snum - 1; j++) {
            sbond1(sfaces[j], sfaces[j+1]);
          }
          sbond1(sfaces[snum - 1], sfaces[0]);
          enextself(neightet);
          senextself(faceloop);
        } // k
      }
      faceloop.sh = shellfacetraverse(subfaces);
    }
  } // if (duplicated_facets_count > 0l)


  if (st_facref_count > 0) {
    // Try to remove the Steiner points added in facets.
    if (b->verbose) {
      printf("  Suppressing %ld Steiner points in facets.\n", st_facref_count);
    }
    int bak_fliplinklevel = b->fliplinklevel;
    b->fliplinklevel = 30; // limit this value
    
    bak_facref_count = st_facref_count;
    for (i = 0; i < subvertstack->objects; i++) {
      // Get the Steiner point.
      parypt = (point *) fastlookup(subvertstack, i);
      rempt = *parypt;
      if (!removevertexbyflips(*parypt)) {
        // Save it in list.
        bdrysteinerptlist->newindex((void **) &parypt);
        *parypt = rempt;
      }
    }
    if (b->verbose) {
      if (st_facref_count < bak_facref_count) {
        printf("  Removed %ld Steiner points in facets.\n", 
               bak_facref_count - st_facref_count);
      }
    }
    
    b->fliplinklevel = bak_fliplinklevel;
    subvertstack->restart();
  }


  // There may be missing segments and subfaces.
  if (misseglist->objects > 0) {
    triface adjtet;
    face checkseg;
    for (i = 0; i < misseglist->objects; i++) {
      checkseg = * (face *) fastlookup(misseglist, i);
      // A saved missing segment might be split or recovered.
      if ((checkseg.sh == NULL) || (checkseg.sh[3] == NULL)) {
        continue; // it is split.
      }
      sstpivot1(checkseg, adjtet);
      if (adjtet.tet != NULL) {
        continue; // it is recovered.
      }
      // This is a missing segmemt.
      subsegstack->newindex((void **) &paryseg);
      *paryseg = checkseg;
    }
    if (subsegstack->objects > 0) {
      if (!b->quiet && !b->nowarning) {
        printf("Warning:  %ld segments are not recovered.\n", subsegstack->objects);
      }
      //assert(0); // to do...
      subsegstack->restart();
    }
  }


  if (bdrysteinerptlist->objects > 0) {
    if (b->verbose) {
      printf("  %ld Steiner points remained in boundary.\n",
             bdrysteinerptlist->objects);
    }
  } // if


  boundary_recovery_flag = 0;

  // Accumulate the dynamic memory.
  totalworkmemory += (misseglist->totalmemory + misshlist->totalmemory +
                      bdrysteinerptlist->totalmemory);

  delete bdrysteinerptlist;
  delete misseglist;
  delete misshlist;
}

//                                                                            //
//                                                                            //
//== steiner_cxx =============================================================//

} // namespace sqmesh::mesh::tet::detail
