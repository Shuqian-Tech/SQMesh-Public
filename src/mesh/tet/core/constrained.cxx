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

//== constrained_cxx =========================================================//
//                                                                            //
//                                                                            //

//============================================================================//
//                                                                            //
// finddirection()    Find the tet on the path from one point to another.     //
//                                                                            //
// The path starts from 'searchtet''s origin and ends at 'endpt'. On finish,  //
// 'searchtet' contains a tet on the path, its origin does not change.        //
//                                                                            //
// The return value indicates one of the following cases (let 'searchtet' be  //
// abcd, a is the origin of the path):                                        //
//   - ACROSSVERT, edge ab is collinear with the path;                        //
//   - ACROSSEDGE, edge bc intersects with the path;                          //
//   - ACROSSFACE, face bcd intersects with the path.                         //
//                                                                            //
// WARNING: This routine is designed for convex triangulations, and will not  //
// generally work after the holes and concavities have been carved.           //
//                                                                            //
//============================================================================//

enum TetMeshCore::interresult 
  TetMeshCore::finddirection(triface* searchtet, point endpt)
{
  triface neightet;
  point pa, pb, pc, pd;
  enum {HMOVE, RMOVE, LMOVE} nextmove;
  double hori, rori, lori;
  int t1ver;
  int s;

  // The origin is fixed.
  pa = org(*searchtet);
  if ((point) searchtet->tet[7] == dummypoint) {
    // A hull tet. Choose the neighbor of its base face.
    decode(searchtet->tet[3], *searchtet);
    // Reset the origin to be pa.
    if ((point) searchtet->tet[4] == pa) {
      searchtet->ver = 11;
    } else if ((point) searchtet->tet[5] == pa) {
      searchtet->ver = 3;
    } else if ((point) searchtet->tet[6] == pa) {
      searchtet->ver = 7;
    } else {
      searchtet->ver = 0;
    }
  }

  pb = dest(*searchtet);
  // Check whether the destination or apex is 'endpt'.
  if (pb == endpt) {
    // pa->pb is the search edge.
    return ACROSSVERT;
  }

  pc = apex(*searchtet);
  if (pc == endpt) {
    // pa->pc is the search edge.
    eprevesymself(*searchtet);
    return ACROSSVERT;
  }

  // Walk through tets around pa until the right one is found.
  while (1) {

    pd = oppo(*searchtet);
    // Check whether the opposite vertex is 'endpt'.
    if (pd == endpt) {
      // pa->pd is the search edge.
      esymself(*searchtet);
      enextself(*searchtet);
      return ACROSSVERT;
    }
    // Check if we have entered outside of the domain.
    if (pd == dummypoint) {
      // This is possible when the mesh is non-convex.
      if (nonconvex) {
        return ACROSSFACE; // return ACROSSSUB; // Hit a bounday.
      } else {
        terminate_tet_core(this, 2);
      }
    }

    // Now assume that the base face abc coincides with the horizon plane,
    //   and d lies above the horizon.  The search point 'endpt' may lie
    //   above or below the horizon.  We test the orientations of 'endpt'
    //   with respect to three planes: abc (horizon), bad (right plane),
    //   and acd (left plane). 
    hori = orient3d(pa, pb, pc, endpt);
    rori = orient3d(pb, pa, pd, endpt);
    lori = orient3d(pa, pc, pd, endpt);

    // Now decide the tet to move.  It is possible there are more than one
    //   tets are viable moves. Is so, randomly choose one. 
    if (hori > 0) {
      if (rori > 0) {
        if (lori > 0) {
          // Any of the three neighbors is a viable move.
          s = randomnation(3); 
          if (s == 0) {
            nextmove = HMOVE;
          } else if (s == 1) {
            nextmove = RMOVE;
          } else {
            nextmove = LMOVE;
          }
        } else {
          // Two tets, below horizon and below right, are viable.
          if (randomnation(2)) {
            nextmove = HMOVE;
          } else {
            nextmove = RMOVE;
          }
        }
      } else {
        if (lori > 0) {
          // Two tets, below horizon and below left, are viable.
          if (randomnation(2)) {
            nextmove = HMOVE;
          } else {
            nextmove = LMOVE;
          }
        } else {
          // The tet below horizon is chosen.
          nextmove = HMOVE;
        }
      }
    } else {
      if (rori > 0) {
        if (lori > 0) {
          // Two tets, below right and below left, are viable.
          if (randomnation(2)) {
            nextmove = RMOVE;
          } else {
            nextmove = LMOVE;
          }
        } else {
          // The tet below right is chosen.
          nextmove = RMOVE;
        }
      } else {
        if (lori > 0) {
          // The tet below left is chosen.
          nextmove = LMOVE;
        } else {
          // 'endpt' lies either on the plane(s) or across face bcd.
          if (hori == 0) {
            if (rori == 0) {
              // pa->'endpt' is COLLINEAR with pa->pb.
              return ACROSSVERT;
            }
            if (lori == 0) {
              // pa->'endpt' is COLLINEAR with pa->pc.
              eprevesymself(*searchtet); // [a,c,d]
              return ACROSSVERT;
            }
            // pa->'endpt' crosses the edge pb->pc.
            return ACROSSEDGE;
          }
          if (rori == 0) {
            if (lori == 0) {
              // pa->'endpt' is COLLINEAR with pa->pd.
              esymself(*searchtet); // face bad.
              enextself(*searchtet); // face [a,d,b]
              return ACROSSVERT;
            }
            // pa->'endpt' crosses the edge pb->pd.
            esymself(*searchtet); // face bad.
            enextself(*searchtet); // face adb
            return ACROSSEDGE;
          }
          if (lori == 0) {
            // pa->'endpt' crosses the edge pc->pd.
            eprevesymself(*searchtet); // [a,c,d]
            return ACROSSEDGE;
          }
          // pa->'endpt' crosses the face bcd.
          return ACROSSFACE;
        }
      }
    }

    // Move to the next tet, fix pa as its origin.
    if (nextmove == RMOVE) {
      fnextself(*searchtet);
    } else if (nextmove == LMOVE) {
      eprevself(*searchtet);
      fnextself(*searchtet);
      enextself(*searchtet);
    } else { // HMOVE
      fsymself(*searchtet);
      enextself(*searchtet);
    }
    if (org(*searchtet) != pa) {
      terminate_tet_core(this, 2);
    }
    pb = dest(*searchtet);
    pc = apex(*searchtet);

  } // while (1)

}

//============================================================================//
//                                                                            //
// scoutsegment()    Search an edge in the tetrahedralization.                //
//                                                                            //
// If the edge is found, it returns SHAREEDGE, and 'searchtet' returns the    //
// edge from startpt to endpt.                                                //
//                                                                            //
// If the edge is missing, it returns either ACROSSEDGE or ACROSSFACE, which  //
// indicates that the edge intersects an edge or a face.  If 'refpt' is NULL, //
// 'searchtet' returns the edge or face. If 'refpt' is not NULL, it returns   //
// a vertex which encroaches upon this edge, and 'searchtet' returns a tet    //
// which containing 'refpt'.                                                  //
//                                                                            //
// The parameter 'sedge' is used to report self-intersection. It is the       //
// whose endpoints are 'startpt' and 'endpt'. It must not be a NULL.          //
//                                                                            //
//============================================================================//

enum TetMeshCore::interresult TetMeshCore::scoutsegment(point startpt,point endpt,
  face *sedge, triface* searchtet, point* refpt, arraypool* intfacelist)
{
  point pd;
  enum interresult dir;
  int t1ver;

  if (b->verbose > 2) {
    printf("      Scout seg (%d, %d).\n",pointmark(startpt),pointmark(endpt));
  }

  point2tetorg(startpt, *searchtet);
  dir = finddirection(searchtet, endpt);

  if (dir == ACROSSVERT) {
    pd = dest(*searchtet);
    if (pd == endpt) {
      if (issubseg(*searchtet)) {
        //report_selfint_edge(startpt, endpt, sedge, searchtet, dir);
        terminate_tet_core(this, 3);
      }
      return SHAREEDGE;
    } else {
      // A point is on the path.
      //report_selfint_edge(startpt, endpt, sedge, searchtet, dir);
      terminate_tet_core(this, 3);
      return ACROSSVERT;
    }
  }

  // dir is either ACROSSEDGE or ACROSSFACE.
  enextesymself(*searchtet); // Go to the opposite face.
  fsymself(*searchtet); // Enter the adjacent tet.

  if (dir == ACROSSEDGE) {
    // Check whether two segments are intersecting.
    if (issubseg(*searchtet)) {
      //report_selfint_edge(startpt, endpt, sedge, searchtet, dir);
      terminate_tet_core(this, 3);
    }
  } else if (dir == ACROSSFACE) {
    if (checksubfaceflag) {
      // Check whether a segment and a subface are intersecting.
      if (issubface(*searchtet)) {
        //report_selfint_edge(startpt, endpt, sedge, searchtet, dir);
        terminate_tet_core(this, 3);
      }
    }
  } else {
    terminate_tet_core(this, 2);
  }

  if (refpt == NULL) {
    // Do not need a reference point. Return.
    return dir;
  }

  triface neightet, reftet;
  point pa, pb, pc;
  double angmax, ang;
  int types[2], poss[4];
  int pos = 0, i, j;

  pa = org(*searchtet);
  angmax = interiorangle(pa, startpt, endpt, NULL);
  *refpt = pa;
  pb = dest(*searchtet);
  ang = interiorangle(pb, startpt, endpt, NULL);
  if (ang > angmax) {
    angmax = ang;
    *refpt = pb;
  }
  pc = apex(*searchtet);
  ang = interiorangle(pc, startpt, endpt, NULL);
  if (ang > angmax) {
    angmax = ang;
    *refpt = pc;
  }
  reftet = *searchtet; // Save the tet containing the refpt.

  // Search intersecting faces along the segment.
  while (1) {


    pd = oppo(*searchtet);


    // Stop if we meet 'endpt'.
    if (pd == endpt) break;

    ang = interiorangle(pd, startpt, endpt, NULL);
    if (ang > angmax) {
      angmax = ang;
      *refpt = pd;
      reftet = *searchtet;
    }

    // Find a face intersecting the segment.
    if (dir == ACROSSFACE) {
      // One of the three oppo faces in 'searchtet' intersects the segment.
      neightet = *searchtet;
      j = (neightet.ver & 3); // j is the current face number.
      for (i = j + 1; i < j + 4; i++) {
        neightet.ver = (i % 4);
        pa = org(neightet);
        pb = dest(neightet);
        pc = apex(neightet);
        pd = oppo(neightet); // The above point.
        if (tri_edge_test(pa, pb, pc, startpt, endpt, pd, 1, types, poss)) {
          dir = (enum interresult) types[0];
          pos = poss[0];
          break;
        } else {
          dir = DISJOINT;
          pos = 0;
        }
      }
    } else if (dir == ACROSSEDGE) {
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
        if (tri_edge_test(pa, pb, pc, startpt, endpt, pd, 1, types, poss)) {
          dir = (enum interresult) types[0];
          pos = poss[0];
          break;
        } else {
          dir = DISJOINT;
          pos = 0;
        }
      }
      if (dir == DISJOINT) {
        // No intersection. Rotate to the next tet at the edge.
        dir = ACROSSEDGE;
        fnextself(*searchtet);
        continue;
      }
    }

    if (dir == ACROSSVERT) {
      // This segment passing a vertex. Choose it and return.
      for (i = 0; i < pos; i++) {
        enextself(neightet);
      }
      eprev(neightet, *searchtet); 
      // dest(*searchtet) lies on the segment.
      //report_selfint_edge(startpt, endpt, sedge, searchtet, dir);
      terminate_tet_core(this, 3);
      return ACROSSVERT;
    } else if (dir == ACROSSEDGE) {
      // Get the edge intersects with the segment.
      for (i = 0; i < pos; i++) {
        enextself(neightet);
      }
    }
    // Go to the next tet.
    fsym(neightet, *searchtet);

    if (dir == ACROSSEDGE) {
      // Check whether two segments are intersecting.
      if (issubseg(*searchtet)) {
        //report_selfint_edge(startpt, endpt, sedge, searchtet, dir);
        terminate_tet_core(this, 3);
      }  
    } else if (dir == ACROSSFACE) {
      if (checksubfaceflag) {
        // Check whether a segment and a subface are intersecting.
        if (issubface(*searchtet)) {
          //report_selfint_edge(startpt, endpt, sedge, searchtet, dir);
          terminate_tet_core(this, 3);
        }
      }
    } else {
      terminate_tet_core(this, 2);
    }

  } // while (1)

  // A valid reference point should inside the diametrial circumsphere of
  //   the missing segment, i.e., it encroaches upon it.
  if (2.0 * angmax < PI) {
    *refpt = NULL;
  }


  *searchtet = reftet;
  return dir;
}

//============================================================================//
//                                                                            //
// getsteinerpointonsegment()    Get a Steiner point on a segment.            //
//                                                                            //
// Return '1' if 'refpt' lies on an adjacent segment of this segment. Other-  //
// wise, return '0'.                                                          //
//                                                                            //
//============================================================================//

int TetMeshCore::getsteinerptonsegment(face* seg, point refpt, point steinpt)
{
  point ei = sorg(*seg);
  point ej = sdest(*seg);
  int adjflag = 0, i;

  if (refpt != NULL) {
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
        // Create a Steiner point at the intersection of the segment
        //   [far_ei, far_ej] and the sphere centered at far_ei with 
        //   radius |far_ei - refpt|.
        L = distance(far_ei, far_ej);
        L1 = distance(far_ei, refpt);
        t = L1 / L;
        for (i = 0; i < 3; i++) {
          steinpt[i] = far_ei[i] + t * (far_ej[i] - far_ei[i]);
        }
        adjflag = 1;
      } else if ((far_pi == far_ej) || (far_pj == far_ej)) {
        L = distance(far_ei, far_ej);
        L1 = distance(far_ej, refpt);
        t = L1 / L;
        for (i = 0; i < 3; i++) {
          steinpt[i] = far_ej[i] + t * (far_ei[i] - far_ej[i]);
        }
        adjflag = 1;
      } else {
        // Cut the segment by the projection point of refpt.
        projpt2edge(refpt, ei, ej, steinpt);
      }
    } else {
      // Cut the segment by the projection point of refpt.
      projpt2edge(refpt, ei, ej, steinpt);
    }

    // Make sure that steinpt is not too close to ei and ej.
    L = distance(ei, ej);
    L1 = distance(steinpt, ei);
    t = L1 / L;
    if ((t < 0.2) || (t > 0.8)) {
      // Split the point at the middle.
      for (i = 0; i < 3; i++) {
        steinpt[i] = ei[i] + 0.5 * (ej[i] - ei[i]);
      }
    }
  } else {
    // Split the point at the middle.
    for (i = 0; i < 3; i++) {
      steinpt[i] = ei[i] + 0.5 * (ej[i] - ei[i]);
    }
  }


  return adjflag;
}



//============================================================================//
//                                                                            //
// delaunizesegments()    Recover segments in a DT.                           //
//                                                                            //
// All segments need to be recovered are in 'subsegstack' (Q).  They will be  //
// be recovered one by one (in a random order).                               //
//                                                                            //
// Given a segment s in the Q, this routine first queries s in the DT, if s   //
// matches an edge in DT, it is 'locked' at the edge. Otherwise, s is split   //
// by inserting a new point p in both the DT and itself. The two new subseg-  //
// ments of s are queued in Q.  The process continues until Q is empty.       //
//                                                                            //
//============================================================================//

void TetMeshCore::delaunizesegments()
{
  triface searchtet, spintet;
  face searchsh;
  face sseg, *psseg;
  point refpt, newpt;
  enum interresult dir;
  insertvertexflags ivf;
  int t1ver; 


  ivf.bowywat = 1; // Use Bowyer-Watson insertion.
  ivf.sloc = (int) ONEDGE; // on 'sseg'.
  ivf.sbowywat = 1; // Use Bowyer-Watson insertion.
  ivf.assignmeshsize = b->metric;
  ivf.smlenflag = useinsertradius; // Return the closet mesh vertex.

  // Loop until 'subsegstack' is empty.
  while (subsegstack->objects > 0l) {
    // seglist is used as a stack.
    subsegstack->objects--;
    psseg = (face *) fastlookup(subsegstack, subsegstack->objects);
    sseg = *psseg;

    // Check if this segment has been recovered.
    sstpivot1(sseg, searchtet);
    if (searchtet.tet != NULL) {
      continue; // Not a missing segment.
    }

    // Search the segment.
    dir = scoutsegment(sorg(sseg), sdest(sseg), &sseg,&searchtet,&refpt,NULL);

    if (dir == SHAREEDGE) {
      // Found this segment, insert it.
      // Let the segment remember an adjacent tet.
      sstbond1(sseg, searchtet);
      // Bond the segment to all tets containing it.
      spintet = searchtet;
      do {
        tssbond1(spintet, sseg);
        fnextself(spintet);
      } while (spintet.tet != searchtet.tet);
    } else {
      if ((dir == ACROSSFACE) || (dir == ACROSSEDGE)) {
        // The segment is missing. Split it.
        // Create a new point.
        makepoint(&newpt, FREESEGVERTEX);
        //setpointtype(newpt, FREESEGVERTEX);
        getsteinerptonsegment(&sseg, refpt, newpt);

        // Start searching from 'searchtet'.
        ivf.iloc = (int) OUTSIDE;
        // Insert the new point into the tetrahedralization T.
        //   Missing segments and subfaces are queued for recovery.
        //   Note that T is convex (nonconvex = 0).
        if (insertpoint(newpt, &searchtet, &searchsh, &sseg, &ivf)) {
          // The new point has been inserted.
          st_segref_count++;
          if (steinerleft > 0) steinerleft--;
          if (useinsertradius) {
            //save_segmentpoint_insradius(newpt, ivf.parentpt, ivf.smlen);
          } 
        } else {
          if (ivf.iloc == (int) NEARVERTEX) {
            // The new point (in the segment) is very close to an existing 
            //   vertex -- a small feature is detected.
            point nearpt = org(searchtet);
            if (pointtype(nearpt) == FREESEGVERTEX) {
              face parentseg;
              sdecode(point2sh(nearpt), parentseg);
              point p1 = farsorg(sseg);
              point p2 = farsdest(sseg);
              point p3 = farsorg(parentseg);
              point p4 = farsdest(parentseg);
              printf("Two segments are very close to each other.\n");
              printf("  Segment 1: [%d, %d] #%d\n", pointmark(p1),
                     pointmark(p2), shellmark(sseg));
              printf("  Segment 2: [%d, %d] #%d\n", pointmark(p3),
                     pointmark(p4), shellmark(parentseg));
              terminate_tet_core(this, 4);
            } else {
              terminate_tet_core(this, 2);
            }
          } else if (ivf.iloc == (int) ONVERTEX) {
            // The new point (in the segment) is coincident with an existing 
            //   vertex -- a self-intersection is detected.
            eprevself(searchtet);
            //report_selfint_edge(sorg(sseg), sdest(sseg), &sseg, &searchtet,
            //                    ACROSSVERT);
            terminate_tet_core(this, 3);
          } else {
            // An unknown case. Report a bug.
            terminate_tet_core(this, 2);
          }
        }
      } else {
        // An unknown case. Report a bug.
        terminate_tet_core(this, 2);
      }
    }
  } // while
}

//============================================================================//
//                                                                            //
// scoutsubface()    Search subface in the tetrahedralization.                //
//                                                                            //
// 'searchsh' is searched in T. If it exists, it is 'locked' at the face in   //
// T. 'searchtet' refers to the face. Otherwise, it is missing.               //
//                                                                            //
// The parameter 'shflag' indicates whether 'searchsh' is a boundary face or  //
// not. It is possible that 'searchsh' is a temporarily subface that is used  //
// as a cavity boundary face.                                                 //
//                                                                            //
//============================================================================//

int TetMeshCore::scoutsubface(face* searchsh, triface* searchtet, int shflag)
{
  point pa = sorg(*searchsh);
  point pb = sdest(*searchsh);

  // Get a tet whose origin is a.
  point2tetorg(pa, *searchtet);
  // Search the edge [a,b].
  enum interresult dir = finddirection(searchtet, pb);
  if (dir == ACROSSVERT) {
    // Check validity of a PLC.
    if (dest(*searchtet) != pb) {
	  if (shflag) {
        // A vertex lies on the search edge. 
	    //report_selfint_edge(pa, pb, searchsh, searchtet, dir);
        terminate_tet_core(this, 3);
	  } else {
	    terminate_tet_core(this, 2);
	  }
    }
	int t1ver;
    // The edge exists. Check if the face exists.
    point pc = sapex(*searchsh);
    // Searchtet holds edge [a,b]. Search a face with apex c.
    triface spintet = *searchtet;
    while (1) {
      if (apex(spintet) == pc) {
        // Found a face matching to 'searchsh'!
        if (!issubface(spintet)) {
          // Insert 'searchsh'.
          tsbond(spintet, *searchsh);
          fsymself(spintet);
          sesymself(*searchsh);
          tsbond(spintet, *searchsh);
          *searchtet = spintet;
          return 1; 
        } else {
          terminate_tet_core(this, 2);
        }
      }
      fnextself(spintet);
      if (spintet.tet == searchtet->tet) break;
    }
  }

  return 0;
}

//============================================================================//
//                                                                            //
// formregion()    Form the missing region of a missing subface.              //
//                                                                            //
// 'missh' is a missing subface. From it we form a missing region R which is  //
// a connected region formed by a set of missing subfaces of a facet.         //
// Comment: There should be no segment inside R.                              //
//                                                                            //
// 'missingshs' returns the list of subfaces in R. All subfaces in this list  //
// are oriented as the 'missh'.  'missingshbds' returns the list of boundary  //
// edges (tetrahedral handles) of R.  'missingshverts' returns all vertices   //
// of R. They are all pmarktested.                                            //
//                                                                            //
// Except the first one (which is 'missh') in 'missingshs', each subface in   //
// this list represents an internal edge of R, i.e., it is missing in the     //
// tetrahedralization. Since R may contain interior vertices, not all miss-   //
// ing edges can be found by this way.                                        //
//============================================================================//

void TetMeshCore::formregion(face* missh, arraypool* missingshs,  
                            arraypool* missingshbds, arraypool* missingshverts)
{
  triface searchtet, spintet;
  face neighsh, *parysh;
  face neighseg, fakeseg;
  point pa, pb, *parypt;
  enum interresult dir;
  int t1ver;
  int i, j;

  smarktest(*missh);
  missingshs->newindex((void **) &parysh);
  *parysh = *missh;

  // Incrementally find other missing subfaces.
  for (i = 0; i < missingshs->objects; i++) {
    missh = (face *) fastlookup(missingshs, i);
    for (j = 0; j < 3; j++) {
      pa = sorg(*missh);
      pb = sdest(*missh);
      point2tetorg(pa, searchtet);
      dir = finddirection(&searchtet, pb);
      if (dir != ACROSSVERT) {
        // This edge is missing. Its neighbor is a missing subface.
        spivot(*missh, neighsh);
        if (!smarktested(neighsh)) {
          // Adjust the face orientation.
          if (sorg(neighsh) != pb) sesymself(neighsh);
          smarktest(neighsh);
          missingshs->newindex((void **) &parysh);
          *parysh = neighsh;
        }
      } else {
        if (dest(searchtet) != pb) {
		  // Report a PLC problem.
		  //report_selfint_edge(pa, pb, missh, &searchtet, dir);
          terminate_tet_core(this, 3);
        }
      }
      // Collect the vertices of R.
      if (!pmarktested(pa)) {
        pmarktest(pa);
        missingshverts->newindex((void **) &parypt);
        *parypt = pa;
      }
      senextself(*missh);
    } // j
  } // i

  // Get the boundary edges of R.
  for (i = 0; i < missingshs->objects; i++) {
    missh = (face *) fastlookup(missingshs, i);
    for (j = 0; j < 3; j++) {
      spivot(*missh, neighsh);
      if ((neighsh.sh == NULL) || !smarktested(neighsh)) {
        // A boundary edge of R.
        // Let the segment point to the adjacent tet.
        point2tetorg(sorg(*missh), searchtet);
        finddirection(&searchtet, sdest(*missh));
        missingshbds->newindex((void **) &parysh);
        *parysh = *missh;
        // Check if this edge is a segment.
        sspivot(*missh, neighseg); 
        if (neighseg.sh == NULL) {
          // Temporarily create a segment at this edge.
          makeshellface(subsegs, &fakeseg);
          setsorg(fakeseg, sorg(*missh));
          setsdest(fakeseg, sdest(*missh));
          sinfect(fakeseg); // Mark it as faked.
          // Connect it to all tets at this edge.
          spintet = searchtet;
          while (1) {
            tssbond1(spintet, fakeseg);
            fnextself(spintet);
            if (spintet.tet == searchtet.tet) break;
          }
          neighseg = fakeseg;
        }
        // Let the segment and the boundary edge point to each other.
        ssbond(*missh, neighseg);
        sstbond1(neighseg, searchtet);
      }
      senextself(*missh);
    } // j
  } // i


  // Unmarktest collected missing subfaces.
  for (i = 0; i < missingshs->objects; i++) {
    parysh = (face *) fastlookup(missingshs, i);
    sunmarktest(*parysh);
  }
}

//============================================================================//
//                                                                            //
// scoutcrossedge()    Search an edge that crosses the missing region.        //
//                                                                            //
// Return 1 if a crossing edge is found. It is returned by 'crosstet'. More-  //
// over, the edge is oriented such that its origin lies below R.  Return 0    //
// if no such edge is found.                                                  //
//                                                                            //
// Assumption: All vertices of the missing region are marktested.             //
//                                                                            //
//============================================================================//

int TetMeshCore::scoutcrossedge(triface& crosstet, arraypool* missingshbds, 
                               arraypool* missingshs)
{
  triface searchtet, spintet, neightet;
  face oldsh, searchsh, *parysh;
  face neighseg;
  point pa, pb, pc, pd, pe;
  double ori;
  int types[2], poss[4];
  int searchflag, interflag;
  int t1ver;
  int i, j;

  searchflag = 0;

  // Search the first new subface to fill the region.
  for (i = 0; i < missingshbds->objects && !searchflag; i++) {
    parysh = (face *) fastlookup(missingshbds, i);
    sspivot(*parysh, neighseg);
    sstpivot1(neighseg, searchtet);
    if (org(searchtet) != sorg(*parysh)) {
      esymself(searchtet);
    }
    spintet = searchtet;
    while (1) {
      if (pmarktested(apex(spintet))) {
        // A possible interior face.
        neightet = spintet;
        oldsh = *parysh;
		// Try to recover an interior edge.
		for (j = 0; j < 2; j++) {
		  enextself(neightet);
          if (!issubseg(neightet)) {
			if (j == 0) {
              senext(oldsh, searchsh);
            } else {
              senext2(oldsh, searchsh);
              sesymself(searchsh);
		      esymself(neightet);
            }
            // Calculate a lifted point. 
			pa = sorg(searchsh);
            pb = sdest(searchsh);
            pc = sapex(searchsh);
            pd = dest(neightet);
            calculateabovepoint4(pa, pb, pc, pd);
            // The lifted point must lie above 'searchsh'.
			ori = orient3d(pa, pb, pc, dummypoint);
			if (ori > 0) {
		      sesymself(searchsh);
		      senextself(searchsh);
		    } else if (ori == 0) {
		      terminate_tet_core(this, 2);
		    }
			if (sscoutsegment(&searchsh,dest(neightet),0,0,1)==SHAREEDGE) {
              // Insert a temp segment to protect the recovered edge.
              face tmpseg;
              makeshellface(subsegs, &tmpseg);
              ssbond(searchsh, tmpseg);
              spivotself(searchsh);
              ssbond(searchsh, tmpseg);
              // Recover locally Delaunay edges.
              lawsonflip();
              // Delete the tmp segment.
              spivot(tmpseg, searchsh);
              ssdissolve(searchsh);
              spivotself(searchsh);
              ssdissolve(searchsh);
              shellfacedealloc(subsegs, tmpseg.sh);
			  searchflag = 1;
            } else {
              // Undo the performed flips.
              if (flipstack != NULL) {
                lawsonflip();
              }
            }
            break;
          } // if (!issubseg(neightet))
		} // j
        if (searchflag) break;
      } // if (pmarktested(apex(spintet)))
      fnextself(spintet);
      if (spintet.tet == searchtet.tet) break;
    }
  } // i

  if (searchflag) {
	// Remove faked segments.
    face checkseg;
    // Remark: We should not use the array 'missingshbds', since the flips may
    //   change the subfaces. We search them from the subfaces in R.
    for (i = 0; i < missingshs->objects; i++) {
      parysh = (face *) fastlookup(missingshs, i);
      oldsh = *parysh;
      for (j = 0; j < 3; j++) {
        if (isshsubseg(oldsh)) {
          sspivot(oldsh, checkseg);
          if (sinfected(checkseg)) {
            // It's a faked segment. Delete it.
            sstpivot1(checkseg, searchtet);
            spintet = searchtet;
            while (1) {
              tssdissolve1(spintet);
              fnextself(spintet);
              if (spintet.tet == searchtet.tet) break;
            }
            shellfacedealloc(subsegs, checkseg.sh);
            ssdissolve(oldsh);
          }
        }
        senextself(oldsh);
      } // j
    }

    fillregioncount++;

    return 0;
  } // if (i < missingshbds->objects)

  searchflag = -1;

  for (j = 0; j < missingshbds->objects && (searchflag == -1); j++) {
    parysh = (face *) fastlookup(missingshbds, j);
    sspivot(*parysh, neighseg);
    sstpivot1(neighseg, searchtet);
    interflag = 0;
    // Let 'spintet' be [#,#,d,e] where [#,#] is the boundary edge of R.
    spintet = searchtet;
    while (1) {
      pd = apex(spintet);
      pe = oppo(spintet);
      // Skip a hull edge.
      if ((pd != dummypoint) && (pe != dummypoint)) {
        // Skip an edge containing a vertex of R.
        if (!pmarktested(pd) && !pmarktested(pe)) {
          // Check if [d,e] intersects R.
          for (i = 0; i < missingshs->objects && !interflag; i++) {
            parysh = (face *) fastlookup(missingshs, i);
            pa = sorg(*parysh);
            pb = sdest(*parysh);
            pc = sapex(*parysh);
            interflag=tri_edge_test(pa, pb, pc, pd, pe, NULL, 1, types, poss);
            if (interflag > 0) { 
              if (interflag == 2) {
                // They intersect at a single point.
				if ((types[0] == (int) ACROSSFACE) ||
				    (types[0] == (int) ACROSSEDGE)) {
                  // Go to the crossing edge [d,e,#,#].
                  edestoppo(spintet, crosstet); // // [d,e,#,#].
                  if (issubseg(crosstet)) {
				    // It is a segment. Report a PLC problem.
				    //report_selfint_face(pa, pb, pc, parysh, &crosstet,
					//                    interflag, types, poss);
                    terminate_tet_core(this, 3);
                  } else {
				    triface chkface = crosstet;
					while (1) {
					  if (issubface(chkface)) break;
					  fsymself(chkface);
					  if (chkface.tet == crosstet.tet) break;
					}
					if (issubface(chkface)) {
					  // Two subfaces are intersecting.
                      //report_selfint_face(pa, pb, pc, parysh, &chkface,
					  //                  interflag, types, poss);
                      terminate_tet_core(this, 3);
					}
				  }
                  // Adjust the edge such that d lies below [a,b,c].
                  ori = orient3d(pa, pb, pc, pd);
                  if (ori < 0) {
                    esymself(crosstet);
                  }
                  searchflag = 1;                  
                } else {
                  // An improper intersection type, ACROSSVERT, TOUCHFACE,
				  //   TOUCHEDGE, SHAREVERT, ...
                  // Maybe it is due to a PLC problem.
				  //report_selfint_face(pa, pb, pc, parysh, &crosstet,
				  //	                  interflag, types, poss);
                  terminate_tet_core(this, 3);
                }
              }
              break;
            } // if (interflag > 0)
          }
        } 
      }
      // Leave search at this bdry edge if an intersection is found.
      if (interflag > 0) break;
      // Go to the next tetrahedron.
      fnextself(spintet);
      if (spintet.tet == searchtet.tet) break; 
    } // while (1)
  } // j

  return searchflag;
}

//============================================================================//
//                                                                            //
// formcavity()    Form the cavity of a missing region.                       //
//                                                                            //
// The missing region R is formed by a set of missing subfaces 'missingshs'.  //
// In the following, we assume R is horizontal and oriented. (All subfaces    //
// of R are oriented in the same way.)  'searchtet' is a tetrahedron [d,e,#,  //
// #] which intersects R in its interior, where the edge [d,e] intersects R,  //
// and d lies below R.                                                        //
//                                                                            //
// 'crosstets' returns the set of crossing tets. Every tet in it has the      //
// form [d,e,#,#] where [d,e] is a crossing edge, and d lies below R.  The    //
// set of tets form the cavity C, which is divided into two parts by R, one   //
// at top and one at bottom. 'topfaces' and 'botfaces' return the upper and   //
// lower boundary faces of C. 'toppoints' contains vertices of 'crosstets'    //
// in the top part of C, and so does 'botpoints'. Both 'toppoints' and        //
// 'botpoints' contain vertices of R.                                         //
//                                                                            //
// Important: This routine assumes all vertices of the facet containing this  //
// subface are marked, i.e., pmarktested(p) returns true.                     //
//                                                                            //
//============================================================================//

bool TetMeshCore::formcavity(triface* searchtet, arraypool* missingshs,
                            arraypool* crosstets, arraypool* topfaces, 
                            arraypool* botfaces, arraypool* toppoints, 
                            arraypool* botpoints)
{
  arraypool *crossedges;
  triface spintet, neightet, chkface, *parytet;
  face *parysh = NULL;
  point pa, pd, pe, *parypt;
  bool testflag, invalidflag;
  int intflag, types[2], poss[4];
  int t1ver;
  int i, j, k;

  // Temporarily re-use 'topfaces' for all crossing edges.
  crossedges = topfaces;

  if (b->verbose > 2) {
    printf("      Form the cavity of a missing region.\n"); 
  }
  // Mark this edge to avoid testing it later.
  markedge(*searchtet);
  crossedges->newindex((void **) &parytet);
  *parytet = *searchtet;

  invalidflag = 0; 
  // Collect all crossing tets.  Each cross tet is saved in the standard
  //   form [d,e,#,#], where [d,e] is a crossing edge, d lies below R.
  //   NEITHER d NOR e is a vertex of R (!pmarktested). 
  for (i = 0; i < crossedges->objects && !invalidflag; i++) {
    // Get a crossing edge [d,e,#,#].
    searchtet = (triface *) fastlookup(crossedges, i);
    // Sort vertices into the bottom and top arrays.
    pd = org(*searchtet);
    if (!pinfected(pd)) {
      pinfect(pd);
      botpoints->newindex((void **) &parypt);
      *parypt = pd;
    }
    pe = dest(*searchtet);
    if (!pinfected(pe)) {
      pinfect(pe);
      toppoints->newindex((void **) &parypt);
      *parypt = pe;
    }

    // All tets sharing this edge are crossing tets.
    spintet = *searchtet;
    while (1) {
      if (!infected(spintet)) {
        infect(spintet);
        crosstets->newindex((void **) &parytet);
        *parytet = spintet;
      }
      // Go to the next crossing tet.
      fnextself(spintet);
      if (spintet.tet == searchtet->tet) break;
    } // while (1)

    // Detect new crossing edges.
    spintet = *searchtet;
    while (1) {
      // spintet is [d,e,a,#], where d lies below R, and e lies above R. 
      pa = apex(spintet);
      if (pa != dummypoint) {
        if (!pmarktested(pa)) {
	      // There exists a crossing edge, either [e,a] or [a,d]. First check
          //   if the crossing edge has already be added, i.e.,to check if one
		  //   of the tetrahedron at this edge has been marked.
          testflag = true;
          for (j = 0; j < 2 && testflag; j++) {
            if (j == 0) {
              enext(spintet, neightet);
            } else {
              eprev(spintet, neightet);
            }
            while (1) {
              if (edgemarked(neightet)) {
                // This crossing edge has already been tested. Skip it.
                testflag = false;
                break;
              }
              fnextself(neightet);
              if (neightet.tet == spintet.tet) break;
            }
          } // j
          if (testflag) {
            // Test if [e,a] or [a,d] intersects R.
            // Do a brute-force search in the set of subfaces of R. Slow!
            //   Need to be improved!
            pd = org(spintet);
            pe = dest(spintet);
            for (k = 0; k < missingshs->objects; k++) {
              parysh = (face *) fastlookup(missingshs, k);
              intflag = tri_edge_test(sorg(*parysh), sdest(*parysh), 
                          sapex(*parysh), pe, pa, NULL, 1, types, poss);
              if (intflag > 0) {
                // Found intersection. 'a' lies below R.
				if (intflag == 2) {
                  enext(spintet, neightet);
				  if ((types[0] == (int) ACROSSFACE) || 
                      (types[0] == (int) ACROSSEDGE)) {
                    // Only this case is valid.
                  } else {
                    // A non-valid intersection. Maybe a PLC problem.
                    invalidflag = 1;
                  }
				} else {
				  // Coplanar intersection. Maybe a PLC problem.
				  invalidflag = 1;
				}
                break;
              }
              intflag = tri_edge_test(sorg(*parysh), sdest(*parysh), 
			              sapex(*parysh), pa, pd, NULL, 1, types, poss);
              if (intflag > 0) {
                // Found intersection. 'a' lies above R.
				if (intflag == 2) {
                  eprev(spintet, neightet);
				  if ((types[0] == (int) ACROSSFACE) || 
                      (types[0] == (int) ACROSSEDGE)) {
					// Only this case is valid.
                  } else {
				    // A non-valid intersection. Maybe a PLC problem.
                    invalidflag = 1;
                  }
				} else {
				  // Coplanar intersection. Maybe a PLC problem.
				  invalidflag = 1;
				}
                break;
              }
            } // k
            if (k < missingshs->objects) {
              // Found a pair of triangle - edge intersection.
              if (invalidflag) {
                break; // the while (1) loop
              }
              // Adjust the edge direction, so that its origin lies below R,
              //   and its destination lies above R.
              esymself(neightet);
              // This edge may be a segment.
              if (issubseg(neightet)) {
				//report_selfint_face(sorg(*parysh), sdest(*parysh),
                //  sapex(*parysh),parysh,&neightet,intflag,types,poss);
                terminate_tet_core(this, 3);
              }
			  // Check if it is an edge of a subface.
			  chkface = neightet;
              while (1) {
                if (issubface(chkface)) break;
                fsymself(chkface);
                if (chkface.tet == neightet.tet) break;
              }
              if (issubface(chkface)) {
                // Two subfaces are intersecting.
                //report_selfint_face(sorg(*parysh), sdest(*parysh),
                //  sapex(*parysh),parysh,&chkface,intflag,types,poss);
                terminate_tet_core(this, 3);
              }
			  
              // Mark this edge to avoid testing it again.
              markedge(neightet);
              crossedges->newindex((void **) &parytet);
              *parytet = neightet;            
            } else {
              // No intersection is found. It may be a PLC problem.
              invalidflag = 1;
              break; // the while (1) loop
            } // if (k == missingshs->objects)
          } // if (testflag)
	    } 
      } // if (pa != dummypoint)
      // Go to the next crossing tet.
      fnextself(spintet);
      if (spintet.tet == searchtet->tet) break;
    } // while (1)
  } // i

  // Unmark all marked edges.
  for (i = 0; i < crossedges->objects; i++) {
    searchtet = (triface *) fastlookup(crossedges, i);
    unmarkedge(*searchtet);
  }
  crossedges->restart();


  if (invalidflag) {
    // Unmark all collected tets.
    for (i = 0; i < crosstets->objects; i++) {
      searchtet = (triface *) fastlookup(crosstets, i);
      uninfect(*searchtet);
    }
    // Unmark all collected vertices.
    for (i = 0; i < botpoints->objects; i++) {
      parypt = (point *) fastlookup(botpoints, i);
      puninfect(*parypt);
    }
    for (i = 0; i < toppoints->objects; i++) {
      parypt = (point *) fastlookup(toppoints, i);
      puninfect(*parypt);
    }
    crosstets->restart();
    botpoints->restart();
    toppoints->restart();

    // Randomly split an interior edge of R.
    i = randomnation(missingshs->objects - 1);
    recentsh = * (face *) fastlookup(missingshs, i);
    return false;
  }

  if (b->verbose > 2) {
    printf("      Formed cavity: %ld (%ld) cross tets (edges).\n", 
           crosstets->objects, crossedges->objects);
  }

  // Collect the top and bottom faces and the middle vertices. Since all top
  //   and bottom vertices have been infected. Uninfected vertices must be
  //   middle vertices (i.e., the vertices of R).
  // NOTE 1: Hull tets may be collected. Process them as a normal one.
  // NOTE 2: Some previously recovered subfaces may be completely inside the
  //   cavity. In such case, we remove these subfaces from the cavity and put
  //   them into 'subfacstack'. They will be recovered later.
  // NOTE 3: Some segments may be completely inside the cavity, e.g., they
  //   attached to a subface which is inside the cavity. Such segments are
  //   put in 'subsegstack'. They will be recovered later. 
  // NOTE4 : The interior subfaces and segments mentioned in NOTE 2 and 3
  //   are identified in the routine "carvecavity()". 

  for (i = 0; i < crosstets->objects; i++) {
    searchtet = (triface *) fastlookup(crosstets, i);
    // searchtet is [d,e,a,b].
    eorgoppo(*searchtet, spintet);
    fsym(spintet, neightet); // neightet is [a,b,e,#]
    if (!infected(neightet)) {
      // A top face.
      topfaces->newindex((void **) &parytet);
      *parytet = neightet;
    }
    edestoppo(*searchtet, spintet);
    fsym(spintet, neightet); // neightet is [b,a,d,#]
    if (!infected(neightet)) {
      // A bottom face.
      botfaces->newindex((void **) &parytet);
      *parytet = neightet;
    }
    // Add middle vertices if there are (skip dummypoint).
    pa = org(neightet);
    if (!pinfected(pa)) {
      if (pa != dummypoint) {
        pinfect(pa);
        botpoints->newindex((void **) &parypt);
        *parypt = pa;
        toppoints->newindex((void **) &parypt);
        *parypt = pa;
      }
    }
    pa = dest(neightet);
    if (!pinfected(pa)) {
      if (pa != dummypoint) {
        pinfect(pa);
        botpoints->newindex((void **) &parypt);
        *parypt = pa;
        toppoints->newindex((void **) &parypt);
        *parypt = pa;
      }
    }
  } // i

  // Uninfect all collected top, bottom, and middle vertices.
  for (i = 0; i < toppoints->objects; i++) {
    parypt = (point *) fastlookup(toppoints, i);
    puninfect(*parypt);
  }
  for (i = 0; i < botpoints->objects; i++) {
    parypt = (point *) fastlookup(botpoints, i);
    puninfect(*parypt);
  }
  cavitycount++;

  return true;
}

//============================================================================//
//                                                                            //
// delaunizecavity()    Fill a cavity by Delaunay tetrahedra.                 //
//                                                                            //
// The cavity C to be tetrahedralized is the top or bottom part of a whole    //
// cavity. 'cavfaces' contains the boundary faces of C. NOTE: faces in 'cav-  //
// faces' do not form a closed polyhedron.  The "open" side are subfaces of   //
// the missing facet. These faces will be recovered later in fillcavity().    //
//                                                                            //
// This routine first constructs the DT of the vertices. Then it identifies   //
// the half boundary faces of the cavity in DT. Possiblely the cavity C will  //
// be enlarged.                                                               //
//                                                                            //
// The DT is returned in 'newtets'.                                           //
//                                                                            //
//============================================================================//

void TetMeshCore::delaunizecavity(arraypool *cavpoints, arraypool *cavfaces, 
                                 arraypool *cavshells, arraypool *newtets, 
                                 arraypool *crosstets, arraypool *misfaces)
{
  triface searchtet, neightet, *parytet, *parytet1;
  face tmpsh, *parysh;
  point pa, pb, pc, pd, pt[3], *parypt;
  insertvertexflags ivf;
  double ori;
  long baknum, bakhullsize;
  int bakchecksubsegflag, bakchecksubfaceflag;
  int t1ver; 
  int i, j;

  if (b->verbose > 2) {
    printf("      Delaunizing cavity: %ld points, %ld faces.\n", 
           cavpoints->objects, cavfaces->objects);
  }
  // Remember the current number of crossing tets. It may be enlarged later.
  baknum = crosstets->objects;
  bakhullsize = hullsize;
  bakchecksubsegflag = checksubsegflag;
  bakchecksubfaceflag = checksubfaceflag;
  hullsize = 0l;
  checksubsegflag = 0;
  checksubfaceflag = 0;
  b->verbose--;  // Suppress informations for creating Delaunay tetra.
  b->plc = 0; // Do not check near vertices.

  ivf.bowywat = 1; // Use Bowyer-Watson algorithm.

  // Get four non-coplanar points (no dummypoint).
  pa = pb = pc = NULL;
  for (i = 0; i < cavfaces->objects; i++) {
    parytet = (triface *) fastlookup(cavfaces, i);
    parytet->ver = epivot[parytet->ver];
    if (apex(*parytet) != dummypoint) {
      pa = org(*parytet);
      pb = dest(*parytet);
      pc = apex(*parytet);
      break;
    }
  }
  pd = NULL;
  for (; i < cavfaces->objects; i++) {
    parytet = (triface *) fastlookup(cavfaces, i);
    pt[0] = org(*parytet);
    pt[1] = dest(*parytet);
    pt[2] = apex(*parytet);
    for (j = 0; j < 3; j++) {
      if (pt[j] != dummypoint) { // Do not include a hull point.
        ori = orient3d(pa, pb, pc, pt[j]);
        if (ori != 0) {
          pd = pt[j];
          if (ori > 0) {  // Swap pa and pb.
            pt[j] = pa; pa = pb; pb = pt[j]; 
          }
          break;
        }
      }
    }
    if (pd != NULL) break;
  }

  // Create an init DT.
  initialdelaunay(pa, pb, pc, pd);

  // Incrementally insert the vertices (duplicated vertices are ignored).
  for (i = 0; i < cavpoints->objects; i++) {
    pt[0] = * (point *) fastlookup(cavpoints, i);
    searchtet = recenttet;
    ivf.iloc = (int) OUTSIDE;
    insertpoint(pt[0], &searchtet, NULL, NULL, &ivf);
  }

  if (b->verbose > 2) {
    printf("      Identifying %ld boundary faces of the cavity.\n", 
           cavfaces->objects);
  }

  while (1) {

    // Identify boundary faces. Mark interior tets. Save missing faces.
    for (i = 0; i < cavfaces->objects; i++) {
      parytet = (triface *) fastlookup(cavfaces, i);
      // Skip an interior face (due to the enlargement of the cavity).
      if (infected(*parytet)) continue;
      parytet->ver = epivot[parytet->ver];
      pt[0] = org(*parytet);
      pt[1] = dest(*parytet);
      pt[2] = apex(*parytet);
      // Create a temp subface.
      makeshellface(subfaces, &tmpsh);
      setshvertices(tmpsh, pt[0], pt[1], pt[2]);
      // Insert tmpsh in DT.
      searchtet.tet = NULL; 
      if (scoutsubface(&tmpsh, &searchtet, 0)) { // shflag = 0
        // Inserted! 'tmpsh' must face toward the inside of the cavity.
        // Remember the boundary tet (outside the cavity) in tmpsh 
        //   (use the adjacent tet slot). 
        tmpsh.sh[0] = (shellface) encode(*parytet);
        // Save this subface.
        cavshells->newindex((void **) &parysh);
        *parysh = tmpsh;
      } 
      else {
        // This boundary face is missing.
        shellfacedealloc(subfaces, tmpsh.sh);
        // Save this face in list.
        misfaces->newindex((void **) &parytet1);
        *parytet1 = *parytet;
      }
    } // i

    if (misfaces->objects > 0) {
      if (b->verbose > 2) {
        printf("      Enlarging the cavity. %ld missing bdry faces\n", 
               misfaces->objects);
      }

      // Removing all temporary subfaces.
      for (i = 0; i < cavshells->objects; i++) {
        parysh = (face *) fastlookup(cavshells, i);
        stpivot(*parysh, neightet);
        tsdissolve(neightet); // Detach it from adj. tets.
        fsymself(neightet);
        tsdissolve(neightet);
        shellfacedealloc(subfaces, parysh->sh);
      }
      cavshells->restart();

      // Infect the points which are of the cavity.
      for (i = 0; i < cavpoints->objects; i++) {
        pt[0] = * (point *) fastlookup(cavpoints, i);
        pinfect(pt[0]); // Mark it as inserted.
      }

      // Enlarge the cavity.
      for (i = 0; i < misfaces->objects; i++) {
        // Get a missing face.
        parytet = (triface *) fastlookup(misfaces, i);
        if (!infected(*parytet)) {
          // Put it into crossing tet list.
          infect(*parytet);
          crosstets->newindex((void **) &parytet1);
          *parytet1 = *parytet;
          // Insert the opposite point if it is not in DT.
          pd = oppo(*parytet);
          if (!pinfected(pd)) {
            searchtet = recenttet;
            ivf.iloc = (int) OUTSIDE;
            insertpoint(pd, &searchtet, NULL, NULL, &ivf);
            pinfect(pd);
            cavpoints->newindex((void **) &parypt);
            *parypt = pd;
          }
          // Add three opposite faces into the boundary list.
          for (j = 0; j < 3; j++) {
            esym(*parytet, neightet);
            fsymself(neightet);
            if (!infected(neightet)) {
              cavfaces->newindex((void **) &parytet1);
              *parytet1 = neightet;
            } 
            enextself(*parytet);
          } // j
        } // if (!infected(parytet))
      } // i

      // Uninfect the points which are of the cavity.
      for (i = 0; i < cavpoints->objects; i++) {
        pt[0] = * (point *) fastlookup(cavpoints, i);
        puninfect(pt[0]);
      }

      misfaces->restart();
      continue;
    } // if (misfaces->objects > 0)

    break;

  } // while (1)

  // Collect all tets of the DT. All new tets are marktested.
  marktest(recenttet);
  newtets->newindex((void **) &parytet);
  *parytet = recenttet;
  for (i = 0; i < newtets->objects; i++) {
    searchtet = * (triface *) fastlookup(newtets, i);
    for (j = 0; j < 4; j++) {
      decode(searchtet.tet[j], neightet);
      if (!marktested(neightet)) {
        marktest(neightet);
        newtets->newindex((void **) &parytet);
        *parytet = neightet;
      }
    }
  }

  cavpoints->restart();
  cavfaces->restart();

  if (crosstets->objects > baknum) {
    // The cavity has been enlarged.
    cavityexpcount++;
  }

  // Restore the original values.
  hullsize = bakhullsize;
  checksubsegflag = bakchecksubsegflag;
  checksubfaceflag = bakchecksubfaceflag;
  b->verbose++;
  b->plc = 1;
}

//============================================================================//
//                                                                            //
// fillcavity()    Fill new tets into the cavity.                             //
//                                                                            //
// The new tets are stored in two disjoint sets(which share the same facet).  //
// 'topfaces' and 'botfaces' are the boundaries of these two sets, respect-   //
// ively. 'midfaces' is empty on input, and will store faces in the facet.    //
//                                                                            //
// Important: This routine assumes all vertices of the missing region R are   //
// marktested, i.e., pmarktested(p) returns true.                             //
//                                                                            //
//============================================================================//

bool TetMeshCore::fillcavity(arraypool* topshells, arraypool* botshells,
                            arraypool* midfaces, arraypool* missingshs,
                            arraypool* topnewtets, arraypool* botnewtets,
                            triface* crossedge)
{
  arraypool *cavshells;
  triface bdrytet, neightet, *parytet;
  triface searchtet, spintet;
  face *parysh;
  face checkseg;
  point pa, pb, pc;
  bool mflag;
  int t1ver;
  int i, j;

  // Connect newtets to tets outside the cavity.  These connections are needed
  //   for identifying the middle faces (which belong to R).
  for (j = 0; j < 2; j++) {
    cavshells = (j == 0 ? topshells : botshells);
    if (cavshells != NULL) {
      for (i = 0; i < cavshells->objects; i++) {
        // Get a temp subface.
        parysh = (face *) fastlookup(cavshells, i);
        // Get the boundary tet outside the cavity (saved in sh[0]).
        decode(parysh->sh[0], bdrytet);
        pa = org(bdrytet);
        pb = dest(bdrytet);
        pc = apex(bdrytet);
        // Get the adjacent new tet inside the cavity.
        stpivot(*parysh, neightet);
        // Mark neightet as an interior tet of this cavity.
        infect(neightet);
        // Connect the two tets (the old connections are replaced).
        bond(bdrytet, neightet); 
        tsdissolve(neightet); // Clear the pointer to tmpsh.
        // Update the point-to-tets map.
        setpoint2tet(pa, (tetrahedron) neightet.tet);
        setpoint2tet(pb, (tetrahedron) neightet.tet);
        setpoint2tet(pc, (tetrahedron) neightet.tet);
      } // i
    } // if (cavshells != NULL)
  } // j

  if (crossedge != NULL) {
    // Glue top and bottom tets at their common facet.
    triface toptet, bottet, spintet, *midface;
    point pd, pe;
    double ori;
    int types[2], poss[4];
    int interflag;
    int bflag;

    mflag = false;
    pd = org(*crossedge);
    pe = dest(*crossedge);

    // Search the first (middle) face in R. 
    // Since R may be non-convex, we must make sure that the face is in the
    //   interior of R.  We search a face in 'topnewtets' whose three vertices
    //   are on R and it intersects 'crossedge' in its interior. Then search
    //   a matching face in 'botnewtets'.
    for (i = 0; i < topnewtets->objects && !mflag; i++) {
      searchtet = * (triface *) fastlookup(topnewtets, i);
      for (searchtet.ver = 0; searchtet.ver < 4 && !mflag; searchtet.ver++) {
        pa = org(searchtet);
        if (pmarktested(pa)) {
          pb = dest(searchtet);
          if (pmarktested(pb)) {
            pc = apex(searchtet);
            if (pmarktested(pc)) {
              // Check if this face intersects [d,e].
              interflag = tri_edge_test(pa,pb,pc,pd,pe,NULL,1,types,poss);
              if (interflag == 2) {
                // They intersect at a single point. Found.
                toptet = searchtet;
                // The face lies in the interior of R.
                // Get the tet (in topnewtets) which lies above R.
                ori = orient3d(pa, pb, pc, pd);
                if (ori < 0) {
                  fsymself(toptet);
                  pa = org(toptet);
                  pb = dest(toptet);
                } else if (ori == 0) {
                  terminate_tet_core(this, 2);
                }
                // Search the face [b,a,c] in 'botnewtets'.
                for (j = 0; j < botnewtets->objects; j++) {
                  neightet = * (triface *) fastlookup(botnewtets, j);                  
                  // Is neightet contains 'b'.
                  if ((point) neightet.tet[4] == pb) {
                    neightet.ver = 11;
                  } else if ((point) neightet.tet[5] == pb) {
                    neightet.ver = 3;
                  } else if ((point) neightet.tet[6] == pb) {
                    neightet.ver = 7;
                  } else if ((point) neightet.tet[7] == pb) {
                    neightet.ver = 0;
                  } else {
                    continue;
                  }
                  // Is the 'neightet' contains edge [b,a].
                  if (dest(neightet) == pa) {
                    // 'neightet' is just the edge.
                  } else if (apex(neightet) == pa) {
                    eprevesymself(neightet);
                  } else if (oppo(neightet) == pa) {
                    esymself(neightet);
                    enextself(neightet);
                  } else {
                    continue;
                  }
                  // Is 'neightet' the face [b,a,c]. 
                  if (apex(neightet) == pc) {
                    bottet = neightet;
                    mflag = true;
                    break;
                  }
                } // j
              } // if (interflag == 2)
            } // pc
          } // pb
        } // pa
      } // toptet.ver
    } // i

    if (mflag) {
      // Found a pair of matched faces in 'toptet' and 'bottet'.
      bond(toptet, bottet);
      // Both are interior tets.
      infect(toptet);
      infect(bottet);
      // Add this face into search list.
      markface(toptet);
      midfaces->newindex((void **) &parytet);
      *parytet = toptet;
    } else {
      // No pair of 'toptet' and 'bottet'.
      toptet.tet = NULL;
      // Randomly split an interior edge of R.
      i = randomnation(missingshs->objects - 1);
      recentsh = * (face *) fastlookup(missingshs, i);
    }

    // Find other middle faces, connect top and bottom tets.
    for (i = 0; i < midfaces->objects && mflag; i++) {
      // Get a matched middle face [a, b, c]
      midface = (triface *) fastlookup(midfaces, i);
      // Check the neighbors at the edges of this face. 
      for (j = 0; j < 3 && mflag; j++) {
        toptet = *midface;
        bflag = false;
        while (1) {
          // Go to the next face in the same tet.
          esymself(toptet);
          pc = apex(toptet);
          if (pmarktested(pc)) {
            break; // Find a subface.
          }
          if (pc == dummypoint) {
            terminate_tet_core(this, 2); // Check this case.
            break; // Find a subface.
          }
          // Go to the adjacent tet.
          fsymself(toptet);
          // Do we walk outside the cavity? 
          if (!marktested(toptet)) {
            // Yes, the adjacent face is not a middle face.
            bflag = true; break; 
          }
        }
        if (!bflag) {
          if (!facemarked(toptet)) {
            fsym(*midface, bottet);
            spintet = bottet;
            while (1) {
              esymself(bottet);
              pd = apex(bottet);
              if (pd == pc) break; // Face matched.
              fsymself(bottet);
              if (bottet.tet == spintet.tet) {
                // Not found a matched bottom face.
                mflag = false;
                break;
              }
            } // while (1)
            if (mflag) {
              if (marktested(bottet)) {
                // Connect two tets together.
                bond(toptet, bottet);
                // Both are interior tets.
                infect(toptet);
                infect(bottet);
                // Add this face into list.
                markface(toptet);
                midfaces->newindex((void **) &parytet);
                *parytet = toptet;
              }
              else {
                // The 'bottet' is not inside the cavity! 
                terminate_tet_core(this, 2); // Check this case
              }
            } else { // mflag == false
              // Adjust 'toptet' and 'bottet' to be the crossing edges.
              fsym(*midface, bottet);
              spintet = bottet;
              while (1) {
                esymself(bottet);
                pd = apex(bottet);
                if (pmarktested(pd)) {
                  // assert(pd != pc);
                  // Let 'toptet' be [a,b,c,#], and 'bottet' be [b,a,d,*].
                  // Adjust 'toptet' and 'bottet' to be the crossing edges.
                  // Test orient3d(b,c,#,d).
                  ori = orient3d(dest(toptet), pc, oppo(toptet), pd);
                  if (ori < 0) {
                    // Edges [a,d] and [b,c] cross each other.
                    enextself(toptet); // [b,c]
                    enextself(bottet); // [a,d]
                  } else if (ori > 0) {
                    // Edges [a,c] and [b,d] cross each other. 
                    eprevself(toptet); // [c,a]
                    eprevself(bottet); // [d,b]
                  } else {
                    // b,c,#,and d are coplanar!.
                    terminate_tet_core(this, 2); //assert(0);
                  }
                  break; // Not matched
                }
                fsymself(bottet);
              }
            } // if (!mflag)
          } // if (!facemarked(toptet))
        } // if (!bflag)
        enextself(*midface);
      } // j
    } // i

    if (mflag) {
      if (b->verbose > 2) {
        printf("      Found %ld middle subfaces.\n", midfaces->objects);
      }
      face oldsh, newsh, casout, casin, neighsh;

      oldsh = * (face *) fastlookup(missingshs, 0);

      // Create new subfaces to fill the region R.
      for (i = 0; i < midfaces->objects; i++) {
        // Get a matched middle face [a, b, c]
        midface = (triface *) fastlookup(midfaces, i);
        unmarkface(*midface);
        makeshellface(subfaces, &newsh);
        setsorg(newsh, org(*midface));
        setsdest(newsh, dest(*midface));
        setsapex(newsh, apex(*midface));
        // The new subface gets its markers from the old one.
        setshellmark(newsh, shellmark(oldsh));
        if (checkconstraints) {
          setareabound(newsh, areabound(oldsh));
        }
        if (useinsertradius) {
          setfacetindex(newsh, getfacetindex(oldsh));
        }
        // Connect the new subface to adjacent tets.
        tsbond(*midface, newsh);
        fsym(*midface, neightet);
        sesymself(newsh);
        tsbond(neightet, newsh);
      }

      // Connect new subfaces together and to the bdry of R.
      // Delete faked segments.
      for (i = 0; i < midfaces->objects; i++) {
        // Get a matched middle face [a, b, c]
        midface = (triface *) fastlookup(midfaces, i);
        for (j = 0; j < 3; j++) {
          tspivot(*midface, newsh);
          spivot(newsh, casout);
          if (casout.sh == NULL) {
            // Search its neighbor.
            fnext(*midface, searchtet);
            while (1) {
              // (1) First check if this side is a bdry edge of R.
              tsspivot1(searchtet, checkseg);
              if (checkseg.sh != NULL) {
                // It's a bdry edge of R.
                // Get the old subface.
                checkseg.shver = 0;
                spivot(checkseg, oldsh);
                if (sinfected(checkseg)) {
                  // It's a faked segment. Delete it.
                  spintet = searchtet;
                  while (1) {
                    tssdissolve1(spintet);
                    fnextself(spintet);
                    if (spintet.tet == searchtet.tet) break;
                  }
                  shellfacedealloc(subsegs, checkseg.sh);
                  ssdissolve(oldsh);
                  checkseg.sh = NULL;
                }
                spivot(oldsh, casout);
                if (casout.sh != NULL) {
                  casin = casout;
                  if (checkseg.sh != NULL) {
                    // Make sure that the subface has the right ori at the 
                    //   segment.
                    checkseg.shver = 0;
                    if (sorg(newsh) != sorg(checkseg)) {
                      sesymself(newsh);
                    }
                    spivot(casin, neighsh);
                    while (neighsh.sh != oldsh.sh) {
                      casin = neighsh;
                      spivot(casin, neighsh);
                    }
                  }
                  sbond1(newsh, casout);
                  sbond1(casin, newsh);
                }
                if (checkseg.sh != NULL) {
                  ssbond(newsh, checkseg);
                }
                break;
              } // if (checkseg.sh != NULL)
              // (2) Second check if this side is an interior edge of R.
              tspivot(searchtet, neighsh);
              if (neighsh.sh != NULL) {
                // Found an adjacent subface of newsh (an interior edge).
                sbond(newsh, neighsh);
                break;
              }
              fnextself(searchtet);
            } // while (1)
          } // if (casout.sh == NULL)
          enextself(*midface);
        } // j
      } // i

      // Delete old subfaces.
      for (i = 0; i < missingshs->objects; i++) {
        parysh = (face *) fastlookup(missingshs, i);
        shellfacedealloc(subfaces, parysh->sh);
      }
    } else {
      if (toptet.tet != NULL) {
        // Faces at top and bottom are not matched. 
        // Choose a Steiner point in R.
        // Split one of the crossing edges.
        pa = org(toptet);
        pb = dest(toptet);
        pc = org(bottet);
        pd = dest(bottet);
        // Search an edge in R which is either [a,b] or [c,d].
        // Reminder:  Subfaces in this list 'missingshs', except the first
        //   one, represents an interior edge of R. 
        parysh = NULL;  // Avoid a warning in MSVC
        for (i = 1; i < missingshs->objects; i++) {
          parysh = (face *) fastlookup(missingshs, i);
          if (((sorg(*parysh) == pa) && (sdest(*parysh) == pb)) ||
              ((sorg(*parysh) == pb) && (sdest(*parysh) == pa))) break;
          if (((sorg(*parysh) == pc) && (sdest(*parysh) == pd)) ||
              ((sorg(*parysh) == pd) && (sdest(*parysh) == pc))) break;
        }
        if (i < missingshs->objects) {
          // Found. Return it.
          recentsh = *parysh;
        } else {
          terminate_tet_core(this, 2); //assert(0);
        }
      } else {
        //terminate_tet_core(this, 2); // Report a bug
      }
    }

    midfaces->restart();
  } else {
    mflag = true;
  } 

  // Delete the temp subfaces.
  for (j = 0; j < 2; j++) {
    cavshells = (j == 0 ? topshells : botshells);
    if (cavshells != NULL) {
      for (i = 0; i < cavshells->objects; i++) {
        parysh = (face *) fastlookup(cavshells, i);
        shellfacedealloc(subfaces, parysh->sh);
      }
    }
  }

  topshells->restart();
  if (botshells != NULL) {
    botshells->restart();
  }

  return mflag;
}

//============================================================================//
//                                                                            //
// carvecavity()    Delete old tets and outer new tets of the cavity.         //
//                                                                            //
//============================================================================//

void TetMeshCore::carvecavity(arraypool *crosstets, arraypool *topnewtets,
                             arraypool *botnewtets)
{
  arraypool *newtets;
  shellface *sptr, *ssptr;
  triface *parytet, *pnewtet, newtet, neightet, spintet;
  face checksh, *parysh;
  face checkseg, *paryseg;
  int t1ver;
  int i, j;

  if (b->verbose > 2) {
    printf("      Carve cavity: %ld old tets.\n", crosstets->objects);
  }

  // First process subfaces and segments which are adjacent to the cavity.
  //   They must be re-connected to new tets in the cavity.
  // Comment: It is possible that some subfaces and segments are completely
  //   inside the cavity. This can happen even if the cavity is not enlarged. 
  //   Before deleting the old tets, find and queue all interior subfaces
  //   and segments. They will be recovered later. 2010-05-06.

  // Collect all subfaces and segments which attached to the old tets.
  for (i = 0; i < crosstets->objects; i++) {
    parytet = (triface *) fastlookup(crosstets, i);
    if ((sptr = (shellface*) parytet->tet[9]) != NULL) {
      for (j = 0; j < 4; j++) {
        if (sptr[j]) {
          sdecode(sptr[j], checksh);
          if (!sinfected(checksh)) {
            sinfect(checksh);
            cavetetshlist->newindex((void **) &parysh);
            *parysh = checksh;
          }
        }
      } // j
    }
    if ((ssptr = (shellface*) parytet->tet[8]) != NULL) {
      for (j = 0; j < 6; j++) {
        if (ssptr[j]) {
          sdecode(ssptr[j], checkseg);
          // Skip a deleted segment (was a faked segment)
          if (checkseg.sh[3] != NULL) {
            if (!sinfected(checkseg)) {
              sinfect(checkseg);
              cavetetseglist->newindex((void **) &paryseg);
              *paryseg = checkseg;
            }
          }
        }
      } // j
    }
  } // i

  // Uninfect collected subfaces.
  for (i = 0; i < cavetetshlist->objects; i++) {
    parysh = (face *) fastlookup(cavetetshlist, i);
    suninfect(*parysh);
  }
  // Uninfect collected segments.
  for (i = 0; i < cavetetseglist->objects; i++) {
    paryseg = (face *) fastlookup(cavetetseglist, i);
    suninfect(*paryseg);
  }

  // Connect subfaces to new tets.
  for (i = 0; i < cavetetshlist->objects; i++) {
    parysh = (face *) fastlookup(cavetetshlist, i);
    // Get an adjacent tet at this subface.
    stpivot(*parysh, neightet);
    // Does this tet lie inside the cavity.
    if (infected(neightet)) {
      // Yes. Get the other adjacent tet at this subface.
      sesymself(*parysh);
      stpivot(*parysh, neightet);
      // Does this tet lie inside the cavity.
      if (infected(neightet)) {
        checksh = *parysh;
        stdissolve(checksh);
        caveencshlist->newindex((void **) &parysh);
        *parysh = checksh;
      }
    }
    if (!infected(neightet)) {
      // Found an outside tet. Re-connect this subface to a new tet.
      fsym(neightet, newtet);
      sesymself(*parysh);
      tsbond(newtet, *parysh);
    }
  } // i


  for (i = 0; i < cavetetseglist->objects; i++) {
    checkseg = * (face *) fastlookup(cavetetseglist, i);
    // Check if the segment is inside the cavity.
    sstpivot1(checkseg, neightet);
    spintet = neightet;
    while (1) {
      if (!infected(spintet)) {
        // This segment is on the boundary of the cavity.
        break;
      }
      fnextself(spintet);
      if (spintet.tet == neightet.tet) {
        sstdissolve1(checkseg);
        caveencseglist->newindex((void **) &paryseg);
        *paryseg = checkseg;
        break;
      }
    }
    if (!infected(spintet)) {
      // A boundary segment. Connect this segment to the new tets.
      sstbond1(checkseg, spintet);
      neightet = spintet;
      while (1) {
        tssbond1(spintet, checkseg);
        fnextself(spintet);
        if (spintet.tet == neightet.tet) break;
      }
    }
  } // i


  cavetetshlist->restart();
  cavetetseglist->restart();

  // Delete the old tets in cavity.
  for (i = 0; i < crosstets->objects; i++) {
    parytet = (triface *) fastlookup(crosstets, i);
    if (ishulltet(*parytet)) {
      hullsize--;
    }
    tetrahedrondealloc(parytet->tet);
  }

  crosstets->restart(); // crosstets will be re-used.

  // Collect new tets in cavity.  Some new tets have already been found 
  //   (and infected) in the fillcavity(). We first collect them.
  for (j = 0; j < 2; j++) {
    newtets = (j == 0 ? topnewtets : botnewtets);
    if (newtets != NULL) {
      for (i = 0; i < newtets->objects; i++) {
        parytet = (triface *) fastlookup(newtets, i);
        if (infected(*parytet)) {
          crosstets->newindex((void **) &pnewtet);
          *pnewtet = *parytet;
        }
      } // i
    }
  } // j

  // Now we collect all new tets in cavity.
  for (i = 0; i < crosstets->objects; i++) {
    parytet = (triface *) fastlookup(crosstets, i);
    for (j = 0; j < 4; j++) {
      decode(parytet->tet[j], neightet);
      if (marktested(neightet)) { // Is it a new tet?
        if (!infected(neightet)) {
          // Find an interior tet.
          //assert((point) neightet.tet[7] != dummypoint); // SELF_CHECK
          infect(neightet);
          crosstets->newindex((void **) &pnewtet);
          *pnewtet = neightet;
        }
      }
    } // j
  } // i

  parytet = (triface *) fastlookup(crosstets, 0);
  recenttet = *parytet; // Remember a live handle.

  // Delete outer new tets.
  for (j = 0; j < 2; j++) {
    newtets = (j == 0 ? topnewtets : botnewtets);
    if (newtets != NULL) {
      for (i = 0; i < newtets->objects; i++) {
        parytet = (triface *) fastlookup(newtets, i);
        if (infected(*parytet)) {
          // This is an interior tet.
          uninfect(*parytet);
          unmarktest(*parytet);
          if (ishulltet(*parytet)) {
            hullsize++;
          }
        } else {
          // An outer tet. Delete it.
          tetrahedrondealloc(parytet->tet);
        }
      }
    }
  }

  crosstets->restart();
  topnewtets->restart();
  if (botnewtets != NULL) {
    botnewtets->restart();
  }
}

//============================================================================//
//                                                                            //
// restorecavity()    Reconnect old tets and delete new tets of the cavity.   //
//                                                                            //
//============================================================================//

void TetMeshCore::restorecavity(arraypool *crosstets, arraypool *topnewtets,
                               arraypool *botnewtets, arraypool *missingshbds)
{
  triface *parytet, neightet, spintet;
  face *parysh;
  face checkseg;
  point *ppt;
  int t1ver;
  int i, j;

  // Reconnect crossing tets to cavity boundary.
  for (i = 0; i < crosstets->objects; i++) {
    parytet = (triface *) fastlookup(crosstets, i);
    parytet->ver = 0;
    for (parytet->ver = 0; parytet->ver < 4; parytet->ver++) {
      fsym(*parytet, neightet);
      if (!infected(neightet)) {
        // Restore the old connections of tets.
        bond(*parytet, neightet);
      }
    }
    // Update the point-to-tet map.
    parytet->ver = 0;
    ppt = (point *) &(parytet->tet[4]);
    for (j = 0; j < 4; j++) {
      setpoint2tet(ppt[j], encode(*parytet));
    }
  }

  // Uninfect all crossing tets.
  for (i = 0; i < crosstets->objects; i++) {
    parytet = (triface *) fastlookup(crosstets, i);
    uninfect(*parytet);
  }

  // Remember a live handle.
  if (crosstets->objects > 0) {
    recenttet = * (triface *) fastlookup(crosstets, 0);
  }

  // Delete faked segments.
  for (i = 0; i < missingshbds->objects; i++) {
    parysh = (face *) fastlookup(missingshbds, i);
    sspivot(*parysh, checkseg);
    if (checkseg.sh[3] != NULL) {
      if (sinfected(checkseg)) {
            // It's a faked segment. Delete it.
        sstpivot1(checkseg, neightet);
        spintet = neightet;
        while (1) {
          tssdissolve1(spintet);
          fnextself(spintet);
          if (spintet.tet == neightet.tet) break;
        }
        shellfacedealloc(subsegs, checkseg.sh);
        ssdissolve(*parysh);
        //checkseg.sh = NULL;
      }
    }
  } // i

  // Delete new tets.
  for (i = 0; i < topnewtets->objects; i++) {
    parytet = (triface *) fastlookup(topnewtets, i);
    tetrahedrondealloc(parytet->tet);
  }

  if (botnewtets != NULL) {
    for (i = 0; i < botnewtets->objects; i++) {
      parytet = (triface *) fastlookup(botnewtets, i);
      tetrahedrondealloc(parytet->tet);
    }
  }

  crosstets->restart();
  topnewtets->restart();
  if (botnewtets != NULL) {
    botnewtets->restart();
  }
}

//============================================================================//
//                                                                            //
// flipcertify()    Insert a crossing face into priority queue.               //
//                                                                            //
// A crossing face of a facet must have at least one top and one bottom ver-  //
// tex of the facet.                                                          //
//                                                                            //
//============================================================================//

void TetMeshCore::flipcertify(triface *chkface,badface **pqueue,point plane_pa,
                             point plane_pb, point plane_pc)
{
  badface *parybf, *prevbf, *nextbf;
  triface neightet;
  face checksh;
  point p[5];
  double w[5];
  double insph, ori4;
  int topi, boti;
  int i;

  // Compute the flip time \tau.
  fsym(*chkface, neightet);

  p[0] = org(*chkface);
  p[1] = dest(*chkface);
  p[2] = apex(*chkface);
  p[3] = oppo(*chkface);
  p[4] = oppo(neightet);

  // Check if the face is a crossing face.
  topi = boti = 0;
  for (i = 0; i < 3; i++) {
    if (pmarktest2ed(p[i])) topi++;
    if (pmarktest3ed(p[i])) boti++;
  }
  if ((topi == 0) || (boti == 0)) {
    // It is not a crossing face.
    // return;
    for (i = 3; i < 5; i++) {
      if (pmarktest2ed(p[i])) topi++;
      if (pmarktest3ed(p[i])) boti++;
    }
    if ((topi == 0) || (boti == 0)) {
      // The two tets sharing at this face are on one side of the facet.
      // Check if this face is locally Delaunay (due to rounding error).
      if ((p[3] != dummypoint) && (p[4] != dummypoint)) {
        // Do not check it if it is a subface.
        tspivot(*chkface, checksh);
        if (checksh.sh == NULL) {
          insph = insphere_s(p[1], p[0], p[2], p[3], p[4]);
          if (insph > 0) {
            // Add the face into queue.
            if (b->verbose > 2) {
              printf("      A locally non-Delanay face (%d, %d, %d)-%d,%d\n", 
                     pointmark(p[0]), pointmark(p[1]), pointmark(p[2]), 
                     pointmark(p[3]), pointmark(p[4]));
            }
            parybf = (badface *) flippool->alloc();
            parybf->key = 0.;  // tau = 0, do immediately.
            parybf->tt = *chkface;
            parybf->forg = p[0];
            parybf->fdest = p[1];
            parybf->fapex = p[2];
            parybf->foppo = p[3];
            parybf->noppo = p[4];
            // Add it at the top of the priority queue.
            if (*pqueue == NULL) {
              *pqueue = parybf;
              parybf->nextitem = NULL;
            } else {
              parybf->nextitem = *pqueue;
              *pqueue = parybf;
            }
          } // if (insph > 0)
        } // if (checksh.sh == NULL)
      }
    }
    return; // Test: omit this face.
  }

  // Decide the "height" for each point.
  for (i = 0; i < 5; i++) {
    if (pmarktest2ed(p[i])) {
      // A top point has a positive weight.
      w[i] = orient3dfast(plane_pa, plane_pb, plane_pc, p[i]);      
      if (w[i] < 0) w[i] = -w[i];
    } else {
      w[i] = 0;
    }
  }

  insph = insphere(p[1], p[0], p[2], p[3], p[4]);
  ori4 = orient4d(p[1], p[0], p[2], p[3], p[4], w[1], w[0], w[2], w[3], w[4]);
  if (ori4 > 0) {
    // Add the face into queue.
    if (b->verbose > 2) {
      printf("      Insert face (%d, %d, %d) - %d, %d\n", pointmark(p[0]),
        pointmark(p[1]), pointmark(p[2]), pointmark(p[3]), pointmark(p[4]));
    }
    
    parybf = (badface *) flippool->alloc();

    parybf->key = -insph / ori4;
    parybf->tt = *chkface;
    parybf->forg = p[0];
    parybf->fdest = p[1];
    parybf->fapex = p[2];
    parybf->foppo = p[3];
    parybf->noppo = p[4];

    // Push the face into priority queue.
    //pq.push(bface);
    if (*pqueue == NULL) {
      *pqueue = parybf;
      parybf->nextitem = NULL;
    } else {
      // Search an item whose key is larger or equal to current key.
      prevbf = NULL;
      nextbf = *pqueue;
      //if (!b->flipinsert_random) { // Default use a priority queue.
        // Insert the item into priority queue.
        while (nextbf != NULL) {
          if (nextbf->key < parybf->key) {
            prevbf = nextbf;
            nextbf = nextbf->nextitem;
          } else {
            break;
          }
        }
      //} // if (!b->flipinsert_random)
      // Insert the new item between prev and next items.
      if (prevbf == NULL) {
        *pqueue = parybf;
      } else {
        prevbf->nextitem = parybf;
      }
      parybf->nextitem = nextbf;
    }
  } else if (ori4 == 0) {
    
  }
}

//============================================================================//
//                                                                            //
// flipinsertfacet()    Insert a facet into a CDT by flips.                   //
//                                                                            //
// The algorithm is described in Shewchuk's paper "Updating and Constructing  //
// Constrained Delaunay and Constrained Regular Triangulations by Flips", in  //
// Proc. 19th Ann. Symp. on Comput. Geom., 86--95, 2003.                      //
//                                                                            //
// 'crosstets' contains the set of crossing tetrahedra (infected) of the      //
// facet.  'toppoints' and 'botpoints' are points lies above and below the    //
// facet, not on the facet.                                                   //
//                                                                            //
//============================================================================//

void TetMeshCore::flipinsertfacet(arraypool *crosstets, arraypool *toppoints,
                                 arraypool *botpoints, arraypool *midpoints)
{
  arraypool *crossfaces, *bfacearray;
  triface fliptets[6], baktets[2], fliptet, newface;
  triface neightet, *parytet;
  badface *pqueue;
  badface *popbf, bface;
  point plane_pa, plane_pb, plane_pc;
  point p1, p2, pd, pe;
  point *parypt;
  flipconstraints fc;
  double ori[3];
  int convcount, copcount;
  int flipflag, fcount;
  int n, i;
  long f23count, f32count, f44count;
  long totalfcount;

  f23count = flip23count;
  f32count = flip32count;
  f44count = flip44count;

  // Get three affinely independent vertices in the missing region R.
  calculateabovepoint(midpoints, &plane_pa, &plane_pb, &plane_pc);

  // Mark top and bottom points. Do not mark midpoints.
  for (i = 0; i < toppoints->objects; i++) {
    parypt = (point *) fastlookup(toppoints, i);
    if (!pmarktested(*parypt)) {
      pmarktest2(*parypt);
    }
  }
  for (i = 0; i < botpoints->objects; i++) {
    parypt = (point *) fastlookup(botpoints, i);
    if (!pmarktested(*parypt)) {
      pmarktest3(*parypt);
    }
  }

  // Collect crossing faces. 
  crossfaces = cavetetlist;  // Re-use array 'cavetetlist'.

  // Each crossing face contains at least one bottom vertex and
  //   one top vertex.
  for (i = 0; i < crosstets->objects; i++) {
    parytet = (triface *) fastlookup(crosstets, i);
    fliptet = *parytet;
    for (fliptet.ver = 0; fliptet.ver < 4; fliptet.ver++) {
      fsym(fliptet, neightet);
      if (infected(neightet)) { // It is an interior face.
        if (!marktested(neightet)) { // It is an unprocessed face.
          crossfaces->newindex((void **) &parytet);
          *parytet = fliptet;
        }
      }
    }
    marktest(fliptet);
  }

  if (b->verbose > 1) {
    printf("    Found %ld crossing faces.\n", crossfaces->objects);
  }

  for (i = 0; i < crosstets->objects; i++) {
    parytet = (triface *) fastlookup(crosstets, i);
    unmarktest(*parytet);
    uninfect(*parytet);
  }

  // Initialize the priority queue.
  pqueue = NULL;

  for (i = 0; i < crossfaces->objects; i++) {
    parytet = (triface *) fastlookup(crossfaces, i);
    flipcertify(parytet, &pqueue, plane_pa, plane_pb, plane_pc);
  }
  crossfaces->restart();

  // The list for temporarily storing unflipable faces.
  bfacearray = new arraypool(sizeof(triface), 4);


  fcount = 0;  // Count the number of flips.

  // Flip insert the facet.
  while (pqueue != NULL) {

    // Pop a face from the priority queue.
    popbf = pqueue;
    bface = *popbf;
    // Update the queue.
    pqueue = pqueue->nextitem;
    // Delete the popped item from the pool.
    flippool->dealloc((void *) popbf);

    if (!isdeadtet(bface.tt)) {
      if ((org(bface.tt) == bface.forg) && (dest(bface.tt) == bface.fdest) &&
          (apex(bface.tt) == bface.fapex) && (oppo(bface.tt) == bface.foppo)) {
        // It is still a crossing face of R.
        fliptet = bface.tt;
        fsym(fliptet, neightet);
        if (oppo(neightet) == bface.noppo) {
          pd = oppo(fliptet);
          pe = oppo(neightet);

          if (b->verbose > 2) {
            printf("      Get face (%d, %d, %d) - %d, %d, tau = %.17g\n",
                   pointmark(bface.forg), pointmark(bface.fdest),
                   pointmark(bface.fapex), pointmark(bface.foppo),
                   pointmark(bface.noppo), bface.key);
          }
          flipflag = 0;

          // Check for which type of flip can we do.
          convcount = 3;
          copcount = 0;
          for (i = 0; i < 3; i++) {
            p1 = org(fliptet);
            p2 = dest(fliptet);
            ori[i] = orient3d(p1, p2, pd, pe);
            if (ori[i] < 0) {
              convcount--;
              //break;
            } else if (ori[i] == 0) {
              convcount--; // Possible 4-to-4 flip.
              copcount++;
              //break;
            }
            enextself(fliptet);
          }

          if (convcount == 3) {
            // A 2-to-3 flip is found.
            fliptets[0] = fliptet; // abcd, d may be the new vertex.
            fliptets[1] = neightet; // bace.
            flip23(fliptets, 1, &fc);
            // Put the link faces into check list.
            for (i = 0; i < 3; i++) {
              eprevesym(fliptets[i], newface);
              crossfaces->newindex((void **) &parytet);
              *parytet = newface;
            }
            for (i = 0; i < 3; i++) {
              enextesym(fliptets[i], newface);
              crossfaces->newindex((void **) &parytet);
              *parytet = newface;
            }
            flipflag = 1;
          } else if (convcount == 2) {
            //if (copcount <= 1) {
            // A 3-to-2 or 4-to-4 may be possible.
            // Get the edge which is locally non-convex or flat. 
            for (i = 0; i < 3; i++) {
              if (ori[i] <= 0) break;
              enextself(fliptet);
            }

            // Collect tets sharing at this edge.
            esym(fliptet, fliptets[0]); // [b,a,d,c]
            n = 0;
            do {
              p1 = apex(fliptets[n]);
              if (!(pmarktested(p1) || pmarktest2ed(p1) || pmarktest3ed(p1))) {
                // This apex is not on the cavity. Hence the face does not
                //   lie inside the cavity. Do not flip this edge.
                n = 1000; break;
              }
              fnext(fliptets[n], fliptets[n + 1]);
              n++;
            } while ((fliptets[n].tet != fliptet.tet) && (n < 5));

            if (n == 3) {
              // Found a 3-to-2 flip.
              flip32(fliptets, 1, &fc);
              // Put the link faces into check list.
              for (i = 0; i < 3; i++) {
                esym(fliptets[0], newface);
                crossfaces->newindex((void **) &parytet);
                *parytet = newface;
                enextself(fliptets[0]);
              }
              for (i = 0; i < 3; i++) {
                esym(fliptets[1], newface);
                crossfaces->newindex((void **) &parytet);
                *parytet = newface;
                enextself(fliptets[1]);
              }
              flipflag = 1;
            } else if (n == 4) {
              if (copcount == 1) {                
                // Found a 4-to-4 flip. 
                // Let the six vertices are: a,b,c,d,e,f, where
                //   fliptets[0] = [b,a,d,c]
                //           [1] = [b,a,c,e]
                //           [2] = [b,a,e,f]
                //           [3] = [b,a,f,d]
                // After the 4-to-4 flip, edge [a,b] is flipped, edge [e,d]
                //   is created.
                // First do a 2-to-3 flip.
                // Comment: This flip temporarily creates a degenerated
                //   tet (whose volume is zero). It will be removed by the 
                //   followed 3-to-2 flip.
                fliptets[0] = fliptet; // = [a,b,c,d], d is the new vertex.
                // fliptets[1];        // = [b,a,c,e].
                baktets[0] = fliptets[2]; // = [b,a,e,f]
                baktets[1] = fliptets[3]; // = [b,a,f,d]
                // The flip may involve hull tets.
                flip23(fliptets, 1, &fc);
                // Put the "outer" link faces into check list.
                //   fliptets[0] = [e,d,a,b] => will be flipped, so 
                //   [a,b,d] and [a,b,e] are not "outer" link faces.
                for (i = 1; i < 3; i++) {
                  eprevesym(fliptets[i], newface);
                  crossfaces->newindex((void **) &parytet);
                  *parytet = newface;
                }
                for (i = 1; i < 3; i++) {
                  enextesym(fliptets[i], newface);
                  crossfaces->newindex((void **) &parytet);
                  *parytet = newface;
                }
                // Then do a 3-to-2 flip.
                enextesymself(fliptets[0]);  // fliptets[0] is [e,d,a,b].
                eprevself(fliptets[0]); // = [b,a,d,c], d is the new vertex.
                fliptets[1] = baktets[0]; // = [b,a,e,f]
                fliptets[2] = baktets[1]; // = [b,a,f,d]
                flip32(fliptets, 1, &fc);
                // Put the "outer" link faces into check list.
                //   fliptets[0] = [d,e,f,a]
                //   fliptets[1] = [e,d,f,b]
                //   Faces [a,b,d] and [a,b,e] are not "outer" link faces.
                enextself(fliptets[0]);
                for (i = 1; i < 3; i++) {
                  esym(fliptets[0], newface);
                  crossfaces->newindex((void **) &parytet);
                  *parytet = newface;
                  enextself(fliptets[0]);
                }
                enextself(fliptets[1]);
                for (i = 1; i < 3; i++) {
                  esym(fliptets[1], newface);
                  crossfaces->newindex((void **) &parytet);
                  *parytet = newface;
                  enextself(fliptets[1]);
                }
                flip23count--;
                flip32count--;
                flip44count++;
                flipflag = 1;
              }
            }
          } else {
            // There are more than 1 non-convex or coplanar cases.
            flipflag = -1; // Ignore this face.
            if (b->verbose > 2) {
              printf("        Ignore face (%d, %d, %d) - %d, %d, tau = %.17g\n",
                     pointmark(bface.forg), pointmark(bface.fdest),
                     pointmark(bface.fapex), pointmark(bface.foppo),
                     pointmark(bface.noppo), bface.key);
            }
          } // if (convcount == 1)

          if (flipflag == 1) {
            // Update the priority queue.
            for (i = 0; i < crossfaces->objects; i++) {
              parytet = (triface *) fastlookup(crossfaces, i);
              flipcertify(parytet, &pqueue, plane_pa, plane_pb, plane_pc);
            }
            crossfaces->restart();
            if (1) { // if (!b->flipinsert_random) {
              // Insert all queued unflipped faces.
              for (i = 0; i < bfacearray->objects; i++) {
                parytet = (triface *) fastlookup(bfacearray, i);
                // This face may be changed.
                if (!isdeadtet(*parytet)) {
                  flipcertify(parytet, &pqueue, plane_pa, plane_pb, plane_pc);
                }
              }
              bfacearray->restart();
            }
            fcount++;
          } else if (flipflag == 0) {
            // Queue an unflippable face. To process it later.
            bfacearray->newindex((void **) &parytet);
            *parytet = fliptet;
          }
        } // if (pe == bface.noppo)  
      } // if ((pa == bface.forg) && ...)
    } // if (bface.tt != NULL)

  } // while (pqueue != NULL)

  if (bfacearray->objects > 0) {
    if (fcount == 0) {
      printf("!! No flip is found in %ld faces.\n", bfacearray->objects);
      terminate_tet_core(this, 2); //assert(0);
    }
  }

  delete bfacearray;

  // Un-mark top and bottom points.
  for (i = 0; i < toppoints->objects; i++) {
    parypt = (point *) fastlookup(toppoints, i);
    punmarktest2(*parypt);
  }
  for (i = 0; i < botpoints->objects; i++) {
    parypt = (point *) fastlookup(botpoints, i);
    punmarktest3(*parypt);
  }

  f23count = flip23count - f23count;
  f32count = flip32count - f32count;
  f44count = flip44count - f44count;
  totalfcount = f23count + f32count + f44count;
  if (b->verbose > 2) {
    printf("      Total %ld flips. f23(%ld), f32(%ld), f44(%ld).\n",
           totalfcount, f23count, f32count, f44count);
  }
}

//============================================================================//
//                                                                            //
// insertpoint_cdt()    Insert a new point into a CDT.                        //
//                                                                            //
//============================================================================//

int TetMeshCore::insertpoint_cdt(point newpt, triface *searchtet, face *splitsh, 
                                face *splitseg, insertvertexflags *ivf,
                                arraypool *cavpoints, arraypool *cavfaces,
                                arraypool *cavshells, arraypool *newtets,
                                arraypool *crosstets, arraypool *misfaces)
{
  triface neightet, *parytet;
  face checksh, *parysh, *parysh1;
  face *paryseg, *paryseg1;
  point *parypt;
  int t1ver;
  int i;

  if (b->verbose > 2) {
    printf("      Insert point %d into CDT\n", pointmark(newpt));
  }

  if (!insertpoint(newpt, searchtet, NULL, NULL, ivf)) {
    // Point is not inserted. Check ivf->iloc for reason.
    return 0;
  }


  for (i = 0; i < cavetetvertlist->objects; i++) {
    cavpoints->newindex((void **) &parypt);
    *parypt = * (point *) fastlookup(cavetetvertlist, i);
  }
  // Add the new point into the point list.
  cavpoints->newindex((void **) &parypt);
  *parypt = newpt;

  for (i = 0; i < cavebdrylist->objects; i++) {
    cavfaces->newindex((void **) &parytet);
    *parytet = * (triface *) fastlookup(cavebdrylist, i);
  }

  for (i = 0; i < caveoldtetlist->objects; i++) {
    crosstets->newindex((void **) &parytet);
    *parytet = * (triface *) fastlookup(caveoldtetlist, i);
  }

  cavetetvertlist->restart();
  cavebdrylist->restart();
  caveoldtetlist->restart();

  // Insert the point using the cavity algorithm.
  delaunizecavity(cavpoints, cavfaces, cavshells, newtets, crosstets, 
                  misfaces);
  fillcavity(cavshells, NULL, NULL, NULL, NULL, NULL, NULL);
  carvecavity(crosstets, newtets, NULL);

  if ((splitsh != NULL) || (splitseg != NULL)) {
    // Insert the point into the surface mesh.
    sinsertvertex(newpt, splitsh, splitseg, ivf->sloc, ivf->sbowywat, 0);

    // Put all new subfaces into stack.
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

    if (splitseg != NULL) {
      // Queue two new subsegments in C(p) for recovery.
      for (i = 0; i < cavesegshlist->objects; i++) {
        paryseg = (face *) fastlookup(cavesegshlist, i);
        subsegstack->newindex((void **) &paryseg1);
        *paryseg1 = *paryseg;
      }
    } // if (splitseg != NULL)

    // Delete the old subfaces in sC(p).
    for (i = 0; i < caveshlist->objects; i++) {
      parysh = (face *) fastlookup(caveshlist, i);
      if (checksubfaceflag) {
        // It is possible that this subface still connects to adjacent
        //   tets which are not in C(p). If so, clear connections in the
        //   adjacent tets at this subface.
        stpivot(*parysh, neightet);
        if (neightet.tet != NULL) {
          if (neightet.tet[4] != NULL) {
            // Found an adjacent tet. It must be not in C(p).
            tsdissolve(neightet);
            fsymself(neightet);
            tsdissolve(neightet);
          }
        }
      }
      shellfacedealloc(subfaces, parysh->sh);
    }
    if (splitseg != NULL) {
      // Delete the old segment in sC(p).
      shellfacedealloc(subsegs, splitseg->sh);
    }

    // Clear working lists.
    caveshlist->restart();
    caveshbdlist->restart();
    cavesegshlist->restart();
  } // if ((splitsh != NULL) || (splitseg != NULL)) 

  // Put all interior subfaces into stack for recovery.
  // They were collected in carvecavity().
  // Note: Some collected subfaces may be deleted by sinsertvertex().
  for (i = 0; i < caveencshlist->objects; i++) {
    parysh = (face *) fastlookup(caveencshlist, i);
    if (parysh->sh[3] != NULL) {
      subfacstack->newindex((void **) &parysh1);
      *parysh1 = *parysh;
    }
  }

  // Put all interior segments into stack for recovery.
  // They were collected in carvecavity().
  // Note: Some collected segments may be deleted by sinsertvertex().
  for (i = 0; i < caveencseglist->objects; i++) {
    paryseg = (face *) fastlookup(caveencseglist, i);
    if (paryseg->sh[3] != NULL) {
      subsegstack->newindex((void **) &paryseg1);
      *paryseg1 = *paryseg;
    }
  }

  caveencshlist->restart();
  caveencseglist->restart();

  return 1;
}

//============================================================================//
//                                                                            //
// refineregion()    Refine a missing region by inserting points.             //
//                                                                            //
// 'splitsh' represents an edge of the facet to be split. It must not be a    //
// segment.                                                                   //
//                                                                            //
// Assumption: The current mesh is a CDT and is convex.                       //
//                                                                            //
//============================================================================//

void TetMeshCore::refineregion(face &splitsh, arraypool *cavpoints, 
                              arraypool *cavfaces, arraypool *cavshells,
                              arraypool *newtets, arraypool *crosstets,
                              arraypool *misfaces)
{
  triface searchtet, spintet;
  face splitseg, *paryseg;
  point steinpt, pa, pb, refpt;
  insertvertexflags ivf;
  enum interresult dir;
  long baknum = points->items;
  int t1ver;
  int i;

  // Do not split a segment.
  for (i = 0; i < 3; i++) {
    sspivot(splitsh, splitseg);
    if (splitseg.sh == NULL) break;
    senextself(splitsh);
  }

  if (b->verbose > 2) {
    printf("      Refining region at edge (%d, %d, %d).\n",
           pointmark(sorg(splitsh)), pointmark(sdest(splitsh)),
           pointmark(sapex(splitsh)));
  }

  // Add the Steiner point at the barycenter of the face.
  pa = sorg(splitsh);
  pb = sdest(splitsh);
  // Create a new point.
  makepoint(&steinpt, FREEFACETVERTEX);
  for (i = 0; i < 3; i++) {
    steinpt[i] = 0.5 * (pa[i] + pb[i]);
  }

  ivf.bowywat = 1; // Use the Bowyer-Watson algorrithm.
  ivf.cdtflag = 1; // Only create the initial cavity.
  ivf.sloc = (int) ONEDGE;
  ivf.sbowywat = 1;
  ivf.assignmeshsize = b->metric;
  ivf.smlenflag = useinsertradius; // Return the closet mesh vertex.

  point2tetorg(pa, searchtet); // Start location from it.
  ivf.iloc = (int) OUTSIDE;

  ivf.rejflag = 1; // Reject it if it encroaches upon any segment.
  if (!insertpoint_cdt(steinpt, &searchtet, &splitsh, NULL, &ivf, cavpoints,
                       cavfaces, cavshells, newtets, crosstets, misfaces)) {
    if (ivf.iloc == (int) ENCSEGMENT) {
      pointdealloc(steinpt);
      // Split an encroached segment.
      i = randomnation(encseglist->objects);
      paryseg = (face *) fastlookup(encseglist, i);
      splitseg = *paryseg;
      encseglist->restart();

      // Split the segment.
      pa = sorg(splitseg);
      pb = sdest(splitseg);
      // Create a new point.
      makepoint(&steinpt, FREESEGVERTEX);
      for (i = 0; i < 3; i++) {
        steinpt[i] = 0.5 * (pa[i] + pb[i]);
      }
      point2tetorg(pa, searchtet);
      ivf.iloc = (int) OUTSIDE;
      ivf.rejflag = 0;
      if (!insertpoint_cdt(steinpt, &searchtet, &splitsh, &splitseg, &ivf,
                           cavpoints, cavfaces, cavshells, newtets, 
                           crosstets, misfaces)) {
        terminate_tet_core(this, 2); 
      }
      if (useinsertradius) {
        //save_segmentpoint_insradius(steinpt, ivf.parentpt, ivf.smlen);
      }
      st_segref_count++;
      if (steinerleft > 0) steinerleft--;
    } else {
      terminate_tet_core(this, 2); // assert(0);
    }
  } else {
    if (useinsertradius) {
      //save_facetpoint_insradius(steinpt, ivf.parentpt, ivf.smlen);
    }
    st_facref_count++;
    if (steinerleft > 0) steinerleft--;
  }

  while (subsegstack->objects > 0l) {
    // seglist is used as a stack.
    subsegstack->objects--;
    paryseg = (face *) fastlookup(subsegstack, subsegstack->objects);
    splitseg = *paryseg;

    // Check if this segment has been recovered.
    sstpivot1(splitseg, searchtet);
    if (searchtet.tet != NULL) continue;

    // Search the segment.
    dir = scoutsegment(sorg(splitseg), sdest(splitseg), &splitseg, &searchtet, 
                       &refpt, NULL);
    if (dir == SHAREEDGE) {
      // Found this segment, insert it.
      // Let the segment remember an adjacent tet.
      sstbond1(splitseg, searchtet);
      // Bond the segment to all tets containing it.
      spintet = searchtet;
      do {
        tssbond1(spintet, splitseg);
        fnextself(spintet);
      } while (spintet.tet != searchtet.tet);
    } else { 
      if ((dir == ACROSSFACE) || (dir == ACROSSEDGE)) {
        // Split the segment.
        makepoint(&steinpt, FREESEGVERTEX);
        getsteinerptonsegment(&splitseg, refpt, steinpt);
        ivf.iloc = (int) OUTSIDE;
        ivf.rejflag = 0;
        if (!insertpoint_cdt(steinpt, &searchtet, &splitsh, &splitseg, &ivf,
                             cavpoints, cavfaces, cavshells, newtets, 
                             crosstets, misfaces)) {
          terminate_tet_core(this, 2);
        }
        if (useinsertradius) {
          //save_segmentpoint_insradius(steinpt, ivf.parentpt, ivf.smlen);
        }
        st_segref_count++;
        if (steinerleft > 0) steinerleft--;
      } else {
        terminate_tet_core(this, 2);
      }
    }
  } // while

  if (b->verbose > 2) {
    printf("      Added %ld Steiner points.\n", points->items - baknum);
  }
}

//============================================================================//
//                                                                            //
// constrainedfacets()    Recover constrained facets in a CDT.                //
//                                                                            //
// All unrecovered subfaces are queued in 'subfacestack'.                     //
//                                                                            //
//============================================================================//

void TetMeshCore::constrainedfacets()
{
  arraypool *tg_crosstets, *tg_topnewtets, *tg_botnewtets;
  arraypool *tg_topfaces, *tg_botfaces, *tg_midfaces;
  arraypool *tg_topshells, *tg_botshells, *tg_facfaces; 
  arraypool *tg_toppoints, *tg_botpoints;
  arraypool *tg_missingshs, *tg_missingshbds, *tg_missingshverts;
  triface searchtet, neightet, crossedge;
  face searchsh, *parysh, *parysh1;
  face *paryseg;
  point *parypt;
  enum interresult dir;
  int facetcount;
  int success;
  int t1ver;
  int i, j;

  // Initialize arrays.
  tg_crosstets      = new arraypool(sizeof(triface), 10);
  tg_topnewtets     = new arraypool(sizeof(triface), 10);
  tg_botnewtets     = new arraypool(sizeof(triface), 10);
  tg_topfaces       = new arraypool(sizeof(triface), 10);
  tg_botfaces       = new arraypool(sizeof(triface), 10);
  tg_midfaces       = new arraypool(sizeof(triface), 10);
  tg_toppoints      = new arraypool(sizeof(point), 8);
  tg_botpoints      = new arraypool(sizeof(point), 8);
  tg_facfaces       = new arraypool(sizeof(face), 10);
  tg_topshells      = new arraypool(sizeof(face), 10);
  tg_botshells      = new arraypool(sizeof(face), 10);
  tg_missingshs     = new arraypool(sizeof(face), 10);
  tg_missingshbds   = new arraypool(sizeof(face), 10);
  tg_missingshverts = new arraypool(sizeof(point), 8);
  // This is a global array used by refineregion().
  encseglist        = new arraypool(sizeof(face), 4);

  facetcount = 0;

  while (subfacstack->objects > 0l) {

    subfacstack->objects--;
    parysh = (face *) fastlookup(subfacstack, subfacstack->objects);
    searchsh = *parysh;

    if (searchsh.sh[3] == NULL) continue; // It is dead.
    if (isshtet(searchsh)) continue; // It is recovered.

    // Collect all unrecovered subfaces which are co-facet.
    smarktest(searchsh);
    tg_facfaces->newindex((void **) &parysh);
    *parysh = searchsh;
    for (i = 0; i < tg_facfaces->objects; i++) {
      parysh = (face *) fastlookup(tg_facfaces, i);
      for (j = 0; j < 3; j++) {
        if (!isshsubseg(*parysh)) {
          spivot(*parysh, searchsh);
          if (!smarktested(searchsh)) {
            if (!isshtet(searchsh)) {
              smarktest(searchsh);
              tg_facfaces->newindex((void **) &parysh1);
              *parysh1 = searchsh;
            }
          }
        }
        senextself(*parysh);
      } // j
    } // i
    // Have found all facet subfaces. Unmark them.
    for (i = 0; i < tg_facfaces->objects; i++) {
      parysh = (face *) fastlookup(tg_facfaces, i);
      sunmarktest(*parysh);
    }


    if (b->verbose > 1) {
      printf("    Recovering facet #%d: %ld subfaces.\n", facetcount + 1, 
             tg_facfaces->objects);      
    }
    facetcount++;

    while (tg_facfaces->objects > 0l) {

      tg_facfaces->objects--;
      parysh = (face *) fastlookup(tg_facfaces, tg_facfaces->objects);
      searchsh = *parysh;

      if (searchsh.sh[3] == NULL) continue; // It is dead.
      if (isshtet(searchsh)) continue; // It is recovered.

      searchtet.tet = NULL;
      if (scoutsubface(&searchsh, &searchtet, 1)) continue;

      // The subface is missing. Form the missing region.
      //   Re-use 'tg_crosstets' for 'adjtets'.
      formregion(&searchsh, tg_missingshs, tg_missingshbds, tg_missingshverts);

      int searchflag = scoutcrossedge(searchtet, tg_missingshbds, tg_missingshs);
      if (searchflag > 0) {
        // Save this crossing edge, will be used by fillcavity().
        crossedge = searchtet;
        // Form a cavity of crossing tets.
        success = formcavity(&searchtet, tg_missingshs, tg_crosstets,
                             tg_topfaces, tg_botfaces, tg_toppoints,
                             tg_botpoints);
        if (success) {
          if (!b->flipinsert) {
            // Tetrahedralize the top part. Re-use 'tg_midfaces'.
            delaunizecavity(tg_toppoints, tg_topfaces, tg_topshells,
                            tg_topnewtets, tg_crosstets, tg_midfaces);
            // Tetrahedralize the bottom part. Re-use 'tg_midfaces'.
            delaunizecavity(tg_botpoints, tg_botfaces, tg_botshells,
                            tg_botnewtets, tg_crosstets, tg_midfaces);
            // Fill the cavity with new tets.
            success = fillcavity(tg_topshells, tg_botshells, tg_midfaces,
                                 tg_missingshs, tg_topnewtets, tg_botnewtets,
                                 &crossedge);
            if (success) {
              // Cavity is remeshed. Delete old tets and outer new tets.
              carvecavity(tg_crosstets, tg_topnewtets, tg_botnewtets);
            } else {
              restorecavity(tg_crosstets, tg_topnewtets, tg_botnewtets,
                            tg_missingshbds);
            }
          } else {
            // Use the flip algorithm of Shewchuk to recover the subfaces.
            flipinsertfacet(tg_crosstets, tg_toppoints, tg_botpoints, 
                            tg_missingshverts);
			// Put all subfaces in R back to tg_facfaces.
		    for (i = 0; i < tg_missingshs->objects; i++) {
		      parysh = (face *) fastlookup(tg_missingshs, i);
		      tg_facfaces->newindex((void **) &parysh1);
              *parysh1 = *parysh;
		    }
			success = 1;
            // Clear working lists.
            tg_crosstets->restart();
            tg_topfaces->restart();
            tg_botfaces->restart();
            tg_toppoints->restart();
            tg_botpoints->restart();
          } // b->flipinsert

          if (success) {
            // Recover interior subfaces.
            for (i = 0; i < caveencshlist->objects; i++) {
              parysh = (face *) fastlookup(caveencshlist, i);
              if (!scoutsubface(parysh, &searchtet, 1)) {
                // Add this face at the end of the list, so it will be
                //   processed immediately.
                tg_facfaces->newindex((void **) &parysh1);
                *parysh1 = *parysh;
              }
            }
            caveencshlist->restart();
            // Recover interior segments. This should always be recovered.
            for (i = 0; i < caveencseglist->objects; i++) {
              paryseg = (face *) fastlookup(caveencseglist, i);
              dir = scoutsegment(sorg(*paryseg), sdest(*paryseg), paryseg,
                                 &searchtet, NULL, NULL);
              if (dir != SHAREEDGE) {
                terminate_tet_core(this, 2);
              }
              // Insert this segment.
              // Let the segment remember an adjacent tet.
              sstbond1(*paryseg, searchtet);
              // Bond the segment to all tets containing it.
              neightet = searchtet;
              do {
                tssbond1(neightet, *paryseg);
                fnextself(neightet);
              } while (neightet.tet != searchtet.tet);
            }
            caveencseglist->restart();
          } // success - remesh cavity
        } // success - form cavity 
        else {
          terminate_tet_core(this, 2); // Report a bug.
        } // Not success - form cavity
      } else {
        // Put all subfaces in R back to tg_facfaces.
		for (i = 0; i < tg_missingshs->objects; i++) {
		  parysh = (face *) fastlookup(tg_missingshs, i);
		  tg_facfaces->newindex((void **) &parysh1);
          *parysh1 = *parysh;
		}
        if (searchflag != -1) {
          // Some edge(s) in the missing regions were flipped.
          success = 1; 
        } else {
          restorecavity(tg_crosstets, tg_topnewtets, tg_botnewtets, 
                        tg_missingshbds); // Only remove fake segments.
          // Choose an edge to split (set in recentsh)
          recentsh = searchsh;
		  success = 0; // Do refineregion();
        }
      } // if (scoutcrossedge)

      // Unmarktest all points of the missing region.
      for (i = 0; i < tg_missingshverts->objects; i++) {
        parypt = (point *) fastlookup(tg_missingshverts, i);
        punmarktest(*parypt);
      }
      tg_missingshverts->restart();
      tg_missingshbds->restart();
      tg_missingshs->restart();

      if (!success) {
        // The missing region can not be recovered. Refine it.
        refineregion(recentsh, tg_toppoints, tg_topfaces, tg_topshells,
                     tg_topnewtets, tg_crosstets, tg_midfaces);
      }
    } // while (tg_facfaces->objects)

  } // while ((subfacstack->objects)

  // Accumulate the dynamic memory.
  totalworkmemory += (tg_crosstets->totalmemory + tg_topnewtets->totalmemory +
                      tg_botnewtets->totalmemory + tg_topfaces->totalmemory +
                      tg_botfaces->totalmemory + tg_midfaces->totalmemory +
                      tg_toppoints->totalmemory + tg_botpoints->totalmemory +
                      tg_facfaces->totalmemory + tg_topshells->totalmemory +
                      tg_botshells->totalmemory + tg_missingshs->totalmemory +
                      tg_missingshbds->totalmemory + 
                      tg_missingshverts->totalmemory + 
                      encseglist->totalmemory);

  // Delete arrays.
  delete tg_crosstets;
  delete tg_topnewtets;
  delete tg_botnewtets;
  delete tg_topfaces;
  delete tg_botfaces;
  delete tg_midfaces;
  delete tg_toppoints;
  delete tg_botpoints;
  delete tg_facfaces;
  delete tg_topshells;
  delete tg_botshells;
  delete tg_missingshs;
  delete tg_missingshbds;
  delete tg_missingshverts;
  delete encseglist;
  encseglist = NULL;
}

//============================================================================//
//                                                                            //
// constraineddelaunay()    Create a constrained Delaunay tetrahedralization. //
//                                                                            //
//============================================================================//

void TetMeshCore::constraineddelaunay(clock_t& tv)
{
  face searchsh, *parysh;
  face searchseg, *paryseg;
  int s, i;

  // Statistics.
  long bakfillregioncount;
  long bakcavitycount, bakcavityexpcount;
  long bakseg_ref_count;

  if (!b->quiet) {
    printf("Constrained Delaunay...\n");
  }

  makesegmentendpointsmap();
  makefacetverticesmap();

  if (b->verbose) {
    printf("  Delaunizing segments.\n");
  }

  checksubsegflag = 1;

  // Put all segments into the list (in random order).
  subsegs->traversalinit();
  for (i = 0; i < subsegs->items; i++) {
    s = randomnation(i + 1);
    // Move the s-th seg to the i-th.
    subsegstack->newindex((void **) &paryseg);
    *paryseg = * (face *) fastlookup(subsegstack, s);
    // Put i-th seg to be the s-th.
    searchseg.sh = shellfacetraverse(subsegs);
    //sinfect(searchseg);  // Only save it once.
    paryseg = (face *) fastlookup(subsegstack, s);
    *paryseg = searchseg;
  }

  // Recover non-Delaunay segments.
  delaunizesegments();

  if (b->verbose) {
    printf("  Inserted %ld Steiner points.\n", st_segref_count); 
  }

  tv = clock();

  if (b->verbose) {
    printf("  Constraining facets.\n");
  }

  // Subfaces will be introduced.
  checksubfaceflag = 1;

  bakfillregioncount = fillregioncount;
  bakcavitycount = cavitycount;
  bakcavityexpcount = cavityexpcount;
  bakseg_ref_count = st_segref_count;

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

  // Recover facets.
  constrainedfacets();

  if (b->verbose) {
    if (fillregioncount > bakfillregioncount) {
      printf("  Remeshed %ld regions.\n", fillregioncount-bakfillregioncount);
    }
    if (cavitycount > bakcavitycount) {
      printf("  Remeshed %ld cavities", cavitycount - bakcavitycount);
      if (cavityexpcount - bakcavityexpcount) {
        printf(" (%ld enlarged)", cavityexpcount - bakcavityexpcount);
      }
      printf(".\n");
    }
    if (st_segref_count + st_facref_count - bakseg_ref_count > 0) {
      printf("  Inserted %ld (%ld, %ld) refine points.\n", 
             st_segref_count + st_facref_count - bakseg_ref_count,
             st_segref_count - bakseg_ref_count, st_facref_count);
    }
  }
}

//                                                                            //
//                                                                            //
//== constrained_cxx =========================================================//

} // namespace sqmesh::mesh::tet::detail
