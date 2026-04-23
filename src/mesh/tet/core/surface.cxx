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

//== surface_cxx =============================================================//
//                                                                            //
//                                                                            //

//============================================================================//
//                                                                            //
// flipshpush()    Push a facet edge into flip stack.                         //
//                                                                            //
//============================================================================//

void TetMeshCore::flipshpush(face* flipedge)
{
  badface *newflipface;

  newflipface = (badface *) flippool->alloc();
  newflipface->ss = *flipedge;
  newflipface->forg = sorg(*flipedge);
  newflipface->fdest = sdest(*flipedge);
  newflipface->nextitem = flipstack;
  flipstack = newflipface;
}

//============================================================================//
//                                                                            //
// flip22()    Perform a 2-to-2 flip in surface mesh.                         //
//                                                                            //
// 'flipfaces' is an array of two subfaces. On input, they are [a,b,c] and    //
// [b,a,d]. On output, they are [c,d,b] and [d,c,a]. As a result, edge [a,b]  //
// is replaced by edge [c,d].                                                 //
//                                                                            //
//============================================================================//

void TetMeshCore::flip22(face* flipfaces, int flipflag, int chkencflag)
{
  face bdedges[4], outfaces[4], infaces[4];
  face bdsegs[4];
  face checkface;
  point pa, pb, pc, pd;
  int i;

  pa = sorg(flipfaces[0]);
  pb = sdest(flipfaces[0]);
  pc = sapex(flipfaces[0]);
  pd = sapex(flipfaces[1]);

  if (sorg(flipfaces[1]) != pb) {
    sesymself(flipfaces[1]);
  }

  flip22count++;

  // Collect the four boundary edges.
  senext(flipfaces[0], bdedges[0]);
  senext2(flipfaces[0], bdedges[1]);
  senext(flipfaces[1], bdedges[2]);
  senext2(flipfaces[1], bdedges[3]);

  // Collect outer boundary faces.
  for (i = 0; i < 4; i++) {
    spivot(bdedges[i], outfaces[i]);
    infaces[i] = outfaces[i];
    sspivot(bdedges[i], bdsegs[i]);
    if (outfaces[i].sh != NULL) {
      if (isshsubseg(bdedges[i])) {
        spivot(infaces[i], checkface);
        while (checkface.sh != bdedges[i].sh) {
          infaces[i] = checkface;
          spivot(infaces[i], checkface);
        }
      }
    }
  }

  // The flags set in these two subfaces do not change.
  // Shellmark does not change.
  // area constraint does not change.

  // Transform [a,b,c] -> [c,d,b].
  setshvertices(flipfaces[0], pc, pd, pb);
  // Transform [b,a,d] -> [d,c,a].
  setshvertices(flipfaces[1], pd, pc, pa);

  // Update the point-to-subface map.
  if (pointtype(pa) == FREEFACETVERTEX) {
    setpoint2sh(pa, sencode(flipfaces[1]));
  }
  if (pointtype(pb) == FREEFACETVERTEX) {
    setpoint2sh(pb, sencode(flipfaces[0]));
  }
  if (pointtype(pc) == FREEFACETVERTEX) {
    setpoint2sh(pc, sencode(flipfaces[0]));
  }
  if (pointtype(pd) == FREEFACETVERTEX) {
    setpoint2sh(pd, sencode(flipfaces[0]));
  }

  // Reconnect boundary edges to outer boundary faces.
  for (i = 0; i < 4; i++) {
    if (outfaces[(3 + i) % 4].sh != NULL) {
      // Make sure that the subface has the ori as the segment.
      if (bdsegs[(3 + i) % 4].sh != NULL) {
        bdsegs[(3 + i) % 4].shver = 0;
        if (sorg(bdedges[i]) != sorg(bdsegs[(3 + i) % 4])) {
          sesymself(bdedges[i]);
        }
      }
      sbond1(bdedges[i], outfaces[(3 + i) % 4]);
      sbond1(infaces[(3 + i) % 4], bdedges[i]);
    } else {
      sdissolve(bdedges[i]);
    }
    if (bdsegs[(3 + i) % 4].sh != NULL) {
      ssbond(bdedges[i], bdsegs[(3 + i) % 4]);
      if (chkencflag & 1) {
        // Queue this segment for encroaching check.
        enqueuesubface(badsubsegs, &(bdsegs[(3 + i) % 4]));
      }
    } else {
      ssdissolve(bdedges[i]);
    }
  }

  if (chkencflag & 2) {
    // Queue the flipped subfaces for quality/encroaching checks.
    for (i = 0; i < 2; i++) {
      enqueuesubface(badsubfacs, &(flipfaces[i]));
    }
  }

  recentsh = flipfaces[0];

  if (flipflag) {
    // Put the boundary edges into flip stack.
    for (i = 0; i < 4; i++) {
      flipshpush(&(bdedges[i]));
    }
  }
}

//============================================================================//
//                                                                            //
// flip31()    Remove a vertex by transforming 3-to-1 subfaces.               //
//                                                                            //
// 'flipfaces' is an array of subfaces. Its length is at least 4.  On input,  //
// the first three faces are: [p,a,b], [p,b,c], and [p,c,a]. This routine     //
// replaces them by one face [a,b,c], it is returned in flipfaces[3].         //
//                                                                            //
// NOTE: The three old subfaces are not deleted within this routine.  They    //
// still hold pointers to their adjacent subfaces. These informations are     //
// needed by the routine 'sremovevertex()' for recovering a segment.          //
// The caller of this routine must delete the old subfaces after their uses.  //
//                                                                            //
//============================================================================//

void TetMeshCore::flip31(face* flipfaces, int flipflag)
{
  face bdedges[3], outfaces[3], infaces[3];
  face bdsegs[3];
  face checkface;
  point pa, pb, pc;
  int i;

  pa = sdest(flipfaces[0]);
  pb = sdest(flipfaces[1]);
  pc = sdest(flipfaces[2]);

  flip31count++;

  // Collect all infos at the three boundary edges.
  for (i = 0; i < 3; i++) {
    senext(flipfaces[i], bdedges[i]);
    spivot(bdedges[i], outfaces[i]);
    infaces[i] = outfaces[i];
    sspivot(bdedges[i], bdsegs[i]);
    if (outfaces[i].sh != NULL) {
      if (isshsubseg(bdedges[i])) {
        spivot(infaces[i], checkface);
        while (checkface.sh != bdedges[i].sh) {
          infaces[i] = checkface;
          spivot(infaces[i], checkface);
        }
      }
    }
  } // i

  // Create a new subface.
  makeshellface(subfaces, &(flipfaces[3]));
  setshvertices(flipfaces[3], pa, pb,pc);
  setshellmark(flipfaces[3], shellmark(flipfaces[0]));
  if (checkconstraints) {
    //area = areabound(flipfaces[0]);
    setareabound(flipfaces[3], areabound(flipfaces[0]));
  }
  if (useinsertradius) {
    setfacetindex(flipfaces[3], getfacetindex(flipfaces[0]));
  }

  // Update the point-to-subface map.
  if (pointtype(pa) == FREEFACETVERTEX) {
    setpoint2sh(pa, sencode(flipfaces[3]));
  }
  if (pointtype(pb) == FREEFACETVERTEX) {
    setpoint2sh(pb, sencode(flipfaces[3]));
  }
  if (pointtype(pc) == FREEFACETVERTEX) {
    setpoint2sh(pc, sencode(flipfaces[3]));
  }

  // Update the three new boundary edges.
  bdedges[0] = flipfaces[3];         // [a,b]
  senext(flipfaces[3], bdedges[1]);  // [b,c]
  senext2(flipfaces[3], bdedges[2]); // [c,a]

  // Reconnect boundary edges to outer boundary faces.
  for (i = 0; i < 3; i++) {
    if (outfaces[i].sh != NULL) {
      // Make sure that the subface has the ori as the segment.
      if (bdsegs[i].sh != NULL) {
        bdsegs[i].shver = 0;
        if (sorg(bdedges[i]) != sorg(bdsegs[i])) {
          sesymself(bdedges[i]);
        }
      }
      sbond1(bdedges[i], outfaces[i]);
      sbond1(infaces[i], bdedges[i]);
    }
    if (bdsegs[i].sh != NULL) {
      ssbond(bdedges[i], bdsegs[i]);
    }
  }

  recentsh = flipfaces[3];

  if (flipflag) {
    // Put the boundary edges into flip stack.
    for (i = 0; i < 3; i++) {
      flipshpush(&(bdedges[i]));
    }
  }
}

//============================================================================//
//                                                                            //
// lawsonflip()    Flip non-locally Delaunay edges.                           //
//                                                                            //
//============================================================================//

long TetMeshCore::lawsonflip()
{
  badface *popface;
  face flipfaces[2];
  point pa, pb, pc, pd;
  double sign;
  long flipcount = 0;

  if (b->verbose > 2) {
    printf("      Lawson flip %ld edges.\n", flippool->items);
  }

  while (flipstack != (badface *) NULL) {

    // Pop an edge from the stack.
    popface = flipstack;
    flipfaces[0] = popface->ss;
    pa = popface->forg;
    pb = popface->fdest;
    flipstack = popface->nextitem; // The next top item in stack.
    flippool->dealloc((void *) popface);

    // Skip it if it is dead.
    if (flipfaces[0].sh[3] == NULL) continue;
    // Skip it if it is not the same edge as we saved.
    if ((sorg(flipfaces[0]) != pa) || (sdest(flipfaces[0]) != pb)) continue;
    // Skip it if it is a subsegment.
    if (isshsubseg(flipfaces[0])) continue;

    // Get the adjacent face.
    spivot(flipfaces[0], flipfaces[1]);
    if (flipfaces[1].sh == NULL) continue; // Skip a hull edge.
    pc = sapex(flipfaces[0]);
    pd = sapex(flipfaces[1]);

    sign = incircle3d(pa, pb, pc, pd);

    if (sign < 0) {
      // It is non-locally Delaunay. Flip it.
      flip22(flipfaces, 1, 0);
      flipcount++;
    }
  }

  if (b->verbose > 2) {
    printf("      Performed %ld flips.\n", flipcount);
  }

  return flipcount;
}

//============================================================================//
//                                                                            //
// sinsertvertex()    Insert a vertex into a triangulation of a facet.        //
//                                                                            //
// This function uses three global arrays: 'caveshlist', 'caveshbdlist', and  //
// 'caveshseglist'. On return, 'caveshlist' contains old subfaces in C(p),    //
// 'caveshbdlist' contains new subfaces in C(p). If the new point lies on a   //
// segment, 'cavesegshlist' returns the two new subsegments.                  //
//                                                                            //
// 'iloc' suggests the location of the point. If it is OUTSIDE, this routine  //
// will first locate the point. It starts searching from 'searchsh' or 'rec-  //
// entsh' if 'searchsh' is NULL.                                              //
//                                                                            //
// If 'bowywat' is set (1), the Bowyer-Watson algorithm is used to insert     //
// the vertex. Otherwise, only insert the vertex in the initial cavity.       //
//                                                                            //
// If 'iloc' is 'INSTAR', this means the cavity of this vertex was already    //
// provided in the list 'caveshlist'.                                         //
//                                                                            //
// If 'splitseg' is not NULL, the new vertex lies on the segment and it will  //
// be split. 'iloc' must be either 'ONEDGE' or 'INSTAR'.                      //
//                                                                            //
// 'rflag' (rounding) is a parameter passed to slocate() function.  If it is  //
// set, after the location of the point is found, either ONEDGE or ONFACE,    //
// round the result using an epsilon.                                         //
//                                                                            //
// NOTE: the old subfaces in C(p) are not deleted. They're needed in case we  //
// want to remove the new point immediately.                                  //
//                                                                            //
//============================================================================//

int TetMeshCore::sinsertvertex(point insertpt, face *searchsh, face *splitseg,
                              int iloc, int bowywat, int rflag)
{
  face cavesh, neighsh, *parysh;
  face newsh, casout, casin;
  face checkseg;
  point pa, pb;
  enum locateresult loc = OUTSIDE;
  double sign, ori;
  int i, j;

  if (b->verbose > 2) {
    printf("      Insert facet point %d.\n", pointmark(insertpt));
  }

  if (bowywat == 3) {
    loc = INSTAR;
  }

  if ((splitseg != NULL) && (splitseg->sh != NULL)) {
    // A segment is going to be split, no point location.
    spivot(*splitseg, *searchsh);
    if (loc != INSTAR) loc = ONEDGE;
  } else {
    if (loc != INSTAR) loc = (enum locateresult) iloc;
    if (loc == OUTSIDE) {
      // Do point location in surface mesh.
      if (searchsh->sh == NULL) {
        *searchsh = recentsh;
      }
      // Search the vertex. An above point must be provided ('aflag' = 1).
      loc = slocate(insertpt, searchsh, 1, 1, rflag);
    }
  }


  // Form the initial sC(p).
  if (loc == ONFACE) {
    // Add the face into list (in B-W cavity).
    smarktest(*searchsh);
    caveshlist->newindex((void **) &parysh);
    *parysh = *searchsh;
  } else if (loc == ONEDGE) {
    if ((splitseg != NULL) && (splitseg->sh != NULL)) {
      splitseg->shver = 0;
      pa = sorg(*splitseg);
    } else {
      pa = sorg(*searchsh);
    }
    if (searchsh->sh != NULL) {
      // Collect all subfaces share at this edge.
      neighsh = *searchsh;
      while (1) {
        // Adjust the origin of its edge to be 'pa'.
        if (sorg(neighsh) != pa) sesymself(neighsh);
        // Add this face into list (in B-W cavity).
        smarktest(neighsh);
        caveshlist->newindex((void **) &parysh);
        *parysh = neighsh;
        // Add this face into face-at-splitedge list.
        cavesegshlist->newindex((void **) &parysh);
        *parysh = neighsh;
        // Go to the next face at the edge.
        spivotself(neighsh);
        // Stop if all faces at the edge have been visited.
        if (neighsh.sh == searchsh->sh) break;
        if (neighsh.sh == NULL) break;
      }
    } // If (not a non-dangling segment).
  } else if (loc == ONVERTEX) {
    return (int) loc;
  } else if (loc == OUTSIDE) {
    // Comment: This should only happen during the surface meshing step.
    // Enlarge the convex hull of the triangulation by including p.
    // An above point of the facet is set in 'dummypoint' to replace
    // orient2d tests by orient3d tests.
    // Imagine that the current edge a->b (in 'searchsh') is horizontal in a
    //   plane, and a->b is directed from left to right, p lies above a->b.  
    //   Find the right-most edge of the triangulation which is visible by p.
    neighsh = *searchsh;
    while (1) {
      senext2self(neighsh);
      spivot(neighsh, casout);
      if (casout.sh == NULL) {
        // A convex hull edge. Is it visible by p.
        ori = orient3d(sorg(neighsh), sdest(neighsh), dummypoint, insertpt);
        if (ori < 0) {
          *searchsh = neighsh; // Visible, update 'searchsh'.
        } else {
          break; // 'searchsh' is the right-most visible edge.
        }
      } else {
        if (sorg(casout) != sdest(neighsh)) sesymself(casout);
        neighsh = casout;
      }
    }
    // Create new triangles for all visible edges of p (from right to left).
    casin.sh = NULL;  // No adjacent face at right.
    pa = sorg(*searchsh);
    pb = sdest(*searchsh);
    while (1) {
      // Create a new subface on top of the (visible) edge.
      makeshellface(subfaces, &newsh); 
      setshvertices(newsh, pb, pa, insertpt);
      setshellmark(newsh, shellmark(*searchsh));
      if (checkconstraints) {
        //area = areabound(*searchsh);
        setareabound(newsh, areabound(*searchsh));
      }
      if (useinsertradius) {
        setfacetindex(newsh, getfacetindex(*searchsh));
      }
      // Connect the new subface to the bottom subfaces.
      sbond1(newsh, *searchsh);
      sbond1(*searchsh, newsh);
      // Connect the new subface to its right-adjacent subface.
      if (casin.sh != NULL) {
        senext(newsh, casout);
        sbond1(casout, casin);
        sbond1(casin, casout);
      }
      // The left-adjacent subface has not been created yet.
      senext2(newsh, casin);
      // Add the new face into list (inside the B-W cavity).
      smarktest(newsh);
      caveshlist->newindex((void **) &parysh);
      *parysh = newsh;
      // Move to the convex hull edge at the left of 'searchsh'.
      neighsh = *searchsh;
      while (1) {
        senextself(neighsh);
        spivot(neighsh, casout);
        if (casout.sh == NULL) {
          *searchsh = neighsh;
          break;
        }
        if (sorg(casout) != sdest(neighsh)) sesymself(casout);
        neighsh = casout;
      }
      // A convex hull edge. Is it visible by p.
      pa = sorg(*searchsh);
      pb = sdest(*searchsh);
      ori = orient3d(pa, pb, dummypoint, insertpt);
      // Finish the process if p is not visible by the hull edge.
      if (ori >= 0) break;
    }
  } else if (loc == INSTAR) {
    // Under this case, the sub-cavity sC(p) has already been formed in
    //   insertvertex().
  }

  // Form the Bowyer-Watson cavity sC(p).
  for (i = 0; i < caveshlist->objects; i++) {
    cavesh = * (face *) fastlookup(caveshlist, i);
    for (j = 0; j < 3; j++) {
      if (!isshsubseg(cavesh)) {
        spivot(cavesh, neighsh);
        if (neighsh.sh != NULL) {
          // The adjacent face exists.
          if (!smarktested(neighsh)) {
            if (bowywat) {
              if (loc == INSTAR) { // if (bowywat > 2) {
                // It must be a boundary edge.
                sign = 1;
              } else {
                // Check if this subface is connected to adjacent tet(s).
                if (!isshtet(neighsh)) {
                  // Check if the subface is non-Delaunay wrt. the new pt.
                  sign = incircle3d(sorg(neighsh), sdest(neighsh), 
                                    sapex(neighsh), insertpt);
                } else {
                  // It is connected to an adjacent tet. A boundary edge.
                  sign = 1;
                }
              }
              if (sign < 0) {
                // Add the adjacent face in list (in B-W cavity).
                smarktest(neighsh);
                caveshlist->newindex((void **) &parysh);
                *parysh = neighsh;
              }
            } else {
              sign = 1; // A boundary edge.
            }
          } else {
            sign = -1; // Not a boundary edge.
          }
        } else {
          // No adjacent face. It is a hull edge.
          if (loc == OUTSIDE) {
            // It is a boundary edge if it does not contain p.
            if ((sorg(cavesh) == insertpt) || (sdest(cavesh) == insertpt)) {
              sign = -1; // Not a boundary edge.
            } else {
              sign = 1; // A boundary edge.
            }
          } else {
            sign = 1; // A boundary edge.
          }
        }
      } else {
        // Do not across a segment. It is a boundary edge.
        sign = 1;
      }
      if (sign >= 0) {
        // Add a boundary edge.
        caveshbdlist->newindex((void **) &parysh);
        *parysh = cavesh;
      }
      senextself(cavesh);
    } // j
  } // i


  // Creating new subfaces.
  for (i = 0; i < caveshbdlist->objects; i++) {
    parysh = (face *) fastlookup(caveshbdlist, i);
    sspivot(*parysh, checkseg);
    if ((parysh->shver & 01) != 0) sesymself(*parysh);
    pa = sorg(*parysh);
    pb = sdest(*parysh);
    // Create a new subface.
    makeshellface(subfaces, &newsh); 
    setshvertices(newsh, pa, pb, insertpt);
    setshellmark(newsh, shellmark(*parysh));
    if (checkconstraints) {
      //area = areabound(*parysh);
      setareabound(newsh, areabound(*parysh));
    }
    if (useinsertradius) {
      setfacetindex(newsh, getfacetindex(*parysh));
    }
    // Update the point-to-subface map.
    if (pointtype(pa) == FREEFACETVERTEX) {
      setpoint2sh(pa, sencode(newsh));
    }
    if (pointtype(pb) == FREEFACETVERTEX) {
      setpoint2sh(pb, sencode(newsh));
    }
    // Connect newsh to outer subfaces.
    spivot(*parysh, casout);
    if (casout.sh != NULL) {
      casin = casout;
      if (checkseg.sh != NULL) {
        // Make sure that newsh has the right ori at this segment.
        checkseg.shver = 0;
        if (sorg(newsh) != sorg(checkseg)) {
          sesymself(newsh);
          sesymself(*parysh); // This side should also be inverse.
        }
        spivot(casin, neighsh);
        while (neighsh.sh != parysh->sh) {
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
    // Connect oldsh <== newsh (for connecting adjacent new subfaces).
    //   *parysh and newsh point to the same edge and the same ori.
    sbond1(*parysh, newsh);
  }

  if (newsh.sh != NULL) {
    // Set a handle for searching.
    recentsh = newsh;
  }

  // Update the point-to-subface map.
  if (pointtype(insertpt) == FREEFACETVERTEX) {
    setpoint2sh(insertpt, sencode(newsh));
  }

  // Connect adjacent new subfaces together.
  for (i = 0; i < caveshbdlist->objects; i++) {
    // Get an old subface at edge [a, b].
    parysh = (face *) fastlookup(caveshbdlist, i);
    spivot(*parysh, newsh); // The new subface [a, b, p].
    senextself(newsh); // At edge [b, p].
    spivot(newsh, neighsh);
    if (neighsh.sh == NULL) {
      // Find the adjacent new subface at edge [b, p].
      pb = sdest(*parysh);
      neighsh = *parysh;
      while (1) {
        senextself(neighsh);
        spivotself(neighsh);
        if (neighsh.sh == NULL) break;
        if (!smarktested(neighsh)) break;
        if (sdest(neighsh) != pb) sesymself(neighsh);
      }
      if (neighsh.sh != NULL) {
        // Now 'neighsh' is a new subface at edge [b, #].
        if (sorg(neighsh) != pb) sesymself(neighsh);
        senext2self(neighsh); // Go to the open edge [p, b].
        sbond(newsh, neighsh);
      }
    }
    spivot(*parysh, newsh); // The new subface [a, b, p].
    senext2self(newsh); // At edge [p, a].
    spivot(newsh, neighsh);
    if (neighsh.sh == NULL) {
      // Find the adjacent new subface at edge [p, a].
      pa = sorg(*parysh);
      neighsh = *parysh;
      while (1) {
        senext2self(neighsh);
        spivotself(neighsh);
        if (neighsh.sh == NULL) break;
        if (!smarktested(neighsh)) break;
        if (sorg(neighsh) != pa) sesymself(neighsh);
      }
      if (neighsh.sh != NULL) {
        // Now 'neighsh' is a new subface at edge [#, a].
        if (sdest(neighsh) != pa) sesymself(neighsh);
        senextself(neighsh); // Go to the open edge [a, p].
        sbond(newsh, neighsh);
      }
    }
  }

  if ((loc == ONEDGE) || ((splitseg != NULL) && (splitseg->sh != NULL))
      || (cavesegshlist->objects > 0l)) {
    // An edge is being split. We distinguish two cases:
    //   (1) the edge is not on the boundary of the cavity;
    //   (2) the edge is on the boundary of the cavity.
    // In case (2), the edge is either a segment or a hull edge. There are
    //   degenerated new faces in the cavity. They must be removed.
    face aseg, bseg, aoutseg, boutseg;

    for (i = 0; i < cavesegshlist->objects; i++) {
      // Get the saved old subface.
      parysh = (face *) fastlookup(cavesegshlist, i);
      // Get a possible new degenerated subface.
      spivot(*parysh, cavesh);
      if (sapex(cavesh) == insertpt) {
        // Found a degenerated new subface, i.e., case (2).
        if (cavesegshlist->objects > 1) {
          // There are more than one subface share at this edge.
          j = (i + 1) % (int) cavesegshlist->objects;
          parysh = (face *) fastlookup(cavesegshlist, j);
          spivot(*parysh, neighsh);
          // Adjust cavesh and neighsh both at edge a->b, and has p as apex.
          if (sorg(neighsh) != sorg(cavesh)) {
            sesymself(neighsh);
          }
          // Connect adjacent faces at two other edges of cavesh and neighsh.
          //   As a result, the two degenerated new faces are squeezed from the
          //   new triangulation of the cavity. Note that the squeezed faces
          //   still hold the adjacent informations which will be used in 
          //   re-connecting subsegments (if they exist). 
          for (j = 0; j < 2; j++) { 
            senextself(cavesh);
            senextself(neighsh);
            spivot(cavesh, newsh);
            spivot(neighsh, casout);
            sbond1(newsh, casout); // newsh <- casout.
          }
        } else {
          // There is only one subface containing this edge [a,b]. Squeeze the
          //   degenerated new face [a,b,c] by disconnecting it from its two 
          //   adjacent subfaces at edges [b,c] and [c,a]. Note that the face
          //   [a,b,c] still hold the connection to them.
          for (j = 0; j < 2; j++) {
            senextself(cavesh);
            spivot(cavesh, newsh);
            sdissolve(newsh);
          }
        }
        //recentsh = newsh;
        // Update the point-to-subface map.
        if (pointtype(insertpt) == FREEFACETVERTEX) {
          setpoint2sh(insertpt, sencode(newsh));
        }
      }
    }

    if ((splitseg != NULL) && (splitseg->sh != NULL)) {
      if (loc != INSTAR) { // if (bowywat < 3) {
        smarktest(*splitseg); // Mark it as being processed.
      }
      
      aseg = *splitseg;
      pa = sorg(*splitseg);
      pb = sdest(*splitseg);

      // Insert the new point p.
      makeshellface(subsegs, &aseg);
      makeshellface(subsegs, &bseg);

      setshvertices(aseg, pa, insertpt, NULL);
      setshvertices(bseg, insertpt, pb, NULL);
      setshellmark(aseg, shellmark(*splitseg));
      setshellmark(bseg, shellmark(*splitseg));
      if (checkconstraints) {
        setareabound(aseg, areabound(*splitseg));
        setareabound(bseg, areabound(*splitseg));
      }
      if (useinsertradius) {
        setfacetindex(aseg, getfacetindex(*splitseg));
        setfacetindex(bseg, getfacetindex(*splitseg));
      }

      // Connect [#, a]<->[a, p].
      senext2(*splitseg, boutseg); // Temporarily use boutseg.
      spivotself(boutseg);
      if (boutseg.sh != NULL) {
        senext2(aseg, aoutseg);
        sbond(boutseg, aoutseg);
      }
      // Connect [p, b]<->[b, #].
      senext(*splitseg, aoutseg);
      spivotself(aoutseg);
      if (aoutseg.sh != NULL) {
        senext(bseg, boutseg);
        sbond(boutseg, aoutseg);
      }
      // Connect [a, p] <-> [p, b].
      senext(aseg, aoutseg);
      senext2(bseg, boutseg);
      sbond(aoutseg, boutseg);

      // Connect subsegs [a, p] and [p, b] to adjacent new subfaces.
      // Although the degenerated new faces have been squeezed. They still
      //   hold the connections to the actual new faces. 
      for (i = 0; i < cavesegshlist->objects; i++) {        
        parysh = (face *) fastlookup(cavesegshlist, i);
        spivot(*parysh, neighsh);
        // neighsh is a degenerated new face.
        if (sorg(neighsh) != pa) {
          sesymself(neighsh);
        }
        senext2(neighsh, newsh);
        spivotself(newsh); // The edge [p, a] in newsh
        ssbond(newsh, aseg);
        senext(neighsh, newsh);
        spivotself(newsh); // The edge [b, p] in newsh
        ssbond(newsh, bseg);
      }


      // Let the point remember the segment it lies on.
      if (pointtype(insertpt) == FREESEGVERTEX) {
        setpoint2sh(insertpt, sencode(aseg));
      }
      // Update the point-to-seg map.
      if (pointtype(pa) == FREESEGVERTEX) {
        setpoint2sh(pa, sencode(aseg));
      }
      if (pointtype(pb) == FREESEGVERTEX) {
        setpoint2sh(pb, sencode(bseg));
      }
    } // if ((splitseg != NULL) && (splitseg->sh != NULL)) 

    // Delete all degenerated new faces.
    for (i = 0; i < cavesegshlist->objects; i++) {
      parysh = (face *) fastlookup(cavesegshlist, i);
      spivotself(*parysh);
      if (sapex(*parysh) == insertpt) {
        shellfacedealloc(subfaces, parysh->sh);
      }
    }
    cavesegshlist->restart();

    if ((splitseg != NULL) && (splitseg->sh != NULL)) {
      // Return the two new subsegments (for further process).
      //   Re-use 'cavesegshlist'.
      cavesegshlist->newindex((void **) &parysh);
      *parysh = aseg;
      cavesegshlist->newindex((void **) &parysh);
      *parysh = bseg;
    }
  } // if (loc == ONEDGE)


  return (int) loc;
}

//============================================================================//
//                                                                            //
// sremovevertex()    Remove a vertex from the surface mesh.                  //
//                                                                            //
// 'delpt' (p) is the vertex to be removed. If 'parentseg' is not NULL, p is  //
// a segment vertex, and the origin of 'parentseg' is p. Otherwise, p is a    //
// facet vertex, and the origin of 'parentsh' is p.                           //
//                                                                            //
// Within each facet, we first use a sequence of 2-to-2 flips to flip any     //
// edge at p, finally use a 3-to-1 flip to remove p.                          //
//                                                                            //
// All new created subfaces are returned in the global array 'caveshbdlist'.  //
// The new segment (when p is on segment) is returned in 'parentseg'.         //
//                                                                            //
// If 'lawson' > 0, the Lawson flip algorithm is used to recover Delaunay-    //
// ness after p is removed.                                                   //
//                                                                            //
//============================================================================//

int TetMeshCore::sremovevertex(point delpt, face* parentsh, face* parentseg,
                              int lawson)
{
  face flipfaces[4], spinsh, *parysh;
  point pa, pb, pc, pd;
  double ori1, ori2;
  int it, i, j;

  if (parentseg != NULL) {
    // 'delpt' (p) should be a Steiner point inserted in a segment [a,b],
    //   where 'parentseg' should be [p,b]. Find the segment [a,p].
    face startsh, neighsh, nextsh;
    face abseg, prevseg, checkseg;
    face adjseg1, adjseg2;
    face fakesh;
    senext2(*parentseg, prevseg);
    spivotself(prevseg);
    prevseg.shver = 0;
    // Restore the original segment [a,b].
    pa = sorg(prevseg);
    pb = sdest(*parentseg);
    if (b->verbose > 2) {
      printf("      Remove vertex %d from segment [%d, %d].\n", 
             pointmark(delpt), pointmark(pa), pointmark(pb));
    }
    makeshellface(subsegs, &abseg);
    setshvertices(abseg, pa, pb, NULL);
    setshellmark(abseg, shellmark(*parentseg));
    if (checkconstraints) {
      setareabound(abseg, areabound(*parentseg));
    }
    if (useinsertradius) {
      setfacetindex(abseg, getfacetindex(*parentseg));
    }
    // Connect [#, a]<->[a, b].
    senext2(prevseg, adjseg1);
    spivotself(adjseg1);
    if (adjseg1.sh != NULL) {
      adjseg1.shver = 0;
      senextself(adjseg1);
      senext2(abseg, adjseg2);
      sbond(adjseg1, adjseg2);
    }
    // Connect [a, b]<->[b, #].
    senext(*parentseg, adjseg1);
    spivotself(adjseg1);
    if (adjseg1.sh != NULL) {
      adjseg1.shver = 0;
      senext2self(adjseg1);
      senext(abseg, adjseg2);
      sbond(adjseg1, adjseg2);
    }
    // Update the point-to-segment map.
    setpoint2sh(pa, sencode(abseg));
    setpoint2sh(pb, sencode(abseg));

    // Get the faces in face ring at segment [p, b].
    //   Re-use array 'caveshlist'.
    spivot(*parentseg, *parentsh);
    if (parentsh->sh != NULL) {
      spinsh = *parentsh;
      while (1) {
        // Save this face in list.
        caveshlist->newindex((void **) &parysh);
        *parysh = spinsh;
        // Go to the next face in the ring.
        spivotself(spinsh);
        if (spinsh.sh == NULL) {
          break; // It is possible there is only one facet.
        }
        if (spinsh.sh == parentsh->sh) break;
      }
    }

    // Create the face ring of the new segment [a,b]. Each face in the ring
    //   is [a,b,p] (degenerated!). It will be removed (automatically).
    for (i = 0; i < caveshlist->objects; i++) {
      parysh = (face *) fastlookup(caveshlist, i);
      startsh = *parysh;
      if (sorg(startsh) != delpt) {
        sesymself(startsh);
      }      
      // startsh is [p, b, #1], find the subface [a, p, #2].
      neighsh = startsh;
      while (1) {
        senext2self(neighsh);
        sspivot(neighsh, checkseg);
        if (checkseg.sh != NULL) {
          // It must be the segment [a, p].
          break;
        }
        spivotself(neighsh);
        if (sorg(neighsh) != delpt) sesymself(neighsh);
      }
      // Now neighsh is [a, p, #2].
      if (neighsh.sh != startsh.sh) {
        // Detach the two subsegments [a,p] and [p,b] from subfaces.
        ssdissolve(startsh);
        ssdissolve(neighsh);
        // Create a degenerated subface [a,b,p]. It is used to: (1) hold the
        //   new segment [a,b]; (2) connect to the two adjacent subfaces
        //   [p,b,#] and [a,p,#].
        makeshellface(subfaces, &fakesh);
        setshvertices(fakesh, pa, pb, delpt);
        setshellmark(fakesh, shellmark(startsh));
        // Connect fakesh to the segment [a,b].
        ssbond(fakesh, abseg);
        // Connect fakesh to adjacent subfaces: [p,b,#1] and [a,p,#2].
        senext(fakesh, nextsh);
        sbond(nextsh, startsh);
        senext2(fakesh, nextsh);
        sbond(nextsh, neighsh);
        smarktest(fakesh); // Mark it as faked.
      } else {
        // Special case. There exists already a degenerated face [a,b,p]!
        //   There is no need to create a faked subface here.
        senext2self(neighsh); // [a,b,p]
        // Since we will re-connect the face ring using the faked subfaces.
        //   We put the adjacent face of [a,b,p] to the list.
        spivot(neighsh, startsh); // The original adjacent subface.
        if (sorg(startsh) != pa) sesymself(startsh);
        sdissolve(startsh);
        // Connect fakesh to the segment [a,b].
        ssbond(startsh, abseg);
        fakesh = startsh; // Do not mark it!
        // Delete the degenerated subface.
        shellfacedealloc(subfaces, neighsh.sh);
      }
      // Save the fakesh in list (for re-creating the face ring).
      cavesegshlist->newindex((void **) &parysh);
      *parysh = fakesh;
    } // i
    caveshlist->restart();

    // Re-create the face ring.
    if (cavesegshlist->objects > 1) {
      for (i = 0; i < cavesegshlist->objects; i++) {
        parysh = (face *) fastlookup(cavesegshlist, i);
        fakesh = *parysh;
        // Get the next face in the ring.
        j = (i + 1) % cavesegshlist->objects;
        parysh = (face *) fastlookup(cavesegshlist, j);
        nextsh = *parysh;
        sbond1(fakesh, nextsh);
      }
    }

    // Delete the two subsegments containing p.
    shellfacedealloc(subsegs, parentseg->sh);
    shellfacedealloc(subsegs, prevseg.sh);
    // Return the new segment.
    *parentseg = abseg;
  } else {
    // p is inside the surface.
    if (b->verbose > 2) {
      printf("      Remove vertex %d from surface.\n", pointmark(delpt));
    }
    // Let 'delpt' be its apex.
    senextself(*parentsh);
    // For unifying the code, we add parentsh to list.
    cavesegshlist->newindex((void **) &parysh);
    *parysh = *parentsh;
  }

  // Remove the point (p).

  for (it = 0; it < cavesegshlist->objects; it++) {
    parentsh = (face *) fastlookup(cavesegshlist, it); // [a,b,p]
    senextself(*parentsh); // [b,p,a].
    spivotself(*parentsh);
    if (sorg(*parentsh) != delpt) sesymself(*parentsh);
    // now parentsh is [p,b,#].
    if (sorg(*parentsh) != delpt) {
      // The vertex has already been removed in above special case.
      continue;
    }

    while (1) {      
      // Initialize the flip edge list. Re-use 'caveshlist'.
      spinsh = *parentsh; // [p, b, #]
      while (1) {
        caveshlist->newindex((void **) &parysh);
        *parysh = spinsh;
        senext2self(spinsh);
        spivotself(spinsh);
        if (spinsh.sh == parentsh->sh) break;
        if (sorg(spinsh) != delpt) sesymself(spinsh);
      } // while (1)

      if (caveshlist->objects == 3) {
        // Delete the point by a 3-to-1 flip.
        for (i = 0; i < 3; i++) {
          parysh = (face *) fastlookup(caveshlist, i);
          flipfaces[i] = *parysh;
        }
        flip31(flipfaces, lawson);
        for (i = 0; i < 3; i++) { 
          shellfacedealloc(subfaces, flipfaces[i].sh);
        }
        caveshlist->restart();
        // Save the new subface.
        caveshbdlist->newindex((void **) &parysh);
        *parysh = flipfaces[3];
        // The vertex is removed.
        break;
      }

      // Search an edge to flip.
      for (i = 0; i < caveshlist->objects; i++) {
        parysh = (face *) fastlookup(caveshlist, i);
        flipfaces[0] = *parysh;
        spivot(flipfaces[0], flipfaces[1]);
        if (sorg(flipfaces[0]) != sdest(flipfaces[1])) 
          sesymself(flipfaces[1]);
        // Skip this edge if it belongs to a faked subface.
        if (!smarktested(flipfaces[0]) && !smarktested(flipfaces[1])) {
          pa = sorg(flipfaces[0]);
          pb = sdest(flipfaces[0]);
          pc = sapex(flipfaces[0]);
          pd = sapex(flipfaces[1]);
          calculateabovepoint4(pa, pb, pc, pd);
          // Check if a 2-to-2 flip is possible.
          ori1 = orient3d(pc, pd, dummypoint, pa);
          ori2 = orient3d(pc, pd, dummypoint, pb);
          if (ori1 * ori2 < 0) {
            // A 2-to-2 flip is found.
            flip22(flipfaces, lawson, 0);
            // The i-th edge is flipped. The i-th and (i-1)-th subfaces are
            //   changed. The 'flipfaces[1]' contains p as its apex.
            senext2(flipfaces[1], *parentsh);
            // Save the new subface.
            caveshbdlist->newindex((void **) &parysh);
            *parysh = flipfaces[0];
            break;
          }
        } //
      } // i

      if (i == caveshlist->objects) {
        // Do a flip22 and a flip31 to remove p.
        parysh = (face *) fastlookup(caveshlist, 0);
        flipfaces[0] = *parysh;
        spivot(flipfaces[0], flipfaces[1]);
        if (sorg(flipfaces[0]) != sdest(flipfaces[1])) {
          sesymself(flipfaces[1]);
        }
        flip22(flipfaces, lawson, 0);
        senext2(flipfaces[1], *parentsh);
        // Save the new subface.
        caveshbdlist->newindex((void **) &parysh);
        *parysh = flipfaces[0];
      }

      // The edge list at p are changed.
      caveshlist->restart();
    } // while (1)

  } // it

  cavesegshlist->restart();

  if (b->verbose > 2) {
    printf("      Created %ld new subfaces.\n", caveshbdlist->objects);
  }


  if (lawson) {
    lawsonflip();
  }

  return 0;
}

//============================================================================//
//                                                                            //
// slocate()    Locate a point in a surface triangulation.                    //
//                                                                            //
// Staring the search from 'searchsh'(it should not be NULL). Perform a line  //
// walk search for a subface containing the point (p).                        //
//                                                                            //
// If 'aflag' is set, the 'dummypoint' is pre-calculated so that it lies      //
// above the 'searchsh' in its current orientation. The test if c is CCW to   //
// the line a->b can be done by the test if c is below the oriented plane     //
// a->b->dummypoint.                                                          //
//                                                                            //
// If 'cflag' is not TRUE, the triangulation may not be convex.  Stop search  //
// when a segment is met and return OUTSIDE.                                  //
//                                                                            //
// If 'rflag' (rounding) is set, after the location of the point is found,    //
// either ONEDGE or ONFACE, round the result using an epsilon.                //
//                                                                            //
// The returned value indicates the following cases:                          //
//   - ONVERTEX, p is the origin of 'searchsh'.                               //
//   - ONEDGE, p lies on the edge of 'searchsh'.                              //
//   - ONFACE, p lies in the interior of 'searchsh'.                          //
//   - OUTSIDE, p lies outside of the triangulation, p is on the left-hand    //
//     side of the edge 'searchsh'(s), i.e., org(s), dest(s), p are CW.       //
//                                                                            //
//============================================================================//

enum TetMeshCore::locateresult TetMeshCore::slocate(point searchpt, 
  face* searchsh, int aflag, int cflag, int rflag)
{
  face neighsh;
  point pa, pb, pc;
  enum locateresult loc;
  enum {MOVE_BC, MOVE_CA} nextmove;
  double ori, ori_bc, ori_ca;
  int i;

  pa = sorg(*searchsh);
  pb = sdest(*searchsh);
  pc = sapex(*searchsh);

  if (!aflag) {
    // No above point is given. Calculate an above point for this facet.
    calculateabovepoint4(pa, pb, pc, searchpt);
  }

  // 'dummypoint' is given. Make sure it is above [a,b,c]
  ori = orient3d(pa, pb, pc, dummypoint);
  if (ori > 0) {
    sesymself(*searchsh); // Reverse the face orientation.
  } else if (ori == 0.0) {
    // This case should not happen theoretically. But... 
    return UNKNOWN; 
  }

  // Find an edge of the face s.t. p lies on its right-hand side (CCW).
  for (i = 0; i < 3; i++) {
    pa = sorg(*searchsh);
    pb = sdest(*searchsh);
    ori = orient3d(pa, pb, dummypoint, searchpt);
    if (ori > 0) break;
    senextself(*searchsh);
  }
  if (i == 3) {
    return UNKNOWN;
  }

  pc = sapex(*searchsh);

  if (pc == searchpt) {
    senext2self(*searchsh);
    return ONVERTEX;
  }

  while (1) {

    ori_bc = orient3d(pb, pc, dummypoint, searchpt);
    ori_ca = orient3d(pc, pa, dummypoint, searchpt);

    if (ori_bc < 0) {
      if (ori_ca < 0) { // (--)
        // Any of the edges is a viable move.
        if (randomnation(2)) {
          nextmove = MOVE_CA;
        } else {
          nextmove = MOVE_BC;
        }
      } else { // (-#)
        // Edge [b, c] is viable.
        nextmove = MOVE_BC;
      }
    } else {
      if (ori_ca < 0) { // (#-)
        // Edge [c, a] is viable.
        nextmove = MOVE_CA;
      } else {
        if (ori_bc > 0) {
          if (ori_ca > 0) { // (++)
            loc = ONFACE;  // Inside [a, b, c].
            break;
          } else { // (+0)
            senext2self(*searchsh); // On edge [c, a].
            loc = ONEDGE;
            break;
          }
        } else { // ori_bc == 0
          if (ori_ca > 0) { // (0+)
            senextself(*searchsh); // On edge [b, c].
            loc = ONEDGE;
            break;
          } else { // (00)
            // p is coincident with vertex c. 
            senext2self(*searchsh);
            return ONVERTEX;
          }
        }
      }
    }

    // Move to the next face.
    if (nextmove == MOVE_BC) {
      senextself(*searchsh);
    } else {
      senext2self(*searchsh);
    }
    if (!cflag) {
      // NON-convex case. Check if we will cross a boundary.
      if (isshsubseg(*searchsh)) {
        return ENCSEGMENT;
      }
    }
    spivot(*searchsh, neighsh);
    if (neighsh.sh == NULL) {
      return OUTSIDE; // A hull edge.
    }
    // Adjust the edge orientation.
    if (sorg(neighsh) != sdest(*searchsh)) {
      sesymself(neighsh);
    }

    // Update the newly discovered face and its endpoints.
    *searchsh = neighsh;
    pa = sorg(*searchsh);
    pb = sdest(*searchsh);
    pc = sapex(*searchsh);

    if (pc == searchpt) {
      senext2self(*searchsh);
      return ONVERTEX;
    }

  } // while (1)

  // assert(loc == ONFACE || loc == ONEDGE);


  if (rflag) {
    // Round the locate result before return.
    double n[3], area_abc, area_abp, area_bcp, area_cap;

    pa = sorg(*searchsh);
    pb = sdest(*searchsh);
    pc = sapex(*searchsh);

    facenormal(pa, pb, pc, n, 1, NULL);
    area_abc = sqrt(dot(n, n));

    facenormal(pb, pc, searchpt, n, 1, NULL);
    area_bcp = sqrt(dot(n, n));
    if ((area_bcp / area_abc) < b->epsilon) {
      area_bcp = 0; // Rounding.
    }

    facenormal(pc, pa, searchpt, n, 1, NULL);
    area_cap = sqrt(dot(n, n));
    if ((area_cap / area_abc) < b->epsilon) {
      area_cap = 0; // Rounding
    }

    if ((loc == ONFACE) || (loc == OUTSIDE)) {
      facenormal(pa, pb, searchpt, n, 1, NULL);
      area_abp = sqrt(dot(n, n));
      if ((area_abp / area_abc) < b->epsilon) {
        area_abp = 0; // Rounding
      }
    } else { // loc == ONEDGE
      area_abp = 0;
    }

    if (area_abp == 0) {
      if (area_bcp == 0) {
        senextself(*searchsh); 
        loc = ONVERTEX; // p is close to b.
      } else {
        if (area_cap == 0) {
          loc = ONVERTEX; // p is close to a.
        } else {
          loc = ONEDGE; // p is on edge [a,b].
        }
      }
    } else if (area_bcp == 0) {
      if (area_cap == 0) {
        senext2self(*searchsh); 
        loc = ONVERTEX; // p is close to c.
      } else {
        senextself(*searchsh);
        loc = ONEDGE; // p is on edge [b,c].
      }
    } else if (area_cap == 0) {
      senext2self(*searchsh);
      loc = ONEDGE; // p is on edge [c,a].
    } else {
      loc = ONFACE; // p is on face [a,b,c].
    }
  } // if (rflag)

  return loc;
}

//============================================================================//
//                                                                            //
// sscoutsegment()    Look for a segment in the surface triangulation.        //
//                                                                            //
// The segment is given by the origin of 'searchsh' and 'endpt'.              //
//                                                                            //
// If an edge in T is found matching this segment, the segment is "locked"    //
// in T at the edge.  Otherwise, flip the first edge in T that the segment    //
// crosses. Continue the search from the flipped face.                        //
//                                                                            //
// This routine uses 'orisent3d' to determine the search direction. It uses   //
// 'dummypoint' as the 'lifted point' in 3d, and it assumes that it (dummy-   //
// point) lies above the 'searchsh' (w.r.t the Right-hand rule).              //
//                                                                            //
//============================================================================//

enum TetMeshCore::interresult TetMeshCore::sscoutsegment(face *searchsh, 
  point endpt, int insertsegflag, int reporterrorflag, int chkencflag)
{
  face flipshs[2], neighsh;
  point startpt, pa, pb, pc, pd;
  enum interresult dir;
  enum {MOVE_AB, MOVE_CA} nextmove;
  double ori_ab, ori_ca, len;

  pc = NULL; // Avoid warnings from MSVC
  // The origin of 'searchsh' is fixed.
  startpt = sorg(*searchsh); 
  nextmove = MOVE_AB; // Avoid compiler warning.

  if (b->verbose > 2) {
    printf("      Scout segment (%d, %d).\n", pointmark(startpt),
           pointmark(endpt));
  }
  len = distance(startpt, endpt);

  // Search an edge in 'searchsh' on the path of this segment.
  while (1) {

    pb = sdest(*searchsh);
    if (pb == endpt) {
      dir = SHAREEDGE; // Found!
      break;
    }

    pc = sapex(*searchsh);
    if (pc == endpt) {
      senext2self(*searchsh);
      sesymself(*searchsh);
      dir = SHAREEDGE; // Found!
      break;
    }


    // Round the results.
    if ((sqrt(triarea(startpt, pb, endpt)) / len) < b->epsilon) {
      ori_ab = 0.0;
    } else {
      ori_ab = orient3d(startpt, pb, dummypoint, endpt);
    }
    if ((sqrt(triarea(pc, startpt, endpt)) / len) < b->epsilon) {
      ori_ca = 0.0;
    } else {
      ori_ca = orient3d(pc, startpt, dummypoint, endpt);
    }

    if (ori_ab < 0) {
      if (ori_ca < 0) { // (--)
        // Both sides are viable moves.
        if (randomnation(2)) {
          nextmove = MOVE_CA;
        } else {
          nextmove = MOVE_AB;
        }
      } else { // (-#)
        nextmove = MOVE_AB;
      }
    } else {
      if (ori_ca < 0) { // (#-)
        nextmove = MOVE_CA;
      } else {
        if (ori_ab > 0) {
          if (ori_ca > 0) { // (++)
            // The segment intersects with edge [b, c].
            dir = ACROSSEDGE;
            break;
          } else { // (+0)
            // The segment collinear with edge [c, a].
            senext2self(*searchsh);
            sesymself(*searchsh); 
            dir = ACROSSVERT;
            break;
          }
        } else {
          if (ori_ca > 0) { // (0+)
            // The segment is collinear with edge [a, b].
            dir = ACROSSVERT;
            break;
          } else { // (00)
            // startpt == endpt. Not possible.
            terminate_tet_core(this, 2);
          }
        }
      }
    }

    // Move 'searchsh' to the next face, keep the origin unchanged.
    if (nextmove == MOVE_AB) {
      if (chkencflag) {
        // Do not cross boundary.
        if (isshsubseg(*searchsh)) {
          return ACROSSEDGE; // ACROSS_SEG
        }
      }
      spivot(*searchsh, neighsh);
      if (neighsh.sh != NULL) {
        if (sorg(neighsh) != pb) sesymself(neighsh);
        senext(neighsh, *searchsh);
      } else {
        // This side (startpt->pb) is outside. It is caused by rounding error.
        // Try the next side, i.e., (pc->startpt).
        senext2(*searchsh, neighsh);
        if (chkencflag) {
          // Do not cross boundary.
          if (isshsubseg(neighsh)) {
            *searchsh = neighsh;
            return ACROSSEDGE; // ACROSS_SEG
          }
        }
        spivotself(neighsh);
        if (sdest(neighsh) != pc) sesymself(neighsh);
        *searchsh = neighsh;
      }
    } else { // MOVE_CA
      senext2(*searchsh, neighsh);
      if (chkencflag) {
        // Do not cross boundary.
        if (isshsubseg(neighsh)) {
          *searchsh = neighsh;
          return ACROSSEDGE; // ACROSS_SEG
        }
      }
      spivotself(neighsh);
      if (neighsh.sh != NULL) {
        if (sdest(neighsh) != pc) sesymself(neighsh);
        *searchsh = neighsh;
      } else {
        // The same reason as above. 
        // Try the next side, i.e., (startpt->pb).
        if (chkencflag) {
          // Do not cross boundary.
          if (isshsubseg(*searchsh)) {
            return ACROSSEDGE; // ACROSS_SEG
          }
        }
        spivot(*searchsh, neighsh);
        if (sorg(neighsh) != pb) sesymself(neighsh);
        senext(neighsh, *searchsh);
      }
    }
  } // while

  if (dir == SHAREEDGE) {
    if (insertsegflag) {
      // Insert the segment into the triangulation.
      face newseg;
      makeshellface(subsegs, &newseg);
      setshvertices(newseg, startpt, endpt, NULL);
      // Set the default segment marker.
      setshellmark(newseg, -1);
      ssbond(*searchsh, newseg);
      spivot(*searchsh, neighsh);
      if (neighsh.sh != NULL) {
        ssbond(neighsh, newseg);
      }
    }
    return dir;
  }

  if (dir == ACROSSVERT) {
    // A point is found collinear with this segment.
    if (reporterrorflag) {
      point pp = sdest(*searchsh);
      printf("PLC Error:  A vertex lies in a segment in facet #%d.\n",
             shellmark(*searchsh));
      printf("  Vertex:  [%d] (%g,%g,%g).\n",pointmark(pp),pp[0],pp[1],pp[2]);
      printf("  Segment: [%d, %d]\n", pointmark(startpt), pointmark(endpt));
    }
    return dir;
  }

  if (dir == ACROSSEDGE) {
    // Edge [b, c] intersects with the segment.
    senext(*searchsh, flipshs[0]);
    if (isshsubseg(flipshs[0])) {
      if (reporterrorflag) {
        double P[3], Q[3], tp = 0, tq = 0;
        linelineint(startpt, endpt, pb, pc, P, Q, &tp, &tq);
        printf("PLC Error:  Two segments intersect at point (%g,%g,%g),", 
               P[0], P[1], P[2]);
        printf(" in facet #%d.\n", shellmark(*searchsh));
        printf("  Segment 1: [%d, %d]\n", pointmark(pb), pointmark(pc));
        printf("  Segment 2: [%d, %d]\n", pointmark(startpt),pointmark(endpt));
      }
      return dir; // ACROSS_SEG
    }
    // Flip edge [b, c], queue unflipped edges (for Delaunay checks).
    spivot(flipshs[0], flipshs[1]);
    if (sorg(flipshs[1]) != sdest(flipshs[0])) sesymself(flipshs[1]);
    flip22(flipshs, 1, 0);
    // The flip may create an inverted triangle, check it.
    pa = sapex(flipshs[1]);
    pb = sapex(flipshs[0]);
    pc = sorg(flipshs[0]);
    pd = sdest(flipshs[0]);
    // Check if pa and pb are on the different sides of [pc, pd]. 
    // Re-use ori_ab, ori_ca for the tests.
    ori_ab = orient3d(pc, pd, dummypoint, pb);
    ori_ca = orient3d(pd, pc, dummypoint, pa);
    if (ori_ab <= 0) {
      flipshpush(&(flipshs[0])); 
    } else if (ori_ca <= 0) {
      flipshpush(&(flipshs[1])); 
    }
    // Set 'searchsh' s.t. its origin is 'startpt'.
    *searchsh = flipshs[0];
  }

  return sscoutsegment(searchsh, endpt, insertsegflag, reporterrorflag, 
                       chkencflag);
}

//============================================================================//
//                                                                            //
// scarveholes()    Remove triangles not in the facet.                        //
//                                                                            //
// This routine re-uses the two global arrays: caveshlist and caveshbdlist.   //
//                                                                            //
//============================================================================//

void TetMeshCore::scarveholes(int holes, double* holelist)
{
  face *parysh, searchsh, neighsh;
  enum locateresult loc;
  int i, j;

  // Get all triangles. Infect unprotected convex hull triangles. 
  smarktest(recentsh);
  caveshlist->newindex((void **) &parysh);
  *parysh = recentsh;
  for (i = 0; i < caveshlist->objects; i++) {
    parysh = (face *) fastlookup(caveshlist, i);
    searchsh = *parysh;
    searchsh.shver = 0;
    for (j = 0; j < 3; j++) {
      spivot(searchsh, neighsh);
      // Is this side on the convex hull?
      if (neighsh.sh != NULL) {
        if (!smarktested(neighsh)) {
          smarktest(neighsh);
          caveshlist->newindex((void **) &parysh);
          *parysh = neighsh;
        }
      } else {
        // A hull side. Check if it is protected by a segment.
        if (!isshsubseg(searchsh)) {
          // Not protected. Save this face.
          if (!sinfected(searchsh)) {
            sinfect(searchsh);
            caveshbdlist->newindex((void **) &parysh);
            *parysh = searchsh;
          }
        }
      }
      senextself(searchsh);
    }
  }

  // Infect the triangles in the holes.
  for (i = 0; i < 3 * holes; i += 3) {
    searchsh = recentsh;
    loc = slocate(&(holelist[i]), &searchsh, 1, 1, 0);
    if (loc != OUTSIDE) {
      sinfect(searchsh);
      caveshbdlist->newindex((void **) &parysh);
      *parysh = searchsh;
    }
  }

  // Find and infect all exterior triangles.
  for (i = 0; i < caveshbdlist->objects; i++) {
    parysh = (face *) fastlookup(caveshbdlist, i);
    searchsh = *parysh;
    searchsh.shver = 0;
    for (j = 0; j < 3; j++) {
      spivot(searchsh, neighsh);
      if (neighsh.sh != NULL) {
        if (!isshsubseg(searchsh)) {
          if (!sinfected(neighsh)) {
            sinfect(neighsh);
            caveshbdlist->newindex((void **) &parysh);
            *parysh = neighsh;
          }
        } else {
          sdissolve(neighsh); // Disconnect a protected face.
        }
      }
      senextself(searchsh);
    }
  }

  // Delete exterior triangles, unmark interior triangles.
  for (i = 0; i < caveshlist->objects; i++) {
    parysh = (face *) fastlookup(caveshlist, i);
    if (sinfected(*parysh)) {
      shellfacedealloc(subfaces, parysh->sh);
    } else {
      sunmarktest(*parysh);
    }
  }

  caveshlist->restart();
  caveshbdlist->restart();
}

//============================================================================//
//                                                                            //
// triangulate()    Create a CDT for the facet.                               //
//                                                                            //
// All vertices of the triangulation have type FACETVERTEX.  The actual type  //
//   of boundary vertices are set by the routine unifysements().              //
//                                                                            //
// All segments created here will have a default marker '-1'.  Some of these  //
//   segments will get their actual marker defined in 'edgemarkerlist'.       //
//                                                                            //
//============================================================================//

int TetMeshCore::triangulate(int shmark, arraypool* ptlist, arraypool* conlist,
                            int holes, double* holelist)
{
  face searchsh, newsh, *parysh; 
  face newseg, *paryseg;
  point pa, pb, pc, *ppt, *cons;
  int iloc;
  int i, j;

  if (b->verbose > 2) {
    printf("      f%d:  %ld vertices, %ld segments", shmark, ptlist->objects,
           conlist->objects);
    if (holes > 0) {
      printf(", %d holes", holes);
    }
    printf(".\n");
  }

  if (ptlist->objects < 2l) {
    // Not a segment or a facet.
    return 1;
  } else if (ptlist->objects == 2l) {
    pa = * (point *) fastlookup(ptlist, 0);
    pb = * (point *) fastlookup(ptlist, 1);
    if (distance(pa, pb) > 0) {
      // It is a single segment.
      makeshellface(subsegs, &newseg);
      setshvertices(newseg, pa, pb, NULL);
      setshellmark(newseg, -1);
    }
    if (pointtype(pa) == VOLVERTEX) {
      setpointtype(pa, FACETVERTEX);
    }
    if (pointtype(pb) == VOLVERTEX) {
      setpointtype(pb, FACETVERTEX);
    }
    return 1;
  } else if (ptlist->objects == 3) {
    pa = * (point *) fastlookup(ptlist, 0);
    pb = * (point *) fastlookup(ptlist, 1);
    pc = * (point *) fastlookup(ptlist, 2);
  } else {
    // Calculate an above point of this facet.
    if (!calculateabovepoint(ptlist, &pa, &pb, &pc)) {
      if (!b->quiet) {
        printf("Warning:  Unable to triangulate facet #%d. Skipped!\n",shmark);
      }
      return 0; // The point set is degenerate.
    }
  }

  // Create an initial triangulation.
  makeshellface(subfaces, &newsh);
  setshvertices(newsh, pa, pb, pc);
  setshellmark(newsh, shmark);
  recentsh = newsh;

  if (pointtype(pa) == VOLVERTEX) {
    setpointtype(pa, FACETVERTEX);
  }
  if (pointtype(pb) == VOLVERTEX) {
    setpointtype(pb, FACETVERTEX);
  }
  if (pointtype(pc) == VOLVERTEX) {
    setpointtype(pc, FACETVERTEX);
  }

  // Are there area constraints?
  if (b->quality && (in->facetconstraintlist != NULL)) {
    for (i = 0; i < in->numberoffacetconstraints; i++) {
      if (shmark == ((int) in->facetconstraintlist[i * 2])) {
        double area = in->facetconstraintlist[i * 2 + 1];
        setareabound(newsh, area);
        break;
      }
    }
  }

  if (ptlist->objects == 3) {
    // The triangulation only has one element.
    for (i = 0; i < 3; i++) {
      makeshellface(subsegs, &newseg);
      setshvertices(newseg, sorg(newsh), sdest(newsh), NULL);
      setshellmark(newseg, -1);
      ssbond(newsh, newseg);
      senextself(newsh);
    }
    return 1;
  }

  // Triangulate the facet. It may not success (due to rounding error, or
  //   incorrect input data), use 'caveencshlist' and 'caveencseglist' are
  //   re-used to store all the newly created subfaces and segments. So we
  //   can clean them if the triangulation is not successful.
  caveencshlist->newindex((void **) &parysh);
  *parysh = newsh;

  // Incrementally build the triangulation.
  pinfect(pa);
  pinfect(pb);
  pinfect(pc);
  for (i = 0; i < ptlist->objects; i++) {
    ppt = (point *) fastlookup(ptlist, i);
    if (!pinfected(*ppt)) {
      searchsh = recentsh; // Start from 'recentsh'.
      iloc = (int) OUTSIDE;
      // Insert the vertex. Use Bowyer-Watson algo. Round the location.
      iloc = sinsertvertex(*ppt, &searchsh, NULL, iloc, 1, 1);
      if (iloc != ((int) ONVERTEX)) { 
        // Point inserted successfully.
        if (pointtype(*ppt) == VOLVERTEX) {
          setpointtype(*ppt, FACETVERTEX);
        }
        // Save the set of new subfaces.
        for (j = 0; j < caveshbdlist->objects; j++) {
          // Get an old subface at edge [a, b].
          parysh = (face *) fastlookup(caveshbdlist, j);
          spivot(*parysh, searchsh); // The new subface [a, b, p].
          // Do not save a deleted new face (degenerated).
          if (searchsh.sh[3] != NULL) {
            caveencshlist->newindex((void **) &parysh);
            *parysh = searchsh;
          }
        }
        // Delete all removed subfaces.
        for (j = 0; j < caveshlist->objects; j++) {
          parysh = (face *) fastlookup(caveshlist, j);
          shellfacedealloc(subfaces, parysh->sh);
        }
        // Clear the global lists.
        caveshbdlist->restart();
        caveshlist->restart();
        cavesegshlist->restart();
      } else {
        // The facet triangulation is failed.
        break;
      }
    } 
  } // i
  puninfect(pa);
  puninfect(pb);
  puninfect(pc);

  if (i < ptlist->objects) {
    //The facet triangulation is failed. Clean the new subfaces.
    //  There is no new segment be created yet.
    if (!b->quiet) {
      printf("Warning:  Fail to triangulate facet #%d. Skipped!\n", shmark);
    }
    for (i = 0; i < caveencshlist->objects; i++) {
      parysh = (face *) fastlookup(caveencshlist, i);
      if (parysh->sh[3] != NULL) {
        shellfacedealloc(subfaces, parysh->sh);
      }
    }
    caveencshlist->restart();
    return 0;
  }

  // Insert the segments.
  for (i = 0; i < conlist->objects; i++) {
    cons = (point *) fastlookup(conlist, i);
    searchsh = recentsh;
    iloc = (int) slocate(cons[0], &searchsh, 1, 1, 0);
    if (iloc != (int) ONVERTEX) {
      // Not found due to roundoff errors. Do a brute-force search.
      bool bflag = false;
      subfaces->traversalinit();
      searchsh.sh = shellfacetraverse(subfaces);
      while (searchsh.sh != NULL) {
        // Only search the subface in the same facet.
        if (shellmark(searchsh) == shmark) {
          if ((point) searchsh.sh[3] == cons[0]) {
            searchsh.shver = 0; bflag = true; //break;
          } else if ((point) searchsh.sh[4] == cons[0]) {
            searchsh.shver = 2; bflag = true; //break;
          } else if ((point) searchsh.sh[5] == cons[0]) {
            searchsh.shver = 4; bflag = true; //break;
          }
        }
        if (bflag) {
          // [2019-12-03] The subface is not guaranteed to be coplanar,
          //    only use "shmark" is not enough.
          point pa = sorg(searchsh);
          point pb = sdest(searchsh);
          point pc = sapex(searchsh);
          double chkori = orient3d(pa, pb, pc, cons[1]);
          if (chkori != 0.0) {
            double len = distance(pa, pb);
            len += distance(pb, pc);
            len += distance(pc, pa);
            len /= 3.0;
            double len3 = len * len * len;
            double eps = fabs(chkori) / len3;
            if (eps < 1e-5) {
              break; // They are almost coplanar.
            }
          } else {
            break;
          }
          bflag = false; // not this subface.
        }
        searchsh.sh = shellfacetraverse(subfaces);
	  }
      //if (searchsh.sh == NULL) {
      //  // Failed to find a subface containing vertex cons[0].
      //}
    }
    if (searchsh.sh != NULL) {
      // Recover the segment. Some edges may be flipped.
      if (sscoutsegment(&searchsh, cons[1], 1, 1, 0) != SHAREEDGE) {
        break; // Fail to recover a segment.
      }
      // Save this newseg.
      sspivot(searchsh, newseg);
      caveencseglist->newindex((void **) &paryseg);
      *paryseg = newseg;
      if (flipstack != NULL) {
        // Recover locally Delaunay edges.
        lawsonflip();
      }
    } else {
      break; // Failed to find a segment.
    }
  } // i

  if (i < conlist->objects) {
    if (!b->quiet) {
      printf("Warning:  Fail to recover a segment in facet #%d. Skipped!\n",
             shmark);
    }
    for (i = 0; i < caveencshlist->objects; i++) {
      parysh = (face *) fastlookup(caveencshlist, i);
      if (parysh->sh[3] != NULL) {
        shellfacedealloc(subfaces, parysh->sh);
      }
    }
    for (i = 0; i < caveencseglist->objects; i++) {
      paryseg = (face *) fastlookup(caveencseglist, i);
      if (paryseg->sh[3] != NULL) {
        shellfacedealloc(subsegs, paryseg->sh);
      }
    }
    caveencshlist->restart();
    caveencseglist->restart();
    return 0;
  }

  // Remove exterior and hole triangles.
  scarveholes(holes, holelist);

  caveencshlist->restart();
  caveencseglist->restart();
  return 1;
}

//============================================================================//
//                                                                            //
// unifysegments()    Remove redundant segments and create face links.        //
//                                                                            //
// After this routine, although segments are unique, but some of them may be  //
// removed later by mergefacet().  All vertices still have type FACETVERTEX.  //
//                                                                            //
//============================================================================//

void TetMeshCore::unifysegments()
{
  badface *facelink = NULL, *newlinkitem, *f1, *f2;
  face *facperverlist, sface;
  face subsegloop, testseg;
  point torg, tdest;
  double ori1, ori2; //, ori3;
  double n1[3], n2[3];
  double cosang, ang, ang_tol;
  int *idx2faclist;
  int idx, k, m;

  if (b->verbose > 1) {
    printf("  Unifying segments.\n");
  }
  // The limit dihedral angle that two facets are not overlapping.
  //ang_tol = b->facet_overlap_ang_tol / 180.0 * PI;
  //if (ang_tol < 0.0) ang_tol = 0.0;

  // Create a mapping from vertices to subfaces.
  makepoint2submap(subfaces, idx2faclist, facperverlist);


  subsegloop.shver = 0;
  subsegs->traversalinit();
  subsegloop.sh = shellfacetraverse(subsegs);
  while (subsegloop.sh != (shellface *) NULL) {
    torg = sorg(subsegloop);
    tdest = sdest(subsegloop);

    idx = pointmark(torg) - in->firstnumber;
    // Loop through the set of subfaces containing 'torg'.  Get all the
    //   subfaces containing the edge (torg, tdest). Save and order them
    //   in 'sfacelist', the ordering is defined by the right-hand rule
    //   with thumb points from torg to tdest.
    for (k = idx2faclist[idx]; k < idx2faclist[idx + 1]; k++) {
      sface = facperverlist[k];
      // The face may be deleted if it is a duplicated face.
      if (sface.sh[3] == NULL) continue;
      // Search the edge torg->tdest.
      if (sdest(sface) != tdest) {
        senext2self(sface);
        sesymself(sface);
      }
      if (sdest(sface) != tdest) continue;

      // Save the face f in facelink.
      if (flippool->items >= 2) {
        f1 = facelink;
        for (m = 0; m < flippool->items - 1; m++) {
          f2 = f1->nextitem;
          ori1 = facedihedral(torg, tdest, sapex(f1->ss), sapex(f2->ss));
          ori2 = facedihedral(torg, tdest, sapex(f1->ss), sapex(sface));
          if (ori1 >= ori2) {
            break; // insert this face between f1 and f2.
          }
          // Go to the next item;
          f1 = f2;
        } // for (m = 0; ...)
        //if (sface.sh[3] != NULL) {
          // Insert sface between f1 and f2.
          newlinkitem = (badface *) flippool->alloc();
          newlinkitem->ss = sface;
          newlinkitem->nextitem = f1->nextitem;
          f1->nextitem = newlinkitem;
        //}
      } else if (flippool->items == 1) {
        f1 = facelink;
        // Add this face to link if it is not deleted.
        //if (sface.sh[3] != NULL) {
          // Add this face into link.
          newlinkitem = (badface *) flippool->alloc();
          newlinkitem->ss = sface;
          newlinkitem->nextitem = NULL;
          f1->nextitem = newlinkitem;
        //}
      } else {
        // The first face.
        newlinkitem = (badface *) flippool->alloc();
        newlinkitem->ss = sface;
        newlinkitem->nextitem = NULL;
        facelink = newlinkitem;
      }
    } // for (k = idx2faclist[idx]; ...)


    // Set the connection between this segment and faces containing it,
    //   at the same time, remove redundant segments.
    f1 = facelink;
    for (k = 0; k < flippool->items; k++) {
      sspivot(f1->ss, testseg);
      // If 'testseg' is not 'subsegloop' and is not dead, it is redundant.
      if ((testseg.sh != subsegloop.sh) && (testseg.sh[3] != NULL)) {
        shellfacedealloc(subsegs, testseg.sh);
      }
      // Bonds the subface and the segment together.
      ssbond(f1->ss, subsegloop);
      f1 = f1->nextitem;
    }

    // Create the face ring at the segment.
    if (flippool->items > 1) {
      f1 = facelink;
      for (k = 1; k <= flippool->items; k++) {
        k < flippool->items ? f2 = f1->nextitem : f2 = facelink;
        // Calculate the dihedral angle between the two facet.
        facenormal(torg, tdest, sapex(f1->ss), n1, 1, NULL);
        facenormal(torg, tdest, sapex(f2->ss), n2, 1, NULL);
        cosang = dot(n1, n2) / (sqrt(dot(n1, n1)) * sqrt(dot(n2, n2)));
        // Rounding.
        if (cosang > 1.0) cosang = 1.0;
        else if (cosang < -1.0) cosang = -1.0;
        ang = acos(cosang);
        //if (ang < ang_tol) {
        //  // Two facets are treated as overlapping each other.
        //  report_overlapping_facets(&(f1->ss), &(f2->ss), ang);
        //} else {
          // Record the smallest input dihedral angle.
          if (ang < minfacetdihed) {
            minfacetdihed = ang;
          }
          sbond1(f1->ss, f2->ss);
        //}
        f1 = f2;
      }
    }

    flippool->restart();

    // Are there length constraints?
    if (b->quality && (in->segmentconstraintlist != (double *) NULL)) {
      int e1, e2;
      double len;
      for (k = 0; k < in->numberofsegmentconstraints; k++) {
        e1 = (int) in->segmentconstraintlist[k * 3];
        e2 = (int) in->segmentconstraintlist[k * 3 + 1];
        if (((pointmark(torg) == e1) && (pointmark(tdest) == e2)) ||
            ((pointmark(torg) == e2) && (pointmark(tdest) == e1))) {
          len = in->segmentconstraintlist[k * 3 + 2];
          setareabound(subsegloop, len);
          break;
        }
      }
    }

    subsegloop.sh = shellfacetraverse(subsegs);
  }

  delete [] idx2faclist;
  delete [] facperverlist;
}

//============================================================================//
//                                                                            //
// identifyinputedges()    Identify input edges.                              //
//                                                                            //
// A set of input edges is provided in the 'in->edgelist'.  We find these     //
// edges in the surface mesh and make them segments of the mesh.              //
//                                                                            //
// It is possible that an input edge is not in any facet, i.e.,it is a float- //
// segment inside the volume.                                                 //
//                                                                            //
//============================================================================//

void TetMeshCore::identifyinputedges(point *idx2verlist)
{
  face* shperverlist;
  int* idx2shlist;
  face searchsh, neighsh;
  face segloop, checkseg, newseg;
  point checkpt, pa = NULL, pb = NULL;
  int *endpts;
  int edgemarker;
  int idx, i, j;

  int e1, e2;
  double len;

  if (!b->quiet) {
    printf("Inserting edges ...\n");
  }

  // Construct a map from points to subfaces.
  makepoint2submap(subfaces, idx2shlist, shperverlist);

  // Process the set of input edges.
  for (i = 0; i < in->numberofedges; i++) {
    endpts = &(in->edgelist[(i << 1)]);
    if (endpts[0] == endpts[1]) {
      if (!b->quiet) {
        printf("Warning:  Edge #%d is degenerated. Skipped.\n", i);
      }
      continue; // Skip a degenerated edge.
    } else if (dupverts > 0l) {
      // Replace duplicated vertices.
      for (j = 0; j < 2; j++) {
        checkpt = idx2verlist[endpts[j]];
        if (pointtype(checkpt) == DUPLICATEDVERTEX) {
          point meshpt = point2ppt(checkpt);
          endpts[j] = pointmark(meshpt);
        }
      }
    }
    // Recall that all existing segments have a default marker '-1'. 
    // We assign all identified segments a default marker '-2'.
    edgemarker = in->edgemarkerlist ? in->edgemarkerlist[i] : -2;

    // Find a face contains the edge.
    newseg.sh = NULL;
    searchsh.sh = NULL;
    idx = endpts[0] - in->firstnumber;
    for (j = idx2shlist[idx]; j < idx2shlist[idx + 1]; j++) {
      checkpt = sdest(shperverlist[j]);
      if (pointmark(checkpt) == endpts[1]) {
        searchsh = shperverlist[j];
        break; // Found.
      } else {
        checkpt = sapex(shperverlist[j]);
        if (pointmark(checkpt) == endpts[1]) {
          senext2(shperverlist[j], searchsh);
          sesymself(searchsh);
          break;
        }
      }
    } // j

    if (searchsh.sh != NULL) {
      // Check if this edge is already a segment of the mesh.
      sspivot(searchsh, checkseg);
      if (checkseg.sh != NULL) {
        // This segment already exist.
        newseg = checkseg;
      } else {
        // Create a new segment at this edge.
        pa = sorg(searchsh);
        pb = sdest(searchsh);
        makeshellface(subsegs, &newseg);
        setshvertices(newseg, pa, pb, NULL);
        ssbond(searchsh, newseg);
        spivot(searchsh, neighsh);
        if (neighsh.sh != NULL) {
          ssbond(neighsh, newseg);
        }
      }
    } else {
      // It is a dangling segment (not belong to any facets).
      // Get the two endpoints of this segment.
      pa = idx2verlist[endpts[0]];
      pb = idx2verlist[endpts[1]];
      if (pa == pb) {
        if (!b->quiet) {
          printf("Warning:  Edge #%d is degenerated. Skipped.\n", i);
        }
        continue;
      }
      // Check if segment [a,b] already exists.
      // TODO: Change the brute-force search. Slow!
      point *ppt;
      subsegs->traversalinit();
      segloop.sh = shellfacetraverse(subsegs);
      while (segloop.sh != NULL) {
        ppt = (point *) &(segloop.sh[3]);
        if (((ppt[0] == pa) && (ppt[1] == pb)) ||
            ((ppt[0] == pb) && (ppt[1] == pa))) {
          // Found!
          newseg = segloop; 
          break;
        }
        segloop.sh = shellfacetraverse(subsegs);
      }
      if (newseg.sh == NULL) {
        makeshellface(subsegs, &newseg);
        setshvertices(newseg, pa, pb, NULL);
      }
    }

    setshellmark(newseg, edgemarker);

    if (b->quality && (in->segmentconstraintlist != (double *) NULL)) {
      for (i = 0; i < in->numberofsegmentconstraints; i++) {
        e1 = (int) in->segmentconstraintlist[i * 3];
        e2 = (int) in->segmentconstraintlist[i * 3 + 1];
        if (((pointmark(pa) == e1) && (pointmark(pb) == e2)) ||
            ((pointmark(pa) == e2) && (pointmark(pb) == e1))) {
          len = in->segmentconstraintlist[i * 3 + 2];
          setareabound(newseg, len);
          break;
        }
      }
    }
  } // i

  delete [] shperverlist;
  delete [] idx2shlist;
}

//============================================================================//
//                                                                            //
// mergefacets()    Merge adjacent facets.                                    //
//                                                                            //
//============================================================================//

void TetMeshCore::mergefacets()
{
  face parentsh, neighsh, neineish;
  face segloop;
  point pa, pb, pc, pd;
  double n1[3], n2[3];
  double cosang, cosang_tol;

  // Allocate an array to save calcaulated dihedral angles at segments.
  arraypool *dihedangarray = new arraypool(sizeof(double), 10);
  double *paryang = NULL;

  // First, remove coplanar segments.
  // The dihedral angle bound for two different facets.
  cosang_tol = cos(b->facet_separate_ang_tol / 180.0 * PI);

  subsegs->traversalinit();
  segloop.sh = shellfacetraverse(subsegs);
  while (segloop.sh != (shellface *) NULL) {
    // Only remove a segment if it has a marker '-1'.
    if (shellmark(segloop) != -1) {
      segloop.sh = shellfacetraverse(subsegs);
      continue;
    }
    spivot(segloop, parentsh);
    if (parentsh.sh != NULL) {
      spivot(parentsh, neighsh);
      if (neighsh.sh != NULL) {
        spivot(neighsh, neineish);
        if (neineish.sh == parentsh.sh) {
          // Exactly two subfaces at this segment.
          // Only merge them if they have the same boundary marker.
          if (shellmark(parentsh) == shellmark(neighsh)) {
            pa = sorg(segloop);
            pb = sdest(segloop);
            pc = sapex(parentsh);
            pd = sapex(neighsh);
            // Calculate the dihedral angle at the segment [a,b].
            facenormal(pa, pb, pc, n1, 1, NULL);
            facenormal(pa, pb, pd, n2, 1, NULL);
            cosang = dot(n1, n2) / (sqrt(dot(n1, n1)) * sqrt(dot(n2, n2)));
            if (cosang < cosang_tol) {
              ssdissolve(parentsh);
              ssdissolve(neighsh);
              shellfacedealloc(subsegs, segloop.sh);
              // Add the edge to flip stack.
              flipshpush(&parentsh);
            } else {
              // Save 'cosang' to avoid re-calculate it.
              // Re-use the pointer at the first segment.
              dihedangarray->newindex((void **) &paryang);
              *paryang = cosang;
              segloop.sh[6] = (shellface) paryang;
            }
          }
        } // if (neineish.sh == parentsh.sh)
      }
    }
    segloop.sh = shellfacetraverse(subsegs);
  }

  // Second, remove ridge segments at small angles.
  // The dihedral angle bound for two different facets.
  cosang_tol = cos(b->facet_small_ang_tol / 180.0 * PI);
  double cosang_sep_tol = cos((b->facet_separate_ang_tol - 5.0) / 180.0 * PI);
  face shloop;
  face seg1, seg2;
  double cosang1, cosang2;
  int i, j;

  subfaces->traversalinit();
  shloop.sh = shellfacetraverse(subfaces);
  while (shloop.sh != (shellface *) NULL) {
    for (i = 0; i < 3; i++) {
      if (isshsubseg(shloop)) {
        senext(shloop, neighsh);
        if (isshsubseg(neighsh)) {
          // Found two segments sharing at one vertex.
          // Check if they form a small angle.
          pa = sorg(shloop);
          pb = sdest(shloop);
          pc = sapex(shloop);
          for (j = 0; j < 3; j++) n1[j] = pa[j] - pb[j];
          for (j = 0; j < 3; j++) n2[j] = pc[j] - pb[j];
          cosang = dot(n1, n2) / (sqrt(dot(n1, n1)) * sqrt(dot(n2, n2)));
          if (cosang > cosang_tol) {
            // Found a small angle.
            segloop.sh = NULL;
            sspivot(shloop, seg1);
            sspivot(neighsh, seg2);
            if (seg1.sh[6] != NULL) {
              paryang = (double *) (seg1.sh[6]);
              cosang1 = *paryang;
            } else {
              cosang1 = 1.0; // 0 degree;
            }
            if (seg2.sh[6] != NULL) {
              paryang = (double *) (seg2.sh[6]);
              cosang2 = *paryang;
            } else {
              cosang2 = 1.0; // 0 degree;
            }
            if (cosang1 < cosang_sep_tol) {
              if (cosang2 < cosang_sep_tol) {
                if (cosang1 < cosang2) {
                  segloop = seg1;
                } else {
                  segloop = seg2;
                }
              } else {
                segloop = seg1;
              }
            } else {
              if (cosang2 < cosang_sep_tol) {
                segloop = seg2;
              }
            }
            if (segloop.sh != NULL) {
              // Remove this segment.
              segloop.shver = 0;
              spivot(segloop, parentsh);
              spivot(parentsh, neighsh);
              ssdissolve(parentsh);
              ssdissolve(neighsh);
              shellfacedealloc(subsegs, segloop.sh);
              // Add the edge to flip stack.
              flipshpush(&parentsh);
              break;
            }
          }
        } // if (isshsubseg)
      } // if (isshsubseg)
      senextself(shloop);
    }
    shloop.sh = shellfacetraverse(subfaces);
  }

  delete dihedangarray;


  if (flipstack != NULL) {
    lawsonflip(); // Recover Delaunayness.
  }
}

//============================================================================//
//                                                                            //
// meshsurface()    Create a surface mesh of the input PLC.                   //
//                                                                            //
//============================================================================//

void TetMeshCore::meshsurface()
{
  arraypool *ptlist, *conlist;
  point *idx2verlist;
  point tstart, tend, *pnewpt, *cons;
  TetMeshData::facet *f;
  TetMeshData::polygon *p;
  int end1, end2;
  int shmark, i, j;

  if (!b->quiet) {
    printf("Creating surface mesh ...\n");
  }

  // Create a map from indices to points.
  makeindex2pointmap(idx2verlist);

  // Initialize arrays (block size: 2^8 = 256).
  ptlist = new arraypool(sizeof(point *), 8);
  conlist = new arraypool(2 * sizeof(point *), 8);

  // Loop the facet list, triangulate each facet.
  for (shmark = 1; shmark <= in->numberoffacets; shmark++) {

    // Get a facet F.
    f = &in->facetlist[shmark - 1];

    // Process the duplicated points first, they are marked with type
    //   DUPLICATEDVERTEX.  If p and q are duplicated, and p'index > q's,
    //   then p is substituted by q.
    if (dupverts > 0l) {
      // Loop all polygons of this facet.
      for (i = 0; i < f->numberofpolygons; i++) {
        p = &(f->polygonlist[i]);
        // Loop other vertices of this polygon.
        for (j = 0; j < p->numberofvertices; j++) {
          end1 = p->vertexlist[j];
          tstart = idx2verlist[end1];
          if (pointtype(tstart) == DUPLICATEDVERTEX) {
            // Reset the index of vertex-j.
            tend = point2ppt(tstart);
            end2 = pointmark(tend);
            p->vertexlist[j] = end2;
          }
        }
      }
    }

    // Loop polygons of F, get the set of vertices and segments.
    for (i = 0; i < f->numberofpolygons; i++) {
      // Get a polygon.
      p = &(f->polygonlist[i]);
      // Get the first vertex.
      end1 = p->vertexlist[0];
      if ((end1 < in->firstnumber) || 
          (end1 >= in->firstnumber + in->numberofpoints)) {
        if (!b->quiet) {
          printf("Warning:  Invalid the 1st vertex %d of polygon", end1);
          printf(" %d in facet %d.\n", i + 1, shmark);
        }
        continue; // Skip this polygon.
      }
      tstart = idx2verlist[end1];
      // Add tstart to V if it haven't been added yet.
      if (!pinfected(tstart)) {
        pinfect(tstart);
        ptlist->newindex((void **) &pnewpt);
        *pnewpt = tstart;
      }
      // Loop other vertices of this polygon.
      for (j = 1; j <= p->numberofvertices; j++) {
        // get a vertex.
        if (j < p->numberofvertices) {
          end2 = p->vertexlist[j];
        } else {
          end2 = p->vertexlist[0];  // Form a loop from last to first.
        }
        if ((end2 < in->firstnumber) ||
            (end2 >= in->firstnumber + in->numberofpoints)) {
          if (!b->quiet) {
            printf("Warning:  Invalid vertex %d in polygon %d", end2, i + 1);
            printf(" in facet %d.\n", shmark);
          }
        } else {
          if (end1 != end2) {
            // 'end1' and 'end2' form a segment.
            tend = idx2verlist[end2];
            // Add tstart to V if it haven't been added yet.
            if (!pinfected(tend)) {
              pinfect(tend);
              ptlist->newindex((void **) &pnewpt);
              *pnewpt = tend;
            }
            // Save the segment in S (conlist).
            conlist->newindex((void **) &cons);
            cons[0] = tstart;
            cons[1] = tend;
            // Set the start for next continuous segment.
            end1 = end2;
            tstart = tend;
          } else {
            // Two identical vertices mean an isolated vertex of F.
            if (p->numberofvertices > 2) {
              // This may be an error in the input, anyway, we can continue
              //   by simply skipping this segment.
              if (!b->quiet) {
                printf("Warning:  Polygon %d has two identical verts", i + 1);
                printf(" in facet %d.\n", shmark);
              }
            } 
            // Ignore this vertex.
          }
        }
        // Is the polygon degenerate (a segment or a vertex)?
        if (p->numberofvertices == 2) break;
      }
    }
    // Unmark vertices.
    for (i = 0; i < ptlist->objects; i++) {
      pnewpt = (point *) fastlookup(ptlist, i);
      puninfect(*pnewpt);
    }

    // Triangulate F into a CDT.
    // If in->facetmarklist is NULL, use the default marker -1. 
    triangulate(in->facetmarkerlist ? in->facetmarkerlist[shmark - 1] : -1,
                ptlist, conlist, f->numberofholes, f->holelist);

    // Clear working lists.
    ptlist->restart();
    conlist->restart();
  }


  // Remove redundant segments and build the face links.
  unifysegments();
  if (in->numberofedges > 0) {
    // There are input segments. Insert them.
    identifyinputedges(idx2verlist);
  }
  if (!b->diagnose && !b->nomergefacet && !b->nobisect) { // No -d -M -Y
    // Merge coplanar facets.
    mergefacets();
  }

  // Mark all segment vertices to be RIDGEVERTEX.
  face segloop;
  point *ppt;
  subsegs->traversalinit();
  segloop.sh = shellfacetraverse(subsegs);
  while (segloop.sh != NULL) {
    ppt = (point *) &(segloop.sh[3]);
    for (i = 0; i < 2; i++) {
      setpointtype(ppt[i], RIDGEVERTEX);
    }
    segloop.sh = shellfacetraverse(subsegs);
  }

  if (b->object == TetMeshBehavior::STL) {
    // Remove redundant vertices (for .stl input mesh).
    jettisonnodes();
    // Update the number of input vertices.
    in->numberofpoints = points->items;
  }

  if (b->verbose) {
    printf("  %ld (%ld) subfaces (segments).\n", subfaces->items, 
           subsegs->items);
  }

  // The total number of iunput segments.
  insegments = subsegs->items;

  delete [] idx2verlist;
  delete ptlist;
  delete conlist;
}

//                                                                            //
//                                                                            //
//== surface_cxx =============================================================//

} // namespace sqmesh::mesh::tet::detail
