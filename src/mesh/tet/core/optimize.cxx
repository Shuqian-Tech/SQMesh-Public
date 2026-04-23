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

//== optimize_cxx ============================================================//
//                                                                            //
//                                                                            //

//============================================================================//
//                                                                            //
// lawsonflip3d()    A three-dimensional Lawson's algorithm.                  //
//                                                                            //
//============================================================================//

long TetMeshCore::lawsonflip3d(flipconstraints *fc)
{
  triface fliptets[5], neightet, hulltet;
  face checksh, casingout;
  badface *popface, *bface;
  point pd, pe, *pts;
  double sign, ori;
  double vol, len3;
  long flipcount, totalcount = 0l;
  long sliver_peels = 0l;
  int t1ver;
  int i;


  while (flippool->items != 0l) {
    if (b->verbose > 2) {
      printf("      Lawson flip %ld faces.\n", flippool->items);
    }
    flipcount = 0l;

    while (flipstack != (badface *) NULL) {
      // Pop a face from the stack.
      popface = flipstack;
      fliptets[0] = popface->tt;
      flipstack = flipstack->nextitem; // The next top item in stack.
      flippool->dealloc((void *) popface);

      // Skip it if it is a dead tet (destroyed by previous flips).
      if (isdeadtet(fliptets[0])) continue;
      // Skip it if it is not the same tet as we saved.
      if (!facemarked(fliptets[0])) continue;

      unmarkface(fliptets[0]);


      if (ishulltet(fliptets[0])) continue;

      fsym(fliptets[0], fliptets[1]);
      if (ishulltet(fliptets[1])) {
        if (nonconvex) {
          // Check if 'fliptets[0]' it is a hull sliver.
          tspivot(fliptets[0], checksh);
          for (i = 0; i < 3; i++) {
            if (!isshsubseg(checksh)) {
              spivot(checksh, casingout);
              //assert(casingout.sh != NULL);
              if (sorg(checksh) != sdest(casingout)) sesymself(casingout);
              stpivot(casingout, neightet);
              if (neightet.tet == fliptets[0].tet) {
                // Found a hull sliver 'neightet'. Let it be [e,d,a,b], where 
                //   [e,d,a] and [d,e,b] are hull faces.
                edestoppo(neightet, hulltet); // [a,b,e,d]
                fsymself(hulltet); // [b,a,e,#]
                if (oppo(hulltet) == dummypoint) {
                  pe = org(neightet);
                  if ((pointtype(pe) == FREEFACETVERTEX) ||
                      (pointtype(pe) == FREESEGVERTEX)) {
                    removevertexbyflips(pe);
                  }
                } else {
                  eorgoppo(neightet, hulltet); // [b,a,d,e]
                  fsymself(hulltet); // [a,b,d,#]
                  if (oppo(hulltet) == dummypoint) {
                    pd = dest(neightet);
                    if ((pointtype(pd) == FREEFACETVERTEX) ||
                        (pointtype(pd) == FREESEGVERTEX)) {
                      removevertexbyflips(pd);
                    }
                  } else {
                    // Perform a 3-to-2 flip to remove the sliver.
                    // To avoid creating an "inverted" subface in the surface
                    //   Check the normals of the two new subfaces, they must
                    //   not be opposite.
                    point chk_pe =  org(neightet);
                    point chk_pd = dest(neightet);
                    point chk_pa = apex(neightet);
                    point chk_pb = oppo(neightet);
                    double n1[3], n2[3];
                    facenormal(chk_pa, chk_pb, chk_pe, n1, 1, NULL);
                    facenormal(chk_pb, chk_pa, chk_pd, n2, 1, NULL);
                    double dot = n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2];
                    if (dot > 0.) {
                      fliptets[0] = neightet;          // [e,d,a,b]
                      fnext(fliptets[0], fliptets[1]); // [e,d,b,c]
                      fnext(fliptets[1], fliptets[2]); // [e,d,c,a]
                      flip32(fliptets, 1, fc);
                      // Update counters.
                      flip32count--;
                      flip22count--;
                      sliver_peels++;
                      if (fc->remove_ndelaunay_edge) {
                        // Update the volume (must be decreased).
                        //assert(fc->tetprism_vol_sum <= 0);
                        tetprism_vol_sum += fc->tetprism_vol_sum;
                        fc->tetprism_vol_sum = 0.0; // Clear it.
                      }
                    } // if (dot. > 0)
                  }
                }
                break;
              } // if (neightet.tet == fliptets[0].tet)
            } // if (!isshsubseg(checksh))
            senextself(checksh);
          } // i
        } // if (nonconvex)
        continue;
      }

      if (checksubfaceflag) {
        // Do not flip if it is a subface.
        if (issubface(fliptets[0])) continue;
      }

      // Test whether the face is locally Delaunay or not.
      pts = (point *) fliptets[1].tet; 
      sign = insphere_s(pts[4], pts[5], pts[6], pts[7], oppo(fliptets[0]));

      if (sign < 0) {
        // A non-Delaunay face. Try to flip it.
        pd = oppo(fliptets[0]);
        pe = oppo(fliptets[1]);

        // Use the length of the edge [d,e] as a reference to determine
        //   a nearly degenerated new tet.
        len3 = distance(pd, pe);
        len3 = (len3 * len3 * len3);
		int round_flag = 0; // [2017-10-20]
        // Check the convexity of its three edges. Stop checking either a
        //   locally non-convex edge (ori < 0) or a flat edge (ori = 0) is
        //   encountered, and 'fliptet' represents that edge.
        for (i = 0; i < 3; i++) {
          ori = orient3d(org(fliptets[0]), dest(fliptets[0]), pd, pe);
          if (ori > 0) {
            // Avoid creating a nearly degenerated new tet at boundary.
            //   Re-use fliptets[2], fliptets[3];
            esym(fliptets[0], fliptets[2]);
            esym(fliptets[1], fliptets[3]);
            if (issubface(fliptets[2]) || issubface(fliptets[3])) {
              vol = orient3dfast(org(fliptets[0]), dest(fliptets[0]), pd, pe);
              if ((fabs(vol) / len3) < b->epsilon) {
                ori = 0.0; // Do rounding.
				round_flag = 1; // [2017-10-20]
              }
            }
          } // Rounding check
          if (ori <= 0) break;
          enextself(fliptets[0]);
          eprevself(fliptets[1]);
        }

        if (ori > 0) {
          // A 2-to-3 flip is found.
          //   [0] [a,b,c,d], 
          //   [1] [b,a,c,e]. no dummypoint.
          flip23(fliptets, 0, fc);
          flipcount++;
          if (fc->remove_ndelaunay_edge) {
            // Update the volume (must be decreased).
            //assert(fc->tetprism_vol_sum <= 0);
            tetprism_vol_sum += fc->tetprism_vol_sum;
            fc->tetprism_vol_sum = 0.0; // Clear it.
          }
          continue;
        } else { // ori <= 0
          // The edge ('fliptets[0]' = [a',b',c',d]) is non-convex or flat,
          //   where the edge [a',b'] is one of [a,b], [b,c], and [c,a].
          if (checksubsegflag) {
            // Do not flip if it is a segment.
            if (issubseg(fliptets[0])) continue;
          }
          // Count the number of interior subfaces for a valid 2-2 flip.
          int scount = 0;
          // Check if there are three or four tets sharing at this edge.        
          esymself(fliptets[0]); // [b,a,d,c]
          for (i = 0; i < 3; i++) {
            if (issubface(fliptets[i])) scount++;
            fnext(fliptets[i], fliptets[i+1]);
          }
          if (fliptets[3].tet == fliptets[0].tet) {
            // A 3-2 flip is found. "scount" must be either 0 or 2.
            if (scount == 1) {
              // This can happen during the boundary recovery. The adjacent
              //   subface is either missing or not recovered yet.
              continue;
            } else if (scount == 2) {
              // Valid if a 2-2 flip is possible.
              for (i = 0; i < 3; i++) {
                if (!issubface(fliptets[i])) break;
              }
              // Assume fliptets[i] is the tet (b,a,c,e). The two subfaces are
              //  fliptets[(i+1)%3] (b,a,e,d) and fliptets[(i+2)%3] (b,a,d,c).
              //  A 2-2 flip is possible if the two faces (d,e,a) and (e,d,b)
              //  are not subfaces.
              triface face1, face2;
              neightet = fliptets[(i+1)%3]; // (b,a,e,d)
              enext(neightet, face1);
              esymself(face1); // (e,a,d)
              eprev(neightet, face2);
              esymself(face2); // (b,e,d)
              if (issubface(face1) || issubface(face2)) {
                continue;
              }
            }
            // A 3-to-2 flip is found. (No hull tet.)
            flip32(fliptets, 0, fc); 
            flipcount++;
            if (fc->remove_ndelaunay_edge) {
              // Update the volume (must be decreased).
              //assert(fc->tetprism_vol_sum <= 0);
              tetprism_vol_sum += fc->tetprism_vol_sum;
              fc->tetprism_vol_sum = 0.0; // Clear it.
            }
            continue;
          } else {
            // There are more than 3 tets at this edge.
            fnext(fliptets[3], fliptets[4]);
            if (fliptets[4].tet == fliptets[0].tet) {
              if (ori != 0.) {
                if (nonconvex) {
                  if (apex(fliptets[3]) == dummypoint) {
                    // This edge is locally non-convex on the hull.
                    // It can be removed by a 4-to-4 flip.
                    ori = 0;
                    round_flag = 1;
                  }
                } // if (nonconvex)
              }
              if (ori == 0) {
                // A 4-to-4 flip is found. (Two hull tets may be involved.)
                // Current tets in 'fliptets':
                //   [0] [b,a,d,c] (d may be newpt)
                //   [1] [b,a,c,e]
                //   [2] [b,a,e,f] (f may be dummypoint)
                //   [3] [b,a,f,d]
                // There are exactly 4 tets at this edge.
                // Moreover, a,b,e,d are coplanar. This 4-4 flip will replace
                //   edge (a,b) to edge (d,e).

                // A valid 2-2 flip is when both faces (a,b,d) and (a,b,e) are
                //   subfaces, and (a,b,c) and (a,b,f) are not subfaces.
                if (issubface(fliptets[0])) { // (a,b,d)
                  if (!issubface(fliptets[2])) { // (a,b,e)
                    continue; // not valid 2-2 flip.
                  }
                  if (issubface(fliptets[1]) ||
                      issubface(fliptets[3])) {
                    continue; // The surface mesh is degnerated.
                  }
                } else {
                  if (issubface(fliptets[1]) ||
                      issubface(fliptets[2]) ||
                      issubface(fliptets[3])) {
                    continue; // not valid 2-2 flip.
                  }
                }
                
                if (round_flag == 1) {
                  //continue; // [2017-10-20]
                  // We want to flip (nearly coplanar) edges [a,b] to [d,e].
                  // Only allow this flip if all new faces are locally Delaunay.
                  // Otherwise, this routine may not terminate.
                  point pb =  org(fliptets[0]);
                  point pa = dest(fliptets[0]);
                  point pc = apex(fliptets[1]);
                  point pf = apex(fliptets[3]); // pf may be dummypoint
                  
                  if (is_collinear_at(pa, pd, pe) ||
                      is_collinear_at(pb, pd, pe)) {
                    continue; // avoid creating a degenerated (sub)face.
                  }

                  // Validate the four new tets (not inverted)
                  double o1, o2;
                  o1 = orient3d(pe, pd, pc, pa);
                  o2 = orient3d(pe, pd, pb, pc);
                  if ((o1 >= 0.) || (o2 >= 0.)) {
                    //assert(0); // to debug...
                    continue; // inverted new tets
                  }
                  if (pf != dummypoint) {
                    double o3, o4;
                    o3 = orient3d(pe, pd, pa, pf);
                    o4 = orient3d(pe, pd, pf, pb);
                    if ((o3 >= 0.) || (o4 >= 0.)) {
                      continue; // inverted new tets
                    }
                  }
                  // Validate locally Delaunay properties of new faces.
                  double test_sign = insphere_s(pe, pd, pc, pa, pb);
                  if (test_sign < 0) {
                    // Locally non-Delaunay. Do not perform the 4-4 flip.
                    continue;
                  }
                  if (pf != dummypoint) {
                    test_sign = insphere_s(pe, pd, pf, pb, pa);
                    if (test_sign < 0) {
                      // Locally non-Delaunay. Do not perform the 4-4 flip.
                      continue;
                    }
                  }
                } // if (round_flag == 1)
                esymself(fliptets[0]); // [a,b,c,d] 
                // A 2-to-3 flip replaces face [a,b,c] by edge [e,d].
                //   This creates a degenerate tet [e,d,a,b] (tmpfliptets[0]).
                //   It will be removed by the followed 3-to-2 flip.
                flip23(fliptets, 0, fc); // No hull tet.
                fnext(fliptets[3], fliptets[1]);
                fnext(fliptets[1], fliptets[2]);
                // Current tets in 'fliptets':
                //   [0] [...]
                //   [1] [b,a,d,e] (degenerated, d may be new point).
                //   [2] [b,a,e,f] (f may be dummypoint)
                //   [3] [b,a,f,d]
                // A 3-to-2 flip replaces edge [b,a] by face [d,e,f].
                //   Hull tets may be involved (f may be dummypoint).
                flip32(&(fliptets[1]), (apex(fliptets[3]) == dummypoint), fc);
                flipcount++;
                flip23count--;
                flip32count--;
                flip44count++;
                if (fc->remove_ndelaunay_edge) {
                  // Update the volume (must be decreased).
                  //assert(fc->tetprism_vol_sum <= 0);
                  tetprism_vol_sum += fc->tetprism_vol_sum;
                  fc->tetprism_vol_sum = 0.0; // Clear it.
                }
                continue;
              } // if (ori == 0)
            }
          }
          // This non-Delaunay face is unflippable. Save it.
          // unflipqueue->newindex((void **) &bface);
          bface = (badface *) flippool->alloc();
          bface->init();
          esymself(fliptets[0]); // *** The original non-Delaunay face ****
          bface->tt = fliptets[0];
          bface->forg  = org(fliptets[0]);
          bface->fdest = dest(fliptets[0]);
          bface->fapex = apex(fliptets[0]);
          // Add it into the unflip queue.
          if (unflip_queue_front == NULL) {
            unflip_queue_front = bface;
          } else {
            unflip_queue_tail->nextitem = bface;
          }
          unflip_queue_tail = bface;
        } // if (ori <= 0)
      } // if (sign < 0)
    } // while (flipstack)

    if (b->verbose > 2) {
      if (flipcount > 0) {
        printf("      Performed %ld flips.\n", flipcount);
      }
      if (flippool->items > 0) {
        printf("      Saved %ld unflippbale faces.\n", flippool->items);
      }
    }
    // Accumulate the counter of flips.
    totalcount += flipcount;

    // Return if no unflippable faces left.
    //if (unflipqueue->objects == 0l) break;
    if (flippool->items == 0l) break;
    // Return if no flip has been performed.
    if (flipcount == 0l) break;

    // Try to flip the unflippable faces.
    while (unflip_queue_front != NULL) {
      bface = unflip_queue_front;
      if (!isdeadtet(bface->tt) &&
        (org(bface->tt) == bface->forg) &&
        (dest(bface->tt) == bface->fdest) &&
        (apex(bface->tt) == bface->fapex)) {
        flippush(flipstack, &(bface->tt));
      }
      unflip_queue_front = bface->nextitem;
      flippool->dealloc((void *) bface);
    }
    unflip_queue_tail = NULL;

  } // while (flippool->items != 0l)

  if (flippool->items > 0l) {
    // Save the unflippable faces to flip them later.
    badface *bf;
    while (unflip_queue_front != NULL) {
      bface = unflip_queue_front;
      if (!isdeadtet(bface->tt) &&
        (org(bface->tt) == bface->forg) &&
        (dest(bface->tt) == bface->fdest) &&
        (apex(bface->tt) == bface->fapex)) {
        //flippush(flipstack, &(bface->tt));
        later_unflip_queue->newindex((void **) &bf);
        *bf = *bface;
      }
      unflip_queue_front = bface->nextitem;
      //flippool->dealloc((void *) bface);
    }
    //unflip_queue_tail = NULL;
    flippool->restart(); // Clear the pool.
  }

  if (b->verbose > 2) {
    if (totalcount > 0) {
      printf("      Performed %ld flips.\n", totalcount);
    }
    if (sliver_peels > 0) {
      printf("      Removed %ld hull slivers.\n", sliver_peels);
    }
    //if (unflipqueue->objects > 0l) {
    //  printf("      %ld unflippable edges remained.\n", unflipqueue->objects);
    //}
  }

  return totalcount + sliver_peels;
}

//============================================================================//
//                                                                            //
// recoverdelaunay()    Recovery the locally Delaunay property.               //
//                                                                            //
//============================================================================//

void TetMeshCore::recoverdelaunay()
{
  badface *bface, *parybface;
  flipconstraints fc;
  int i, j;

  if (b->verbose > 2) {
    printf("    Recovering Delaunayness...\n");
  }
  tetprism_vol_sum = 0.0; // Initialize it.

  if (later_unflip_queue->objects > 0) {
    // Flip the saved unflippable faces.
    for (i = 0; i < later_unflip_queue->objects; i++) {
      bface = (badface *) fastlookup(later_unflip_queue, i);
      if (!isdeadtet(bface->tt) &&
          (org(bface->tt) == bface->forg) &&
          (dest(bface->tt) == bface->fdest) &&
          (apex(bface->tt) == bface->fapex)) {
        flippush(flipstack, &(bface->tt));
      }
    }
    later_unflip_queue->restart(); // clean it.
    if (flippool->items == 0l) {
      return;
    }
  } else {
    if (flippool->items == 0l) {
      // Flip all locally non-Delaunay faces of the tetrahedralisation.
      triface tetloop, neightet; //, *parytet;
      tetrahedrons->traversalinit();
      tetloop.tet = tetrahedrontraverse();
      while (tetloop.tet != NULL) {
        for (tetloop.ver = 0; tetloop.ver < 4; tetloop.ver++) {
          decode(tetloop.tet[tetloop.ver], neightet);
          if (!facemarked(neightet)) {
            flippush(flipstack, &tetloop);
          }
        }
        point *ppt = (point *) &(tetloop.tet[4]);
        tetprism_vol_sum += tetprismvol(ppt[0], ppt[1], ppt[2], ppt[3]);
        tetloop.tet = tetrahedrontraverse();
      }
    }
  }
  
  recover_delaunay_count++;

  // Calulate a relatively lower bound for small improvement.
  //   Used to avoid rounding error in volume calculation.
  fc.bak_tetprism_vol = tetprism_vol_sum * b->epsilon * 1e-3;

  if (b->verbose > 2) {
    printf("    Initial obj = %.17g\n", tetprism_vol_sum);
  }

  if (b->verbose > 2) {
    printf("    Recover Delaunay [Lawson] : %ld\n", flippool->items);
  }

  // First only use the basic Lawson's flip.
  fc.remove_ndelaunay_edge = 1;
  fc.enqflag = 2;

  lawsonflip3d(&fc);

  if (b->verbose > 2) {
    printf("    obj (after Lawson) = %.17g\n", tetprism_vol_sum);
  }

  if (later_unflip_queue->objects == 0l) {
    return;
  }
  
  fc.unflip = 0; // fc.unflip = 1; // Unflip if the edge is not flipped.
  fc.collectnewtets = 1; // new tets are returned in 'cavetetlist'.
  fc.enqflag = 0;

  int bak_autofliplinklevel = autofliplinklevel; 
  int bak_fliplinklevel = b->fliplinklevel;
  autofliplinklevel = 1; // Init level.
  b->fliplinklevel = -1; // No fixed level.

  badface *bfarray = new badface[later_unflip_queue->objects];

  while ((later_unflip_queue->objects > 0) &&
         (autofliplinklevel < 4)) { // level = 1,2,3 //< 10

    int nbf = later_unflip_queue->objects;
    for (i = 0; i < nbf; i++) {
      bfarray[i] = * (badface *) fastlookup(later_unflip_queue, i);
    }
    later_unflip_queue->restart(); // clean it.

    if (b->verbose > 2) {
      printf("    Recover Delaunay [level = %2d] #:  %d.\n",
             autofliplinklevel, nbf);
    }

    for (i = 0; i < nbf; i++) {
      bface = &(bfarray[i]);
      if (getedge(bface->forg, bface->fdest, &bface->tt)) {
        if (removeedgebyflips(&(bface->tt), &fc) == 2) {
          tetprism_vol_sum += fc.tetprism_vol_sum;
        } else {
          // This edge is not removed. Save it in later_flip_queue.
          later_unflip_queue->newindex((void **) &parybface);
          *parybface = bfarray[i]; // *bface;
        }
        fc.tetprism_vol_sum = 0.0; // Clear it.
        if (cavetetlist->objects > 0) {
          // Queue new faces for flips.
          triface neightet, *parytet;
          for (j = 0; j < cavetetlist->objects; j++) {
            parytet = (triface *) fastlookup(cavetetlist, j);
            // A queued new tet may be dead.
            if (!isdeadtet(*parytet)) {
              for (parytet->ver = 0; parytet->ver < 4; parytet->ver++) {
                // Avoid queue a face twice.
                decode(parytet->tet[parytet->ver], neightet);
                if (!facemarked(neightet)) {
                  flippush(flipstack, parytet);
                }
              } // parytet->ver
            }
          } // j
          cavetetlist->restart();
        } // if (cavetetlist->objects > 0)
      }
    } // i

    autofliplinklevel++; // =b->fliplinklevelinc;
  } // while (later_unflip_queue->objects > 0)

  delete [] bfarray;

  if (b->verbose > 2) {
    if (later_unflip_queue->objects > 0l) {
      printf("    %ld non-Delaunay edges remained.\n", later_unflip_queue->objects);
    }
  }
  
  if (flippool->items > 0l) {
    // Flip locally non-Delaunay faces. Unflippable faces are queued
    //   in later_flip_queue.
    fc.remove_ndelaunay_edge = 1;
    fc.enqflag = 2; // queue exteior faces of a flip.
    lawsonflip3d(&fc);
    //fc.enqflag = 0; // for removedgebyflips().
  }

  if (b->verbose > 2) {
    printf("  Final obj  = %.17g\n", tetprism_vol_sum);
  }

  if (later_unflip_queue->objects > 0l) {
    if (b->verbose > 2) {
      printf("    %ld non-Delaunay edges remained.\n", later_unflip_queue->objects);
    }
    later_unflip_queue->restart();
  }

  autofliplinklevel = bak_autofliplinklevel; // Restore this value.
  b->fliplinklevel = bak_fliplinklevel;
}

//============================================================================//
//                                                                            //
// get_seg_laplacian_center()    Get the Laplcian center of a mesh vertex.    //
//                                                                            //
// "mesh_vert" must be a Steiner vertex (FREESEGVERTEX) in a segment.         //
//                                                                            //
//============================================================================//

int TetMeshCore::get_seg_laplacian_center(point mesh_vert, double target[3])
{
  if (pointtype(mesh_vert) == UNUSEDVERTEX) {
    return 0;
  }

  face leftseg, rightseg;

  sdecode(point2sh(mesh_vert), leftseg);
  leftseg.shver = 0;
  if (sdest(leftseg) == mesh_vert) {
    senext(leftseg, rightseg);
    spivotself(rightseg);
    rightseg.shver = 0;
    if (sorg(rightseg) != mesh_vert) {
      sesymself(rightseg);
    }
    if (sorg(rightseg) != mesh_vert) {
      terminate_tet_core(this, 2);
    }
  } else {
    rightseg = leftseg;
    senext2(rightseg, leftseg);
    spivotself(leftseg);
    leftseg.shver = 0;
    if (sdest(leftseg) != mesh_vert) {
      sesymself(leftseg);
    }
    if (sdest(leftseg) != mesh_vert) {
      terminate_tet_core(this, 2);
    }
  }
  point lpt = sorg(leftseg);
  point rpt = sdest(rightseg);

  int j;
  
  for (j = 0; j < 3; j++) {
    target[j] = 0.5 * (lpt[j] + rpt[j]);
  }
  
  return 1;
}

//============================================================================//
//                                                                            //
// get_surf_laplacian_center()    Get the Laplcian center of a mesh vertex.   //
//                                                                            //
// "mesh_vert" must be a Steiner vertex (FREEFACETVERTEX) in a facet.         //
//                                                                            //
//============================================================================//

int TetMeshCore::get_surf_laplacian_center(point mesh_vert, double target[3])
{
  if (pointtype(mesh_vert) == UNUSEDVERTEX) {
    return 0;
  }

  getvertexstar(1, mesh_vert, caveoldtetlist, NULL, caveshlist);

  // The number of vertices is the same as the number of edges.
  int npt = (int) caveshlist->objects;
  int i, j;
  
  for (j = 0; j < 3; j++) {
    target[j] = 0.;
  }

  for (i = 0; i < npt; i++) {
    face *cavesh = (face *) fastlookup(caveshlist, i);
    point e1 =  sorg(*cavesh);
    point e2 = sdest(*cavesh);
    for (j = 0; j < 3; j++) {
      target[j] += e1[j];
    }
    for (j = 0; j < 3; j++) {
      target[j] += e2[j];
    }
  }

  // We added every link vertex twice.
  int npt2 = npt * 2;

  for (j = 0; j < 3; j++) {
    target[j] /= (double) npt2;
  }

  caveoldtetlist->restart();
  caveshlist->restart();
  return 1;
}

//============================================================================//
//                                                                            //
// get_laplacian_center()    Get the Laplcian center of a mesh vertex.        //
//                                                                            //
// "mesh_vert" must be a Steiner vertex (FREEVOLVERTEX) in volume.            //
//                                                                            //
//============================================================================//

int TetMeshCore::get_laplacian_center(point mesh_vert, double target[3])
{
  if (pointtype(mesh_vert) == UNUSEDVERTEX) {
    return 0;
  }
  getvertexstar(1, mesh_vert, caveoldtetlist, cavetetvertlist, NULL);

  // Calculate the laplacian center.
  int npt = (int) cavetetvertlist->objects;
  int i, j;

  for (j = 0; j < 3; j++) {
    target[j] = 0.;
  }

  for (i = 0; i < npt; i++) {
    point *pt = (point *) fastlookup(cavetetvertlist, i);
    for (j = 0; j < 3; j++) {
      target[j] += (*pt)[j];
    }
  }

  for (j = 0; j < 3; j++) {
    target[j] /= (double) npt;
  }

  cavetetvertlist->restart();
  return 1;
}

//============================================================================//
//                                                                            //
// move_vertex()    Try to move a given vertex towards the target position.   //
//                                                                            //
//============================================================================//

bool TetMeshCore::move_vertex(point mesh_vert, double target[3])
{
  if (pointtype(mesh_vert) == UNUSEDVERTEX) {
    if (caveoldtetlist->objects > 0l) {
      caveoldtetlist->restart();
    }
    return 0;
  }
  // Do not move if the target is already very close the vertex.
  if (distance(mesh_vert, target) < minedgelength) {
    if (caveoldtetlist->objects > 0l) {
      caveoldtetlist->restart();
    }
    return 0;
  }
  triface* cavetet;
  point pa, pb, pc;
  double ori;
  int i, j;

  double dir[3], newpos[3];
  double alpha = b->smooth_alpha; // 0.3;

  for (j = 0; j < 3; j++) {
    dir[j] = target[j] - mesh_vert[j];
    newpos[j] = mesh_vert[j] + alpha * dir[j];
  }

  if (caveoldtetlist->objects == 0l) {
    getvertexstar(1, mesh_vert, caveoldtetlist, NULL, NULL);
  }

  bool moveflag = true;
  int iter = 0;

  while (iter < 3) {
    for (i = 0; i < caveoldtetlist->objects; i++) {
      cavetet = (triface *) fastlookup(caveoldtetlist, i);
      if (ishulltet(*cavetet)) continue; // Skip a hull face.

      pa =  org(*cavetet);
      pb = dest(*cavetet);
      pc = apex(*cavetet);
      ori = orient3d(pa, pb, pc, newpos);
      if (ori >= 0) {
        moveflag = false;
        break; // This tet becomes invalid.
      }
    }
    if (moveflag) {
      break;
    } else {
      alpha = (alpha / 2.);
      for (j = 0; j < 3; j++) {
        newpos[j] = mesh_vert[j] + alpha * dir[j];
      }
      iter++;
    }
  } // while (iter < 3)

  if (moveflag) {
    for (j = 0; j < 3; j++) {
      mesh_vert[j] = newpos[j];
    }

    triface checkface, neightet;
    //int j;

    // Push all faces of this vertex star and link into queue.
    for (i = 0; i < caveoldtetlist->objects; i++) {
      cavetet = (triface *) fastlookup(caveoldtetlist, i);
      if (ishulltet(*cavetet)) continue; // Skip a hull face.
      flippush(flipstack, cavetet);
      for (j = 0; j < 3; j++) {
        esym(*cavetet, checkface);
        fsym(checkface, neightet);
        if (!facemarked(neightet)) {
          flippush(flipstack, &checkface);
        }
        enextself(*cavetet);
      }
    }

    if (badtetrahedrons != NULL) {
      // queue all cavity tets for quality check.
      for (i = 0; i < caveoldtetlist->objects; i++) {
        cavetet = (triface *) fastlookup(caveoldtetlist, i);
        if (ishulltet(*cavetet)) continue; // Skip a hull face.
        enqueuetetrahedron(cavetet);
      }
    }

    flipconstraints fc;
    fc.enqflag = 2; // queue all exterior faces of a flip.
    if (badtetrahedrons != NULL) {
      fc.chkencflag = 4; // queue new tets for quality check.
    }
    lawsonflip3d(&fc);


  } // if (moveflag)

  caveoldtetlist->restart();
  return moveflag;
}

//============================================================================//
//                                                                            //
// smooth_vertices()    Smooth vertices.                                      //
//                                                                            //
//============================================================================//

void TetMeshCore::smooth_vertices()
{
  if (!b->quiet) {
    printf("Smoothing vertices...\n");
  }

  if (b->verbose) {
    printf("  Smooth criterion   = %d\n", b->smooth_cirterion);
    printf("  Smooth iterations  = %d\n", b->smooth_maxiter);
    printf("  Smooth relax-alpha = %g\n", b->smooth_alpha);
  }

  point *smpt_list = NULL;
  point *surf_smpt_list = NULL;
  point *seg_smpt_list = NULL;
  int volcount = 0, faccount = 0, segcount = 0;

  // Only use it when we have Steiner points.
  if (st_segref_count > 0) {
    seg_smpt_list = new point[st_segref_count];
  }
  if (st_volref_count > 0) {
    smpt_list = new point[st_volref_count];
  }
  if (st_facref_count > 0) {
    surf_smpt_list = new point[st_facref_count];
  }

  points->traversalinit();
  point ptloop = pointtraverse();
  while (ptloop != NULL) {
    enum verttype vt = pointtype(ptloop);
    if (vt == FREEVOLVERTEX) {
      smpt_list[volcount++] = ptloop;
    } else if (vt == FREEFACETVERTEX) {
      surf_smpt_list[faccount++] = ptloop;
    } else if (vt == FREESEGVERTEX) {
      seg_smpt_list[segcount++] = ptloop;
    }
    ptloop = pointtraverse();
  }

  if ((volcount != st_volref_count) ||
      (faccount != st_facref_count) ||
      (segcount != st_segref_count)) {
    terminate_tet_core(this, 2);
  }

  if (b->verbose > 1) {
    printf("  Smoothing (%ld, %ld) %ld vertices.\n",
           st_segref_count, st_facref_count, st_volref_count);
  }

  // Allocate a list of target points.
  double *target_list = NULL;
  double *surf_target_list = NULL;
  double *seg_target_list = NULL;

  if (st_volref_count > 0) {
    target_list = new double[st_volref_count * 3];
  }
  if (st_facref_count > 0) {
    surf_target_list = new double[st_facref_count * 3];
  }
  if (st_segref_count > 0) {
    seg_target_list = new double[st_segref_count * 3];
  }

  long bak_flipcount = flip23count + flip32count + flip44count;
  int  movedcount = 0, total_movecount = 0;
  int  unmovedcount = 0, total_unmovedcount = 0;
  int  iter = 0, maxiter = b->smooth_maxiter;
  int  i;

  for (iter = 0; iter < maxiter; iter++) {

    movedcount = unmovedcount = 0;

    if (((b->smooth_cirterion & 4) > 0)) { // -s4, -s5, -s6, -s7, default -s3
      //if (st_segref_count > 0) {
      for (i = 0; i < st_segref_count; i++) {
        get_seg_laplacian_center(seg_smpt_list[i], &(seg_target_list[i*3]));
      }
      for (i = 0; i < st_segref_count; i++) {
        if (move_vertex(seg_smpt_list[i], &(seg_target_list[i*3]))) {
          if (later_unflip_queue->objects > b->unflip_queue_limit) {
            recoverdelaunay();
          }
          movedcount++;
        } else {
          unmovedcount++;
        }
      }
      //} // if (st_segref_count > 0)
    }

    if (((b->smooth_cirterion & 2) > 0)) { // default -s3
      //if (st_facref_count > 0) {
      for (i = 0; i < st_facref_count; i++) {
        get_surf_laplacian_center(surf_smpt_list[i], &(surf_target_list[i*3]));
      }
      for (i = 0; i < st_facref_count; i++) {
        if (move_vertex(surf_smpt_list[i], &(surf_target_list[i*3]))) {
          if (later_unflip_queue->objects > b->unflip_queue_limit) {
            recoverdelaunay();
          }
          movedcount++;
        } else {
          unmovedcount++;
        }
      }
      //} // if (st_facref_count > 0)
    }

    if (((b->smooth_cirterion & 1) > 0)) { // default -s3
      //if (st_volref_count > 0) {
      for (i = 0; i < st_volref_count; i++) {
        get_laplacian_center(smpt_list[i], &(target_list[i*3]));
        caveoldtetlist->restart();
      }
      for (i = 0; i < st_volref_count; i++) {
        if (move_vertex(smpt_list[i], &(target_list[i*3]))) {
          if (later_unflip_queue->objects > b->unflip_queue_limit) {
            recoverdelaunay();
          }
          movedcount++;
        } else {
          unmovedcount++;
        }
      }
      //} // if (st_volref_count > 0)
    }

    if (movedcount == 0) break;

    if (b->verbose > 1) {
      printf("  iter=%d, smoothed %d vertices, %d unsmoothed\n", iter,
             movedcount, unmovedcount);
    }


    total_movecount += movedcount;
    total_unmovedcount += unmovedcount;
    
    if (later_unflip_queue->objects > 0) {
      recoverdelaunay();
    }
  } // iter

  if (b->verbose > 1) {
    printf("  Smoothed %d (%d) times, flipped %ld faces.\n",
           total_movecount, total_unmovedcount,
           flip23count + flip32count + flip44count - bak_flipcount);
  }

  if (st_segref_count > 0) {
    delete [] seg_smpt_list;
    delete [] seg_target_list;
  }
  if (st_facref_count > 0) {
    delete [] surf_target_list;
    delete [] surf_smpt_list;
  }
  if (st_volref_count > 0) {
    delete [] target_list;
    delete [] smpt_list;
  }
}

//============================================================================//
//                                                                            //
// get_tet()    Get the tetrahedron with the given vertices.                  //
//                                                                            //
//============================================================================//

bool TetMeshCore::get_tet(point pa, point pb, point pc, point pd, triface *searchtet)
{
  if (getedge(pa, pb, searchtet)) {
    int t1ver;
    triface spintet = *searchtet;
    while (1) {
      if (apex(spintet) == pc) {
        *searchtet = spintet;
        break;
      }
      fnextself(spintet);
      if (spintet.tet == searchtet->tet) break;
    }
    if (apex(*searchtet) == pc) {
      if (oppo(*searchtet) == pd) {
        return true;
      } else {
        fsymself(*searchtet);
        if (oppo(*searchtet) == pd) {
          return true;
        }
      }
    }
  }

  return false;
}

//============================================================================//
//                                                                            //
// get_tetqual()    Calculate various quality measures of a given tetrahedron.//
//                                                                            //
// Calculate the aspect ratio (Lmax / hmin), edge ratio (Lmax/Lmin), maximal  //
// and minimal dihedral angles of this tetrahedron.                           //
//                                                                            //
// These values are returned by:                                              //
//   bf->key,  aspect ratio                                                   //
//   bf->cent[0], cosine of maximal dihedral angle                            //
//   bf->cent[1], cosine of minimal dihedral angle                            //
//   bf->cent[2], edge ratio                                                  //
//   bf->cent[3], minimal edge length                                         //
//   bf->cent[4], volume (used to validate whether it is modified or not).    //
//   bf->cent[5], (no use).                                                   //
//   bf->tet, the edge with maximal dihedral angle.                           //
//   bf->ss.shver, (re-used) count the number of dihedrals > 165 degree.      //
//                                                                            //
//============================================================================//

bool TetMeshCore::get_tetqual(triface *chktet, point oppo_pt, badface *bf)
{
  if (chktet != NULL) {
    bf->init();
    if (oppo_pt == NULL) {
      point *ppt = (point *) &(chktet->tet[4]);
      bf->forg  = ppt[0]; // pa
      bf->fdest = ppt[1]; // pb
      bf->fapex = ppt[2]; // pc
      bf->foppo = ppt[3]; // pd
    } else {
      bf->forg  =  org(*chktet);
      bf->fdest = dest(*chktet);
      bf->fapex = apex(*chktet);
      bf->foppo = oppo_pt;
    }
  }

  double A[4][4], rhs[4], D;
  int indx[4];
  int i, j;

  // get the entries of A[3][3].
  for (i = 0; i < 3; i++) A[0][i] =  bf->forg[i] - bf->foppo[i];  // d->a vec
  for (i = 0; i < 3; i++) A[1][i] = bf->fdest[i] - bf->foppo[i];  // d->b vec
  for (i = 0; i < 3; i++) A[2][i] = bf->fapex[i] - bf->foppo[i];  // d->c vec

  // Get the max-min edge length
  double L[6], Lmax, Lmin;
  double Vab[3], Vbc[3], Vca[3];

  for (i = 0; i < 3; i++) Vab[i] = bf->fdest[i] - bf->forg[i];   // a->b vec
  for (i = 0; i < 3; i++) Vbc[i] = bf->fapex[i] - bf->fdest[i];  // b->c vec
  for (i = 0; i < 3; i++) Vca[i] = bf->forg[i]  - bf->fapex[i];  // c->a vec

  // Use the idx2edge
  L[0] = dot(A[2], A[2]); // edge c,d
  L[1] = dot(A[0], A[0]); // edge a,d
  L[2] = dot(Vab, Vab);   // edge a,b
  L[3] = dot(Vbc, Vbc);   // edge b,c
  L[4] = dot(A[1], A[1]); // edge b,d
  L[5] = dot(Vca, Vca);   // edge a,c

  Lmax = Lmin = L[0];
  //int idx = 0;
  for (i = 1; i < 6; i++) {
    Lmax = (Lmax < L[i] ? L[i] : Lmax);
    Lmin = (Lmin > L[i] ? L[i] : Lmin);
    //if (Lmin > L[i]) {
    //  Lmin = L[i]; idx = i;
    //}
  }

  Lmax = sqrt(Lmax);
  Lmin = sqrt(Lmin);

  // Caluclate the Lmax / Lmin edge ratio (to detect very short edge).
  bf->cent[2] = Lmax / Lmin;
  bf->cent[3] = Lmin;

  // Calculate the normals and heights.
  double N[4][3]; // The normals of the four faces.
  double H[4]; // H[i] is the inverse of the height of its corresponding face.
  bool flat_flag = false;

  if (lu_decmp(A, 3, indx, &D, 0)) {
    // Get the volume of this tet.
    bf->cent[4] = fabs((A[indx[0]][0] * A[indx[1]][1] * A[indx[2]][2]));
    if (bf->cent[4] > 0.0) {
      // Compute the inverse of matrix A, to get 3 normals of the 4 faces.
      for (j = 0; j < 3; j++) {
        for (i = 0; i < 3; i++) rhs[i] = 0.0;
        rhs[j] = 1.0;  // Positive means the inside direction
        lu_solve(A, 3, indx, rhs, 0);
        for (i = 0; i < 3; i++) N[j][i] = rhs[i];
      }
      // Get the fourth normal by summing up the first three.
      for (i = 0; i < 3; i++) N[3][i] = - N[0][i] - N[1][i] - N[2][i];
    } else {
      // This is a very flat tet.
      flat_flag = true;
    }
  } else {
    flat_flag = true;
  }
  
  if (flat_flag) {
    // This tet is nearly degenerate.
    bf->cent[4] = orient3d(bf->fdest, bf->forg, bf->fapex, bf->foppo);
    if (bf->cent[4] <= 0.0) {
      return false; // degenerated or inverted.
    }
    // Calculate the normals of the four faces.
    facenormal(bf->fapex, bf->fdest, bf->foppo, N[0], 1, NULL); // face [c,b,d]
    facenormal(bf->forg,  bf->fapex, bf->foppo, N[1], 1, NULL); // face [a,c,d]
    facenormal(bf->fdest, bf->forg,  bf->foppo, N[2], 1, NULL); // face [b,a,d]
    facenormal(bf->forg,  bf->fdest, bf->fapex, N[3], 1, NULL); // face [a,b,c]
  } // if (!success)

  // Normalized the normals.
  for (i = 0; i < 4; i++) {
    H[i] = sqrt(dot(N[i], N[i]));
    if (H[i] > 0.0) {
      for (j = 0; j < 3; j++) N[i][j] /= H[i];
    } else {
      return false; // H[i] == 0.0;
    }
  }

  if (!flat_flag) {
    // Get the biggest H[i] (corresponding to the smallest height).
    double minheightinv = H[0];
    for (i = 1; i < 4; i++) {
      if (H[i] > minheightinv) minheightinv = H[i];
    }
    // Calulcate the aspect ratio = L_max / h_min.
    bf->key = Lmax * minheightinv;
  } else {
    // A very flat tet.
    //if (bf->key <= 0.0) {
      bf->key = 1.e+30; // infinity.
    //}
  }

  // Calculate the cosine of the dihedral angles of the edges.
  double cosmaxd = 1.0, cosmind = -1.0, cosd;
  int f1, f2, idx = 0;
  bf->ss.shver = 0; // // Count the number of large dihedrals.
  for (i = 0; i < 6; i++) {
    switch (i) {
    case 0: f1 = 0; f2 = 1; break; // [c,d].
    case 1: f1 = 1; f2 = 2; break; // [a,d].
    case 2: f1 = 2; f2 = 3; break; // [a,b].
    case 3: f1 = 0; f2 = 3; break; // [b,c].
    case 4: f1 = 2; f2 = 0; break; // [b,d].
    case 5: f1 = 1; f2 = 3; break; // [a,c].
    }
    cosd = -dot(N[f1], N[f2]);
    if (cosd < -1.0) cosd = -1.0; // Rounding.
    if (cosd >  1.0) cosd =  1.0; // Rounding.
    // cosmaxd = cosd < cosmaxd ? cosd : cosmaxd;
    if (cosd < cosmaxd) {cosmaxd = cosd; idx = i;}
    cosmind = (cosd > cosmind ? cosd : cosmind);
    // Count the number of large dihedrals.
    if (cosd < cos_large_dihed) bf->ss.shver++;
  } // i

  bf->cent[0] = cosmaxd;
  bf->cent[1] = cosmind;

  // Remember the edge with largest dihedral angle.
  if (chktet) bf->tt.tet = chktet->tet;
  bf->tt.ver = edge2ver[idx];

  bf->cent[5] = 0.0;

  return true;
}

bool TetMeshCore::get_tetqual(point pa, point pb, point pc, point pd, badface *bf)
{
  bf->init();

  bf->forg  = pa;
  bf->fdest = pb;
  bf->fapex = pc;
  bf->foppo = pd;

  return get_tetqual(NULL, NULL, bf);
}

//============================================================================//
//                                                                            //
// enqueue_badtet()    Push a bad-quality tet into the proority queue.        //
//                                                                            //
//============================================================================//

void TetMeshCore:: enqueue_badtet(badface *bf)
{
  badface *bt = (badface *) badqual_tets_pool->alloc();

  *bt = *bf;

  // The following vertices are used by get_tet(...) to identify whether
  //   the saved tet is still alive or not.
  //bt->forg  =  org(bf->tt);
  //bt->fdest = dest(bf->tt);
  //bt->fapex = apex(bf->tt);
  //bt->foppo = oppo(bf->tt);

  bt->nextitem = NULL; // important, this pointer is used to recongise the last
                       //   item in each queue.

  // Push it into the priority queue.
  double qual = 1.0 / log(bf->key);

  // Determine the appropriate queue to put the bad subface into.
  int queuenumber = 0;
  if (qual < 1.0) {
    queuenumber = (int) (64.0 * (1.0 - qual));
    if (queuenumber > 63) {
      queuenumber = 63;
    }
  } else {
    // It's not a bad shape; put the subface in the lowest-priority queue.
    queuenumber = 0;
  }
  
  // Are we inserting into an empty queue?
  if (bt_queuefront[queuenumber] == (badface *) NULL) {
    // Yes, we are inserting into an empty queue.
    //   Will this become the highest-priority queue?
    if (queuenumber > bt_firstnonemptyq) {
      // Yes, this is the highest-priority queue.
      bt_nextnonemptyq[queuenumber] = bt_firstnonemptyq;
      bt_firstnonemptyq = queuenumber;
    } else {
      // No, this is not the highest-priority queue.
      //   Find the queue with next higher priority.
      int i = queuenumber + 1;
      while (bt_queuefront[i] == (badface *) NULL) {
        i++;
      }
      // Mark the newly nonempty queue as following a higher-priority queue.
      bt_nextnonemptyq[queuenumber] = bt_nextnonemptyq[i];
      bt_nextnonemptyq[i] = queuenumber;
    }
    // Put the bad subface at the beginning of the (empty) queue.
    bt_queuefront[queuenumber] = bt;
  } else {
    // Add the bad tetrahedron to the end of an already nonempty queue.
    bt_queuetail[queuenumber]->nextitem = bt;
  }
  // Maintain a pointer to the last subface of the queue.
  bt_queuetail[queuenumber] = bt;
}

//============================================================================//
//                                                                            //
// top_badtet()    Get a bad-quality tet from the priority queue.           //
//                                                                            //
//============================================================================//

TetMeshCore::badface* TetMeshCore::top_badtet()
{
  // Keep a record of which queue was accessed in case dequeuebadtetra()
  //   is called later.
  bt_recentq = bt_firstnonemptyq;
  // If no queues are nonempty, return NULL.
  if (bt_firstnonemptyq < 0) {
    return (badface *) NULL;
  } else {
    // Return the first tetrahedron of the highest-priority queue.
    return bt_queuefront[bt_firstnonemptyq];
  }
}

//============================================================================//
//                                                                            //
// dequeue_badtet()    Popup a bad-quality tet from the priority queue.       //
//                                                                            //
//============================================================================//
 
void TetMeshCore::dequeue_badtet()
{
  badface *bt;
  int i;

  // If queues were empty last time topbadtetra() was called, do nothing.
  if (bt_recentq >= 0) {
    // Find the tetrahedron last returned by topbadtetra().
    bt = bt_queuefront[bt_recentq];
    // Remove the tetrahedron from the queue.
    bt_queuefront[bt_recentq] = bt->nextitem;
    // If this queue is now empty, update the list of nonempty queues.
    if (bt == bt_queuetail[bt_recentq]) {
      // Was this the highest-priority queue?
      if (bt_firstnonemptyq == bt_recentq) {
        // Yes; find the queue with next lower priority.
        bt_firstnonemptyq = bt_nextnonemptyq[bt_firstnonemptyq];
      } else {
        // No; find the queue with next higher priority.
        i = bt_recentq + 1;
        while (bt_queuefront[i] == (badface *) NULL) {
          i++;
        }
        bt_nextnonemptyq[i] = bt_nextnonemptyq[bt_recentq];
      }
    }
    // Return the badface to the pool.
    badqual_tets_pool->dealloc((void *) bt);
  }
}


//============================================================================//
//                                                                            //
// add_steinerpt_to_repair()    Add Steiner to repair a bad-qaulity tet.      //
//                                                                            //
//============================================================================//

bool TetMeshCore::add_steinerpt_to_repair(badface *bf, bool bSmooth)
{
  double cosmaxd = bf->cent[0];
  double eta = bf->cent[2];
  int lcount = bf->ss.shver; // the number of large dihedrals.

  triface splittet;
  splittet.tet = NULL;

  if (cosmaxd < cosslidihed) { // cossmtdihed
    // It is a sliver (flat) (might contain a short edge -- skinny).
    triface sliver_edge;
    char shape = 0;

    // Determine the outer shape of this sliver, i.e., a square of a triangle?
    if (lcount == 2) {
      // It is a square. Try to remove the edge [a,b]
      shape = 'S';
      sliver_edge = bf->tt;
    } else if (lcount == 3) {
      // It is a triangle. Try to remove the edge [c,d]
      shape = 'T';
      edestoppo(bf->tt, sliver_edge); // face [c,d,a]
    }

    // Determine a Steiner point according to the shape of this sliver.
    if (shape == 'S') {
      double vol, max_vol = 0.0;

      triface check_sliver = sliver_edge;
      for (int i = 0; i < 2; i++) {
        bool is_bdry = false;
        if (issubseg(check_sliver)) {
          is_bdry = true;
        } else {
          triface spintet = check_sliver;
          int t1ver;
          do {
            if (issubface(spintet)) {
              is_bdry = true; break;
            }
            fnextself(spintet);
          } while (spintet.tet != check_sliver.tet);
        }
        
        if (!is_bdry) {
          triface spintet = check_sliver;
          int t1ver;
          do {
            point *ppt = (point *) &(spintet.tet[4]);
            vol = orient3d(ppt[1], ppt[0], ppt[2], ppt[3]);
            if (vol > max_vol) {
              max_vol = vol;
              splittet = spintet;
            }
            fnextself(spintet);
          } while (spintet.tet != check_sliver.tet);
        }
          
        // Check the opposite edge.
        edestoppoself(check_sliver);
      } // i
    } else if (shape == 'T') {
    }
  } else if (eta > b->opt_max_edge_ratio) {
    // It is a skinny tet.
    // This tet contains a relatively short edge. Check if it can be collapsed.
    double Lmin = bf->cent[3];

    // Get the shortest edge of this tet.
    triface short_edge = bf->tt;
    int i;
    for (i = 0; i < 6; i++) {
      short_edge.ver = edge2ver[i];
      double dd = distance(org(short_edge), dest(short_edge));
      if ((fabs(Lmin - dd) / Lmin) < 1e-4) break;
    }
    if (i == 6) {
      terminate_tet_core(this, 2);
    }

    if (Lmin <= minedgelength) {
      // A very short edge. Check if it was correctly created.
      point e1 = org(short_edge);
      point e2 = dest(short_edge);
      if (issteinerpoint(e1)) {
        if (!create_a_shorter_edge(e1, e2)) {
          terminate_tet_core(this, 2);
        }
      } else if (issteinerpoint(e2)) {
        if (!create_a_shorter_edge(e2, e1)) {
          terminate_tet_core(this, 2);
        }
      }
    }
  }

  if (splittet.tet == NULL) {
    return false; // not added.
  }

  // Do not add if the splittet is also a bad qual tet.
  badface tmpbf;
  if (get_tetqual(&splittet, NULL, &tmpbf)) {
    if (tmpbf.cent[0] < cosslidihed) {
      return false;
    }
  } else {
    return false;
  }

  point steinerpt;
  makepoint(&steinerpt, FREEVOLVERTEX);
  point *ppt = (point *) &(splittet.tet[4]);
  for (int j = 0; j < 3; j++) {
    steinerpt[j] = (ppt[0][j]+ppt[1][j]+ppt[2][j]+ppt[3][j]) / 4.0;
  }

  insertvertexflags ivf;

  //triface searchtet = splittet;
  ivf.iloc = (int) OUTSIDE;
  ivf.bowywat = 3;
  ivf.lawson = 2;
  ivf.rejflag = 0;
  if (badtetrahedrons != NULL) {
    ivf.chkencflag = 4; // queue new tets.
  }
  ivf.sloc = ivf.sbowywat = 0; // No use.
  ivf.splitbdflag = 0; // No use (its an interior vertex).
  ivf.validflag = 1;
  ivf.respectbdflag = 1;

  ivf.smlenflag = 1; // avoid creating very short edges
  ivf.parentpt = NULL; // init.

  if (insertpoint(steinerpt, &splittet, NULL, NULL, &ivf)) {
    st_volref_count++;
    //if (steinerleft > 0) steinerleft--;

    if (flipstack != NULL) {
      flipconstraints fc;
      fc.enqflag = 2;
      if (badtetrahedrons != NULL) {
        fc.chkencflag = 4;
      }
      lawsonflip3d(&fc);
    }

    if (later_unflip_queue->objects > b->unflip_queue_limit) {
      //recoverdelaunay();
      later_unflip_queue->restart(); // clean it.
    }
  } else {
    // Point is not inserted.
    pointdealloc(steinerpt);
    return false;
  }

  if (bSmooth) {
    double ccent[3];
    get_laplacian_center(steinerpt, ccent);
    if (move_vertex(steinerpt, ccent)) {
      opt_smooth_count++;
    }
  } // if (bSmooth)

  if (badtetrahedrons->items > 0) {
    // Push new bad quality tetrahedron into queue.
    badface bf;
    double max_asp = 0., cosmaxd = 1.;
    badtetrahedrons->traversalinit();
    triface *bface = (triface *) badtetrahedrons->traverse();
    while (bface != NULL) {
      if (!isdeadtet(*bface)) {
        // A queued tet may have been processed.
        if (marktest2ed(*bface)) {
          unmarktest2(*bface);
          if (!ishulltet(*bface)) {
            get_tetqual(bface, NULL, &bf);
            // Save the worst quality.
            max_asp = (max_asp > bf.key ? max_asp : bf.key);
            cosmaxd = (cosmaxd < bf.cent[0] ? cosmaxd : bf.cent[0]);
            if ((bf.key > b->opt_max_asp_ratio) || (bf.cent[0] < cosmaxdihed)) {
              bf.forg  =  org(bf.tt);
              bf.fdest = dest(bf.tt);
              bf.fapex = apex(bf.tt);
              bf.foppo = oppo(bf.tt);
              enqueue_badtet(&bf);
            }
          } // if (!ishulltet(*bface))
        }
      }
      bface = (triface *) badtetrahedrons->traverse();
    }
    badtetrahedrons->restart();
  }

  // Check if the bad quality tet is removed or not.
  if (get_tet(bf->forg, bf->fdest, bf->fapex, bf->foppo, &(bf->tt))) {
    // Try to remove it.
    if (repair_tet(bf, true, false, false)) {
      return true;
    }
  } else {
    // This tet is removed.
    return true;
  }

  return false; // not added.
}






//============================================================================//
//                                                                            //
// flip_edge_to_improve()    Flip an edge of a bad-quality tet.               //
//                                                                            //
//============================================================================//

bool TetMeshCore::flip_edge_to_improve(triface *sliver_edge, double& improved_cosmaxd)
{
  if (issubseg(*sliver_edge)) {
    return false;
  }

  flipconstraints fc;

  //fc.noflip_in_surface = 1; // do not flip in surface.
  fc.noflip_in_surface = ((b->nobisect > 0) || ((b->cdtrefine & 2) == 0));
  fc.remove_large_angle = 1;
  fc.unflip = 1;
  fc.collectnewtets = 1;
  fc.checkflipeligibility = 1;
  fc.cosdihed_in = improved_cosmaxd; // cosmaxd;
  fc.cosdihed_out = 0.0; // 90 degree.
  fc.max_asp_out = 0.0;

  if (removeedgebyflips(sliver_edge, &fc) == 2) {
    // This sliver is removed by flips.
    if ((fc.cosdihed_out < cosmaxdihed) || (fc.max_asp_out > b->opt_max_asp_ratio)) {
      // Queue new bad tets for further improvements.
      badface bf;
      for (int j = 0; j < cavetetlist->objects; j++) {
        triface *parytet = (triface *) fastlookup(cavetetlist, j);
        if (!isdeadtet(*parytet) && !ishulltet(*parytet)) {
          if (get_tetqual(parytet, NULL, &bf)) {
            if ((bf.key > b->opt_max_asp_ratio) || (bf.cent[0] < cosmaxdihed)) {
              bf.forg  =  org(bf.tt);
              bf.fdest = dest(bf.tt);
              bf.fapex = apex(bf.tt);
              bf.foppo = oppo(bf.tt);
              enqueue_badtet(&bf);
            }
          } else {
            terminate_tet_core(this, 2);
          }
        }
      } // j
    }
    cavetetlist->restart();
    return true;
  }

  return false;
}

//============================================================================//
//                                                                            //
// repair_tet()    Repair a bad-qaulity tet.                                  //
//                                                                            //
//============================================================================//

bool TetMeshCore::repair_tet(badface *bf, bool bFlips, bool bSmooth, bool bSteiners)
{
  double cosmaxd = bf->cent[0];
  double eta = bf->cent[2];
  int lcount = bf->ss.shver; // the number of large dihedrals.

  if (cosmaxd < cossmtdihed) {
    // It is a sliver (flat) (it might contain a short edge -- skinny).
    //triface sliver_edge;
    char shape = '0';

    // Determine the outer shape of this sliver, i.e., a square of a triangle?
    if (lcount == 2) {
      // It is a square. Try to remove the edge [a,b]
      shape = 'S';
    } else if (lcount == 3) {
      // It is a triangle. Try to remove the edge [c,d]
      shape = 'T';
      //edestoppo(bf->tt, sliver_edge); // face [c,d,a]
    }

    if (bFlips) {
      if (shape == 'S') {
        triface sliver_edge = bf->tt;
        if (flip_edge_to_improve(&sliver_edge, cosmaxd)) {
          opt_flips_count++;
          return true;
        }
        // Due to 'unflip', the flip function may modify the sliver.
        if (get_tet(bf->forg, bf->fdest, bf->fapex, bf->foppo, &(bf->tt))) {
          // Try to flip the opposite edge of this sliver.
          edestoppo(bf->tt, sliver_edge); // face [c,d,a]
          if (flip_edge_to_improve(&sliver_edge, cosmaxd)) {
            opt_flips_count++;
            return true;
          }
        }
      } else if (shape == 'T') {
        triface sliver_edge;
        // flip_face_to_improve(...)
        edestoppo(bf->tt, sliver_edge); // face [c,d,a]
        if (flip_edge_to_improve(&sliver_edge, cosmaxd)) {
          opt_flips_count++;
          return true;
        }
      }
    }
  } else if (eta > b->opt_max_edge_ratio) {
    // It is a skinny tet.
    // This tet contains a relatively short edge. Check if it can be collapsed.
    double Lmin = bf->cent[3];

    // Get the shortest edge of this tet.
    triface short_edge = bf->tt;
    int i;
    for (i = 0; i < 6; i++) {
      short_edge.ver = edge2ver[i];
      double dd = distance(org(short_edge), dest(short_edge));
      //if (fabs(Lmin - dd) < 1e-8) break;
      if ((fabs(Lmin - dd) / Lmin) < 1e-4) break;
    }
    if (i == 6) {
      terminate_tet_core(this, 2);
    }


    if (Lmin <= minedgelength) {
      // A very short edge. Check if it was correctly created.
      point e1 = org(short_edge);
      point e2 = dest(short_edge);
      if (issteinerpoint(e1)) {
        if (!create_a_shorter_edge(e1, e2)) {
          terminate_tet_core(this, 2);
        }
      } else if (issteinerpoint(e2)) {
        if (!create_a_shorter_edge(e2, e1)) {
          terminate_tet_core(this, 2);
        }
      }
    }

  } else {
    // It is neither a flat nor skinny tet. While it has a large asp.

  }


  if (bSteiners &&
      ((bf->key > opt_max_sliver_asp_ratio) || (cosmaxd < cosslidihed))) {
    // This sliver is not removed. Due to 'unflip', the flip function may
    //   modify the sliver.
    if (get_tet(bf->forg, bf->fdest, bf->fapex, bf->foppo, &(bf->tt))) {
      if (add_steinerpt_to_repair(bf, bSmooth)) {
        return true;
      }
    }
  } // if (bSteiners)

  return false; // not repaired
}

//============================================================================//
//                                                                            //
// repair_badqual_tets()    Repair all queued bad quality tet.                //
//                                                                            //
//============================================================================//

long TetMeshCore::repair_badqual_tets(bool bFlips, bool bSmooth, bool bSteiners)
{
  if (b->verbose > 1) {
    printf("  Repairing %ld bad quality tets.\n", badqual_tets_pool->items);
  }
  long repaired_count = 0l;

  while (badqual_tets_pool->items > 0) {
  
    // Get a badtet of highest priority.
    badface *bt = top_badtet();

    if (get_tet(bt->forg, bt->fdest, bt->fapex, bt->foppo, &(bt->tt))) {
      if (repair_tet(bt, bFlips, bSmooth, bSteiners)) {
        repaired_count++;
      } else {
        // Failed to repair this tet. Save it.
        badface *bf = NULL;
        unsplit_badtets->newindex((void **) &bf);
        *bf = *bt;
      }
    } // if (get_tet(...))

    // Return the badtet to the pool.
    dequeue_badtet();
  } // while (badqual_tets_pool->items > 0)

  if (unsplit_badtets->objects > 0l) {
    // Re-initialise the priority queue
    for (int i = 0; i < 64; i++) {
      bt_queuefront[i] = bt_queuetail[i] = NULL;
    }
    bt_firstnonemptyq = -1;
    bt_recentq = -1;
        
    for (int i = 0; i < unsplit_badtets->objects; i++) {
      badface *bt = (badface *) fastlookup(unsplit_badtets, i);
      enqueue_badtet(bt);
    }
    unsplit_badtets->restart();
  } // if (unsplit_badtets->objects > 0l)

  return repaired_count;
}

//============================================================================//
//                                                                            //
// improve_mesh()    Mesh improvement.                                        //
//                                                                            //
//============================================================================//

void TetMeshCore::improve_mesh()
{
  if (!b->quiet) {
    printf("Improving mesh...\n");
  }

  if (b->verbose) {
    printf("  Target maximum aspect ratio = %g.\n", b->opt_max_asp_ratio);
    printf("  Target maximum dihedral angle = %g.\n", b->optmaxdihedral);
    printf("  Maximum flip level   = %d.\n", b->opt_max_flip_level); // -O#
    printf("  Number of iterations = %d.\n", b->opt_iterations); // -O///#
  }

  long blt = b->tetrahedraperblock;
  badqual_tets_pool = new memorypool(sizeof(badface), blt, sizeof(void *), 0);
  badtetrahedrons = new memorypool(sizeof(triface), blt, sizeof(void *), 0);
  unsplit_badtets = new arraypool(sizeof(badface), 10);

  for (int i = 0; i < 64; i++) {
    bt_queuefront[i] = NULL;
  }
  bt_firstnonemptyq = -1;
  bt_recentq = -1;

  cos_large_dihed = cos(135. / 180. * PI); // used in get_tetqual

  cosmaxdihed = cos(b->optmaxdihedral / 180.0 * PI); // set by -o/#

  // The smallest dihedral angle to identify slivers.
  double sliver_ang_tol = b->optmaxdihedral - 5.0;
  if (sliver_ang_tol < 172.0) {
    sliver_ang_tol = 172.;
  }
  cossmtdihed = cos(sliver_ang_tol / 180.0 * PI);

  // The smallest dihedral angle to split slivers.
  double split_sliver_ang_tol = b->optmaxdihedral + 10.0;
  if (split_sliver_ang_tol < 179.0) {
    split_sliver_ang_tol = 179.0;
  } else if (split_sliver_ang_tol > 180.0) {
    split_sliver_ang_tol = 179.9;
  }
  cosslidihed = cos(split_sliver_ang_tol / 180.0 * PI);

  opt_max_sliver_asp_ratio = b->opt_max_asp_ratio * 10.; // set by -o//#

  int attrnum = numelemattrib - 1;
  triface checktet; badface bf;

  // Put all bad tetrahedra into array.
  tetrahedrons->traversalinit();
  checktet.tet = tetrahedrontraverse();
  while (checktet.tet != NULL) {
    if (b->convex) { // -c
      // Skip this tet if it lies in the exterior.
      if (elemattribute(checktet.tet, attrnum) == -1.0) {
        checktet.tet = tetrahedrontraverse();
        continue;
      }
    }
    if (get_tetqual(&checktet, NULL, &bf)) {
      if ((bf.key > b->opt_max_asp_ratio) || (bf.cent[0] < cosmaxdihed)) {
        bf.forg  =  org(bf.tt);
        bf.fdest = dest(bf.tt);
        bf.fapex = apex(bf.tt);
        bf.foppo = oppo(bf.tt);
        enqueue_badtet(&bf);
      }
    } else {
      terminate_tet_core(this, 2); // a degenerated tet.
    }
    checktet.tet = tetrahedrontraverse();
  }

  // Backup flip edge options.
  int bakautofliplinklevel = autofliplinklevel;
  int bakfliplinklevel = b->fliplinklevel;
  int bakmaxflipstarsize = b->flipstarsize;

  b->fliplinklevel = 1; // initial (<= b->opt_max_flip_level, -O#)
  b->flipstarsize = 10;  // b->optmaxflipstarsize;

  long total_repaired_count = 0l;
  long bak_pt_count = points->items;

  // Only using flips.
  while (badqual_tets_pool->items > 0) {
    long repaired_count = repair_badqual_tets(true, false, false);
    total_repaired_count += repaired_count;
    if (b->fliplinklevel < b->opt_max_flip_level) {
      b->fliplinklevel++;
    } else {
      break; // maximal flip level is reached.
    }
  } // while (badqual_tets_pool->items > 0)
  
  if (b->verbose > 1) {
    printf("  Repaired %ld tetrahedra by flips.\n", total_repaired_count);
    printf("  %ld badqual tets remained.\n", badqual_tets_pool->items);
  }

  int iter = 0;
  long bak_st_count = st_volref_count;
  while ((badqual_tets_pool->items > 0) && (iter < b->opt_iterations)) {
    //b->fliplinklevel++;
    long repaired_count = repair_badqual_tets(true, true, true);
    // Break if no repair and no new Steiner point.
    if ((repaired_count == 0l) && (bak_st_count == st_volref_count)) {
      break;
    }
    total_repaired_count += repaired_count;
    bak_st_count = st_volref_count;
    iter++;
  } // while (badqual_tets_pool->items > 0)

  // Do last flips.
  if (badqual_tets_pool->items > 0) {
    long repaired_count = repair_badqual_tets(true, false, false);
    total_repaired_count += repaired_count;
  }

  if (b->verbose > 1) {
    printf("  Repaired %ld tetrahedra.\n", total_repaired_count);
    printf("  %ld badqual tets remained.\n", badqual_tets_pool->items);
  }

  if (later_unflip_queue->objects > b->unflip_queue_limit) {
    //recoverdelaunay();
    later_unflip_queue->restart(); // clean it.
  }

  if (b->verbose) {
    if (opt_flips_count > 0l) {
      printf("  Removed %ld edges/faces.\n", opt_flips_count);
    }
    if (opt_collapse_count > 0l) {
      printf("  Collapsed %ld edges/faces.\n", opt_collapse_count);
    }
    if (opt_smooth_count > 0l) {
      printf("  Smoothed %ld vertices.\n", opt_smooth_count);
    }
    if ((points->items - bak_pt_count) > 0l) {
      printf("  Added %ld Steiner points.\n", points->items - bak_pt_count);
    }
  }


  // Restore original flip edge options.
  autofliplinklevel = bakautofliplinklevel;
  b->fliplinklevel = bakfliplinklevel;
  b->flipstarsize = bakmaxflipstarsize;

  delete badtetrahedrons;
  badtetrahedrons = NULL;
  delete badqual_tets_pool;
  badqual_tets_pool = NULL;
  delete unsplit_badtets;
  unsplit_badtets = NULL;
}

//                                                                            //
//                                                                            //
//== optimize_cxx ============================================================//

} // namespace sqmesh::mesh::tet::detail
