// Copyright (c) 2026 Suzhou AI Lab & Shuqian Tech.
// Licensed under the GNU Affero General Public License v3.0 or later.
// See LICENSE file in the project root for details.
//
// Forward declarations of Jonathan Richard Shewchuk's robust geometric
// predicates. The implementation lives in src/core/predicates.cxx and is
// Shewchuk's public-domain code (CMU 1996); see THIRD_PARTY_LICENSES.md §2.
//
// These predicates are shared between the tetrahedral mesher (src/mesh/tet)
// and the auto-CFD surface mesher (src/mesh/auto_cfd). exactinit() must be
// called once before any predicate is invoked.
//
#ifndef SQMESH_CORE_PREDICATES_HPP
#define SQMESH_CORE_PREDICATES_HPP

void exactinit(int verbose, int noexact, int nofilter,
               double maxx, double maxy, double maxz);

double orient2d(double *pa, double *pb, double *pc);
double orient3d(double *pa, double *pb, double *pc, double *pd);
double orient4d(double *pa, double *pb, double *pc, double *pd, double *pe,
                double ah, double bh, double ch, double dh, double eh);

double incircle(double *pa, double *pb, double *pc, double *pd);
double insphere(double *pa, double *pb, double *pc, double *pd, double *pe);

double orient2dexact(double *pa, double *pb, double *pc);
double orient3dexact(double *pa, double *pb, double *pc, double *pd);
double orient4dexact(double *pa, double *pb, double *pc, double *pd, double *pe,
                     double ah, double bh, double ch, double dh, double eh);

#endif // SQMESH_CORE_PREDICATES_HPP
