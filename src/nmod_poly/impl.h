/*
    Copyright (C) 2025 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_POLY_IMPL_H
#define NMOD_POLY_IMPL_H

#include "nmod_types.h"
#include "nmod_poly.h"  /* TODO is this ok? */

void _nmod_poly_inv_series_basecase_preinv1(nn_ptr Qinv, nn_srcptr Q, slong Qlen, slong n, ulong q, nmod_t mod);
void _nmod_poly_div_series_basecase_preinv1(nn_ptr Qinv, nn_srcptr P, slong Plen, nn_srcptr Q, slong Qlen, slong n, ulong q, nmod_t mod);
int nmod_poly_irreducible_binomial(nmod_poly_t res, ulong n);
int nmod_poly_irreducible_trinomial(nmod_poly_t res, ulong n);
int nmod_poly_irreducible_tetranomial(nmod_poly_t res, ulong n);
int nmod_poly_irreducible_pentanomial(nmod_poly_t res, ulong n);
void _nmod_poly_divrem_q0_preinv1(nn_ptr Q, nn_ptr R, nn_srcptr A, nn_srcptr B, slong lenA, ulong invL, nmod_t mod);

void _nmod_geometric_progression_evaluate_init_nonfullword(nmod_geometric_progression_t G,
                                                           ulong r, slong len, nmod_t mod,
                                                           ulong q, ulong q_pr_quo, ulong q_pr_rem,
                                                           ulong inv_r,
                                                           ulong inv_q, ulong inv_q_pr_quo, ulong inv_q_pr_rem);
void _nmod_geometric_progression_evaluate_init(nmod_geometric_progression_t G,
                                               ulong r, slong len, nmod_t mod,
                                               ulong q, ulong inv_r, ulong inv_q);
void _nmod_geometric_progression_interpolate_init_nonfullword(nmod_geometric_progression_t G,
                                                              slong len, nmod_t mod,
                                                              ulong q, ulong q_pr_quo, ulong q_pr_rem,
                                                              ulong inv_q, ulong inv_q_pr_quo, ulong inv_q_pr_rem);
void _nmod_geometric_progression_interpolate_init(nmod_geometric_progression_t G, slong len, nmod_t mod,
                                                  ulong q, ulong inv_q);


#endif
