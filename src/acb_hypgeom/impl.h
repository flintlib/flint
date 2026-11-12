/*
    Copyright (C) 2025 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ACB_HYPGEOM_IMPL_H
#define ACB_HYPGEOM_IMPL_H

#include "acb_types.h"

int _mag_gt_norm_ui(const mag_t a, const mag_t b, const mag_t c, ulong n);

void acb_hypgeom_gamma_stirling_choose_param(int * reflect, slong * r, slong * n, const acb_t z, int use_reflect, int digamma, slong prec);
void acb_hypgeom_gamma_stirling_inner(acb_t s, const acb_t z, slong N, slong prec);
slong acb_hypgeom_pfq_choose_n_max(acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, slong prec, slong n_max);
int acb_hypgeom_u_asymp_determine_region(const mag_t r, const mag_t zlo, const acb_t z);
void acb_hypgeom_mag_chi(mag_t chi, ulong n);
void _acb_hypgeom_legendre_q_single(acb_t res, const acb_t n, const acb_t m, const acb_t z, slong prec);
void _acb_hypgeom_legendre_q_double(acb_t res, const acb_t n, const acb_t m, const acb_t z, slong prec);

void acb_poly_reciprocal_majorant(arb_poly_t res, const acb_poly_t poly, slong prec);

/* Functions selected by an integer, for the Bessel integration code */
#define ACB_HYPGEOM_KIND_BESSEL_J 0
#define ACB_HYPGEOM_KIND_BESSEL_Y 1
#define ACB_HYPGEOM_KIND_BESSEL_JY 2
#define ACB_HYPGEOM_KIND_HANKEL_1 3
#define ACB_HYPGEOM_KIND_HANKEL_2 4
#define ACB_HYPGEOM_KIND_BESSEL_I 5
#define ACB_HYPGEOM_KIND_BESSEL_K 6

/* Series and asymptotic expansions, optionally scaled */
int _acb_hypgeom_bessel_asymp_converges(const acb_t nu, const acb_t z, slong prec);
double _acb_hypgeom_bessel_debye_re(double nr, double ni, double wr, double wi);

/* 1F1 by the series and asymptotic expansions only (acb_hypgeom_m falls back
   on numerical integration for real parameters if these are inaccurate) */
void _acb_hypgeom_m_nointegration(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, slong prec);
void _acb_hypgeom_bessel_i_nointegration(acb_t res, const acb_t nu, const acb_t z, int scaled, slong prec);
void _acb_hypgeom_bessel_k_nointegration(acb_t res, const acb_t nu, const acb_t z, int scaled, slong prec);

/* Numerical integration along steepest-descent contours (see sd.h):
   the Bessel function of the given kind (with ACB_HYPGEOM_KIND_BESSEL_JY,
   J and Y in res and res2; otherwise res2 is unused), I and K optionally
   scaled by exp(-+z), and the upper or lower incomplete gamma function.
   The result encloses the function on the input balls.  These return 0
   (leaving res undefined) if the parameters are unsupported or no
   contour is found. */
int _acb_hypgeom_bessel_integration(acb_t res, acb_t res2, int kind, const acb_t nu, const acb_t z, int scaled, slong prec);
int _acb_hypgeom_gamma_incomplete_integration(acb_t res, const acb_t s, const acb_t z, int upper, slong prec);

/* Bounds for sup |F|, sup |dF/dz| and sup |dF/dnu| (F = J, Y, I or K),
   and for sup |F| and sup |dF/ds| (F an incomplete gamma function) over
   the input balls, from the L1 norms of the integrands along the same
   contours; these return 0 on failure (the bounds set to infinity). */
int _acb_hypgeom_bessel_integration_bounds(mag_t M, mag_t Dz, mag_t Dnu, int kind, const acb_t nu, const acb_t z);
int _acb_hypgeom_gamma_incomplete_integration_bounds(mag_t M, mag_t Ds, const acb_t s, const acb_t z, int upper);

/* res (and res2 unless NULL) = the union of f over the parts Im >= 0 and
   Im <= 0 (side = 1, -1) of a ball z crossing the branch cut (-inf, 0) */
int _acb_hypgeom_split_at_cut(acb_t res, acb_t res2, const acb_t nu, const acb_t z, slong prec,
    int (*f)(acb_t, acb_t, const acb_t, const acb_t, int, slong, void *), void * param);

#endif
