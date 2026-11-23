/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.2

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ACB_CALC_IMPL_H
#define ACB_CALC_IMPL_H

#include "acb_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* The cache of Gauss-Legendre rules (per thread) used by
   acb_calc_integrate, and shared with other quadrature code: the rules
   with _acb_calc_gl_steps[i] nodes, 0 <= i < ACB_CALC_GL_STEPS, in
   steps of about 2^(1/2) (1, 2, 4, 6, 8, 12, 16, 22, 32, ...). */
#define ACB_CALC_GL_STEPS 38

extern const slong _acb_calc_gl_steps[ACB_CALC_GL_STEPS];

/* Sets *x and *w to the (n + 1) / 2 nonnegative nodes of the rule with
   n = _acb_calc_gl_steps[i] nodes, in decreasing order, and to their
   weights, with a precision of at least prec bits (the entries can have
   a higher precision and should be rounded by the user).  The pointers
   remain valid until the next call for the same i with a higher
   precision, or until flint_cleanup. */
void _acb_calc_gl_nodes(arb_srcptr * x, arb_srcptr * w, slong i, slong prec);

/* Building blocks for the quadrature of exp(f(w)) along polygons
   (gl_exp.c).

   The rules i < ACB_CALC_GL_DFLOAT_STEPS (up to 256 nodes) as arrays of
   dNb balls (N = 1, ..., 4; see dfloat.h): all the n nodes in
   decreasing order, and their weights.  The pointers remain valid
   until flint_cleanup. */
#define ACB_CALC_GL_DFLOAT_STEPS 15
#define ACB_CALC_GL_DFLOAT_NMAX 256
void _acb_calc_gl_nodes_dfloat(const void ** x, const void ** w, int N, slong i);

/* err = (64 / 15) M ah / ((rho - 1) rho^(2 n - 1)), rho = X + sqrt(X^2 - 1):
   a bound for the error of the rule with n nodes on a chord with
   half-length at most ah for an integrand bounded by M on the Bernstein
   ellipse with parameter rho of the chord (with the semi-axes ah X
   along the chord and ah sqrt(X^2 - 1)), or on a set which contains it
   such as the disk of radius ah X around the midpoint; X > 1. */
void _acb_calc_gl_error_bound(mag_t err, const mag_t M, const mag_t ah, double X, slong n);

/* res = the value of the rule i for the integral of
   exp(sum_{j < nc} c_j w^j) over the chord between two complex doubles
   in dNb arithmetic; scratch must
   have room for ACB_CALC_GL_DFLOAT_SCRATCH(N) bytes.  The error of the
   rule is not included. */
#define ACB_CALC_GL_DFLOAT_SCRATCH(N) (ACB_CALC_GL_DFLOAT_NMAX * 12 * ((N) + 1) * sizeof(double))
void _acb_calc_gl_exp_poly_dfloat(acb_t res, acb_srcptr coef, slong nc, double are, double aim,
    double bre, double bim, slong i, int N, void * scratch);

#ifdef __cplusplus
}
#endif

#endif
