/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ACB_MZETA_H
#define ACB_MZETA_H

#include "acb_types.h"
#include "fmpq_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Multivariate zeta functions */

/* Multiple zeta values *****************************************************/

/* zeta(s_1, ..., s_r) = sum_{n_1 > ... > n_r >= 1} n_1^-s_1 ... n_r^-s_r,
   s_1 >= 2 (the convention of PARI/GP zetamult) */
void acb_mzeta_mzv(arb_t res, const ulong * s, slong r, slong prec);

/* star values (n_1 >= ... >= n_r) and Yamamoto's interpolation
   sum_p t^(r - depth(p)) zeta(p) over the merges p of s */
void acb_mzeta_mzv_star(arb_t res, const ulong * s, slong r, slong prec);
void acb_mzeta_mzv_t(acb_t res, const ulong * s, slong r, const acb_t t, slong prec);

/* multiple polylogarithms
   sum_{n_1 > ... > n_r >= 1} x_1^n_1 ... x_r^n_r / (n_1^s_1 ... n_r^s_r)
   (PARI/GP polylogmult), and colored MZVs x_i = exp(2 pi i e_i / N) */
void acb_mzeta_polylog(acb_t res, const ulong * s, acb_srcptr x, slong r, slong prec);
void acb_mzeta_mzv_colored(acb_t res, const ulong * s, const ulong * e, slong r, ulong N, slong prec);

/* all values of weight <= k (or = k), ordered by the index m of PARI/GP:
   the word 0 y 1 of an index has m = (1 y)_2 */
#define ACB_MZETA_MZV_STAR 1
#define ACB_MZETA_MZV_WEIGHT 4

void acb_mzeta_mzv_all(arb_ptr res, ulong k, int flags, slong prec);
ulong acb_mzeta_mzv_index(const ulong * s, slong r);
slong acb_mzeta_mzv_index_composition(ulong * s, ulong m);

/* Cache ********************************************************************/

/* An opaque store of precomputed data (reduced lattices, coefficient
   tables, ...) that algorithms may reuse between calls; any function
   taking a cache accepts NULL.  Not thread-safe: one cache per thread. */
typedef struct
{
    void * entries;
    slong length;
    slong alloc;
}
acb_mzeta_cache_struct;

typedef acb_mzeta_cache_struct acb_mzeta_cache_t[1];

void acb_mzeta_cache_init(acb_mzeta_cache_t cache);
void acb_mzeta_cache_clear(acb_mzeta_cache_t cache);

/* Epstein zeta function ****************************************************/

/* Z(x; y) = sum'_{z in A Z^d} exp(-2 pi i y.z) |x - z|^(-nu)
   (basis vectors = columns of A); with ACB_MZETA_EPSTEIN_REGULARIZED,
   exp(2 pi i x.y) Z(x; y) - s_nu(y) / |det A|. */
#define ACB_MZETA_EPSTEIN_REGULARIZED 1

/* All variants share this interface.  They return 1 when res has been
   computed (it may still be indeterminate, e.g. at a pole), and 0 (with
   res indeterminate) when the algorithm does not apply to the input. */
int acb_mzeta_epstein(acb_t res, const acb_t nu, const arb_mat_t A,
    arb_srcptr x, arb_srcptr y, int flags, slong prec);
int acb_mzeta_epstein_generic(acb_t res, const acb_t nu, const arb_mat_t A,
    arb_srcptr x, arb_srcptr y, int flags, acb_mzeta_cache_t cache, slong prec);
int acb_mzeta_epstein_crandall(acb_t res, const acb_t nu, const arb_mat_t A,
    arb_srcptr x, arb_srcptr y, int flags, acb_mzeta_cache_t cache, slong prec);
int acb_mzeta_epstein_crandall_dfloat(acb_t res, const acb_t nu, const arb_mat_t A,
    arb_srcptr x, arb_srcptr y, int flags, acb_mzeta_cache_t cache, slong prec);
int acb_mzeta_epstein_theta(acb_t res, const acb_t nu, const arb_mat_t A,
    arb_srcptr x, arb_srcptr y, int flags, acb_mzeta_cache_t cache, slong prec);
int acb_mzeta_epstein_shell(acb_t res, const acb_t nu, const arb_mat_t A,
    arb_srcptr x, arb_srcptr y, int flags, acb_mzeta_cache_t cache, slong prec);

/* The same for the lattice with rational Gram matrix G, x = A u and
   y = A^-T v (u, v the coordinates in the basis and the dual basis);
   the algorithm by shells, which needs exact rational input. */
int acb_mzeta_epstein_gram_shell(acb_t res, const acb_t nu, const fmpq_mat_t G,
    const fmpq * u, const fmpq * v, int flags, slong prec);

/* The Madelung constant of NaCl, -Z_{Z^3,1}(0; (1/2,1/2,1/2)) =
   -sum' (-1)^(i+j+k) / sqrt(i^2+j^2+k^2) = 1.7475645946... */
void acb_mzeta_madelung_nacl(arb_t res, slong prec);

/* Generalized incomplete gamma function E_{1-a}(t) = t^(-a) Gamma(a, t)
   for real t > 0, evaluated by this module's own choice between the
   asymptotic series and the convergent series; ga = Gamma(a) or NULL. */
void _acb_mzeta_gamma_upper_scaled(acb_t res, const acb_t a, const arb_t t,
    const acb_t ga, slong prec);

#ifdef __cplusplus
}
#endif

#endif
