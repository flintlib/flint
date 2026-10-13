/*
    Copyright (C) 2026 Maël Hostettler
    Developed using Claude Opus 5 and Opus 5.5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef GR_EC_IMPL_H
#define GR_EC_IMPL_H

#include "gr_ec.h"
#include "gr_poly.h"

/*
    Variants of the Jacobian doubling and mixed addition taking a caller
    supplied workspace of GR_EC_JAC_SCRATCH base ring elements. Allocating
    and initializing the temporaries costs a noticeable fraction of an
    operation over small moduli, so the scalar multiplication ladder hoists
    the workspace out of its loop and calls these directly.
*/

#define GR_EC_JAC_SCRATCH 9

/*
    Over a ring that is not an integral domain, a nonzero value is not
    enough to make a branch of the group law right: it has to be a unit.
    The functions taking a witness w multiply it by every quantity a branch
    treated as nonzero, so that a non-unit w (over Z/n: gcd(w, n) != 1)
    flags the computation as untrustworthy. A NULL witness is ignored.
*/
static inline WARN_UNUSED_RESULT int
_gr_ec_witness(gr_ptr w, gr_srcptr x, gr_ctx_t R)
{
    return (w == NULL) ? GR_SUCCESS : gr_mul(w, w, x, R);
}

/* Model-specific group law */

WARN_UNUSED_RESULT int _gr_ec_aff_point_add_long_weierstrass(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_aff_point_add_short_weierstrass(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_aff_point_dbl_long_weierstrass(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_aff_point_dbl_short_weierstrass(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_aff_point_mul_fmpz_binary(gr_ec_aff_point_t res, const gr_ec_aff_point_t P, const fmpz_t n, gr_ec_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_ec_jac_point_add_long_weierstrass(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const gr_ec_jac_point_t Q, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_jac_point_add_short_weierstrass(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const gr_ec_jac_point_t Q, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_jac_point_dbl_long_weierstrass(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_jac_point_dbl_short_weierstrass(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_jac_point_add_aff_point_long_weierstrass(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_jac_point_add_aff_point_short_weierstrass(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const gr_ec_aff_point_t Q, gr_ec_ctx_t ctx);

/* Same, with a workspace of GR_EC_JAC_SCRATCH elements and a witness */

WARN_UNUSED_RESULT int _gr_ec_jac_point_dbl_long_weierstrass_ws(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, gr_ptr t, gr_ptr w, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_jac_point_dbl_short_weierstrass_ws(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, gr_ptr t, gr_ptr w, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_jac_point_add_aff_point_long_weierstrass_ws(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const gr_ec_aff_point_t Q, gr_ptr t, gr_ptr w, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_jac_point_add_aff_point_short_weierstrass_ws(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const gr_ec_aff_point_t Q, gr_ptr t, gr_ptr w, gr_ec_ctx_t ctx);

/* Scalar multiplication and batch normalisation */

WARN_UNUSED_RESULT int _gr_ec_jac_point_mul_fmpz_binary(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const fmpz_t n, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_jac_point_mul_fmpz_naf(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const fmpz_t n, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_jac_point_mul_fmpz_naf_witness(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const fmpz_t n, gr_ptr w, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_jac_point_vec_get_aff_point_vec(gr_ec_aff_point_struct * res, const gr_ec_jac_point_struct * P, slong len, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_jac_point_vec_get_aff_point_vec_witness(gr_ec_aff_point_struct * res, const gr_ec_jac_point_struct * P, slong len, gr_ptr w, gr_ec_ctx_t ctx);

/* Context: generic method table (generic.c), printing (ctx.c) and the
   uncached point count (cardinality.c) */

void _gr_ec_ctx_init_methods(gr_ec_ctx_t ctx);
int _gr_ec_ctx_write(gr_stream_t out, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_ctx_cardinality(fmpz_t res, gr_ec_ctx_t ctx);

/*
    The largest field the naive count will walk over. The caller asked for
    the O(q) algorithm, but not for it never to return; the dispatcher
    also reads it, to know when the walk is still an affordable fallback.
*/
#define GR_EC_NAIVE_MAX_Q WORD(100000000)

/*
    Arithmetic on the l-torsion in F_q[x]/(m), shared by the point counting
    algorithms that work there. Defined in torsion.c; see the notes there.
*/

/* a point of the l-torsion as (u, v y), or the identity */
typedef struct
{
    gr_poly_t u;
    gr_poly_t v;
    int is_inf;
}
gr_ec_tors_struct;

typedef struct
{
    gr_ctx_struct * R;
    gr_poly_preinv_t P;         /* preconditioned psi_l */
    gr_poly_t psi;              /* psi_l itself, for the xgcd in _gr_ec_tors_inv */
    gr_poly_t f;                /* x^3 + a4 x + a6 mod psi_l */
    gr_poly_t a4;
    gr_poly_t factor;           /* a proper factor of psi, when one shows up */
    int have_factor;
}
gr_ec_tors_ctx_struct;

void _gr_ec_tors_init(gr_ec_tors_struct * T, gr_ec_tors_ctx_struct * C);
void _gr_ec_tors_clear(gr_ec_tors_struct * T, gr_ec_tors_ctx_struct * C);
WARN_UNUSED_RESULT int _gr_ec_tors_set(gr_ec_tors_struct * D, const gr_ec_tors_struct * S, gr_ec_tors_ctx_struct * C);
WARN_UNUSED_RESULT int _gr_ec_tors_mulmod(gr_poly_t res, const gr_poly_t a, const gr_poly_t b, gr_ec_tors_ctx_struct * C);
WARN_UNUSED_RESULT int _gr_ec_tors_inv(gr_poly_t res, const gr_poly_t a, gr_ec_tors_ctx_struct * C);
truth_t _gr_ec_tors_equal(const gr_ec_tors_struct * A, const gr_ec_tors_struct * B, gr_ec_tors_ctx_struct * C);
WARN_UNUSED_RESULT int _gr_ec_tors_neg(gr_ec_tors_struct * D, const gr_ec_tors_struct * S, gr_ec_tors_ctx_struct * C);
WARN_UNUSED_RESULT int _gr_ec_tors_add(gr_ec_tors_struct * D, const gr_ec_tors_struct * A, const gr_ec_tors_struct * B, gr_ec_tors_ctx_struct * C);
WARN_UNUSED_RESULT int _gr_ec_tors_mul_ui(gr_ec_tors_struct * D, const gr_ec_tors_struct * S, ulong k, gr_ec_tors_ctx_struct * C);
WARN_UNUSED_RESULT int _gr_ec_tors_ctx_init(gr_ec_tors_ctx_struct * C, const gr_poly_t m, gr_ec_ctx_t ctx);
void _gr_ec_tors_ctx_clear(gr_ec_tors_ctx_struct * C);
WARN_UNUSED_RESULT int _gr_ec_tors_frobenius(gr_ec_tors_struct * P, gr_ec_tors_struct * phiP, const fmpz_t q, gr_ec_tors_ctx_struct * C);

/* a question asked of the points of the subgroup C's modulus describes */
typedef int (*_gr_ec_tors_step_t)(ulong * res, ulong l, const fmpz_t q, gr_ec_tors_ctx_struct * C, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_tors_solve(ulong * res, const gr_poly_t m, ulong l, const fmpz_t q, _gr_ec_tors_step_t step, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_cardinality_bsgs_progression(fmpz_t res, const fmpz_t t0, const fmpz_t M, slong max_baby, gr_ec_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_ec_tors_discrete_log(ulong * k, int * found, const gr_ec_tors_struct * target, const gr_ec_tors_struct * base, ulong l, gr_ec_tors_ctx_struct * C);

/*
    Phi_l(X, Y) reduced into R, as l + 2 polynomials in Y with
    Phi_l = sum_a X^a phi[a](Y). l prime, characteristic of R above l + 1.
    Defined in modular_poly.c.
*/
WARN_UNUSED_RESULT int _gr_ec_modular_polynomial(gr_poly_struct * phi, ulong l, gr_ctx_t R);
WARN_UNUSED_RESULT int _gr_ec_modular_polynomial_canonical(gr_poly_struct * phi, ulong * s_out, ulong l, gr_ctx_t R);
WARN_UNUSED_RESULT int _gr_ec_j_qexp(gr_poly_t J, slong len, gr_ctx_t R);

/*
    For an Elkies prime l, the kernel polynomial of a rational l-isogeny,
    of degree (l - 1)/2; *elkies says whether l was one. Defined in
    elkies.c.
*/
WARN_UNUSED_RESULT int _gr_ec_elkies_kernel(gr_poly_t h, int * elkies, ulong * atkin_r, ulong l, const fmpz_t q, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_elkies_trace(ulong * tl, int * elkies, ulong * atkin_r, ulong l, const fmpz_t q, gr_ec_ctx_t ctx);
slong _gr_ec_atkin_candidates(ulong * T, ulong l, ulong r, ulong ql);

/* t mod l is one of T[0], ..., T[n-1] */
typedef struct
{
    ulong l;
    slong n;
    ulong * T;
}
gr_ec_atkin_struct;

double _gr_ec_atkin_plan(int * side, slong * Z1, int * use_atkin, const gr_ec_atkin_struct * A, slong nA, const fmpz_t m3, const fmpz_t q);
WARN_UNUSED_RESULT int _gr_ec_cardinality_mod_p(fmpz_t res, const fmpz_t p, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_cardinality_match_sort(fmpz_t res, const fmpz_t t3, const fmpz_t m3, const gr_ec_atkin_struct * A, slong nA, gr_ec_ctx_t ctx);

#endif
