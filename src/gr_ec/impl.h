/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef GR_EC_IMPL_H
#define GR_EC_IMPL_H

#include "gr_ec.h"

/*
    Variants of the Jacobian doubling and mixed addition taking a caller
    supplied workspace of GR_EC_JAC_SCRATCH base ring elements. Allocating
    and initializing the temporaries costs a noticeable fraction of an
    operation over small moduli, so the scalar multiplication ladder hoists
    the workspace out of its loop and calls these directly.
*/

#define GR_EC_JAC_SCRATCH 9

/*
    The group law branches on whether certain quantities vanish, and over a
    ring that is not an integral domain a nonzero value is not enough to
    make the branch right: it has to be a unit, or the formula is correct
    modulo one factor of the modulus and wrong modulo another. The
    functions below therefore take an optional witness w, which they
    multiply by every quantity a branch has just treated as nonzero. Where
    the caller can test it -- gcd(w, n) over Z/n -- a non-unit witness says
    the computation is not to be trusted, and over Z/n it is a factor.

    Passing NULL asks for none of this and costs nothing.
*/
#define GR_EC_WITNESS(w, x) \
    do { if ((w) != NULL) status |= gr_mul((w), (w), (x), R); } while (0)

WARN_UNUSED_RESULT int _gr_ec_jac_point_dbl_long_weierstrass_ws(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, gr_ptr t, gr_ptr w, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_jac_point_dbl_short_weierstrass_ws(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, gr_ptr t, gr_ptr w, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_jac_point_add_aff_point_long_weierstrass_ws(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const gr_ec_aff_point_t Q, gr_ptr t, gr_ptr w, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_jac_point_add_aff_point_short_weierstrass_ws(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const gr_ec_aff_point_t Q, gr_ptr t, gr_ptr w, gr_ec_ctx_t ctx);

/* the ladder and the batch normalisation, both taking a witness */
WARN_UNUSED_RESULT int _gr_ec_jac_point_mul_fmpz_naf_witness(gr_ec_jac_point_t res, const gr_ec_jac_point_t P, const fmpz_t n, gr_ptr w, gr_ec_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_ec_jac_point_vec_get_aff_point_vec_witness(gr_ec_aff_point_struct * res, const gr_ec_jac_point_struct * P, slong len, gr_ptr w, gr_ec_ctx_t ctx);

/* Points the context at the method table and element size of a
   representation. Defined in generic.c; called when a context is created. */
void _gr_ec_ctx_init_methods(gr_ec_ctx_t ctx, gr_ec_repr_t repr);

/*
    Scalars taken from another ring. Defined in mul_other.c and wired into
    the method tables as MUL_OTHER and OTHER_MUL; see that file for which
    rings are allowed to act and why.
*/

WARN_UNUSED_RESULT int _gr_ec_scalar_of_other(fmpz_t k, gr_srcptr y, gr_ctx_t y_ctx, gr_ec_ctx_t ctx);

#define GR_EC_MUL_OTHER_DECL(kind) \
WARN_UNUSED_RESULT int _ ## kind ## _mul_other(kind ## _t res, const kind ## _t P, gr_srcptr y, gr_ctx_t y_ctx, gr_ec_ctx_t ctx); \
WARN_UNUSED_RESULT int _ ## kind ## _other_mul(kind ## _t res, gr_srcptr y, gr_ctx_t y_ctx, const kind ## _t P, gr_ec_ctx_t ctx);

GR_EC_MUL_OTHER_DECL(gr_ec_point)
GR_EC_MUL_OTHER_DECL(gr_ec_aff_point)
GR_EC_MUL_OTHER_DECL(gr_ec_jac_point)

#endif
