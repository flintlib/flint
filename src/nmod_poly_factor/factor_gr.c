/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "nmod_poly_factor.h"
#include "nmod_poly_factor_gr.h"

/*
    Wrappers around the gr_poly implementations. The factors are moved
    (not copied) out of the gr_poly_vec into the nmod_poly_factor.
*/
void
_nmod_poly_factor_gr(nmod_poly_factor_t res, const nmod_poly_t f, int algorithm)
{
    gr_ctx_t ctx;
    gr_poly_t P;
    gr_poly_vec_t fac;
    fmpz_vec_t exp;
    nmod_t mod = f->mod;
    ulong lc;

    _gr_ctx_init_nmod(ctx, &mod);
    GR_MUST_SUCCEED(gr_ctx_set_is_field(ctx, T_TRUE));

    NMOD_POLY_AS_GR(P, f);
    gr_poly_vec_init(fac, 0, ctx);
    fmpz_vec_init(exp, 0);

    GR_MUST_SUCCEED(_gr_poly_factor_finite_field(&lc, fac, exp, P, algorithm, ctx));

    _nmod_poly_factor_set_gr(res, fac, exp, mod, ctx);

    gr_poly_vec_clear(fac, ctx);
    fmpz_vec_clear(exp);
    gr_ctx_clear(ctx);
}

void
nmod_poly_factor_cantor_zassenhaus(nmod_poly_factor_t res, const nmod_poly_t f)
{
    _nmod_poly_factor_gr(res, f, GR_POLY_FACTOR_ALGORITHM_CANTOR_ZASSENHAUS);
}

void
nmod_poly_factor_berlekamp(nmod_poly_factor_t res, const nmod_poly_t f)
{
    _nmod_poly_factor_gr(res, f, GR_POLY_FACTOR_ALGORITHM_BERLEKAMP);
}

void
nmod_poly_factor_kaltofen_shoup(nmod_poly_factor_t res, const nmod_poly_t f)
{
    _nmod_poly_factor_gr(res, f, GR_POLY_FACTOR_ALGORITHM_KALTOFEN_SHOUP);
}
