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

void
nmod_poly_factor_equal_deg(nmod_poly_factor_t factors, const nmod_poly_t pol,
    slong d)
{
    gr_ctx_t ctx;
    gr_poly_t P;
    gr_poly_vec_t fac;
    fmpz_vec_t exp;
    nmod_t mod = pol->mod;
    slong i;

    _gr_ctx_init_nmod(ctx, &mod);
    GR_MUST_SUCCEED(gr_ctx_set_is_field(ctx, T_TRUE));

    NMOD_POLY_AS_GR(P, pol);
    gr_poly_vec_init(fac, 0, ctx);
    fmpz_vec_init(exp, 0);

    GR_MUST_SUCCEED(gr_poly_factor_equal_deg(fac, P, d, ctx));

    for (i = 0; i < fac->length; i++)
        fmpz_vec_append_ui(exp, 1);

    _nmod_poly_factor_set_gr(factors, fac, exp, mod, ctx);

    gr_poly_vec_clear(fac, ctx);
    fmpz_vec_clear(exp);
    gr_ctx_clear(ctx);
}
