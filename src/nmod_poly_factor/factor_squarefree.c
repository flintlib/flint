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
nmod_poly_factor_squarefree(nmod_poly_factor_t res, const nmod_poly_t f)
{
    gr_ctx_t ctx;
    gr_poly_t P;
    gr_poly_vec_t fac;
    fmpz_vec_t exp;
    nmod_t mod = f->mod;
    ulong lc;

    if (f->length <= 1)
    {
        res->num = 0;
        return;
    }

    _gr_ctx_init_nmod(ctx, &mod);
    GR_MUST_SUCCEED(gr_ctx_set_is_field(ctx, T_TRUE));

    NMOD_POLY_AS_GR(P, f);
    gr_poly_vec_init(fac, 0, ctx);
    fmpz_vec_init(exp, 0);

    GR_MUST_SUCCEED(gr_poly_factor_squarefree(&lc, fac, exp, P, ctx));

    res->num = 0;
    _nmod_poly_factor_set_gr(res, fac, exp, mod, ctx);

    gr_poly_vec_clear(fac, ctx);
    fmpz_vec_clear(exp);
    gr_ctx_clear(ctx);
}
