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
nmod_poly_factor_distinct_deg(nmod_poly_factor_t res, const nmod_poly_t poly,
    slong * const * degs)
{
    gr_ctx_t ctx;
    gr_poly_t P;
    gr_poly_vec_t fac;
    fmpz_vec_t degrees;
    nmod_t mod = poly->mod;
    slong i, num;

    _gr_ctx_init_nmod(ctx, &mod);
    GR_MUST_SUCCEED(gr_ctx_set_is_field(ctx, T_TRUE));

    NMOD_POLY_AS_GR(P, poly);
    gr_poly_vec_init(fac, 0, ctx);
    fmpz_vec_init(degrees, 0);

    GR_MUST_SUCCEED(gr_poly_factor_distinct_deg(fac, degrees, P, ctx));

    num = fac->length;
    nmod_poly_factor_fit_length(res, res->num + num);

    for (i = 0; i < num; i++)
    {
        nmod_poly_struct * p = res->p + res->num + i;

        nmod_poly_clear(p);
        NMOD_POLY_FROM_GR(p, fac->entries + i);
        p->mod = mod;
        gr_poly_init(fac->entries + i, ctx);

        res->exp[res->num + i] = 1;
        (*degs)[i] = fmpz_get_si(degrees->entries + i);
    }

    res->num += num;

    gr_poly_vec_clear(fac, ctx);
    fmpz_vec_clear(degrees);
    gr_ctx_clear(ctx);
}
