/*
    Copyright (C) 2011, 2012 Sebastian Pancratz 
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod.h"

void fq_nmod_randtest(fq_nmod_t rop, flint_rand_t state, const fq_nmod_ctx_t ctx)
{
    const slong d = fq_nmod_ctx_degree(ctx);
    slong i, sparse;

    nmod_poly_fit_length(rop, d);

    if (n_randint(state, 2))
    {
        for (i = 0; i < d; i++)
            rop->coeffs[i] = n_randint(state, ctx->mod.n);
    }
    else
    {
        sparse = 1 + n_randint(state, FLINT_MAX(2, d));

        for (i = 0; i < d; i++)
            if (n_randint(state, sparse))
                rop->coeffs[i] = WORD(0);
            else
                rop->coeffs[i] = n_randint(state, ctx->mod.n);
    }

    _nmod_poly_set_length(rop, d);
    _nmod_poly_normalise(rop);
}

void fq_nmod_randtest_dense(fq_nmod_t rop, flint_rand_t state, const fq_nmod_ctx_t ctx)
{
    const slong d = fq_nmod_ctx_degree(ctx);
    slong i;

    nmod_poly_fit_length(rop, d);

    for (i = 0; i < d - 1; i++)
        rop->coeffs[i] = n_randint(state, ctx->mod.n);

    rop->coeffs[d - 1] = 1;

    _nmod_poly_set_length(rop, d);
    _nmod_poly_normalise(rop);
}

void fq_nmod_randtest_not_zero(fq_nmod_t rop, flint_rand_t state, const fq_nmod_ctx_t ctx)
{
    slong i;

    fq_nmod_randtest(rop, state, ctx);
    for (i = 0; fq_nmod_is_zero(rop, ctx) && (i < 10); i++)
        fq_nmod_randtest(rop, state, ctx);

    if (fq_nmod_is_zero(rop, ctx))
        fq_nmod_one(rop, ctx);
}

