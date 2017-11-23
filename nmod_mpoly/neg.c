/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void _nmod_mpoly_neg(ulong * poly1, ulong * exps1,
                     const ulong * poly2, const ulong * exps2, slong len,
                                               slong N, const nmodf_ctx_t fctx)
{
    slong i, D;

    D = fctx->deg;

    for (i = 0; i < len; i++)
        nmodf_neg(poly1 + i*D, poly2 + i*D, fctx);

    if (exps1 != exps2)
    {
        for (i = 0; i < len; i++)
            mpoly_monomial_set(exps1 + i*N, exps2 + i*N, N);
    }
}

void nmod_mpoly_neg(nmod_mpoly_t poly1, const nmod_mpoly_t poly2,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong N;
    N = words_per_exp(ctx->n, poly2->bits);

    nmod_mpoly_fit_length(poly1, poly2->length, ctx);
    nmod_mpoly_fit_bits(poly1, poly2->bits, ctx);

    _nmod_mpoly_neg(poly1->coeffs, poly1->exps,
                   poly2->coeffs, poly2->exps, poly2->length, N, ctx->ffinfo);

    _nmod_mpoly_set_length(poly1, poly2->length, ctx);
    poly1->bits = poly2->bits;
}
