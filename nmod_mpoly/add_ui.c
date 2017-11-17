/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void nmod_mpoly_add_ui(nmod_mpoly_t poly1, const nmod_mpoly_t poly2,
                                           ulong c, const nmod_mpoly_ctx_t ctx)
{
    ulong cr;
    slong i, N;
    slong len2 = poly2->length;
    const nmodf_ctx_struct * fctx = ctx->ffinfo;

    if (len2 == 0)
    {
        nmod_mpoly_set_ui(poly1, c, ctx);
        return;
    }

    NMOD_RED(cr, c, ctx->ffinfo->mod);
    if (cr != 0)
    {
        N = words_per_exp(ctx->n, poly2->bits);

        if (mpoly_monomial_is_zero(poly2->exps + (len2 - 1)*N, N))
        {
            if (poly1 != poly2)
            {
                nmod_mpoly_fit_length(poly1, poly2->length, ctx);
                nmod_mpoly_fit_bits(poly1, poly2->bits, ctx);

                for (i = 0; i < len2 - 1; i++)
                    nmodf_set(poly1->coeffs + i*fctx->deg, poly2->coeffs + i*fctx->deg, fctx);

                for (i = 0; i < len2; i++)
                    mpoly_monomial_set(poly1->exps + i*N, poly2->exps + i*N, N);

                _nmod_mpoly_set_length(poly1, poly2->length, ctx);
                poly1->bits = poly2->bits;
            }

            nmodf_sadd_nmod(poly1->coeffs + (len2 - 1)*fctx->deg, poly2->coeffs + (len2 - 1)*fctx->deg, cr, fctx);

            if (nmodf_is_zero(poly1->coeffs + (len2 - 1)*fctx->deg, fctx))
                _nmod_mpoly_set_length(poly1, len2 - 1, ctx);
        } else
        {
            nmod_mpoly_fit_length(poly1, len2 + 1, ctx);

            if (poly1 != poly2)
            {
                nmod_mpoly_fit_bits(poly1, poly2->bits, ctx);
                poly1->bits = poly2->bits;

                for (i = 0; i < len2; i++)
                    nmodf_set(poly1->coeffs + i*fctx->deg, poly2->coeffs + i*fctx->deg, fctx);

                for (i = 0; i < len2; i++)
                    mpoly_monomial_set(poly1->exps + i*N, poly2->exps + i*N, N);
            } 

            mpoly_monomial_zero(poly1->exps + len2*N, N);

            nmodf_sset_nmod(poly1->coeffs + len2*fctx->deg, cr, fctx);
            _nmod_mpoly_set_length(poly1, poly2->length + 1, ctx);
        }
    } else if (poly1 != poly2)
    {
        nmod_mpoly_set(poly1, poly2, ctx);
    }
}
