/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

slong _nmod_mpoly_scalar_mul_nmod(ulong * poly1, ulong * exps1,
        const ulong * poly2, const ulong * exps2, slong len2, slong N, ulong c,
                                                        const nmodf_ctx_t fctx)
{
    slong i;

    if (exps1 != exps2)
        mpn_copyi(exps1, exps2, N*len2);

    for (i = 0; i < len2; i++)
        nmodf_smul_nmod(poly1 + i*fctx->deg, poly2 + i*fctx->deg, c, fctx);

    while (i > 0 && nmodf_is_zero(poly1 + (i-1)*fctx->deg, fctx))
        i--;

    return i;  
}

void nmod_mpoly_scalar_mul_ui(nmod_mpoly_t poly1, const nmod_mpoly_t poly2,
                                          ulong c, const nmod_mpoly_ctx_t ctx)
{
    slong N, len1;
    ulong cr;

    N = words_per_exp(ctx->n, poly2->bits);

    nmod_mpoly_fit_length(poly1, poly2->length, ctx);
    nmod_mpoly_fit_bits(poly1, poly2->bits, ctx);

    NMOD_RED(cr, c, ctx->ffinfo->mod);
    len1 = _nmod_mpoly_scalar_mul_nmod(poly1->coeffs, poly1->exps, 
                poly2->coeffs, poly2->exps, poly2->length, N, cr, ctx->ffinfo);
      
    _nmod_mpoly_set_length(poly1, len1, ctx);
    poly1->bits = poly2->bits;
}
