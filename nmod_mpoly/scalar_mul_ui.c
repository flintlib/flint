/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

slong _nmod_mpoly_scalar_mul_nmod(mp_limb_t * coeff1, ulong * exp1,
                const mp_limb_t * coeff2, const ulong * exp2, slong len2,
                                      slong N, ulong c, const nmodf_ctx_t fctx)
{
    slong i;

    if (exp1 != exp2)
    {
        for (i = 0; i < len2; i++)
            mpoly_monomial_set(exp1 + N*i, exp2 + N*i, N);
    }

    for (i = 0; i < len2; i++)
        coeff1[i] = nmod_mul(coeff2[i], c, fctx->mod);

    while (i > 0 && coeff1[i - 1] == 0)
        i--;

    return i;  
}

void nmod_mpoly_scalar_mul_ui(nmod_mpoly_t poly1, const nmod_mpoly_t poly2,
                                          ulong c, const nmod_mpoly_ctx_t ctx)
{
    slong N, len1;
    mp_limb_t cr;

    N = mpoly_words_per_exp(poly2->bits, ctx->minfo);

    nmod_mpoly_fit_length(poly1, poly2->length, ctx);
    nmod_mpoly_fit_bits(poly1, poly2->bits, ctx);

    NMOD_RED(cr, c, ctx->ffinfo->mod);
    len1 = _nmod_mpoly_scalar_mul_nmod(poly1->coeffs, poly1->exps, 
                poly2->coeffs, poly2->exps, poly2->length, N, cr, ctx->ffinfo);
      
    _nmod_mpoly_set_length(poly1, len1, ctx);
    poly1->bits = poly2->bits;
}
