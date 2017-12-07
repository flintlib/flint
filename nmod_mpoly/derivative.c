/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

slong _nmod_mpoly_derivative(mp_limb_t * coeff1,       ulong * exp1,
                       const mp_limb_t * coeff2, const ulong * exp2, slong len2,
     slong var, slong bits, slong N, slong offset, slong shift, ulong * oneexp,
                                                        const nmodf_ctx_t fctx)
{
    slong i, len1;

    /* x^c -> c*x^(c-1) */
    len1 = 0;
    for (i = 0; i < len2; i++)
    {
        ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
        ulong c = (exp2[N*i + offset] >> shift) & mask;
        mp_limb_t cr;
        if (c != 0)
        {
            mpoly_monomial_sub(exp1 + N*len1, exp2 + N*i, oneexp, N);
            NMOD_RED(cr, c, fctx->mod);
            coeff1[len1] = nmod_mul(coeff2[i], cr, fctx->mod);
            len1 += (coeff1[len1] != 0);
        }
    }
    return len1;
}

void nmod_mpoly_derivative(nmod_mpoly_t poly1, const nmod_mpoly_t poly2,
                                         slong var, const nmod_mpoly_ctx_t ctx)
{
    slong bits, N, offset, shift;
    ulong * oneexp;
    slong len1;
    TMP_INIT;

    TMP_START;
    bits = poly2->bits;
    nmod_mpoly_fit_length(poly1, poly2->length, ctx);
    nmod_mpoly_fit_bits(poly1, bits, ctx);
    poly1->bits = bits;

    N = mpoly_words_per_exp(bits, ctx->minfo);
    oneexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_gen_oneexp_offset_shift(oneexp, &offset, &shift,
                                                     var, N, bits, ctx->minfo);

    len1 = _nmod_mpoly_derivative(poly1->coeffs, poly1->exps,
                                  poly2->coeffs, poly2->exps, poly2->length,
                             var, bits, N, offset, shift, oneexp, ctx->ffinfo);

    _nmod_mpoly_set_length(poly1, len1, ctx);
    TMP_END;
}
