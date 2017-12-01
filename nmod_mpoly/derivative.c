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
                slong var, int deg, int rev, slong nfields, slong bits, slong N,
                                                        const nmodf_ctx_t fctx)
{
    slong off, shift, fpw, i, len1;
    ulong c, mask;
    mp_limb_t cr;
    ulong * one;
    TMP_INIT;

    TMP_START;

    fpw = FLINT_BITS/bits;
    mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    mpoly_off_shift(&off, &shift, var, deg, rev, fpw, nfields, bits);
    one = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_univar_exp(one, var, deg, N, off, shift, fpw, bits);

    /* x^c -> c*x^(c-1) */
    len1 = 0;
    for (i = 0; i < len2; i++)
    {
        c = (exp2[N*i + off] >> shift) & mask;
        if (c != 0)
        {
            mpoly_monomial_sub(exp1 + N*len1, exp2 + N*i, one, N);
            NMOD_RED(cr, c, fctx->mod);
            coeff1[len1] = nmod_mul(coeff2[i], cr, fctx->mod);
            len1 += (coeff1[len1] != 0);
        }
    }
    TMP_END;
    return len1;
}

void nmod_mpoly_derivative(nmod_mpoly_t poly1, const nmod_mpoly_t poly2,
                                         slong var, const nmod_mpoly_ctx_t ctx)
{
    int deg, rev;
    slong N, len1;

    N = words_per_exp(ctx->n, poly2->bits);
    degrev_from_ord(deg, rev, ctx->ord);

    nmod_mpoly_fit_length(poly1, poly2->length, ctx);
    nmod_mpoly_fit_bits(poly1, poly2->bits, ctx);
    poly1->bits = poly2->bits;

    len1 = _nmod_mpoly_derivative(poly1->coeffs, poly1->exps,
                   poly2->coeffs, poly2->exps, poly2->length,
                        var, deg, rev, ctx->n, poly2->bits, N, ctx->ffinfo);

    _nmod_mpoly_set_length(poly1, len1, ctx);
}
