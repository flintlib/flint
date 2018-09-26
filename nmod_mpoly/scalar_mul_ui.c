/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

slong _nmod_mpoly_scalar_mul_nmod(mp_limb_t * Acoeff, ulong * Aexp,
                const mp_limb_t * Bcoeff, const ulong * Bexp, slong Blen,
                                      slong N, ulong c, const nmodf_ctx_t fctx)
{
    slong i, Alen;

    Alen = 0;
    for (i = 0; i < Blen; i++)
    {
        mpoly_monomial_set(Aexp + N*Alen, Bexp + N*i, N);
        Acoeff[Alen] = nmod_mul(Bcoeff[i], c, fctx->mod);
        Alen += (Acoeff[Alen] != 0);
    }

    return Alen;
}

void nmod_mpoly_scalar_mul_ui(nmod_mpoly_t A, const nmod_mpoly_t B,
                                           ulong c, const nmod_mpoly_ctx_t ctx)
{
    slong N, len1;
    mp_limb_t cr;

    N = mpoly_words_per_exp(B->bits, ctx->minfo);

    nmod_mpoly_fit_length(A, B->length, ctx);
    nmod_mpoly_fit_bits(A, B->bits, ctx);

    NMOD_RED(cr, c, ctx->ffinfo->mod);
    len1 = _nmod_mpoly_scalar_mul_nmod(A->coeffs, A->exps, 
                            B->coeffs, B->exps, B->length, N, cr, ctx->ffinfo);
      
    _nmod_mpoly_set_length(A, len1, ctx);
    A->bits = B->bits;
}
