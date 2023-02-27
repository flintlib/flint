/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void nmod_mpoly_add_ui(nmod_mpoly_t A, const nmod_mpoly_t B,
                                           ulong c, const nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(B->bits, ctx->minfo);
    slong Blen = B->length;

    if (c >= ctx->mod.n)
        NMOD_RED(c, c, ctx->mod);

    if (c == 0)
    {
        nmod_mpoly_set(A, B, ctx);
        return;
    }

    if (Blen < 1)
    {
        nmod_mpoly_set_ui(A, c, ctx);
        return;
    }

    if (mpoly_monomial_is_zero(B->exps + (Blen - 1)*N, N))
    {
        if (A != B)
        {
            nmod_mpoly_fit_length_reset_bits(A, B->length, B->bits, ctx);
            _nmod_vec_set(A->coeffs, B->coeffs, Blen - 1);
            mpoly_copy_monomials(A->exps, B->exps, Blen, N);
            _nmod_mpoly_set_length(A, B->length, ctx);
        }

        A->coeffs[Blen - 1] = nmod_add(B->coeffs[Blen - 1], c, ctx->mod);
        if (A->coeffs[Blen - 1] == 0)
            _nmod_mpoly_set_length(A, Blen - 1, ctx);
    }
    else
    {
        if (A != B)
        {
            nmod_mpoly_fit_length_reset_bits(A, Blen + 1, B->bits, ctx);
            _nmod_vec_set(A->coeffs, B->coeffs, Blen);
            mpoly_copy_monomials(A->exps, B->exps, Blen, N);
        }
        else
        {
            nmod_mpoly_fit_length(A, Blen + 1, ctx);
        }

        mpoly_monomial_zero(A->exps + N*Blen, N);
        A->coeffs[Blen] = c;
        _nmod_mpoly_set_length(A, Blen + 1, ctx);
    }
}
