/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void nmod_mpoly_add_ui(nmod_mpoly_t A, const nmod_mpoly_t B,
                                           ulong c, const nmod_mpoly_ctx_t ctx)
{
    slong i, N;
    slong Blen = B->length;
    const nmodf_ctx_struct * fctx = ctx->ffinfo;

    if (Blen == 0)
    {
        nmod_mpoly_set_ui(A, c, ctx);
        return;
    }

    if (c >= ctx->ffinfo->mod.n)
    {
        NMOD_RED(c, c, ctx->ffinfo->mod);
    }
    if (c != 0)
    {
        N = mpoly_words_per_exp(B->bits, ctx->minfo);

        if (mpoly_monomial_is_zero(B->exps + (Blen - 1)*N, N))
        {
            if (A != B)
            {
                nmod_mpoly_fit_length(A, B->length, ctx);
                nmod_mpoly_fit_bits(A, B->bits, ctx);
                A->bits = B->bits;

                for (i = 0; i < Blen - 1; i++)
                    A->coeffs[i] = B->coeffs[i];

                for (i = 0; i < Blen; i++)
                    mpoly_monomial_set(A->exps + i*N, B->exps + i*N, N);

                _nmod_mpoly_set_length(A, B->length, ctx);
            }

            A->coeffs[Blen - 1] = nmod_add(B->coeffs[Blen - 1], c, fctx->mod);
            if (A->coeffs[Blen - 1] == 0)
                _nmod_mpoly_set_length(A, Blen - 1, ctx);
        }
        else
        {
            nmod_mpoly_fit_length(A, Blen + 1, ctx);

            if (A != B)
            {
                nmod_mpoly_fit_bits(A, B->bits, ctx);
                A->bits = B->bits;

                for (i = 0; i < Blen; i++)
                    A->coeffs[i] = B->coeffs[i];

                for (i = 0; i < Blen; i++)
                    mpoly_monomial_set(A->exps + i*N, B->exps + i*N, N);
            }

            mpoly_monomial_zero(A->exps + Blen*N, N);

            A->coeffs[Blen] = c;
            _nmod_mpoly_set_length(A, Blen + 1, ctx);
        }
    }
    else if (A != B)
    {
        nmod_mpoly_set(A, B, ctx);
    }
}
