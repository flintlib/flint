/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

void fmpz_mod_mpoly_deflate(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz * shift,
    const fmpz * stride,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong bits = B->bits;
    slong NA = mpoly_words_per_exp(bits, ctx->minfo);

    if (fmpz_mod_mpoly_is_zero(B, ctx))
    {
        fmpz_mod_mpoly_zero(A, ctx);
        return;
    }

    if (A == B)
    {
        slong new_alloc = NA*A->length;
        ulong * texps = flint_malloc(new_alloc*sizeof(ulong));
        mpoly_monomials_deflate(texps, bits, B->exps, B->bits, B->length,
                                                    shift, stride, ctx->minfo);
        flint_free(A->exps);
        A->exps = texps;
        A->bits = bits;
        A->exps_alloc = new_alloc;
    }
    else
    {
        fmpz_mod_mpoly_fit_length_reset_bits(A, B->length, bits, ctx);
        _fmpz_vec_set(A->coeffs, B->coeffs, B->length);
        mpoly_monomials_deflate(A->exps, bits, B->exps, B->bits, B->length,
                                                    shift, stride, ctx->minfo);
        _fmpz_mod_mpoly_set_length(A, B->length, ctx);
    }

    if (ctx->minfo->ord != ORD_LEX)
        fmpz_mod_mpoly_sort_terms(A, ctx);
}
