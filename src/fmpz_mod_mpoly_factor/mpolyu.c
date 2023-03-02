/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"


/** 0 variables *************************************/


void fmpz_mod_mpoly_to_mpolyl_perm_deflate(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_ctx_t lctx,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride)
{
    slong j, k, l;
    slong m = lctx->minfo->nvars;
    slong n = ctx->minfo->nvars;
    slong NA, NB;
    ulong * lexps;
    ulong * Bexps;
    TMP_INIT;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(m <= n);
    TMP_START;

    fmpz_mod_mpoly_fit_length(A, B->length, ctx);
    A->length = B->length;

    lexps = (ulong *) TMP_ALLOC(m*sizeof(ulong));
    Bexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));
    NA = mpoly_words_per_exp(A->bits, lctx->minfo);
    NB = mpoly_words_per_exp(B->bits, ctx->minfo);

    for (j = 0; j < B->length; j++)
    {
        fmpz_set(A->coeffs + j, B->coeffs + j);

        mpoly_get_monomial_ui(Bexps, B->exps + NB*j, B->bits, ctx->minfo);
        for (k = 0; k < m; k++)
        {
            l = perm[k];
            if (stride[l] == 1)
            {
                lexps[k] = (Bexps[l] - shift[l]);
            }
            else
            {
                FLINT_ASSERT(stride[l] != 0);
                FLINT_ASSERT(((Bexps[l] - shift[l]) % stride[l]) == 0);
                lexps[k] = (Bexps[l] - shift[l]) / stride[l];
            }
        }

        mpoly_set_monomial_ui(A->exps + NA*j, lexps, A->bits, lctx->minfo);
    }

    TMP_END;

    fmpz_mod_mpoly_sort_terms(A, lctx);

    FLINT_ASSERT(fmpz_mod_mpoly_is_canonical(A, lctx));
}


/*
    Convert B to A using the variable permutation vector perm.
    A must be constructed with bits = Abits.

    operation on each term:

        for 0 <= l < n
            Aexp[l] = shift[l]

        for 0 <= k < m + 2
            l = perm[k]
            Aexp[l] += scale[l]*Bexp[k]
*/
void fmpz_mod_mpoly_from_mpolyl_perm_inflate(
    fmpz_mod_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz_mod_mpoly_ctx_t ctx,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t lctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride)
{
    slong n = ctx->minfo->nvars;
    slong m = lctx->minfo->nvars;
    slong i, k, l;
    slong NA, NB;
    ulong * Bexps;
    ulong * Aexps;
    TMP_INIT;

    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(Abits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(m <= n);
    TMP_START;

    Bexps = (ulong *) TMP_ALLOC(m*sizeof(ulong));
    Aexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    NA = mpoly_words_per_exp(Abits, ctx->minfo);
    NB = mpoly_words_per_exp(B->bits, lctx->minfo);

    fmpz_mod_mpoly_fit_length_reset_bits(A, B->length, Abits, ctx);
    A->length = B->length;

    for (i = 0; i < B->length; i++)
    {
	    fmpz_set(A->coeffs + i, B->coeffs + i);
	    mpoly_get_monomial_ui(Bexps, B->exps + NB*i, B->bits, lctx->minfo);

        for (l = 0; l < n; l++)
        {
            Aexps[l] = shift[l];
        }
        for (k = 0; k < m; k++)
        {
            l = perm[k];
            Aexps[l] += stride[l]*Bexps[k];
        }

        mpoly_set_monomial_ui(A->exps + NA*i, Aexps, Abits, ctx->minfo);
    }

    TMP_END;

    fmpz_mod_mpoly_sort_terms(A, ctx);

    FLINT_ASSERT(fmpz_mod_mpoly_is_canonical(A, ctx));
}
