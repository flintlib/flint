/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"


void fmpz_mpoly_compression_do(
    fmpz_mpoly_t L,
    const fmpz_mpoly_ctx_t Lctx,
    fmpz * Acoeffs,
    slong Alen,
    mpoly_compression_t M)
{    
    slong i, max_deg;
    flint_bitcnt_t Lbits;
    slong mvars = Lctx->minfo->nvars;
    slong nvars = M->nvars;
    slong LN;

    FLINT_ASSERT(mvars == M->mvars);

    max_deg = M->degs[0];
    for (i = 1; i < mvars; i++)
        max_deg = FLINT_MAX(max_deg, M->degs[i]);
    Lbits = mpoly_fix_bits(1 + FLINT_BIT_COUNT(max_deg), Lctx->minfo);

    fmpz_mpoly_fit_length_reset_bits(L, Alen, Lbits, Lctx);

    LN = mpoly_words_per_exp_sp(Lbits, Lctx->minfo);

    L->length = Alen;
    for (i = 0; i < Alen; i++)
    {
        fmpz_swap(L->coeffs + i, Acoeffs + i);
        mpoly_set_monomial_ui(L->exps + LN*i, (ulong *)M->exps + nvars*i,
                                                           Lbits, Lctx->minfo);
    }

    fmpz_mpoly_sort_terms(L, Lctx);
    fmpz_mpoly_unit_normalize(L, Lctx);
}


void fmpz_mpoly_compression_undo(
    fmpz_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz_mpoly_ctx_t Actx,
    fmpz_mpoly_t L,
    const fmpz_mpoly_ctx_t Lctx,
    mpoly_compression_t M)
{
    slong i, k, l;
    slong nvars = Actx->minfo->nvars;
    slong NA = mpoly_words_per_exp(Abits, Actx->minfo);
    slong mvars = Lctx->minfo->nvars;
    flint_bitcnt_t Lbits = L->bits;
    slong NL = mpoly_words_per_exp(Lbits, Lctx->minfo);
    slong * mins, * texps;
    TMP_INIT;

    FLINT_ASSERT(mvars == M->mvars);
    FLINT_ASSERT(fmpz_mpoly_degrees_fit_si(L, Lctx));

    TMP_START;
    texps = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    mins = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    for (k = 0; k < nvars; k++)
        mins[k] = WORD_MAX;

    _slong_array_fit_length(&M->exps, &M->exps_alloc, L->length*nvars);
    fmpz_mpoly_fit_length_reset_bits(A, L->length, Abits, Actx);
    _fmpz_mpoly_set_length(A, L->length, Actx);
    for (i = 0; i < L->length; i++)
    {
        fmpz_swap(A->coeffs + i, L->coeffs + i);
        mpoly_get_monomial_ui((ulong *)texps, L->exps + NL*i, Lbits, Lctx->minfo);
        for (k = 0; k < nvars; k++)
        {
            slong tot = M->deltas[k];
            for (l = 0; l < mvars; l++)
                tot += M->umat[k*nvars + l]*texps[l];
            M->exps[i*nvars + k] = tot;
            mins[k] = FLINT_MIN(mins[k], tot);
        }
    }

    for (i = 0; i < L->length; i++)
    {
        for (k = 0; k < nvars; k++)
            M->exps[i*nvars + k] -= mins[k];
        mpoly_set_monomial_ui(A->exps + NA*i, (ulong *)M->exps + i*nvars,
                                                           Abits, Actx->minfo);
    }

    TMP_END;

    fmpz_mpoly_sort_terms(A, Actx);
    fmpz_mpoly_unit_normalize(A, Actx);
}

