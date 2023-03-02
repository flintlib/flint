/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

/* evaluate B(xbar) at xbar = C */
int fq_nmod_mpoly_compose_fq_nmod_mpoly_geobucket(fq_nmod_mpoly_t A,
               const fq_nmod_mpoly_t B, fq_nmod_mpoly_struct * const * C,
               const fq_nmod_mpoly_ctx_t ctxB, const fq_nmod_mpoly_ctx_t ctxAC)
{
    slong d = fq_nmod_ctx_degree(ctxAC->fqctx);
    int success = 1;
    slong i, j;
    slong Blen = B->length;
    const ulong * Bexp = B->exps;
    flint_bitcnt_t Bbits = B->bits;
    slong BN = mpoly_words_per_exp(Bbits, ctxB->minfo);
    fq_nmod_mpoly_t U, V, W;
    fq_nmod_mpoly_geobucket_t T;
    fmpz * e;

    fq_nmod_mpoly_init(U, ctxAC);
    fq_nmod_mpoly_init(V, ctxAC);
    fq_nmod_mpoly_init(W, ctxAC);
    fq_nmod_mpoly_geobucket_init(T, ctxAC);
    e = _fmpz_vec_init(ctxB->minfo->nvars);

    for (i = 0; success && i < Blen; i++)
    {
        fq_nmod_mpoly_set_n_fq(U, B->coeffs + d*i, ctxAC);
        mpoly_get_monomial_ffmpz(e, Bexp + BN*i, Bbits, ctxB->minfo);
        for (j = 0; j < ctxB->minfo->nvars; j++)
        {
            success = success && fq_nmod_mpoly_pow_fmpz(V, C[j], e + j, ctxAC);
            fq_nmod_mpoly_mul(W, U, V, ctxAC);
            fq_nmod_mpoly_swap(U, W, ctxAC);
        }
        fq_nmod_mpoly_geobucket_add(T, U, ctxAC);
    }

    if (success)
        fq_nmod_mpoly_geobucket_empty(A, T, ctxAC);

    fq_nmod_mpoly_clear(U, ctxAC);
    fq_nmod_mpoly_clear(V, ctxAC);
    fq_nmod_mpoly_clear(W, ctxAC);
    fq_nmod_mpoly_geobucket_clear(T, ctxAC);
    _fmpz_vec_clear(e, ctxB->minfo->nvars);

    return success;
}

