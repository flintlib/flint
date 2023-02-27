/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

/* evaluate B(xbar) at xbar = C */
int fmpz_mod_mpoly_compose_fmpz_mod_mpoly_geobucket(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    fmpz_mod_mpoly_struct * const * C,
    const fmpz_mod_mpoly_ctx_t ctxB,
    const fmpz_mod_mpoly_ctx_t ctxAC)
{
    int success = 1;
    slong i, j;
    slong Blen = B->length;
    const fmpz * Bcoeff = B->coeffs;
    const ulong * Bexp = B->exps;
    flint_bitcnt_t Bbits = B->bits;
    slong BN = mpoly_words_per_exp(Bbits, ctxB->minfo);
    fmpz_mod_mpoly_t U, V, W;
    fmpz_mod_mpoly_geobucket_t T;
    fmpz * e;

    fmpz_mod_mpoly_init(U, ctxAC);
    fmpz_mod_mpoly_init(V, ctxAC);
    fmpz_mod_mpoly_init(W, ctxAC);
    fmpz_mod_mpoly_geobucket_init(T, ctxAC);
    e = _fmpz_vec_init(ctxB->minfo->nvars);

    for (i = 0; success && i < Blen; i++)
    {
        fmpz_mod_mpoly_set_fmpz(U, Bcoeff + i, ctxAC);
        mpoly_get_monomial_ffmpz(e, Bexp + BN*i, Bbits, ctxB->minfo);
        for (j = 0; j < ctxB->minfo->nvars; j++)
        {
            success = success && fmpz_mod_mpoly_pow_fmpz(V, C[j], e + j, ctxAC);
            fmpz_mod_mpoly_mul(W, U, V, ctxAC);
            fmpz_mod_mpoly_swap(U, W, ctxAC);
        }
        fmpz_mod_mpoly_geobucket_add(T, U, ctxAC);
    }

    if (success)
        fmpz_mod_mpoly_geobucket_empty(A, T, ctxAC);

    fmpz_mod_mpoly_clear(U, ctxAC);
    fmpz_mod_mpoly_clear(V, ctxAC);
    fmpz_mod_mpoly_clear(W, ctxAC);
    fmpz_mod_mpoly_geobucket_clear(T, ctxAC);
    _fmpz_vec_clear(e, ctxB->minfo->nvars);

    return success;
}
