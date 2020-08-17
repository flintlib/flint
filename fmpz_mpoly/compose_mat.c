/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

/* essentially exps(A) = M*exps(B) */
void _fmpz_mpoly_compose_mat(fmpz_mpoly_t A,
                            const fmpz_mpoly_t B, const fmpz_mat_t M,
                     const fmpz_mpoly_ctx_t ctxB, const fmpz_mpoly_ctx_t ctxAC)
{
    slong i;
    fmpz * u, * v;
    flint_bitcnt_t vbits;
    slong Blen = B->length;
    flint_bitcnt_t Bbits = B->bits;
    slong BN = mpoly_words_per_exp(Bbits, ctxB->minfo);
    const ulong * Bexp = B->exps;
    const fmpz * Bcoeffs = B->coeffs;
    slong Alen_old = A->length;
    slong AN;

    FLINT_ASSERT(A != B);

    FLINT_ASSERT(fmpz_mat_nrows(M) == ctxAC->minfo->nfields + 1);
    FLINT_ASSERT(fmpz_mat_ncols(M) == ctxB->minfo->nfields);

    u = _fmpz_vec_init(ctxB->minfo->nfields);
    v = _fmpz_vec_init(ctxAC->minfo->nfields + 1);

    fmpz_mpoly_fit_length(A, Blen, ctxAC);
    A->length = 0;
    fmpz_mpoly_fit_bits(A, MPOLY_MIN_BITS, ctxAC);
    A->bits = MPOLY_MIN_BITS;
    for (i = 0; i < Blen; i++)
    {
        mpoly_unpack_vec_fmpz(u, Bexp + BN*i, Bbits, ctxB->minfo->nfields, 1);
        fmpz_mat_mul_vec(v, M, u);
        if (!fmpz_is_zero(v + ctxAC->minfo->nfields))
            continue;
        vbits = _fmpz_vec_max_bits(v, ctxAC->minfo->nfields);
        FLINT_ASSERT(vbits >= 0);
        fmpz_mpoly_fit_bits(A, mpoly_fix_bits(vbits + 1, ctxAC->minfo), ctxAC);
        fmpz_set(A->coeffs + A->length, Bcoeffs + i);
        AN = mpoly_words_per_exp(A->bits, ctxAC->minfo);
        mpoly_pack_vec_fmpz(A->exps + AN*A->length, v, A->bits, ctxAC->minfo->nfields, 1);
        A->length++;
    }

    while (--Alen_old >= A->length)
        _fmpz_demote(A->coeffs + Alen_old);

    _fmpz_vec_clear(u, ctxB->minfo->nfields);
    _fmpz_vec_clear(v, ctxAC->minfo->nfields + 1);

    fmpz_mpoly_sort_terms(A, ctxAC);
    fmpz_mpoly_combine_like_terms(A, ctxAC);
    return;
}

