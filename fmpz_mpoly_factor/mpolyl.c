/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"

void fmpz_mpoly_to_mpolyl_perm_deflate(
    fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t lctx,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride)
{
    fmpz_mpoly_fit_length(A, B->length, ctx);
    mpoly_to_mpolyl_perm_deflate(A->exps, A->bits, lctx->minfo,
                                 B->exps, B->bits, ctx->minfo,
                                 B->length, perm, shift, stride);
    _fmpz_vec_set(A->coeffs, B->coeffs, B->length);
    _fmpz_mpoly_set_length(A, B->length, lctx);
    fmpz_mpoly_sort_terms(A, lctx);
    FLINT_ASSERT(fmpz_mpoly_is_canonical(A, lctx));
}

void fmpz_mpoly_from_mpolyl_perm_inflate(
    fmpz_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t lctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride)
{
    fmpz_mpoly_fit_length_reset_bits(A, B->length, Abits, ctx);
    _fmpz_vec_set(A->coeffs, B->coeffs, B->length);
    mpoly_from_mpolyl_perm_inflate(A->exps, Abits, ctx->minfo,
                                   B->exps, B->bits, lctx->minfo,
                                   B->length, perm, shift, stride);
    _fmpz_mpoly_set_length(A, B->length, ctx);
    fmpz_mpoly_sort_terms(A, ctx);
    FLINT_ASSERT(fmpz_mpoly_is_canonical(A, ctx));
}
