/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

int fmpq_mpoly_is_fmpq_poly(
    const fmpq_mpoly_t A,
    slong var,
    const fmpq_mpoly_ctx_t ctx)
{
    return fmpz_mpoly_is_fmpz_poly(A->zpoly, var, ctx->zctx);
}

int fmpq_mpoly_get_fmpq_poly(
    fmpq_poly_t A,
    const fmpq_mpoly_t B,
    slong var,
    const fmpq_mpoly_ctx_t ctx)
{
    int success;
    fmpz_poly_t a;
    fmpz_poly_init(a);
    success = fmpz_mpoly_get_fmpz_poly(a, B->zpoly, var, ctx->zctx);
    if (success)
    {
        fmpq_poly_set_fmpz_poly(A, a);
        fmpq_poly_scalar_mul_fmpq(A, A, B->content);
    }
    fmpz_poly_clear(a);
    return success;
}

void fmpq_mpoly_set_fmpq_poly(
    fmpq_mpoly_t A,
    const fmpq_poly_t B,
    slong v,
    const fmpq_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bits;

    if (B->length < 1)
    {
        fmpq_mpoly_zero(A, ctx);
        return;
    }

    bits = mpoly_gen_pow_exp_bits_required(v, B->length - 1, ctx->zctx->minfo);
    bits = mpoly_fix_bits(bits, ctx->zctx->minfo);
    _fmpz_mpoly_set_fmpz_poly(A->zpoly, bits, B->coeffs, B->length, v, ctx->zctx);

    fmpz_one(fmpq_numref(A->content));
    fmpz_set(fmpq_denref(A->content), B->den);

    fmpq_mpoly_reduce(A, ctx);
}
