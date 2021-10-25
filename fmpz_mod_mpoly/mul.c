/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"
#include "ulong_extras.h"
#include "long_extras.h"


static int _try_dense(
    const fmpz * maxBfields,
    const fmpz * maxCfields,
    slong Blen,
    slong Clen,
    const mpoly_ctx_t mctx)
{
    int ret;
    const int max_bit_size = FLINT_MIN(FLINT_BITS/3 + 16, FLINT_BITS - 4);
    slong i, product_count;
    ulong * Bdegs, * Cdegs;
    ulong t, dense_size;
    TMP_INIT;

    if (mctx->nvars < 1)
        return 0;

    TMP_START;

    Bdegs = TMP_ARRAY_ALLOC(2*mctx->nvars, ulong);
    Cdegs = Bdegs + mctx->nvars;
    mpoly_get_monomial_ui_unpacked_ffmpz(Bdegs, maxBfields, mctx);
    mpoly_get_monomial_ui_unpacked_ffmpz(Cdegs, maxCfields, mctx);

    FLINT_ASSERT(Blen > 0);
    FLINT_ASSERT(Clen > 0);

    dense_size = 1;
    for (i = 0; i < mctx->nvars; i++)
    {
        if (n_add_checked(&t, Bdegs[i], Cdegs[i] + 1) || 
            n_mul_checked(&dense_size, dense_size, t))
        {
            ret = 0;
            goto cleanup;
        }
    }

    if (dense_size >= (UWORD(1) << max_bit_size))
    {
        ret = 0;
        goto cleanup;
    }

    if (z_mul_checked(&product_count, Blen, Clen))
    {
        ret = 1;
        goto cleanup;
    }

    ret = dense_size < product_count/32;

cleanup:

    TMP_END;
    return ret;
}

void fmpz_mod_mpoly_mul(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_t C,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;
    fmpz * maxBfields, * maxCfields;
    slong min_length, max_length;
    TMP_INIT;

    if (B->length < 1 || C->length < 1)
    {
        fmpz_mod_mpoly_zero(A, ctx);
        return;
    }

    TMP_START;

    maxBfields = TMP_ARRAY_ALLOC(2*ctx->minfo->nfields, fmpz);
    maxCfields = maxBfields + ctx->minfo->nfields;
    for (i = 0; i < 2*ctx->minfo->nfields; i++)
        fmpz_init(maxBfields + i);

    mpoly_max_fields_fmpz(maxBfields, B->exps, B->length, B->bits, ctx->minfo);
    mpoly_max_fields_fmpz(maxCfields, C->exps, C->length, C->bits, ctx->minfo);

    min_length = FLINT_MIN(B->length, C->length);
    max_length = FLINT_MAX(B->length, C->length);

    if (min_length < 20 || max_length < 50 ||
        B->bits > FLINT_BITS || C->bits > FLINT_BITS)
    {
        goto do_heap;
    }

    if (_try_dense(maxBfields, maxCfields, B->length, C->length, ctx->minfo))
    {
        if (_fmpz_mod_mpoly_mul_dense_maxfields(A, B, maxBfields, C, maxCfields, ctx))
            goto cleanup;
    }

do_heap:

    _fmpz_mod_mpoly_mul_johnson_maxfields(A, B, maxBfields, C, maxCfields, ctx);

cleanup:

    for (i = 0; i < 2*ctx->minfo->nfields; i++)
        fmpz_clear(maxBfields + i);

    TMP_END;
}

