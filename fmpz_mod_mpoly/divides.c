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

static int _fmpz_mod_mpoly_divides_try_dense(
    const fmpz * maxAfields,
    const fmpz * maxBfields,
    slong Alen,
    slong Blen,
    const mpoly_ctx_t mctx)
{
    const int max_bit_size = FLINT_MIN(FLINT_BITS/3 + 16, FLINT_BITS - 4);
    slong i;
    ulong * Adegs;
    ulong total_dense_size;
    TMP_INIT;

    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(Blen > 0);

    if (mctx->nvars < 1 || mctx->nvars > FLINT_BITS)
        return 0;

    TMP_START;

    Adegs = TMP_ARRAY_ALLOC(mctx->nvars, ulong);
    mpoly_get_monomial_ui_unpacked_ffmpz(Adegs, maxAfields, mctx);

    total_dense_size = 1;
    for (i = 0; i < mctx->nvars; i++)
    {
        if (n_mul_checked(&total_dense_size, total_dense_size, Adegs[i] + 1))
        {
            TMP_END;
            return 0;
        }
    }

    TMP_END;

    return total_dense_size < (UWORD(1) << max_bit_size) &&
           total_dense_size/16 < Alen;
}


int fmpz_mod_mpoly_divides(
    fmpz_mod_mpoly_t Q,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int success;
    slong i;
    fmpz * maxAfields, * maxBfields;
    TMP_INIT;

    if (fmpz_mod_mpoly_is_zero(B, ctx))
    {
        if (!fmpz_mod_mpoly_is_zero(A, ctx) &&
            !fmpz_is_one(fmpz_mod_mpoly_ctx_modulus(ctx)))
        {
            flint_throw(FLINT_DIVZERO, "fmpz_mod_mpoly_divides: divide by zero");
        }

        fmpz_mod_mpoly_zero(Q, ctx);
        return 1;
    }

    if (fmpz_mod_mpoly_is_zero(A, ctx))
    {
        fmpz_mod_mpoly_zero(Q, ctx);
        return 1;
    }

    TMP_START;

    maxAfields = TMP_ARRAY_ALLOC(2*ctx->minfo->nfields, fmpz);
    maxBfields = maxAfields + ctx->minfo->nfields;
    for (i = 0; i < 2*ctx->minfo->nfields; i++)
        fmpz_init(maxAfields + i);

    mpoly_max_fields_fmpz(maxAfields, A->exps, A->length, A->bits, ctx->minfo);
    mpoly_max_fields_fmpz(maxBfields, B->exps, B->length, B->bits, ctx->minfo);

    /* quick degree check */
    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        if (fmpz_cmp(maxAfields + i, maxBfields + i) < 0)
        {
            fmpz_mod_mpoly_zero(Q, ctx);
            success = 0;
            goto cleanup;
        }
    }

    if (A->length < 30 || A->bits > FLINT_BITS || B->bits > FLINT_BITS)
    {
        goto do_heap;
    }

    if (_fmpz_mod_mpoly_divides_try_dense(maxAfields, maxBfields,
                                             A->length, B->length, ctx->minfo))
    {
        success = _fmpz_mod_mpoly_divides_dense_maxfields(Q,
                                            A, maxAfields, B, maxBfields, ctx);
        if (success >= 0)
            goto cleanup;
    }

do_heap:

    success = _fmpz_mod_mpoly_divides_monagan_pearce_maxfields(Q,
                                            A, maxAfields, B, maxBfields, ctx);
cleanup:

    for (i = 0; i < 2*ctx->minfo->nfields; i++)
        fmpz_clear(maxAfields + i);

    TMP_END;
    return success;
}
