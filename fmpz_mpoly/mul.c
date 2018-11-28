/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

static int _try_dense(slong * Adegs, slong * Bdegs,
                                           slong Alen, slong Blen, slong nvars)
{
    slong i, total_dense_size;
    ulong hi;

    FLINT_ASSERT(Alen > WORD(0));
    FLINT_ASSERT(Blen > WORD(0));

    total_dense_size = WORD(1);
    for (i = 0; i < nvars; i++)
    {
        if (Adegs[i] + Bdegs[i] + 1 < WORD(0))
            return 0;

        umul_ppmm(hi, total_dense_size, total_dense_size, Adegs[i] + Bdegs[i] + 1);

        if (hi != WORD(0) || total_dense_size <= WORD(0))
            return 0;
    }

    return total_dense_size < WORD(5000000)
            && total_dense_size/Alen/Blen < WORD(10);
}

static int _try_array_LEX(slong * Adegs, slong * Bdegs,
                                           slong Alen, slong Blen, slong nvars)
{
    slong i, total_dense_size;
    ulong hi;

    FLINT_ASSERT(Alen > WORD(0));
    FLINT_ASSERT(Blen > WORD(0));

    /* accept array method if the array is probably at least 10% full */

    total_dense_size = WORD(1);
    for (i = 0; i < nvars; i++)
    {
        umul_ppmm(hi, total_dense_size, total_dense_size, Adegs[i] + Bdegs[i] + 1);
        if (hi != WORD(0) || total_dense_size <= WORD(0))
            return 0;
    }
    return total_dense_size <= WORD(5000000)
            && total_dense_size/Alen/Blen < WORD(10);
}

static int _try_array_DEG(slong * Adegs, slong * Bdegs,
                                           slong Alen, slong Blen, slong nvars)
{
    slong i, total_dense_size, total_degree;
    ulong hi;

    FLINT_ASSERT(Alen > WORD(0));
    FLINT_ASSERT(Blen > WORD(0));

    total_degree = 0;
    for (i = 0; i < nvars; i++)
    {
        total_degree += Adegs[i];
        if (total_degree < 0)
            return 0;
        total_degree += Bdegs[i];
        if (total_degree < 0)
            return 0;
    }

    /* the relevent portion of the array has approx size d^nvars/nvars!*/
    total_dense_size = WORD(1);
    for (i = 0; i < nvars; i++)
    {
        umul_ppmm(hi, total_dense_size, total_dense_size, total_degree);
        if (hi != WORD(0) || total_dense_size < WORD(0))
            return 0;
    }
    for (i = 0; i < nvars; i++)
    {
        total_dense_size /= i+1;
    }

    return total_dense_size <= WORD(5000000)
            && total_dense_size/Alen/Blen < WORD(10);
}


void fmpz_mpoly_mul(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                             const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    int success;
    slong * Bdegs, * Cdegs;
    fmpz * maxBfields, * maxCfields;
    TMP_INIT;

    if (B->length == 0 || C->length == 0)
    {
        fmpz_mpoly_zero(A, ctx);
        return;
    }

    TMP_START;

    /*
        All methods require a linear scan of the exponents.
        Do it here once and for all.
    */
    maxBfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    maxCfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_init(maxBfields + i);
        fmpz_init(maxCfields + i);
    }
    mpoly_max_fields_fmpz(maxBfields, B->exps, B->length, B->bits, ctx->minfo);
    mpoly_max_fields_fmpz(maxCfields, C->exps, C->length, C->bits, ctx->minfo);


    if (B->bits > FLINT_BITS || C->bits > FLINT_BITS
             || B->length < 50 || C->length < 50)
    {
        _fmpz_mpoly_mul_johnson_maxfields(A, B, maxBfields, C, maxCfields, ctx);
        goto done;
    }

    /*
        The multiplication is not trivial and each packed field fits
        into one word. In particular, the degs must fit an slong.
    */

    Bdegs = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    Cdegs = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    mpoly_get_monomial_ui_unpacked_ffmpz((ulong *)Bdegs, maxBfields, ctx->minfo);
    mpoly_get_monomial_ui_unpacked_ffmpz((ulong *)Cdegs, maxCfields, ctx->minfo);

    success = 0;
    if (_try_dense(Bdegs, Cdegs, B->length, C->length, ctx->minfo->nvars))
    {
        success = _fmpz_mpoly_mul_dense(A, B, maxBfields, C, maxCfields, ctx);
        if (success)
        {
            goto done;
        }
    }

    if (   ctx->minfo->nvars <= WORD(1)
        || ctx->minfo->nvars >= WORD(8)
        || 1 != mpoly_words_per_exp(B->bits, ctx->minfo)
        || 1 != mpoly_words_per_exp(C->bits, ctx->minfo)
       )
    {
        goto do_heap;
    }

    switch (ctx->minfo->ord)
    {
        case ORD_LEX:
        {
            if (_try_array_LEX(Bdegs, Cdegs, B->length, C->length,
                                                            ctx->minfo->nvars))
            {
                success = _fmpz_mpoly_mul_array_threaded_LEX(A,
                                            B, maxBfields, C, maxCfields, ctx);
            }
            break;
        }
        case ORD_DEGREVLEX:
        case ORD_DEGLEX:
        {
            if (_try_array_DEG(Bdegs, Cdegs, B->length, C->length,
                                                            ctx->minfo->nvars))
            {
                success = _fmpz_mpoly_mul_array_threaded_DEG(A,
                                            B, maxBfields, C, maxCfields, ctx);
            }
            break;
        }
        default:
        {
            success = 0;
            break;
        }
    }

    if (success)
    {
        goto done;
    }

do_heap:

    _fmpz_mpoly_mul_heap_threaded_maxfields(A, B, maxBfields, C, maxCfields, ctx);

done:

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(maxBfields + i);
        fmpz_clear(maxCfields + i);
    }

    TMP_END;
}
