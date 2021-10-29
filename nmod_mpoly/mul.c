/*
    Copyright (C) 2018, 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

static int _try_dense(int try_array, slong * Bdegs, slong * Cdegs,
                                           slong Blen, slong Clen, slong nvars)
{
    const int max_bit_size = FLINT_MIN(FLINT_BITS/3 + 16, FLINT_BITS - 3);
    slong i, product_count, dense_size;
    ulong hi;

    FLINT_ASSERT(Blen > 0);
    FLINT_ASSERT(Clen > 0);

    dense_size = WORD(1);
    for (i = 0; i < nvars; i++)
    {
        umul_ppmm(hi, dense_size, dense_size, Bdegs[i] + Cdegs[i] + 1);

        if (hi != 0 || dense_size <= 0)
            return 0;
    }

    if (dense_size >= WORD(1) << max_bit_size)
        return 0;

    umul_ppmm(hi, product_count, Blen, Clen);

    if (hi != 0 || product_count < 0)
        return 1;

    /*
        Assume that the running time of the dense method is linear
        in "dense_size" and that the running time of the array|heap
        method is linear in "product_count".
        Assume further that the array method is 4x faster than heap.
    */

    if (try_array)
        return dense_size < product_count/128;
    else
        return dense_size < product_count/32;
}


static int _try_array_LEX(slong * Bdegs, slong * Cdegs,
                                           slong Blen, slong Clen, slong nvars)
{
    slong i, dense_size;
    ulong hi;

    FLINT_ASSERT(Blen > 0);
    FLINT_ASSERT(Clen > 0);

    /* accept array method if the array is probably at least 10% full */

    dense_size = WORD(1);
    for (i = 0; i < nvars; i++)
    {
        umul_ppmm(hi, dense_size, dense_size, Bdegs[i] + Cdegs[i] + 1);
        if (hi != 0 || dense_size <= 0)
            return 0;
    }

    return dense_size <= WORD(50000000) &&
           dense_size/Blen/Clen < WORD(10);
}


static int _try_array_DEG(slong Btotaldeg, slong Ctotaldeg,
                                           slong Blen, slong Clen, slong nvars)
{
    slong i, dense_size, total_degree;
    ulong hi;

    FLINT_ASSERT(Blen > 0);
    FLINT_ASSERT(Clen > 0);

    total_degree = Btotaldeg + Btotaldeg;
    if (total_degree <= 0)
        return 0;

    /* the relevent portion of the array has approx size d^nvars/nvars!*/
    dense_size = WORD(1);
    for (i = 0; i < nvars; i++)
    {
        umul_ppmm(hi, dense_size, dense_size, total_degree);
        if (hi != WORD(0) || dense_size < 0)
            return 0;
    }
    for (i = 0; i < nvars; i++)
    {
        dense_size /= i + 1;
    }

    return dense_size <= WORD(5000000) &&
           dense_size/Blen/Clen < WORD(10);
}


void nmod_mpoly_mul(
    nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_t C,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong nvars = ctx->minfo->nvars;
    int success, try_array;
    slong * Bdegs, * Cdegs;
    fmpz * maxBfields, * maxCfields;
    thread_pool_handle * handles;
    slong num_handles;
    slong min_length, max_length;
    slong thread_limit;
    TMP_INIT;

    if (B->length == 0 || C->length == 0)
    {
        nmod_mpoly_zero(A, ctx);
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

    min_length = FLINT_MIN(B->length, C->length);
    max_length = FLINT_MAX(B->length, C->length);
    thread_limit = min_length/512;

    /*
        If one polynomial is tiny or if both polynomials are small,
        heap method with operational complexity O(B->length*C->length) is fine.
    */
    if (nvars < 1 || min_length < 20 || max_length < 50)
    {
        _nmod_mpoly_mul_johnson_maxfields(A, B, maxBfields, C, maxCfields, ctx);
        goto cleanup;
    }

    /*
        If either polynomial has multi-word fields, only heap will do.
    */
    if (B->bits > FLINT_BITS || C->bits > FLINT_BITS)
    {
        num_handles = flint_request_threads(&handles, thread_limit);
        goto do_heap;
    }

    /*
        The multiplication is not trivial and each packed field fits
        into one word. In particular, the degrees must fit an slong.
    */
    Bdegs = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    Cdegs = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    mpoly_get_monomial_ui_unpacked_ffmpz((ulong *)Bdegs, maxBfields, ctx->minfo);
    mpoly_get_monomial_ui_unpacked_ffmpz((ulong *)Cdegs, maxCfields, ctx->minfo);

    /*
        See if array method is applicable.
        If so, it should be about 4x faster than heap.
    */
    try_array = 0;
    if (nvars > WORD(1) &&
        nvars < WORD(8) &&
        1 == mpoly_words_per_exp(B->bits, ctx->minfo) &&
        1 == mpoly_words_per_exp(C->bits, ctx->minfo))
    {
        if (ctx->minfo->ord == ORD_LEX)
        {
            try_array = _try_array_LEX(Bdegs, Cdegs, B->length, C->length, nvars);
        }
        else if (ctx->minfo->ord == ORD_DEGLEX || ctx->minfo->ord == ORD_DEGREVLEX)
        {
            slong Btdeg = fmpz_get_si(maxBfields + nvars);
            slong Ctdeg = fmpz_get_si(maxCfields + nvars);
            try_array = _try_array_DEG(Btdeg, Ctdeg, B->length, C->length, nvars);
        }
    }

    success = 0;
    if (_try_dense(try_array, Bdegs, Cdegs, B->length, C->length, nvars))
    {
        success = _nmod_mpoly_mul_dense(A, B, maxBfields, C, maxCfields, ctx);
        if (success)
        {
            goto cleanup;
        }
    }

    num_handles = flint_request_threads(&handles, thread_limit);

    if (!try_array)
    {
        goto do_heap;
    }

    if (ctx->minfo->ord == ORD_LEX)
    {
        success = (num_handles > 0)
                ? _nmod_mpoly_mul_array_threaded_pool_LEX(
                                    A, B, maxBfields, C, maxCfields, ctx,
                                                         handles, num_handles)
                : _nmod_mpoly_mul_array_LEX(
                                    A, B, maxBfields, C, maxCfields, ctx);
    }
    else if (ctx->minfo->ord == ORD_DEGLEX || ctx->minfo->ord == ORD_DEGREVLEX)
    {
        success = (num_handles > 0)
                ? _nmod_mpoly_mul_array_threaded_pool_DEG(
                                    A, B, maxBfields, C, maxCfields, ctx,
                                                         handles, num_handles)
                : _nmod_mpoly_mul_array_DEG(
                                    A, B, maxBfields, C, maxCfields, ctx);
    }

    if (success)
    {
        goto cleanup_threads;
    }

do_heap:

    if (num_handles > 0)
    {
        _nmod_mpoly_mul_heap_threaded_pool_maxfields(A,
                      B, maxBfields, C, maxCfields, ctx, handles, num_handles);
    }
    else
    {
        _nmod_mpoly_mul_johnson_maxfields(A, B, maxBfields, C, maxCfields, ctx);
    }

cleanup_threads:

    flint_give_back_threads(handles, num_handles);

cleanup:

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(maxBfields + i);
        fmpz_clear(maxCfields + i);
    }

    TMP_END;
}

