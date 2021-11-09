/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"
#include "long_extras.h"


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

/* !!! this function DOES need to change with new orderings */
static int _try_dense_univar(
    fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_t C,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    ulong maskB = (-UWORD(1)) >> (FLINT_BITS - B->bits);
    ulong maskC = (-UWORD(1)) >> (FLINT_BITS - C->bits);
    slong NB = mpoly_words_per_exp_sp(B->bits, ctx->minfo);
    slong NC = mpoly_words_per_exp_sp(C->bits, ctx->minfo);
    slong Blen = B->length;
    slong Clen = C->length;
    slong BClen;
    const ulong * Bexps = B->exps;
    const ulong * Cexps = C->exps;
    fmpz * Acoeffs, * Bcoeffs, * Ccoeffs;
    slong Adeg, Bdeg, Cdeg;
    TMP_INIT;

    /* the variable power is stored at offset 0 */
    Bdeg = Bexps[NB*0 + 0] & maskB;
    Cdeg = Cexps[NC*0 + 0] & maskC;

    if (z_mul_checked(&BClen, Blen, Clen) || z_add_checked(&Adeg, Bdeg, Cdeg))
        return 0;

    if (Adeg > WORD_MAX/FLINT_BITS)
        return 0;

    if (Adeg > BClen)
        return 0;

    {
        slong Bcoeffbits = _fmpz_vec_max_bits(B->coeffs, Blen);
        slong Ccoeffbits = _fmpz_vec_max_bits(C->coeffs, Clen);
        slong t = FLINT_ABS(Bcoeffbits) + FLINT_ABS(Ccoeffbits);

        if (t > FLINT_BITS && Adeg > BClen/4)
            return 0;
    }

    TMP_START;

    Acoeffs = TMP_ARRAY_ALLOC(Adeg + 1 + Bdeg + 1 + Cdeg + 1, fmpz);
    Bcoeffs = Acoeffs + Adeg + 1;
    Ccoeffs = Bcoeffs + Bdeg + 1;

    /* we own the fmpz's in Acoeffs */
    for (i = 0; i < Adeg + 1; i++)
        fmpz_init(Acoeffs + i);

    if (A != B && A != C)
    {
        for (i = FLINT_MIN(A->length - 1, Adeg); i >= 0; i--)
            fmpz_swap(Acoeffs + i, A->coeffs + i);
    }

    /* Bcoeffs and Ccoeffs are shallow copies */
    for (i = 0; i < Bdeg + 1 + Cdeg + 1; i++)
        Bcoeffs[i] = 0;

    for (i = 0; i < Blen; i++)
        Bcoeffs[Bexps[NB*i + 0] & maskB] = B->coeffs[i];

    for (i = 0; i < Clen; i++)
        Ccoeffs[Cexps[NC*i + 0] & maskC] = C->coeffs[i];

    if (Bdeg >= Cdeg)
        _fmpz_poly_mul(Acoeffs, Bcoeffs, Bdeg + 1, Ccoeffs, Cdeg + 1);
    else
        _fmpz_poly_mul(Acoeffs, Ccoeffs, Cdeg + 1, Bcoeffs, Bdeg + 1);

    /* Acoeffs cleared */
    _fmpz_mpoly_set_fmpz_poly_one_var(A, FLINT_MAX(B->bits, C->bits),
                                                           Acoeffs, Adeg, ctx);
    TMP_END;
    return 1;
}


void fmpz_mpoly_mul(
    fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_t C,
    const fmpz_mpoly_ctx_t ctx)
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

    if (B->length < 1 || C->length < 1)
    {
        fmpz_mpoly_zero(A, ctx);
        return;
    }
    else if (B->length == 1)
    {
        fmpz_mpoly_mul_monomial(A, C, B, ctx);
        return;
    }
    else if (C->length == 1)
    {
        fmpz_mpoly_mul_monomial(A, B, C, ctx);
        return;
    }

    FLINT_ASSERT(nvars > 0);

    if (nvars == 1 && B->bits <= FLINT_BITS && C->bits <= FLINT_BITS)
    {
        if (_try_dense_univar(A, B, C, ctx))
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
    if (min_length < 20 || max_length < 50)
    {
        _fmpz_mpoly_mul_johnson_maxfields(A, B, maxBfields, C, maxCfields, ctx);
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
        success = _fmpz_mpoly_mul_dense(A, B, maxBfields, C, maxCfields, ctx);
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
                ? _fmpz_mpoly_mul_array_threaded_pool_LEX(
                                    A, B, maxBfields, C, maxCfields, ctx,
                                                         handles, num_handles)
                : _fmpz_mpoly_mul_array_LEX(
                                    A, B, maxBfields, C, maxCfields, ctx);
    }
    else if (ctx->minfo->ord == ORD_DEGLEX || ctx->minfo->ord == ORD_DEGREVLEX)
    {
        success = (num_handles > 0)
                ? _fmpz_mpoly_mul_array_threaded_pool_DEG(
                                    A, B, maxBfields, C, maxCfields, ctx,
                                                         handles, num_handles)
                : _fmpz_mpoly_mul_array_DEG(
                                    A, B, maxBfields, C, maxCfields, ctx);
    }

    if (success)
    {
        goto cleanup_threads;
    }

do_heap:

    if (num_handles > 0)
    {
        _fmpz_mpoly_mul_heap_threaded_pool_maxfields(A,
                      B, maxBfields, C, maxCfields, ctx, handles, num_handles);
    }
    else
    {
        _fmpz_mpoly_mul_johnson_maxfields(A, B, maxBfields, C, maxCfields, ctx);
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

