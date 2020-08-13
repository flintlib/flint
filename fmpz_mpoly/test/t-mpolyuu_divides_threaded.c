/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void bivar_divides_check(
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx,
    slong ii,
    slong jj,
    const char * name,
    flint_rand_t randstate)
{
    int divides, uudivides;
    fmpz_mpoly_ctx_t uuctx;
    fmpz_mpolyu_t Auu, Buu, Quu;
    fmpz_mpoly_t Q, Qcheck;
    flint_bitcnt_t ABbits;
    slong * Adegs, * Bdegs;
    ulong * shift, * stride;
    slong * perm;
    slong i, j, k;
    thread_pool_handle * handles;
    slong num_workers;

    if (   A->bits > FLINT_BITS
        || B->bits > FLINT_BITS
        || A->length == 0
        || B->length == 0
        || ctx->minfo->nvars < 3)
    {
        return;
    }

    Adegs = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    Bdegs = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    perm = (slong *) flint_malloc((ctx->minfo->nvars)*sizeof(slong));
    shift = (ulong *) flint_malloc((ctx->minfo->nvars)*sizeof(ulong));
    stride = (ulong *) flint_malloc((ctx->minfo->nvars)*sizeof(ulong));

    fmpz_mpoly_init(Q, ctx);
    fmpz_mpoly_init(Qcheck, ctx);

    mpoly_degrees_si(Adegs, A->exps, A->length, A->bits, ctx->minfo);
    mpoly_degrees_si(Bdegs, B->exps, B->length, B->bits, ctx->minfo);

    fmpz_mpoly_ctx_init(uuctx, ctx->minfo->nvars - 2, ORD_LEX);

    ABbits = FLINT_MAX(A->bits, B->bits);

    fmpz_mpolyu_init(Auu, ABbits, uuctx);
    fmpz_mpolyu_init(Buu, ABbits, uuctx);
    fmpz_mpolyu_init(Quu, ABbits, uuctx);

    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        perm[i] = i;
        shift[i] = 0;
        stride[i] = 1;
    }

    /* randomize perm */
    for (k = 0; k < ctx->minfo->nvars; k++)
    {
        slong t1, t2;
        i = n_randint(randstate, ctx->minfo->nvars - 1);
        j = i + n_randint(randstate, ctx->minfo->nvars - i);
        t1 = perm[i];
        t2 = perm[j];
        perm[i] = t2;
        perm[j] = t1;
    }

    /* ensure that main two variables can be packed into FLINT_BITS/2 */
    if (   FLINT_BIT_COUNT(Adegs[perm[0]]) >= FLINT_BITS/2
        || FLINT_BIT_COUNT(Adegs[perm[1]]) >= FLINT_BITS/2
        || FLINT_BIT_COUNT(Bdegs[perm[0]]) >= FLINT_BITS/2
        || FLINT_BIT_COUNT(Bdegs[perm[1]]) >= FLINT_BITS/2)
    {
        goto cleanup;
    }

    num_workers = flint_request_threads(&handles, 9999);

    fmpz_mpoly_to_mpolyuu_perm_deflate_threaded_pool(Auu, uuctx, A, ctx,
                              perm, shift, stride, NULL, handles, num_workers);
    fmpz_mpoly_to_mpolyuu_perm_deflate_threaded_pool(Buu, uuctx, B, ctx,
                              perm, shift, stride, NULL, handles, num_workers);

    uudivides = fmpz_mpolyuu_divides_threaded_pool(Quu, Auu, Buu, 2, uuctx,
                                                         handles, num_workers);
    flint_give_back_threads(handles, num_workers);

    divides = fmpz_mpoly_divides(Q, A, B, ctx);

    if (divides != uudivides)
    {
        flint_printf("check bivariate divisibility\n"
                                       "i = %wd, j = %wd, %s\n", ii, jj, name);
        flint_abort();
    }

    if (!divides)
    {
        goto cleanup;
    }

    fmpz_mpoly_from_mpolyuu_perm_inflate(Qcheck, ABbits, ctx, Quu, uuctx,
                                                          perm, shift, stride);
    if (!fmpz_mpoly_equal(Q, Qcheck, ctx))
    {
        flint_printf("check bivariate quotient\n"
                                       "i = %wd, j = %wd, %s\n", ii, jj, name);
        flint_abort();
    }

cleanup:

    fmpz_mpolyu_clear(Auu, uuctx);
    fmpz_mpolyu_clear(Buu, uuctx);
    fmpz_mpolyu_clear(Quu, uuctx);

    fmpz_mpoly_ctx_clear(uuctx);

    fmpz_mpoly_clear(Q, ctx);
    fmpz_mpoly_clear(Qcheck, ctx);

    flint_free(Adegs);
    flint_free(Bdegs);
    flint_free(perm);
    flint_free(shift);
    flint_free(stride);
}

void univar_divides_check(
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx,
    slong ii,
    slong jj,
    const char * name,
    flint_rand_t randstate)
{
    int divides, udivides;
    fmpz_mpoly_ctx_t uctx;
    fmpz_mpolyu_t Au, Bu, Qu;
    fmpz_mpoly_t Q, Qcheck;
    flint_bitcnt_t ABbits;
    ulong * shift, * stride;
    slong * perm;
    slong i, j, k;
    thread_pool_handle * handles;
    slong num_workers;

    if (   A->bits > FLINT_BITS
        || B->bits > FLINT_BITS
        || A->length == 0
        || B->length == 0
        || ctx->minfo->nvars < 2)
    {
        return;
    }

    perm = (slong *) flint_malloc((ctx->minfo->nvars)*sizeof(slong));
    shift = (ulong *) flint_malloc((ctx->minfo->nvars)*sizeof(ulong));
    stride = (ulong *) flint_malloc((ctx->minfo->nvars)*sizeof(ulong));

    fmpz_mpoly_init(Q, ctx);
    fmpz_mpoly_init(Qcheck, ctx);

    fmpz_mpoly_ctx_init(uctx, ctx->minfo->nvars - 1, ORD_LEX);

    ABbits = FLINT_MAX(A->bits, B->bits);

    fmpz_mpolyu_init(Au, ABbits, uctx);
    fmpz_mpolyu_init(Bu, ABbits, uctx);
    fmpz_mpolyu_init(Qu, ABbits, uctx);

    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        perm[i] = i;
        shift[i] = 0;
        stride[i] = 1;
    }

    /* randomize perm */
    for (k = 0; k < ctx->minfo->nvars; k++)
    {
        slong t1, t2;
        i = n_randint(randstate, ctx->minfo->nvars - 1);
        j = i + n_randint(randstate, ctx->minfo->nvars - i);
        t1 = perm[i];
        t2 = perm[j];
        perm[i] = t2;
        perm[j] = t1;
    }

    num_workers = flint_request_threads(&handles, 9999);

    fmpz_mpoly_to_mpolyu_perm_deflate_threaded_pool(Au, uctx, A, ctx,
                              perm, shift, stride, NULL, handles, num_workers);
    fmpz_mpoly_to_mpolyu_perm_deflate_threaded_pool(Bu, uctx, B, ctx,
                              perm, shift, stride, NULL, handles, num_workers);

    udivides = fmpz_mpolyuu_divides_threaded_pool(Qu, Au, Bu, 1, uctx,
                                                         handles, num_workers);

    flint_give_back_threads(handles, num_workers);

    divides = fmpz_mpoly_divides(Q, A, B, ctx);

    if (divides != udivides)
    {
        flint_printf("check univariate divisibility\n"
                                       "i = %wd, j = %wd, %s\n", ii, jj, name);
        flint_abort();
    }

    if (!divides)
    {
        goto cleanup;
    }

    fmpz_mpoly_from_mpolyu_perm_inflate(Qcheck, ABbits, ctx, Qu, uctx,
                                                          perm, shift, stride);

    if (!fmpz_mpoly_equal(Q, Qcheck, ctx))
    {
        flint_printf("check univariate quotient\n"
                                       "i = %wd, j = %wd, %s\n", ii, jj, name);
        flint_abort();
    }

cleanup:

    fmpz_mpolyu_clear(Au, uctx);
    fmpz_mpolyu_clear(Bu, uctx);
    fmpz_mpolyu_clear(Qu, uctx);

    fmpz_mpoly_ctx_clear(uctx);

    fmpz_mpoly_clear(Q, ctx);
    fmpz_mpoly_clear(Qcheck, ctx);

    flint_free(perm);
    flint_free(shift);
    flint_free(stride);
}


int
main(void)
{
    slong i, j, max_threads = 5, tmul = 20;
    FLINT_TEST_INIT(state);
#ifdef _WIN32
    tmul = 2;
#endif

    flint_printf("mpolyuu_divides_threaded....");
    fflush(stdout);

    /* Check (a*b)/b */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b;
        flint_bitcnt_t coeff_bits;
        slong len1, len2;
        flint_bitcnt_t exp_bits;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);

        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);

        len1 = n_randint(state, 30);
        len2 = n_randint(state, 30);

        exp_bits = n_randint(state, FLINT_BITS/2 - 2) + 2;

        coeff_bits = n_randint(state, 100);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bits(a, state, len1, coeff_bits, exp_bits, ctx);
            fmpz_mpoly_randtest_bits(b, state, len2, coeff_bits, exp_bits, ctx);
            fmpz_mpoly_mul(a, a, b, ctx);
            univar_divides_check(a, b, ctx, i, j, "univar Check (a*b)/b", state);
            bivar_divides_check(a, b, ctx, i, j, "bivar Check (a*b)/b", state);
        }

        flint_set_num_threads(n_randint(state, max_threads) + 1);

        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check (a*b + c)/b */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, c;
        flint_bitcnt_t coeff_bits;
        slong len1, len2, len3;
        flint_bitcnt_t exp_bits;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);

        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(c, ctx);

        len1 = n_randint(state, 30);
        len2 = n_randint(state, 30);
        len3 = n_randint(state, 40);

        exp_bits = n_randint(state, 7) + 3;

        coeff_bits = n_randint(state, 100);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bits(a, state, len1, coeff_bits, exp_bits, ctx);
            fmpz_mpoly_randtest_bits(b, state, len2, coeff_bits, exp_bits, ctx);
            fmpz_mpoly_randtest_bits(c, state, len3, coeff_bits, exp_bits, ctx);
            fmpz_mpoly_mul(a, a, b, ctx);
            fmpz_mpoly_add(a, a, c, ctx);
            flint_set_num_threads(n_randint(state, max_threads) + 1);
            univar_divides_check(a, b, ctx, i, j, "univar Check (a*b + c)/b", state);
            bivar_divides_check(a, b, ctx, i, j, "bivar Check (a*b + c)/b", state);
        }

        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(c, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}
