/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "thread_support.h"
#include "nmod_mpoly.h"

#if defined(nmod_mpolyn_divides_threaded_pool)
void _divides_check(
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx,
    slong ii,
    slong jj,
    const char * name,
    flint_rand_t randstate)
{
    int divides, ndivides;
    nmod_mpoly_ctx_t nctx;
    nmod_mpolyn_t An, Bn, Qn;
    nmod_mpoly_t Q, Qcheck;
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
        || ctx->minfo->nvars <= 2)
    {
        return;
    }

    perm = (slong *) flint_malloc((ctx->minfo->nvars)*sizeof(slong));
    shift = (ulong *) flint_malloc((ctx->minfo->nvars)*sizeof(ulong));
    stride = (ulong *) flint_malloc((ctx->minfo->nvars)*sizeof(ulong));

    nmod_mpoly_init(Q, ctx);
    nmod_mpoly_init(Qcheck, ctx);

    nmod_mpoly_ctx_init(nctx, ctx->minfo->nvars, ORD_LEX, ctx->mod.n);

    ABbits = FLINT_MAX(A->bits, B->bits);

    nmod_mpolyn_init(An, ABbits, nctx);
    nmod_mpolyn_init(Bn, ABbits, nctx);
    nmod_mpolyn_init(Qn, ABbits, nctx);

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

    nmod_mpoly_to_mpolyn_perm_deflate_threaded_pool(An, nctx, A, ctx,
                                           perm, shift, stride, NULL, 0);
    nmod_mpoly_to_mpolyn_perm_deflate_threaded_pool(Bn, nctx, B, ctx,
                                           perm, shift, stride, NULL, 0);

    num_workers = flint_request_threads(&handles, 9999);

    ndivides = nmod_mpolyn_divides_threaded_pool(Qn, An, Bn, nctx, handles, num_workers);

    flint_give_back_threads(handles, num_workers);

    divides = nmod_mpoly_divides(Q, A, B, ctx);

    if (divides != ndivides)
    {
        flint_printf("check divisibility i = %wd, j = %wd, %s\n", ii, jj, name);
        fflush(stdout);
        flint_abort();
    }

    if (!divides)
    {
        goto cleanup;
    }

    nmod_mpoly_from_mpolyn_perm_inflate(Qcheck, ABbits, ctx, Qn, nctx,
                                                          perm, shift, stride);

    if (!nmod_mpoly_equal(Q, Qcheck, ctx))
    {
        flint_printf("check quotient i = %wd, j = %wd, %s\n", ii, jj, name);
        fflush(stdout);
        flint_abort();
    }

cleanup:

    nmod_mpolyn_clear(An, nctx);
    nmod_mpolyn_clear(Bn, nctx);
    nmod_mpolyn_clear(Qn, nctx);

    nmod_mpoly_ctx_clear(nctx);

    nmod_mpoly_clear(Q, ctx);
    nmod_mpoly_clear(Qcheck, ctx);

    flint_free(perm);
    flint_free(shift);
    flint_free(stride);
}

TEST_FUNCTION_START(nmod_mpolyn_divides_threaded_pool, state)
{
    slong i, j, max_threads = 5, tmul = 50;
#ifdef _WIN32
    tmul = 2;
#endif

    /* Check (a*b)/b = a */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b;
        slong len1, len2;
        flint_bitcnt_t exp_bits;
        mp_limb_t modulus;

        modulus = n_randint(state, (i % 10 == 0) ? 4: FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 20, modulus);

        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);

        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20);

        exp_bits = n_randint(state, 8) + 3;

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bits(a, state, len1, exp_bits, ctx);
            nmod_mpoly_randtest_bits(b, state, len2, exp_bits, ctx);
            nmod_mpoly_mul(a, a, b, ctx);
            _divides_check(a, b, ctx, i, j, "check (a*b)/b = a", state);
        }

        flint_set_num_threads(n_randint(state, max_threads) + 1);

        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check (a*b + c)/b */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, c;
        slong len1, len2, len3;
        ulong exp_bound, exp_bound2;
        mp_limb_t modulus;

        modulus = n_randint(state, (i % 10 == 0) ? 4: FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 20, modulus);
        if (ctx->minfo->nvars < 1)
        {
            nmod_mpoly_ctx_clear(ctx);
            continue;
        }

        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(c, ctx);

        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20);
        len3 = n_randint(state, 30);

        exp_bound = 1 + 200/ctx->minfo->nvars;
        exp_bound2 = 1 + 100/ctx->minfo->nvars;

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bound(a, state, len1, exp_bound, ctx);
            nmod_mpoly_randtest_bound(b, state, len2, exp_bound, ctx);
            nmod_mpoly_randtest_bound(c, state, len3, exp_bound2, ctx);
            nmod_mpoly_mul(a, a, b, ctx);
            nmod_mpoly_add(a, a, c, ctx);
            _divides_check(a, b, ctx, i, j, "check (a*b + c)/b", state);
        }

        flint_set_num_threads(n_randint(state, max_threads) + 1);

        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(c, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
#else
TEST_FUNCTION_START(nmod_mpolyn_divides_threaded_pool, state)
{
    TEST_FUNCTION_END_SKIPPED(state);
}
#endif
