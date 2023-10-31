/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_mpoly.h"

void univar_divides_check(
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx,
    slong ii,
    slong jj,
    const char * name,
    flint_rand_t randstate)
{
    int divides, udivides;
    nmod_mpoly_ctx_t uctx;
    nmod_mpolyu_t Au, Bu, Qu;
    nmod_mpoly_t Q, Qcheck;
    flint_bitcnt_t ABbits;
    ulong * shift, * stride;
    slong * perm;
    slong i, j, k;

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

    nmod_mpoly_init(Q, ctx);
    nmod_mpoly_init(Qcheck, ctx);

    nmod_mpoly_ctx_init(uctx, ctx->minfo->nvars - 1, ORD_LEX, ctx->mod.n);

    ABbits = FLINT_MAX(A->bits, B->bits);

    nmod_mpolyu_init(Au, ABbits, uctx);
    nmod_mpolyu_init(Bu, ABbits, uctx);
    nmod_mpolyu_init(Qu, ABbits, uctx);

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

    nmod_mpoly_to_mpolyu_perm_deflate_threaded_pool(Au, uctx, A, ctx,
                                           perm, shift, stride, NULL, 0);
    nmod_mpoly_to_mpolyu_perm_deflate_threaded_pool(Bu, uctx, B, ctx,
                                           perm, shift, stride, NULL, 0);

    udivides = nmod_mpolyuu_divides(Qu, Au, Bu, 1, uctx);
    divides = nmod_mpoly_divides(Q, A, B, ctx);

    if (divides != udivides)
    {
        flint_printf("check univariate divisibility\n"
                                       "i = %wd, j = %wd, %s\n", ii, jj, name);
        fflush(stdout);
        flint_abort();
    }

    if (!divides)
    {
        goto cleanup;
    }

    nmod_mpoly_from_mpolyu_perm_inflate(Qcheck, ABbits, ctx, Qu, uctx,
                                                          perm, shift, stride);

    if (!nmod_mpoly_equal(Q, Qcheck, ctx))
    {
        flint_printf("check univariate quotient\n"
                                       "i = %wd, j = %wd, %s\n", ii, jj, name);
        fflush(stdout);
        flint_abort();
    }

cleanup:

    nmod_mpolyu_clear(Au, uctx);
    nmod_mpolyu_clear(Bu, uctx);
    nmod_mpolyu_clear(Qu, uctx);

    nmod_mpoly_ctx_clear(uctx);

    nmod_mpoly_clear(Q, ctx);
    nmod_mpoly_clear(Qcheck, ctx);

    flint_free(perm);
    flint_free(shift);
    flint_free(stride);
}

TEST_FUNCTION_START(nmod_mpoly_mpolyuu_divides, state)
{
    slong i, j, tmul = 50;

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

        exp_bits = n_randint(state, FLINT_BITS/2 - 2) + 2;

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bits(a, state, len1, exp_bits, ctx);
            nmod_mpoly_randtest_bits(b, state, len2, exp_bits, ctx);
            nmod_mpoly_mul(a, a, b, ctx);
            univar_divides_check(a, b, ctx, i, j, "univar Check (a*b)/b = a", state);
        }

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
        flint_bitcnt_t exp_bits;
        mp_limb_t modulus;

        modulus = n_randint(state, (i % 10 == 0) ? 4: FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 15, modulus);

        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(c, ctx);

        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20);
        len3 = n_randint(state, 30);

        exp_bits = n_randint(state, 8) + 3;

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bits(a, state, len1, exp_bits, ctx);
            nmod_mpoly_randtest_bits(b, state, len2, exp_bits, ctx);
            nmod_mpoly_randtest_bits(c, state, len3, exp_bits, ctx);
            nmod_mpoly_mul(a, a, b, ctx);
            nmod_mpoly_add(a, a, c, ctx);
            univar_divides_check(a, b, ctx, i, j, "univar Check (a*b + c)/b", state);
        }

        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(c, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
