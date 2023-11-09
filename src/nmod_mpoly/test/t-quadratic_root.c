/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_mpoly.h"

void _test_root(
    nmod_mpoly_t x,
    const nmod_mpoly_t a,
    const nmod_mpoly_t b,
    const nmod_mpoly_ctx_t ctx,
    int sol_exists)
{
    int success, success2;
    nmod_mpoly_t s, t;

    nmod_mpoly_init(s, ctx);
    nmod_mpoly_init(t, ctx);

    success = nmod_mpoly_quadratic_root(x, a, b, ctx);

    if (sol_exists && !success)
    {
        flint_printf("FAIL: solution exists but root failed\n");
        fflush(stdout);
        flint_abort();
    }

    if (success)
    {
        nmod_mpoly_add(t, x, a, ctx);
        nmod_mpoly_mul(s, t, x, ctx);
        if (!nmod_mpoly_equal(s, b, ctx))
        {

            flint_printf("FAIL: reported solution is not a solution\n");
            fflush(stdout);
            flint_abort();
        }
    }

    nmod_mpoly_set(t, a, ctx);
    success2 = nmod_mpoly_quadratic_root(t, t, b, ctx);
    if (success != success2 || (success && !nmod_mpoly_equal(x, t, ctx)))
    {
        flint_printf("FAIL: Check aliasing first argument\n");
        fflush(stdout);
        flint_abort();
    }

    nmod_mpoly_set(t, b, ctx);
    success2 = nmod_mpoly_quadratic_root(t, a, t, ctx);
    if (success != success2 || (success && !nmod_mpoly_equal(x, t, ctx)))
    {
        flint_printf("FAIL: Check aliasing second argument\n");
        fflush(stdout);
        flint_abort();
    }

    nmod_mpoly_clear(s, ctx);
    nmod_mpoly_clear(t, ctx);
}

TEST_FUNCTION_START(nmod_mpoly_quadratic_root, state)
{
    slong i, j, tmul = 20;

    /* Check sqrt(f^2) = +-f */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, a, b, x;
        slong len, len1;
        flint_bitcnt_t exp_bits, exp_bits1;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);
        if (i % 2)
            modulus = 2;

        nmod_mpoly_ctx_init_rand(ctx, state, 10, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(x, ctx);

        for (j = 0; j < 5; j++)
        {
            len = n_randint(state, 100);
            len1 = n_randint(state, 100) + 1;
            exp_bits =  n_randint(state, 100) + 1;
            exp_bits1 = n_randint(state, 100) + 1;
            nmod_mpoly_randtest_bits(a, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);
            nmod_mpoly_add(b, f, a, ctx);
            nmod_mpoly_mul(b, b, f, ctx);
            _test_root(x, a, b, ctx, 1);

            len = n_randint(state, 50);
            len1 = n_randint(state, 50) + 1;
            exp_bits =  n_randint(state, 20) + 1;
            exp_bits1 = n_randint(state, 20) + 1;
            nmod_mpoly_randtest_bits(a, state, len1, 10, ctx);
            nmod_mpoly_randtest_bits(b, state, len, 10, ctx);
            _test_root(x, a, b, ctx, 0);
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(x, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
