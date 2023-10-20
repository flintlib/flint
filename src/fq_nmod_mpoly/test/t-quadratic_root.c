/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_nmod_mpoly.h"

void _test_root(
    fq_nmod_mpoly_t x,
    const fq_nmod_mpoly_t a,
    const fq_nmod_mpoly_t b,
    const fq_nmod_mpoly_ctx_t ctx,
    int sol_exists)
{
    int success, success2;
    fq_nmod_mpoly_t s, t;

    fq_nmod_mpoly_init(s, ctx);
    fq_nmod_mpoly_init(t, ctx);

    success = fq_nmod_mpoly_quadratic_root(x, a, b, ctx);

    if (sol_exists && !success)
    {
        flint_printf("FAIL: solution exists but root failed\n");
        fflush(stdout);
        flint_abort();
    }

    if (success)
    {
        fq_nmod_mpoly_add(t, x, a, ctx);
        fq_nmod_mpoly_mul(s, t, x, ctx);
        if (!fq_nmod_mpoly_equal(s, b, ctx))
        {
            flint_printf("FAIL: reported solution is not a solution\n");
            fflush(stdout);
            flint_abort();
        }
    }

    fq_nmod_mpoly_set(t, a, ctx);
    success2 = fq_nmod_mpoly_quadratic_root(t, t, b, ctx);
    if (success != success2 || (success && !fq_nmod_mpoly_equal(x, t, ctx)))
    {
        flint_printf("FAIL: Check aliasing first argument\n");
        fflush(stdout);
        flint_abort();
    }

    fq_nmod_mpoly_set(t, b, ctx);
    success2 = fq_nmod_mpoly_quadratic_root(t, a, t, ctx);
    if (success != success2 || (success && !fq_nmod_mpoly_equal(x, t, ctx)))
    {
        flint_printf("FAIL: Check aliasing second argument\n");
        fflush(stdout);
        flint_abort();
    }

    fq_nmod_mpoly_clear(s, ctx);
    fq_nmod_mpoly_clear(t, ctx);
}

TEST_FUNCTION_START(fq_nmod_mpoly_quadratic_root, state)
{
    slong i, j, tmul = 20;

    /* Check sqrt(f^2) = +-f */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, a, b, x;
        slong len, len1;
        flint_bitcnt_t exp_bits, exp_bits1;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 9, (i % 2) ? FLINT_BITS : 1, 5);

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(b, ctx);
        fq_nmod_mpoly_init(x, ctx);

        for (j = 0; j < 5; j++)
        {
            len = n_randint(state, 100);
            len1 = n_randint(state, 100) + 1;
            exp_bits =  n_randint(state, 100) + 1;
            exp_bits1 = n_randint(state, 100) + 1;
            fq_nmod_mpoly_randtest_bits(a, state, len1, exp_bits1, ctx);
            fq_nmod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);
            fq_nmod_mpoly_add(b, f, a, ctx);
            fq_nmod_mpoly_mul(b, b, f, ctx);
            _test_root(x, a, b, ctx, 1);

            len = n_randint(state, 50);
            len1 = n_randint(state, 50) + 1;
            exp_bits =  n_randint(state, 20) + 1;
            exp_bits1 = n_randint(state, 20) + 1;
            fq_nmod_mpoly_randtest_bits(a, state, len1, 10, ctx);
            fq_nmod_mpoly_randtest_bits(b, state, len, 10, ctx);
            _test_root(x, a, b, ctx, 0);
        }

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_clear(b, ctx);
        fq_nmod_mpoly_clear(x, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
