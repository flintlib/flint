/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mod_mpoly.h"

void _test_root(
    fmpz_mod_mpoly_t x,
    const fmpz_mod_mpoly_t a,
    const fmpz_mod_mpoly_t b,
    const fmpz_mod_mpoly_ctx_t ctx,
    int sol_exists)
{
    int success, success2;
    fmpz_mod_mpoly_t s, t;

    fmpz_mod_mpoly_init(s, ctx);
    fmpz_mod_mpoly_init(t, ctx);

    success = fmpz_mod_mpoly_quadratic_root(x, a, b, ctx);

    if (sol_exists && !success)
    {
        flint_printf("FAIL: solution exists but root failed\n");
        fflush(stdout);
        flint_abort();
    }

    if (success)
    {
        fmpz_mod_mpoly_add(t, x, a, ctx);
        fmpz_mod_mpoly_mul(s, t, x, ctx);
        if (!fmpz_mod_mpoly_equal(s, b, ctx))
        {

            flint_printf("FAIL: reported solution is not a solution\n");
            fflush(stdout);
            flint_abort();
        }
    }

    fmpz_mod_mpoly_set(t, a, ctx);
    success2 = fmpz_mod_mpoly_quadratic_root(t, t, b, ctx);
    if (success != success2 || (success && !fmpz_mod_mpoly_equal(x, t, ctx)))
    {
        flint_printf("FAIL: Check aliasing first argument\n");
        fflush(stdout);
        flint_abort();
    }

    fmpz_mod_mpoly_set(t, b, ctx);
    success2 = fmpz_mod_mpoly_quadratic_root(t, a, t, ctx);
    if (success != success2 || (success && !fmpz_mod_mpoly_equal(x, t, ctx)))
    {
        flint_printf("FAIL: Check aliasing second argument\n");
        fflush(stdout);
        flint_abort();
    }

    fmpz_mod_mpoly_clear(s, ctx);
    fmpz_mod_mpoly_clear(t, ctx);
}

TEST_FUNCTION_START(fmpz_mod_mpoly_quadratic_root, state)
{
    slong i, j, tmul = 20;

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f, a, b, x;
        slong len, len1;
        flint_bitcnt_t exp_bits, exp_bits1;

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 10, 100);

        fmpz_mod_mpoly_init(f, ctx);
        fmpz_mod_mpoly_init(a, ctx);
        fmpz_mod_mpoly_init(b, ctx);
        fmpz_mod_mpoly_init(x, ctx);

        for (j = 0; j < 5; j++)
        {
            len = n_randint(state, 100);
            len1 = n_randint(state, 100) + 1;
            exp_bits =  n_randint(state, 100) + 1;
            exp_bits1 = n_randint(state, 100) + 1;
            fmpz_mod_mpoly_randtest_bits(a, state, len1, exp_bits1, ctx);
            fmpz_mod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);
            fmpz_mod_mpoly_add(b, f, a, ctx);
            fmpz_mod_mpoly_mul(b, b, f, ctx);
            _test_root(x, a, b, ctx, 1);

            len = n_randint(state, 50);
            len1 = n_randint(state, 50) + 1;
            exp_bits =  n_randint(state, 20) + 1;
            exp_bits1 = n_randint(state, 20) + 1;
            fmpz_mod_mpoly_randtest_bits(a, state, len1, 10, ctx);
            fmpz_mod_mpoly_randtest_bits(b, state, len, 10, ctx);
            _test_root(x, a, b, ctx, 0);
        }

        fmpz_mod_mpoly_clear(f, ctx);
        fmpz_mod_mpoly_clear(a, ctx);
        fmpz_mod_mpoly_clear(b, ctx);
        fmpz_mod_mpoly_clear(x, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
