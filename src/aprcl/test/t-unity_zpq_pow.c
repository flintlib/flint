/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "aprcl.h"

TEST_FUNCTION_START(aprcl_unity_zpq_pow, state)
{
    int i, j;

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        ulong p, q, pow;
        fmpz_t n;
        unity_zpq f, g, h1, h2;

        p = n_randprime(state, 2 + n_randint(state, 5), 0);
        q = n_randprime(state, 2 + n_randint(state, 5), 0);

        pow =  n_randint(state, 32);

        fmpz_init(n);
        fmpz_randtest_unsigned(n, state, 100);
        while (fmpz_equal_ui(n, 0) != 0)
            fmpz_randtest_unsigned(n, state, 100);

        unity_zpq_init(f, q, p, n);
        unity_zpq_init(g, q, p, n);
        unity_zpq_init(h1, q, p, n);
        unity_zpq_init(h2, q, p, n);

        for (j = 0; j < 50; j++)
        {
            ulong x, y;
            fmpz_t val;

            fmpz_init(val);

            x = n_randint(state, p);
            y = n_randint(state, q);

            fmpz_randtest_not_zero(val, state, 100);

            unity_zpq_coeff_set_fmpz(h1, y, x, val);

            fmpz_clear(val);
        }

        unity_zpq_copy(h2, h1);

        unity_zpq_pow_ui(f, h2, pow);

        if (pow == 0)
        {
            unity_zpq_coeff_set_ui(g, 0, 0, 1);
        }
        else
        {
            for (j = 0; j < pow; j++)
            {
                unity_zpq_mul(g, h1, h2);
                unity_zpq_swap(h2, g);
            }
        }

        if (unity_zpq_equal(f, g) == 0)
        {
            flint_printf("FAIL\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(n);
        unity_zpq_clear(f);
        unity_zpq_clear(g);
        unity_zpq_clear(h1);
        unity_zpq_clear(h2);
    }

    TEST_FUNCTION_END(state);
}
