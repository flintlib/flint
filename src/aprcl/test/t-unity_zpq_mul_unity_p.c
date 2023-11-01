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

TEST_FUNCTION_START(aprcl_unity_zpq_mul_unity_p, state)
{
    int i, j;

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        ulong p, q;
        fmpz_t n;
        unity_zpq f, g;

        fmpz_init(n);

        p = n_randprime(state, 2 + n_randint(state, 6), 0);
        q = n_randprime(state, 2 + n_randint(state, 6), 0);

        fmpz_randtest_unsigned(n, state, 200);
        while (fmpz_equal_ui(n, 0) != 0)
            fmpz_randtest_unsigned(n, state, 200);

        unity_zpq_init(f, q, p, n);
        unity_zpq_init(g, q, p, n);

        for (j = 0; j < 100; j++)
        {
            ulong x, y;
            fmpz_t val;

            fmpz_init(val);

            x = n_randint(state, p);
            y = n_randint(state, q);

            fmpz_randtest_unsigned(val, state, 200);

            unity_zpq_coeff_set_fmpz(g, y, x, val);

            fmpz_clear(val);
        }

        for (j = 0; j < p; j++)
        {
            unity_zpq h1, h2;
            unity_zpq_init(h1, q, p, n);
            unity_zpq_init(h2, q, p, n);

            unity_zpq_coeff_set_ui(h1, 0, j, 1);
            unity_zpq_mul(h2, g, h1);
            unity_zpq_mul_unity_p_pow(f, g, j);

            if (unity_zpq_equal(f, h2) == 0)
            {
                flint_printf("FAIL\n");
                fflush(stdout);
                flint_abort();
            }

            unity_zpq_clear(h1);
            unity_zpq_clear(h2);
        }

        fmpz_clear(n);
        unity_zpq_clear(f);
        unity_zpq_clear(g);
    }

    TEST_FUNCTION_END(state);
}
