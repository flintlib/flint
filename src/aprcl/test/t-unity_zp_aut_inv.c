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

TEST_FUNCTION_START(aprcl_unity_zp_aut_inv, state)
{
    ulong i, j;

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        ulong ind, q, p, k, x;
        fmpz_t n;
        unity_zp f, g, h;
        n_factor_t q_factors;

        n_factor_init(&q_factors);

        q = n_randprime(state, 2 + n_randint(state, 6), 0);
        while (q < 3)
            q = n_randprime(state, 2 + n_randint(state, 6), 0);

        n_factor(&q_factors, q - 1, 1);
        ind = n_randint(state, q_factors.num);
        p = q_factors.p[ind];
        k = q_factors.exp[ind];

        x = n_randint(state, n_pow(p, k));
        while (n_gcd(p, x) != 1 || x == 0)
            x = n_randint(state, n_pow(p, k));

        fmpz_init(n);
        fmpz_randtest_unsigned(n, state, 200);
        while (fmpz_equal_ui(n, 0) != 0)
            fmpz_randtest_unsigned(n, state, 200);

        unity_zp_init(f, p, k, n);
        unity_zp_init(g, p, k, n);
        unity_zp_init(h, p, k, n);

        for (j = 0; j < 100; j++)
        {
            ulong ind;
            fmpz_t val;

            fmpz_init(val);

            ind = n_randint(state, n_pow(p, k));
            fmpz_randtest_unsigned(val, state, 200);
            unity_zp_coeff_set_fmpz(g, ind, val);

            fmpz_clear(val);
        }

        /* reduce random element h */
        unity_zp_reduce_cyclotomic(h, g);
        /* \sigma_x(f) == h now */
        unity_zp_aut_inv(f, h, x);
        /* g = \sigma_x(f) */
        unity_zp_aut(g, f, x);

        if (unity_zp_equal(h, g) == 0)
        {
            flint_printf("FAIL\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(n);
        unity_zp_clear(f);
        unity_zp_clear(g);
        unity_zp_clear(h);
    }

    TEST_FUNCTION_END(state);
}
