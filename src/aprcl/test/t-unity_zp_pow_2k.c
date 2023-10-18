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

TEST_FUNCTION_START(aprcl_unity_zp_pow_2k, state)
{
    int i, j;

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        ulong p;
        fmpz_t n, pow;
        unity_zp f, g, h;

        p = n_randprime(state, 2 + n_randint(state, 6), 0);

        fmpz_init(n);
        fmpz_randtest_unsigned(n, state, 200);
        while (fmpz_equal_ui(n, 0) != 0)
            fmpz_randtest_unsigned(n, state, 200);
        fmpz_init(pow);
        fmpz_randtest_unsigned(pow, state, 200);
        while (fmpz_equal_ui(pow, 0) != 0)
            fmpz_randtest_unsigned(pow, state, 200);

        unity_zp_init(f, p, 1, n);
        unity_zp_init(g, p, 1, n);
        unity_zp_init(h, p, 1, n);

        for (j = 0; j < 100; j++)
        {
            ulong ind;
            fmpz_t val;

            fmpz_init(val);

            ind = n_randint(state, p);
            fmpz_randtest_unsigned(val, state, 200);
            unity_zp_coeff_set_fmpz(h, ind, val);

            fmpz_clear(val);
        }

        unity_zp_pow_2k_fmpz(f, h, pow);
        unity_zp_pow_fmpz(g, h, pow);

        if (unity_zp_equal(f, g) == 0)
        {
            flint_printf("FAIL\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(n);
        fmpz_clear(pow);
        unity_zp_clear(f);
        unity_zp_clear(g);
        unity_zp_clear(h);
    }

    TEST_FUNCTION_END(state);
}
