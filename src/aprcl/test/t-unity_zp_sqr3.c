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

TEST_FUNCTION_START(aprcl_unity_zp_sqr3, state)
{
    int i, j;
    fmpz_t * t;

    t = (fmpz_t*) flint_malloc(sizeof(fmpz_t) * (SQUARING_SPACE));
    for (i = 0; i < SQUARING_SPACE; i++)
        fmpz_init(t[i]);

    /* test squaring in Z[\zeta_3] */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        ulong p;
        fmpz_t n;
        unity_zp f, g, temp;

        p = 3;

        fmpz_init(n);
        fmpz_randtest_unsigned(n, state, 200);
        while (fmpz_equal_ui(n, 0) != 0)
            fmpz_randtest_unsigned(n, state, 200);

        unity_zp_init(f, p, 1, n);
        unity_zp_init(g, p, 1, n);
        unity_zp_init(temp, p, 1, n);

        for (j = 0; j < 100; j++)
        {
            ulong ind;
            fmpz_t val;

            fmpz_init(val);

            ind = n_randint(state, p);

            fmpz_randtest_unsigned(val, state, 200);

            unity_zp_coeff_set_fmpz(temp, ind, val);

            fmpz_clear(val);
        }

        _unity_zp_reduce_cyclotomic(temp);
        unity_zp_sqr3(f, temp, t);
        unity_zp_sqr(g, temp);

        if (unity_zp_equal(f, g) == 0)
        {
            flint_printf("FAIL\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(n);
        unity_zp_clear(f);
        unity_zp_clear(g);
        unity_zp_clear(temp);
    }

    /* test squaring in Z[\zeta_9] */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        ulong p, k;
        fmpz_t n;
        unity_zp f, g, temp;

        p = 3;
        k = 2;

        fmpz_init(n);
        fmpz_randtest_unsigned(n, state, 200);
        while (fmpz_equal_ui(n, 0) != 0)
            fmpz_randtest_unsigned(n, state, 200);

        unity_zp_init(f, p, k, n);
        unity_zp_init(g, p, k, n);
        unity_zp_init(temp, p, k, n);

        for (j = 0; j < 100; j++)
        {
            ulong ind;
            fmpz_t val;

            fmpz_init(val);

            ind = n_randint(state, n_pow(p, k));

            fmpz_randtest_unsigned(val, state, 200);

            unity_zp_coeff_set_fmpz(temp, ind, val);

            fmpz_clear(val);
        }

        _unity_zp_reduce_cyclotomic(temp);
        unity_zp_sqr9(f, temp, t);
        unity_zp_sqr(g, temp);

        if (unity_zp_equal(f, g) == 0)
        {
            flint_printf("FAIL\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(n);
        unity_zp_clear(f);
        unity_zp_clear(g);
        unity_zp_clear(temp);
    }

    for (i = 0; i < SQUARING_SPACE; i++)
        fmpz_clear(t[i]);
    flint_free(t);

    TEST_FUNCTION_END(state);
}
