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

TEST_FUNCTION_START(aprcl_unity_zp_is_unity, state)
{
    int i;

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        ulong p, exp;
        fmpz_t n;
        unity_zp f;
        ulong ind;

        p = n_randprime(state, 2 + n_randint(state, 4), 0);
        exp =  n_randint(state, 5);
        while (exp == 0)
            exp = n_randint(state, 5);

        fmpz_init(n);
        fmpz_randtest_unsigned(n, state, 200);
        while (fmpz_equal_ui(n, 0) != 0)
            fmpz_randtest_unsigned(n, state, 200);

        unity_zp_init(f, p, 1, n);

        ind = n_randint(state, n_pow(p, exp));

        unity_zp_coeff_set_ui(f, ind, 1);

        if (unity_zp_is_unity(f) < 0)
        {
            flint_printf("FAIL\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(n);
        unity_zp_clear(f);
    }

    TEST_FUNCTION_END(state);
}
