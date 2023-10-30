/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(n_CRT, state)
{
    slong i, j;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t m1, m2, m1m2, a1, a2, r;
        flint_bitcnt_t b1 = n_randint(state, FLINT_BITS) + 1;
        flint_bitcnt_t b2 = n_randint(state, FLINT_BITS + 1 - b1) + 1;

        fmpz_init(m1);
        fmpz_init(m2);
        fmpz_init(m1m2);
        fmpz_init(a1);
        fmpz_init(a2);
        fmpz_init(r);

        fmpz_randtest_not_zero(m1, state, b1);
        fmpz_randtest_not_zero(m2, state, b2);
        fmpz_abs(m1, m1);
        fmpz_abs(m2, m2);
        fmpz_mul(m1m2, m1, m2);

        fmpz_gcd(r, m1, m2);

        if (fmpz_is_one(r) &&
            fmpz_abs_fits_ui(m1) &&
            fmpz_abs_fits_ui(m2) &&
            fmpz_abs_fits_ui(m1m2))
        {
            for (j = 0; j < 50; j++)
            {
                fmpz_randm(a1, state, m1);
                fmpz_randm(a2, state, m2);
                fmpz_CRT(r, a1, m1, a2, m2, 0);
                if (fmpz_get_ui(r) != n_CRT(fmpz_get_ui(a1), fmpz_get_ui(m1),
                                            fmpz_get_ui(a2), fmpz_get_ui(m2)))
                {
                    flint_printf("FAIL\n:");
                    flint_printf("a1: "); fmpz_print(a1); flint_printf("\n");
                    flint_printf("m1: "); fmpz_print(m1); flint_printf("\n");
                    flint_printf("a2: "); fmpz_print(a2); flint_printf("\n");
                    flint_printf("m2: "); fmpz_print(m2); flint_printf("\n");
                }
            }
        }

        fmpz_clear(m1);
        fmpz_clear(m2);
        fmpz_clear(m1m2);
        fmpz_clear(a1);
        fmpz_clear(a2);
        fmpz_clear(r);
    }

    TEST_FUNCTION_END(state);
}
