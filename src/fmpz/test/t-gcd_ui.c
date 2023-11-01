/*
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_gcd_ui, state)
{
    int i, result;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t r1, r2, a, fb;
        ulong b;

        fmpz_init(r1);
        fmpz_init(r2);
        fmpz_init(a);
        fmpz_init(fb);

        fmpz_randtest(a, state, 200);
        b = n_randtest(state);
        fmpz_set_ui(fb, b);

        fmpz_gcd(r1, a, fb);

        if (n_randint(state, 2))
        {
            fmpz_gcd_ui(r2, a, b);
        }
        else
        {
            fmpz_set(r2, a);
            fmpz_gcd_ui(r2, r2, b);
        }

        result = (fmpz_cmp(r1, r2) == 0) && _fmpz_is_canonical(r2);

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = ");
            fmpz_print(a);
            flint_printf("\nb = %wu", b);
            flint_printf("\nr1 = ");
            fmpz_print(r1);
            flint_printf("\nr2 = ");
            fmpz_print(r2);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(r1);
        fmpz_clear(r2);
        fmpz_clear(a);
        fmpz_clear(fb);
    }

    TEST_FUNCTION_END(state);
}
