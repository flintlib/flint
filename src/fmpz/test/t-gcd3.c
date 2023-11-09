/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_gcd3, state)
{
    slong i;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c, d, e;
        int aliasing;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(d);
        fmpz_init(e);

        fmpz_randtest(a, state, 200);
        fmpz_randtest(b, state, 200);
        fmpz_randtest(c, state, 200);
        fmpz_randtest(d, state, 200);
        fmpz_randtest(e, state, 200);

        fmpz_mul(a, a, d);
        fmpz_mul(b, b, d);
        fmpz_mul(c, c, d);

        fmpz_gcd(d, a, b);
        fmpz_gcd(d, d, c);

        aliasing = n_randint(state, 4);

        if (aliasing == 0)
        {
            fmpz_gcd3(e, a, b, c);
        }
        else if (aliasing == 1)
        {
            fmpz_set(e, a);
            fmpz_gcd3(e, e, b, c);
        }
        else if (aliasing == 2)
        {
            fmpz_set(e, b);
            fmpz_gcd3(e, a, e, c);
        }
        else
        {
            fmpz_set(e, c);
            fmpz_gcd3(e, a, b, e);
        }

        if (!fmpz_equal(d, e))
        {
            flint_printf("FAIL:\n");
            flint_printf("aliasing = %d\n", aliasing);
            fmpz_print(a); printf("\n");
            fmpz_print(b); printf("\n");
            fmpz_print(c); printf("\n");
            fmpz_print(d); printf("\n");
            fmpz_print(e); printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(d);
        fmpz_clear(e);
    }

    TEST_FUNCTION_END(state);
}
