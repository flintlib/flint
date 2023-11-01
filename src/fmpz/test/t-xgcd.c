/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_xgcd, state)
{
    int i, result;

    /* Test a f  + b g == d and d >= 0 */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t d, a, b, f, g, t1, t2;
        int aliasing;

        fmpz_init(d);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(f);
        fmpz_init(g);
        fmpz_init(t1);
        fmpz_init(t2);

        fmpz_randtest_unsigned(g, state, 200);
        fmpz_add_ui(g, g, 1);
        fmpz_randm(f, state, g);
        if (n_randint(state, 2)) fmpz_neg(g, g);
        if (n_randint(state, 2)) fmpz_neg(f, f);

        aliasing = n_randint(state, 5);

        if (aliasing == 0)
        {
            fmpz_xgcd(d, a, b, f, g);
        }
        else if (aliasing == 1)
        {
            /* Test aliasing of d and f, a and g */
            fmpz_set(d, f);
            fmpz_set(a, g);
            fmpz_xgcd(d, a, b, d, a);
        }
        else if (aliasing == 2)
        {
            /* Test aliasing of a and f, d and g */
            fmpz_set(a, f);
            fmpz_set(d, g);
            fmpz_xgcd(d, a, b, a, d);
        }
        else if (aliasing == 3)
        {
            /* Test aliasing of d and f, b and g */
            fmpz_set(d, f);
            fmpz_set(b, g);
            fmpz_xgcd(d, a, b, d, b);
        }
        else
        {
            /* Test aliasing of b and f, d and g */
            fmpz_set(b, f);
            fmpz_set(d, g);
            fmpz_xgcd(d, a, b, b, d);
        }

        fmpz_mul(t1, a, f);
        fmpz_mul(t2, b, g);
        fmpz_add(t1, t1, t2);

        result = fmpz_equal(t1, d) && fmpz_sgn(d) >= 0 &&
            _fmpz_is_canonical(d) && _fmpz_is_canonical(a) && _fmpz_is_canonical(b);
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("aliasing = %d\n", aliasing);
            flint_printf("d = "), fmpz_print(d), flint_printf("\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("f = "), fmpz_print(f), flint_printf("\n");
            flint_printf("g = "), fmpz_print(g), flint_printf("\n");
            flint_printf("t1 = "), fmpz_print(t1), flint_printf("\n");
            flint_printf("t2 = "), fmpz_print(t2), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(d);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(f);
        fmpz_clear(g);
        fmpz_clear(t1);
        fmpz_clear(t2);
    }

    TEST_FUNCTION_END(state);
}
