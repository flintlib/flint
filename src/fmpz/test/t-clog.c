/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_clog, state)
{
    int i, result;

    /* Check aliasing */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a;
        slong l;

        fmpz_init(a);

        while (fmpz_cmp_ui(a, 2) < 0)
            fmpz_randtest(a, state, 200);

        l = fmpz_clog(a, a);

        result = (l == 1);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("l = %wd\n", l);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
    }

    /* Check correctness */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, x, y;
        slong k;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(x);
        fmpz_init(y);

        while (fmpz_cmp_ui(a, 1) < 0)
            fmpz_randtest(a, state, 200);
        while (fmpz_cmp_ui(b, 2) < 0)
            fmpz_randtest(b, state, 200);

        k = fmpz_clog(a, b);  /* p^{k-1} < a <= p^k*/

        if (k > 0)
            fmpz_pow_ui(x, b, k - 1);
        else
            fmpz_zero(x);
        fmpz_pow_ui(y, b, k);

        result = (fmpz_cmp(x, a) < 0 && fmpz_cmp(a, y) <= 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("x = "), fmpz_print(x), flint_printf("\n");
            flint_printf("y = "), fmpz_print(y), flint_printf("\n");
            flint_printf("k = %wd\n", k);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(x);
        fmpz_clear(y);
    }

    /* Check correctness:  exact powers */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        slong k, l;

        fmpz_init(a);
        fmpz_init(b);

        while (fmpz_cmp_ui(b, 2) < 0)
            fmpz_randtest(b, state, 200);
        l = n_randint(state, 20);
        fmpz_pow_ui(a, b, l);

        k = fmpz_clog(a, b);

        result = (k == l);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("k = %wd\n", k);
            flint_printf("l = %wd\n", l);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
    }

    TEST_FUNCTION_END(state);
}
