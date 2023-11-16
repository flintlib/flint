/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_extras.h"

TEST_FUNCTION_START(fmpz_add_inline, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpz_t a, b, c, d;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(d);

        fmpz_randtest(a, state, 1 + n_randint(state, 200));
        fmpz_randtest(b, state, 1 + n_randint(state, 200));
        fmpz_randtest(c, state, 1 + n_randint(state, 200));
        fmpz_randtest(d, state, 1 + n_randint(state, 200));

        fmpz_add(c, a, b);
        fmpz_add_inline(d, a, b);
        if (!fmpz_equal(c, d))
        {
            flint_printf("FAIL\n");
            fmpz_print(a); flint_printf("\n\n");
            fmpz_print(b); flint_printf("\n\n");
            fmpz_print(c); flint_printf("\n\n");
            fmpz_print(d); flint_printf("\n\n");
            flint_abort();
        }

        fmpz_add_inline(a, a, b);
        if (!fmpz_equal(c, a))
        {
            flint_printf("FAIL (aliasing 1)\n");
            fmpz_print(a); flint_printf("\n\n");
            fmpz_print(b); flint_printf("\n\n");
            fmpz_print(c); flint_printf("\n\n");
            fmpz_print(d); flint_printf("\n\n");
            flint_abort();
        }
        fmpz_randtest(a, state, 1 + n_randint(state, 200));

        fmpz_add(c, a, b);
        fmpz_add_inline(b, a, b);
        if (!fmpz_equal(c, b))
        {
            flint_printf("FAIL (aliasing 2)\n");
            fmpz_print(a); flint_printf("\n\n");
            fmpz_print(b); flint_printf("\n\n");
            fmpz_print(c); flint_printf("\n\n");
            fmpz_print(d); flint_printf("\n\n");
            flint_abort();
        }

        fmpz_add(d, a, a);
        fmpz_add_inline(c, a, a);
        if (!fmpz_equal(c, d))
        {
            flint_printf("FAIL (aliasing 3)\n");
            fmpz_print(a); flint_printf("\n\n");
            fmpz_print(b); flint_printf("\n\n");
            fmpz_print(c); flint_printf("\n\n");
            fmpz_print(d); flint_printf("\n\n");
            flint_abort();
        }

        fmpz_add(d, a, a);
        fmpz_add_inline(a, a, a);
        if (!fmpz_equal(d, a))
        {
            flint_printf("FAIL (aliasing 4)\n");
            fmpz_print(a); flint_printf("\n\n");
            fmpz_print(b); flint_printf("\n\n");
            fmpz_print(c); flint_printf("\n\n");
            fmpz_print(d); flint_printf("\n\n");
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(d);
    }

    TEST_FUNCTION_END(state);
}
