/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "qqbar.h"

TEST_FUNCTION_START(qqbar_get_quadratic, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        qqbar_t x, y;
        fmpz_t a, b, c, q, d, n;

        qqbar_init(x);
        qqbar_init(y);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(q);
        fmpz_init(d);
        fmpz_init(n);

        switch (n_randint(state, 3))
        {
            case 0:
                qqbar_randtest(x, state, 2, 2 + n_randint(state, 20));
                qqbar_get_quadratic(a, b, c, q, x, 1);
                break;
            case 1:
                qqbar_randtest(x, state, 2, 2 + n_randint(state, 40));
                qqbar_get_quadratic(a, b, c, q, x, 2);
                break;
            default:
                qqbar_randtest(x, state, 2, 2 + n_randint(state, 80));
                qqbar_get_quadratic(a, b, c, q, x, 0);
                break;
        }

        qqbar_set_fmpz(y, c);
        qqbar_sqrt(y, y);
        qqbar_mul_fmpz(y, y, b);
        qqbar_add_fmpz(y, y, a);
        qqbar_div_fmpz(y, y, q);

        if (!qqbar_equal(x, y))
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
            flint_printf("a = "); fmpz_print(a); flint_printf("\n\n");
            flint_printf("b = "); fmpz_print(b); flint_printf("\n\n");
            flint_printf("c = "); fmpz_print(c); flint_printf("\n\n");
            flint_printf("q = "); fmpz_print(q); flint_printf("\n\n");
            flint_abort();
        }

        fmpz_gcd(d, a, b);
        fmpz_gcd(d, d, q);

        if (!fmpz_is_one(d) || fmpz_sgn(q) < 0 || (qqbar_degree(x) == 2 && (fmpz_is_square(c) || fmpz_val2(c) >= 2)))
        {
            flint_printf("FAIL (normalization)\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("a = "); fmpz_print(a); flint_printf("\n\n");
            flint_printf("b = "); fmpz_print(b); flint_printf("\n\n");
            flint_printf("c = "); fmpz_print(c); flint_printf("\n\n");
            flint_printf("q = "); fmpz_print(q); flint_printf("\n\n");
            flint_abort();
        }

        qqbar_clear(x);
        qqbar_clear(y);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(q);
        fmpz_clear(d);
        fmpz_clear(n);
    }

    TEST_FUNCTION_END(state);
}
