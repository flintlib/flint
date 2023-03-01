/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("get_interval_fmpz_2exp....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 100000 * arb_test_multiplier(); iter++)
    {
        arb_t x;
        arf_t y;
        fmpz_t a, b, exp;

        arb_init(x);
        arf_init(y);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(exp);

        arb_randtest(x, state, 200, 10);

        arb_get_interval_fmpz_2exp(a, b, exp, x);

        arf_set_fmpz_2exp(y, a, exp);
        if (!arb_contains_arf(x, y))
        {
            flint_printf("FAIL (a):\n\n");
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("a = "); fmpz_print(a);
            flint_printf(" exp = "); fmpz_print(exp); flint_printf("\n\n");
            flint_abort();
        }

        arf_set_fmpz_2exp(y, b, exp);
        if (!arb_contains_arf(x, y))
        {
            flint_printf("FAIL (b):\n\n");
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("b = "); fmpz_print(b);
            flint_printf(" exp = "); fmpz_print(exp); flint_printf("\n\n");
            flint_abort();
        }

        if (fmpz_is_even(a) && fmpz_is_even(b) &&
            !(fmpz_is_zero(a) && fmpz_is_zero(b)))
        {
            flint_printf("FAIL (odd):\n\n");
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("a = "); fmpz_print(a); flint_printf("\n\n");
            flint_printf("b = "); fmpz_print(b); flint_printf("\n\n");
            flint_printf(" exp = "); fmpz_print(exp); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        arf_clear(y);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(exp);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
