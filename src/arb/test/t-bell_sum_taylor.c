/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_bell_sum_taylor, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t s1, s2;
        fmpz_t a, b, n;
        slong prec;

        arb_init(s1);
        arb_init(s2);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(n);

        prec = 2 + n_randint(state, 300);
        fmpz_randtest_unsigned(n, state, 1 + n_randint(state, 100));
        fmpz_randtest_unsigned(a, state, 1 + n_randint(state, 100));
        fmpz_add_ui(b, a, n_randint(state, 100));

        arb_bell_sum_bsplit(s1, n, a, b, NULL, prec);
        arb_bell_sum_taylor(s2, n, a, b, NULL, prec);

        if (!arb_overlaps(s1, s2) || (arb_rel_accuracy_bits(s1) < prec - 4)
            || (arb_rel_accuracy_bits(s2) < prec - 4))
        {
            flint_printf("FAIL: overlap or accuracy\n\n");
            flint_printf("prec = %wd\n\n", prec);
            flint_printf("n = "); fmpz_print(n); flint_printf("\n\n");
            flint_printf("a = "); fmpz_print(a); flint_printf("\n\n");
            flint_printf("b = "); fmpz_print(b); flint_printf("\n\n");
            flint_printf("s1 = "); arb_printn(s1, 100, 0); flint_printf("\n\n");
            flint_printf("s2 = "); arb_printn(s2, 100, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(s1);
        arb_clear(s2);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(n);
    }

    TEST_FUNCTION_END(state);
}
