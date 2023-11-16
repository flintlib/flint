/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "arb.h"

TEST_FUNCTION_START(arb_contains_arf, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a;
        arf_t b;
        fmpq_t am, ar, bm, t;
        int c1, c2;

        arb_init(a);
        arf_init(b);

        fmpq_init(am);
        fmpq_init(ar);
        fmpq_init(bm);
        fmpq_init(t);

        arb_randtest(a, state, 1 + n_randint(state, 500), 14);
        arf_randtest(b, state, 1 + n_randint(state, 500), 14);

        arf_get_fmpq(am, arb_midref(a));
        mag_get_fmpq(ar, arb_radref(a));
        arf_get_fmpq(bm, b);

        c1 = arb_contains_arf(a, b);

        fmpq_sub(t, am, ar);
        c2 = fmpq_cmp(t, bm) <= 0;

        fmpq_add(t, am, ar);
        c2 = c2 && (fmpq_cmp(t, bm) >= 0);

        if (c1 != c2)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arf_print(b); flint_printf("\n\n");
            flint_printf("am = "); fmpq_print(am); flint_printf("\n\n");
            flint_printf("ar = "); fmpq_print(ar); flint_printf("\n\n");
            flint_printf("bm = "); fmpq_print(bm); flint_printf("\n\n");
            flint_printf("t = "); fmpq_print(t); flint_printf("\n\n");
            flint_printf("c1 = %d, c2 = %d\n\n", c1, c2);
            flint_abort();
        }

        arb_clear(a);
        arf_clear(b);

        fmpq_clear(am);
        fmpq_clear(ar);
        fmpq_clear(bm);
        fmpq_clear(t);
    }

    TEST_FUNCTION_END(state);
}
