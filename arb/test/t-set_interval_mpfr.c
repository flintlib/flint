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

    flint_printf("set_interval_mpfr....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 100000 * arb_test_multiplier(); iter++)
    {
        arb_t x;
        arf_t a, b;
        mpfr_t aa, bb;

        arb_init(x);
        arf_init(a);
        arf_init(b);
        mpfr_init2(aa, 200);
        mpfr_init2(bb, 200);

        arf_randtest_special(a, state, 200, 10);
        arf_randtest_special(b, state, 200, 10);
        if (arf_cmp(a, b) > 0)
            arf_swap(a, b);

        arf_get_mpfr(aa, a, MPFR_RNDD);
        arf_get_mpfr(bb, b, MPFR_RNDU);

        arb_set_interval_mpfr(x, aa, bb, 2 + n_randint(state, 200));

        if (!arb_contains_arf(x, a) || !arb_contains_arf(x, b))
        {
            flint_printf("FAIL:\n\n");
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("a = "); arf_print(a); flint_printf("\n\n");
            flint_printf("b = "); arf_print(b); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        arf_clear(a);
        arf_clear(b);
        mpfr_clear(aa);
        mpfr_clear(bb);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

