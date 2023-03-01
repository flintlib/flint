/*
    Copyright (C) 2014 Fredrik Johansson

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

    flint_printf("get_mag_lower_nonnegative....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 100000 * arb_test_multiplier(); iter++)
    {
        arb_t a;
        mag_t m;
        int result;
        arf_struct t[3];
        arf_t s;

        arb_init(a);
        mag_init(m);
        arf_init(s);

        arb_randtest_special(a, state, 200, 1 + n_randint(state, 100));
        arb_get_mag_lower_nonnegative(m, a);
        MAG_CHECK_BITS(m)

        if (arb_contains_nonpositive(a))
        {
            result = mag_is_zero(m);
        }
        else
        {
            arf_init_set_shallow(t + 0, arb_midref(a));
            arf_init_neg_mag_shallow(t + 1, arb_radref(a));
            arf_init_neg_mag_shallow(t + 2, m);
            arf_sum(s, t, 3, 16, ARF_RND_DOWN);
            result = (arf_sgn(s) >= 0);
        }

        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("m = "); mag_print(m); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        mag_clear(m);
        arf_clear(s);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

