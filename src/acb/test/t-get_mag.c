/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("get_mag....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        acb_t a;
        arf_t m2, x, y, s;
        mag_t m;

        acb_init(a);
        mag_init(m);
        arf_init(m2);
        arf_init(x);
        arf_init(y);
        arf_init(s);

        acb_randtest_special(a, state, 200, 10);
        acb_get_mag(m, a);
        MAG_CHECK_BITS(m)

        /* check m^2 >= x^2 + y^2 */
        arf_set_mag(m2, m);
        arf_mul(m2, m2, m2, ARF_PREC_EXACT, ARF_RND_DOWN);

        arb_get_abs_ubound_arf(x, acb_realref(a), ARF_PREC_EXACT);
        arb_get_abs_ubound_arf(y, acb_imagref(a), ARF_PREC_EXACT);
        arf_sosq(s, x, y, ARF_PREC_EXACT, ARF_RND_DOWN);

        if (arf_cmp(m2, s) < 0)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("m = "); mag_print(m); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(a);
        mag_clear(m);
        arf_clear(m2);
        arf_clear(x);
        arf_clear(y);
        arf_clear(s);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

