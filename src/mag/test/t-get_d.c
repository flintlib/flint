/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arf.h"
#include "double_extras.h"
#include "mag.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("get_d....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000 * arb_test_multiplier(); iter++)
    {
        mag_t x;
        arf_t a, b;
        double y;

        mag_init(x);
        arf_init(a);
        arf_init(b);

        mag_randtest_special(x, state, 1 + n_randint(state, 100));

        y = mag_get_d(x);
        arf_set_d(a, y);

        arf_set_mag(b, x);

        if (mag_cmp_2exp_si(x, 1000) < 0 && mag_cmp_2exp_si(x, -1000) > 0)
        {
            if (!arf_equal(a, b))
            {
                flint_printf("FAIL (equality)\n\n");
                flint_printf("x = "); mag_print(x); flint_printf("\n\n");
                flint_printf("a = "); arf_print(a); flint_printf("\n\n");
                flint_printf("b = "); arf_print(b); flint_printf("\n\n");
                flint_abort();
            }
        }
        else
        {
            if (arf_cmp(a, b) < 0)
            {
                flint_printf("FAIL (bound)\n\n");
                flint_printf("x = "); mag_print(x); flint_printf("\n\n");
                flint_printf("a = "); arf_print(a); flint_printf("\n\n");
                flint_printf("b = "); arf_print(b); flint_printf("\n\n");
                flint_abort();
            }
        }

        arf_clear(a);
        arf_clear(b);
        mag_clear(x);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}

