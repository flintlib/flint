/*
    Copyright (C) 2017 Fredrik Johansson

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

    flint_printf("sinc_pi....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
        acb_t x, y, z;
        slong prec;

        acb_init(x);
        acb_init(y);
        acb_init(z);

        prec = 2 + n_randint(state, 200);

        acb_randtest(x, state, 1 + n_randint(state, 200), 1 + n_randint(state, 100));
        acb_randtest(y, state, 1 + n_randint(state, 200), 1 + n_randint(state, 100));

        if (n_randint(state, 2))
        {
            acb_set(y, x);
            acb_sinc_pi(y, y, prec);
        }
        else
            acb_sinc_pi(y, x, prec);

        acb_const_pi(z, prec);
        acb_mul(z, z, x, prec);
        acb_sinc(z, z, prec);

        if (!acb_overlaps(y, z))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("x = "); acb_printd(x, 30); flint_printf("\n\n");
            flint_printf("y = "); acb_printd(y, 30); flint_printf("\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(x);
        acb_clear(y);
        acb_clear(z);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
