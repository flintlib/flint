/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_modular.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("psl2z_inv....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        psl2z_t f, g, h, i;

        psl2z_init(f);
        psl2z_init(g);
        psl2z_init(h);
        psl2z_init(i);

        psl2z_randtest(f, state, n_randint(state, 100));
        psl2z_randtest(g, state, n_randint(state, 100));
        psl2z_randtest(h, state, n_randint(state, 100));

        psl2z_inv(g, f);
        psl2z_mul(h, f, g);
        psl2z_one(i);

        if (!psl2z_equal(h, i) || !psl2z_is_correct(g))
        {
            flint_printf("FAIL\n");
            flint_printf("f = "); psl2z_print(f); flint_printf("\n");
            flint_printf("g = "); psl2z_print(g); flint_printf("\n");
            flint_printf("h = "); psl2z_print(h); flint_printf("\n");
            flint_printf("i = "); psl2z_print(i); flint_printf("\n");
            flint_abort();
        }

        psl2z_inv(f, f);

        if (!psl2z_equal(f, g) || !psl2z_is_correct(f))
        {
            flint_printf("FAIL (aliasing)\n");
            flint_printf("f = "); psl2z_print(f); flint_printf("\n");
            flint_printf("g = "); psl2z_print(g); flint_printf("\n");
            flint_abort();
        }

        psl2z_clear(f);
        psl2z_clear(g);
        psl2z_clear(h);
        psl2z_clear(i);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

