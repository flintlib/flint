/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);

    flint_printf("taylor_shift....");
    fflush(stdout);

    

    /* Check aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t f, g;
        mp_limb_t c, mod;

        mod = n_randtest_prime(state, 0);

        nmod_poly_init(f, mod);
        nmod_poly_init(g, mod);

        nmod_poly_randtest(f, state, 1 + n_randint(state, 50));
        c = n_randtest(state) % mod;

        nmod_poly_taylor_shift(g, f, c);
        nmod_poly_taylor_shift(f, f, c);

        if (!nmod_poly_equal(g, f))
        {
            flint_printf("FAIL\n");
            nmod_poly_print(f); flint_printf("\n");
            nmod_poly_print(g); flint_printf("\n");
            abort();
        }

        nmod_poly_clear(f);
        nmod_poly_clear(g);
    }

    /* Compare with composition */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t f, g, h1, h2;
        mp_limb_t mod, c;

        mod = n_randtest_prime(state, 0);

        nmod_poly_init(f, mod);
        nmod_poly_init(g, mod);
        nmod_poly_init(h1, mod);
        nmod_poly_init(h2, mod);

        nmod_poly_randtest(f, state, 1 + n_randint(state, 50));
        c = n_randtest(state) % mod;

        nmod_poly_set_coeff_ui(g, 1, 1);
        nmod_poly_set_coeff_ui(g, 0, c);

        nmod_poly_taylor_shift(h1, f, c);
        nmod_poly_compose(h2, f, g);

        if (!nmod_poly_equal(h1, h2))
        {
            flint_printf("FAIL\n");
            nmod_poly_print(f); flint_printf("\n");
            nmod_poly_print(g); flint_printf("\n");
            nmod_poly_print(h1); flint_printf("\n");
            nmod_poly_print(h2); flint_printf("\n");
            abort();
        }

        nmod_poly_clear(f);
        nmod_poly_clear(g);
        nmod_poly_clear(h1);
        nmod_poly_clear(h2);
    }

    /* Check some large cases */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_poly_t f, g, h1, h2;
        mp_limb_t mod, c;

        mod = n_randtest_prime(state, 0);

        nmod_poly_init(f, mod);
        nmod_poly_init(g, mod);
        nmod_poly_init(h1, mod);
        nmod_poly_init(h2, mod);

        nmod_poly_randtest(f, state, 1 + n_randint(state, 2000));
        c = n_randtest(state) % mod;

        nmod_poly_set_coeff_ui(g, 1, 1);
        nmod_poly_set_coeff_ui(g, 0, c);

        nmod_poly_taylor_shift(h1, f, c);
        nmod_poly_compose(h2, f, g);

        if (!nmod_poly_equal(h1, h2))
        {
            flint_printf("FAIL\n");
            nmod_poly_print(f); flint_printf("\n");
            nmod_poly_print(g); flint_printf("\n");
            nmod_poly_print(h1); flint_printf("\n");
            nmod_poly_print(h2); flint_printf("\n");
            abort();
        }

        nmod_poly_clear(f);
        nmod_poly_clear(g);
        nmod_poly_clear(h1);
        nmod_poly_clear(h2);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
