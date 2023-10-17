/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_taylor_shift_convolution, state)
{
    int i;

    /* Check aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t f, g;
        slong n;
        mp_limb_t c, mod;

        n = n_randint(state, 100);
        do {
            mod = n_randtest_prime(state, 0);
        } while (mod <= n);

        nmod_poly_init(f, mod);
        nmod_poly_init(g, mod);

        nmod_poly_randtest(f, state, n);
        c = n_randtest(state) % mod;

        nmod_poly_taylor_shift_convolution(g, f, c);
        nmod_poly_taylor_shift_convolution(f, f, c);

        if (!nmod_poly_equal(g, f))
        {
            flint_printf("FAIL\n");
            nmod_poly_print(f); flint_printf("\n");
            nmod_poly_print(g); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(f);
        nmod_poly_clear(g);
    }

    /* Compare with composition */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t f, g, h1, h2;
        mp_limb_t mod, c;
        slong n;

        n = n_randint(state, 100);
        do {
            mod = n_randtest_prime(state, 0);
        } while (mod <= n);

        nmod_poly_init(f, mod);
        nmod_poly_init(g, mod);
        nmod_poly_init(h1, mod);
        nmod_poly_init(h2, mod);

        nmod_poly_randtest(f, state, n);
        c = n_randtest(state) % mod;

        nmod_poly_set_coeff_ui(g, 1, 1);
        nmod_poly_set_coeff_ui(g, 0, c);

        nmod_poly_taylor_shift_convolution(h1, f, c);
        nmod_poly_compose(h2, f, g);

        if (!nmod_poly_equal(h1, h2))
        {
            flint_printf("FAIL\n");
            nmod_poly_print(f); flint_printf("\n");
            nmod_poly_print(g); flint_printf("\n");
            nmod_poly_print(h1); flint_printf("\n");
            nmod_poly_print(h2); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(f);
        nmod_poly_clear(g);
        nmod_poly_clear(h1);
        nmod_poly_clear(h2);
    }

    TEST_FUNCTION_END(state);
}
