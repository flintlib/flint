/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_compose_series, state)
{
    int i, result;

    /* Check aliasing of the first argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_poly_t f, g, h;
        mp_limb_t m;
        slong n;

        m = n_randtest_prime(state, 0);
        nmod_poly_init(f, m);
        nmod_poly_init(g, m);
        nmod_poly_init(h, m);
        nmod_poly_randtest(g, state, n_randint(state, 50));
        nmod_poly_randtest(h, state, n_randint(state, 30));
        nmod_poly_set_coeff_ui(h, 0, 0);
        n = n_randint(state, 40);

        nmod_poly_compose_series(f, g, h, n);
        nmod_poly_compose_series(g, g, h, n);

        result = (nmod_poly_equal(f, g));
        if (!result)
        {
            flint_printf("FAIL (aliasing 1):\n");
            nmod_poly_print(f), flint_printf("\n\n");
            nmod_poly_print(g), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(f);
        nmod_poly_clear(g);
        nmod_poly_clear(h);
    }

    /* Check aliasing of the second argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_poly_t f, g, h;
        mp_limb_t m;
        slong n;

        m = n_randtest_prime(state, 0);
        nmod_poly_init(f, m);
        nmod_poly_init(g, m);
        nmod_poly_init(h, m);
        nmod_poly_randtest(g, state, n_randint(state, 50));
        nmod_poly_randtest(h, state, n_randint(state, 30));
        nmod_poly_set_coeff_ui(h, 0, 0);
        n = n_randint(state, 40);

        nmod_poly_compose_series(f, g, h, n);
        nmod_poly_compose_series(h, g, h, n);

        result = (nmod_poly_equal(f, h));
        if (!result)
        {
            flint_printf("FAIL (aliasing 2):\n");
            nmod_poly_print(f), flint_printf("\n\n");
            nmod_poly_print(h), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(f);
        nmod_poly_clear(g);
        nmod_poly_clear(h);
    }

    /* Compare with compose */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_poly_t f, g, h, s, t;
        mp_limb_t m;
        slong n;

        m = n_randtest_prime(state, 0);
        nmod_poly_init(f, m);
        nmod_poly_init(g, m);
        nmod_poly_init(h, m);
        nmod_poly_init(s, m);
        nmod_poly_init(t, m);
        nmod_poly_randtest(g, state, n_randint(state, 50));
        nmod_poly_randtest(h, state, n_randint(state, 30));
        nmod_poly_set_coeff_ui(h, 0, 0);
        n = n_randint(state, 40);

        nmod_poly_compose(s, g, h);
        nmod_poly_truncate(s, n);
        nmod_poly_compose_series(f, g, h, n);

        result = (nmod_poly_equal(f, s));
        if (!result)
        {
            flint_printf("FAIL (comparison):\n");
            flint_printf("n = %wd\n", n);
            flint_printf("g = "), nmod_poly_print(g), flint_printf("\n\n");
            flint_printf("h = "), nmod_poly_print(h), flint_printf("\n\n");
            flint_printf("f = "), nmod_poly_print(f), flint_printf("\n\n");
            flint_printf("s = "), nmod_poly_print(s), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(f);
        nmod_poly_clear(g);
        nmod_poly_clear(h);
        nmod_poly_clear(s);
        nmod_poly_clear(t);
    }

    TEST_FUNCTION_END(state);
}
