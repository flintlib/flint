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
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_revert_series, state)
{
    int i, result;

    /* Check aliasing */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_poly_t f, g;
        mp_limb_t m;
        slong n;

        m = n_randtest_prime(state, 0);
        nmod_poly_init(f, m);
        nmod_poly_init(g, m);
        do {
            nmod_poly_randtest(g, state, n_randint(state, 100));
        } while (nmod_poly_get_coeff_ui(g, 1) == 0);
        nmod_poly_set_coeff_ui(g, 0, 0);
        do {
            n = n_randint(state, 100);
        } while (n >= m);

        nmod_poly_revert_series(f, g, n);
        nmod_poly_revert_series(g, g, n);

        result = (nmod_poly_equal(f, g));
        if (!result)
        {
            flint_printf("FAIL (aliasing):\n");
            nmod_poly_print(f), flint_printf("\n\n");
            nmod_poly_print(g), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(f);
        nmod_poly_clear(g);
    }

    /* Check f(f^(-1)) = id */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t f, g, h;
        mp_limb_t m;
        slong n;

        m = n_randtest_prime(state, 0);
        nmod_poly_init(f, m);
        nmod_poly_init(g, m);
        nmod_poly_init(h, m);
        do {
            nmod_poly_randtest(g, state, n_randint(state, 100));
        } while (nmod_poly_get_coeff_ui(g, 1) == 0);
        nmod_poly_set_coeff_ui(g, 0, 0);
        do {
            n = n_randint(state, 100);
        } while (n >= m);

        nmod_poly_revert_series(f, g, n);
        nmod_poly_compose_series(h, g, f, n);

        result = ((n <= 1 && nmod_poly_is_zero(h)) ||
            (h->length == 2 && h->coeffs[0] == 0 && h->coeffs[1] == 1));
        if (!result)
        {
            flint_printf("FAIL (comparison):\n");
            nmod_poly_print(g), flint_printf("\n\n");
            nmod_poly_print(f), flint_printf("\n\n");
            nmod_poly_print(h), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(f);
        nmod_poly_clear(g);
        nmod_poly_clear(h);
    }

    TEST_FUNCTION_END(state);
}
