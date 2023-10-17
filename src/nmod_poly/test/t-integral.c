/*
    Copyright (C) 2010 William Hart
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

TEST_FUNCTION_START(nmod_poly_integral, state)
{
    int i, result = 1;

    /* Check with derivative */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c;
        ulong len;
        mp_limb_t c0;
        mp_limb_t n = n_randtest_prime(state, 0);

        len = n_randint(state, 100);
        len = FLINT_MIN(len, n - 1);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);

        nmod_poly_randtest(a, state, len);

        nmod_poly_integral(b, a);
        c0 = nmod_poly_get_coeff_ui(b, 0);
        nmod_poly_derivative(c, b);

        result = (c0 == UWORD(0)) && nmod_poly_equal(a, c);

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a->length = %wd, n = %wu\n", a->length, a->mod.n);
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(c), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
    }

    /* Check aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b;
        mp_limb_t n = n_randtest_prime(state, 0);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_randtest(a, state, n_randint(state, FLINT_MIN(100, n)));

        nmod_poly_integral(b, a);
        nmod_poly_integral(a, a);

        result = nmod_poly_equal(a, b);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a->length = %wd, n = %wu\n", a->length, a->mod.n);
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
    }

    TEST_FUNCTION_END(state);
}
