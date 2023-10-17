/*
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

TEST_FUNCTION_START(nmod_poly_inflate, state)
{
    int iter;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        nmod_poly_t poly1, poly2, poly3, xp;
        mp_limb_t modulus;
        ulong inflation;

        modulus = n_randtest_prime(state, 0);

        nmod_poly_init(poly1, modulus);
        nmod_poly_init(poly2, modulus);
        nmod_poly_init(poly3, modulus);
        nmod_poly_init(xp, modulus);

        nmod_poly_randtest(poly1, state, n_randint(state, 20));
        inflation = n_randint(state, 10);

        nmod_poly_inflate(poly2, poly1, inflation);

        nmod_poly_set_coeff_ui(xp, inflation, 1);
        nmod_poly_compose(poly3, poly1, xp);

        if (!nmod_poly_equal(poly2, poly3))
        {
            flint_printf("FAIL: not equal to compose (inflation = %wu)\n", inflation);
            flint_printf("poly1:\n"); nmod_poly_print(poly1); flint_printf("\n\n");
            flint_printf("poly2:\n"); nmod_poly_print(poly2); flint_printf("\n\n");
            flint_printf("poly3:\n"); nmod_poly_print(poly3); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_inflate(poly1, poly1, inflation);
        if (!nmod_poly_equal(poly1, poly2))
        {
            flint_printf("FAIL: aliasing (inflation = %wu)\n", inflation);
            flint_printf("poly1:\n"); nmod_poly_print(poly1); flint_printf("\n\n");
            flint_printf("poly2:\n"); nmod_poly_print(poly2); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(poly1);
        nmod_poly_clear(poly2);
        nmod_poly_clear(poly3);
        nmod_poly_clear(xp);
    }

    TEST_FUNCTION_END(state);
}
