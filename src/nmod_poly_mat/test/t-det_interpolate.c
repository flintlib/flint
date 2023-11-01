/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_poly_mat.h"

TEST_FUNCTION_START(nmod_poly_mat_det_interpolate, state)
{
    slong i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_mat_t A;
        nmod_poly_t a, b;
        slong n, deg;
        mp_limb_t mod;

        mod = n_randtest_prime(state, 0);
        n = n_randint(state, 10);
        deg = 1 + n_randint(state, 5);

        nmod_poly_mat_init(A, n, n, mod);

        nmod_poly_init(a, mod);
        nmod_poly_init(b, mod);

        nmod_poly_mat_randtest(A, state, deg);

        nmod_poly_mat_det(a, A);
        nmod_poly_mat_det_interpolate(b, A);

        if (!nmod_poly_equal(a, b))
        {
            flint_printf("FAIL:\n");
            flint_printf("determinants don't agree!\n");
            flint_printf("A:\n");
            nmod_poly_mat_print(A, "x");
            flint_printf("det(A):\n");
            nmod_poly_print(a);
            flint_printf("\ndet_interpolate(A):\n");
            nmod_poly_print(b);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);

        nmod_poly_mat_clear(A);
    }

    TEST_FUNCTION_END(state);
}
