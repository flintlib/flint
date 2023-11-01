/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_mat.h"
#include "nmod_poly_mat.h"

TEST_FUNCTION_START(nmod_poly_mat_set_nmod_mat, state)
{
    slong i;

    for (i = 0; i < 400 * flint_test_multiplier(); i++)
    {
        nmod_poly_mat_t pmat;
        nmod_mat_t cmat;
        mp_limb_t mod;
        slong m, n, deg;

        mod = n_randtest_prime(state, 0);
        m = n_randint(state, 20);
        n = n_randint(state, 20);
        deg = 1 + n_randint(state, 10);

        nmod_poly_mat_init(pmat, m, n, mod);
        nmod_mat_init(cmat, m, n, mod);

        nmod_poly_mat_set_nmod_mat(pmat, cmat);
        if (! nmod_poly_mat_equal_nmod_mat(pmat, cmat))
        {
            flint_printf("FAIL:\n");
            flint_printf("pmat:\n");
            nmod_poly_mat_print(pmat, "x");
            flint_printf("cmat:\n");
            nmod_mat_print(cmat);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_mat_randtest(pmat, state, deg);
        nmod_poly_mat_set_nmod_mat(pmat, cmat);
        if (! nmod_poly_mat_equal_nmod_mat(pmat, cmat))
        {
            flint_printf("FAIL:\n");
            flint_printf("pmat:\n");
            nmod_poly_mat_print(pmat, "x");
            flint_printf("cmat:\n");
            nmod_mat_print(cmat);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_mat_clear(pmat);
        nmod_mat_clear(cmat);
    }

    TEST_FUNCTION_END(state);
}
