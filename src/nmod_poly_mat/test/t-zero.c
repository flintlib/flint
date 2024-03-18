/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly_mat.h"

TEST_FUNCTION_START(nmod_poly_mat_zero, state)
{
    int iter;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        nmod_poly_mat_t A;
        slong m, n;
        mp_limb_t mod;

        mod = n_randtest_prime(state, 0);
        m = n_randint(state, 10);
        n = n_randint(state, 10);

        nmod_poly_mat_init(A, m, n, mod);
        nmod_poly_mat_randtest(A, state, n_randint(state, 5));
        nmod_poly_mat_zero(A);

        if (!nmod_poly_mat_is_zero(A))
        {
            flint_printf("FAIL: expected matrix to be zero\n");
            fflush(stdout);
            flint_abort();
        }

        if (m > 0 && n > 0)
        {
            m = n_randint(state, m);
            n = n_randint(state, n);
            nmod_poly_randtest_not_zero(nmod_poly_mat_entry(A, m, n), state, 5);

            if (nmod_poly_mat_is_zero(A))
            {
                flint_printf("FAIL: expected matrix not to be zero\n");
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_poly_mat_clear(A);
    }

    TEST_FUNCTION_END(state);
}
