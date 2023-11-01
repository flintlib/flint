/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_mat.h"

TEST_FUNCTION_START(nmod_mat_rank, state)
{
    nmod_mat_t A;
    slong i, m, n, d, r;
    mp_limb_t mod;

    /* Maximally sparse matrices of given rank */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 20);
        n = n_randint(state, 20);
        mod = n_randtest_prime(state, 0);

        for (r = 0; r <= FLINT_MIN(m,n); r++)
        {
            nmod_mat_init(A, m, n, mod);
            nmod_mat_randrank(A, state, r);
            /* flint_printf("SPARSE %wd\n", r);
            nmod_mat_print_pretty(A); */
            if (r != nmod_mat_rank(A))
            {
                flint_printf("FAIL:\n");
                flint_printf("wrong rank!\n");
                fflush(stdout);
                flint_abort();
            }
            nmod_mat_clear(A);
        }
    }

    /* Dense */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 20);
        n = n_randint(state, 20);
        mod = n_randtest_prime(state, 0);

        for (r = 0; r <= FLINT_MIN(m,n); r++)
        {
            d = n_randint(state, 2*m*n + 1);
            nmod_mat_init(A, m, n, mod);
            nmod_mat_randrank(A, state, r);
            nmod_mat_randops(A, d, state);
            /*
            flint_printf("DENSE %wd %wd\n", r, d);
            nmod_mat_print_pretty(A); */
            if (r != nmod_mat_rank(A))
            {
                flint_printf("FAIL:\n");
                flint_printf("wrong rank!\n");
                fflush(stdout);
                flint_abort();
            }
            nmod_mat_clear(A);
        }
    }

    TEST_FUNCTION_END(state);
}
