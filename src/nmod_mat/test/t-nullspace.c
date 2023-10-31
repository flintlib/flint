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
#include "nmod_mat.h"

TEST_FUNCTION_START(nmod_mat_nullspace, state)
{
    slong i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_mat_t A, B, ker;
        mp_limb_t mod;
        slong m, n, d, r, nullity, nulrank;

        m = n_randint(state, 30);
        n = n_randint(state, 30);

        for (r = 0; r <= FLINT_MIN(m,n); r++)
        {
            mod = n_randtest_prime(state, 0);
            d = n_randint(state, 2*m*n + 1);

            nmod_mat_init(A, m, n, mod);
            nmod_mat_init(ker, n, n, mod);
            nmod_mat_init(B, m, n, mod);

            nmod_mat_randrank(A, state, r);
            /* Densify */
            if (n_randlimb(state) % 2)
                nmod_mat_randops(A, d, state);

            nullity = nmod_mat_nullspace(ker, A);
            nulrank = nmod_mat_rank(ker);

            if (nullity != nulrank)
            {
                flint_printf("FAIL:\n");
                flint_printf("rank(ker) != nullity!\n");
                nmod_mat_print_pretty(A);
                flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            if (nullity + r != n)
            {
                flint_printf("FAIL:\n");
                flint_printf("nullity + rank != n\n");
                nmod_mat_print_pretty(A);
                flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            nmod_mat_mul(B, A, ker);

            if (nmod_mat_rank(B) != 0)
            {
                flint_printf("FAIL:\n");
                flint_printf("A * ker != 0\n");
                nmod_mat_print_pretty(A);
                flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            nmod_mat_clear(A);
            nmod_mat_clear(ker);
            nmod_mat_clear(B);
        }
    }

    TEST_FUNCTION_END(state);
}
