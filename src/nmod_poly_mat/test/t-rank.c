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

TEST_FUNCTION_START(nmod_poly_mat_rank, state)
{
    slong i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_mat_t A;
        nmod_mat_t Ax;
        mp_limb_t mod, x;
        slong j, m, n, deg, rank, zrank;
        float density;

        /* Don't pick a too small modulus, to avoid failure in
           the probabilistic rank computation (todo: test
           for small moduli) */
        do {
            mod = n_randtest_prime(state, 0);
        } while (mod < 20);

        m = n_randint(state, 15);
        n = n_randint(state, 15);
        deg = 1 + n_randint(state, 5);
        density = n_randint(state, 100) * 0.01;

        nmod_poly_mat_init(A, m, n, mod);
        nmod_mat_init(Ax, m, n, mod);

        nmod_poly_mat_randtest_sparse(A, state, deg, density);

        /* Probabilistic rank computation */
        zrank = 0;
        for (j = 0; j < 5; j++)
        {
            slong r;
            x = n_randint(state, mod);
            nmod_poly_mat_evaluate_nmod(Ax, A, x);
            r = nmod_mat_rank(Ax);
            zrank = FLINT_MAX(zrank, r);
        }

        rank = nmod_poly_mat_rank(A);

        if (rank != zrank)
        {
            flint_printf("FAIL:\n");
            flint_printf("A:\n");
            nmod_poly_mat_print(A, "x");
            flint_printf("Computed rank: %wd (zrank = %wd)\n", rank, zrank);
            fflush(stdout);
            flint_abort();
        }

        nmod_mat_clear(Ax);
        nmod_poly_mat_clear(A);
    }

    TEST_FUNCTION_END(state);
}
