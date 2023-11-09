/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpz_poly_mat.h"

TEST_FUNCTION_START(fmpz_poly_mat_rank, state)
{
    slong i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_mat_t A;
        fmpz_mat_t Ax;
        fmpz_t x;
        slong j, m, n, bits, deg, rank, zrank;
        float density;

        m = n_randint(state, 15);
        n = n_randint(state, 15);
        deg = 1 + n_randint(state, 5);
        bits = 1 + n_randint(state, 100);
        density = n_randint(state, 100) * 0.01;

        fmpz_poly_mat_init(A, m, n);
        fmpz_mat_init(Ax, m, n);
        fmpz_init(x);

        fmpz_poly_mat_randtest_sparse(A, state, deg, bits, density);

        /* Probabilistic rank computation */
        zrank = 0;
        for (j = 0; j < 5; j++)
        {
            slong r;
            fmpz_randbits(x, state, 15);
            fmpz_poly_mat_evaluate_fmpz(Ax, A, x);
            r = fmpz_mat_rank(Ax);
            zrank = FLINT_MAX(zrank, r);
        }

        rank = fmpz_poly_mat_rank(A);

        if (rank != zrank)
        {
            flint_printf("FAIL:\n");
            flint_printf("A:\n");
            fmpz_poly_mat_print(A, "x");
            flint_printf("Computed rank: %wd (zrank = %wd)\n", rank, zrank);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(x);
        fmpz_mat_clear(Ax);
        fmpz_poly_mat_clear(A);
    }

    TEST_FUNCTION_END(state);
}
