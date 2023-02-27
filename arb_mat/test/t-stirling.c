/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"
#include "flint/arith.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("stirling....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100 * arb_test_multiplier(); iter++)
    {
        arb_mat_t A;
        fmpz_mat_t B;
        slong n, m, prec;

        n = n_randint(state, 10);
        m = n_randint(state, 10);
        prec = 2 + n_randint(state, 200);

        arb_mat_init(A, n, m);
        fmpz_mat_init(B, n, m);
        arb_mat_randtest(A, state, 100, 10);

        arb_mat_stirling(A, 0, prec);
        arith_stirling_matrix_1u(B);
        if (!arb_mat_contains_fmpz_mat(A, B))
        {
            flint_printf("FAIL: containment (0)\n");
            flint_abort();
        }

        arb_mat_stirling(A, 1, prec);
        arith_stirling_matrix_1(B);
        if (!arb_mat_contains_fmpz_mat(A, B))
        {
            flint_printf("FAIL: containment (1)\n");
            flint_abort();
        }

        arb_mat_stirling(A, 2, prec);
        arith_stirling_matrix_2(B);
        if (!arb_mat_contains_fmpz_mat(A, B))
        {
            flint_printf("FAIL: containment (2)\n");
            flint_abort();
        }

        arb_mat_clear(A);
        fmpz_mat_clear(B);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
