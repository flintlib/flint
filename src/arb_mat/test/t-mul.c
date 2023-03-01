/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"


int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("mul....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        slong m, n, k, qbits1, qbits2, rbits1, rbits2, rbits3;
        fmpq_mat_t A, B, C;
        arb_mat_t a, b, c, d;

        qbits1 = 2 + n_randint(state, 200);
        qbits2 = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);
        rbits3 = 2 + n_randint(state, 200);

        m = n_randint(state, 10);
        n = n_randint(state, 10);
        k = n_randint(state, 10);

        fmpq_mat_init(A, m, n);
        fmpq_mat_init(B, n, k);
        fmpq_mat_init(C, m, k);

        arb_mat_init(a, m, n);
        arb_mat_init(b, n, k);
        arb_mat_init(c, m, k);
        arb_mat_init(d, m, k);

        fmpq_mat_randtest(A, state, qbits1);
        fmpq_mat_randtest(B, state, qbits2);
        fmpq_mat_mul(C, A, B);

        arb_mat_set_fmpq_mat(a, A, rbits1);
        arb_mat_set_fmpq_mat(b, B, rbits2);
        arb_mat_mul(c, a, b, rbits3);

        if (!arb_mat_contains_fmpq_mat(c, C))
        {
            flint_printf("FAIL\n\n");
            flint_printf("m = %wd, n = %wd, k = %wd, bits3 = %wd\n", m, n, k, rbits3);

            flint_printf("A = "); fmpq_mat_print(A); flint_printf("\n\n");
            flint_printf("B = "); fmpq_mat_print(B); flint_printf("\n\n");
            flint_printf("C = "); fmpq_mat_print(C); flint_printf("\n\n");

            flint_printf("a = "); arb_mat_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); arb_mat_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); arb_mat_printd(c, 15); flint_printf("\n\n");

            flint_abort();
        }

        /* test aliasing with a */
        if (arb_mat_nrows(a) == arb_mat_nrows(c) &&
            arb_mat_ncols(a) == arb_mat_ncols(c))
        {
            arb_mat_set(d, a);
            arb_mat_mul(d, d, b, rbits3);
            if (!arb_mat_equal(d, c))
            {
                flint_printf("FAIL (aliasing 1)\n\n");
                flint_abort();
            }
        }

        /* test aliasing with b */
        if (arb_mat_nrows(b) == arb_mat_nrows(c) &&
            arb_mat_ncols(b) == arb_mat_ncols(c))
        {
            arb_mat_set(d, b);
            arb_mat_mul(d, a, d, rbits3);
            if (!arb_mat_equal(d, c))
            {
                flint_printf("FAIL (aliasing 2)\n\n");
                flint_abort();
            }
        }

        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
        fmpq_mat_clear(C);

        arb_mat_clear(a);
        arb_mat_clear(b);
        arb_mat_clear(c);
        arb_mat_clear(d);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
