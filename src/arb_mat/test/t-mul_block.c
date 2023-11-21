/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_mat.h"
#include "arb_mat.h"

FLINT_DLL extern slong arb_mat_mul_block_min_block_size;

TEST_FUNCTION_START(arb_mat_mul_block, state)
{
    slong iter;

    for (iter = 0; iter < 5000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n, k, qbits1, qbits2, rbits1, rbits2, rbits3;
        fmpq_mat_t A, B, C;
        arb_mat_t a, b, c, d;

        qbits1 = 2 + n_randint(state, 200);
        qbits2 = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);
        rbits3 = 2 + n_randint(state, 200);

        arb_mat_mul_block_min_block_size = n_randint(state, 10);

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
        arb_mat_mul_block(c, a, b, rbits3);

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
            arb_mat_mul_block(d, d, b, rbits3);
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
            arb_mat_mul_block(d, a, d, rbits3);
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

    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_mat_t A, B, C, D;
        slong m, n, p, bits1, bits2, exp1, exp2, prec1, prec2;

        m = n_randint(state, 40);
        n = n_randint(state, 40);
        p = n_randint(state, 40);

        arb_mat_mul_block_min_block_size = n_randint(state, 10);

        if (n_randint(state, 4) == 0)
        {
            exp1 = 4 + n_randint(state, FLINT_BITS);
            exp2 = 4 + n_randint(state, FLINT_BITS);
        }
        else
        {
            exp1 = exp2 = 20;
        }

        bits1 = 2 + n_randint(state, 200);
        bits2 = 2 + n_randint(state, 200);
        prec1 = 2 + n_randint(state, 200);
        prec2 = 2 + n_randint(state, 200);

        arb_mat_init(A, m, n);
        arb_mat_init(B, n, p);
        arb_mat_init(C, m, p);
        arb_mat_init(D, m, p);

        arb_mat_randtest(A, state, bits1, exp1);
        arb_mat_randtest(B, state, bits2, exp2);
        arb_mat_randtest(C, state, bits2, exp2);

        arb_mat_mul_block(C, A, B, prec1);
        arb_mat_mul_classical(D, A, B, prec2);

        if (!arb_mat_overlaps(C, D))
        {
            flint_printf("FAIL\n");
            flint_printf("m = %wd, n = %wd, p = %wd\n", m, n, p);
            flint_abort();
        }

        arb_mat_clear(A);
        arb_mat_clear(B);
        arb_mat_clear(C);
        arb_mat_clear(D);
    }

    TEST_FUNCTION_END(state);
}
