/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "fmpq_mat.h"
#include "acb_mat.h"

void _fmpq_mat_randtest_for_exp(fmpq_mat_t mat, flint_rand_t state, flint_bitcnt_t bits)
{
    slong i, j;
    slong d, l, u;
    d = n_randint(state, 5);
    l = n_randint(state, 5);
    u = n_randint(state, 5);
    fmpq_mat_zero(mat);
    for (i = 0; i < fmpq_mat_nrows(mat); i++)
    {
        for (j = 0; j < fmpq_mat_ncols(mat); j++)
        {
            if ((i == j && d) || (i < j && u) || (i > j && l))
            {
                fmpq_randtest(fmpq_mat_entry(mat, i, j), state, bits);
            }
        }
    }
}

TEST_FUNCTION_START(acb_mat_exp, state)
{
    slong iter;

    /* check exp(A)*exp(c*A) = exp((1+c)*A) */
    for (iter = 0; iter < 500 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_mat_t A, E, F, EF, G;
        fmpq_mat_t Q;
        acb_t c, d;
        slong n, qbits, prec;

        n = n_randint(state, 5);
        qbits = 2 + n_randint(state, 300);
        prec = 2 + n_randint(state, 300);

        acb_init(c);
        acb_init(d);
        fmpq_mat_init(Q, n, n);
        acb_mat_init(A, n, n);
        acb_mat_init(E, n, n);
        acb_mat_init(F, n, n);
        acb_mat_init(EF, n, n);
        acb_mat_init(G, n, n);

        _fmpq_mat_randtest_for_exp(Q, state, qbits);
        acb_mat_set_fmpq_mat(A, Q, prec);

        acb_mat_exp(E, A, prec);

        acb_randtest(c, state, prec, 10);
        acb_mat_scalar_mul_acb(F, A, c, prec);
        acb_mat_exp(F, F, prec);

        acb_add_ui(d, c, 1, prec);
        acb_mat_scalar_mul_acb(G, A, d, prec);
        acb_mat_exp(G, G, prec);

        acb_mat_mul(EF, E, F, prec);

        if (!acb_mat_overlaps(EF, G))
        {
            flint_printf("FAIL\n\n");
            flint_printf("n = %wd, prec = %wd\n", n, prec);

            flint_printf("c = \n"); acb_printd(c, 15); flint_printf("\n\n");

            flint_printf("A = \n"); acb_mat_printd(A, 15); flint_printf("\n\n");
            flint_printf("E   = \n"); acb_mat_printd(E, 15); flint_printf("\n\n");
            flint_printf("F   = \n"); acb_mat_printd(F, 15); flint_printf("\n\n");
            flint_printf("E*F = \n"); acb_mat_printd(EF, 15); flint_printf("\n\n");
            flint_printf("G   = \n"); acb_mat_printd(G, 15); flint_printf("\n\n");

            flint_abort();
        }

        acb_clear(c);
        acb_clear(d);
        fmpq_mat_clear(Q);
        acb_mat_clear(A);
        acb_mat_clear(E);
        acb_mat_clear(F);
        acb_mat_clear(EF);
        acb_mat_clear(G);
    }

    TEST_FUNCTION_END(state);
}
