/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_mat.h"

TEST_FUNCTION_START(arb_mat_det_precond, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_mat_t A, B, AB;
        arb_t detA, detB, detAB, detAb, detBb, detABb, t;
        slong n, prec1, prec2, prec3;

        n = n_randint(state, 12);
        prec1 = 2 + n_randint(state, 200);
        prec2 = 2 + n_randint(state, 200);
        prec3 = 2 + n_randint(state, 200);

        arb_mat_init(A, n, n);
        arb_mat_init(B, n, n);
        arb_mat_init(AB, n, n);
        arb_init(detA);
        arb_init(detB);
        arb_init(detAB);
        arb_init(t);
        arb_init(detAb);
        arb_init(detBb);
        arb_init(detABb);

        arb_mat_randtest(A, state, 2 + n_randint(state, 200), 2 + n_randint(state, 20));
        arb_mat_randtest(B, state, 2 + n_randint(state, 200), 2 + n_randint(state, 20));
        arb_mat_mul(AB, A, B, prec3);

        arb_mat_det_precond(detA, A, prec1);
        arb_mat_det_precond(detB, B, prec2);
        arb_mat_det_precond(detAB, AB, prec3);

        arb_mul(t, detA, detB, 1000);

        if (!arb_overlaps(t, detAB))
        {
            flint_printf("FAIL (overlap, iter = %wd)\n", iter);
            flint_printf("n = %wd, prec1 = %wd, prec2 = %wd, prec3 = %wd\n", n, prec1, prec2, prec3);
            flint_printf("\n");

            flint_printf("A = \n"); arb_mat_printd(A, 15); flint_printf("\n\n");
            flint_printf("detA = \n"); arb_printn(detA, 50, 0); flint_printf("\n\n");

            flint_printf("B = \n"); arb_mat_printd(B, 15); flint_printf("\n\n");
            flint_printf("detB = \n"); arb_printn(detB, 50, 0); flint_printf("\n\n");

            flint_printf("A = \n"); arb_mat_printd(AB, 15); flint_printf("\n\n");
            flint_printf("detAB = \n"); arb_printn(detAB, 50, 0); flint_printf("\n\n");

            flint_printf("detA*detB = \n"); arb_printn(t, 50, 0); flint_printf("\n\n");

            flint_abort();
        }

        arb_mat_det_lu(detAb, A, prec1);
        arb_mat_det_lu(detBb, B, prec2);
        arb_mat_det_lu(detABb, AB, prec3);

        if (!arb_overlaps(detA, detAb) || !arb_overlaps(detB, detBb) || !arb_overlaps(detAB, detABb))
        {
            flint_printf("FAIL (overlap, iter = %wd)\n", iter);
            flint_printf("n = %wd, prec1 = %wd, prec2 = %wd, prec3 = %wd\n", n, prec1, prec2, prec3);
            flint_printf("\n");

            flint_printf("A = \n"); arb_mat_printd(A, 15); flint_printf("\n\n");
            flint_printf("detA = \n"); arb_printn(detA, 50, 0); flint_printf("\n\n");
            flint_printf("detAb = \n"); arb_printn(detAb, 50, 0); flint_printf("\n\n");

            flint_printf("B = \n"); arb_mat_printd(B, 15); flint_printf("\n\n");
            flint_printf("detB = \n"); arb_printn(detB, 50, 0); flint_printf("\n\n");
            flint_printf("detBb = \n"); arb_printn(detBb, 50, 0); flint_printf("\n\n");

            flint_printf("A = \n"); arb_mat_printd(AB, 15); flint_printf("\n\n");
            flint_printf("detAB = \n"); arb_printn(detAB, 50, 0); flint_printf("\n\n");
            flint_printf("detABb = \n"); arb_printn(detABb, 50, 0); flint_printf("\n\n");

            flint_abort();
        }

        arb_mat_clear(A);
        arb_mat_clear(B);
        arb_mat_clear(AB);
        arb_clear(detA);
        arb_clear(detB);
        arb_clear(detAB);
        arb_clear(t);
        arb_clear(detAb);
        arb_clear(detBb);
        arb_clear(detABb);
    }

    TEST_FUNCTION_END(state);
}
