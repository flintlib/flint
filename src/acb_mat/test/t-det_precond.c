/*
    Copyright (C) 2018 Fredrik Johansson

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

TEST_FUNCTION_START(acb_mat_det_precond, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpq_mat_t Q;
        fmpq_t Qdet;
        acb_mat_t A;
        acb_t Adet;
        slong n, qbits, prec;

        n = n_randint(state, 12);
        qbits = 1 + n_randint(state, 100);
        prec = 2 + n_randint(state, 200);

        fmpq_mat_init(Q, n, n);
        fmpq_init(Qdet);

        acb_mat_init(A, n, n);
        acb_init(Adet);

        fmpq_mat_randtest(Q, state, qbits);
        fmpq_mat_det(Qdet, Q);

        acb_mat_set_fmpq_mat(A, Q, prec);
        acb_mat_det_precond(Adet, A, prec);

        if (!acb_contains_fmpq(Adet, Qdet))
        {
            flint_printf("FAIL (containment, iter = %wd)\n", iter);
            flint_printf("n = %wd, prec = %wd\n", n, prec);
            flint_printf("\n");

            flint_printf("Q = \n"); fmpq_mat_print(Q); flint_printf("\n\n");
            flint_printf("Qdet = \n"); fmpq_print(Qdet); flint_printf("\n\n");

            flint_printf("A = \n"); acb_mat_printd(A, 15); flint_printf("\n\n");
            flint_printf("Adet = \n"); acb_printd(Adet, 15); flint_printf("\n\n");
            flint_printf("Adet = \n"); acb_print(Adet); flint_printf("\n\n");

            flint_abort();
        }

        fmpq_mat_clear(Q);
        fmpq_clear(Qdet);
        acb_mat_clear(A);
        acb_clear(Adet);
    }

    for (iter = 0; iter < 5000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_mat_t A, B, AB;
        acb_t detA, detB, detAB, detAb, detBb, detABb, t;
        slong n, prec1, prec2, prec3;

        n = n_randint(state, 12);
        prec1 = 2 + n_randint(state, 200);
        prec2 = 2 + n_randint(state, 200);
        prec3 = 2 + n_randint(state, 200);

        acb_mat_init(A, n, n);
        acb_mat_init(B, n, n);
        acb_mat_init(AB, n, n);
        acb_init(detA);
        acb_init(detB);
        acb_init(detAB);
        acb_init(t);
        acb_init(detAb);
        acb_init(detBb);
        acb_init(detABb);

        acb_mat_randtest(A, state, 2 + n_randint(state, 200), 2 + n_randint(state, 20));
        acb_mat_randtest(B, state, 2 + n_randint(state, 200), 2 + n_randint(state, 20));
        acb_mat_mul(AB, A, B, prec3);

        acb_mat_det_precond(detA, A, prec1);
        acb_mat_det_precond(detB, B, prec2);
        acb_mat_det_precond(detAB, AB, prec3);

        acb_mul(t, detA, detB, 1000);

        if (!acb_overlaps(t, detAB))
        {
            flint_printf("FAIL (overlap, iter = %wd)\n", iter);
            flint_printf("n = %wd, prec1 = %wd, prec2 = %wd, prec3 = %wd\n", n, prec1, prec2, prec3);
            flint_printf("\n");

            flint_printf("A = \n"); acb_mat_printd(A, 15); flint_printf("\n\n");
            flint_printf("detA = \n"); acb_printn(detA, 50, 0); flint_printf("\n\n");

            flint_printf("B = \n"); acb_mat_printd(B, 15); flint_printf("\n\n");
            flint_printf("detB = \n"); acb_printn(detB, 50, 0); flint_printf("\n\n");

            flint_printf("A = \n"); acb_mat_printd(AB, 15); flint_printf("\n\n");
            flint_printf("detAB = \n"); acb_printn(detAB, 50, 0); flint_printf("\n\n");

            flint_printf("detA*detB = \n"); acb_printn(t, 50, 0); flint_printf("\n\n");

            flint_abort();
        }

        acb_mat_det_lu(detAb, A, prec1);
        acb_mat_det_lu(detBb, B, prec2);
        acb_mat_det_lu(detABb, AB, prec3);

        if (!acb_overlaps(detA, detAb) || !acb_overlaps(detB, detBb) || !acb_overlaps(detAB, detABb))
        {
            flint_printf("FAIL (overlap, iter = %wd)\n", iter);
            flint_printf("n = %wd, prec1 = %wd, prec2 = %wd, prec3 = %wd\n", n, prec1, prec2, prec3);
            flint_printf("\n");

            flint_printf("A = \n"); acb_mat_printd(A, 15); flint_printf("\n\n");
            flint_printf("detA = \n"); acb_printn(detA, 50, 0); flint_printf("\n\n");
            flint_printf("detAb = \n"); acb_printn(detAb, 50, 0); flint_printf("\n\n");

            flint_printf("B = \n"); acb_mat_printd(B, 15); flint_printf("\n\n");
            flint_printf("detB = \n"); acb_printn(detB, 50, 0); flint_printf("\n\n");
            flint_printf("detBb = \n"); acb_printn(detBb, 50, 0); flint_printf("\n\n");

            flint_printf("A = \n"); acb_mat_printd(AB, 15); flint_printf("\n\n");
            flint_printf("detAB = \n"); acb_printn(detAB, 50, 0); flint_printf("\n\n");
            flint_printf("detABb = \n"); acb_printn(detABb, 50, 0); flint_printf("\n\n");

            flint_abort();
        }

        acb_mat_clear(A);
        acb_mat_clear(B);
        acb_mat_clear(AB);
        acb_clear(detA);
        acb_clear(detB);
        acb_clear(detAB);
        acb_clear(t);
        acb_clear(detAb);
        acb_clear(detBb);
        acb_clear(detABb);
    }

    TEST_FUNCTION_END(state);
}
