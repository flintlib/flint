/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ca_mat.h"

TEST_FUNCTION_START(ca_mat_det, state)
{
    slong iter;

    for (iter = 0; iter < 100 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_mat_t A;
        ca_t detA, detB;
        slong n;

        ca_ctx_init(ctx);

        n = n_randint(state, 6);

        ca_mat_init(A, n, n, ctx);
        ca_init(detA, ctx);
        ca_init(detB, ctx);

        ca_mat_randtest(A, state, 3, 5, ctx);

        switch (n_randint(state, 5))
        {
            case 0: ca_mat_det_berkowitz(detA, A, ctx); break;
            case 1: ca_mat_det_lu(detA, A, ctx); break;
            case 2: ca_mat_det_bareiss(detA, A, ctx); break;
            default: ca_mat_det(detA, A, ctx); break;
        }

        switch (n_randint(state, 5))
        {
            case 0: ca_mat_det_berkowitz(detB, A, ctx); break;
            case 1: ca_mat_det_lu(detB, A, ctx); break;
            case 2: ca_mat_det_bareiss(detB, A, ctx); break;
            default: ca_mat_det(detB, A, ctx); break;
        }

        if (ca_check_equal(detA, detB, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n\n");
            flint_printf("A = "); ca_mat_print(A, ctx); flint_printf("\n");
            flint_printf("detA = "); ca_print(detA, ctx); flint_printf("\n");
            flint_printf("detB = "); ca_print(detB, ctx); flint_printf("\n");
            flint_abort();
        }

        ca_mat_clear(A, ctx);
        ca_clear(detA, ctx);
        ca_clear(detB, ctx);

        ca_ctx_clear(ctx);
    }

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_mat_t A, B, AB;
        ca_t detA, detB, detAB, detAdetB;
        slong n;

        ca_ctx_init(ctx);

        n = n_randint(state, 5);

        ca_mat_init(A, n, n, ctx);
        ca_mat_init(B, n, n, ctx);
        ca_mat_init(AB, n, n, ctx);
        ca_init(detA, ctx);
        ca_init(detB, ctx);
        ca_init(detAB, ctx);
        ca_init(detAdetB, ctx);

        if (n_randint(state, 2))
            ca_mat_randtest_rational(A, state, 100, ctx);
        else
            ca_mat_randtest(A, state, 3, 5, ctx);

        if (n_randint(state, 2))
            ca_mat_randtest_rational(B, state, 10, ctx);
        else
            ca_mat_randtest(B, state, 3, 5, ctx);

        ca_mat_mul(AB, A, B, ctx);

        switch (n_randint(state, 5))
        {
            case 0: ca_mat_det_berkowitz(detA, A, ctx); break;
            case 1: ca_mat_det_lu(detA, A, ctx); break;
            case 2: ca_mat_det_bareiss(detA, A, ctx); break;
            default: ca_mat_det(detA, A, ctx); break;
        }

        switch (n_randint(state, 5))
        {
            case 0: ca_mat_det_berkowitz(detB, B, ctx); break;
            case 1: ca_mat_det_lu(detB, B, ctx); break;
            case 2: ca_mat_det_bareiss(detB, B, ctx); break;
            default: ca_mat_det(detB, B, ctx); break;
        }

        switch (n_randint(state, 5))
        {
            case 0: ca_mat_det_berkowitz(detAB, AB, ctx); break;
            case 1: ca_mat_det_lu(detAB, AB, ctx); break;
            case 2: ca_mat_det_bareiss(detAB, AB, ctx); break;
            default: ca_mat_det(detAB, AB, ctx); break;
        }

        ca_mul(detAdetB, detA, detB, ctx);

        if (ca_check_equal(detAB, detAdetB, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n\n");
            flint_printf("A = "); ca_mat_print(A, ctx); flint_printf("\n");
            flint_printf("B = "); ca_mat_print(B, ctx); flint_printf("\n");
            flint_printf("AB = "); ca_mat_print(AB, ctx); flint_printf("\n");
            flint_printf("detA = "); ca_print(detA, ctx); flint_printf("\n");
            flint_printf("detB = "); ca_print(detB, ctx); flint_printf("\n");
            flint_printf("detAB = "); ca_print(detAB, ctx); flint_printf("\n");
            flint_printf("detAdetB = "); ca_print(detAdetB, ctx); flint_printf("\n");
            flint_abort();
        }

        ca_mat_clear(A, ctx);
        ca_mat_clear(B, ctx);
        ca_mat_clear(AB, ctx);
        ca_clear(detA, ctx);
        ca_clear(detB, ctx);
        ca_clear(detAB, ctx);
        ca_clear(detAdetB, ctx);

        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
