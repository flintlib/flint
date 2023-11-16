/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "fmpq_vec.h"
#include "fmpq_poly.h"
#include "arb_poly.h"

TEST_FUNCTION_START(arb_poly_evaluate_vec_fast, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong i, n, qbits1, qbits2, rbits1, rbits2, rbits3;
        fmpq_poly_t F;
        fmpq * X, * Y;
        arb_poly_t f;
        arb_ptr x, y;

        qbits1 = 2 + n_randint(state, 100);
        qbits2 = 2 + n_randint(state, 100);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);
        rbits3 = 2 + n_randint(state, 200);

        n = n_randint(state, 10);

        fmpq_poly_init(F);
        X = _fmpq_vec_init(n);
        Y = _fmpq_vec_init(n);

        arb_poly_init(f);
        x = _arb_vec_init(n);
        y = _arb_vec_init(n);

        fmpq_poly_randtest(F, state, 1 + n_randint(state, 20), qbits1);
        for (i = 0; i < n; i++)
            fmpq_randtest(X + i, state, qbits2);
        for (i = 0; i < n; i++)
            fmpq_poly_evaluate_fmpq(Y + i, F, X + i);

        arb_poly_set_fmpq_poly(f, F, rbits1);
        for (i = 0; i < n; i++)
            arb_set_fmpq(x + i, X + i, rbits2);
        arb_poly_evaluate_vec_fast(y, f, x, n, rbits3);

        for (i = 0; i < n; i++)
        {
            if (!arb_contains_fmpq(y + i, Y + i))
            {
                flint_printf("FAIL (%wd of %wd)\n\n", i, n);

                flint_printf("F = "); fmpq_poly_print(F); flint_printf("\n\n");
                flint_printf("X = "); fmpq_print(X + i); flint_printf("\n\n");
                flint_printf("Y = "); fmpq_print(Y + i); flint_printf("\n\n");

                flint_printf("f = "); arb_poly_printd(f, 15); flint_printf("\n\n");
                flint_printf("x = "); arb_printd(x + i, 15); flint_printf("\n\n");
                flint_printf("y = "); arb_printd(y + i, 15); flint_printf("\n\n");

                flint_abort();
            }
        }

        fmpq_poly_clear(F);
        _fmpq_vec_clear(X, n);
        _fmpq_vec_clear(Y, n);

        arb_poly_clear(f);
        _arb_vec_clear(x, n);
        _arb_vec_clear(y, n);
    }

    TEST_FUNCTION_END(state);
}
