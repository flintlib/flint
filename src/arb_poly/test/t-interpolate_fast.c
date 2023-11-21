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
#include "fmpq_poly.h"
#include "arb_poly.h"

TEST_FUNCTION_START(arb_poly_interpolate_fast, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong i, n, qbits1, qbits2, rbits1, rbits2, rbits3;
        fmpq_poly_t P;
        arb_poly_t R, S;
        fmpq_t t, u;
        arb_ptr xs, ys;

        fmpq_poly_init(P);
        arb_poly_init(R);
        arb_poly_init(S);
        fmpq_init(t);
        fmpq_init(u);

        qbits1 = 2 + n_randint(state, 200);
        qbits2 = 2 + n_randint(state, 5);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);
        rbits3 = 2 + n_randint(state, 200);

        fmpq_poly_randtest(P, state, 1 + n_randint(state, 20), qbits1);
        n = P->length;

        xs = _arb_vec_init(n);
        ys = _arb_vec_init(n);

        arb_poly_set_fmpq_poly(R, P, rbits1);

        if (n > 0)
        {
            fmpq_randtest(t, state, qbits2);
            arb_set_fmpq(xs, t, rbits2);

            for (i = 1; i < n; i++)
            {
                fmpq_randtest_not_zero(u, state, qbits2);
                fmpq_abs(u, u);
                fmpq_add(t, t, u);
                arb_set_fmpq(xs + i, t, rbits2);
            }
        }

        for (i = 0; i < n; i++)
            arb_poly_evaluate(ys + i, R, xs + i, rbits2);

        arb_poly_interpolate_fast(S, xs, ys, n, rbits3);

        if (!arb_poly_contains_fmpq_poly(S, P))
        {
            flint_printf("FAIL:\n");
            flint_printf("P = "); fmpq_poly_print(P); flint_printf("\n\n");
            flint_printf("R = "); arb_poly_printd(R, 15); flint_printf("\n\n");
            flint_printf("S = "); arb_poly_printd(S, 15); flint_printf("\n\n");
            flint_abort();
        }

        fmpq_poly_clear(P);
        arb_poly_clear(R);
        arb_poly_clear(S);
        fmpq_clear(t);
        fmpq_clear(u);
        _arb_vec_clear(xs, n);
        _arb_vec_clear(ys, n);
    }

    TEST_FUNCTION_END(state);
}
