/*
    Copyright (C) 2025 Fredrik Johansson, Rémi Prébet

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_vec.h"
#include "fmpq.h"
#include "fmpq_poly.h"

TEST_FUNCTION_START(fmpq_poly_interpolate_fast, state)
{
    int i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t P, Q;
        fmpq *x, *y;
        slong j, npoints, n, bits;

        npoints = n_randint(state, 200);
        n = n_randint(state, npoints + 1);
        bits = n_randint(state, 400) + 1;

        x = _fmpq_vec_init(npoints);
        y = _fmpq_vec_init(npoints);

        fmpq_poly_init(P);
        fmpq_poly_init(Q);

        fmpq_poly_randtest(P, state, n, bits);

        for (j = 0; j < npoints; j++)
            fmpq_set_si(x + j, -npoints/2 + j, 1);
        for (j = 0; j < npoints; j++)
            fmpq_poly_evaluate_fmpq(y + j, P, x + j);

        fmpq_poly_interpolate_fast(Q, x, y, npoints);

        if (!fmpq_poly_equal(P, Q))
        {
            flint_printf("FAIL (P != Q):\n");
            fmpq_poly_print(P), flint_printf("\n\n");
            fmpq_poly_print(Q), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(P);
        fmpq_poly_clear(Q);
        _fmpq_vec_clear(x, npoints);
        _fmpq_vec_clear(y, npoints);
    }

    TEST_FUNCTION_END(state);
}
