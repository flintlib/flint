/*
    Copyright (C) 2011 Fredrik Johansson, Rémi Prébet

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

TEST_FUNCTION_START(fmpq_poly_interpolate_fmpq_vec, state)
{
    int i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t P, Q;
        fmpq *x, *y;
        fmpq_t z;
        slong j, npoints, n, bits;
        int result;

        npoints = n_randint(state, 200);
        n = n_randint(state, npoints + 1);
        bits = n_randint(state, 400) + 1;

        x = _fmpq_vec_init(npoints);
        y = _fmpq_vec_init(npoints);
        fmpq_init(z);

        fmpq_poly_init(P);
        fmpq_poly_init(Q);

        fmpq_poly_randtest(P, state, n, bits);

        for (j = 0; j < npoints; j++)
            fmpq_set_si(x + j, -npoints/2 + j, 1);
        for (j = 0; j < npoints; j++)
            fmpq_poly_evaluate_fmpq(y + j, P, x + j);

        result = fmpq_poly_interpolate_fmpq_vec(Q, x, y, npoints);

        if (!result)
        {
            flint_printf("FAIL (exact):\n");
            flint_printf("P  %{fmpq_poly}\n\n", P);
            flint_printf("x  %{fmpq*}\n\n", x, npoints);
            flint_printf("y  %{fmpq*}\n\n", y, npoints);
            fflush(stdout);
            flint_abort();
        }

        result = fmpq_poly_equal(P, Q);
        if (!result)
        {
            flint_printf("FAIL (P != Q):\n");
            fmpq_poly_print(P), flint_printf("\n\n");
            fmpq_poly_print(Q), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        /* Test arbitrary x and y */
        if (n_randint(state, 2))
            _fmpq_vec_randtest(x, state, npoints, 1 + n_randint(state, 10));
        _fmpq_vec_randtest(y, state, npoints, 1 + n_randint(state, 10));

        result = fmpq_poly_interpolate_fmpq_vec(Q, x, y, npoints);

        if (result)
        {
            for (j = 0; j < npoints && result; j++)
            {
                fmpq_poly_evaluate_fmpq(z, Q, x + j);
                result = fmpq_equal(y + j, z);
            }
            if (!result)
            {
                flint_printf("FAIL (P != Q, 2):\n");
                flint_printf("x  %{fmpq*}\n\n", x, npoints);
                flint_printf("y  %{fmpq*}\n\n", y, npoints);
                flint_printf("Q  %{fmpq_poly}\n\n", Q);
                fflush(stdout);
                flint_abort();
            }
        }


        fmpq_poly_clear(P);
        fmpq_poly_clear(Q);
        _fmpq_vec_clear(x, npoints);
        _fmpq_vec_clear(y, npoints);
        fmpq_clear(z);
    }

    TEST_FUNCTION_END(state);
}
