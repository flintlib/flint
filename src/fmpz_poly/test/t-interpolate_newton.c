/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

TEST_FUNCTION_START(fmpz_poly_interpolate_newton, state)
{
    int i, result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t P, Q;
        fmpz *x, *y, *z;
        slong j, n, npoints, bits;

        npoints = n_randint(state, 10);
        n = n_randint(state, npoints + 1);
        bits = n_randint(state, 400);

        x = _fmpz_vec_init(npoints);
        y = _fmpz_vec_init(npoints);
        z = _fmpz_vec_init(npoints);

        fmpz_poly_init(P);
        fmpz_poly_init(Q);

        fmpz_poly_randtest(P, state, n, bits);

        for (j = 0; j < npoints; j++)
            fmpz_set_si(x + j, -npoints/2 + j);

        fmpz_poly_evaluate_fmpz_vec(y, P, x, npoints);

        if (n_randint(state, 2))
        {
            fmpz_poly_interpolate_exact_newton(Q, x, y, npoints);
            result = 1;
        }
        else
            result = fmpz_poly_interpolate_newton(Q, x, y, npoints);

        if (!result)
        {
            flint_printf("FAIL (exact):\n");
            flint_printf("P  %{fmpz_poly}\n\n", P);
            flint_printf("x  %{fmpz*}\n\n", x, npoints);
            flint_printf("y  %{fmpz*}\n\n", y, npoints);
            fflush(stdout);
            flint_abort();
        }

        result = (fmpz_poly_equal(P, Q));
        if (!result)
        {
            flint_printf("FAIL (P != Q):\n");
            flint_printf("P  %{fmpz_poly}\n\n", P);
            flint_printf("x  %{fmpz*}\n\n", x, npoints);
            flint_printf("y  %{fmpz*}\n\n", y, npoints);
            flint_printf("Q  %{fmpz_poly}\n\n", Q);
            fflush(stdout);
            flint_abort();
        }

        /* Test arbitrary x and y */
        if (n_randint(state, 2))
            _fmpz_vec_randtest(x, state, npoints, 1 + n_randint(state, 10));
        _fmpz_vec_randtest(y, state, npoints, 1 + n_randint(state, 10));

        result = fmpz_poly_interpolate_newton(Q, x, y, npoints);

        if (result)
        {
            fmpz_poly_evaluate_fmpz_vec(z, Q, x, npoints);

            result = (_fmpz_vec_equal(z, y, npoints));

            if (!result)
            {
                flint_printf("FAIL (P != Q, 2):\n");
                flint_printf("x  %{fmpz*}\n\n", x, npoints);
                flint_printf("y  %{fmpz*}\n\n", y, npoints);
                flint_printf("Q  %{fmpz_poly}\n\n", Q);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_poly_clear(P);
        fmpz_poly_clear(Q);
        _fmpz_vec_clear(x, npoints);
        _fmpz_vec_clear(y, npoints);
        _fmpz_vec_clear(z, npoints);
    }

    TEST_FUNCTION_END(state);
}
