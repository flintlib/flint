/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

TEST_FUNCTION_START(fmpz_poly_interpolate_fmpz_vec, state)
{
    int i, result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t P, Q;
        fmpz *x, *y;
        slong j, n, npoints, bits;

        npoints = n_randint(state, 50);
        n = n_randint(state, npoints + 1);
        bits = n_randint(state, 100);

        x = _fmpz_vec_init(npoints);
        y = _fmpz_vec_init(npoints);

        fmpz_poly_init(P);
        fmpz_poly_init(Q);

        fmpz_poly_randtest(P, state, n, bits);

        for (j = 0; j < npoints; j++)
            fmpz_set_si(x + j, -npoints/2 + j);

        fmpz_poly_evaluate_fmpz_vec(y, P, x, npoints);
        fmpz_poly_interpolate_fmpz_vec(Q, x, y, npoints);

        result = (fmpz_poly_equal(P, Q));
        if (!result)
        {
            flint_printf("FAIL (P != Q):\n");
            fmpz_poly_print(P), flint_printf("\n\n");
            fmpz_poly_print(Q), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(P);
        fmpz_poly_clear(Q);
        _fmpz_vec_clear(x, npoints);
        _fmpz_vec_clear(y, npoints);
    }

    TEST_FUNCTION_END(state);
}
