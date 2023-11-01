/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_vec.h"
#include "fmpq.h"
#include "fmpq_poly.h"

TEST_FUNCTION_START(fmpq_poly_interpolate_fmpz_vec, state)
{
    int i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t P;
        fmpz *x, *y, *z;
        fmpq_t q;
        slong j, n, bits;

        n = n_randint(state, 50);
        bits = n_randint(state, 100);

        x = _fmpz_vec_init(n);
        y = _fmpz_vec_init(n);
        z = _fmpz_vec_init(n);

        fmpq_poly_init(P);

        for (j = 0; j < n; j++)
            fmpz_set_si(x + j, -n/2 + j);

        _fmpz_vec_randtest(y, state, n, bits);

        fmpq_poly_interpolate_fmpz_vec(P, x, y, n);

        fmpq_init(q);
        for (j = 0; j < n; j++)
        {
            fmpq_poly_evaluate_fmpz(q, P, x + j);
            fmpz_set(z + j, fmpq_numref(q));

            if (!fmpz_equal(z + j, y + j) || !fmpz_is_one(fmpq_denref(q)))
            {
                flint_printf("FAIL:\n");
                flint_printf("x:\n"); _fmpz_vec_print(x, n); flint_printf("\n\n");
                flint_printf("y:\n"); _fmpz_vec_print(y, n); flint_printf("\n\n");
                flint_printf("P:\n"); fmpq_poly_print(P), flint_printf("\n\n");
                fflush(stdout);
                flint_abort();
            }
        }
        fmpq_clear(q);

        fmpq_poly_clear(P);
        _fmpz_vec_clear(x, n);
        _fmpz_vec_clear(y, n);
        _fmpz_vec_clear(z, n);
    }

    TEST_FUNCTION_END(state);
}
