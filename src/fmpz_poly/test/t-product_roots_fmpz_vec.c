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
#include "fmpz_poly.h"

TEST_FUNCTION_START(fmpz_poly_product_roots_fmpz_vec, state)
{
    int i, result;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t P, Q, tmp;
        fmpz * x;
        slong j, n, bits;

        n = n_randint(state, 100);
        bits = n_randint(state, 10);

        x = _fmpz_vec_init(n);
        _fmpz_vec_randtest(x, state, n, bits);

        fmpz_poly_init(P);
        fmpz_poly_init(Q);
        fmpz_poly_init(tmp);

        fmpz_poly_product_roots_fmpz_vec(P, x, n);

        fmpz_poly_set_ui(Q, UWORD(1));
        for (j = 0; j < n; j++)
        {
            fmpz_poly_zero(tmp);
            fmpz_poly_set_coeff_si(tmp, 1, WORD(-1));
            fmpz_poly_set_coeff_fmpz(tmp, 0, x + j);
            fmpz_poly_neg(tmp, tmp);
            fmpz_poly_mul(Q, Q, tmp);
        }

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
        fmpz_poly_clear(tmp);
        _fmpz_vec_clear(x, n);
    }

    TEST_FUNCTION_END(state);
}
