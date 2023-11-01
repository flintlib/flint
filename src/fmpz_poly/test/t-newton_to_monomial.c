/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

TEST_FUNCTION_START(fmpz_poly_newton_to_monomial, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f, g;
        fmpz * r;
        slong k, n;

        fmpz_poly_init(f);
        fmpz_poly_init(g);

        fmpz_poly_randtest(f, state, 1 + n_randint(state, 20),
                                     1 + n_randint(state, 200));

        n = fmpz_poly_length(f);
        r = _fmpz_vec_init(n);

        for (k = 0; k < n; k++)
            fmpz_randtest(r + k, state, n_randint(state, 200));

        fmpz_poly_set(g, f);

        _fmpz_poly_newton_to_monomial(g->coeffs, r, n);
        _fmpz_poly_monomial_to_newton(g->coeffs, r, n);

        if (!fmpz_poly_equal(f, g))
        {
            flint_printf("FAIL: roundtrip\n");
            fmpz_poly_print(f); flint_printf("\n");
            fmpz_poly_print(g); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        _fmpz_vec_clear(r, n);
    }

    TEST_FUNCTION_END(state);
}
