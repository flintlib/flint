/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"

TEST_FUNCTION_START(fmpz_poly_theta_qexp, state)
{
    int i;

    for (i = 0; i < 2000; i++)
    {
        fmpz_poly_t a, b;
        slong e, n;

        fmpz_poly_init(a);
        fmpz_poly_init(b);

        e = n_randint(state, 100) - 50;
        n = n_randint(state, 250);

        fmpz_poly_randtest(a, state, n_randint(state, 250),
            1 + n_randint(state, 100));

        fmpz_poly_theta_qexp(a, e, n);

        fmpz_poly_theta_qexp(b, 1, n + n_randint(state, 10));

        if (n == 0)
        {
            fmpz_poly_zero(b);
        }
        else
        {
            if (e >= 0)
            {
                fmpz_poly_pow_trunc(b, b, e, n);
            }
            else
            {
                fmpz_poly_inv_series(b, b, n);
                fmpz_poly_pow_trunc(b, b, -e, n);
            }
        }

        if (!fmpz_poly_equal(a, b))
        {
            flint_printf("FAIL (powering):\n");
            flint_printf("e = %wd, n = %wd\n\n", e, n);
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
    }

    TEST_FUNCTION_END(state);
}
