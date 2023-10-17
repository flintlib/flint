/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_poly.h"

TEST_FUNCTION_START(fmpq_poly_exp_expinv_series, state)
{
    int i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, c, d;
        slong n = n_randint(state, 50) + 1;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(c);
        fmpq_poly_init(d);

        fmpq_poly_randtest_not_zero(a, state, n_randint(state, 50) + 1, 10);
        fmpq_poly_set_coeff_ui(a, 0, UWORD(0));

        fmpq_poly_canonicalise(a);

        if (n_randint(state, 2))
        {
            fmpq_poly_set(c, a);
            fmpq_poly_exp_expinv_series(b, c, c, n);
        }
        else if (n_randint(state, 2))
        {
            fmpq_poly_set(b, a);
            fmpq_poly_exp_expinv_series(b, c, b, n);
        }
        else
        {
            fmpq_poly_exp_expinv_series(b, c, a, n);
        }

        fmpq_poly_mullow(d, b, c, n);

        if (!fmpq_poly_is_one(d) || !fmpq_poly_is_canonical(b) || !fmpq_poly_is_canonical(c))
        {
            flint_printf("FAIL:\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            fmpq_poly_debug(c), flint_printf("\n\n");
            fmpq_poly_debug(d), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(c);
        fmpq_poly_clear(d);
    }

    TEST_FUNCTION_END(state);
}
