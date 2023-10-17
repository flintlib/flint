/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpq_poly.h"

TEST_FUNCTION_START(fmpq_poly_revert_series_lagrange_fast, state)
{
    int i, result;

    /* Check aliasing */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g;
        slong n;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        do {
            fmpq_poly_randtest(g, state, n_randint(state, 50), 1+n_randint(state,100));
        } while (fmpq_poly_length(g) < 2 || fmpz_is_zero(g->coeffs + 1));
        fmpq_poly_set_coeff_ui(g, 0, 0);
        n = n_randint(state, 50);

        fmpq_poly_revert_series_lagrange_fast(f, g, n);
        fmpq_poly_revert_series_lagrange_fast(g, g, n);

        result = (fmpq_poly_equal(f, g));
        if (!result)
        {
            flint_printf("FAIL (aliasing):\n");
            fmpq_poly_print(f), flint_printf("\n\n");
            fmpq_poly_print(g), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
    }

    /* Check f(f^(-1)) = id */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g, h;
        slong n;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_init(h);
        do {
            if (n_randint(state, 20) == 0)
                fmpq_poly_randtest(g, state,
                    n_randint(state, 50), 1);
            else
                fmpq_poly_randtest(g, state,
                    n_randint(state, 50), 1+n_randint(state,100));
        } while (fmpq_poly_length(g) < 2 || fmpz_is_zero(g->coeffs + 1));
        fmpq_poly_set_coeff_ui(g, 0, 0);
        n = n_randint(state, 50);

        fmpq_poly_revert_series_lagrange_fast(f, g, n);
        fmpq_poly_compose_series(h, g, f, n);

        result = ((n <= 1 && fmpq_poly_is_zero(h)) ||
            (h->length == 2 && fmpz_is_zero(h->coeffs + 0) &&
                fmpz_is_one(h->coeffs + 1)));
        if (!result)
        {
            flint_printf("FAIL (comparison):\n");
            fmpq_poly_print(f), flint_printf("\n\n");
            fmpq_poly_print(g), flint_printf("\n\n");
            fmpq_poly_print(h), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_poly_clear(h);
    }

    TEST_FUNCTION_END(state);
}
