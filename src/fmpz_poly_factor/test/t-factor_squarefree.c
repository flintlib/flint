/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_poly_factor.h"

TEST_FUNCTION_START(fmpz_poly_factor_squarefree, state)
{
    int i, result;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f, g[10], h, t;
        fmpz_poly_factor_t fac;
        slong k, l, n = n_randint(state, 10) + 1;

        fmpz_poly_init(f);
        fmpz_poly_init(h);
        fmpz_poly_init(t);
        fmpz_poly_factor_init(fac);

        fmpz_poly_one(f);
        for (k = 0; k < n; k++)
        {
            fmpz_poly_init(g[k]);
            fmpz_poly_randtest_not_zero(g[k], state, n_randint(state, 40)+1, 20);
            l = n_randint(state, 2) + 1;
            while (l--)
                fmpz_poly_mul(f, f, g[k]);
        }
        fmpz_poly_factor_squarefree(fac, f);

        /* Squarefree? */
        result = 1;
        for (k = 0; k < fac->num && result; k++)
        {
            fmpz_poly_derivative(h, fac->p + k);
            fmpz_poly_gcd(t, h, fac->p + k);
            result &= fmpz_poly_is_one(t);
        }

        /* Product? */
        fmpz_poly_set_fmpz(h, &(fac->c));
        for (k = 0; k < fac->num; k++)
        {
            if (fac->exp[k] == 1)
                fmpz_poly_mul(h, h, fac->p + k);
            else
            {
                fmpz_poly_pow(t, fac->p + k, fac->exp[k]);
                fmpz_poly_mul(h, h, t);
            }
        }

        result &= fmpz_poly_equal(f, h);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("f = "), fmpz_poly_print_pretty(f, "x"), flint_printf("\n\n");
            for (k = 0; k < n; k++)
            {
                flint_printf("g[%wd] = ", k), fmpz_poly_print_pretty(g[k], "x"), flint_printf("\n\n");
            }
            flint_printf("h = "), fmpz_poly_print_pretty(h, "x"), flint_printf("\n\n");
            fmpz_poly_factor_print(fac);
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(f);
        for (k = 0; k < n; k++)
            fmpz_poly_clear(g[k]);
        fmpz_poly_clear(h);
        fmpz_poly_clear(t);
        fmpz_poly_factor_clear(fac);
    }

    /* check that the zero polynomial doesn't crash */
    {
        fmpz_poly_factor_t fac;
        fmpz_poly_t f;

        fmpz_poly_init(f);
        fmpz_poly_factor_init(fac);

        fmpz_poly_factor_squarefree(fac, f);

        if (fac->num != 0 || !fmpz_is_zero(&(fac->c)))
        {
            flint_printf("FAIL: 0\n");
            flint_abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_factor_clear(fac);
    }

    TEST_FUNCTION_END(state);
}
