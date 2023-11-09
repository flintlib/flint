/*
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"

TEST_FUNCTION_START(fmpz_poly_eulerian_polynomial, state)
{
    ulong n, ix, mx;
    fmpz_t sum, fac;
    fmpz_poly_t poly;

    fmpz_poly_init(poly);
    fmpz_init(sum);
    fmpz_init(fac);

    for (ix = 0; ix < 100; ix++)
    {
        n = n_randint(state, 1000) + 1; /* Don't want n = 0. */
        fmpz_poly_eulerian_polynomial(poly, n);

        if (!fmpz_is_one(poly->coeffs))
        {
            flint_printf("The first coefficient is not 1 for n = %u. Received:\n", n);
            fmpz_poly_print_pretty(poly, "x");
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        if (!fmpz_is_one(poly->coeffs + (n - 1)))
        {
            flint_printf("The last coefficient is not 1 for n = %u. Received:\n", n);
            fmpz_poly_print_pretty(poly, "x");
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        for (mx = 0; mx <= n / 2; mx++)
        {
            if (!fmpz_equal(poly->coeffs + mx, poly->coeffs + (n - mx - 1)))
            {
                flint_printf("A(n, m) is not equal to A(n, n - m - 1)"
                        "for n = %u and m = %u.\n", n, mx);
                fflush(stdout);
                flint_abort();
            }
            if (fmpz_cmp_si(poly->coeffs + mx, 0) <= 0)
            {
                flint_printf("Negative coefficients for A(n, m)"
                        "where n = %u and m = %u.\n", n, mx);
                fflush(stdout);
                flint_abort();
            }
        }
        fmpz_zero(sum);
        for (mx = 0; mx < n; mx++)
            fmpz_add(sum, sum, poly->coeffs + mx);
        fmpz_fac_ui(fac, n);
        if (!fmpz_equal(sum, fac))
        {
            flint_printf("The sum of the coefficients of the %u'th polynomial"
                    "was not equal to %u!.\n", n, n);
            flint_printf("Expected: "), fmpz_print(fac), flint_printf("\n");
            flint_printf("Got:      "), fmpz_print(sum), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }
    }

    fmpz_poly_clear(poly);
    fmpz_clear(sum);
    fmpz_clear(fac);

    TEST_FUNCTION_END(state);
}
