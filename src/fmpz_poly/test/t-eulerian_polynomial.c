/*
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"

TEST_FUNCTION_START(fmpz_poly_eulerian_polynomial, state)
{
    ulong n, mx;
    fmpz_t sum, fac;
    fmpz_poly_t poly;

    fmpz_poly_init(poly);
    fmpz_init(sum);
    fmpz_init(fac);

    for (n = 0; n < FLINT_MIN(100 * flint_test_multiplier(), 500 + 10 * flint_test_multiplier()); n++)
    {
        fmpz_poly_eulerian_polynomial(poly, n);

        if (!fmpz_is_one(poly->coeffs + 0))
            TEST_FUNCTION_FAIL(
                    "First coefficient is not 1.\n"
                    "n = %wu\n"
                    "Got: %{fmpz}\n",
                    n, poly->coeffs + 0);

        if (n == UWORD(0))
            continue;

        if (!fmpz_is_one(poly->coeffs + (n - 1)))
            TEST_FUNCTION_FAIL(
                    "Last coefficient is not 1.\n"
                    "n = %wu\n"
                    "Got: %{fmpz}\n",
                    n, poly->coeffs + (n - 1));

        for (mx = 0; mx <= n / 2; mx++)
        {
            if (!fmpz_equal(poly->coeffs + mx, poly->coeffs + (n - mx - 1)))
                TEST_FUNCTION_FAIL(
                        "A(n, m) != A(n, n - m - 1)\n"
                        "n = %wu\n"
                        "m = %wu\n"
                        "A(n, m) = %{fmpz}\n",
                        "A(n, n - m - 1) = %{fmpz}\n",
                        n, mx, poly->coeffs + mx, poly->coeffs + (n - mx - 1));

            if (fmpz_cmp_si(poly->coeffs + mx, 0) <= 0)
                TEST_FUNCTION_FAIL(
                        "Negative coefficient(s)\n"
                        "n = %wu\n",
                        n);
        }

        fmpz_zero(sum);
        for (mx = 0; mx < n; mx++)
            fmpz_add(sum, sum, poly->coeffs + mx);

        fmpz_fac_ui(fac, n);
        if (!fmpz_equal(sum, fac))
            TEST_FUNCTION_FAIL(
                    "Sum of coefficients is not equal to n-factorial\n"
                    "n = %wu\n"
                    "%wu! = %{fmpz}\n"
                    "sum = %{fmpz}\n",
                    n, fac, sum);
    }

    fmpz_poly_clear(poly);
    fmpz_clear(sum);
    fmpz_clear(fac);

    TEST_FUNCTION_END(state);
}
