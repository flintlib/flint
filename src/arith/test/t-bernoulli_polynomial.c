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
#include "fmpq_poly.h"
#include "arith.h"

TEST_FUNCTION_START(arith_bernoulli_polynomial, state)
{
    fmpq_poly_t P, Q;
    fmpz_t t;

    slong k, n;


    for (n = 0; n <= 100; n++)
    {
        fmpq_poly_init(P);
        fmpq_poly_init(Q);

        fmpz_init(t);

        for (k = 0; k <= n; k++)
        {
            arith_bernoulli_polynomial(P, k);
            fmpz_bin_uiui(t, n+1, k);
            fmpq_poly_scalar_mul_fmpz(P, P, t);
            fmpq_poly_add(Q, Q, P);
        }

        fmpq_poly_scalar_div_ui(Q, Q, n+1);
        fmpz_clear(t);

        fmpq_poly_zero(P);
        fmpq_poly_set_coeff_ui(P, n, UWORD(1));

        if (!fmpq_poly_equal(P, Q))
        {
            flint_printf("ERROR: sum up to n = %wd did not add to x^n\n", n);
            flint_printf("Sum: ");
            fmpq_poly_print_pretty(Q, "x");
            flint_printf("\nExpected: ");
            fmpq_poly_print_pretty(P, "x");
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(P);
        fmpq_poly_clear(Q);
    }

    TEST_FUNCTION_END(state);
}
