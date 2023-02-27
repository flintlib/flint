/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
#include "flint.h"
#include "arith.h"
#include "profiler.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpq_poly.h"


int main()
{
    fmpq_poly_t P, Q;
    mpz_t t;

    slong k, n;

    FLINT_TEST_INIT(state);

    flint_printf("euler_polynomial....");
    fflush(stdout);

    for (n = 0; n <= 100; n++)
    {
        fmpq_poly_init(P);
        fmpq_poly_init(Q);

        mpz_init(t);

        for (k = 0; k < n; k++)
        {
            arith_euler_polynomial(P, k);
            flint_mpz_bin_uiui(t, n, k);
            fmpq_poly_scalar_mul_mpz(P, P, t);
            fmpq_poly_add(Q, Q, P);
        }

        fmpq_poly_scalar_div_ui(Q, Q, 2);

        arith_euler_polynomial(P, n);
        fmpq_poly_add(Q, Q, P);

        mpz_clear(t);

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

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
