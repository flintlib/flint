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
#include "flint.h"
#include "arith.h"
#include "fmpz.h"
#include "fmpz_poly.h"


int main()
{
    fmpz_poly_t T, U;

    slong n;

    FLINT_TEST_INIT(state);

    flint_printf("chebyshev_u_polynomial....");
    fflush(stdout);

    fmpz_poly_init(T);
    fmpz_poly_init(U);

    for (n = 0; n <= 500; n++)
    {
        arith_chebyshev_u_polynomial(U, n);
        arith_chebyshev_t_polynomial(T, n + 1);
        fmpz_poly_derivative(T, T);
        fmpz_poly_scalar_divexact_ui(T, T, n + 1);

        if (!fmpz_poly_equal(T, U))
        {
            flint_printf("FAIL: n = %wd\n", n);
            flint_printf("T: "); fmpz_poly_print_pretty(T, "x"); flint_printf("\n");
            flint_printf("U: "); fmpz_poly_print_pretty(U, "x"); flint_printf("\n");
            abort();
        }

    }

    fmpz_poly_clear(T);
    fmpz_poly_clear(U);

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
