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
#include "flint.h"
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"


int
main(void)
{
    slong i;

    FLINT_TEST_INIT(state);

    flint_printf("det_interpolate....");
    fflush(stdout);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_mat_t A;
        fmpz_poly_t a, b;
        slong n, bits, deg;

        n = n_randint(state, 10);
        deg = 1 + n_randint(state, 5);
        bits = 1 + n_randint(state, 100);

        fmpz_poly_mat_init(A, n, n);

        fmpz_poly_init(a);
        fmpz_poly_init(b);

        fmpz_poly_mat_randtest(A, state, deg, bits);

        fmpz_poly_mat_det(a, A);
        fmpz_poly_mat_det_interpolate(b, A);

        if (!fmpz_poly_equal(a, b))
        {
            flint_printf("FAIL:\n");
            flint_printf("determinants don't agree!\n");
            flint_printf("A:\n");
            fmpz_poly_mat_print(A, "x");
            flint_printf("det(A):\n");
            fmpz_poly_print_pretty(a, "x");
            flint_printf("\ndet_interpolate(A):\n");
            fmpz_poly_print_pretty(b, "x");
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);

        fmpz_poly_mat_clear(A);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
