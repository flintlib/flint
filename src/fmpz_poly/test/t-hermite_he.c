/*
    Copyright (C) 2016  Ralf Stephan

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_poly.h"

TEST_FUNCTION_START(fmpz_poly_hermite_he, state)
{
    fmpz_poly_t T0, T1, t1, t2;
    slong n;

    fmpz_poly_init(T0);
    fmpz_poly_init(T1);
    fmpz_poly_init(t1);
    fmpz_poly_init(t2);

    fmpz_poly_hermite_he(T0, 0);
    fmpz_poly_hermite_he(t1, 0);

    for (n = 1; n <= 500; n++)
    {
        fmpz_poly_hermite_he(T1, n);

        /* Verify H_{n+1} = 2 x H_n - diff(H_n) */
        fmpz_poly_shift_left(t1, t1, 1);
        fmpz_poly_set(t2, T0);
        fmpz_poly_derivative(t2, t2);
        fmpz_poly_sub(t1, t1, t2);

        if (!fmpz_poly_equal(t1, T1))
        {
            flint_printf("FAIL: n = %wd\n", n);
            flint_printf("t1: ");
            fmpz_poly_print_pretty(t1, "x");
            flint_printf("\n");
            flint_printf("T1: ");
            fmpz_poly_print_pretty(T1, "x");
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_swap(T0, T1);
    }

    fmpz_poly_clear(T0);
    fmpz_poly_clear(T1);
    fmpz_poly_clear(t1);
    fmpz_poly_clear(t2);

    TEST_FUNCTION_END(state);
}
