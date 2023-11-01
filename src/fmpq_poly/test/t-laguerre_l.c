/*
    Copyright (C) 2016  Ralf Stephan

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_poly.h"

TEST_FUNCTION_START(fmpq_poly_laguerre_l, state)
{
    fmpq_poly_t T0, T1, T2, t, tt;
    slong n;

    fmpq_poly_init(T0);
    fmpq_poly_init(T1);
    fmpq_poly_init(T2);
    fmpq_poly_init(t);
    fmpq_poly_init(tt);

    fmpq_poly_laguerre_l(T0, 0);
    fmpq_poly_laguerre_l(T1, 1);

    for (n = 1; n <= 500; n++)
    {
        fmpq_poly_laguerre_l(T2, n+1);
        fmpq_poly_set(t, T1);

        /* Verify (n+1)P_{n+1} = (2n+1-x) P_n - nP_{n-1} */
        fmpq_poly_shift_left(tt, t, 1);
        fmpq_poly_scalar_mul_ui(t, t, 2*n+1);
        fmpq_poly_sub(t, t, tt);
        fmpq_poly_scalar_mul_ui(tt, T0, n);
        fmpq_poly_sub(t, t, tt);
        fmpq_poly_scalar_mul_ui(tt, T2, n+1);

        if (!fmpq_poly_equal(t, tt))
        {
            flint_printf("FAIL: n = %wd\n", n);
            flint_printf("t: "); fmpq_poly_print_pretty(t, "x"); flint_printf("\n");
            flint_printf("tt: "); fmpq_poly_print_pretty(tt, "x"); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_swap(T0, T1);
        fmpq_poly_swap(T1, T2);
    }

    fmpq_poly_clear(T0);
    fmpq_poly_clear(T1);
    fmpq_poly_clear(T2);
    fmpq_poly_clear(t);
    fmpq_poly_clear(tt);

    TEST_FUNCTION_END(state);
}
