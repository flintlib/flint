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
    fmpz_poly_t T0, T1, T2, t;
    slong n;

    FLINT_TEST_INIT(state);

    flint_printf("chebyshev_t_polynomial....");
    fflush(stdout);

    fmpz_poly_init(T0);
    fmpz_poly_init(T1);
    fmpz_poly_init(T2);
    fmpz_poly_init(t);

    arith_chebyshev_t_polynomial(T0, 0);
    arith_chebyshev_t_polynomial(T1, 1);

    for (n = 2; n <= 500; n++)
    {
        arith_chebyshev_t_polynomial(T2, n);

        /* Verify T_{n+1} = 2 x T_n - T_{n-1} */
        fmpz_poly_scalar_mul_ui(t, T1, UWORD(2));
        fmpz_poly_shift_left(t, t, 1);
        fmpz_poly_sub(t, t, T0);

        if (!fmpz_poly_equal(t, T2))
        {
            flint_printf("FAIL: n = %wd\n", n);
            flint_printf("t: "); fmpz_poly_print_pretty(t, "x"); flint_printf("\n");
            flint_printf("T2: "); fmpz_poly_print_pretty(T2, "x"); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_swap(T0, T1);
        fmpz_poly_swap(T1, T2);
    }

    fmpz_poly_clear(T0);
    fmpz_poly_clear(T1);
    fmpz_poly_clear(T2);
    fmpz_poly_clear(t);

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
