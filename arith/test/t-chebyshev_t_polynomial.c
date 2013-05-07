/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

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
    long n;

    printf("chebyshev_t_polynomial....");
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
        fmpz_poly_scalar_mul_ui(t, T1, 2UL);
        fmpz_poly_shift_left(t, t, 1);
        fmpz_poly_sub(t, t, T0);

        if (!fmpz_poly_equal(t, T2))
        {
            printf("FAIL: n = %ld\n", n);
            printf("t: "); fmpz_poly_print_pretty(t, "x"); printf("\n");
            printf("T2: "); fmpz_poly_print_pretty(T2, "x"); printf("\n");
            abort();
        }

        fmpz_poly_swap(T0, T1);
        fmpz_poly_swap(T1, T2);
    }

    fmpz_poly_clear(T0);
    fmpz_poly_clear(T1);
    fmpz_poly_clear(T2);
    fmpz_poly_clear(t);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
