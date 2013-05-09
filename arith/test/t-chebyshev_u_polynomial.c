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
    fmpz_poly_t T, U;

    len_t n;

    printf("chebyshev_u_polynomial....");
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
            printf("FAIL: n = %ld\n", n);
            printf("T: "); fmpz_poly_print_pretty(T, "x"); printf("\n");
            printf("U: "); fmpz_poly_print_pretty(U, "x"); printf("\n");
            abort();
        }

    }

    fmpz_poly_clear(T);
    fmpz_poly_clear(U);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
