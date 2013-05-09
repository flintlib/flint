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

    len_t k, n;

    printf("euler_polynomial....");
    fflush(stdout);

    for (n = 0; n <= 100; n++)
    {
        fmpq_poly_init(P);
        fmpq_poly_init(Q);

        mpz_init(t);

        for (k = 0; k < n; k++)
        {
            arith_euler_polynomial(P, k);
            mpz_bin_uiui(t, n, k);
            fmpq_poly_scalar_mul_mpz(P, P, t);
            fmpq_poly_add(Q, Q, P);
        }

        fmpq_poly_scalar_div_ui(Q, Q, 2);

        arith_euler_polynomial(P, n);
        fmpq_poly_add(Q, Q, P);

        mpz_clear(t);

        fmpq_poly_zero(P);
        fmpq_poly_set_coeff_ui(P, n, 1UL);

        if (!fmpq_poly_equal(P, Q))
        {
            printf("ERROR: sum up to n = %ld did not add to x^n\n", n);
            printf("Sum: ");
            fmpq_poly_print_pretty(Q, "x");
            printf("\nExpected: ");
            fmpq_poly_print_pretty(P, "x");
            printf("\n");
            abort();
        }

        fmpq_poly_clear(P);
        fmpq_poly_clear(Q);
    }

    mpfr_free_cache();
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
