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
#include "profiler.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpq_poly.h"


int main()
{
    fmpq_poly_t Pn, Pn1, Pn2, R;

    long n;

    printf("legendre_polynomial....");
    fflush(stdout);

    fmpq_poly_init(Pn);
    fmpq_poly_init(Pn1);
    fmpq_poly_init(Pn2);
    fmpq_poly_init(R);

    fmpq_poly_set_ui(Pn, 1UL);
    fmpq_poly_set_coeff_ui(Pn1, 1, 1UL);

    for (n = 0; n <= 500; n++)
    {
        arith_legendre_polynomial(R, n);

        if (!fmpq_poly_equal(Pn, R))
        {
            printf("FAIL: n = %ld\n", n);
            printf("Direct: "); fmpq_poly_print_pretty(R, "x"); printf("\n");
            printf("Recur.: "); fmpq_poly_print_pretty(Pn, "x"); printf("\n");
            abort();
        }

        fmpq_poly_shift_left(Pn2, Pn1, 1);
        fmpq_poly_scalar_mul_ui(Pn2, Pn2, 2*n + 3);
        fmpq_poly_scalar_mul_si(Pn, Pn, -(n+1));
        fmpq_poly_add(Pn2, Pn2, Pn);
        fmpq_poly_scalar_div_ui(Pn2, Pn2, n+2);

        fmpq_poly_swap(Pn, Pn1);
        fmpq_poly_swap(Pn1, Pn2);
    }

    fmpq_poly_clear(Pn);
    fmpq_poly_clear(Pn1);
    fmpq_poly_clear(Pn2);
    fmpq_poly_clear(R);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
