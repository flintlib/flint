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
#include "profiler.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpq_poly.h"


int main()
{
    fmpq_poly_t Pn, Pn1, Pn2, R;

    slong n;

    FLINT_TEST_INIT(state);

    flint_printf("legendre_polynomial....");
    fflush(stdout);

    fmpq_poly_init(Pn);
    fmpq_poly_init(Pn1);
    fmpq_poly_init(Pn2);
    fmpq_poly_init(R);

    fmpq_poly_set_ui(Pn, UWORD(1));
    fmpq_poly_set_coeff_ui(Pn1, 1, UWORD(1));

    for (n = 0; n <= 500; n++)
    {
        fmpq_poly_legendre_p(R, n);

        if (!fmpq_poly_equal(Pn, R))
        {
            flint_printf("FAIL: n = %wd\n", n);
            flint_printf("Direct: "); fmpq_poly_print_pretty(R, "x"); flint_printf("\n");
            flint_printf("Recur.: "); fmpq_poly_print_pretty(Pn, "x"); flint_printf("\n");
            fflush(stdout);
            flint_abort();
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

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
