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
#include "fmpz.h"
#include "fmpz_poly.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("div_root....");
    fflush(stdout);

    flint_randinit(state);

    /* Compare with standard divrem */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t P, Q, D, DQ;
        fmpz_t c;
        long n, b;

        n = n_randint(state, 100);
        b = n_randint(state, 200);

        fmpz_init(c);
        fmpz_poly_init(P);
        fmpz_poly_init(Q);
        fmpz_poly_init(D);
        fmpz_poly_init(DQ);

        fmpz_poly_randtest(P, state, n, b);
        fmpz_randtest(c, state, b);

        fmpz_poly_div_root(Q, P, c);

        fmpz_poly_set_coeff_fmpz(D, 0, c);
        fmpz_poly_neg(D, D);
        fmpz_poly_set_coeff_ui(D, 1, 1UL);

        fmpz_poly_div_basecase(DQ, P, D);

        result = fmpz_poly_equal(Q, DQ);

        if (!result)
        {
            printf("FAIL!\n");
            printf("P:\n"); fmpz_poly_print(P); printf("\n\n");
            printf("Q:\n"); fmpz_poly_print(Q); printf("\n\n");
            printf("D:\n"); fmpz_poly_print(D); printf("\n\n");
            printf("DQ:\n"); fmpz_poly_print(DQ); printf("\n\n");
            abort();
        }

        fmpz_clear(c);
        fmpz_poly_clear(P);
        fmpz_poly_clear(Q);
        fmpz_poly_clear(D);
        fmpz_poly_clear(DQ);
    }

    /* Check aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t P, Q1, Q2;
        fmpz_t c;
        long n, b;

        n = n_randint(state, 100);
        b = n_randint(state, 200);

        fmpz_init(c);
        fmpz_poly_init(P);
        fmpz_poly_init(Q1);
        fmpz_poly_init(Q2);

        fmpz_randtest(c, state, b);
        fmpz_poly_randtest(P, state, n, b);
        fmpz_poly_set(Q2, P);

        fmpz_poly_div_root(Q1, P, c);
        fmpz_poly_div_root(Q2, Q2, c);

        result = fmpz_poly_equal(Q1, Q2);

        if (!result)
        {
            printf("FAIL (aliasing)!\n");
            printf("P:\n"); fmpz_poly_print(P); printf("\n\n");
            printf("Q1:\n"); fmpz_poly_print(Q1); printf("\n\n");
            printf("Q2:\n"); fmpz_poly_print(Q2); printf("\n\n");
            abort();
        }

        fmpz_clear(c);
        fmpz_poly_clear(P);
        fmpz_poly_clear(Q1);
        fmpz_poly_clear(Q2);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
