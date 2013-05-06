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
    Copyright (C) 2012 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    long m, n, rep;
    flint_rand_t state;

    printf("charpoly....");
    fflush(stdout);

    flint_randinit(state);

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        fmpz_mat_t A, B, C, D;
        fmpz_poly_t f, g;

        m = n_randint(state, 4);
        n = m;

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(B, m, n);
        fmpz_mat_init(C, m, m);
        fmpz_mat_init(D, n, n);
        fmpz_poly_init(f);
        fmpz_poly_init(g);

        fmpz_mat_randtest(A, state, 10);
        fmpz_mat_randtest(B, state, 10);

        fmpz_mat_mul(C, A, B);
        fmpz_mat_mul(D, B, A);

        fmpz_mat_charpoly(f, C);
        fmpz_mat_charpoly(g, D);

        if (!fmpz_poly_equal(f, g))
        {
            printf("FAIL: charpoly(AB) != charpoly(BA).\n");
            printf("Matrix A:\n"), fmpz_mat_print(A), printf("\n");
            printf("Matrix B:\n"), fmpz_mat_print(B), printf("\n");
            printf("cp(AB) = "), fmpz_poly_print_pretty(f, "X"), printf("\n");
            printf("cp(BA) = "), fmpz_poly_print_pretty(g, "X"), printf("\n");
            abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_mat_clear(D);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
