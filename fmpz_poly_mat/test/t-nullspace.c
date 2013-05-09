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
#include "flint.h"
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"


int
main(void)
{
    flint_rand_t state;
    len_t i;

    printf("nullspace....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_mat_t A, N, AN;
        len_t n, m, bits, deg, rank, nullity;
        float density;

        m = n_randint(state, 13);
        n = n_randint(state, 13);
        deg = 1 + n_randint(state, 5);
        bits = 1 + n_randint(state, 100);
        density = n_randint(state, 100) * 0.01;

        fmpz_poly_mat_init(A, m, n);
        fmpz_poly_mat_init(N, n, n);
        fmpz_poly_mat_init(AN, m, n);

        fmpz_poly_mat_randtest_sparse(A, state, deg, bits, density);

        rank = fmpz_poly_mat_rank(A);
        nullity = fmpz_poly_mat_nullspace(N, A);

        if (nullity + rank != n)
        {
            printf("FAIL: wrong nullity!\n");
            printf("rank = %ld\n", rank);
            printf("nullity = %ld\n", nullity);
            fmpz_poly_mat_print(A, "x");
            printf("\n");
            fmpz_poly_mat_print(N, "x");
            printf("\n");
            abort();
        }

        if (fmpz_poly_mat_rank(N) != nullity)
        {
            printf("FAIL: wrong rank(N) != nullity!\n");
            abort();
        }

        fmpz_poly_mat_mul(AN, A, N);

        if (!fmpz_poly_mat_is_zero(AN))
        {
            printf("FAIL: A * N != 0\n");
            abort();
        }

        fmpz_poly_mat_clear(A);
        fmpz_poly_mat_clear(N);
        fmpz_poly_mat_clear(AN);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
