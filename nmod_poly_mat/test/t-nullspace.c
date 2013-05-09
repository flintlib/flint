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
#include "nmod_poly.h"
#include "nmod_poly_mat.h"
#include "fmpz.h"

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
        nmod_poly_mat_t A, N, AN;
        len_t n, m, deg, rank, nullity;
        float density;
        mp_limb_t mod;

        mod = n_randtest_prime(state, 0);
        m = n_randint(state, 13);
        n = n_randint(state, 13);
        deg = 1 + n_randint(state, 5);
        density = n_randint(state, 100) * 0.01;

        nmod_poly_mat_init(A, m, n, mod);
        nmod_poly_mat_init(N, n, n, mod);
        nmod_poly_mat_init(AN, m, n, mod);

        nmod_poly_mat_randtest_sparse(A, state, deg, density);

        rank = nmod_poly_mat_rank(A);
        nullity = nmod_poly_mat_nullspace(N, A);

        if (nullity + rank != n)
        {
            printf("FAIL: wrong nullity!\n");
            printf("rank = %ld\n", rank);
            printf("nullity = %ld\n", nullity);
            nmod_poly_mat_print(A, "x");
            printf("\n");
            nmod_poly_mat_print(N, "x");
            printf("\n");
            abort();
        }

        if (nmod_poly_mat_rank(N) != nullity)
        {
            printf("FAIL: wrong rank(N) != nullity!\n");
            abort();
        }

        nmod_poly_mat_mul(AN, A, N);

        if (!nmod_poly_mat_is_zero(AN))
        {
            printf("FAIL: A * N != 0\n");
            abort();
        }

        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(N);
        nmod_poly_mat_clear(AN);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
