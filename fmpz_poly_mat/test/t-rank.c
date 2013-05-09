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
#include "fmpz_mat.h"
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"


int
main(void)
{
    flint_rand_t state;
    len_t i;

    printf("rank....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_mat_t A;
        fmpz_mat_t Ax;
        fmpz_t x;
        len_t j, m, n, bits, deg, rank, zrank;
        float density;

        m = n_randint(state, 15);
        n = n_randint(state, 15);
        deg = 1 + n_randint(state, 5);
        bits = 1 + n_randint(state, 100);
        density = n_randint(state, 100) * 0.01;

        fmpz_poly_mat_init(A, m, n);
        fmpz_mat_init(Ax, m, n);
        fmpz_init(x);

        fmpz_poly_mat_randtest_sparse(A, state, deg, bits, density);

        /* Probabilistic rank computation */
        zrank = 0;
        for (j = 0; j < 5; j++)
        {
            len_t r;
            fmpz_randbits(x, state, 15);
            fmpz_poly_mat_evaluate_fmpz(Ax, A, x);
            r = fmpz_mat_rank(Ax);
            zrank = FLINT_MAX(zrank, r);
        }

        rank = fmpz_poly_mat_rank(A);

        if (rank != zrank)
        {
            printf("FAIL:\n");
            printf("A:\n");
            fmpz_poly_mat_print(A, "x");
            printf("Computed rank: %ld (zrank = %ld)\n", rank, zrank);
            abort();
        }

        fmpz_clear(x);
        fmpz_mat_clear(Ax);
        fmpz_poly_mat_clear(A);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
