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
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"

int
main(void)
{
    flint_rand_t state;
    len_t i;

    printf("prod....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_poly_mat_t A, B, *V;
        len_t m, j, count, bits, deg;
        float density;

        m = n_randint(state, 6);
        deg = 1 + n_randint(state, 6);
        bits = 1 + n_randint(state, 100);
        count = n_randint(state, 20);
        density = n_randint(state, 100) * 0.01;

        fmpz_poly_mat_init(A, m, m);
        fmpz_poly_mat_init(B, m, m);

        V = flint_malloc(sizeof(fmpz_poly_mat_t) * count);
        for (j = 0; j < count; j++)
        {
            fmpz_poly_mat_init(V[j], m, m);
            fmpz_poly_mat_randtest_sparse(V[j], state, deg, bits, density);
        }

        fmpz_poly_mat_prod(A, V, count);

        fmpz_poly_mat_one(B);
        for (j = 0; j < count; j++)
            fmpz_poly_mat_mul(B, B, V[j]);

        if (!fmpz_poly_mat_equal(A, B))
        {
            printf("FAIL:\n");
            printf("count = %ld\n", count);
            printf("A:\n");
            fmpz_poly_mat_print(A, "x");
            printf("B:\n");
            fmpz_poly_mat_print(B, "x");
            abort();
        }

        fmpz_poly_mat_clear(A);
        fmpz_poly_mat_clear(B);
        for (j = 0; j < count; j++)
            fmpz_poly_mat_clear(V[j]);
        flint_free(V);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
