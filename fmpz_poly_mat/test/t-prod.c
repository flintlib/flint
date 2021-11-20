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
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"

int
main(void)
{
    slong i;

    FLINT_TEST_INIT(state);

    flint_printf("prod....");
    fflush(stdout);    

    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_poly_mat_t A, B, *V;
        slong m, j, count, bits, deg;
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
            flint_printf("FAIL:\n");
            flint_printf("count = %wd\n", count);
            flint_printf("A:\n");
            fmpz_poly_mat_print(A, "x");
            flint_printf("B:\n");
            fmpz_poly_mat_print(B, "x");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_mat_clear(A);
        fmpz_poly_mat_clear(B);
        for (j = 0; j < count; j++)
            fmpz_poly_mat_clear(V[j]);
        flint_free(V);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
