/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

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
#include "d_mat.h"
#include "ulong_extras.h"

#define D_MAT_MUL_CLASSICAL_EPS (1e-11)

int
main(void)
{
    d_mat_t A, B, C, D, E, F, G;
    slong i;
    FLINT_TEST_INIT(state);

    flint_printf("mul....");
    fflush(stdout);

    /* check associative law */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        slong m, n, k, l;

        m = n_randint(state, 50);
        n = n_randint(state, 50);
        k = n_randint(state, 50);
        l = n_randint(state, 50);

        d_mat_init(A, m, n);
        d_mat_init(B, n, k);
        d_mat_init(C, k, l);
        d_mat_init(D, n, l);
        d_mat_init(E, m, k);
        d_mat_init(F, m, l);
        d_mat_init(G, m, l);

        d_mat_randtest(A, state, 0, 0);
        d_mat_randtest(B, state, 0, 0);
        d_mat_randtest(C, state, 0, 0);

        d_mat_mul_classical(D, B, C);
        d_mat_mul_classical(E, A, B);
        d_mat_mul_classical(F, A, D);
        d_mat_mul_classical(G, E, C);

        if (!d_mat_approx_equal(F, G, D_MAT_MUL_CLASSICAL_EPS))
        {
            flint_printf("FAIL: results not equal\n");
            d_mat_print(F);
            d_mat_print(G);
            abort();
        }

        d_mat_mul_classical(A, A, B);

        if (!d_mat_equal(A, E))
        {
            flint_printf("FAIL: aliasing failed\n");
            abort();
        }

        d_mat_clear(A);
        d_mat_clear(B);
        d_mat_clear(C);
        d_mat_clear(D);
        d_mat_clear(E);
        d_mat_clear(F);
        d_mat_clear(G);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
