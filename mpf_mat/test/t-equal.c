/*
    Copyright (C) 2010 William Hart
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
#include "mpf_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);

    flint_printf("equal....");
    fflush(stdout);



    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        mpf_mat_t A, B, D, E;
        slong m, n, j;

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        mpf_mat_init(A, m, n, 200);
        mpf_mat_init(B, m, n, 200);
        mpf_mat_init(D, m + 1, n, 200);
        mpf_mat_init(E, m, n + 1, 200);

        if (mpf_mat_equal(A, D) || mpf_mat_equal(A, E))
        {
            flint_printf("FAIL: different dimensions should not be equal\n");
            abort();
        }

        mpf_mat_randtest(A, state, 200);
        mpf_mat_set(B, A);

        if (!mpf_mat_equal(A, B))
        {
            flint_printf("FAIL: copied matrices should be equal\n");
            abort();
        }

        if (m && n)
        {
            j = n_randint(state, m * n);
            mpf_add_ui(A->entries + j, A->entries + j, 1);

            if (mpf_mat_equal(A, B))
            {
                flint_printf("FAIL: modified matrices should not be equal\n");
                abort();
            }
        }

        mpf_mat_clear(A);
        mpf_mat_clear(B);
        mpf_mat_clear(D);
        mpf_mat_clear(E);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
