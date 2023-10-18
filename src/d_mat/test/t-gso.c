/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "double_extras.h"
#include "d_mat.h"
#include "ulong_extras.h"

#define D_MAT_GSO_NORM_EPS (4 * D_EPS)
#define D_MAT_GSO_ORTHO_EPS (2 * D_EPS)

TEST_FUNCTION_START(d_mat_gso, state)
{
    int i, tmul = 100;
#ifdef _WIN32
    tmul = 1;
#endif

    /* check norm(column(gso)) = 1 or 0
     * check dot product of columns of gso is zero */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        double dot;
        int j, k, l;
        d_mat_t A;

        slong m, n;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        d_mat_init(A, m, n);

        d_mat_randtest(A, state, 0, 0);

        d_mat_gso(A, A);

        for (j = 0; j < n; j++)
        {
            double norm = 0;
            for (l = 0; l < m; l++)
            {
                norm += d_mat_entry(A, l, j) * d_mat_entry(A, l, j);
            }
            if (norm != 0 && fabs(norm - 1) > D_MAT_GSO_NORM_EPS)
            {
                flint_printf("FAIL:\n");
                flint_printf("A:\n");
                d_mat_print(A);
                flint_printf("%g\n", norm);
                flint_printf("%d\n", j);
                fflush(stdout);
                flint_abort();
            }
            for (k = j + 1; k < n; k++)
            {

                dot = 0;
                for (l = 0; l < m; l++)
                {
                    dot += d_mat_entry(A, l, j) * d_mat_entry(A, l, k);
                }

                if (fabs(dot) > D_MAT_GSO_ORTHO_EPS)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("A:\n");
                    d_mat_print(A);
                    flint_printf("%g\n", dot);
                    flint_printf("%d %d\n", j, k);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        d_mat_clear(A);
    }

    TEST_FUNCTION_END(state);
}
