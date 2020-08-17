/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "d_mat.h"
#include "ulong_extras.h"

#define D_MAT_QR_EQ_EPS (3 * D_EPS)
#define D_MAT_QR_NORM_EPS (4 * D_EPS)
#define D_MAT_QR_ORTHO_EPS (2 * D_EPS)

int
main(void)
{
    int i;
    int tmul = 100;
#ifdef _WIN32
    tmul = 1;
#endif
    FLINT_TEST_INIT(state);

    flint_printf("qr....");
    fflush(stdout);

    /* check QR = A
     * check norm(column(Q)) = 1 or 0
     * check dot product of columns of Q is zero */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        double dot;
        int j, k, l;
        d_mat_t A, Q, R, B;

        slong m, n;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        d_mat_init(A, m, n);
        d_mat_init(Q, m, n);
        d_mat_init(R, n, n);
        d_mat_init(B, m, n);

        d_mat_randtest(A, state, 0, 0);
        d_mat_zero(R);

        d_mat_qr(Q, R, A);

        d_mat_mul_classical(B, Q, R);

        if (!d_mat_approx_equal(A, B, D_MAT_QR_EQ_EPS))
        {
            flint_printf("FAIL:\n");
            flint_printf("A:\n");
            d_mat_print(A);
            flint_printf("Q:\n");
            d_mat_print(Q);
            flint_printf("R:\n");
            d_mat_print(R);
            flint_printf("B:\n");
            d_mat_print(B);
            abort();
        }

        for (j = 0; j < n; j++)
        {
            double norm = 0;
            for (l = 0; l < m; l++)
            {
                norm += d_mat_entry(Q, l, j) * d_mat_entry(Q, l, j);
            }
            if (norm != 0 && fabs(norm - 1) > D_MAT_QR_NORM_EPS)
            {
                flint_printf("FAIL:\n");
                flint_printf("Q:\n");
                d_mat_print(Q);
                flint_printf("%g\n", norm);
                flint_printf("%d\n", j);
                abort();
            }
            for (k = j + 1; k < n; k++)
            {

                dot = 0;
                for (l = 0; l < m; l++)
                {
                    dot += d_mat_entry(Q, l, j) * d_mat_entry(Q, l, k);
                }

                if (fabs(dot) > D_MAT_QR_ORTHO_EPS)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("Q:\n");
                    d_mat_print(Q);
                    flint_printf("%g\n", dot);
                    flint_printf("%d %d\n", j, k);
                    abort();
                }
            }
        }

        d_mat_clear(A);
        d_mat_clear(Q);
        d_mat_clear(R);
        d_mat_clear(B);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
