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

    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#include "mpf_mat.h"
#include "ulong_extras.h"

#define MPF_MAT_QR_EQ_BITS (42)
#define MPF_MAT_QR_NORM_EPS (2e-16)
#define MPF_MAT_QR_ORTHO_EPS (2e-20)

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);

    flint_printf("qr....");
    fflush(stdout);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        mpf_t dot, tmp;
        int j, k, l;
        mpf_mat_t A, Q, R, B;

        slong m, n;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        mpf_mat_init(A, m, n, 200);
        mpf_mat_init(Q, m, n, 200);
        mpf_mat_init(R, n, n, 200);
        mpf_mat_init(B, m, n, 200);
        mpf_inits(dot, tmp, '\0');

        mpf_mat_randtest(A, state, 200);

        mpf_mat_qr(Q, R, A);

        mpf_mat_mul(B, Q, R);

        if (!mpf_mat_approx_equal(A, B, MPF_MAT_QR_EQ_BITS))
        {
            flint_printf("FAIL:\n");
            flint_printf("A:\n");
            mpf_mat_print(A);
            flint_printf("Q:\n");
            mpf_mat_print(Q);
            flint_printf("R:\n");
            mpf_mat_print(R);
            flint_printf("B:\n");
            mpf_mat_print(B);
            abort();
        }

        for (j = 0; j < n; j++)
        {
            double d;
            mpf_t norm;
            mpf_init_set_ui(norm, 0);
            for (l = 0; l < m; l++)
            {
                mpf_mul(tmp, mpf_mat_entry(Q, l, j), mpf_mat_entry(Q, l, j));
                mpf_add(norm, norm, tmp);
            }

            d = mpf_get_d(norm);
            if (mpf_cmp_ui(norm, 0) != 0 && fabs(d - 1) > MPF_MAT_QR_NORM_EPS)
            {
                flint_printf("FAIL:\n");
                flint_printf("Q:\n");
                mpf_mat_print(Q);
                flint_printf("%g\n", norm);
                flint_printf("%d\n", j);
                abort();
            }
            mpf_clear(norm);
            for (k = j + 1; k < n; k++)
            {
                mpf_set_ui(dot, 0);
                for (l = 0; l < m; l++)
                {
                    mpf_mul(tmp, mpf_mat_entry(Q, l, j),
                            mpf_mat_entry(Q, l, k));
                    mpf_add(dot, dot, tmp);
                }

                d = mpf_get_d(dot);
                if (fabs(d) > MPF_MAT_QR_ORTHO_EPS)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("Q:\n");
                    mpf_mat_print(Q);
                    flint_printf("%g\n", dot);
                    flint_printf("%d %d\n", j, k);
                    abort();
                }
            }
        }

        mpf_mat_clear(A);
        mpf_mat_clear(Q);
        mpf_mat_clear(R);
        mpf_mat_clear(B);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
