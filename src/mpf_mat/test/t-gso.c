/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gmpcompat.h"
#include "mpf_mat.h"

TEST_FUNCTION_START(mpf_mat_gso, state)
{
    int i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        mpf_t dot, tmp;
        int j, k, l;
        mpf_mat_t A, B;
        flint_bitcnt_t prec;

        slong m, n;

        m = n_randint(state, 10);
        n = n_randint(state, 10);
        prec = n_randint(state, 200) + 3;

        mpf_mat_init(A, m, n, prec);
        mpf_mat_init(B, m, n, prec);
        mpf_init2(dot, prec);
        mpf_init2(tmp, prec);

        mpf_mat_randtest(A, state, prec);

        mpf_mat_set(B, A);

        mpf_mat_gso(A, A);

        for (j = 0; j < n; j++)
        {
            mpf_t norm;
            mpf_init2(norm, prec);
            for (l = 0; l < m; l++)
            {
                mpf_mul(tmp, mpf_mat_entry(A, l, j), mpf_mat_entry(A, l, j));
                mpf_add(norm, norm, tmp);
            }

            mpf_sub_ui(dot, norm, 1);
            mpf_abs(dot, dot);
            flint_mpf_set_ui(tmp, 1);
            mpf_div_2exp(tmp, tmp, prec - 3);
            if (flint_mpf_cmp_ui(norm, 0) != 0 && mpf_cmp(dot, tmp) > 0)
            {
                flint_printf("FAIL:\n");
                flint_printf("A:\n");
                mpf_mat_print(A);
                mpf_out_str(stdout, 10, 0, norm);
                flint_printf("\n");
                mpf_out_str(stdout, 10, 0, dot);
                flint_printf("\n");
                mpf_out_str(stdout, 10, 0, tmp);
                flint_printf("\n");
                flint_printf("%d\n", prec);
                flint_printf("%d\n", j);
                fflush(stdout);
                flint_abort();
            }
            mpf_clear(norm);
            for (k = j + 1; k < n; k++)
            {
                flint_mpf_set_ui(dot, 0);
                for (l = 0; l < m; l++)
                {
                    mpf_mul(tmp, mpf_mat_entry(A, l, j),
                            mpf_mat_entry(A, l, k));
                    mpf_add(dot, dot, tmp);
                }

                flint_mpf_set_ui(tmp, 1);
                mpf_div_2exp(tmp, tmp, prec);
                mpf_abs(dot, dot);
                if (mpf_cmp(dot, tmp) > 0)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("A:\n");
                    mpf_mat_print(A);
                    flint_printf("B:\n");
                    mpf_mat_print(B);
                    mpf_out_str(stdout, 10, 0, dot);
                    flint_printf("\n");
                    flint_printf("%d %d\n", j, k);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        mpf_mat_clear(A);
        mpf_mat_clear(B);
        mpf_clears(dot, tmp, NULL);
    }

    TEST_FUNCTION_END(state);
}
