/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "arith.h"

TEST_FUNCTION_START(arith_stirling, state)
{
    fmpz_mat_t mat, mat2, mat3;
    fmpz * row;
    fmpz_t s;

    slong n, k, mm, nn;

    const slong maxn = 40;


    fmpz_init(s);

    /* test a few large Stirling numbers */
    {
        slong iter;
        fmpz_t a, b, c, d;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(d);

        for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
        {
            n = n_randint(state, 200);
            k = 1 + n_randint(state, 200);

            arith_stirling_number_1u(a, n + 1, k);
            arith_stirling_number_1u(b, n, k);
            arith_stirling_number_1u(c, n, k - 1);

            fmpz_set(d, c);
            fmpz_addmul_ui(d, b, n);

            if (!fmpz_equal(a, d))
            {
                flint_printf("stirling1u");
                flint_printf("n = %wd, k = %wd\n", n, k);
                flint_printf("a: "); fmpz_print(a); flint_printf("\n");
                flint_printf("b: "); fmpz_print(b); flint_printf("\n");
                flint_printf("c: "); fmpz_print(c); flint_printf("\n");
                flint_printf("d: "); fmpz_print(d); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
        {
            n = n_randint(state, 1000);
            k = 1 + n_randint(state, 1000);

            arith_stirling_number_2(a, n + 1, k);
            arith_stirling_number_2(b, n, k);
            arith_stirling_number_2(c, n, k - 1);

            fmpz_set(d, c);
            fmpz_addmul_ui(d, b, k);

            if (!fmpz_equal(a, d))
            {
                flint_printf("stirling2");
                flint_printf("n = %wd, k = %wd\n", n, k);
                flint_printf("a: "); fmpz_print(a); flint_printf("\n");
                flint_printf("b: "); fmpz_print(b); flint_printf("\n");
                flint_printf("c: "); fmpz_print(c); flint_printf("\n");
                flint_printf("d: "); fmpz_print(d); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
        {
            n = n_randint(state, 200);
            k = 1 + n_randint(state, 200);

            row = _fmpz_vec_init(k + 1);

            arith_stirling_number_1u_vec(row, n, k + 1);
            arith_stirling_number_1u(b, n, k);

            if (!fmpz_equal(row + k, b))
            {
                flint_printf("stirling1u (2)");
                flint_printf("n = %wd, k = %wd\n", n, k);
                flint_printf("a: "); fmpz_print(row + k); flint_printf("\n");
                flint_printf("b: "); fmpz_print(b); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            _fmpz_vec_clear(row, k + 1);
        }

        for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
        {
            n = n_randint(state, 400);
            k = 1 + n_randint(state, 400);

            row = _fmpz_vec_init(k + 1);

            arith_stirling_number_2_vec(row, n, k + 1);
            arith_stirling_number_2(b, n, k);

            if (!fmpz_equal(row + k, b))
            {
                flint_printf("stirling2 (2)");
                flint_printf("n = %wd, k = %wd\n", n, k);
                flint_printf("a: "); fmpz_print(row + k); flint_printf("\n");
                flint_printf("b: "); fmpz_print(b); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            _fmpz_vec_clear(row, k + 1);
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(d);
    }

    for (mm = 0; mm < maxn / 2; mm++)
    {
        /* Consistency test for stirling1u */

        for (nn = 0; nn < maxn; nn++)
        {
            fmpz_mat_init(mat, mm, nn);
            arith_stirling_matrix_1u(mat);

            for (n = 0; n < mm; n++)
            {
                for (k = 0; k < nn; k++)
                {
                    row = _fmpz_vec_init(k);
                    arith_stirling_number_1u_vec(row, n, k);
                    if (!_fmpz_vec_equal(row, mat->rows[n], k))
                    {
                        flint_printf("stirling1u mat != vec ");
                        flint_printf("nn,n,k=%wd,%wd,%wd\n", nn, n, k);
                        flint_printf("mat: ");
                        _fmpz_vec_print(mat->rows[n], k);
                        flint_printf("\nvec: ");
                        _fmpz_vec_print(row, k);
                        fflush(stdout);
                        flint_abort();
                    }
                    _fmpz_vec_clear(row, k);

                    arith_stirling_number_1u(s, n, k);
                    if (!fmpz_equal(mat->rows[n]+k, s))
                    {
                        flint_printf("stirling1u mat != single ");
                        flint_printf("nn,n,k=%wd,%wd,%wd\n", nn, n, k);
                        flint_printf("mat: ");
                        fmpz_print(mat->rows[n]+k);
                        flint_printf("\nsingle: ");
                        fmpz_print(s);
                        fflush(stdout);
                        flint_abort();
                    }
                }
            }

            fmpz_mat_clear(mat);
        }

        /* Consistency test for stirling1 */
        for (nn = 0; nn < maxn; nn++)
        {
            fmpz_mat_init(mat, mm, nn);
            arith_stirling_matrix_1(mat);

            for (n = 0; n < mm; n++)
            {
                for (k = 0; k < nn; k++)
                {
                    row = _fmpz_vec_init(k);
                    arith_stirling_number_1_vec(row, n, k);
                    if (!_fmpz_vec_equal(row, mat->rows[n], k))
                    {
                        flint_printf("stirling1 mat != vec ");
                        flint_printf("nn,n,k=%wd,%wd,%wd\n", nn, n, k);
                        flint_printf("mat: ");
                        _fmpz_vec_print(mat->rows[n], k);
                        flint_printf("\nvec: ");
                        _fmpz_vec_print(row, k);
                        fflush(stdout);
                        flint_abort();
                    }
                    _fmpz_vec_clear(row, k);

                    arith_stirling_number_1(s, n, k);
                    if (!fmpz_equal(mat->rows[n]+k, s))
                    {
                        flint_printf("stirling1 mat != single ");
                        flint_printf("nn,n,k=%wd,%wd,%wd\n", nn, n, k);
                        flint_printf("mat: ");
                        fmpz_print(mat->rows[n]+k);
                        flint_printf("\nsingle: ");
                        fmpz_print(s);
                        fflush(stdout);
                        flint_abort();
                    }
                }
            }

            fmpz_mat_clear(mat);
        }

        /* Consistency test for stirling2 */
        for (nn = 0; nn < maxn; nn++)
        {
            fmpz_mat_init(mat, mm, nn);
            arith_stirling_matrix_2(mat);

            for (n = 0; n < mm; n++)
            {
                for (k = 0; k < nn; k++)
                {
                    row = _fmpz_vec_init(k);
                    arith_stirling_number_2_vec(row, n, k);
                    if (!_fmpz_vec_equal(row, mat->rows[n], k))
                    {
                        flint_printf("stirling2 mat != vec ");
                        flint_printf("nn,n,k=%wd,%wd,%wd\n", nn, n, k);
                        flint_printf("mat: ");
                        _fmpz_vec_print(mat->rows[n], k);
                        flint_printf("\nvec: ");
                        _fmpz_vec_print(row, k);
                        fflush(stdout);
                        flint_abort();
                    }
                    _fmpz_vec_clear(row, k);

                    arith_stirling_number_2(s, n, k);
                    if (!fmpz_equal(mat->rows[n]+k, s))
                    {
                        flint_printf("stirling2 mat != single ");
                        flint_printf("nn,n,k=%wd,%wd,%wd\n", nn, n, k);
                        flint_printf("mat: ");
                        fmpz_print(mat->rows[n]+k);
                        flint_printf("\nsingle: ");
                        fmpz_print(s);
                        fflush(stdout);
                        flint_abort();
                    }
                }
            }

            fmpz_mat_clear(mat);
        }

    }

    /* Matrix inverse test */
    for (nn = 1; nn < 50; nn++)
    {
        fmpz_mat_init(mat, nn, nn);
        fmpz_mat_init(mat2, nn, nn);
        fmpz_mat_init(mat3, nn, nn);

        arith_stirling_matrix_1(mat);
        arith_stirling_matrix_2(mat2);
        fmpz_mat_mul(mat3, mat, mat2);

        for (n = 0; n < nn; n++)
        {
            for (k = 0; k < nn; k++)
            {
                if (fmpz_get_ui(mat3->rows[n]+k) != (n == k))
                {
                    flint_printf("not identity matrix: %wd, %wd, %wd\n", nn, n, k);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }
        fmpz_mat_clear(mat);
        fmpz_mat_clear(mat2);
        fmpz_mat_clear(mat3);
    }

    fmpz_clear(s);

    TEST_FUNCTION_END(state);
}
