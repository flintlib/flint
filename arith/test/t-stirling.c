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

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <gmp.h>
#include "flint.h"
#include "arith.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "ulong_extras.h"
#include "profiler.h"

int main(void)
{
    fmpz_mat_t mat, mat2, mat3;

    fmpz * row;
    fmpz_t s;

    slong n, k, mm, nn;

    const slong maxn = 40;

    FLINT_TEST_INIT(state);

    flint_printf("stirling....");
    fflush(stdout);

    fmpz_init(s);

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
                        abort();
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
                        abort();
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
                        abort();
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
                        abort();
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
                        abort();
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
                        abort();
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
                    abort();
                }
            }
        }
        fmpz_mat_clear(mat);
        fmpz_mat_clear(mat2);
        fmpz_mat_clear(mat3);
    }

    fmpz_clear(s);

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
