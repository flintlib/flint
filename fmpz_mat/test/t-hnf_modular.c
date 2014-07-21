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

    Copyright (C) 2014 Alex J. Best

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"

/* checks that the input matrix is in Hermite normal form */
int in_hnf(const fmpz_mat_t A)
{
    slong i, last_i, j, prev_j;

    /* find last non-zero row */
    for (last_i = A->r - 1; last_i >= 0; last_i--)
    {
        for (j = 0; j < A->c; j++)
        {
            if (!fmpz_is_zero(fmpz_mat_entry(A, last_i, j)))
                break;
        }
        if (j < A->c)
            break;
    }

    /* hermite form structure */
    prev_j = -1;
    for (i = 0; i <= last_i; i++)
    {
        slong i2;

        for (j = 0; j < A->c; j++)
        {
            if (!fmpz_is_zero(fmpz_mat_entry(A, i, j)))
            {
                if (fmpz_sgn(fmpz_mat_entry(A, i, j)) < 0)
                    return 0;
                break;
            }
        }
        if (j == A->c || j <= prev_j)
            return 0;
        prev_j = j;
        for (i2 = 0; i2 < i; i2++)
        {
            if (fmpz_cmp(fmpz_mat_entry(A, i2, j), fmpz_mat_entry(A, i, j)) >= 0)
                return 0;
            if (fmpz_sgn(fmpz_mat_entry(A, i2, j)) < 0)
                return 0;
        }
    }

    return 1;
}

int
main(void)
{
    slong iter;
    FLINT_TEST_INIT(state);

    flint_printf("hnf_modular....");
    fflush(stdout);

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        fmpz_t det;
        fmpz_mat_t A, B, H, H2;
        slong m, n, b, c, d, i, j;
        int equal;

        n = n_randint(state, 10);
        m = n + n_randint(state, 10);

        fmpz_init(det);

        fmpz_mat_init(A, n, n);
        fmpz_mat_init(B, m, n);
        fmpz_mat_init(H, m, n);
        fmpz_mat_init(H2, m, n);

        /* sparse */
        b = 1 + n_randint(state, 10) * n_randint(state, 10);
        fmpz_mat_randrank(A, state, n, b);

        fmpz_mat_det(det, A);
        c = 1 + n_randint(state, 10);
        fmpz_abs(det, det);
        fmpz_mul_ui(det, det, c);

        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                fmpz_set(fmpz_mat_entry(B, i, j), fmpz_mat_entry(A, i, j));

        /* dense */
        d = n_randint(state, 2*m*n + 1);
        if (n_randint(state, 2))
            fmpz_mat_randops(B, state, d);

        fmpz_mat_hnf_modular(H, B, det);

        if (!in_hnf(H))
        {
            flint_printf("FAIL:\n");
            flint_printf("matrix not in hnf!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(B); flint_printf("\n\n");
            fmpz_mat_print_pretty(H); flint_printf("\n\n");
            abort();
        }

        fmpz_mat_hnf_xgcd(H2, B);
        equal = fmpz_mat_equal(H, H2);

        if (!equal)
        {
            flint_printf("FAIL:\n");
            flint_printf("hnfs produced by different methods should be the same!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(B); flint_printf("\n\n");
            fmpz_mat_print_pretty(H); flint_printf("\n\n");
            fmpz_mat_print_pretty(H2); flint_printf("\n\n");
            fmpz_print(det); flint_printf("\n\n");
            abort();
        }

        fmpz_mat_hnf_modular(H2, H, det);
        equal = fmpz_mat_equal(H, H2);

        if (!equal)
        {
            flint_printf("FAIL:\n");
            flint_printf("hnf of a matrix in hnf should be the same!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(B); flint_printf("\n\n");
            fmpz_mat_print_pretty(H); flint_printf("\n\n");
            fmpz_mat_print_pretty(H2); flint_printf("\n\n");
            abort();
        }

        fmpz_mat_clear(H2);
        fmpz_mat_clear(H);
        fmpz_mat_clear(B);
        fmpz_mat_clear(A);
        fmpz_clear(det);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}

