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
    slong i, j, prev_f;
    fmpz_t zero;

    fmpz_init(zero);
    fmpz_zero(zero);

    for (j = 0; j < A->c; j++)
    {
        for (i = 0; i < A->r; i++)
        {
            if (!fmpz_is_zero(fmpz_mat_entry(A, i, j)))
            {
                break;
            }
        }
        if (i < A->r)
        {
            if (!fmpz_is_zero(fmpz_mat_entry(A, i, j)))
            {
                break;
            }
        }
    }

    prev_f = -1;
    for (; j < A->c; j++)
    {
        fmpz_t last_nonzero;
        fmpz_init(last_nonzero);
        for (i = A->r - 1; i >= prev_f + 1; i--)
        {
            if (!fmpz_is_zero(fmpz_mat_entry(A, i, j)))
            {
                if (fmpz_cmp(fmpz_mat_entry(A, i, j), zero) < 0)
                    return 0;
                fmpz_set(last_nonzero, fmpz_mat_entry(A, i, j));
                break;
            }
        }
        if (i == prev_f)
        {
            return 0;
        }
        prev_f = i;
        i--;
        for (; i >= 0; i--)
        {
            if (fmpz_cmp(fmpz_mat_entry(A, i, j), last_nonzero) >= 0)
                return 0;
            if (fmpz_cmp(fmpz_mat_entry(A, i, j), zero) < 0)
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

    flint_printf("hnf....");
    fflush(stdout);

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        fmpz_mat_t A, H, H2;
        fmpz_t c;
        slong m, n, b, d, r;
        int equal;

        m = n_randint(state, 7);
        n = n_randint(state, 7);
        r = n_randint(state, FLINT_MIN(m, n) + 1);

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(H, m, n);
        fmpz_mat_init(H2, m, n);

        fmpz_init(c);

        /* sparse */
        b = 1 + n_randint(state, 10) * n_randint(state, 10);
        d = n_randint(state, 2*m*n + 1);
        fmpz_mat_randrank(A, state, r, b);

        /* dense */
        if (n_randint(state, 2))
            fmpz_mat_randops(A, state, d);

        fmpz_mat_hnf_classical(H, A);

        if (!in_hnf(H))
        {
            flint_printf("FAIL:\n");
            flint_printf("matrix not in hnf!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(H); flint_printf("\n\n");
            abort();
        }

        fmpz_randtest_not_zero(c, state, 5);

        fmpz_mat_hnf_classical(H2, H);
        equal = fmpz_mat_equal(H,H2);

        if (!equal)
        {
            flint_printf("FAIL:\n");
            flint_printf("hnf of a matrix in hnf should be the same!\n");
            fmpz_mat_print_pretty(H); flint_printf("\n\n");
            fmpz_mat_print_pretty(H2); flint_printf("\n\n");
            abort();
        }

        fmpz_clear(c);

        fmpz_mat_clear(A);
        fmpz_mat_clear(H);
        fmpz_mat_clear(H2);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}

