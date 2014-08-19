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

/* checks that the input matrix is in Smith normal form */
int in_snf(const fmpz_mat_t A)
{
    slong i, j;
    int snf = 1;
    for (i = 0; i < A->r && snf; i++)
    {
        for (j = 0; j < A->c && snf; j++)
        {
            if (i == j)
            {
                snf = (fmpz_sgn(fmpz_mat_entry(A, i, i)) >= 0);
                if (i > 0)
                    snf &= fmpz_divisible(fmpz_mat_entry(A, i, i),
                            fmpz_mat_entry(A, i - 1, i - 1));
            }
            else
            {
                snf = fmpz_is_zero(fmpz_mat_entry(A, i, j));
            }
        }
    }

    return snf;
}

int
main(void)
{
    slong iter, fails;
    FLINT_TEST_INIT(state);

    flint_printf("snf_saunders_wan....");
    fflush(stdout);
    fails = 0;

    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        fmpz_mat_t A, S, S2;
        slong m, n, b, d, r;
        int equal;

        m = 2 + n_randint(state, 10);
        n = 2 + n_randint(state, 10);
        r = n_randint(state, FLINT_MIN(m, n) + 1);

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(S, m, n);
        fmpz_mat_init(S2, m, n);

        /* sparse */
        b = 1 + n_randint(state, 10) * n_randint(state, 10);
        d = n_randint(state, 2*m*n + 1);
        fmpz_mat_randrank(A, state, r, b);

        /* dense */
        if (n_randint(state, 2))
            fmpz_mat_randops(A, state, d);

        if (iter==45)
            fmpz_mat_print_pretty(A);

        fmpz_mat_snf_saunders_wan(S, A);

        if (!in_snf(S))
        {
            flint_printf("FAIL:\n");
            flint_printf("matrix not in snf!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(S); flint_printf("\n\n");
            abort();
        }

        fmpz_mat_snf_saunders_wan(S2, S);
        equal = fmpz_mat_equal(S, S2);

        if (!equal)
        {
            flint_printf("WRONG ANSWER:\n");
            flint_printf("snf of a matrix in snf should be the same!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(S); flint_printf("\n\n");
            fmpz_mat_print_pretty(S2); flint_printf("\n\n");
            fails++;
            fmpz_mat_clear(S2);
            fmpz_mat_clear(S);
            fmpz_mat_clear(A);
            continue;
        }

        fmpz_mat_snf_kannan_bachem(S2, A);
        equal = fmpz_mat_equal(S, S2);

        if (!equal)
        {
            flint_printf("WRONG ANSWER:\n");
            flint_printf("different methods should agree!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(S); flint_printf("\n\n");
            fmpz_mat_print_pretty(S2); flint_printf("\n\n");
            fails++;
        }
        flint_printf("%wd\n", iter);

        fmpz_mat_clear(S2);
        fmpz_mat_clear(S);
        fmpz_mat_clear(A);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("%wd / %wd correct (%wd failures)\n",
            iter - fails, iter, fails);
    if (fails > iter/1000)
    {
        flint_printf("FAIL:\n");
        flint_printf("too many wrong answers!\n");
        abort();
    }

    flint_printf("PASS\n");
    return 0;
}

