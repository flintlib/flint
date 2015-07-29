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

/* 
   Though fmpz_mat_hnf_pernet_stein may fail to give the correct result, so
   that this test code could in theory fail with low probability, we always
   generate the same random values, so it should always pass.
*/

int
main(void)
{
    slong iter;
    FLINT_TEST_INIT(state);

    flint_printf("hnf_pernet_stein....");
    fflush(stdout);

    /* matrices of random rank */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        fmpz_mat_t A, B, H, H2;
        slong m, n, r, b, d;
        int equal, ok;

        n = 1 + n_randint(state, 10);
        m = 1 + n_randint(state, 10);
        r = n_randint(state, FLINT_MIN(m, n) + 1);

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(B, m, n);
        fmpz_mat_init(H, m, n);
        fmpz_mat_init(H2, m, n);

        /* sparse */
        b = 1 + n_randint(state, 10) * n_randint(state, 10);
        d = n_randint(state, 2*m*n + 1);
        fmpz_mat_randrank(A, state, r, b);

        /* dense */
        if (n_randint(state, 2))
            fmpz_mat_randops(A, state, d);

        ok = fmpz_mat_hnf_pernet_stein(H, A, state);

        if (ok && !fmpz_mat_is_in_hnf(H))
        {
            flint_printf("FAIL:\n");
            flint_printf("matrix not in hnf!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(H); flint_printf("\n\n");
            abort();
        }

        fmpz_mat_hnf_classical(H2, A);
        equal = fmpz_mat_equal(H, H2);

        if (ok && !equal)
        {
            flint_printf("FAIL:\n");
            flint_printf("hnfs produced by different methods should be the same!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(H); flint_printf("\n\n");
            fmpz_mat_print_pretty(H2); flint_printf("\n\n");
            abort();
        }

        ok = fmpz_mat_hnf_pernet_stein(H2, H, state);
        equal = fmpz_mat_equal(H, H2);

        if (ok && !equal)
        {
            flint_printf("FAIL:\n");
            flint_printf("hnf of a matrix in hnf should be the same!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(H); flint_printf("\n\n");
            fmpz_mat_print_pretty(H2); flint_printf("\n\n");
            abort();
        }

        fmpz_mat_clear(H2);
        fmpz_mat_clear(H);
        fmpz_mat_clear(B);
        fmpz_mat_clear(A);
    }

    /* matrices with random entries */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        fmpz_mat_t A, B, H, H2;
        slong m, n, b;
        int equal, ok;

        n = 1 + n_randint(state, 10);
        m = 1 + n_randint(state, 10);

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(B, m, n);
        fmpz_mat_init(H, m, n);
        fmpz_mat_init(H2, m, n);

        b = 1 + n_randint(state, 8) * n_randint(state, 8);
        fmpz_mat_randtest(A, state, b);

        ok = fmpz_mat_hnf_pernet_stein(H, A, state);

        if (ok && !fmpz_mat_is_in_hnf(H))
        {
            flint_printf("FAIL:\n");
            flint_printf("matrix not in hnf!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(H); flint_printf("\n\n");
            abort();
        }

        fmpz_mat_hnf_classical(H2, A);
        equal = fmpz_mat_equal(H, H2);

        if (ok && !equal)
        {
            flint_printf("FAIL:\n");
            flint_printf("hnfs produced by different methods should be the same!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(H); flint_printf("\n\n");
            fmpz_mat_print_pretty(H2); flint_printf("\n\n");
            abort();
        }

        ok = fmpz_mat_hnf_pernet_stein(H2, H, state);
        equal = fmpz_mat_equal(H, H2);

        if (ok && !equal)
        {
            flint_printf("FAIL:\n");
            flint_printf("hnf of a matrix in hnf should be the same!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(H); flint_printf("\n\n");
            fmpz_mat_print_pretty(H2); flint_printf("\n\n");
            abort();
        }

        fmpz_mat_clear(H2);
        fmpz_mat_clear(H);
        fmpz_mat_clear(B);
        fmpz_mat_clear(A);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
