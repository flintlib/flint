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

int
main(void)
{
    slong iter;
    FLINT_TEST_INIT(state);

    flint_printf("hnf_transform....");
    fflush(stdout);

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        fmpz_mat_t A, H, H2, U;
        fmpz_t det;
        slong m, n, b, d, r;
        int equal;

        m = n_randint(state, 10);
        n = n_randint(state, 10);
        r = n_randint(state, FLINT_MIN(m, n) + 1);

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(H, m, n);
        fmpz_mat_init(H2, m, n);
        fmpz_mat_init(U, m, m);

        /* sparse */
        b = 1 + n_randint(state, 10) * n_randint(state, 10);
        d = n_randint(state, 2*m*n + 1);
        fmpz_mat_randrank(A, state, r, b);

        /* dense */
        if (n_randint(state, 2))
            fmpz_mat_randops(A, state, d);

        fmpz_mat_hnf_transform(H, U, A);

        if (!fmpz_mat_is_in_hnf(H))
        {
            flint_printf("FAIL:\n");
            flint_printf("matrix not in hnf!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(H); flint_printf("\n\n");
            abort();
        }

        fmpz_init(det);

        fmpz_mat_det(det, U);
        if (!fmpz_is_pm1(det))
        {
            flint_printf("FAIL:\n");
            flint_printf("transformation matrices should have determinant +-1, U does not!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(U); flint_printf("\n\n");
            fmpz_mat_print_pretty(H); flint_printf("\n\n");
            abort();
        }

        fmpz_clear(det);

        fmpz_mat_mul(H2, U, A);
        equal = fmpz_mat_equal(H, H2);

        if (!equal)
        {
            flint_printf("FAIL:\n");
            flint_printf("multiplying by the transformation matrix should give the same HNF!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(U); flint_printf("\n\n");
            fmpz_mat_print_pretty(H); flint_printf("\n\n");
            fmpz_mat_print_pretty(H2); flint_printf("\n\n");
            abort();
        }

        fmpz_mat_hnf(H2, A);
        equal = fmpz_mat_equal(H, H2);

        if (!equal)
        {
            flint_printf("FAIL:\n");
            flint_printf("hnfs produced by different methods should be the same!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(H); flint_printf("\n\n");
            fmpz_mat_print_pretty(H2); flint_printf("\n\n");
            abort();
        }

        fmpz_mat_clear(U);
        fmpz_mat_clear(H2);
        fmpz_mat_clear(H);
        fmpz_mat_clear(A);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}

