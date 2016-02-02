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

    Copyright (C) 2016 Aaditya Thakkar

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz_mat.h"
#include "fmpz.h"

int
main(void)
{
    slong i;
    FLINT_TEST_INIT(state);
    
    flint_printf("strassen_mul....");
    fflush(stdout);

    for (i = 0; i < 20 * 10; i++)
    {
        fmpz_mat_t A, B, C, D;

        slong m, k, n;

        m = n_randint(state, 400);
        k = n_randint(state, 400);
        n = n_randint(state, 400);

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(B, n, k);
        fmpz_mat_init(C, m, k);
        fmpz_mat_init(D, m, k);

        fmpz_mat_randtest(A, state, n_randint(state, 200)+1);
        fmpz_mat_randtest(B, state, n_randint(state, 200)+1);

        fmpz_mat_mul_classical(C, A, B);
        fmpz_mat_mul_strassen(D, A, B);

        if (!fmpz_mat_equal(C, D))
        {
            flint_printf("FAIL: results not equal\n");
            fmpz_mat_print_pretty(A);
            fmpz_mat_print_pretty(B);
            fmpz_mat_print_pretty(C);
            fmpz_mat_print_pretty(D);
            abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_mat_clear(D);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
