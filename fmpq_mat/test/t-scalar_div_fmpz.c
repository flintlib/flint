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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpq.h"
#include "fmpq_mat.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;
    flint_randinit(state);

    printf("scalar_div_fmpz....");
    fflush(stdout);

    /* Aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t A, B;
        fmpz_t x;

        long m, n, bits;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        bits = 1 + n_randint(state, 100);

        fmpq_mat_init(A, m, n);
        fmpq_mat_init(B, m, n);
        fmpz_init(x);

        fmpq_mat_randtest(B, state, bits);
        fmpz_randtest_not_zero(x, state, bits);

        fmpq_mat_scalar_div_fmpz(A, B, x);
        fmpq_mat_scalar_div_fmpz(B, B, x);

        result = fmpq_mat_equal(A, B);
        if (!result)
        {
            printf("FAIL:\n");
            printf("A:\n"), fmpq_mat_print(A);
            printf("B:\n"), fmpq_mat_print(B);
            abort();
        }

        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
        fmpz_clear(x);
    }

    /* (A + B) / x == A / x + B / x */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t A, B, C, D;
        fmpz_t x;

        long m, n, bits;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        bits = 1 + n_randint(state, 100);

        fmpq_mat_init(A, m, n);
        fmpq_mat_init(B, m, n);
        fmpq_mat_init(C, m, n);
        fmpq_mat_init(D, m, n);
        fmpz_init(x);

        fmpq_mat_randtest(A, state, bits);
        fmpq_mat_randtest(B, state, bits);
        fmpz_randtest_not_zero(x, state, bits);

        fmpq_mat_scalar_div_fmpz(C, A, x);
        fmpq_mat_scalar_div_fmpz(D, B, x);
        fmpq_mat_add(D, C, D);

        fmpq_mat_add(C, A, B);
        fmpq_mat_scalar_div_fmpz(C, C, x);

        result = fmpq_mat_equal(C, D);
        if (!result)
        {
            printf("FAIL:\n");
            printf("A:\n"), fmpq_mat_print(A);
            printf("B:\n"), fmpq_mat_print(B);
            abort();
        }

        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
        fmpq_mat_clear(C);
        fmpq_mat_clear(D);
        fmpz_clear(x);
    }

    flint_randclear(state);

    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
