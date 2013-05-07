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

    Copyright (C) 2011 Fredrik Johansson

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
    int i;
    flint_rand_t state;
    flint_randinit(state);

    printf("solve_fraction_free....");
    fflush(stdout);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t A, B, X, AX;
        fmpq_t d;
        int success;

        long n, m, bits;

        n = n_randint(state, 10);
        m = n_randint(state, 10);
        bits = 1 + n_randint(state, 100);

        fmpq_mat_init(A, n, n);
        fmpq_mat_init(B, n, m);
        fmpq_mat_init(X, n, m);
        fmpq_mat_init(AX, n, m);

        fmpq_init(d);
        /* XXX: replace with a randtest function */
        do {
            fmpq_mat_randtest(A, state, bits);
            fmpq_mat_det(d, A);
        } while (fmpq_is_zero(d));
        fmpq_clear(d);

        fmpq_mat_randtest(B, state, bits);

        success = fmpq_mat_solve_fraction_free(X, A, B);
        fmpq_mat_mul(AX, A, X);

        if (!fmpq_mat_equal(AX, B) || !success)
        {
            printf("FAIL!\n");
            printf("success: %d\n", success);
            printf("A:\n");
            fmpq_mat_print(A);
            printf("B:\n");
            fmpq_mat_print(B);
            printf("X:\n");
            fmpq_mat_print(X);
            printf("AX:\n");
            fmpq_mat_print(AX);
            abort();
        }

        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
        fmpq_mat_clear(X);
        fmpq_mat_clear(AX);
    }

    /* Check singular systems */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t A, B, X;
        fmpz_mat_t M;
        fmpz_t den;
        long n, m, bits;
        int success;

        n = 1 + n_randint(state, 10);
        m = 1 + n_randint(state, 10);
        bits = 1 + n_randint(state, 100);

        fmpz_init(den);
        fmpz_mat_init(M, n, n);
        fmpz_mat_randrank(M, state, n_randint(state, n), bits);
        if (i % 2)
            fmpz_mat_randops(M, state, n_randint(state, 2*m*n + 1));
        fmpz_randtest_not_zero(den, state, bits);
        fmpq_mat_init(A, n, n);
        fmpq_mat_set_fmpz_mat_div_fmpz(A, M, den);

        fmpq_mat_init(B, n, m);
        fmpq_mat_randtest(B, state, bits);
        fmpq_mat_init(X, n, m);

        success = fmpq_mat_solve_fraction_free(X, A, B);

        if (success != 0)
        {
            printf("FAIL!\n");
            printf("Expected success = 0\n");
            fmpq_mat_print(A);
            printf("\n");
            abort();
        }

        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
        fmpq_mat_clear(X);
        fmpz_mat_clear(M);
        fmpz_clear(den);
    }

    flint_randclear(state);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
