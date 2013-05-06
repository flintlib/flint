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
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "ulong_extras.h"


int
main(void)
{
    flint_rand_t state;
    long i;

    printf("solve_bound....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t A, B, X;
        fmpz_t N, D, den;
        long m, n, b1, b2;
        long j, k;

        b1 = 1 + n_randint(state, 100);
        b2 = 1 + n_randint(state, 100);
        m = n_randint(state, 20);
        n = n_randint(state, 20);

        fmpz_init(den);
        fmpz_init(N);
        fmpz_init(D);
        fmpz_mat_init(A, m, m);
        fmpz_mat_init(B, m, n);
        fmpz_mat_init(X, m, n);

        fmpz_mat_randrank(A, state, m, b1);
        fmpz_mat_randops(A, state, n_randint(state, m)*n_randint(state, m));
        fmpz_mat_randtest(B, state, b2);

        fmpz_mat_solve_bound(N, D, A, B);
        fmpz_mat_solve(X, den, A, B);

        if (fmpz_cmpabs(D, den) < 0)
        {
            printf("FAIL:\n");
            printf("denominator bound:\n");
            fmpz_print(D);
            printf("\ndenominator:\n");
            fmpz_print(den);
            printf("\n");
            printf("A:\n");
            fmpz_mat_print_pretty(A);
            printf("B:\n");
            fmpz_mat_print_pretty(B);
            printf("\n");
            abort();
        }

        for (j = 0; j < m; j++)
        {
            for (k = 0; k < n; k++)
            {
                if (fmpz_cmpabs(N, fmpz_mat_entry(X, j, k)) < 0)
                {
                    printf("FAIL:\n");
                    printf("numerator bound:\n");
                    fmpz_print(N);
                    printf("\nnumerator:\n");
                    fmpz_print(fmpz_mat_entry(X, j, k));
                    printf("\n");
                    printf("A:\n");
                    fmpz_mat_print_pretty(A);
                    printf("B:\n");
                    fmpz_mat_print_pretty(B);
                    printf("\n");
                    abort();
                }
            }
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(X);

        fmpz_clear(den);
        fmpz_clear(N);
        fmpz_clear(D);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
