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
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    fmpz_mat_t A, X, B, AX, AXm, Bm;
    fmpz_t mod;
    flint_rand_t state;
    long i, m, n, r;
    int success;

    printf("solve_dixon....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 20);
        n = n_randint(state, 20);

        fmpz_mat_init(A, m, m);
        fmpz_mat_init(B, m, n);
        fmpz_mat_init(Bm, m, n);
        fmpz_mat_init(X, m, n);
        fmpz_mat_init(AX, m, n);
        fmpz_mat_init(AXm, m, n);
        fmpz_init(mod);

        fmpz_mat_randrank(A, state, m, 1+n_randint(state, 2)*n_randint(state, 100));
        fmpz_mat_randtest(B, state, 1+n_randint(state, 2)*n_randint(state, 100));

        /* Dense */
        if (n_randint(state, 2))
            fmpz_mat_randops(A, state, 1+n_randint(state, 1 + m*m));

        success = fmpz_mat_solve_dixon(X, mod, A, B);

        fmpz_mat_set(AXm, X);

        fmpz_mat_mul(AX, A, AXm);

        fmpz_mat_scalar_mod_fmpz(AXm, AX, mod);
        fmpz_mat_scalar_mod_fmpz(Bm, B, mod);

        if (!fmpz_mat_equal(AXm, Bm) || !success)
        {
            printf("FAIL:\n");
            printf("AX != B!\n");
            printf("A:\n"),      fmpz_mat_print_pretty(A),  printf("\n");
            printf("B:\n"),      fmpz_mat_print_pretty(B),  printf("\n");
            printf("X:\n"),      fmpz_mat_print_pretty(X),  printf("\n");
            printf("mod = "),    fmpz_print(mod),           printf("\n");
            printf("AX:\n"),     fmpz_mat_print_pretty(AX), printf("\n");
            abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(Bm);
        fmpz_mat_clear(X);
        fmpz_mat_clear(AX);
        fmpz_mat_clear(AXm);
        fmpz_clear(mod);
    }

    /* Test singular systems */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        m = 1 + n_randint(state, 10);
        n = 1 + n_randint(state, 10);
        r = n_randint(state, m);

        fmpz_mat_init(A, m, m);
        fmpz_mat_init(B, m, n);
        fmpz_mat_init(X, m, n);
        fmpz_init(mod);

        fmpz_mat_randrank(A, state, r, 1+n_randint(state, 2)*n_randint(state, 100));
        fmpz_mat_randtest(B, state, 1+n_randint(state, 2)*n_randint(state, 100));

        /* Dense */
        if (n_randint(state, 2))
            fmpz_mat_randops(A, state, 1+n_randint(state, 1 + m*m));

        if (fmpz_mat_solve_dixon(X, mod, A, B) != 0)
        {
            printf("FAIL:\n");
            printf("singular system, returned nonzero\n");
            abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(X);
        fmpz_clear(mod);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
