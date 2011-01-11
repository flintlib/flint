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
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "ulong_extras.h"


int
main(void)
{
    fmpz_mat_t A, x, b, Ax;
    fmpz_t den;
    flint_rand_t state;
    long i, m, r;

    printf("solve....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 20000; i++)
    {
        m = n_randint(state, 10);

        fmpz_mat_init(A, m, m);
        fmpz_mat_init(b, m, 1);
        fmpz_mat_init(x, m, 1);
        fmpz_mat_init(Ax, m, 1);
        fmpz_init(den);

        fmpz_mat_randrank(A, state, m, 1+n_randint(state, 2)*n_randint(state, 100));
        fmpz_mat_randtest(b, state, 1+n_randint(state, 2)*n_randint(state, 100));

        /* Dense */
        if (n_randint(state, 2))
            fmpz_mat_randops(A, state, 1+n_randint(state, 1 + m*m));

        fmpz_mat_solve(x->entries, den, A, b->entries);

        fmpz_mat_mul(Ax, A, x);
        _fmpz_vec_scalar_divexact_fmpz(Ax->entries, Ax->entries, m, den);

        if (!fmpz_mat_equal(Ax, b))
        {
            printf("FAIL:\n");
            printf("Ax != b!\n");
            printf("A:\n"),      fmpz_mat_print_pretty(A),  printf("\n");
            printf("b:\n"),      fmpz_mat_print_pretty(b),  printf("\n");
            printf("x:\n"),      fmpz_mat_print_pretty(x),  printf("\n");
            printf("den(x) = "), fmpz_print(den),           printf("\n");
            printf("Ax:\n"),     fmpz_mat_print_pretty(Ax), printf("\n");
            abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(b);
        fmpz_mat_clear(x);
        fmpz_mat_clear(Ax);
        fmpz_clear(den);
    }

    /* Test singular systems */
    for (i = 0; i < 10000; i++)
    {
        m = 1 + n_randint(state, 10);
        r = n_randint(state, m);

        fmpz_mat_init(A, m, m);
        fmpz_mat_init(b, m, 1);
        fmpz_mat_init(x, m, 1);
        fmpz_mat_init(Ax, m, 1);
        fmpz_init(den);

        fmpz_mat_randrank(A, state, r, 1+n_randint(state, 2)*n_randint(state, 100));
        fmpz_mat_randtest(b, state, 1+n_randint(state, 2)*n_randint(state, 100));

        /* Dense */
        if (n_randint(state, 2))
            fmpz_mat_randops(A, state, 1+n_randint(state, 1 + m*m));

        fmpz_mat_solve(x->entries, den, A, b->entries);

        if (!fmpz_is_zero(den))
        {
            printf("FAIL:\n");
            printf("singular system gave nonzero determinant\n");
            abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(b);
        fmpz_mat_clear(x);
        fmpz_mat_clear(Ax);
        fmpz_clear(den);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
