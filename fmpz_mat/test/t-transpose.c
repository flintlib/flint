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
#include "fmpz_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    long m, n, rep;
    flint_rand_t state;

    printf("transpose....");
    fflush(stdout);

    flint_randinit(state);

    /* Rectangular transpose */
    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        fmpz_mat_t A, B, C;

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(B, n, m);
        fmpz_mat_init(C, m, n);

        fmpz_mat_randtest(A, state, 1+n_randint(state, 100));
        fmpz_mat_randtest(B, state, 1+n_randint(state, 100));

        fmpz_mat_transpose(B, A);
        fmpz_mat_transpose(C, B);

        if (!fmpz_mat_equal(C, A))
        {
            printf("FAIL: C != A\n");
            abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
    }

    /* Self-transpose */
    for (rep = 0; rep < 1000; rep++)
    {
        fmpz_mat_t A, B;

        m = n_randint(state, 20);

        fmpz_mat_init(A, m, m);
        fmpz_mat_init(B, m, m);

        fmpz_mat_randtest(A, state, 1+n_randint(state, 100));
        fmpz_mat_set(B, A);
        fmpz_mat_transpose(B, B);
        fmpz_mat_transpose(B, B);

        if (!fmpz_mat_equal(B, A))
        {
            printf("FAIL: B != A\n");
            abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
