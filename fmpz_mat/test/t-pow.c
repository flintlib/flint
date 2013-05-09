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

int main(void)
{
    len_t i;
    flint_rand_t state;

    printf("pow....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t A, B, C;
        len_t i, n;
        ulong e;

        n = n_randint(state, 10);
        e = n_randint(state, 20);

        fmpz_mat_init(A, n, n);
        fmpz_mat_init(B, n, n);
        fmpz_mat_init(C, n, n);

        fmpz_mat_randtest(A, state, n_randint(state, 200) + 1);
        fmpz_mat_randtest(B, state, n_randint(state, 200) + 1);

        /* Make sure noise in the output is ok */
        fmpz_mat_randtest(B, state, n_randint(state, 200) + 1);

        fmpz_mat_pow(B, A, e);

        fmpz_mat_one(C);
        for (i = 0; i < e; i++)
            fmpz_mat_mul(C, C, A);

        if (!fmpz_mat_equal(C, B))
        {
            printf("FAIL: results not equal\n");
            abort();
        }

        fmpz_mat_pow(A, A, e);

        if (!fmpz_mat_equal(A, B))
        {
            printf("FAIL: aliasing failed\n");
            abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
