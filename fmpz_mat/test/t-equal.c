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

    Copyright (C) 2010 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    fmpz_randstate_t rnd;

    printf("equal....");
    fflush(stdout);

    fmpz_randinit(rnd);

    for (i = 0; i < 10000; i++)
    {
        fmpz_mat_t A, B, C, D, E;
        long m, n, j;

        m = n_randint(20);
        n = n_randint(20);

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(B, m, n);
        fmpz_mat_init(C, m, n);
        fmpz_mat_init(D, m+1, n);
        fmpz_mat_init(E, m, n+1);

        if (fmpz_mat_equal(A, D) || fmpz_mat_equal(A, E))
        {
            printf("FAIL: different dimensions should not be equal\n");
            abort();
        }

        fmpz_mat_randtest(A, rnd, 1 + n_randint(100));
        fmpz_mat_copy(B, A);

        if (!fmpz_mat_equal(A, B))
        {
            printf("FAIL: copied matrices should be equal\n");
            abort();
        }

        if (m && n)
        {
            j = n_randint(m * n);
            fmpz_add_ui(A->entries + j, A->entries + j, 1);

            if (fmpz_mat_equal(A, B))
            {
                printf("FAIL: modified matrices should not be equal\n");
                abort();
            }
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_mat_clear(D);
        fmpz_mat_clear(E);
    }

    fmpz_randclear(rnd);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
