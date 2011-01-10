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
    fmpz_mat_t A;
    flint_rand_t rnd;
    long i, m, n, b, d, r;

    printf("rank....");
    fflush(stdout);

    fmpz_randinit(rnd);

    /* Maximally sparse matrices of given rank */
    for (i = 0; i < 1000; i++)
    {
        m = n_randint(10);
        n = n_randint(10);

        for (r = 0; r <= FLINT_MIN(m,n); r++)
        {
            b = 1 + n_randint(10) * n_randint(10);
            d = n_randint(2*m*n + 1);
            fmpz_mat_init(A, m, n);
            fmpz_mat_randrank(A, rnd, r, b);
            if (r != fmpz_mat_rank(A))
            {
                printf("FAIL:\n");
                printf("wrong rank!\n");
                abort();
            }
            fmpz_mat_clear(A);
        }
    }

    /* Dense */
    for (i = 0; i < 1000; i++)
    {
        m = n_randint(10);
        n = n_randint(10);

        for (r = 0; r <= FLINT_MIN(m,n); r++)
        {
            b = 1 + n_randint(10) * n_randint(10);
            d = n_randint(2*m*n + 1);
            fmpz_mat_init(A, m, n);
            fmpz_mat_randrank(A, rnd, r, b);
            fmpz_mat_randops(A, rnd, d);
            if (r != fmpz_mat_rank(A))
            {
                printf("FAIL:\n");
                printf("wrong rank!\n");
                abort();
            }
            fmpz_mat_clear(A);
        }
    }

    fmpz_mat_randclear(rnd);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
