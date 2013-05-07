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
    fmpz_mat_t A, B, ker;
    flint_rand_t state;
    long i, m, n, b, d, r, nullity, nulrank;

    printf("nullspace....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 10);
        n = n_randint(state, 10);

        for (r = 0; r <= FLINT_MIN(m,n); r++)
        {
            b = 1 + n_randint(state, 10) * n_randint(state, 10);
            d = n_randint(state, 2*m*n + 1);

            fmpz_mat_init(A, m, n);
            fmpz_mat_init(ker, n, n);
            fmpz_mat_init(B, m, n);

            fmpz_mat_randrank(A, state, r, b);
            /* Densify */
            if (n_randlimb(state) % 2)
                fmpz_mat_randops(A, state, d);

            nullity = fmpz_mat_nullspace(ker, A);
            nulrank = fmpz_mat_rank(ker);

            if (nullity != nulrank)
            {
                printf("FAIL:\n");
                printf("rank(ker) != nullity!\n");
                fmpz_mat_print_pretty(A);
                printf("\n");
                abort();
            }

            if (nullity + r != n)
            {
                printf("FAIL:\n");
                printf("nullity + rank != n\n");
                fmpz_mat_print_pretty(A);
                printf("\n");
                abort();
            }

            fmpz_mat_mul(B, A, ker);

            if (fmpz_mat_rank(B) != 0)
            {
                printf("FAIL:\n");
                printf("A * ker != 0\n");
                fmpz_mat_print_pretty(A);
                printf("\n");
                abort();
            }

            fmpz_mat_clear(A);
            fmpz_mat_clear(ker);
            fmpz_mat_clear(B);
        }
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
