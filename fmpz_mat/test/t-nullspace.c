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
    slong i, m, n, b, d, r, nullity, nulrank;

    FLINT_TEST_INIT(state);

    flint_printf("nullspace....");
    fflush(stdout);    

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
                flint_printf("FAIL:\n");
                flint_printf("rank(ker) != nullity!\n");
                fmpz_mat_print_pretty(A);
                flint_printf("\n");
                abort();
            }

            if (nullity + r != n)
            {
                flint_printf("FAIL:\n");
                flint_printf("nullity + rank != n\n");
                fmpz_mat_print_pretty(A);
                flint_printf("\n");
                abort();
            }

            fmpz_mat_mul(B, A, ker);

            if (fmpz_mat_rank(B) != 0)
            {
                flint_printf("FAIL:\n");
                flint_printf("A * ker != 0\n");
                fmpz_mat_print_pretty(A);
                flint_printf("\n");
                abort();
            }

            fmpz_mat_clear(A);
            fmpz_mat_clear(ker);
            fmpz_mat_clear(B);
        }
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
