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
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "ulong_extras.h"


int
main(void)
{
    flint_rand_t state;
    len_t i;

    printf("nullspace....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_mat_t A, B, ker;
        mp_limb_t mod;
        len_t m, n, d, r, nullity, nulrank;

        m = n_randint(state, 30);
        n = n_randint(state, 30);

        for (r = 0; r <= FLINT_MIN(m,n); r++)
        {
            mod = n_randtest_prime(state, 0);
            d = n_randint(state, 2*m*n + 1);

            nmod_mat_init(A, m, n, mod);
            nmod_mat_init(ker, n, n, mod);
            nmod_mat_init(B, m, n, mod);

            nmod_mat_randrank(A, state, r);
            /* Densify */
            if (n_randlimb(state) % 2)
                nmod_mat_randops(A, d, state);

            nullity = nmod_mat_nullspace(ker, A);
            nulrank = nmod_mat_rank(ker);

            if (nullity != nulrank)
            {
                printf("FAIL:\n");
                printf("rank(ker) != nullity!\n");
                nmod_mat_print_pretty(A);
                printf("\n");
                abort();
            }

            if (nullity + r != n)
            {
                printf("FAIL:\n");
                printf("nullity + rank != n\n");
                nmod_mat_print_pretty(A);
                printf("\n");
                abort();
            }

            nmod_mat_mul(B, A, ker);

            if (nmod_mat_rank(B) != 0)
            {
                printf("FAIL:\n");
                printf("A * ker != 0\n");
                nmod_mat_print_pretty(A);
                printf("\n");
                abort();
            }

            nmod_mat_clear(A);
            nmod_mat_clear(ker);
            nmod_mat_clear(B);
        }
    }

    flint_randclear(state);
    printf("PASS\n");
    return 0;
}
