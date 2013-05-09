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
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "ulong_extras.h"


int
main(void)
{
    nmod_mat_t A;
    len_t i, m, n, d, r;
    mp_limb_t mod;
    flint_rand_t state;
    flint_randinit(state);

    printf("rank....");
    fflush(stdout);

    /* Maximally sparse matrices of given rank */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 20);
        n = n_randint(state, 20);
        mod = n_randtest_prime(state, 0);

        for (r = 0; r <= FLINT_MIN(m,n); r++)
        {
            nmod_mat_init(A, m, n, mod);
            nmod_mat_randrank(A, state, r);
            /* printf("SPARSE %ld\n", r);
            nmod_mat_print_pretty(A); */
            if (r != nmod_mat_rank(A))
            {
                printf("FAIL:\n");
                printf("wrong rank!\n");
                abort();
            }
            nmod_mat_clear(A);
        }
    }

    /* Dense */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 20);
        n = n_randint(state, 20);
        mod = n_randtest_prime(state, 0);

        for (r = 0; r <= FLINT_MIN(m,n); r++)
        {
            d = n_randint(state, 2*m*n + 1);
            nmod_mat_init(A, m, n, mod);
            nmod_mat_randrank(A, state, r);
            nmod_mat_randops(A, d, state);
            /*
            printf("DENSE %ld %ld\n", r, d);
            nmod_mat_print_pretty(A); */
            if (r != nmod_mat_rank(A))
            {
                printf("FAIL:\n");
                printf("wrong rank!\n");
                abort();
            }
            nmod_mat_clear(A);
        }
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
