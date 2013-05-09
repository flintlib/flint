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
#include "flint.h"
#include "nmod_mat.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"

int
main(void)
{
    flint_rand_t state;
    int iter;

    printf("zero/is_zero....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        nmod_poly_mat_t A;
        len_t m, n;
        mp_limb_t mod;

        mod = n_randtest_prime(state, 0);
        m = n_randint(state, 10);
        n = n_randint(state, 10);

        nmod_poly_mat_init(A, m, n, mod);
        nmod_poly_mat_randtest(A, state, n_randint(state, 5));
        nmod_poly_mat_zero(A);

        if (!nmod_poly_mat_is_zero(A))
        {
            printf("FAIL: expected matrix to be zero\n");
            abort();
        }

        if (m > 0 && n > 0)
        {
            m = n_randint(state, m);
            n = n_randint(state, n);
            nmod_poly_randtest_not_zero(nmod_poly_mat_entry(A, m, n), state, 5);

            if (nmod_poly_mat_is_zero(A))
            {
                printf("FAIL: expected matrix not to be zero\n");
                abort();
            }
        }

        nmod_poly_mat_clear(A);
    }

    flint_randclear(state);
    printf("PASS\n");
    return 0;
}
