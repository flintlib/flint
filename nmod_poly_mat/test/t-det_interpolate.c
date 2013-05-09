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
#include "nmod_poly.h"
#include "nmod_poly_mat.h"


int
main(void)
{
    flint_rand_t state;
    len_t i;

    printf("det_interpolate....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_mat_t A;
        nmod_poly_t a, b;
        len_t n, deg;
        mp_limb_t mod;

        mod = n_randtest_prime(state, 0);
        n = n_randint(state, 10);
        deg = 1 + n_randint(state, 5);

        nmod_poly_mat_init(A, n, n, mod);

        nmod_poly_init(a, mod);
        nmod_poly_init(b, mod);

        nmod_poly_mat_randtest(A, state, deg);

        nmod_poly_mat_det(a, A);
        nmod_poly_mat_det_interpolate(b, A);

        if (!nmod_poly_equal(a, b))
        {
            printf("FAIL:\n");
            printf("determinants don't agree!\n");
            printf("A:\n");
            nmod_poly_mat_print(A, "x");
            printf("det(A):\n");
            nmod_poly_print(a);
            printf("\ndet_interpolate(A):\n");
            nmod_poly_print(b);
            printf("\n");
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);

        nmod_poly_mat_clear(A);
    }

    flint_randclear(state);
    printf("PASS\n");
    return 0;
}
