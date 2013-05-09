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
    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result = 1;
    flint_rand_t state;
    flint_randinit(state);

    printf("atan_series....");
    fflush(stdout);

    /* Check 2*atan(A) = atan(2*A/(1-A^2)) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t A, B, atanA, atanB;
        len_t n;
        mp_limb_t mod;

        mod = n_randtest_prime(state, 0);
        n = 1 + n_randtest(state) % 100;
        n = FLINT_MIN(n, mod);

        nmod_poly_init(A, mod);
        nmod_poly_init(B, mod);
        nmod_poly_init(atanA, mod);
        nmod_poly_init(atanB, mod);

        nmod_poly_randtest(A, state, n_randint(state, 100));
        nmod_poly_set_coeff_ui(A, 0, 0UL);

        nmod_poly_mullow(B, A, A, n);
        nmod_poly_neg(B, B);
        nmod_poly_set_coeff_ui(B, 0, 1UL);
        nmod_poly_div_series(B, A, B, n);
        nmod_poly_add(B, B, B);

        nmod_poly_atan_series(atanA, A, n);
        nmod_poly_atan_series(atanB, B, n);
        nmod_poly_add(atanA, atanA, atanA);

        result = nmod_poly_equal(atanA, atanB);

        if (!result)
        {
            printf("FAIL:\n");
            printf("n = %ld, mod = %lu\n", n, mod);
            printf("A: "); nmod_poly_print(A), printf("\n\n");
            printf("B: "); nmod_poly_print(B), printf("\n\n");
            printf("2*atan(A): "); nmod_poly_print(atanA), printf("\n\n");
            printf("atan(B): "); nmod_poly_print(atanB), printf("\n\n");
            abort();
        }

        nmod_poly_clear(A);
        nmod_poly_clear(B);
        nmod_poly_clear(atanA);
        nmod_poly_clear(atanB);
    }

    /* Check aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t A, B;
        len_t n;
        mp_limb_t mod;
        mod = n_randtest_prime(state, 0);
        n = n_randtest(state) % 50;
        n = FLINT_MIN(n, mod);

        nmod_poly_init(A, mod);
        nmod_poly_init(B, mod);
        nmod_poly_randtest(A, state, n_randint(state, 50));
        nmod_poly_set_coeff_ui(A, 0, 0UL);

        nmod_poly_atan_series(B, A, n);
        nmod_poly_atan_series(A, A, n);

        result = nmod_poly_equal(A, B);
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(A), printf("\n\n");
            nmod_poly_print(B), printf("\n\n");
            abort();
        }

        nmod_poly_clear(A);
        nmod_poly_clear(B);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
