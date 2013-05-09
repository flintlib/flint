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

    printf("tanh_series....");
    fflush(stdout);

    /* Check atanh(tanh(A)) = A */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t A, tanhA, B;
        len_t n;
        mp_limb_t mod;

        do { mod = n_randtest_prime(state, 0); } while (mod == 2);
        n = 1 + n_randtest(state) % 100;
        n = FLINT_MIN(n, mod);

        nmod_poly_init(A, mod);
        nmod_poly_init(tanhA, mod);
        nmod_poly_init(B, mod);

        nmod_poly_randtest(A, state, n_randint(state, 100));
        nmod_poly_set_coeff_ui(A, 0, 0UL);

        nmod_poly_tanh_series(tanhA, A, n);
        nmod_poly_atanh_series(B, tanhA, n);

        nmod_poly_truncate(A, n);

        result = nmod_poly_equal(A, B);

        if (!result)
        {
            printf("FAIL:\n");
            printf("n = %ld, mod = %lu\n", n, mod);
            printf("A: "); nmod_poly_print(A), printf("\n\n");
            printf("tanh(A): "); nmod_poly_print(tanhA), printf("\n\n");
            printf("B: "); nmod_poly_print(B), printf("\n\n");
            abort();
        }

        nmod_poly_clear(A);
        nmod_poly_clear(tanhA);
        nmod_poly_clear(B);
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

        nmod_poly_tanh_series(B, A, n);
        nmod_poly_tanh_series(A, A, n);

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
