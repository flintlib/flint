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

    printf("exp_series_basecase....");
    fflush(stdout);

    /* Check exp(A+B) = exp(A) * exp(B) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t A, B, AB, expA, expB, expAB, S;
        long n;
        mp_limb_t mod;

        mod = n_randtest_prime(state, 0);
        n = n_randtest(state) % 100;
        n = FLINT_MIN(n, mod);

        nmod_poly_init(A, mod);
        nmod_poly_init(B, mod);
        nmod_poly_init(AB, mod);
        nmod_poly_init(expA, mod);
        nmod_poly_init(expB, mod);
        nmod_poly_init(expAB, mod);
        nmod_poly_init(S, mod);

        nmod_poly_randtest(A, state, n_randint(state, 100));
        nmod_poly_set_coeff_ui(A, 0, 0UL);
        nmod_poly_randtest(B, state, n_randint(state, 100));
        nmod_poly_set_coeff_ui(B, 0, 0UL);

        if (n_randlimb(state) % 100 == 0)
        {
            nmod_poly_zero(A);
            nmod_poly_set_coeff_ui(A, n_randlimb(state) % (n+5), \
                n_randtest_not_zero(state) % mod);
            nmod_poly_set_coeff_ui(A, 0, 0UL);
        }

        nmod_poly_exp_series_basecase(expA, A, n);
        nmod_poly_exp_series_basecase(expB, B, n);
        nmod_poly_add(AB, A, B);
        nmod_poly_exp_series(expAB, AB, n);
        nmod_poly_mullow(S, expA, expB, n);

        result = nmod_poly_equal(S, expAB);

        if (!result)
        {
            printf("FAIL:\n");
            printf("n = %ld, mod = %lu\n", n, mod);
            printf("A: "); nmod_poly_print(A), printf("\n\n");
            printf("B: "); nmod_poly_print(B), printf("\n\n");
            printf("exp(A): "); nmod_poly_print(expA), printf("\n\n");
            printf("exp(B): "); nmod_poly_print(expB), printf("\n\n");
            printf("exp(A+B):       "); nmod_poly_print(expAB), printf("\n\n");
            printf("exp(A)*exp(B): "); nmod_poly_print(S), printf("\n\n");
            abort();
        }

        nmod_poly_clear(A);
        nmod_poly_clear(B);
        nmod_poly_clear(AB);
        nmod_poly_clear(expA);
        nmod_poly_clear(expB);
        nmod_poly_clear(expAB);
        nmod_poly_clear(S);
    }

    /* Check aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t A, B;
        long n;
        mp_limb_t mod;
        mod = n_randtest_prime(state, 0);
        n = n_randtest(state) % 50;
        n = FLINT_MIN(n, mod);

        nmod_poly_init(A, mod);
        nmod_poly_init(B, mod);
        nmod_poly_randtest(A, state, n_randint(state, 50));
        nmod_poly_set_coeff_ui(A, 0, 0UL);

        nmod_poly_exp_series_basecase(B, A, n);
        nmod_poly_exp_series_basecase(A, A, n);

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
