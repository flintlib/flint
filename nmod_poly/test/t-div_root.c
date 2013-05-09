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
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("div_root....");
    fflush(stdout);

    flint_randinit(state);

    /* Compare with standard divrem */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t P, Q, D, DQ, DR;
        mp_limb_t mod, r, rem;
        len_t n;

        mod = n_randtest_prime(state, 0);
        n = n_randint(state, 100);
        r = n_randint(state, mod);

        nmod_poly_init(P, mod);
        nmod_poly_init(Q, mod);
        nmod_poly_init(D, mod);
        nmod_poly_init(DQ, mod);
        nmod_poly_init(DR, mod);

        nmod_poly_randtest(P, state, n);

        rem = nmod_poly_div_root(Q, P, r);

        nmod_poly_set_coeff_ui(D, 0, n_negmod(r, mod));
        nmod_poly_set_coeff_ui(D, 1, 1UL);

        nmod_poly_divrem(DQ, DR, P, D);

        result = nmod_poly_equal(Q, DQ) &&
            (rem == nmod_poly_get_coeff_ui(DR, 0));

        if (!result)
        {
            printf("FAIL!\n");
            printf("P:\n"); nmod_poly_print(P); printf("\n\n");
            printf("Q:\n"); nmod_poly_print(Q); printf("\n\n");
            printf("D:\n"); nmod_poly_print(D); printf("\n\n");
            printf("DQ:\n"); nmod_poly_print(DQ); printf("\n\n");
            printf("DR:\n"); nmod_poly_print(DR); printf("\n\n");
            abort();
        }

        nmod_poly_clear(P);
        nmod_poly_clear(Q);
        nmod_poly_clear(D);
        nmod_poly_clear(DQ);
        nmod_poly_clear(DR);
    }

    /* Check aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t P, Q1, Q2;
        mp_limb_t mod, r, rem1, rem2;
        len_t n;

        mod = n_randtest_prime(state, 0);
        n = n_randint(state, 100);
        r = n_randint(state, mod);

        nmod_poly_init(P, mod);
        nmod_poly_init(Q1, mod);
        nmod_poly_init(Q2, mod);

        nmod_poly_randtest(P, state, n);
        nmod_poly_set(Q2, P);

        rem1 = nmod_poly_div_root(Q1, P, r);
        rem2 = nmod_poly_div_root(Q2, Q2, r);

        result = nmod_poly_equal(Q1, Q2) && (rem1 == rem2);

        if (!result)
        {
            printf("FAIL (aliasing)!\n");
            printf("P:\n"); nmod_poly_print(P); printf("\n\n");
            printf("Q1:\n"); nmod_poly_print(Q1); printf("\n\n");
            printf("Q2:\n"); nmod_poly_print(Q2); printf("\n\n");
            abort();
        }

        nmod_poly_clear(P);
        nmod_poly_clear(Q1);
        nmod_poly_clear(Q2);
    }

    flint_randclear(state);
    printf("PASS\n");
    return 0;
}
