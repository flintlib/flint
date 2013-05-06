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

    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;
    flint_randinit(state);

    printf("powmod_ui_binexp....");
    fflush(stdout);

    /* Aliasing of res and a */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, res1, t, f;
        mp_limb_t n;
        ulong exp;

        n = n_randtest_prime(state, 0);
        exp = n_randlimb(state) % 32;

        nmod_poly_init(a, n);
        nmod_poly_init(f, n);
        nmod_poly_init(res1, n);
        nmod_poly_init(t, n);

        nmod_poly_randtest(a, state, n_randint(state, 50));
        do {
            nmod_poly_randtest(f, state, n_randint(state, 50));
        } while (nmod_poly_is_zero(f));

        nmod_poly_powmod_ui_binexp(res1, a, exp, f);
        nmod_poly_powmod_ui_binexp(a, a, exp, f);

        result = (nmod_poly_equal(res1, a));
        if (!result)
        {
            printf("FAIL:\n");
            printf("exp: %lu\n\n", exp);
            printf("a:\n"); nmod_poly_print(a), printf("\n\n");
            printf("f:\n"); nmod_poly_print(f), printf("\n\n");
            printf("res1:\n"); nmod_poly_print(res1), printf("\n\n");
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(f);
        nmod_poly_clear(res1);
        nmod_poly_clear(t);
    }

    /* Aliasing of res and f */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, res1, t, f;
        mp_limb_t n;
        ulong exp;

        n = n_randtest_prime(state, 0);
        exp = n_randlimb(state) % 32;

        nmod_poly_init(a, n);
        nmod_poly_init(f, n);
        nmod_poly_init(res1, n);
        nmod_poly_init(t, n);

        nmod_poly_randtest(a, state, n_randint(state, 50));
        do {
            nmod_poly_randtest(f, state, n_randint(state, 50));
        } while (nmod_poly_is_zero(f));

        nmod_poly_powmod_ui_binexp(res1, a, exp, f);
        nmod_poly_powmod_ui_binexp(f, a, exp, f);

        result = (nmod_poly_equal(res1, f));
        if (!result)
        {
            printf("FAIL:\n");
            printf("exp: %lu\n\n", exp);
            printf("a:\n"); nmod_poly_print(a), printf("\n\n");
            printf("f:\n"); nmod_poly_print(f), printf("\n\n");
            printf("res1:\n"); nmod_poly_print(res1), printf("\n\n");
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(f);
        nmod_poly_clear(res1);
        nmod_poly_clear(t);
    }

    /* No aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, res1, res2, t, f;
        mp_limb_t n;
        ulong exp;
        int j;

        n = n_randtest_prime(state, 0);
        exp = n_randlimb(state) % 32;

        nmod_poly_init(a, n);
        nmod_poly_init(f, n);
        nmod_poly_init(res1, n);
        nmod_poly_init(res2, n);
        nmod_poly_init(t, n);

        nmod_poly_randtest(a, state, n_randint(state, 50));
        do {
            nmod_poly_randtest(f, state, n_randint(state, 50));
        } while (nmod_poly_is_zero(f));

        nmod_poly_powmod_ui_binexp(res1, a, exp, f);

        nmod_poly_zero(res2);
        nmod_poly_set_coeff_ui(res2, 0, 1);
        for (j = 1; j <= exp; j++)
            nmod_poly_mulmod(res2, res2, a, f);

        result = (nmod_poly_equal(res1, res2));
        if (!result)
        {
            printf("FAIL:\n");
            printf("exp: %lu\n\n", exp);
            printf("a:\n"); nmod_poly_print(a), printf("\n\n");
            printf("f:\n"); nmod_poly_print(f), printf("\n\n");
            printf("res1:\n"); nmod_poly_print(res1), printf("\n\n");
            printf("res2:\n"); nmod_poly_print(res2), printf("\n\n");
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(f);
        nmod_poly_clear(res1);
        nmod_poly_clear(res2);
        nmod_poly_clear(t);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
