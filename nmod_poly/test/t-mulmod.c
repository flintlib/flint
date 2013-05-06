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

    printf("mulmod....");
    fflush(stdout);

    /* Aliasing res and a */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, res, t, f;

        mp_limb_t n = n_randtest_prime(state, 0);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(f, n);
        nmod_poly_init(res, n);
        nmod_poly_init(t, n);

        nmod_poly_randtest(a, state, n_randint(state, 50));
        nmod_poly_randtest(b, state, n_randint(state, 50));
        do {
            nmod_poly_randtest(f, state, n_randint(state, 50));
        } while (nmod_poly_is_zero(f));

        nmod_poly_mulmod(res, a, b, f);
        nmod_poly_mulmod(a, a, b, f);

        result = (nmod_poly_equal(res, a));
        if (!result)
        {
            printf("FAIL:\n");
            printf("a:\n"); nmod_poly_print(a), printf("\n\n");
            printf("b:\n"); nmod_poly_print(b), printf("\n\n");
            printf("f:\n"); nmod_poly_print(f), printf("\n\n");
            printf("res1:\n"); nmod_poly_print(res), printf("\n\n");
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(f);
        nmod_poly_clear(res);
        nmod_poly_clear(t);
    }

    /* Aliasing res and b */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, res, t, f;

        mp_limb_t n = n_randtest_prime(state, 0);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(f, n);
        nmod_poly_init(res, n);
        nmod_poly_init(t, n);

        nmod_poly_randtest(a, state, n_randint(state, 50));
        nmod_poly_randtest(b, state, n_randint(state, 50));
        do {
            nmod_poly_randtest(f, state, n_randint(state, 50));
        } while (nmod_poly_is_zero(f));

        nmod_poly_mulmod(res, a, b, f);
        nmod_poly_mulmod(b, a, b, f);

        result = (nmod_poly_equal(res, b));
        if (!result)
        {
            printf("FAIL:\n");
            printf("a:\n"); nmod_poly_print(a), printf("\n\n");
            printf("b:\n"); nmod_poly_print(b), printf("\n\n");
            printf("f:\n"); nmod_poly_print(f), printf("\n\n");
            printf("res1:\n"); nmod_poly_print(res), printf("\n\n");
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(f);
        nmod_poly_clear(res);
        nmod_poly_clear(t);
    }

    /* Aliasing res and f */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, res, t, f;

        mp_limb_t n = n_randtest_prime(state, 0);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(f, n);
        nmod_poly_init(res, n);
        nmod_poly_init(t, n);

        nmod_poly_randtest(a, state, n_randint(state, 50));
        nmod_poly_randtest(b, state, n_randint(state, 50));
        do {
            nmod_poly_randtest(f, state, n_randint(state, 50));
        } while (nmod_poly_is_zero(f));

        nmod_poly_mulmod(res, a, b, f);
        nmod_poly_mulmod(f, a, b, f);

        result = (nmod_poly_equal(res, f));
        if (!result)
        {
            printf("FAIL:\n");
            printf("a:\n"); nmod_poly_print(a), printf("\n\n");
            printf("b:\n"); nmod_poly_print(b), printf("\n\n");
            printf("f:\n"); nmod_poly_print(f), printf("\n\n");
            printf("res1:\n"); nmod_poly_print(res), printf("\n\n");
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(f);
        nmod_poly_clear(res);
        nmod_poly_clear(t);
    }

    /* No aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, res1, res2, t, f;

        mp_limb_t n = n_randtest_prime(state, 0);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(f, n);
        nmod_poly_init(res1, n);
        nmod_poly_init(res2, n);
        nmod_poly_init(t, n);

        nmod_poly_randtest(a, state, n_randint(state, 50));
        nmod_poly_randtest(b, state, n_randint(state, 50));
        do {
            nmod_poly_randtest(f, state, n_randint(state, 50));
        } while (nmod_poly_is_zero(f));

        nmod_poly_mulmod(res1, a, b, f);
        nmod_poly_mul(res2, a, b);
        nmod_poly_divrem(t, res2, res2, f);

        result = (nmod_poly_equal(res1, res2));
        if (!result)
        {
            printf("FAIL:\n");
            printf("a:\n"); nmod_poly_print(a), printf("\n\n");
            printf("b:\n"); nmod_poly_print(b), printf("\n\n");
            printf("f:\n"); nmod_poly_print(f), printf("\n\n");
            printf("res1:\n"); nmod_poly_print(res1), printf("\n\n");
            printf("res2:\n"); nmod_poly_print(res2), printf("\n\n");
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(f);
        nmod_poly_clear(res1);
        nmod_poly_clear(res2);
        nmod_poly_clear(t);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
