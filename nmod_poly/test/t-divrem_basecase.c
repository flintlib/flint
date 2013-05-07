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
    int i, result;
    flint_rand_t state;
    flint_randinit(state);

    printf("divrem_basecase....");
    fflush(stdout);

    /* Check result of divrem */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, q, r, prod;

        mp_limb_t n;
        do
        {
            n = n_randtest_not_zero(state);
        } while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(q, n);
        nmod_poly_init(r, n);
        nmod_poly_init(prod, n);
        
        nmod_poly_randtest(a, state, n_randint(state, 200));
        do
        {
            nmod_poly_randtest(b, state, n_randint(state, 200));
        } while (b->length == 0);

        nmod_poly_divrem_basecase(q, r, a, b);
        nmod_poly_mul(prod, q, b);
        nmod_poly_add(prod, prod, r);

        result = (nmod_poly_equal(a, prod));
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(prod), printf("\n\n");
            nmod_poly_print(q), printf("\n\n");
            nmod_poly_print(r), printf("\n\n");
            printf("n = %ld\n", n);
            abort();
        }
        
        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(q);
        nmod_poly_clear(r);
        nmod_poly_clear(prod);
    }

    /* Check aliasing of a and q */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, q, r;

        mp_limb_t n;
        do
        {
            n = n_randtest(state);
        } while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(q, n);
        nmod_poly_init(r, n);
        nmod_poly_randtest(a, state, n_randint(state, 200));
        do
        {
            nmod_poly_randtest(b, state, n_randint(state, 200));
        } while (b->length == 0);

        nmod_poly_divrem_basecase(q, r, a, b);
        nmod_poly_divrem_basecase(a, r, a, b);

        result = (nmod_poly_equal(a, q));
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(b), printf("\n\n");
            nmod_poly_print(q), printf("\n\n");
            nmod_poly_print(r), printf("\n\n");
            printf("n = %ld\n", n);
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(q);
        nmod_poly_clear(r);
    }

    /* Check aliasing of b and q */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, q, r;

        mp_limb_t n;
        do
        {
            n = n_randtest(state);
        } while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(q, n);
        nmod_poly_init(r, n);
        nmod_poly_randtest(a, state, n_randint(state, 200));
        do
        {
            nmod_poly_randtest(b, state, n_randint(state, 200));
        } while (b->length == 0);

        nmod_poly_divrem_basecase(q, r, a, b);
        nmod_poly_divrem_basecase(b, r, a, b);

        result = (nmod_poly_equal(b, q));
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(b), printf("\n\n");
            nmod_poly_print(q), printf("\n\n");
            nmod_poly_print(r), printf("\n\n");
            printf("n = %ld\n", n);
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(q);
        nmod_poly_clear(r);
    }

    /* Check aliasing of a and r */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, q, r;

        mp_limb_t n;
        do
        {
            n = n_randtest(state);
        } while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(q, n);
        nmod_poly_init(r, n);
        nmod_poly_randtest(a, state, n_randint(state, 200));
        do
        {
            nmod_poly_randtest(b, state, n_randint(state, 200));
        } while (b->length == 0);

        nmod_poly_divrem_basecase(q, r, a, b);
        nmod_poly_divrem_basecase(q, a, a, b);

        result = (nmod_poly_equal(a, r));
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(b), printf("\n\n");
            nmod_poly_print(q), printf("\n\n");
            nmod_poly_print(r), printf("\n\n");
            printf("n = %ld\n", n);
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(q);
        nmod_poly_clear(r);
    }

    /* Check aliasing of b and r */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, q, r;

        mp_limb_t n;
        do
        {
            n = n_randtest(state);
        } while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(q, n);
        nmod_poly_init(r, n);
        nmod_poly_randtest(a, state, n_randint(state, 200));
        do
        {
            nmod_poly_randtest(b, state, n_randint(state, 200));
        } while (b->length == 0);

        nmod_poly_divrem_basecase(q, r, a, b);
        nmod_poly_divrem_basecase(q, b, a, b);

        result = (nmod_poly_equal(b, r));
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(b), printf("\n\n");
            nmod_poly_print(q), printf("\n\n");
            nmod_poly_print(r), printf("\n\n");
            printf("n = %ld\n", n);
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(q);
        nmod_poly_clear(r);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
