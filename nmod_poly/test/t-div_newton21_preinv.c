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
    Copyright (C) 2013 Martin Lee

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

    printf("div_newton21_preinv....");
    fflush(stdout);

    /* Check result of divrem */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, binv, q, r, test;

        mp_limb_t n;
        do n = n_randtest_not_zero(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(binv, n);
        nmod_poly_init(q, n);
        nmod_poly_init(r, n);
        nmod_poly_init(test, n);

        do
        nmod_poly_randtest(b, state, n_randint(state, 200));
        while (b->length <= 2);
        nmod_poly_randtest(a, state, n_randint(state, 200));
        if (a->length > 2*(b->length)-3)
          nmod_poly_truncate (a, 2*(b->length)-3);

        nmod_poly_reverse (binv, b, b->length);
        nmod_poly_inv_series (binv, binv, b->length);
        nmod_poly_div_newton21_preinv(q, a, b, binv);
        nmod_poly_divrem (test, r, a, b);

        result = (nmod_poly_equal(q, test));
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(test), printf("\n\n");
            nmod_poly_print(q), printf("\n\n");
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(b), printf("\n\n");
            printf("n = %ld\n", n);
            abort();
        }
        
        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(binv);
        nmod_poly_clear(q);
        nmod_poly_clear(r);
        nmod_poly_clear(test);
    }

    /* Check aliasing of a and q */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, binv, q;

        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(binv, n);
        nmod_poly_init(q, n);
        do
        nmod_poly_randtest(b, state, n_randint(state, 200));
        while (b->length <= 2);
        nmod_poly_randtest(a, state, n_randint(state, 200));
        if (a->length > 2*(b->length)-3)
          nmod_poly_truncate (a, 2*(b->length)-3);

        nmod_poly_reverse (binv, b, b->length);
        nmod_poly_inv_series (binv, binv, b->length);

        nmod_poly_div_newton21_preinv(q, a, b, binv);
        nmod_poly_div_newton21_preinv(a, a, b, binv);

        result = (nmod_poly_equal(a, q));
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(b), printf("\n\n");
            nmod_poly_print(q), printf("\n\n");
            printf("n = %ld\n", n);
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(binv);
        nmod_poly_clear(q);
    }

    /* Check aliasing of b and q */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, binv, q;

        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(binv, n);
        nmod_poly_init(q, n);
        do
        nmod_poly_randtest(b, state, n_randint(state, 200));
        while (b->length <= 2);
        nmod_poly_randtest(a, state, n_randint(state, 200));
        if (a->length > 2*(b->length)-3)
          nmod_poly_truncate (a, 2*(b->length)-3);

        nmod_poly_reverse (binv, b, b->length);
        nmod_poly_inv_series (binv, binv, b->length);

        nmod_poly_div_newton21_preinv(q, a, b, binv);
        nmod_poly_div_newton21_preinv(b, a, b, binv);

        result = (nmod_poly_equal(b, q));
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(b), printf("\n\n");
            nmod_poly_print(q), printf("\n\n");
            printf("n = %ld\n", n);
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(binv);
        nmod_poly_clear(q);
    }

    /* Check aliasing of binv and q */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, binv, q;

        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(binv, n);
        nmod_poly_init(q, n);
        do
        nmod_poly_randtest(b, state, n_randint(state, 200));
        while (b->length <= 2);
        nmod_poly_randtest(a, state, n_randint(state, 200));
        if (a->length > 2*(b->length)-3)
          nmod_poly_truncate (a, 2*(b->length)-3);

        nmod_poly_reverse (binv, b, b->length);
        nmod_poly_inv_series (binv, binv, b->length);

        nmod_poly_div_newton21_preinv(q, a, b, binv);
        nmod_poly_div_newton21_preinv(binv, a, b, binv);

        result = (nmod_poly_equal(binv, q));
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(b), printf("\n\n");
            nmod_poly_print(binv), printf("\n\n");
            nmod_poly_print(q), printf("\n\n");
            printf("n = %ld\n", n);
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(binv);
        nmod_poly_clear(q);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
