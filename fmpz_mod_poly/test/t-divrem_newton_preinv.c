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
    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;
    flint_randinit(state);

    printf("divrem_newton_preinv....");
    fflush(stdout);

    /* Check result of divrem */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, b, binv, q, r, test;

        fmpz_t p;
        mp_limb_t n;
        do n = n_randtest_not_zero(state);
        while (!n_is_probabprime(n));

        fmpz_init(p);
        fmpz_set_ui(p, n);

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(binv, p);
        fmpz_mod_poly_init(q, p);
        fmpz_mod_poly_init(r, p);
        fmpz_mod_poly_init(test, p);

        do
            fmpz_mod_poly_randtest(b, state, n_randint(state, 200));
        while (b->length <= 2);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 200));
        if (a->length > 2*(b->length)-3)
          fmpz_mod_poly_truncate (a, 2*(b->length)-3);

        fmpz_mod_poly_reverse (binv, b, b->length);
        fmpz_mod_poly_inv_series_newton (binv, binv, b->length);
        fmpz_mod_poly_divrem_newton_preinv(q, r, a, b, binv);
        fmpz_mod_poly_mul(test, q, b);
        fmpz_mod_poly_add(test, test, r);

        result = (fmpz_mod_poly_equal(a, test));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_mod_poly_print(a), printf("\n\n");
            fmpz_mod_poly_print(test), printf("\n\n");
            fmpz_mod_poly_print(q), printf("\n\n");
            fmpz_mod_poly_print(r), printf("\n\n");
            fmpz_mod_poly_print(b), printf("\n\n");
            printf("n = %ld\n", n);
            abort();
        }
        
        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(binv);
        fmpz_mod_poly_clear(q);
        fmpz_mod_poly_clear(r);
        fmpz_mod_poly_clear(test);
        fmpz_clear(p);
    }

    /* Check aliasing of a and q */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, b, binv, q, r;

        fmpz_t p;
        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        fmpz_init(p);
        fmpz_set_ui(p, n);

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(binv, p);
        fmpz_mod_poly_init(q, p);
        fmpz_mod_poly_init(r, p);
        do
        fmpz_mod_poly_randtest(b, state, n_randint(state, 200));
        while (b->length <= 2);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 200));
        if (a->length > 2*(b->length)-3)
          fmpz_mod_poly_truncate (a, 2*(b->length)-3);

        fmpz_mod_poly_reverse (binv, b, b->length);
        fmpz_mod_poly_inv_series_newton (binv, binv, b->length);

        fmpz_mod_poly_divrem_newton_preinv(q, r, a, b, binv);
        fmpz_mod_poly_divrem_newton_preinv(a, r, a, b, binv);

        result = (fmpz_mod_poly_equal(a, q));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_mod_poly_print(a), printf("\n\n");
            fmpz_mod_poly_print(b), printf("\n\n");
            fmpz_mod_poly_print(q), printf("\n\n");
            fmpz_mod_poly_print(r), printf("\n\n");
            printf("n = %ld\n", n);
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(binv);
        fmpz_mod_poly_clear(q);
        fmpz_mod_poly_clear(r);
        fmpz_clear(p);
    }

    /* Check aliasing of b and q */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, b, binv, q, r;

        fmpz_t p;
        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        fmpz_init(p);
        fmpz_set_ui(p, n);

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(binv, p);
        fmpz_mod_poly_init(q, p);
        fmpz_mod_poly_init(r, p);
        do
        fmpz_mod_poly_randtest(b, state, n_randint(state, 200));
        while (b->length <= 2);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 200));
        if (a->length > 2*(b->length)-3)
          fmpz_mod_poly_truncate (a, 2*(b->length)-3);

        fmpz_mod_poly_reverse (binv, b, b->length);
        fmpz_mod_poly_inv_series_newton (binv, binv, b->length);

        fmpz_mod_poly_divrem_newton_preinv(q, r, a, b, binv);
        fmpz_mod_poly_divrem_newton_preinv(b, r, a, b, binv);

        result = (fmpz_mod_poly_equal(b, q));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_mod_poly_print(a), printf("\n\n");
            fmpz_mod_poly_print(b), printf("\n\n");
            fmpz_mod_poly_print(q), printf("\n\n");
            fmpz_mod_poly_print(r), printf("\n\n");
            printf("n = %ld\n", n);
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(binv);
        fmpz_mod_poly_clear(q);
        fmpz_mod_poly_clear(r);
        fmpz_clear(p);
    }

    /* Check aliasing of binv and q */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, b, binv, q, r;

        fmpz_t p;
        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        fmpz_init(p);
        fmpz_set_ui(p, n);

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(binv, p);
        fmpz_mod_poly_init(q, p);
        fmpz_mod_poly_init(r, p);
        do
        fmpz_mod_poly_randtest(b, state, n_randint(state, 200));
        while (b->length <= 2);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 200));
        if (a->length > 2*(b->length)-3)
          fmpz_mod_poly_truncate (a, 2*(b->length)-3);

        fmpz_mod_poly_reverse (binv, b, b->length);
        fmpz_mod_poly_inv_series_newton (binv, binv, b->length);

        fmpz_mod_poly_divrem_newton_preinv(q, r, a, b, binv);
        fmpz_mod_poly_divrem_newton_preinv(binv, r, a, b, binv);

        result = (fmpz_mod_poly_equal(binv, q));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_mod_poly_print(a), printf("\n\n");
            fmpz_mod_poly_print(b), printf("\n\n");
            fmpz_mod_poly_print(binv), printf("\n\n");
            fmpz_mod_poly_print(q), printf("\n\n");
            fmpz_mod_poly_print(r), printf("\n\n");
            printf("n = %ld\n", n);
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(binv);
        fmpz_mod_poly_clear(q);
        fmpz_mod_poly_clear(r);
        fmpz_clear(p);
    }

    /* Check aliasing of a and r */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, b, binv, q, r;

        fmpz_t p;
        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        fmpz_init(p);
        fmpz_set_ui(p, n);

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(binv, p);
        fmpz_mod_poly_init(q, p);
        fmpz_mod_poly_init(r, p);
        do
        fmpz_mod_poly_randtest(b, state, n_randint(state, 200));
        while (b->length <= 2);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 200));
        if (a->length > 2*(b->length)-3)
          fmpz_mod_poly_truncate (a, 2*(b->length)-3);

        fmpz_mod_poly_reverse (binv, b, b->length);
        fmpz_mod_poly_inv_series_newton (binv, binv, b->length);

        fmpz_mod_poly_divrem_newton_preinv(q, r, a, b, binv);
        fmpz_mod_poly_divrem_newton_preinv(q, a, a, b, binv);

        result = (fmpz_mod_poly_equal(a, r));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_mod_poly_print(a), printf("\n\n");
            fmpz_mod_poly_print(b), printf("\n\n");
            fmpz_mod_poly_print(q), printf("\n\n");
            fmpz_mod_poly_print(r), printf("\n\n");
            printf("n = %ld\n", n);
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(binv);
        fmpz_mod_poly_clear(q);
        fmpz_mod_poly_clear(r);
        fmpz_clear(p);
    }

    /* Check aliasing of b and r */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, b, binv, q, r;

        fmpz_t p;
        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        fmpz_init(p);
        fmpz_set_ui(p, n);

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(binv, p);
        fmpz_mod_poly_init(q, p);
        fmpz_mod_poly_init(r, p);
        do
        fmpz_mod_poly_randtest(b, state, n_randint(state, 200));
        while (b->length <= 2);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 200));
        if (a->length > 2*(b->length)-3)
          fmpz_mod_poly_truncate (a, 2*(b->length)-3);

        fmpz_mod_poly_reverse (binv, b, b->length);
        fmpz_mod_poly_inv_series_newton (binv, binv, b->length);

        fmpz_mod_poly_divrem_newton_preinv(q, r, a, b, binv);
        fmpz_mod_poly_divrem_newton_preinv(q, b, a, b, binv);

        result = (fmpz_mod_poly_equal(b, r));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_mod_poly_print(a), printf("\n\n");
            fmpz_mod_poly_print(b), printf("\n\n");
            fmpz_mod_poly_print(q), printf("\n\n");
            fmpz_mod_poly_print(r), printf("\n\n");
            printf("n = %ld\n", n);
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(binv);
        fmpz_mod_poly_clear(q);
        fmpz_mod_poly_clear(r);
        fmpz_clear(p);
    }

    /* Check aliasing of binv and r */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, b, binv, q, r;

        fmpz_t p;
        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        fmpz_init(p);
        fmpz_set_ui(p, n);

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(binv, p);
        fmpz_mod_poly_init(q, p);
        fmpz_mod_poly_init(r, p);
        do
        fmpz_mod_poly_randtest(b, state, n_randint(state, 200));
        while (b->length <= 2);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 200));
        if (a->length > 2*(b->length)-3)
          fmpz_mod_poly_truncate (a, 2*(b->length)-3);

        fmpz_mod_poly_reverse (binv, b, b->length);
        fmpz_mod_poly_inv_series_newton (binv, binv, b->length);

        fmpz_mod_poly_divrem_newton_preinv(q, r, a, b, binv);
        fmpz_mod_poly_divrem_newton_preinv(q, binv, a, b, binv);

        result = (fmpz_mod_poly_equal(binv, r));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_mod_poly_print(a), printf("\n\n");
            fmpz_mod_poly_print(b), printf("\n\n");
            fmpz_mod_poly_print(binv), printf("\n\n");
            fmpz_mod_poly_print(q), printf("\n\n");
            fmpz_mod_poly_print(r), printf("\n\n");
            printf("n = %ld\n", n);
            abort();
        }

        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(binv);
        fmpz_mod_poly_clear(q);
        fmpz_mod_poly_clear(r);
        fmpz_clear(p);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
