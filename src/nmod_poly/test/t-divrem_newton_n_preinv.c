/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2013 Martin Lee

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_divrem_newton_n_preinv, state)
{
    int i, result;

    /* Check result of divrem */
    for (i = 0; i < 2000 * flint_test_multiplier(); i++)
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
        nmod_poly_randtest(a, state, 2*b->length - 1);

        nmod_poly_reverse (binv, b, b->length);
        nmod_poly_inv_series (binv, binv, b->length);
        nmod_poly_divrem_newton_n_preinv(q, r, a, b, binv);
        nmod_poly_mul(test, q, b);
        nmod_poly_add(test, test, r);

        result = (nmod_poly_equal(a, test));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a.length: %wd\n", a->length);
            flint_printf("b.length: %wd\n", b->length);
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(test), flint_printf("\n\n");
            nmod_poly_print(q), flint_printf("\n\n");
            nmod_poly_print(r), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(binv);
        nmod_poly_clear(q);
        nmod_poly_clear(r);
        nmod_poly_clear(test);
    }
#if 0
    /* Check aliasing of a and q */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, binv, q, r;

        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(binv, n);
        nmod_poly_init(q, n);
        nmod_poly_init(r, n);
        do
        nmod_poly_randtest(b, state, n_randint(state, 200));
        while (b->length <= 2);
        nmod_poly_randtest(a, state, n_randint(state, 200));
        if (a->length > 2*(b->length)-3)
          nmod_poly_truncate (a, 2*(b->length)-3);

        nmod_poly_reverse (binv, b, b->length);
        nmod_poly_inv_series (binv, binv, b->length);

        nmod_poly_divrem_newton_n_preinv(q, r, a, b, binv);
        nmod_poly_divrem_newton_n_preinv(a, r, a, b, binv);

        result = (nmod_poly_equal(a, q));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(q), flint_printf("\n\n");
            nmod_poly_print(r), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(binv);
        nmod_poly_clear(q);
        nmod_poly_clear(r);
    }

    /* Check aliasing of b and q */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, binv, q, r;

        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(binv, n);
        nmod_poly_init(q, n);
        nmod_poly_init(r, n);
        do
        nmod_poly_randtest(b, state, n_randint(state, 200));
        while (b->length <= 2);
        nmod_poly_randtest(a, state, n_randint(state, 200));
        if (a->length > 2*(b->length)-3)
          nmod_poly_truncate (a, 2*(b->length)-3);

        nmod_poly_reverse (binv, b, b->length);
        nmod_poly_inv_series (binv, binv, b->length);

        nmod_poly_divrem_newton_n_preinv(q, r, a, b, binv);
        nmod_poly_divrem_newton_n_preinv(b, r, a, b, binv);

        result = (nmod_poly_equal(b, q));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(q), flint_printf("\n\n");
            nmod_poly_print(r), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(binv);
        nmod_poly_clear(q);
        nmod_poly_clear(r);
    }

    /* Check aliasing of binv and q */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, binv, q, r;

        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(binv, n);
        nmod_poly_init(q, n);
        nmod_poly_init(r, n);
        do
        nmod_poly_randtest(b, state, n_randint(state, 200));
        while (b->length <= 2);
        nmod_poly_randtest(a, state, n_randint(state, 200));
        if (a->length > 2*(b->length)-3)
          nmod_poly_truncate (a, 2*(b->length)-3);

        nmod_poly_reverse (binv, b, b->length);
        nmod_poly_inv_series (binv, binv, b->length);

        nmod_poly_divrem_newton_n_preinv(q, r, a, b, binv);
        nmod_poly_divrem_newton_n_preinv(binv, r, a, b, binv);

        result = (nmod_poly_equal(binv, q));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(binv), flint_printf("\n\n");
            nmod_poly_print(q), flint_printf("\n\n");
            nmod_poly_print(r), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(binv);
        nmod_poly_clear(q);
        nmod_poly_clear(r);
    }

    /* Check aliasing of a and r */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, binv, q, r;

        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(binv, n);
        nmod_poly_init(q, n);
        nmod_poly_init(r, n);
        do
        nmod_poly_randtest(b, state, n_randint(state, 200));
        while (b->length <= 2);
        nmod_poly_randtest(a, state, n_randint(state, 200));
        if (a->length > 2*(b->length)-3)
          nmod_poly_truncate (a, 2*(b->length)-3);

        nmod_poly_reverse (binv, b, b->length);
        nmod_poly_inv_series (binv, binv, b->length);

        nmod_poly_divrem_newton_n_preinv(q, r, a, b, binv);
        nmod_poly_divrem_newton_n_preinv(q, a, a, b, binv);

        result = (nmod_poly_equal(a, r));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(q), flint_printf("\n\n");
            nmod_poly_print(r), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(binv);
        nmod_poly_clear(q);
        nmod_poly_clear(r);
    }

    /* Check aliasing of b and r */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, binv, q, r;

        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(binv, n);
        nmod_poly_init(q, n);
        nmod_poly_init(r, n);
        do
        nmod_poly_randtest(b, state, n_randint(state, 200));
        while (b->length <= 2);
        nmod_poly_randtest(a, state, n_randint(state, 200));
        if (a->length > 2*(b->length)-3)
          nmod_poly_truncate (a, 2*(b->length)-3);

        nmod_poly_reverse (binv, b, b->length);
        nmod_poly_inv_series (binv, binv, b->length);

        nmod_poly_divrem_newton_n_preinv(q, r, a, b, binv);
        nmod_poly_divrem_newton_n_preinv(q, b, a, b, binv);

        result = (nmod_poly_equal(b, r));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(q), flint_printf("\n\n");
            nmod_poly_print(r), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(binv);
        nmod_poly_clear(q);
        nmod_poly_clear(r);
    }

    /* Check aliasing of binv and r */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, binv, q, r;

        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(binv, n);
        nmod_poly_init(q, n);
        nmod_poly_init(r, n);
        do
        nmod_poly_randtest(b, state, n_randint(state, 200));
        while (b->length <= 2);
        nmod_poly_randtest(a, state, n_randint(state, 200));
        if (a->length > 2*(b->length)-3)
          nmod_poly_truncate (a, 2*(b->length)-3);

        nmod_poly_reverse (binv, b, b->length);
        nmod_poly_inv_series (binv, binv, b->length);

        nmod_poly_divrem_newton_n_preinv(q, r, a, b, binv);
        nmod_poly_divrem_newton_n_preinv(q, binv, a, b, binv);

        result = (nmod_poly_equal(binv, r));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(binv), flint_printf("\n\n");
            nmod_poly_print(q), flint_printf("\n\n");
            nmod_poly_print(r), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(binv);
        nmod_poly_clear(q);
        nmod_poly_clear(r);
    }
#endif

    TEST_FUNCTION_END(state);
}
