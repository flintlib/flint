/*
    Copyright (C) 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_xgcd_euclidean, state)
{
    int i, result;

    /*
       Compare with result from gcd and check a*s + b*t = g
    */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, g1, s, t, g2;

        mp_limb_t n;
        do n = n_randtest_not_zero(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_init(g1, n);
        nmod_poly_init(g2, n);
        nmod_poly_init(s, n);
        nmod_poly_init(t, n);

        nmod_poly_randtest(a, state, n_randtest(state) % 200);
        nmod_poly_randtest(b, state, n_randtest(state) % 200);
        nmod_poly_randtest(c, state, n_randtest(state) % 200);

        nmod_poly_mul(a, a, c);
        nmod_poly_mul(b, b, c);

        nmod_poly_gcd_euclidean(g1, a, b);
        nmod_poly_xgcd_euclidean(g2, s, t, a, b);

        nmod_poly_mul(s, s, a);
        nmod_poly_mul(t, t, b);
        nmod_poly_add(s, s, t);

        result = (nmod_poly_equal(g1, g2) && nmod_poly_equal(s, g1)
            && (g1->length == 0 || g1->coeffs[g1->length - 1] == 1));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(c), flint_printf("\n\n");
            nmod_poly_print(g1), flint_printf("\n\n");
            nmod_poly_print(g2), flint_printf("\n\n");
            nmod_poly_print(s), flint_printf("\n\n");
            nmod_poly_print(t), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(g1);
        nmod_poly_clear(g2);
        nmod_poly_clear(s);
        nmod_poly_clear(t);
    }

    /* Check aliasing of a and g */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, g, s, t;

        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(g, n);
        nmod_poly_init(s, n);
        nmod_poly_init(t, n);
        nmod_poly_randtest(a, state, n_randtest(state) % 200);
        nmod_poly_randtest(b, state, n_randtest(state) % 200);

        nmod_poly_xgcd_euclidean(g, s, t, a, b);
        nmod_poly_xgcd_euclidean(a, s, t, a, b);

        result = (nmod_poly_equal(a, g));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(g), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(g);
        nmod_poly_clear(s);
        nmod_poly_clear(t);
    }

    /* Check aliasing of b and g */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, g, s, t;

        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(g, n);
        nmod_poly_init(s, n);
        nmod_poly_init(t, n);
        nmod_poly_randtest(a, state, n_randtest(state) % 200);
        nmod_poly_randtest(b, state, n_randtest(state) % 200);

        nmod_poly_xgcd_euclidean(g, s, t, a, b);
        nmod_poly_xgcd_euclidean(b, s, t, a, b);

        result = (nmod_poly_equal(b, g));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(g), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(g);
        nmod_poly_clear(s);
        nmod_poly_clear(t);
    }

    /* Check aliasing of s and a */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, g, s, t;

        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(g, n);
        nmod_poly_init(s, n);
        nmod_poly_init(t, n);
        nmod_poly_randtest(a, state, n_randtest(state) % 200);
        nmod_poly_randtest(b, state, n_randtest(state) % 200);

        nmod_poly_xgcd_euclidean(g, s, t, a, b);
        nmod_poly_xgcd_euclidean(g, a, t, a, b);

        result = (nmod_poly_equal(s, a));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(s), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(g);
        nmod_poly_clear(s);
        nmod_poly_clear(t);
    }

    /* Check aliasing of s and b */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, g, s, t;

        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(g, n);
        nmod_poly_init(s, n);
        nmod_poly_init(t, n);
        nmod_poly_randtest(a, state, n_randtest(state) % 200);
        nmod_poly_randtest(b, state, n_randtest(state) % 200);

        nmod_poly_xgcd_euclidean(g, s, t, a, b);
        nmod_poly_xgcd_euclidean(g, b, t, a, b);

        result = (nmod_poly_equal(s, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(s), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(g);
        nmod_poly_clear(s);
        nmod_poly_clear(t);
    }

    /* Check aliasing of t and a */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, g, s, t;

        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(g, n);
        nmod_poly_init(s, n);
        nmod_poly_init(t, n);
        nmod_poly_randtest(a, state, n_randtest(state) % 200);
        nmod_poly_randtest(b, state, n_randtest(state) % 200);

        nmod_poly_xgcd_euclidean(g, s, t, a, b);
        nmod_poly_xgcd_euclidean(g, s, a, a, b);

        result = (nmod_poly_equal(t, a));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(t), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(g);
        nmod_poly_clear(s);
        nmod_poly_clear(t);
    }

    /* Check aliasing of t and b */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, g, s, t;

        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(g, n);
        nmod_poly_init(s, n);
        nmod_poly_init(t, n);
        nmod_poly_randtest(a, state, n_randtest(state) % 200);
        nmod_poly_randtest(b, state, n_randtest(state) % 200);

        nmod_poly_xgcd_euclidean(g, s, t, a, b);
        nmod_poly_xgcd_euclidean(g, s, b, a, b);

        result = (nmod_poly_equal(t, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(t), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(g);
        nmod_poly_clear(s);
        nmod_poly_clear(t);
    }

    TEST_FUNCTION_END(state);
}
