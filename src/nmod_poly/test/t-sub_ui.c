/*
    Copyright (C) 2026 Pavlos Athinaios

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_sub_ui, state)
{
    int i, result;

    /* Check a - x == a + (-x), with c possibly >= n */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, c, d, e;
        ulong n = n_randtest_not_zero(state);
        ulong x = n_randtest(state);

        nmod_poly_init(a, n);
        nmod_poly_init(c, n);
        nmod_poly_init(d, n);
        nmod_poly_init(e, n);
        nmod_poly_randtest(a, state, n_randint(state, 100));

        /* Reference: constant polynomial (n - (x mod n)) mod n */
        nmod_poly_set_coeff_ui(c, 0, x);   /* set_coeff_ui reduces x mod n */
        nmod_poly_neg(c, c);              /* negation gives -x */

        nmod_poly_add(d, a, c);          /* a + (-x) */
        nmod_poly_sub_ui(e, a, x);       /* a - x */

        result = (nmod_poly_equal(d, e));
        if (!result)
        {
            flint_printf("FAIL (compare against a + (-x)):\n");
            flint_printf("n = %wu, x = %wu\n", n, x);
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(c), flint_printf("\n\n");
            nmod_poly_print(d), flint_printf("\n\n");
            nmod_poly_print(e), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(c);
        nmod_poly_clear(d);
        nmod_poly_clear(e);
    }

    /* Check aliasing of res and poly */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c;
        ulong n = n_randtest_not_zero(state);
        ulong x = (n == 1) ? 0 : n_randint(state, n);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_randtest(a, state, n_randint(state, 100));
        nmod_poly_set(b, a);

        nmod_poly_sub_ui(c, a, x);
        nmod_poly_sub_ui(a, a, x);

        result = (nmod_poly_equal(a, c));
        if (!result)
        {
            flint_printf("FAIL (aliasing res == poly):\n");
            flint_printf("n = %wu, x = %wu\n", n, x);
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(c), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
    }

    /* Check empty poly with x >= n reduces to (n - (x mod n)) */
    {
        nmod_poly_t a, b, c;
        ulong n = 1000;
        ulong x = 2 * n + 7;
        ulong expected;

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);

        /* x reduces to 7, so (n - 7) = 993 */
        expected = n - (x % n);

        nmod_poly_sub_ui(b, a, x);
        nmod_poly_set_coeff_ui(c, 0, expected);

        result = (nmod_poly_equal(b, c));
        if (!result)
        {
            flint_printf("FAIL (empty poly, x = %wu mod n = %wu):\n", x, n);
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(c), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
    }

    /* Check x == n reduces to zero (so result is a copy of a) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c;
        ulong n = n_randtest_not_zero(state);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_randtest(a, state, n_randint(state, 100));
        nmod_poly_set(c, a);

        nmod_poly_sub_ui(b, a, n);

        result = (nmod_poly_equal(b, c));
        if (!result)
        {
            flint_printf("FAIL (x == n should reduce to 0):\n");
            flint_printf("n = %wu\n", n);
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(c), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
    }

    /* Check a - a[0] produces a zero polynomial when a is constant */
    {
        nmod_poly_t a, b, zero;

        nmod_poly_init(a, 100);
        nmod_poly_init(b, 100);
        nmod_poly_init(zero, 100);

        nmod_poly_set_coeff_ui(a, 0, 13);
        nmod_poly_zero(zero);

        nmod_poly_sub_ui(b, a, 13);

        result = (nmod_poly_equal(b, zero));
        if (!result)
        {
            flint_printf("FAIL (constant x - constant x should be zero):\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(zero);
    }

    TEST_FUNCTION_END(state);
}
