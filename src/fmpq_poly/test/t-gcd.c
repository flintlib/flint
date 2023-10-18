/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_poly.h"

TEST_FUNCTION_START(fmpq_poly_gcd, state)
{
    int cflags = 0, i, result;

    /* Check aliasing of a and c */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, c;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(c);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 100);
        fmpq_poly_randtest(b, state, n_randint(state, 100), 100);

        fmpq_poly_gcd(c, a, b);
        fmpq_poly_gcd(a, a, b);

        cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(c) ? 0 : 2;
        result = (fmpq_poly_equal(a, c) && !cflags);
        if (!result)
        {
            flint_printf("FAIL (aliasing a, c):\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            fmpq_poly_debug(c), flint_printf("\n\n");
            flint_printf("cflags = %d\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(c);
    }

    /* Check aliasing of b and c */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, c;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(c);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 100);
        fmpq_poly_randtest(b, state, n_randint(state, 100), 100);

        fmpq_poly_gcd(c, a, b);
        fmpq_poly_gcd(b, a, b);

        cflags |= fmpq_poly_is_canonical(b) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(c) ? 0 : 2;
        result = (fmpq_poly_equal(b, c) && !cflags);
        if (!result)
        {
            flint_printf("FAIL (aliasing b, c):\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            fmpq_poly_debug(c), flint_printf("\n\n");
            flint_printf("cflags = %d\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(c);
    }

    /* Generic case when a, b are most likely co-prime ***********************/

    /* Verify commutativity and that c is monic */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, c;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(c);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 100);
        fmpq_poly_randtest(b, state, n_randint(state, 100), 100);

        fmpq_poly_gcd(c, a, b);
        fmpq_poly_gcd(a, b, a);

        cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(c) ? 0 : 2;
        result = (fmpq_poly_equal(a, c) && !cflags
                  && (fmpq_poly_is_zero(c) || fmpq_poly_is_monic(c)));
        if (!result)
        {
            flint_printf("FAIL (commutativity #1):\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            fmpq_poly_debug(c), flint_printf("\n\n");
            flint_printf("cflags = %d\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(c);
    }

    /* Verify that GCD(a, b) divides a, b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, c, r1, r2;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(c);
        fmpq_poly_init(r1);
        fmpq_poly_init(r2);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 100);
        fmpq_poly_randtest(b, state, n_randint(state, 100), 100);

        fmpq_poly_gcd(c, a, b);
        if (!fmpq_poly_is_zero(c))
        {
            fmpq_poly_rem(r1, a, c);
            fmpq_poly_rem(r2, b, c);
        }

        result = fmpq_poly_is_zero(r1) && fmpq_poly_is_zero(r2);
        if (!result)
        {
            flint_printf("FAIL (division #1):\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            fmpq_poly_debug(c), flint_printf("\n\n");
            flint_printf("cflags = %d\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(c);
        fmpq_poly_clear(r1);
        fmpq_poly_clear(r2);
    }

    /* Case when a, b are not co-prime ***************************************/

    /* Verify commutativity and that c is monic */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, c, t;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(c);
        fmpq_poly_init(t);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 100);
        fmpq_poly_randtest(b, state, n_randint(state, 100), 100);
        fmpq_poly_randtest(t, state, n_randint(state, 50), 20);
        fmpq_poly_mul(a, a, t);
        fmpq_poly_mul(b, b, t);

        fmpq_poly_gcd(c, a, b);
        fmpq_poly_gcd(a, b, a);

        cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(c) ? 0 : 2;
        result = (fmpq_poly_equal(a, c) && !cflags
                  && (fmpq_poly_is_zero(c) || fmpq_poly_is_monic(c)));
        if (!result)
        {
            flint_printf("FAIL (commutativity #2):\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            fmpq_poly_debug(c), flint_printf("\n\n");
            flint_printf("cflags = %d\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(c);
        fmpq_poly_clear(t);
    }

    /* Verify that GCD(a, b) divides a, b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, c, r1, r2, t;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(c);
        fmpq_poly_init(r1);
        fmpq_poly_init(r2);
        fmpq_poly_init(t);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 100);
        fmpq_poly_randtest(b, state, n_randint(state, 100), 100);
        fmpq_poly_randtest(t, state, n_randint(state, 50), 20);
        fmpq_poly_mul(a, a, t);
        fmpq_poly_mul(b, b, t);

        fmpq_poly_gcd(c, a, b);
        if (!fmpq_poly_is_zero(c))
        {
            fmpq_poly_rem(r1, a, c);
            fmpq_poly_rem(r2, b, c);
        }

        result = fmpq_poly_is_zero(r1) && fmpq_poly_is_zero(r2);
        if (!result)
        {
            flint_printf("FAIL (division #2):\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            fmpq_poly_debug(c), flint_printf("\n\n");
            flint_printf("cflags = %d\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(c);
        fmpq_poly_clear(r1);
        fmpq_poly_clear(r2);
        fmpq_poly_clear(t);
    }

    TEST_FUNCTION_END(state);
}
