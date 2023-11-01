/*
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_poly.h"

TEST_FUNCTION_START(fmpq_poly_remove, state)
{
    int i, result;
    ulong cflags = UWORD(0);

    /* Check aliasing of q and a */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, p, q;
        slong e1, e2;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(p);
        fmpq_poly_init(q);

        e1 = n_randint(state, 5);

        do {
            fmpq_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, 100);
        } while (b->length < 2);
        do {
            fmpq_poly_randtest(a, state, n_randint(state, 50), 100);
            fmpq_poly_gcd(q, a, b);
        } while (!fmpq_poly_is_one(q));

        fmpq_poly_pow(p, b, e1);
        fmpq_poly_mul(p, p, a);

        e1 = fmpq_poly_remove(q, p, b);
        e2 = fmpq_poly_remove(p, p, b);

        cflags |= fmpq_poly_is_canonical(q) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(p) ? 0 : 2;
        result = (e1 == e2 && fmpq_poly_equal(q, p) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("q = "), fmpq_poly_debug(q), flint_printf("\n\n");
            flint_printf("a = "), fmpq_poly_debug(a), flint_printf("\n\n");
            flint_printf("b = "), fmpq_poly_debug(b), flint_printf("\n\n");
            flint_printf("p = "), fmpq_poly_debug(p), flint_printf("\n\n");
            flint_printf("cflags = %wu, e1 = %wd, e2 = %wd\n\n", cflags, e1, e2);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(p);
        fmpq_poly_clear(q);
    }

    /* Check aliasing of q and b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, p, q;
        slong e1, e2;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(p);
        fmpq_poly_init(q);

        e1 = n_randint(state, 5);

        do {
            fmpq_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, 100);
        } while (b->length < 2);
        do {
            fmpq_poly_randtest(a, state, n_randint(state, 50), 100);
            fmpq_poly_gcd(q, a, b);
        } while (!fmpq_poly_is_one(q));

        fmpq_poly_pow(p, b, e1);
        fmpq_poly_mul(p, p, a);

        e1 = fmpq_poly_remove(q, p, b);
        e2 = fmpq_poly_remove(b, p, b);

        cflags |= fmpq_poly_is_canonical(q) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(b) ? 0 : 2;
        result = (e1 == e2 && fmpq_poly_equal(q, b) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("q = "), fmpq_poly_debug(q), flint_printf("\n\n");
            flint_printf("a = "), fmpq_poly_debug(a), flint_printf("\n\n");
            flint_printf("b = "), fmpq_poly_debug(b), flint_printf("\n\n");
            flint_printf("p = "), fmpq_poly_debug(p), flint_printf("\n\n");
            flint_printf("cflags = %wu, e1 = %wd, e2 = %wd\n\n", cflags, e1, e2);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(p);
        fmpq_poly_clear(q);
    }

    /* Compare with construction */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, p, q;
        slong e1, e2;
        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(p);
        fmpq_poly_init(q);

        e1 = n_randint(state, 5);

        do {
            fmpq_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, 100);
        } while (b->length < 2);
        do {
            fmpq_poly_randtest(a, state, n_randint(state, 50), 100);
            fmpq_poly_gcd(q, a, b);
        } while (!fmpq_poly_is_one(q));

        fmpq_poly_pow(p, b, e1);
        fmpq_poly_mul(p, p, a);

        e2 = fmpq_poly_remove(q, p, b);

        cflags |= fmpq_poly_is_canonical(q) ? 0 : 1;
        result = (e1 == e2 && fmpq_poly_equal(q, a) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("q = "), fmpq_poly_debug(q), flint_printf("\n\n");
            flint_printf("a = "), fmpq_poly_debug(a), flint_printf("\n\n");
            flint_printf("b = "), fmpq_poly_debug(b), flint_printf("\n\n");
            flint_printf("p = "), fmpq_poly_debug(p), flint_printf("\n\n");
            flint_printf("cflags = %wu, e1 = %wd, e2 = %wd\n\n", cflags, e1, e2);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(p);
        fmpq_poly_clear(q);
    }

    TEST_FUNCTION_END(state);
}
