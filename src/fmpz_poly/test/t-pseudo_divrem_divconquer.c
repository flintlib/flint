/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpz_poly_pseudo_divrem_divconquer, state)
{
    int i, result;

    /* Check q*b + r = l(b)^d a, no aliasing */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, q, r, prod;
        fmpz_t p;
        ulong d;

        fmpz_init(p);
        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(q);
        fmpz_poly_init(r);
        fmpz_poly_init(prod);
        fmpz_poly_randtest(a, state, n_randint(state, 200), 50);
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, 50);

        fmpz_poly_pseudo_divrem_divconquer(q, r, &d, a, b);
        fmpz_poly_mul(prod, q, b);
        fmpz_poly_add(prod, prod, r);
        fmpz_pow_ui(p, b->coeffs + b->length - 1, d);
        fmpz_poly_scalar_mul_fmpz(a, a, p);

        result = (fmpz_poly_equal(a, prod));
        if (!result)
        {
            flint_printf("FAIL (check qb + r = l(b)^d a):\n");
            flint_printf("l^d a = "), fmpz_poly_print(a), flint_printf("\n\n");
            flint_printf("qb + r = "), fmpz_poly_print(prod), flint_printf("\n\n");
            flint_printf("q = "), fmpz_poly_print(q), flint_printf("\n\n");
            flint_printf("r = "), fmpz_poly_print(r), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(p);
        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(q);
        fmpz_poly_clear(r);
        fmpz_poly_clear(prod);
    }

    /* Check r and a alias */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, q, r;
        ulong d;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(q);
        fmpz_poly_init(r);
        fmpz_poly_randtest(a, state, n_randint(state, 100), 50);
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, 50);

        fmpz_poly_pseudo_divrem_divconquer(q, r, &d, a, b);
        fmpz_poly_pseudo_divrem_divconquer(q, a, &d, a, b);

        result = (fmpz_poly_equal(a, r));
        if (!result)
        {
            flint_printf("FAIL (alias r and a):\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(q), flint_printf("\n\n");
            fmpz_poly_print(r), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(q);
        fmpz_poly_clear(r);
    }

    /* Check r and b alias */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, q, r;
        ulong d;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(q);
        fmpz_poly_init(r);
        fmpz_poly_randtest(a, state, n_randint(state, 100), 50);
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, 50);

        fmpz_poly_pseudo_divrem_divconquer(q, r, &d, a, b);
        fmpz_poly_pseudo_divrem_divconquer(q, b, &d, a, b);

        result = (fmpz_poly_equal(b, r));
        if (!result)
        {
            flint_printf("FAIL (alias r and b):\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(q), flint_printf("\n\n");
            fmpz_poly_print(r), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(q);
        fmpz_poly_clear(r);
    }

    /* Check q and a alias */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, q, r;
        ulong d;
        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(q);
        fmpz_poly_init(r);
        fmpz_poly_randtest(a, state, n_randint(state, 100), 50);
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, 50);

        fmpz_poly_pseudo_divrem_divconquer(q, r, &d, a, b);
        fmpz_poly_pseudo_divrem_divconquer(a, r, &d, a, b);

        result = (fmpz_poly_equal(a, q));
        if (!result)
        {
            flint_printf("FAIL (alias q and a):\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(q), flint_printf("\n\n");
            fmpz_poly_print(r), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(q);
        fmpz_poly_clear(r);
    }

    /* Check q and b alias */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, q, r;
        ulong d;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(q);
        fmpz_poly_init(r);
        fmpz_poly_randtest(a, state, n_randint(state, 100), 50);
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, 50);

        fmpz_poly_pseudo_divrem_divconquer(q, r, &d, a, b);
        fmpz_poly_pseudo_divrem_divconquer(b, r, &d, a, b);

        result = (fmpz_poly_equal(b, q));
        if (!result)
        {
            flint_printf("FAIL (alias q and b):\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(q), flint_printf("\n\n");
            fmpz_poly_print(r), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(q);
        fmpz_poly_clear(r);
    }

    TEST_FUNCTION_END(state);
}
