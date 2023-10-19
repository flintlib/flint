/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2013 Martin Lee

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

TEST_FUNCTION_START(fmpz_mod_poly_divrem_newton_n_preinv, state)
{
    int i, result;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctx, 2);

    /* Check q*b + r = a, no aliasing */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, binv, q, r, t, q2, r2;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);
        fmpz_mod_ctx_set_modulus(ctx, p);
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(binv, ctx);
        fmpz_mod_poly_init(q, ctx);
        fmpz_mod_poly_init(q2, ctx);
        fmpz_mod_poly_init(r, ctx);
        fmpz_mod_poly_init(r2, ctx);
        fmpz_mod_poly_init(t, ctx);
        do
        fmpz_mod_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, ctx);
        while (b->length <= 2);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);
        if (a->length > 2*(b->length)-3)
          fmpz_mod_poly_truncate (a, 2*(b->length)-3, ctx);

        {
            fmpz_t d;
            fmpz *leadB = fmpz_mod_poly_lead(b, ctx);

            fmpz_init(d);
            fmpz_gcd(d, p, leadB);
            while (!fmpz_is_one(d))
            {
                fmpz_divexact(leadB, leadB, d);
                fmpz_gcd(d, p, leadB);
            }
            fmpz_clear(d);
        }

        fmpz_mod_poly_reverse (binv, b, b->length, ctx);
        fmpz_mod_poly_inv_series (binv, binv, b->length, ctx);

        fmpz_mod_poly_divrem_newton_n_preinv (q, r, a, b, binv, ctx);
        fmpz_mod_poly_mul(t, q, b, ctx);
        fmpz_mod_poly_add(t, t, r, ctx);

        result = (fmpz_mod_poly_equal(a, t, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("p = "), fmpz_print(p), flint_printf("\n\n");
            flint_printf("a = "), fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b = "), fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("q = "), fmpz_mod_poly_print(q, ctx), flint_printf("\n\n");
            flint_printf("r = "), fmpz_mod_poly_print(r, ctx), flint_printf("\n\n");
            flint_printf("t = "), fmpz_mod_poly_print(t, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_divrem_basecase(q2, r2, a, b, ctx);
        result = (fmpz_mod_poly_equal(q2, q, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("p = "), fmpz_print(p), flint_printf("\n\n");
            flint_printf("a = "), fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b = "), fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("q = "), fmpz_mod_poly_print(q, ctx), flint_printf("\n\n");
            flint_printf("q2 = "), fmpz_mod_poly_print(q2, ctx), flint_printf("\n\n");
            flint_printf("r = "), fmpz_mod_poly_print(r, ctx), flint_printf("\n\n");
            flint_printf("t = "), fmpz_mod_poly_print(t, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(binv, ctx);
        fmpz_mod_poly_clear(q, ctx);
        fmpz_mod_poly_clear(q2, ctx);
        fmpz_mod_poly_clear(r, ctx);
        fmpz_mod_poly_clear(r2, ctx);
        fmpz_mod_poly_clear(t, ctx);
        fmpz_clear(p);
    }

    /* Alias a and q, b and r */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, binv, q, r;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(binv, ctx);
        fmpz_mod_poly_init(q, ctx);
        fmpz_mod_poly_init(r, ctx);
        do
        fmpz_mod_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, ctx);
        while (b->length <= 2);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);
        if (a->length > 2*(b->length)-3)
          fmpz_mod_poly_truncate (a, 2*(b->length)-3, ctx);

        {
            fmpz_t d;
            fmpz *leadB = fmpz_mod_poly_lead(b, ctx);

            fmpz_init(d);
            fmpz_gcd(d, p, leadB);
            while (!fmpz_is_one(d))
            {
                fmpz_divexact(leadB, leadB, d);
                fmpz_gcd(d, p, leadB);
            }
            fmpz_clear(d);
        }

        fmpz_mod_poly_reverse (binv, b, b->length, ctx);
        fmpz_mod_poly_inv_series (binv, binv, b->length, ctx);

        fmpz_mod_poly_divrem_newton_n_preinv(q, r, a, b, binv, ctx);
        fmpz_mod_poly_divrem_newton_n_preinv(a, b, a, b, binv, ctx);

        result = (fmpz_mod_poly_equal(q, a, ctx) && fmpz_mod_poly_equal(r, b, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("p = "), fmpz_print(p), flint_printf("\n\n");
            flint_printf("a = "), fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b = "), fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("q = "), fmpz_mod_poly_print(q, ctx), flint_printf("\n\n");
            flint_printf("r = "), fmpz_mod_poly_print(r, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(binv, ctx);
        fmpz_mod_poly_clear(q, ctx);
        fmpz_mod_poly_clear(r, ctx);
        fmpz_clear(p);
    }

    /* Alias b and q, a and r */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, binv, q, r;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(binv, ctx);
        fmpz_mod_poly_init(q, ctx);
        fmpz_mod_poly_init(r, ctx);
        do
        fmpz_mod_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, ctx);
        while (b->length <= 2);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);
        if (a->length > 2*(b->length)-3)
          fmpz_mod_poly_truncate (a, 2*(b->length)-3, ctx);

        {
            fmpz_t d;
            fmpz *leadB = fmpz_mod_poly_lead(b, ctx);

            fmpz_init(d);
            fmpz_gcd(d, p, leadB);
            while (!fmpz_is_one(d))
            {
                fmpz_divexact(leadB, leadB, d);
                fmpz_gcd(d, p, leadB);
            }
            fmpz_clear(d);
        }

        fmpz_mod_poly_reverse (binv, b, b->length, ctx);
        fmpz_mod_poly_inv_series (binv, binv, b->length, ctx);

        fmpz_mod_poly_divrem_newton_n_preinv(q, r, a, b, binv, ctx);
        fmpz_mod_poly_divrem_newton_n_preinv(b, a, a, b, binv, ctx);

        result = (fmpz_mod_poly_equal(q, b, ctx) && fmpz_mod_poly_equal(r, a, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("p = "), fmpz_print(p), flint_printf("\n\n");
            flint_printf("a = "), fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b = "), fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("q = "), fmpz_mod_poly_print(q, ctx), flint_printf("\n\n");
            flint_printf("r = "), fmpz_mod_poly_print(r, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(binv, ctx);
        fmpz_mod_poly_clear(q, ctx);
        fmpz_mod_poly_clear(r, ctx);
        fmpz_clear(p);
    }

    /* Alias binv and q, a and r*/
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, binv, q, r;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(binv, ctx);
        fmpz_mod_poly_init(q, ctx);
        fmpz_mod_poly_init(r, ctx);
        do
        fmpz_mod_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, ctx);
        while (b->length <= 2);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);
        if (a->length > 2*(b->length)-3)
          fmpz_mod_poly_truncate (a, 2*(b->length)-3, ctx);

        {
            fmpz_t d;
            fmpz *leadB = fmpz_mod_poly_lead(b, ctx);

            fmpz_init(d);
            fmpz_gcd(d, p, leadB);
            while (!fmpz_is_one(d))
            {
                fmpz_divexact(leadB, leadB, d);
                fmpz_gcd(d, p, leadB);
            }
            fmpz_clear(d);
        }

        fmpz_mod_poly_reverse (binv, b, b->length, ctx);
        fmpz_mod_poly_inv_series (binv, binv, b->length, ctx);

        fmpz_mod_poly_divrem_newton_n_preinv(q, r, a, b, binv, ctx);
        fmpz_mod_poly_divrem_newton_n_preinv(binv, a, a, b, binv, ctx);

        result = (fmpz_mod_poly_equal(q, binv, ctx) && fmpz_mod_poly_equal(r, a, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("p = "), fmpz_print(p), flint_printf("\n\n");
            flint_printf("a = "), fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b = "), fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("binv = "), fmpz_mod_poly_print(binv, ctx), flint_printf("\n\n");
            flint_printf("q = "), fmpz_mod_poly_print(q, ctx), flint_printf("\n\n");
            flint_printf("r = "), fmpz_mod_poly_print(r, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(binv, ctx);
        fmpz_mod_poly_clear(q, ctx);
        fmpz_mod_poly_clear(r, ctx);
        fmpz_clear(p);
    }

    /* Alias binv and q, b and r*/
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, binv, q, r;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(binv, ctx);
        fmpz_mod_poly_init(q, ctx);
        fmpz_mod_poly_init(r, ctx);
        do
        fmpz_mod_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, ctx);
        while (b->length <= 2);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);
        if (a->length > 2*(b->length)-3)
          fmpz_mod_poly_truncate (a, 2*(b->length)-3, ctx);

        {
            fmpz_t d;
            fmpz *leadB = fmpz_mod_poly_lead(b, ctx);

            fmpz_init(d);
            fmpz_gcd(d, p, leadB);
            while (!fmpz_is_one(d))
            {
                fmpz_divexact(leadB, leadB, d);
                fmpz_gcd(d, p, leadB);
            }
            fmpz_clear(d);
        }

        fmpz_mod_poly_reverse (binv, b, b->length, ctx);
        fmpz_mod_poly_inv_series (binv, binv, b->length, ctx);

        fmpz_mod_poly_divrem_newton_n_preinv(q, r, a, b, binv, ctx);
        fmpz_mod_poly_divrem_newton_n_preinv(binv, b, a, b, binv, ctx);

        result = (fmpz_mod_poly_equal(q, binv, ctx) && fmpz_mod_poly_equal(r, b, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("p = "), fmpz_print(p), flint_printf("\n\n");
            flint_printf("a = "), fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b = "), fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("binv = "), fmpz_mod_poly_print(binv, ctx), flint_printf("\n\n");
            flint_printf("q = "), fmpz_mod_poly_print(q, ctx), flint_printf("\n\n");
            flint_printf("r = "), fmpz_mod_poly_print(r, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(binv, ctx);
        fmpz_mod_poly_clear(q, ctx);
        fmpz_mod_poly_clear(r, ctx);
        fmpz_clear(p);
    }

    /* Alias a and q, binv and r*/
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, binv, q, r;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(binv, ctx);
        fmpz_mod_poly_init(q, ctx);
        fmpz_mod_poly_init(r, ctx);
        do
        fmpz_mod_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, ctx);
        while (b->length <= 2);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);
        if (a->length > 2*(b->length)-3)
          fmpz_mod_poly_truncate (a, 2*(b->length)-3, ctx);

        {
            fmpz_t d;
            fmpz *leadB = fmpz_mod_poly_lead(b, ctx);

            fmpz_init(d);
            fmpz_gcd(d, p, leadB);
            while (!fmpz_is_one(d))
            {
                fmpz_divexact(leadB, leadB, d);
                fmpz_gcd(d, p, leadB);
            }
            fmpz_clear(d);
        }

        fmpz_mod_poly_reverse (binv, b, b->length, ctx);
        fmpz_mod_poly_inv_series (binv, binv, b->length, ctx);

        fmpz_mod_poly_divrem_newton_n_preinv(q, r, a, b, binv, ctx);
        fmpz_mod_poly_divrem_newton_n_preinv(a, binv, a, b, binv, ctx);

        result = (fmpz_mod_poly_equal(q, a, ctx) && fmpz_mod_poly_equal(r, binv, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("p = "), fmpz_print(p), flint_printf("\n\n");
            flint_printf("a = "), fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b = "), fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("binv = "), fmpz_mod_poly_print(binv, ctx), flint_printf("\n\n");
            flint_printf("q = "), fmpz_mod_poly_print(q, ctx), flint_printf("\n\n");
            flint_printf("r = "), fmpz_mod_poly_print(r, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(binv, ctx);
        fmpz_mod_poly_clear(q, ctx);
        fmpz_mod_poly_clear(r, ctx);
        fmpz_clear(p);
    }

    /* Alias b and q, binv and r*/
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, binv, q, r;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(binv, ctx);
        fmpz_mod_poly_init(q, ctx);
        fmpz_mod_poly_init(r, ctx);
        do
        fmpz_mod_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, ctx);
        while (b->length <= 2);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);
        if (a->length > 2*(b->length)-3)
          fmpz_mod_poly_truncate (a, 2*(b->length)-3, ctx);

        {
            fmpz_t d;
            fmpz *leadB = fmpz_mod_poly_lead(b, ctx);

            fmpz_init(d);
            fmpz_gcd(d, p, leadB);
            while (!fmpz_is_one(d))
            {
                fmpz_divexact(leadB, leadB, d);
                fmpz_gcd(d, p, leadB);
            }
            fmpz_clear(d);
        }

        fmpz_mod_poly_reverse (binv, b, b->length, ctx);
        fmpz_mod_poly_inv_series (binv, binv, b->length, ctx);

        fmpz_mod_poly_divrem_newton_n_preinv(q, r, a, b, binv, ctx);
        fmpz_mod_poly_divrem_newton_n_preinv(b, binv, a, b, binv, ctx);

        result = (fmpz_mod_poly_equal(q, b, ctx) && fmpz_mod_poly_equal(r, binv, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("p = "), fmpz_print(p), flint_printf("\n\n");
            flint_printf("a = "), fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b = "), fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("binv = "), fmpz_mod_poly_print(binv, ctx), flint_printf("\n\n");
            flint_printf("q = "), fmpz_mod_poly_print(q, ctx), flint_printf("\n\n");
            flint_printf("r = "), fmpz_mod_poly_print(r, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(binv, ctx);
        fmpz_mod_poly_clear(q, ctx);
        fmpz_mod_poly_clear(r, ctx);
        fmpz_clear(p);
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}
