/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

TEST_FUNCTION_START(fmpz_mod_poly_gcd, state)
{
    int i, result;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctx, 2);

    /* Generic case, most likely co-prime arguments ******************************/

    /* Check aliasing of a and c */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, c;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(c, ctx);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);
        fmpz_mod_poly_randtest(b, state, n_randint(state, 100), ctx);

        fmpz_mod_poly_gcd(c, a, b, ctx);
        fmpz_mod_poly_gcd(a, a, b, ctx);

        result = (fmpz_mod_poly_equal(a, c, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(c, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(c, ctx);
        fmpz_clear(p);
    }

    /* Check aliasing of b and c */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, c;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(c, ctx);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);
        fmpz_mod_poly_randtest(b, state, n_randint(state, 100), ctx);

        fmpz_mod_poly_gcd(c, a, b, ctx);
        fmpz_mod_poly_gcd(b, a, b, ctx);

        result = (fmpz_mod_poly_equal(b, c, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(c, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(c, ctx);
        fmpz_clear(p);
    }

    /*
        Check that g = GCD(a,b) divides a and b,
        and that 1 == GCD(a/g, b/g)
     */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, c, d, g, h, s, t;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(c, ctx);
        fmpz_mod_poly_init(d, ctx);
        fmpz_mod_poly_init(g, ctx);
        fmpz_mod_poly_init(h, ctx);
        fmpz_mod_poly_init(s, ctx);
        fmpz_mod_poly_init(t, ctx);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);
        fmpz_mod_poly_randtest(b, state, n_randint(state, 100), ctx);

        fmpz_mod_poly_gcd(g, a, b, ctx);

        if (fmpz_mod_poly_is_zero(g, ctx))
        {
            result = 1;
        }
        else
        {
            fmpz_mod_poly_divrem_basecase(c, s, a, g, ctx);
            fmpz_mod_poly_divrem_basecase(d, t, b, g, ctx);
            fmpz_mod_poly_gcd(h, c, d, ctx);

            result = (fmpz_mod_poly_is_zero(s, ctx) && fmpz_mod_poly_is_zero(t, ctx)
                      && (h->length == 1) && (fmpz_is_one(h->coeffs)));
        }

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("p = "), fmpz_print(p), flint_printf("\n\n");
            flint_printf("a = "), fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b = "), fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("c = "), fmpz_mod_poly_print(c, ctx), flint_printf("\n\n");
            flint_printf("d = "), fmpz_mod_poly_print(d, ctx), flint_printf("\n\n");
            flint_printf("g = "), fmpz_mod_poly_print(g, ctx), flint_printf("\n\n");
            flint_printf("h = "), fmpz_mod_poly_print(h, ctx), flint_printf("\n\n");
            flint_printf("s = "), fmpz_mod_poly_print(s, ctx), flint_printf("\n\n");
            flint_printf("t = "), fmpz_mod_poly_print(t, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(c, ctx);
        fmpz_mod_poly_clear(d, ctx);
        fmpz_mod_poly_clear(g, ctx);
        fmpz_mod_poly_clear(h, ctx);
        fmpz_mod_poly_clear(s, ctx);
        fmpz_mod_poly_clear(t, ctx);
        fmpz_clear(p);
    }

    /* Special case, arguments share a factor ********************************/

    /* Check aliasing of a and c */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, c, f;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(c, ctx);
        fmpz_mod_poly_init(f, ctx);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);
        fmpz_mod_poly_randtest(b, state, n_randint(state, 100), ctx);
        fmpz_mod_poly_randtest(f, state, n_randint(state, 20), ctx);
        fmpz_mod_poly_mul(a, a, f, ctx);
        fmpz_mod_poly_mul(b, b, f, ctx);

        fmpz_mod_poly_gcd(c, a, b, ctx);
        fmpz_mod_poly_gcd(a, a, b, ctx);

        result = (fmpz_mod_poly_equal(a, c, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(c, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(c, ctx);
        fmpz_mod_poly_clear(f, ctx);
        fmpz_clear(p);
    }

    /* Check aliasing of b and c */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, c, f;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(c, ctx);
        fmpz_mod_poly_init(f, ctx);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);
        fmpz_mod_poly_randtest(b, state, n_randint(state, 100), ctx);
        fmpz_mod_poly_randtest(f, state, n_randint(state, 20), ctx);
        fmpz_mod_poly_mul(a, a, f, ctx);
        fmpz_mod_poly_mul(b, b, f, ctx);

        fmpz_mod_poly_gcd(c, a, b, ctx);
        fmpz_mod_poly_gcd(b, a, b, ctx);

        result = (fmpz_mod_poly_equal(b, c, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(c, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(c, ctx);
        fmpz_mod_poly_clear(f, ctx);
        fmpz_clear(p);
    }

    /*
        Check that g = GCD(a,b) divides a and b,
        and that 1 == GCD(a/g, b/g)
     */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, c, d, f, g, h, s, t;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(c, ctx);
        fmpz_mod_poly_init(d, ctx);
        fmpz_mod_poly_init(f, ctx);
        fmpz_mod_poly_init(g, ctx);
        fmpz_mod_poly_init(h, ctx);
        fmpz_mod_poly_init(s, ctx);
        fmpz_mod_poly_init(t, ctx);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);
        fmpz_mod_poly_randtest(b, state, n_randint(state, 100), ctx);
        fmpz_mod_poly_randtest(f, state, n_randint(state, 20), ctx);
        fmpz_mod_poly_mul(a, a, f, ctx);
        fmpz_mod_poly_mul(b, b, f, ctx);

        fmpz_mod_poly_gcd(g, a, b, ctx);

        if (fmpz_mod_poly_is_zero(g, ctx))
        {
            result = 1;
        }
        else
        {
            fmpz_mod_poly_divrem_basecase(c, s, a, g, ctx);
            fmpz_mod_poly_divrem_basecase(d, t, b, g, ctx);
            fmpz_mod_poly_gcd(h, c, d, ctx);

            result = (fmpz_mod_poly_is_zero(s, ctx) && fmpz_mod_poly_is_zero(t, ctx)
                      && (h->length == 1) && (fmpz_is_one(h->coeffs)));
        }

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("p = "), fmpz_print(p), flint_printf("\n\n");
            flint_printf("a = "), fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b = "), fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("c = "), fmpz_mod_poly_print(c, ctx), flint_printf("\n\n");
            flint_printf("d = "), fmpz_mod_poly_print(d, ctx), flint_printf("\n\n");
            flint_printf("g = "), fmpz_mod_poly_print(g, ctx), flint_printf("\n\n");
            flint_printf("h = "), fmpz_mod_poly_print(h, ctx), flint_printf("\n\n");
            flint_printf("s = "), fmpz_mod_poly_print(s, ctx), flint_printf("\n\n");
            flint_printf("t = "), fmpz_mod_poly_print(t, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(c, ctx);
        fmpz_mod_poly_clear(d, ctx);
        fmpz_mod_poly_clear(f, ctx);
        fmpz_mod_poly_clear(g, ctx);
        fmpz_mod_poly_clear(h, ctx);
        fmpz_mod_poly_clear(s, ctx);
        fmpz_mod_poly_clear(t, ctx);
        fmpz_clear(p);
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}
