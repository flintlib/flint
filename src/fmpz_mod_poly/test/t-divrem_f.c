/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz

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

TEST_FUNCTION_START(fmpz_mod_poly_divrem_f, state)
{
    int i, result;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctx, 2);

    /* Check q*b + r = a when gcd(lead(B),p) = 1, no aliasing */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        fmpz_t f, p;
        fmpz_mod_poly_t a, b, q, r, t;

        fmpz_init(f);
        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(q, ctx);
        fmpz_mod_poly_init(r, ctx);
        fmpz_mod_poly_init(t, ctx);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);
        fmpz_mod_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, ctx);

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

        fmpz_mod_poly_divrem_f(f, q, r, a, b, ctx);
        fmpz_mod_poly_mul(t, q, b, ctx);
        fmpz_mod_poly_add(t, t, r, ctx);

        result = (fmpz_is_one(f) && fmpz_mod_poly_equal(a, t, ctx));
        if (!result)
        {
            flint_printf("FAIL (divrem):\n");
            flint_printf("p = "), fmpz_print(p), flint_printf("\n\n");
            flint_printf("f = "), fmpz_print(f), flint_printf("\n\n");
            flint_printf("a = "), fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b = "), fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("q = "), fmpz_mod_poly_print(q, ctx), flint_printf("\n\n");
            flint_printf("r = "), fmpz_mod_poly_print(r, ctx), flint_printf("\n\n");
            flint_printf("t = "), fmpz_mod_poly_print(t, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(q, ctx);
        fmpz_mod_poly_clear(r, ctx);
        fmpz_mod_poly_clear(t, ctx);
        fmpz_clear(f);
        fmpz_clear(p);
    }

    /* Check f | p when gcd(lead(B),p) > 1 */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        fmpz_t f, p, q1, q2;
        fmpz_mod_poly_t a, b, q, r, t;

        fmpz_init(f);
        fmpz_init(p);
        fmpz_init(q1);
        fmpz_init(q2);
        fmpz_randtest_unsigned(q1, state, 2 * FLINT_BITS);
        fmpz_randtest_unsigned(q2, state, 2 * FLINT_BITS);
        fmpz_add_ui(q1, q1, 2);
        fmpz_add_ui(q2, q2, 2);
        fmpz_mul(p, q1, q2);
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(q, ctx);
        fmpz_mod_poly_init(r, ctx);
        fmpz_mod_poly_init(t, ctx);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);
        fmpz_mod_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, ctx);

        {
            fmpz_t d;
            fmpz *leadB = fmpz_mod_poly_lead(b, ctx);

            fmpz_init(d);
            fmpz_gcd(d, p, leadB);
            if (fmpz_is_one(d))
                fmpz_set(leadB, q1);
            fmpz_clear(d);
        }

        fmpz_mod_poly_divrem_f(f, q, r, a, b, ctx);
        fmpz_mod_poly_mul(t, q, b, ctx);
        fmpz_mod_poly_add(t, t, r, ctx);

        result = (fmpz_cmp_ui(f, 1) > 0 && fmpz_cmp(f, p) < 0 && fmpz_divisible(p, f));
        if (!result)
        {
            flint_printf("FAIL (factor):\n");
            flint_printf("p = "), fmpz_print(p), flint_printf("\n\n");
            flint_printf("f = "), fmpz_print(f), flint_printf("\n\n");
            flint_printf("q1 = "), fmpz_print(q1), flint_printf("\n\n");
            flint_printf("q2 = "), fmpz_print(q2), flint_printf("\n\n");
            flint_printf("a = "), fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b = "), fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("q = "), fmpz_mod_poly_print(q, ctx), flint_printf("\n\n");
            flint_printf("r = "), fmpz_mod_poly_print(r, ctx), flint_printf("\n\n");
            flint_printf("t = "), fmpz_mod_poly_print(t, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(q, ctx);
        fmpz_mod_poly_clear(r, ctx);
        fmpz_mod_poly_clear(t, ctx);
        fmpz_clear(f);
        fmpz_clear(p);
        fmpz_clear(q1);
        fmpz_clear(q2);
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}
