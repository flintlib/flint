/*
    Copyright (C) 2012 Sebastian Pancratz

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

TEST_FUNCTION_START(fmpz_mod_poly_invmod, state)
{
    int i, result;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctx, 2);

    /* Test aliasing *************************************************************/

    /* Aliasing c and a */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, c;
        int ans1, ans2;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(c, ctx);

        do
            fmpz_mod_poly_randtest(b, state, n_randint(state, 100), ctx);
        while (b->length < 3);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);

        ans1 = fmpz_mod_poly_invmod(c, a, b, ctx);
        ans2 = fmpz_mod_poly_invmod(a, a, b, ctx);

        result = (ans1 == ans2 && fmpz_mod_poly_equal(a, c, ctx));
        if (!result)
        {
            flint_printf("FAIL (alias a and c):\n");
            fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(c, ctx), flint_printf("\n\n");
            flint_printf("ans1 = %d\n\n", ans1);
            flint_printf("ans2 = %d\n\n", ans2);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(c, ctx);
        fmpz_clear(p);
    }

    /* Aliasing c and b */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, c;
        int ans1, ans2;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(c, ctx);

        do
            fmpz_mod_poly_randtest(b, state, n_randint(state, 100), ctx);
        while (b->length < 3);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);

        ans1 = fmpz_mod_poly_invmod(c, a, b, ctx);
        ans2 = fmpz_mod_poly_invmod(b, a, b, ctx);

        result = ((ans1 == ans2) && fmpz_mod_poly_equal(b, c, ctx));
        if (!result)
        {
            flint_printf("FAIL (alias b and c):\n");
            fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(c, ctx), flint_printf("\n\n");
            flint_printf("ans1 = %d\n\n", ans1);
            flint_printf("ans2 = %d\n\n", ans2);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(c, ctx);
        fmpz_clear(p);
    }

    /* Compare with result from XGCD */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, g, s, t, u;
        int ans;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(g, ctx);
        fmpz_mod_poly_init(s, ctx);
        fmpz_mod_poly_init(t, ctx);
        fmpz_mod_poly_init(u, ctx);

        do
            fmpz_mod_poly_randtest(b, state, n_randint(state, 100), ctx);
        while (b->length < 3);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);

        ans = fmpz_mod_poly_invmod(u, a, b, ctx);
        fmpz_mod_poly_xgcd(g, s, t, a, b, ctx);

        result = (((ans) && g->length == 1
                        && fmpz_is_one(g->coeffs) && fmpz_mod_poly_equal(s, u, ctx))
                 || (!(ans) && g->length > 1));

        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(g, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(s, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(t, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(u, ctx), flint_printf("\n\n");
            flint_printf("ans = %d\n\n", ans);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(g, ctx);
        fmpz_mod_poly_clear(s, ctx);
        fmpz_mod_poly_clear(t, ctx);
        fmpz_mod_poly_clear(u, ctx);
        fmpz_clear(p);
    }

    /* Special case, arguments share a factor ********************************/

    /* Check correctness */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t a, b, f, u;
        int ans;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(f, ctx);
        fmpz_mod_poly_init(u, ctx);

        do
            fmpz_mod_poly_randtest(b, state, n_randint(state, 100), ctx);
        while (b->length < 2);
        fmpz_mod_poly_randtest(a, state, n_randint(state, 100), ctx);
        do
            fmpz_mod_poly_randtest_not_zero(f, state, n_randint(state, 20) + 1, ctx);
        while (f->length < 2);
        fmpz_mod_poly_mul(a, f, a, ctx);
        fmpz_mod_poly_mul(b, f, b, ctx);

        ans = fmpz_mod_poly_invmod(u, a, b, ctx);

        result = (!ans);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(b, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(f, ctx), flint_printf("\n\n");
            fmpz_mod_poly_print(u, ctx), flint_printf("\n\n");
            flint_printf("ans = %d\n\n", ans);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(f, ctx);
        fmpz_mod_poly_clear(u, ctx);
        fmpz_clear(p);
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}
