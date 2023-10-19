/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

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

TEST_FUNCTION_START(fmpz_mod_poly_powmod_ui_binexp, state)
{
    int i, result;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctx, 2);

    /* Aliasing of res and a */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, res1, t, f;
        fmpz_t p;
        ulong exp;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);

        exp = n_randint(state, 50);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(f, ctx);
        fmpz_mod_poly_init(res1, ctx);
        fmpz_mod_poly_init(t, ctx);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 50), ctx);
        fmpz_mod_poly_randtest_not_zero(f, state, n_randint(state, 50) + 1, ctx);

        fmpz_mod_poly_powmod_ui_binexp(res1, a, exp, f, ctx);
        fmpz_mod_poly_powmod_ui_binexp(a, a, exp, f, ctx);

        result = (fmpz_mod_poly_equal(res1, a, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("exp: %wu\n\n", exp);
            flint_printf("a:\n"); fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("f:\n"); fmpz_mod_poly_print(f, ctx), flint_printf("\n\n");
            flint_printf("res:\n"); fmpz_mod_poly_print(res1, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(p);
        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(f, ctx);
        fmpz_mod_poly_clear(res1, ctx);
        fmpz_mod_poly_clear(t, ctx);
    }

    /* Aliasing of res and f */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, res1, t, f;
        fmpz_t p;
        ulong exp;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);

        exp = n_randint(state, 50);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(f, ctx);
        fmpz_mod_poly_init(res1, ctx);
        fmpz_mod_poly_init(t, ctx);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 50), ctx);
        fmpz_mod_poly_randtest_not_zero(f, state, n_randint(state, 50) + 1, ctx);

        fmpz_mod_poly_powmod_ui_binexp(res1, a, exp, f, ctx);
        fmpz_mod_poly_powmod_ui_binexp(f, a, exp, f, ctx);

        result = (fmpz_mod_poly_equal(res1, f, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("exp: %wu\n\n", exp);
            flint_printf("a:\n"); fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("f:\n"); fmpz_mod_poly_print(f, ctx), flint_printf("\n\n");
            flint_printf("res1:\n"); fmpz_mod_poly_print(res1, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(p);
        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(f, ctx);
        fmpz_mod_poly_clear(res1, ctx);
        fmpz_mod_poly_clear(t, ctx);
    }

    /* No aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, res1, res2, t, f;
        fmpz_t p;
        ulong exp;
        int j;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);

        exp = n_randint(state, 50);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(f, ctx);
        fmpz_mod_poly_init(res1, ctx);
        fmpz_mod_poly_init(res2, ctx);
        fmpz_mod_poly_init(t, ctx);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 50), ctx);
        fmpz_mod_poly_randtest_not_zero(f, state, n_randint(state, 50) + 1, ctx);

        fmpz_mod_poly_powmod_ui_binexp(res1, a, exp, f, ctx);

        fmpz_mod_poly_zero(res2, ctx);
        if (fmpz_mod_poly_length(f, ctx) > 1)
            fmpz_mod_poly_set_coeff_ui(res2, 0, 1, ctx);
        for (j = 1; j <= exp; j++)
            fmpz_mod_poly_mulmod(res2, res2, a, f, ctx);

        result = (fmpz_mod_poly_equal(res1, res2, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("exp: %wu\n\n", exp);
            flint_printf("a:\n"); fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("f:\n"); fmpz_mod_poly_print(f, ctx), flint_printf("\n\n");
            flint_printf("res1:\n"); fmpz_mod_poly_print(res1, ctx), flint_printf("\n\n");
            flint_printf("res2:\n"); fmpz_mod_poly_print(res2, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(p);
        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(f, ctx);
        fmpz_mod_poly_clear(res1, ctx);
        fmpz_mod_poly_clear(res2, ctx);
        fmpz_mod_poly_clear(t, ctx);
    }

    /* Check that a^(b+c) = a^b * a^c */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, res1, res2, res3, res4, t, f;
        fmpz_t p;
        ulong exp1, exp2, exp3;

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);

        exp1 = n_randint(state, 50);
        exp2 = n_randint(state, 50);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(f, ctx);
        fmpz_mod_poly_init(res1, ctx);
        fmpz_mod_poly_init(res2, ctx);
        fmpz_mod_poly_init(res3, ctx);
        fmpz_mod_poly_init(res4, ctx);
        fmpz_mod_poly_init(t, ctx);

        fmpz_mod_poly_randtest(a, state, n_randint(state, 50), ctx);
        fmpz_mod_poly_randtest_not_zero(f, state, n_randint(state, 50) + 1, ctx);

        fmpz_mod_poly_powmod_ui_binexp(res1, a, exp1, f, ctx);
        fmpz_mod_poly_powmod_ui_binexp(res2, a, exp2, f, ctx);
        fmpz_mod_poly_mulmod(res4, res1, res2, f, ctx);
        exp3 = exp1 + exp2;
        fmpz_mod_poly_powmod_ui_binexp(res3, a, exp3, f, ctx);

        result = (fmpz_mod_poly_equal(res4, res3, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a:\n"); fmpz_mod_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("f:\n"); fmpz_mod_poly_print(f, ctx), flint_printf("\n\n");
            flint_printf("res3:\n"); fmpz_mod_poly_print(res3, ctx), flint_printf("\n\n");
            flint_printf("res4:\n"); fmpz_mod_poly_print(res4, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(p);
        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(f, ctx);
        fmpz_mod_poly_clear(res1, ctx);
        fmpz_mod_poly_clear(res2, ctx);
        fmpz_mod_poly_clear(res3, ctx);
        fmpz_mod_poly_clear(res4, ctx);
        fmpz_mod_poly_clear(t, ctx);
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}
