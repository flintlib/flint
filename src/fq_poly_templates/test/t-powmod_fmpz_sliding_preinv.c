/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "test_helpers.h"
#include "templates.h"
#include "fmpz.h"

TEST_TEMPLATE_FUNCTION_START(T, poly_powmod_fmpz_sliding_preinv, state)
{
    int i, result;

    /* Aliasing of res and a */
    for (i = 0; i < 2.5 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, res, t, f, finv;
        ulong exp;
        fmpz_t expz;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        exp = n_randint(state, 50);
        fmpz_init_set_ui(expz, exp);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (f, ctx);
        TEMPLATE(T, poly_init) (finv, ctx);
        TEMPLATE(T, poly_init) (res, ctx);
        TEMPLATE(T, poly_init) (t, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 50), ctx);
        TEMPLATE(T, poly_randtest_not_zero) (f, state,
                                             n_randint(state, 50) + 1, ctx);

        TEMPLATE(T, poly_reverse) (finv, f, f->length, ctx);
        TEMPLATE(T, poly_inv_series_newton) (finv, finv, f->length, ctx);

        TEMPLATE(T, poly_powmod_fmpz_sliding_preinv) (res, a, expz, 0, f, finv,
                                                      ctx);
        TEMPLATE(T, poly_powmod_fmpz_sliding_preinv) (a, a, expz, 0, f, finv,
                                                      ctx);

        result = (TEMPLATE(T, poly_equal) (res, a, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a:\n");
            TEMPLATE(T, poly_print) (a, ctx), flint_printf("\n\n");
            flint_printf("f:\n");
            TEMPLATE(T, poly_print) (f, ctx), flint_printf("\n\n");
            flint_printf("res:\n");
            TEMPLATE(T, poly_print) (res, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (f, ctx);
        TEMPLATE(T, poly_clear) (finv, ctx);
        TEMPLATE(T, poly_clear) (res, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);
        fmpz_clear(expz);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Aliasing of res and f */
    for (i = 0; i < 2.5 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, res, t, f, finv;
        ulong exp;
        fmpz_t expz;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        exp = n_randint(state, 50);
        fmpz_init_set_ui(expz, exp);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (f, ctx);
        TEMPLATE(T, poly_init) (finv, ctx);
        TEMPLATE(T, poly_init) (res, ctx);
        TEMPLATE(T, poly_init) (t, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 50), ctx);
        TEMPLATE(T, poly_randtest_not_zero) (f, state,
                                             n_randint(state, 50) + 1, ctx);

        TEMPLATE(T, poly_reverse) (finv, f, f->length, ctx);
        TEMPLATE(T, poly_inv_series_newton) (finv, finv, f->length, ctx);

        TEMPLATE(T, poly_powmod_fmpz_sliding_preinv) (res, a, expz, 0, f, finv,
                                                      ctx);
        TEMPLATE(T, poly_powmod_fmpz_sliding_preinv) (f, a, expz, 0, f, finv,
                                                      ctx);

        result = (TEMPLATE(T, poly_equal) (res, f, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a:\n");
            TEMPLATE(T, poly_print) (a, ctx), flint_printf("\n\n");
            flint_printf("f:\n");
            TEMPLATE(T, poly_print) (f, ctx), flint_printf("\n\n");
            flint_printf("res:\n");
            TEMPLATE(T, poly_print) (res, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (f, ctx);
        TEMPLATE(T, poly_clear) (finv, ctx);
        TEMPLATE(T, poly_clear) (res, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);
        fmpz_clear(expz);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* No aliasing */
    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, res1, res2, t, f, finv;
        ulong exp;
        fmpz_t expz;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        exp = n_randint(state, 50);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (f, ctx);
        TEMPLATE(T, poly_init) (finv, ctx);
        TEMPLATE(T, poly_init) (res1, ctx);
        TEMPLATE(T, poly_init) (res2, ctx);
        TEMPLATE(T, poly_init) (t, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 50), ctx);
        TEMPLATE(T, poly_randtest_not_zero) (f, state,
                                             n_randint(state, 50) + 1, ctx);
        fmpz_init_set_ui(expz, exp);

        TEMPLATE(T, poly_reverse) (finv, f, f->length, ctx);
        TEMPLATE(T, poly_inv_series_newton) (finv, finv, f->length, ctx);

        TEMPLATE(T, poly_powmod_fmpz_sliding_preinv) (res1, a, expz, 0, f,
                                                      finv, ctx);
        TEMPLATE(T, poly_powmod_fmpz_sliding_preinv) (res2, a, expz, 0, f, finv, ctx);  /* TODO: Fix */

        result = (TEMPLATE(T, poly_equal) (res1, res2, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a:\n");
            TEMPLATE(T, poly_print) (a, ctx), flint_printf("\n\n");
            flint_printf("f:\n");
            TEMPLATE(T, poly_print) (f, ctx), flint_printf("\n\n");
            flint_printf("res1:\n");
            TEMPLATE(T, poly_print) (res1, ctx), flint_printf("\n\n");
            flint_printf("res2:\n");
            TEMPLATE(T, poly_print) (res2, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (f, ctx);
        TEMPLATE(T, poly_clear) (finv, ctx);
        TEMPLATE(T, poly_clear) (res1, ctx);
        TEMPLATE(T, poly_clear) (res2, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);
        fmpz_clear(expz);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Check that a^(b+c) = a^b * a^c */
    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, res1, res2, res3, res4, t, f, finv;
        fmpz_t exp1, exp2, exp3;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        fmpz_init(exp1);
        fmpz_init(exp2);
        fmpz_randtest(exp1, state, 200);
        if (fmpz_sgn(exp1) == -1)
            fmpz_neg(exp1, exp1);
        fmpz_randtest(exp2, state, 200);
        if (fmpz_sgn(exp2) == -1)
            fmpz_neg(exp2, exp2);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (f, ctx);
        TEMPLATE(T, poly_init) (finv, ctx);
        TEMPLATE(T, poly_init) (res1, ctx);
        TEMPLATE(T, poly_init) (res2, ctx);
        TEMPLATE(T, poly_init) (res3, ctx);
        TEMPLATE(T, poly_init) (res4, ctx);
        TEMPLATE(T, poly_init) (t, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 50), ctx);
        TEMPLATE(T, poly_randtest_not_zero) (f, state,
                                             n_randint(state, 50) + 1, ctx);

        TEMPLATE(T, poly_reverse) (finv, f, f->length, ctx);
        TEMPLATE(T, poly_inv_series_newton) (finv, finv, f->length, ctx);

        TEMPLATE(T, poly_powmod_fmpz_sliding_preinv) (res1, a, exp1, 0, f,
                                                      finv, ctx);
        TEMPLATE(T, poly_powmod_fmpz_sliding_preinv) (res2, a, exp2, 0, f,
                                                      finv, ctx);
        TEMPLATE(T, poly_mulmod) (res4, res1, res2, f, ctx);
        fmpz_init(exp3);
        fmpz_add(exp3, exp1, exp2);
        TEMPLATE(T, poly_powmod_fmpz_sliding_preinv) (res3, a, exp3, 0, f,
                                                      finv, ctx);

        result = (TEMPLATE(T, poly_equal) (res4, res3, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            TEMPLATE(T, ctx_print) (ctx);
            flint_printf("exp1:\n");
            fmpz_print(exp1);
            flint_printf("\n\n");
            flint_printf("exp2:\n");
            fmpz_print(exp2);
            flint_printf("\n\n");
            flint_printf("a:\n");
            TEMPLATE(T, poly_print_pretty) (a, "x", ctx), flint_printf("\n\n");
            flint_printf("f:\n");
            TEMPLATE(T, poly_print_pretty) (f, "x", ctx), flint_printf("\n\n");
            flint_printf("res1:\n");
            TEMPLATE(T, poly_print_pretty) (res1, "x", ctx),
                flint_printf("\n\n");
            flint_printf("res2:\n");
            TEMPLATE(T, poly_print_pretty) (res2, "x", ctx),
                flint_printf("\n\n");
            flint_printf("res3:\n");
            TEMPLATE(T, poly_print_pretty) (res3, "x", ctx),
                flint_printf("\n\n");
            flint_printf("res4:\n");
            TEMPLATE(T, poly_print_pretty) (res4, "x", ctx),
                flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (f, ctx);
        TEMPLATE(T, poly_clear) (finv, ctx);
        TEMPLATE(T, poly_clear) (res1, ctx);
        TEMPLATE(T, poly_clear) (res2, ctx);
        TEMPLATE(T, poly_clear) (res3, ctx);
        TEMPLATE(T, poly_clear) (res4, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);
        fmpz_clear(exp1);
        fmpz_clear(exp2);
        fmpz_clear(exp3);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
