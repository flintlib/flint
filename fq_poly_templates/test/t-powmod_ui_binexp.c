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

#include "templates.h"

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("powmod_ui_binexp....");
    fflush(stdout);

    /* Aliasing of res and a */
    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, res1, t, f;
        ulong exp;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        exp = n_randint(state, 50);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (f, ctx);
        TEMPLATE(T, poly_init) (res1, ctx);
        TEMPLATE(T, poly_init) (t, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 20), ctx);
        TEMPLATE(T, poly_randtest_not_zero) (f, state,
                                             n_randint(state, 20) + 1, ctx);

        TEMPLATE(T, poly_powmod_ui_binexp) (res1, a, exp, f, ctx);
        TEMPLATE(T, poly_powmod_ui_binexp) (a, a, exp, f, ctx);

        result = (TEMPLATE(T, poly_equal) (res1, a, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("exp: %wu\n\n", exp);
            flint_printf("a:\n");
            TEMPLATE(T, poly_print) (a, ctx), flint_printf("\n\n");
            flint_printf("f:\n");
            TEMPLATE(T, poly_print) (f, ctx), flint_printf("\n\n");
            flint_printf("res:\n");
            TEMPLATE(T, poly_print) (res1, ctx), flint_printf("\n\n");
            abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (f, ctx);
        TEMPLATE(T, poly_clear) (res1, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Aliasing of res and f */
    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, res1, t, f;
        ulong exp;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        exp = n_randint(state, 50);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (f, ctx);
        TEMPLATE(T, poly_init) (res1, ctx);
        TEMPLATE(T, poly_init) (t, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 20), ctx);
        TEMPLATE(T, poly_randtest_not_zero) (f, state,
                                             n_randint(state, 20) + 1, ctx);

        TEMPLATE(T, poly_powmod_ui_binexp) (res1, a, exp, f, ctx);
        TEMPLATE(T, poly_powmod_ui_binexp) (f, a, exp, f, ctx);

        result = (TEMPLATE(T, poly_equal) (res1, f, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("exp: %wu\n\n", exp);
            flint_printf("a:\n");
            TEMPLATE(T, poly_print) (a, ctx), flint_printf("\n\n");
            flint_printf("f:\n");
            TEMPLATE(T, poly_print) (f, ctx), flint_printf("\n\n");
            flint_printf("res1:\n");
            TEMPLATE(T, poly_print) (res1, ctx), flint_printf("\n\n");
            abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (f, ctx);
        TEMPLATE(T, poly_clear) (res1, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* No aliasing */
    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, res1, res2, t, f;
        ulong exp;
        int j;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        exp = n_randint(state, 50);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (f, ctx);
        TEMPLATE(T, poly_init) (res1, ctx);
        TEMPLATE(T, poly_init) (res2, ctx);
        TEMPLATE(T, poly_init) (t, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 20), ctx);
        TEMPLATE(T, poly_randtest_not_zero) (f, state,
                                             n_randint(state, 20) + 1, ctx);

        TEMPLATE(T, poly_powmod_ui_binexp) (res1, a, exp, f, ctx);

        TEMPLATE(T, poly_zero) (res2, ctx);
        TEMPLATE(T, poly_one) (res2, ctx);

        for (j = 1; j <= exp; j++)
            TEMPLATE(T, poly_mulmod) (res2, res2, a, f, ctx);

        result = (TEMPLATE(T, poly_equal) (res1, res2, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("exp: %wu\n\n", exp);
            flint_printf("a:\n");
            TEMPLATE(T, poly_print) (a, ctx), flint_printf("\n\n");
            flint_printf("f:\n");
            TEMPLATE(T, poly_print) (f, ctx), flint_printf("\n\n");
            flint_printf("res1:\n");
            TEMPLATE(T, poly_print) (res1, ctx), flint_printf("\n\n");
            flint_printf("res2:\n");
            TEMPLATE(T, poly_print) (res2, ctx), flint_printf("\n\n");
            abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (f, ctx);
        TEMPLATE(T, poly_clear) (res1, ctx);
        TEMPLATE(T, poly_clear) (res2, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Check that a^(b+c) = a^b * a^c */
    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, res1, res2, res3, res4, t, f;

        ulong exp1, exp2, exp3;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        exp1 = n_randint(state, 50);
        exp2 = n_randint(state, 50);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (f, ctx);
        TEMPLATE(T, poly_init) (res1, ctx);
        TEMPLATE(T, poly_init) (res2, ctx);
        TEMPLATE(T, poly_init) (res3, ctx);
        TEMPLATE(T, poly_init) (res4, ctx);
        TEMPLATE(T, poly_init) (t, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 20), ctx);
        TEMPLATE(T, poly_randtest_not_zero) (f, state,
                                             n_randint(state, 20) + 1, ctx);

        TEMPLATE(T, poly_powmod_ui_binexp) (res1, a, exp1, f, ctx);
        TEMPLATE(T, poly_powmod_ui_binexp) (res2, a, exp2, f, ctx);
        TEMPLATE(T, poly_mulmod) (res4, res1, res2, f, ctx);
        exp3 = exp1 + exp2;
        TEMPLATE(T, poly_powmod_ui_binexp) (res3, a, exp3, f, ctx);

        result = (TEMPLATE(T, poly_equal) (res4, res3, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a:\n");
            TEMPLATE(T, poly_print) (a, ctx), flint_printf("\n\n");
            flint_printf("f:\n");
            TEMPLATE(T, poly_print) (f, ctx), flint_printf("\n\n");
            flint_printf("res3:\n");
            TEMPLATE(T, poly_print) (res3, ctx), flint_printf("\n\n");
            flint_printf("res4:\n");
            TEMPLATE(T, poly_print) (res4, ctx), flint_printf("\n\n");
            abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (f, ctx);
        TEMPLATE(T, poly_clear) (res1, ctx);
        TEMPLATE(T, poly_clear) (res2, ctx);
        TEMPLATE(T, poly_clear) (res3, ctx);
        TEMPLATE(T, poly_clear) (res4, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}


#endif
