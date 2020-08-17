/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2013 Martin Lee
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("div_newton_n_preinv....");
    fflush(stdout);

    /* Check result of divrem */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, poly_t) a, b, binv, q, r, test;

        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (binv, ctx);
        TEMPLATE(T, poly_init) (q, ctx);
        TEMPLATE(T, poly_init) (r, ctx);
        TEMPLATE(T, poly_init) (test, ctx);

        do
            TEMPLATE(T, poly_randtest) (b, state, n_randint(state, 200), ctx);
        while (b->length <= 2);
        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 200), ctx);
        if (a->length > 2 * (b->length) - 3)
            TEMPLATE(T, poly_truncate) (a, 2 * (b->length) - 3, ctx);

        TEMPLATE(T, poly_reverse) (binv, b, b->length, ctx);
        TEMPLATE(T, poly_inv_series_newton) (binv, binv, b->length, ctx);
        TEMPLATE(T, poly_div_newton_n_preinv) (q, a, b, binv, ctx);
        TEMPLATE(T, poly_divrem) (test, r, a, b, ctx);

        result = (TEMPLATE(T, poly_equal) (q, test, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            TEMPLATE(T, poly_print) (test, ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print) (q, ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print) (a, ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print) (b, ctx), flint_printf("\n\n");
            abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (binv, ctx);
        TEMPLATE(T, poly_clear) (q, ctx);
        TEMPLATE(T, poly_clear) (r, ctx);
        TEMPLATE(T, poly_clear) (test, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Check aliasing of a and q */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, poly_t) a, b, binv, q;

        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (binv, ctx);
        TEMPLATE(T, poly_init) (q, ctx);

        do
            TEMPLATE(T, poly_randtest) (b, state, n_randint(state, 200), ctx);
        while (b->length <= 2);
        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 200), ctx);
        if (a->length > 2 * (b->length) - 3)
            TEMPLATE(T, poly_truncate) (a, 2 * (b->length) - 3, ctx);

        TEMPLATE(T, poly_reverse) (binv, b, b->length, ctx);
        TEMPLATE(T, poly_inv_series_newton) (binv, binv, b->length, ctx);

        TEMPLATE(T, poly_div_newton_n_preinv) (q, a, b, binv, ctx);
        TEMPLATE(T, poly_div_newton_n_preinv) (a, a, b, binv, ctx);

        result = (TEMPLATE(T, poly_equal) (a, q, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            TEMPLATE(T, poly_print) (a, ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print) (b, ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print) (q, ctx), flint_printf("\n\n");
            abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (binv, ctx);
        TEMPLATE(T, poly_clear) (q, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Check aliasing of b and q */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, poly_t) a, b, binv, q;

        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (binv, ctx);
        TEMPLATE(T, poly_init) (q, ctx);

        do
            TEMPLATE(T, poly_randtest) (b, state, n_randint(state, 200), ctx);
        while (b->length <= 2);
        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 200), ctx);
        if (a->length > 2 * (b->length) - 3)
            TEMPLATE(T, poly_truncate) (a, 2 * (b->length) - 3, ctx);

        TEMPLATE(T, poly_reverse) (binv, b, b->length, ctx);
        TEMPLATE(T, poly_inv_series_newton) (binv, binv, b->length, ctx);

        TEMPLATE(T, poly_div_newton_n_preinv) (q, a, b, binv, ctx);
        TEMPLATE(T, poly_div_newton_n_preinv) (b, a, b, binv, ctx);

        result = (TEMPLATE(T, poly_equal) (b, q, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            TEMPLATE(T, poly_print) (a, ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print) (b, ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print) (q, ctx), flint_printf("\n\n");
            abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (binv, ctx);
        TEMPLATE(T, poly_clear) (q, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Check aliasing of binv and q */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, poly_t) a, b, binv, q;

        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (binv, ctx);
        TEMPLATE(T, poly_init) (q, ctx);


        do
            TEMPLATE(T, poly_randtest) (b, state, n_randint(state, 200), ctx);
        while (b->length <= 2);
        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 200), ctx);
        if (a->length > 2 * (b->length) - 3)
            TEMPLATE(T, poly_truncate) (a, 2 * (b->length) - 3, ctx);

        TEMPLATE(T, poly_reverse) (binv, b, b->length, ctx);
        TEMPLATE(T, poly_inv_series_newton) (binv, binv, b->length, ctx);

        TEMPLATE(T, poly_div_newton_n_preinv) (q, a, b, binv, ctx);
        TEMPLATE(T, poly_div_newton_n_preinv) (binv, a, b, binv, ctx);

        result = (TEMPLATE(T, poly_equal) (binv, q, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            TEMPLATE(T, poly_print) (a, ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print) (b, ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print) (binv, ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print) (q, ctx), flint_printf("\n\n");
            abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (binv, ctx);
        TEMPLATE(T, poly_clear) (q, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}


#endif
