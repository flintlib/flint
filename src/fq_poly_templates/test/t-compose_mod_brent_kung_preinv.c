/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Martin Lee
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "test_helpers.h"
#include "templates.h"

TEST_TEMPLATE_FUNCTION_START(T, poly_compose_mod_brent_kung_preinv, state)
{
    int i;

    for (i = 0; i < 30 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, b, c, cinv, d, e;

        TEMPLATE(T, ctx_init_randtest)(ctx, state, 3);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (c, ctx);
        TEMPLATE(T, poly_init) (cinv, ctx);
        TEMPLATE(T, poly_init) (d, ctx);
        TEMPLATE(T, poly_init) (e, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 20) + 1, ctx);
        TEMPLATE(T, poly_randtest) (b, state, n_randint(state, 20) + 1, ctx);
        TEMPLATE(T, poly_randtest_not_zero) (c, state,
                                             n_randint(state, 20) + 1, ctx);

        TEMPLATE(T, poly_reverse) (cinv, c, c->length, ctx);
        TEMPLATE(T, poly_inv_series_newton) (cinv, cinv, c->length, ctx);

        TEMPLATE(T, poly_rem) (a, a, c, ctx);
        TEMPLATE(T, poly_compose_mod_brent_kung_preinv) (d, a, b, c, cinv,
                                                         ctx);
        TEMPLATE(T, poly_compose) (e, a, b, ctx);
        TEMPLATE(T, poly_rem) (e, e, c, ctx);

        if (!TEMPLATE(T, poly_equal) (d, e, ctx))
        {
            flint_printf("FAIL (composition):\n");
            flint_printf("a:\n");
            TEMPLATE(T, poly_print) (a, ctx);
            flint_printf("\n");
            flint_printf("b:\n");
            TEMPLATE(T, poly_print) (b, ctx);
            flint_printf("\n");
            flint_printf("c:\n");
            TEMPLATE(T, poly_print) (c, ctx);
            flint_printf("\n");
            flint_printf("d:\n");
            TEMPLATE(T, poly_print) (d, ctx);
            flint_printf("\n");
            flint_printf("e:\n");
            TEMPLATE(T, poly_print) (e, ctx);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (c, ctx);
        TEMPLATE(T, poly_clear) (cinv, ctx);
        TEMPLATE(T, poly_clear) (d, ctx);
        TEMPLATE(T, poly_clear) (e, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Test aliasing of res and a */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, b, c, cinv, d;

        TEMPLATE(T, ctx_init_randtest)(ctx, state, 3);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (c, ctx);
        TEMPLATE(T, poly_init) (cinv, ctx);
        TEMPLATE(T, poly_init) (d, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 20) + 1, ctx);
        TEMPLATE(T, poly_randtest) (b, state, n_randint(state, 20) + 1, ctx);
        TEMPLATE(T, poly_randtest_not_zero) (c, state,
                                             n_randint(state, 20) + 1, ctx);

        TEMPLATE(T, poly_reverse) (cinv, c, c->length, ctx);
        TEMPLATE(T, poly_inv_series_newton) (cinv, cinv, c->length, ctx);

        TEMPLATE(T, poly_rem) (a, a, c, ctx);
        TEMPLATE(T, poly_compose_mod_brent_kung_preinv) (d, a, b, c, cinv,
                                                         ctx);
        TEMPLATE(T, poly_compose_mod_brent_kung_preinv) (a, a, b, c, cinv,
                                                         ctx);

        if (!TEMPLATE(T, poly_equal) (d, a, ctx))
        {
            flint_printf("FAIL (aliasing a):\n");
            flint_printf("a:\n");
            TEMPLATE(T, poly_print) (a, ctx);
            flint_printf("\n");
            flint_printf("b:\n");
            TEMPLATE(T, poly_print) (b, ctx);
            flint_printf("\n");
            flint_printf("c:\n");
            TEMPLATE(T, poly_print) (c, ctx);
            flint_printf("\n");
            flint_printf("d:\n");
            TEMPLATE(T, poly_print) (d, ctx);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (c, ctx);
        TEMPLATE(T, poly_clear) (cinv, ctx);
        TEMPLATE(T, poly_clear) (d, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Test aliasing of res and b */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, b, c, cinv, d;

        TEMPLATE(T, ctx_init_randtest)(ctx, state, 3);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (c, ctx);
        TEMPLATE(T, poly_init) (cinv, ctx);
        TEMPLATE(T, poly_init) (d, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 20) + 1, ctx);
        TEMPLATE(T, poly_randtest) (b, state, n_randint(state, 20) + 1, ctx);
        TEMPLATE(T, poly_randtest_not_zero) (c, state,
                                             n_randint(state, 20) + 1, ctx);

        TEMPLATE(T, poly_reverse) (cinv, c, c->length, ctx);
        TEMPLATE(T, poly_inv_series_newton) (cinv, cinv, c->length, ctx);

        TEMPLATE(T, poly_rem) (a, a, c, ctx);
        TEMPLATE(T, poly_compose_mod_brent_kung_preinv) (d, a, b, c, cinv,
                                                         ctx);
        TEMPLATE(T, poly_compose_mod_brent_kung_preinv) (b, a, b, c, cinv,
                                                         ctx);

        if (!TEMPLATE(T, poly_equal) (d, b, ctx))
        {
            flint_printf("FAIL (aliasing b)\n");
            flint_printf("a:\n");
            TEMPLATE(T, poly_print) (a, ctx);
            flint_printf("\n");
            flint_printf("b:\n");
            TEMPLATE(T, poly_print) (b, ctx);
            flint_printf("\n");
            flint_printf("c:\n");
            TEMPLATE(T, poly_print) (c, ctx);
            flint_printf("\n");
            flint_printf("d:\n");
            TEMPLATE(T, poly_print) (d, ctx);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (c, ctx);
        TEMPLATE(T, poly_clear) (cinv, ctx);
        TEMPLATE(T, poly_clear) (d, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Test aliasing of res and c */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, b, c, cinv, d;

        TEMPLATE(T, ctx_init_randtest)(ctx, state, 3);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (c, ctx);
        TEMPLATE(T, poly_init) (cinv, ctx);
        TEMPLATE(T, poly_init) (d, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 20) + 1, ctx);
        TEMPLATE(T, poly_randtest) (b, state, n_randint(state, 20) + 1, ctx);
        TEMPLATE(T, poly_randtest_not_zero) (c, state,
                                             n_randint(state, 20) + 1, ctx);

        TEMPLATE(T, poly_reverse) (cinv, c, c->length, ctx);
        TEMPLATE(T, poly_inv_series_newton) (cinv, cinv, c->length, ctx);

        TEMPLATE(T, poly_rem) (a, a, c, ctx);
        TEMPLATE(T, poly_compose_mod_brent_kung_preinv) (d, a, b, c, cinv,
                                                         ctx);
        TEMPLATE(T, poly_compose_mod_brent_kung_preinv) (c, a, b, c, cinv,
                                                         ctx);

        if (!TEMPLATE(T, poly_equal) (d, c, ctx))
        {
            flint_printf("FAIL (aliasing c)\n");
            flint_printf("a:\n");
            TEMPLATE(T, poly_print) (a, ctx);
            flint_printf("\n");
            flint_printf("b:\n");
            TEMPLATE(T, poly_print) (b, ctx);
            flint_printf("\n");
            flint_printf("c:\n");
            TEMPLATE(T, poly_print) (c, ctx);
            flint_printf("\n");
            flint_printf("d:\n");
            TEMPLATE(T, poly_print) (d, ctx);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (c, ctx);
        TEMPLATE(T, poly_clear) (cinv, ctx);
        TEMPLATE(T, poly_clear) (d, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
