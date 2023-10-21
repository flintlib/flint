/*
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

TEST_TEMPLATE_FUNCTION_START(T, poly_compose_mod_horner, state)
{
    int i;

    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, b, c, d, e;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (c, ctx);
        TEMPLATE(T, poly_init) (d, ctx);
        TEMPLATE(T, poly_init) (e, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 20) + 1, ctx);
        TEMPLATE(T, poly_randtest) (b, state, n_randint(state, 20) + 1, ctx);
        TEMPLATE(T, poly_randtest_not_zero) (c, state,
                                             n_randint(state, 20) + 1, ctx);

        TEMPLATE(T, poly_compose_mod_horner) (d, a, b, c, ctx);
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
        TEMPLATE(T, poly_clear) (d, ctx);
        TEMPLATE(T, poly_clear) (e, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Test aliasing of res and a */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, b, c, d;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (c, ctx);
        TEMPLATE(T, poly_init) (d, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 20) + 1, ctx);
        TEMPLATE(T, poly_randtest) (b, state, n_randint(state, 20) + 1, ctx);
        TEMPLATE(T, poly_randtest_not_zero) (c, state,
                                             n_randint(state, 20) + 1, ctx);

        TEMPLATE(T, poly_compose_mod_horner) (d, a, b, c, ctx);
        TEMPLATE(T, poly_compose_mod_horner) (a, a, b, c, ctx);

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
        TEMPLATE(T, poly_clear) (d, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Test aliasing of res and b */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, b, c, d;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (c, ctx);
        TEMPLATE(T, poly_init) (d, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 20) + 1, ctx);
        TEMPLATE(T, poly_randtest) (b, state, n_randint(state, 20) + 1, ctx);
        TEMPLATE(T, poly_randtest_not_zero) (c, state,
                                             n_randint(state, 20) + 1, ctx);

        TEMPLATE(T, poly_compose_mod_horner) (d, a, b, c, ctx);
        TEMPLATE(T, poly_compose_mod_horner) (b, a, b, c, ctx);

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
        TEMPLATE(T, poly_clear) (d, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Test aliasing of res and c */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, b, c, d;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (c, ctx);
        TEMPLATE(T, poly_init) (d, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 20) + 1, ctx);
        TEMPLATE(T, poly_randtest) (b, state, n_randint(state, 20) + 1, ctx);
        TEMPLATE(T, poly_randtest_not_zero) (c, state,
                                             n_randint(state, 20) + 1, ctx);

        TEMPLATE(T, poly_compose_mod_horner) (d, a, b, c, ctx);
        TEMPLATE(T, poly_compose_mod_horner) (c, a, b, c, ctx);

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
        TEMPLATE(T, poly_clear) (d, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
