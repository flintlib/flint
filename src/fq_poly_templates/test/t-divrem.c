/*
    Copyright (C) 2012 Sebastian Pancratz
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

TEST_TEMPLATE_FUNCTION_START(T, poly_divrem, state)
{
    int i, result;

    /* Check q*b + r == a */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) a, b, c, q, r;

        TEMPLATE(T, ctx_init_randtest)(ctx, state, 3);
        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (c, ctx);
        TEMPLATE(T, poly_init) (q, ctx);
        TEMPLATE(T, poly_init) (r, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 50), ctx);
        TEMPLATE(T, poly_randtest_not_zero) (b, state,
                                             n_randint(state, 50) + 1, ctx);

        TEMPLATE(T, poly_divrem) (q, r, a, b, ctx);
        TEMPLATE(T, poly_mul) (c, q, b, ctx);
        TEMPLATE(T, poly_add) (c, c, r, ctx);

        result = (TEMPLATE(T, poly_equal) (a, c, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), TEMPLATE(T, poly_print_pretty) (a, "X", ctx),
                flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, poly_print_pretty) (b, "X", ctx),
                flint_printf("\n");
            flint_printf("c = "), TEMPLATE(T, poly_print_pretty) (c, "X", ctx),
                flint_printf("\n");
            flint_printf("q = "), TEMPLATE(T, poly_print_pretty) (q, "X", ctx),
                flint_printf("\n");
            flint_printf("r = "), TEMPLATE(T, poly_print_pretty) (r, "X", ctx),
                flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (c, ctx);
        TEMPLATE(T, poly_clear) (q, ctx);
        TEMPLATE(T, poly_clear) (r, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Check aliasing: a and r */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) a, b, q, r;

        TEMPLATE(T, ctx_init_randtest)(ctx, state, 3);
        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (q, ctx);
        TEMPLATE(T, poly_init) (r, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 50), ctx);
        TEMPLATE(T, poly_randtest_not_zero) (b, state,
                                             n_randint(state, 50) + 1, ctx);

        TEMPLATE(T, poly_divrem) (q, r, a, b, ctx);
        TEMPLATE(T, poly_divrem) (q, a, a, b, ctx);

        result = (TEMPLATE(T, poly_equal) (a, r, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), TEMPLATE(T, poly_print_pretty) (a, "X", ctx),
                flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, poly_print_pretty) (b, "X", ctx),
                flint_printf("\n");
            flint_printf("q = "), TEMPLATE(T, poly_print_pretty) (q, "X", ctx),
                flint_printf("\n");
            flint_printf("r = "), TEMPLATE(T, poly_print_pretty) (r, "X", ctx),
                flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (q, ctx);
        TEMPLATE(T, poly_clear) (r, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Check aliasing: b and r */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) a, b, q, r;

        TEMPLATE(T, ctx_init_randtest)(ctx, state, 3);
        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (q, ctx);
        TEMPLATE(T, poly_init) (r, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 50), ctx);
        TEMPLATE(T, poly_randtest_not_zero) (b, state,
                                             n_randint(state, 50) + 1, ctx);

        TEMPLATE(T, poly_divrem) (q, r, a, b, ctx);
        TEMPLATE(T, poly_divrem) (q, b, a, b, ctx);

        result = (TEMPLATE(T, poly_equal) (b, r, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), TEMPLATE(T, poly_print_pretty) (a, "X", ctx),
                flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, poly_print_pretty) (b, "X", ctx),
                flint_printf("\n");
            flint_printf("q = "), TEMPLATE(T, poly_print_pretty) (q, "X", ctx),
                flint_printf("\n");
            flint_printf("r = "), TEMPLATE(T, poly_print_pretty) (r, "X", ctx),
                flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (q, ctx);
        TEMPLATE(T, poly_clear) (r, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Check aliasing: a and q */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) a, b, q, r;

        TEMPLATE(T, ctx_init_randtest)(ctx, state, 3);
        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (q, ctx);
        TEMPLATE(T, poly_init) (r, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 50), ctx);
        TEMPLATE(T, poly_randtest_not_zero) (b, state,
                                             n_randint(state, 50) + 1, ctx);

        TEMPLATE(T, poly_divrem) (q, r, a, b, ctx);
        TEMPLATE(T, poly_divrem) (a, r, a, b, ctx);

        result = (TEMPLATE(T, poly_equal) (a, q, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), TEMPLATE(T, poly_print_pretty) (a, "X", ctx),
                flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, poly_print_pretty) (b, "X", ctx),
                flint_printf("\n");
            flint_printf("q = "), TEMPLATE(T, poly_print_pretty) (q, "X", ctx),
                flint_printf("\n");
            flint_printf("r = "), TEMPLATE(T, poly_print_pretty) (r, "X", ctx),
                flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (q, ctx);
        TEMPLATE(T, poly_clear) (r, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Check aliasing: b and q */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) a, b, q, r;

        TEMPLATE(T, ctx_init_randtest)(ctx, state, 3);
        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (q, ctx);
        TEMPLATE(T, poly_init) (r, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 50), ctx);
        TEMPLATE(T, poly_randtest_not_zero) (b, state,
                                             n_randint(state, 50) + 1, ctx);

        TEMPLATE(T, poly_divrem) (q, r, a, b, ctx);
        TEMPLATE(T, poly_divrem) (b, r, a, b, ctx);

        result = (TEMPLATE(T, poly_equal) (b, q, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), TEMPLATE(T, poly_print_pretty) (a, "X", ctx),
                flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, poly_print_pretty) (b, "X", ctx),
                flint_printf("\n");
            flint_printf("q = "), TEMPLATE(T, poly_print_pretty) (q, "X", ctx),
                flint_printf("\n");
            flint_printf("r = "), TEMPLATE(T, poly_print_pretty) (r, "X", ctx),
                flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (q, ctx);
        TEMPLATE(T, poly_clear) (r, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Check _poly_rem against divrem, including lenA < lenB */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, b, q, r;
        TEMPLATE(T, struct) * R;
        TEMPLATE(T, t) invB;
        slong lenA, lenB, rlen;

        TEMPLATE(T, ctx_init_randtest)(ctx, state, 3);
        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (q, ctx);
        TEMPLATE(T, poly_init) (r, ctx);
        TEMPLATE(T, init) (invB, ctx);

        TEMPLATE(T, poly_randtest_not_zero) (b, state,
                                             n_randint(state, 20) + 1, ctx);
        lenB = b->length;
        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, lenB + 10), ctx);
        lenA = a->length;

        TEMPLATE(T, inv) (invB, b->coeffs + lenB - 1, ctx);

        rlen = FLINT_MAX(lenB - 1, 1);
        R = _TEMPLATE(T, vec_init) (rlen, ctx);

        /* garbage in the output */
        TEMPLATE(T, randtest_not_zero) (R + 0, state, ctx);
        if (lenB >= 2)
            TEMPLATE(T, randtest_not_zero) (R + lenB - 2, state, ctx);

        _TEMPLATE(T, poly_rem) (R, a->coeffs, lenA, b->coeffs, lenB, invB, ctx);

        TEMPLATE(T, poly_divrem) (q, r, a, b, ctx);

        result = (r->length <= lenB - 1)
            && _TEMPLATE(T, vec_equal) (R, r->coeffs, r->length, ctx)
            && _TEMPLATE(T, vec_is_zero) (R + r->length, lenB - 1 - r->length, ctx);

        if (!result)
        {
            flint_printf("FAIL (_poly_rem):\n\n");
            flint_printf("lenA = %wd, lenB = %wd\n", lenA, lenB);
            flint_printf("a = "), TEMPLATE(T, poly_print_pretty) (a, "X", ctx),
                flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, poly_print_pretty) (b, "X", ctx),
                flint_printf("\n");
            flint_printf("r = "), TEMPLATE(T, poly_print_pretty) (r, "X", ctx),
                flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        _TEMPLATE(T, vec_clear) (R, rlen, ctx);
        TEMPLATE(T, clear) (invB, ctx);
        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (q, ctx);
        TEMPLATE(T, poly_clear) (r, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
