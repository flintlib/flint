/*
    Copyright (C) 2012 Sebastian Pancratz
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

    flint_printf("xgcd_euclidean....");
    fflush(stdout);

    /* Generic case, most likely co-prime arguments ***************************** */

    /* Compare with result from GCD and check correctness */
    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, b, d, g, s, t, v, w;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (d, ctx);
        TEMPLATE(T, poly_init) (g, ctx);
        TEMPLATE(T, poly_init) (s, ctx);
        TEMPLATE(T, poly_init) (t, ctx);
        TEMPLATE(T, poly_init) (v, ctx);
        TEMPLATE(T, poly_init) (w, ctx);
        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 100), ctx);
        TEMPLATE(T, poly_randtest) (b, state, n_randint(state, 100), ctx);

        TEMPLATE(T, poly_gcd_euclidean) (d, a, b, ctx);
        TEMPLATE(T, poly_xgcd_euclidean) (g, s, t, a, b, ctx);

        TEMPLATE(T, poly_mul) (v, s, a, ctx);
        TEMPLATE(T, poly_mul) (w, t, b, ctx);
        TEMPLATE(T, poly_add) (w, v, w, ctx);

        result = (TEMPLATE(T, poly_equal) (d, g, ctx) &&
                  TEMPLATE(T, poly_equal) (g, w, ctx));

        if (!result)
        {
            flint_printf("FAIL:\n");
            TEMPLATE(T, poly_print_pretty) (a, "x", ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (b, "x", ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (d, "x", ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (g, "x", ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (s, "x", ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (t, "x", ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (v, "x", ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (w, "x", ctx), flint_printf("\n\n");
            abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (d, ctx);
        TEMPLATE(T, poly_clear) (g, ctx);
        TEMPLATE(T, poly_clear) (s, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);
        TEMPLATE(T, poly_clear) (v, ctx);
        TEMPLATE(T, poly_clear) (w, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Special case, arguments share a factor ******************************* */

    /* Compare with result from GCD and check correctness */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, b, d, f, g, s, t, v, w;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (d, ctx);
        TEMPLATE(T, poly_init) (f, ctx);
        TEMPLATE(T, poly_init) (g, ctx);
        TEMPLATE(T, poly_init) (s, ctx);
        TEMPLATE(T, poly_init) (t, ctx);
        TEMPLATE(T, poly_init) (v, ctx);
        TEMPLATE(T, poly_init) (w, ctx);
        TEMPLATE(T, poly_randtest) (s, state, n_randint(state, 10), ctx);
        TEMPLATE(T, poly_randtest) (t, state, n_randint(state, 10), ctx);
        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 10), ctx);
        TEMPLATE(T, poly_randtest) (b, state, n_randint(state, 10), ctx);
        TEMPLATE(T, poly_randtest) (f, state, n_randint(state, 5), ctx);
        TEMPLATE(T, poly_mul) (a, a, f, ctx);
        TEMPLATE(T, poly_mul) (b, b, f, ctx);

        TEMPLATE(T, poly_gcd_euclidean) (d, a, b, ctx);
        TEMPLATE(T, poly_xgcd_euclidean) (g, s, t, a, b, ctx);

        TEMPLATE(T, poly_mul) (v, s, a, ctx);
        TEMPLATE(T, poly_mul) (w, t, b, ctx);
        TEMPLATE(T, poly_add) (w, v, w, ctx);

        result = (TEMPLATE(T, poly_equal) (d, g, ctx) &&
                  TEMPLATE(T, poly_equal) (g, w, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            TEMPLATE(T, poly_print_pretty) (a, "x", ctx);
            flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (b, "x", ctx);
            flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (d, "x", ctx);
            flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (f, "x", ctx);
            flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (g, "x", ctx);
            flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (s, "x", ctx);
            flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (t, "x", ctx);
            flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (v, "x", ctx);
            flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (w, "x", ctx);
            flint_printf("\n\n");
            abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (d, ctx);
        TEMPLATE(T, poly_clear) (f, ctx);
        TEMPLATE(T, poly_clear) (g, ctx);
        TEMPLATE(T, poly_clear) (s, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);
        TEMPLATE(T, poly_clear) (v, ctx);
        TEMPLATE(T, poly_clear) (w, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}

#endif
