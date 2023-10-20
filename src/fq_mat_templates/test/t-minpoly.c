/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "test_helpers.h"
#include "templates.h"

TEST_TEMPLATE_FUNCTION_START(T, mat_minpoly, state)
{
    TEMPLATE(T, ctx_t) ctx;
    TEMPLATE(T, t) t;
    TEMPLATE(T, mat_t) A, B;
    TEMPLATE(T, poly_t) p1, p2, q, r;
    slong i, j, k, m, n;

    /* minpoly(A) divides charpoly(A) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 10);
        n = m;
        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (p1, ctx);
        TEMPLATE(T, poly_init) (p2, ctx);
        TEMPLATE(T, poly_init) (q, ctx);
        TEMPLATE(T, poly_init) (r, ctx);
        TEMPLATE(T, mat_init) (A, m, n, ctx);
        TEMPLATE(T, mat_randtest) (A, state, ctx);

        TEMPLATE(T, mat_charpoly) (p1, A, ctx);
        TEMPLATE(T, mat_minpoly) (p2, A, ctx);

        TEMPLATE(T, poly_divrem) (q, r, p1, p2, ctx);

        if (!TEMPLATE(T, poly_is_zero) (r, ctx))
        {
            flint_printf("FAIL:\n");
            flint_printf("minpoly(A) doesn't divide charpoly(A).\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, mat_clear) (A, ctx);
        TEMPLATE(T, poly_clear) (p1, ctx);
        TEMPLATE(T, poly_clear) (p2, ctx);
        TEMPLATE(T, poly_clear) (q, ctx);
        TEMPLATE(T, poly_clear) (r, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* minpoly(P^{-1}AP) == minpoly(A) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 10);
        n = m;
        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, init) (t, ctx);
        TEMPLATE(T, poly_init) (p1, ctx);
        TEMPLATE(T, poly_init) (p2, ctx);
        TEMPLATE(T, mat_init) (A, m, n, ctx);
        TEMPLATE(T, mat_init) (B, m, n, ctx);
        TEMPLATE(T, mat_randtest) (A, state, ctx);

        for (j = 0; j < n/2; j++)
        {
           for (k = 0; k < n/2; k++)
           {
              TEMPLATE(T, zero) (TEMPLATE(T, mat_entry) (A, j + n/2, k), ctx);
              TEMPLATE(T, zero) (TEMPLATE(T, mat_entry) (A, j, k + n/2), ctx);
              TEMPLATE(T, set) (TEMPLATE(T, mat_entry) (A, j + n/2, k + n/2),
                                TEMPLATE(T, mat_entry) (A, j, k), ctx);
           }
        }

        TEMPLATE(T, mat_set) (B, A, ctx);

        TEMPLATE(T, mat_minpoly) (p1, A, ctx);

        for (j = 0; j < n; j++)
        {
           TEMPLATE(T, set_ui) (t, n_randint(state, 6) - 3, ctx);
           TEMPLATE(T, mat_similarity) (B, n_randint(state, n), t, ctx);
        }

        TEMPLATE(T, mat_minpoly) (p2, B, ctx);

        if (!TEMPLATE(T, poly_equal) (p1, p2, ctx))
        {
            flint_printf("FAIL:\n");
            flint_printf("minpoly(P^{-1}AP) != minpoly(A).\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, clear) (t, ctx);
        TEMPLATE(T, mat_clear) (A, ctx);
        TEMPLATE(T, mat_clear) (B, ctx);
        TEMPLATE(T, poly_clear) (p1, ctx);
        TEMPLATE(T, poly_clear) (p2, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
