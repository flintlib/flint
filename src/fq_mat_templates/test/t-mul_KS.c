/*
    Copyright (C) 2011 Fredrik Johansson
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

TEST_TEMPLATE_FUNCTION_START(T, mat_mul_KS, state)
{
    slong i;

    /* Check aliasing C and A */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, mat_t) A, B, C;
        slong m, n;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        m = n_randint(state, 50);
        n = n_randint(state, 50);

        TEMPLATE(T, mat_init) (A, m, n, ctx);
        TEMPLATE(T, mat_init) (B, n, n, ctx);
        TEMPLATE(T, mat_init) (C, m, n, ctx);

        TEMPLATE(T, mat_randtest) (A, state, ctx);
        TEMPLATE(T, mat_randtest) (B, state, ctx);
        TEMPLATE(T, mat_randtest) (C, state, ctx);  /* noise in output */

        TEMPLATE(T, mat_mul_KS) (C, A, B, ctx);
        TEMPLATE(T, mat_mul_KS) (A, A, B, ctx);

        if (!TEMPLATE(T, mat_equal) (C, A, ctx))
        {
            printf("FAIL:\n");
            printf("A:\n");
            TEMPLATE(T, mat_print) (A, ctx);
            printf("B:\n");
            TEMPLATE(T, mat_print) (B, ctx);
            printf("C:\n");
            TEMPLATE(T, mat_print) (C, ctx);
            printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, mat_clear) (A, ctx);
        TEMPLATE(T, mat_clear) (B, ctx);
        TEMPLATE(T, mat_clear) (C, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Check aliasing C and B */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, mat_t) A, B, C;
        slong m, n;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        m = n_randint(state, 50);
        n = n_randint(state, 50);

        TEMPLATE(T, mat_init) (A, m, m, ctx);
        TEMPLATE(T, mat_init) (B, m, n, ctx);
        TEMPLATE(T, mat_init) (C, m, n, ctx);

        TEMPLATE(T, mat_randtest) (A, state, ctx);
        TEMPLATE(T, mat_randtest) (B, state, ctx);
        TEMPLATE(T, mat_randtest) (C, state, ctx);  /* noise in output */

        TEMPLATE(T, mat_mul_KS) (C, A, B, ctx);
        TEMPLATE(T, mat_mul_KS) (B, A, B, ctx);

        if (!TEMPLATE(T, mat_equal) (C, B, ctx))
        {
            printf("FAIL:\n");
            printf("A:\n");
            TEMPLATE(T, mat_print) (A, ctx);
            printf("B:\n");
            TEMPLATE(T, mat_print) (B, ctx);
            printf("C:\n");
            TEMPLATE(T, mat_print) (C, ctx);
            printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, mat_clear) (A, ctx);
        TEMPLATE(T, mat_clear) (B, ctx);
        TEMPLATE(T, mat_clear) (C, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
