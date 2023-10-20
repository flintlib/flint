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

TEST_TEMPLATE_FUNCTION_START(T, mat_charpoly, state)
{
    TEMPLATE(T, ctx_t) ctx;
    TEMPLATE(T, mat_t) A, B, C;
    TEMPLATE(T, poly_t) p1, p2;
    slong i, m, n;

    /* charpoly(AB) == charpoly(BA) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 10);
        n = m;
        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (p1, ctx);
        TEMPLATE(T, poly_init) (p2, ctx);
        TEMPLATE(T, mat_init) (A, m, n, ctx);
        TEMPLATE(T, mat_init) (B, m, n, ctx);
        TEMPLATE(T, mat_init) (C, m, n, ctx);
        TEMPLATE(T, mat_randtest) (A, state, ctx);
        TEMPLATE(T, mat_randtest) (B, state, ctx);

        TEMPLATE(T, mat_mul) (C, A, B, ctx);
        TEMPLATE(T, mat_charpoly) (p1, C, ctx);

        TEMPLATE(T, mat_mul) (C, B, A, ctx);
        TEMPLATE(T, mat_charpoly) (p2, C, ctx);

        if (!TEMPLATE(T, poly_equal) (p1, p2, ctx))
        {
            flint_printf("FAIL:\n");
            flint_printf("charpoly(AB) != charpoly(BA).\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, mat_clear) (A, ctx);
        TEMPLATE(T, mat_clear) (B, ctx);
        TEMPLATE(T, mat_clear) (C, ctx);
        TEMPLATE(T, poly_clear) (p1, ctx);
        TEMPLATE(T, poly_clear) (p2, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
