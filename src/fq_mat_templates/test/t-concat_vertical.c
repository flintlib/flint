/*
    Copyright (C) 2015 Elena Sergeicheva

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "test_helpers.h"
#include "templates.h"

TEST_TEMPLATE_FUNCTION_START(T, mat_concat_vertical, state)
{
    slong i;

    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
    	TEMPLATE(T, ctx_t) ctx;
    	TEMPLATE(T, mat_t) A, B, C;
    	TEMPLATE(T, mat_t) window1, window2;
        slong r1, r2, c1;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        r1 = n_randint(state, 10);
        r2 = n_randint(state, 10);
        c1 = n_randint(state, 10);

        TEMPLATE(T, mat_init) (A, r1, c1, ctx);
        TEMPLATE(T, mat_init) (B, r2, c1, ctx);
        TEMPLATE(T, mat_init) (C, (r1 + r2), c1, ctx);

        TEMPLATE(T, mat_randtest) (A, state, ctx);
        TEMPLATE(T, mat_randtest) (B, state, ctx);

        TEMPLATE(T, mat_randtest) (C, state, ctx);

        TEMPLATE(T, mat_concat_vertical) (C, A, B, ctx);

        TEMPLATE(T, mat_window_init) (window1, C, 0, 0, r1, c1, ctx);
        TEMPLATE(T, mat_window_init) (window2, C, r1, 0, (r1+r2), c1, ctx);

        if (!(TEMPLATE(T, mat_equal) (window1, A, ctx) && TEMPLATE(T, mat_equal) (window2, B, ctx)))
        {
            flint_printf("FAIL: results not equal\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, mat_clear) (A, ctx);
        TEMPLATE(T, mat_clear) (B, ctx);
        TEMPLATE(T, mat_clear) (C, ctx);

        TEMPLATE(T, mat_window_clear) (window1, ctx);
        TEMPLATE(T, mat_window_clear) (window2, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
