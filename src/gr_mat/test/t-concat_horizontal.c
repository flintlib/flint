/*
    Copyright (C) 2015 Elena Sergeicheva

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_mat.h"

TEST_FUNCTION_START(gr_mat_concat_horizontal, state)
{
    slong iter;

    for (iter = 0; iter < 100; iter++)
    {
        int status = GR_SUCCESS;
    	gr_ctx_t ctx;
    	gr_mat_t A, B, C;
    	gr_mat_t window1, window2;
        slong c1, c2, r1;

        gr_ctx_init_random(ctx, state);

        c1 = n_randint(state, 5);
        c2 = n_randint(state, 5);
        r1 = n_randint(state, 5);

        gr_mat_init(A, r1, c1, ctx);
        gr_mat_init(B, r1, c2, ctx);
        gr_mat_init(C, r1, (c1 + c2), ctx);

        status |= gr_mat_randtest(A, state, ctx);
        status |= gr_mat_randtest(B, state, ctx);
        status |= gr_mat_randtest(C, state, ctx);

        status |= gr_mat_concat_horizontal(C, A, B, ctx);

        gr_mat_window_init(window1, C, 0, 0, r1, c1, ctx);
        gr_mat_window_init(window2, C, 0, c1, r1, (c1 + c2), ctx);

        if (status == GR_SUCCESS)
        {
            if (gr_mat_equal(window1, A, ctx) == T_FALSE || gr_mat_equal(window2, B, ctx) == T_FALSE)
            {
                flint_printf("FAIL: results not equal\n");
                fflush(stdout);
                flint_abort();
            }
        }

        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);
        gr_mat_clear(C, ctx);

        gr_mat_window_clear(window1, ctx);
        gr_mat_window_clear(window2, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
