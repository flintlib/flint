/*
    Copyright (C) 2015 Elena Sergeicheva

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_mat.h"

TEST_FUNCTION_START(gr_mat_concat_vertical, state)
{
    slong iter;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        int status = GR_SUCCESS;
        gr_ctx_t ctx;
        gr_mat_t A, B, C, BA, BB, BC;
        gr_mat_t window1, window2;
        gr_mat_struct *a, *b, *c;
        slong r1, r2, c1, ea, eb, ec;

        gr_ctx_init_random(ctx, state);

        r1 = n_randint(state, 5);
        r2 = n_randint(state, 5);
        c1 = n_randint(state, 5);

        /* optionally use windows of larger matrices, so that the
           rows are not stored contiguously */
        ea = n_randint(state, 2) ? n_randint(state, 3) : 0;
        eb = n_randint(state, 2) ? n_randint(state, 3) : 0;
        ec = n_randint(state, 2) ? n_randint(state, 3) : 0;

        gr_mat_init(BA, r1, c1 + ea, ctx);
        gr_mat_init(BB, r2, c1 + eb, ctx);
        gr_mat_init(BC, r1 + r2, c1 + ec, ctx);
        gr_mat_window_init(A, BA, 0, 0, r1, c1, ctx);
        gr_mat_window_init(B, BB, 0, 0, r2, c1, ctx);
        gr_mat_window_init(C, BC, 0, 0, r1 + r2, c1, ctx);
        a = (ea == 0) ? BA : A;
        b = (eb == 0) ? BB : B;
        c = (ec == 0) ? BC : C;

        status |= gr_mat_randtest(BA, state, ctx);
        status |= gr_mat_randtest(BB, state, ctx);
        status |= gr_mat_randtest(BC, state, ctx);

        status |= gr_mat_concat_vertical(c, a, b, ctx);

        gr_mat_window_init(window1, c, 0, 0, r1, c1, ctx);
        gr_mat_window_init(window2, c, r1, 0, r1 + r2, c1, ctx);

        if (status == GR_SUCCESS)
        {
            if (gr_mat_equal(window1, a, ctx) == T_FALSE || gr_mat_equal(window2, b, ctx) == T_FALSE)
            {
                flint_printf("FAIL: results not equal\n");
                fflush(stdout);
                flint_abort();
            }
        }

        gr_mat_window_clear(window1, ctx);
        gr_mat_window_clear(window2, ctx);
        gr_mat_window_clear(A, ctx);
        gr_mat_window_clear(B, ctx);
        gr_mat_window_clear(C, ctx);
        gr_mat_clear(BA, ctx);
        gr_mat_clear(BB, ctx);
        gr_mat_clear(BC, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
