/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr.h"

TEST_FUNCTION_START(gr_dirichlet, state)
{
    gr_ctx_t G;
    int flags = GR_TEST_ALWAYS_ABLE;
    ulong q;

    for (q = 1; q <= 100; q++)
    {
        GR_MUST_SUCCEED(gr_ctx_init_dirichlet_group(G, q));
        gr_test_multiplicative_group(G, 100, flags);
        gr_ctx_clear(G);
    }

    {
        slong iter;

        for (iter = 0; iter < 100; iter++)
        {
            q = n_randtest(state);
            if (gr_ctx_init_dirichlet_group(G, q) == GR_SUCCESS)
            {
                gr_test_multiplicative_group(G, 100, flags);
                gr_ctx_clear(G);
            }
        }
    }

    TEST_FUNCTION_END(state);
}
