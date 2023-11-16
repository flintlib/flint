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

TEST_FUNCTION_START(gr_nmod32, state)
{
    gr_ctx_t ZZn;
    int flags = GR_TEST_ALWAYS_ABLE;
    ulong n;
    int i;

    for (i = 0; i < 100; i++)
    {
        n = (unsigned int) n_randtest(state);
        if (n == 0)
            n = 1;
        gr_ctx_init_nmod32(ZZn, n);
        gr_test_ring(ZZn, 100, flags);
        gr_ctx_clear(ZZn);
    }

    gr_ctx_init_nmod32(ZZn, UWORD(4294967291));
    gr_test_ring(ZZn, 10000, flags);
    gr_ctx_clear(ZZn);

    gr_ctx_init_nmod32(ZZn, UWORD(4294967295));
    gr_test_ring(ZZn, 10000, flags);
    gr_ctx_clear(ZZn);

    TEST_FUNCTION_END(state);
}
