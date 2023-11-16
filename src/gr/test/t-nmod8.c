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

TEST_FUNCTION_START(gr_nmod8, state)
{
    gr_ctx_t ZZn;
    int flags = GR_TEST_ALWAYS_ABLE;
    ulong n;

    for (n = 1; n < 256; n++)
    {
        gr_ctx_init_nmod8(ZZn, n);
        gr_test_ring(ZZn, 100, flags);
        gr_ctx_clear(ZZn);
    }

    gr_ctx_init_nmod8(ZZn, 107);
    gr_test_ring(ZZn, 10000, flags);
    gr_ctx_clear(ZZn);

    gr_ctx_init_nmod8(ZZn, 10);
    gr_test_ring(ZZn, 10000, flags);
    gr_ctx_clear(ZZn);

    TEST_FUNCTION_END(state);
}
