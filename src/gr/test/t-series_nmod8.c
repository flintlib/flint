/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr.h"

TEST_FUNCTION_START(gr_series_nmod8, state)
{
    gr_ctx_t ZZn, ZZnx;
    int flags = 0;
    slong i;

    for (i = 0; i < 5; i++)
    {
        gr_ctx_init_nmod8(ZZn, 1 + n_randtest(state) % 255);
        gr_ctx_init_gr_series(ZZnx, ZZn, i);
        gr_test_ring(ZZnx, 100, flags);
        gr_ctx_clear(ZZnx);
        gr_ctx_clear(ZZn);
    }

    for (i = 0; i < 5; i++)
    {
        gr_ctx_init_nmod8(ZZn, 1 + n_randtest(state) % 255);
        gr_ctx_init_gr_series_mod(ZZnx, ZZn, i);
        gr_test_ring(ZZnx, 100, flags);
        gr_ctx_clear(ZZnx);
        gr_ctx_clear(ZZn);
    }

    TEST_FUNCTION_END(state);
}
