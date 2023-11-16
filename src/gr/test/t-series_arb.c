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

TEST_FUNCTION_START(gr_series_arb, state)
{
    gr_ctx_t RRn, RRnx;
    int flags = 0;
    slong i;

    for (i = 0; i < 5; i++)
    {
        gr_ctx_init_real_arb(RRn, 2 + n_randint(state, 200));
        gr_ctx_init_gr_series(RRnx, RRn, i);
        gr_test_ring(RRnx, 100, flags);
        gr_ctx_clear(RRnx);
        gr_ctx_clear(RRn);
    }

    for (i = 0; i < 5; i++)
    {
        gr_ctx_init_real_arb(RRn, 2 + n_randint(state, 200));
        gr_ctx_init_gr_series_mod(RRnx, RRn, i);
        gr_test_ring(RRnx, 100, flags);
        gr_ctx_clear(RRnx);
        gr_ctx_clear(RRn);
    }

    gr_ctx_clear(RRn);

    TEST_FUNCTION_END(state);
}
