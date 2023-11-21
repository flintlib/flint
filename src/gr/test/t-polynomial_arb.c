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

TEST_FUNCTION_START(gr_polynomial_arb, state)
{
    gr_ctx_t RR, RRx;
    int flags = 0;

    gr_ctx_init_real_arb(RR, 64);
    gr_ctx_init_gr_poly(RRx, RR);
    RRx->size_limit = 50;
    gr_test_ring(RRx, 1000, flags);
    gr_ctx_clear(RRx);
    gr_ctx_clear(RR);

    TEST_FUNCTION_END(state);
}
