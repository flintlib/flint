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

TEST_FUNCTION_START(gr_vector_arb, state)
{
    gr_ctx_t RR, VecRR;
    slong n;
    int flags = 0;

    gr_ctx_init_real_arb(RR, 64);

    for (n = 0; n <= 4; n++)
    {
        gr_ctx_init_vector_space_gr_vec(VecRR, RR, n);
        gr_test_ring(VecRR, 100, flags);
        gr_ctx_clear(VecRR);
    }

    gr_ctx_clear(RR);

    TEST_FUNCTION_END(state);
}
