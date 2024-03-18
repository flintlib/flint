/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr.h"

TEST_FUNCTION_START(gr_matrix_arb, state)
{
    gr_ctx_t RR, MatRR;
    slong n;
    int flags = 0;

    gr_ctx_init_real_arb(RR, 64);

    for (n = 0; n <= 5; n++)
    {
        gr_ctx_init_matrix_ring(MatRR, RR, n);
        gr_test_ring(MatRR, n <= 2 ? 100 : 10, flags);
        gr_ctx_clear(MatRR);
    }

    gr_ctx_clear(RR);

    TEST_FUNCTION_END(state);
}
