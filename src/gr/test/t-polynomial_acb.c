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

TEST_FUNCTION_START(gr_polynomial_acb, state)
{
    gr_ctx_t CC, CCx;
    int flags = 0;

    gr_ctx_init_complex_acb(CC, 64);
    gr_ctx_init_gr_poly(CCx, CC);
    CCx->size_limit = 50;
    gr_test_ring(CCx, 1000, flags);
    gr_ctx_clear(CCx);
    gr_ctx_clear(CC);

    TEST_FUNCTION_END(state);
}
