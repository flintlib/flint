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

TEST_FUNCTION_START(gr_fmpz_poly, state)
{
    gr_ctx_t ZZx;
    int flags = GR_TEST_ALWAYS_ABLE;

    gr_ctx_init_fmpz_poly(ZZx);
    gr_test_ring(ZZx, 1000, flags);

    gr_ctx_clear(ZZx);

    TEST_FUNCTION_END(state);
}
