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

TEST_FUNCTION_START(gr_polynomial_fmpz, state)
{
    gr_ctx_t ZZ, ZZx;
    int flags = 0;

    gr_ctx_init_fmpz(ZZ);
    ZZ->size_limit = 200;
    gr_ctx_init_gr_poly(ZZx, ZZ);
    ZZx->size_limit = 30;
    gr_test_ring(ZZx, 1000, flags);
    gr_ctx_clear(ZZx);
    gr_ctx_clear(ZZ);

    TEST_FUNCTION_END(state);
}
