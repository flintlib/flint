/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr.h"

TEST_FUNCTION_START(gr_polynomial_nmod8, state)
{
    gr_ctx_t ZZn, ZZnx;
    slong iter;
    int flags = 0;

    for (iter = 0; iter < 50; iter++)
    {
        gr_ctx_init_nmod8(ZZn, 1 + n_randtest(state) % 255);

        gr_ctx_init_gr_poly(ZZnx, ZZn);
        ZZnx->size_limit = 30;
        gr_test_ring(ZZnx, 1000, flags);
        gr_ctx_clear(ZZnx);

        gr_ctx_clear(ZZn);
    }

    TEST_FUNCTION_END(state);
}
