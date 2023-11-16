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

TEST_FUNCTION_START(gr_acb, state)
{
    gr_ctx_t CC;
    int flags = 0;
    slong prec;

    for (prec = 64; prec <= 256; prec *= 2)
    {
        gr_ctx_init_complex_acb(CC, prec);
        gr_test_ring(CC, 1000, flags);
        gr_ctx_clear(CC);
    }

    TEST_FUNCTION_END(state);
}
