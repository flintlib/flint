/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr.h"

TEST_FUNCTION_START(gr_acf, state)
{
    gr_ctx_t CC;
    int flags = 0;
    slong prec;

    for (prec = 64; prec <= 256; prec *= 2)
    {
        gr_ctx_init_complex_float_acf(CC, prec);
        gr_test_floating_point(CC, 100, flags);
        gr_ctx_clear(CC);
    }

    TEST_FUNCTION_END(state);
}
