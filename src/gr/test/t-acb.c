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

TEST_FUNCTION_START(gr_acb, state)
{
    gr_ctx_t CC;
    int flags = 0;
    slong prec;

    /* the precision is clamped to [2, WORD_MAX / 8] */
    gr_ctx_init_complex_acb(CC, WORD_MAX);
    if (_gr_ctx_get_real_prec(CC) != WORD_MAX / 8)
        flint_abort();
    gr_ctx_clear(CC);
    gr_ctx_init_complex_acb(CC, -1);
    if (_gr_ctx_get_real_prec(CC) != 2)
        flint_abort();
    gr_ctx_clear(CC);

    for (prec = 64; prec <= 256; prec *= 2)
    {
        gr_ctx_init_complex_acb(CC, prec);
        gr_test_ring(CC, 1000, flags);
        gr_ctx_clear(CC);
    }

    TEST_FUNCTION_END(state);
}
