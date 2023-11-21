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

TEST_FUNCTION_START(gr_series_fmpq, state)
{
    gr_ctx_t QQ, QQx;
    int flags = 0;
    slong i;

    gr_ctx_init_fmpq(QQ);

    for (i = 0; i < 5; i++)
    {
        gr_ctx_init_gr_series(QQx, QQ, i);
        gr_test_ring(QQx, 100, flags);
        gr_ctx_clear(QQx);
    }

    for (i = 0; i < 5; i++)
    {
        gr_ctx_init_gr_series_mod(QQx, QQ, i);
        gr_test_ring(QQx, 100, flags);
        gr_ctx_clear(QQx);
    }

    gr_ctx_clear(QQ);

    TEST_FUNCTION_END(state);
}
