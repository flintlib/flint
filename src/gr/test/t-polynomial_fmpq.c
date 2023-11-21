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

TEST_FUNCTION_START(gr_polynomial_fmpq, state)
{
    gr_ctx_t QQ, QQx;
    int flags = 0;

    gr_ctx_init_fmpq(QQ);
    QQ->size_limit = 100;
    gr_ctx_init_gr_poly(QQx, QQ);
    QQx->size_limit = 30;
    gr_test_ring(QQx, 1000, flags);
    gr_ctx_clear(QQx);
    gr_ctx_clear(QQ);

    TEST_FUNCTION_END(state);
}
