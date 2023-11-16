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

TEST_FUNCTION_START(gr_qqbar, state)
{
    gr_ctx_t QQbar_real, QQbar;
    int flags = GR_TEST_ALWAYS_ABLE;

    gr_ctx_init_real_qqbar(QQbar_real);
    gr_test_ring(QQbar_real, 100, flags);
    gr_ctx_clear(QQbar_real);

    gr_ctx_init_complex_qqbar(QQbar);
    gr_test_ring(QQbar, 100, flags);
    gr_ctx_clear(QQbar);

    TEST_FUNCTION_END(state);
}
