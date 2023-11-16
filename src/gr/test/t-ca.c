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

TEST_FUNCTION_START(gr_ca, state)
{
    gr_ctx_t RR, CC, QQbar_real, QQbar;
    int flags = 0;

    gr_ctx_init_real_ca(RR);
    gr_test_ring(RR, 100, flags);
    gr_ctx_clear(RR);

    gr_ctx_init_complex_ca(CC);
    gr_test_ring(CC, 100, flags);
    gr_ctx_clear(CC);

    gr_ctx_init_real_algebraic_ca(QQbar_real);
    gr_test_ring(QQbar_real, 100, flags);
    gr_ctx_clear(QQbar_real);

    gr_ctx_init_complex_algebraic_ca(QQbar);
    gr_test_ring(QQbar, 100, flags);
    gr_ctx_clear(QQbar);

    TEST_FUNCTION_END(state);
}
