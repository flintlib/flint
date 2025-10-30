/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr.h"

TEST_FUNCTION_START(gr_debug, state)
{
    gr_ctx_t R, Rdebug;

    gr_ctx_init_fmpz(R);

    gr_ctx_init_debug(Rdebug, R, 0, 0.0);
    gr_test_ring(Rdebug, 1000, 0);
    gr_ctx_clear(Rdebug);

    gr_ctx_init_debug(Rdebug, R, 0, 0.2);
    gr_test_ring(Rdebug, 1000, 0);
    gr_ctx_clear(Rdebug);

    gr_ctx_init_debug(Rdebug, R, GR_DEBUG_WRAP, 0.2);
    gr_test_ring(Rdebug, 1000, 0);
    gr_ctx_clear(Rdebug);

    gr_ctx_clear(R);

    TEST_FUNCTION_END(state);
}
