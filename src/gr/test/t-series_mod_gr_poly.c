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

TEST_FUNCTION_START(gr_series_mod_gr_poly, state)
{
    gr_ctx_t R, Rx;
    int flags = 0;
    slong i;
    slong iter;

    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        for (i = 0; i < 4; i++)
        {
            gr_ctx_init_random(R, state);
            gr_ctx_init_series_mod_gr_poly(Rx, R, i);
            /* hack: avoid collisions with parent ring generators */
            GR_MUST_SUCCEED(gr_ctx_set_gen_name(Rx, "u"));
            gr_test_ring(Rx, 10, flags);
            gr_ctx_clear(Rx);
            gr_ctx_clear(R);
        }
    }

    TEST_FUNCTION_END(state);
}
