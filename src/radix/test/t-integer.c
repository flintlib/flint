/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "radix.h"
#include "gr.h"

TEST_FUNCTION_START(radix_integer, state)
{
    slong iter;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        radix_t radix;
        gr_ctx_t ctx;

        radix_init_randtest(radix, state);
        gr_ctx_init_radix_integer(ctx, DIGIT_RADIX(radix), radix->exp);

        gr_test_ring(ctx, 10, 0);

        gr_ctx_clear(ctx);
        radix_clear(radix);
    }

    TEST_FUNCTION_END(state);
}
