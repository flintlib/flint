/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr.h"

TEST_FUNCTION_START(gr_nmod, state)
{
    gr_ctx_t ZZn;
    int flags = GR_TEST_ALWAYS_ABLE;
    slong iter;
    ulong n;

    for (iter = 0; iter < 100; iter++)
    {
        n = n_randtest_not_zero(state);
        gr_ctx_init_nmod(ZZn, n);
        if (n_randint(state, 2))
            GR_MUST_SUCCEED(gr_ctx_set_is_field(ZZn, n_is_prime(n) ? T_TRUE : T_FALSE));
        gr_test_ring(ZZn, 100, flags);
        gr_ctx_clear(ZZn);
    }

    TEST_FUNCTION_END(state);
}
