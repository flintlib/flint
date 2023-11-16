/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr.h"

TEST_FUNCTION_START(gr_vector_nmod, state)
{
    gr_ctx_t ZZn, VecZZn;
    slong iter, n;
    int flags = 0;

    for (iter = 0; iter < 100; iter++)
    {
        gr_ctx_init_nmod(ZZn, n_randtest_not_zero(state));

        for (n = 0; n <= 4; n++)
        {
            gr_ctx_init_vector_space_gr_vec(VecZZn, ZZn, n);
            gr_test_ring(VecZZn, 10, flags);
            gr_ctx_clear(VecZZn);
        }

        gr_ctx_clear(ZZn);
    }

    TEST_FUNCTION_END(state);
}
