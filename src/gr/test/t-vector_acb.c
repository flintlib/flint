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

TEST_FUNCTION_START(gr_vector_acb, state)
{
    gr_ctx_t CC, VecCC;
    slong n;
    int flags = 0;

    gr_ctx_init_complex_acb(CC, 64);

    for (n = 0; n <= 4; n++)
    {
        gr_ctx_init_vector_space_gr_vec(VecCC, CC, n);
        gr_test_ring(VecCC, 100, flags);
        gr_ctx_clear(VecCC);
    }

    gr_ctx_clear(CC);

    TEST_FUNCTION_END(state);
}
