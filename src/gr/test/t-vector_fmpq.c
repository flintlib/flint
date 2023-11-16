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

TEST_FUNCTION_START(gr_vector_fmpq, state)
{
    gr_ctx_t QQ, VecQQ;
    slong n;
    int flags = 0;

    gr_ctx_init_fmpz(QQ);
    QQ->size_limit = 1000;

    for (n = 0; n <= 4; n++)
    {
        gr_ctx_init_vector_space_gr_vec(VecQQ, QQ, n);
        gr_test_ring(VecQQ, 100, flags);
        gr_ctx_clear(VecQQ);
    }

    gr_ctx_clear(QQ);

    TEST_FUNCTION_END(state);
}
