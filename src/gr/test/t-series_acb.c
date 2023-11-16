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

TEST_FUNCTION_START(gr_series_acb, state)
{
    gr_ctx_t CCn, CCnx;
    int flags = 0;
    slong i;

    for (i = 0; i < 5; i++)
    {
        gr_ctx_init_complex_acb(CCn, 2 + n_randint(state, 200));
        gr_ctx_init_gr_series(CCnx, CCn, i);
        gr_test_ring(CCnx, 100, flags);
        gr_ctx_clear(CCnx);
        gr_ctx_clear(CCn);
    }

    for (i = 0; i < 5; i++)
    {
        gr_ctx_init_complex_acb(CCn, 2 + n_randint(state, 200));
        gr_ctx_init_gr_series_mod(CCnx, CCn, i);
        gr_test_ring(CCnx, 100, flags);
        gr_ctx_clear(CCnx);
        gr_ctx_clear(CCn);
    }

    gr_ctx_clear(CCn);

    TEST_FUNCTION_END(state);
}
