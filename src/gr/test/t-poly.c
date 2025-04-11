/*
    Copyright (C) 2023 Fredrik Johansson
    Copyright (C) 2025 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr.h"

int gr_test_factor(gr_ctx_t R, flint_rand_t state, int test_flags);

TEST_FUNCTION_START(gr_poly, state)
{
    gr_ctx_t cctx, pctx;
    int flags = 0;

    for (slong iter = 0; iter < 100; iter++)
    {
        gr_ctx_init_random(cctx, state);
        cctx->size_limit = 30;
        gr_ctx_init_gr_poly(pctx, cctx);
        pctx->size_limit = 30;
        gr_test_ring(pctx, 10, flags);
        gr_ctx_clear(pctx);
        gr_ctx_clear(cctx);
    }

    /* selected tests for rings currently not covered by gr_ctx_init_random */

    {
        gr_ctx_init_complex_qqbar(cctx);
        gr_ctx_init_gr_poly(pctx, cctx);
        cctx->size_limit = 30;
        pctx->size_limit = 30;
        gr_test_factor(pctx, state, flags);
        gr_ctx_clear(pctx);
        gr_ctx_clear(cctx);
    }

    TEST_FUNCTION_END(state);
}
