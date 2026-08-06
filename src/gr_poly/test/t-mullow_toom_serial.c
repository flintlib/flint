/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr_poly.h"

slong toom_num_points = 0;
slong toom_splitting_degree = 0;

static int
_gr_poly_mullow_toom_serial_wrapper(gr_ptr res, gr_srcptr A, slong len1, gr_srcptr B, slong len2, slong trunc, gr_ctx_t ctx)
{
    return _gr_poly_mullow_toom_serial(res, A, len1, B, len2, trunc, toom_num_points, toom_splitting_degree, ctx);
}

TEST_FUNCTION_START(gr_poly_mullow_toom_serial, state)
{
    slong iter;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        slong n;

        if (n_randint(state, 2))
        {
            gr_ctx_init_fmpq(ctx);
            n = 10;
        }
        else
        {
            gr_ctx_init_nmod(ctx, n_randtest_not_zero(state));
            n = 50;
        }

        if (n_randint(state, 2))
        {
            toom_num_points = 1 + n_randint(state, 10);
            toom_splitting_degree = 0;
        }
        else
        {
            toom_num_points = 0;
            toom_splitting_degree = 1 + n_randint(state, 10);
        }

        _gr_poly_test_mullow((gr_method_poly_binary_trunc_op) _gr_poly_mullow_toom_serial_wrapper, NULL, state, 10, n, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}

