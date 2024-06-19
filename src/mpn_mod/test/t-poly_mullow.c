/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_mod.h"
#include "gr_poly.h"

TEST_FUNCTION_START(mpn_mod_poly_mullow, state)
{
    gr_ctx_t ctx;
    slong iter;

    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        gr_ctx_init_mpn_mod_randtest(ctx, state);
        _gr_poly_test_mullow((gr_method_poly_binary_trunc_op) _mpn_mod_poly_mullow, NULL, state, 10, 100, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
