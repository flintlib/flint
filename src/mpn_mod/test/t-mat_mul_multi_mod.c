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
#include "gr_mat.h"
#include "fmpz.h"

TEST_FUNCTION_START(mpn_mod_mat_mul_multi_mod, state)
{
    gr_ctx_t ctx;
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_init_mpn_mod_randtest(ctx, state);
        flint_set_num_threads(1 + n_randint(state, 4));
        gr_mat_test_mul((gr_method_mat_binary_op) mpn_mod_mat_mul_multi_mod, state, 1, 20, ctx);
        flint_set_num_threads(1);
        gr_ctx_clear(ctx);
    }

    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        gr_ctx_init_mpn_mod_randtest(ctx, state);
        flint_set_num_threads(1 + n_randint(state, 4));
        gr_mat_test_mul((gr_method_mat_binary_op) mpn_mod_mat_mul_multi_mod, state, 1, 120, ctx);
        flint_set_num_threads(1);
        gr_ctx_clear(ctx);
    }


    TEST_FUNCTION_END(state);
}
