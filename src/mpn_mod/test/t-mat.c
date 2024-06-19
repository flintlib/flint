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

TEST_FUNCTION_START(mpn_mod_mat, state)
{
    gr_ctx_t ctx;
    slong iter;

    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        gr_ctx_init_mpn_mod_randtest(ctx, state);
        gr_mat_test_det((gr_method_mat_unary_op_get_scalar) mpn_mod_mat_det, state, 5, 10, ctx);
        gr_mat_test_nonsingular_solve_tril((gr_method_mat_binary_op_with_flag) mpn_mod_mat_nonsingular_solve_tril, state, 5, 50, ctx);
        gr_mat_test_nonsingular_solve_triu((gr_method_mat_binary_op_with_flag) mpn_mod_mat_nonsingular_solve_triu, state, 5, 50, ctx);
        gr_mat_test_mul((gr_method_mat_binary_op) mpn_mod_mat_mul, state, 1, 100, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
