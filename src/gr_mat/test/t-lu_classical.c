/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_mat.h"

TEST_FUNCTION_START(gr_mat_lu_classical, state)
{
    gr_ctx_t ctx;

    gr_ctx_init_fmpz(ctx);
    gr_mat_test_lu((gr_method_mat_lu_op) gr_mat_lu_classical, state, 100 * flint_test_multiplier(), 10, ctx);
    gr_ctx_clear(ctx);

    gr_ctx_init_fmpq(ctx);
    gr_mat_test_lu((gr_method_mat_lu_op) gr_mat_lu_classical, state, 100 * flint_test_multiplier(), 10, ctx);
    gr_ctx_clear(ctx);

    gr_ctx_init_nmod(ctx, 2);
    gr_mat_test_lu((gr_method_mat_lu_op) gr_mat_lu_classical, state, 100 * flint_test_multiplier(), 10, ctx);
    gr_ctx_clear(ctx);

    gr_mat_test_lu((gr_method_mat_lu_op) gr_mat_lu_classical, state, 100 * flint_test_multiplier(), 10, NULL);

    TEST_FUNCTION_END(state);
}
