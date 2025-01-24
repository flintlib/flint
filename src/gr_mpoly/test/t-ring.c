/*
    Copyright (C) 2020 Daniel Schultz
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_mpoly.h"

FLINT_DLL extern gr_static_method_table _ca_methods;

TEST_FUNCTION_START(gr_mpoly_ring, state)
{
    slong iter;

    for (iter = 0; iter < 30 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t cctx;
        gr_mpoly_ctx_t ctx;
        slong reps;

        gr_ctx_init_random(cctx, state);

        if (gr_ctx_is_finite(cctx) == T_TRUE ||
            gr_ctx_has_real_prec(cctx) == T_TRUE)
        {
            gr_mpoly_ctx_init_rand(ctx, state, cctx, 4);
            reps = 10;
        }
        else if (cctx->methods == _ca_methods) /* hack: slow */
        {
            gr_mpoly_ctx_init_rand(ctx, state, cctx, 1);
            reps = 1;
        }
        else
        {
            gr_mpoly_ctx_init_rand(ctx, state, cctx, 2);
            reps = 3;
        }

        /* Hack: for string conversion tests, make sure we don't have
           overlapping generator names. */
        gr_vec_t vec;
        gr_vec_init(vec, 0, cctx);
        if (gr_gens_recursive(vec, cctx) == GR_SUCCESS)
        {
            const char * vars[] = { "mv1", "mv2", "mv3", "mv4" };

            GR_MUST_SUCCEED(gr_ctx_set_gen_names(ctx, vars));

        }
        gr_vec_clear(vec, cctx);

        /* gr_ctx_println(ctx); */
        gr_test_ring(ctx, reps, 0 * GR_TEST_VERBOSE);

        gr_mpoly_ctx_clear(ctx);
        gr_ctx_clear(cctx);
    }

    TEST_FUNCTION_END(state);
}
