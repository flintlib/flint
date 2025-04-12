/*
    Copyright (C) 2025 Ricardo Buring

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_vec.h"
#include "gr_ore_poly.h"

/* This is mostly a copy of src/gr_mpoly/test/t-ring.c */

FLINT_DLL extern gr_static_method_table _ca_methods;

TEST_FUNCTION_START(gr_ore_poly_ring, state)
{
    slong iter;

    for (iter = 0; iter < 30 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        gr_ore_poly_ctx_t ore_ctx;
        slong reps;

        gr_ctx_init_random(ctx, state);
        /* TODO: ensure ctx has a generator */
        gr_ore_poly_ctx_init_rand(ore_ctx, state, ctx);

        if (gr_ctx_is_finite(ctx) == T_TRUE ||
            gr_ctx_has_real_prec(ctx) == T_TRUE)
        {
            reps = 10;
        }
        else if (ctx->methods == _ca_methods) /* hack: slow */
        {
            reps = 1;
        }
        else
        {
            reps = 3;
        }

        /* Hack: for string conversion tests, make sure we don't have
           overlapping generator names. */
        gr_vec_t vec;
        gr_vec_init(vec, 0, ctx);
        if (gr_gens_recursive(vec, ctx) == GR_SUCCESS)
        {
            const char * vars[] = { "DD" };

            GR_MUST_SUCCEED(gr_ctx_set_gen_names(ore_ctx, vars));

        }
        gr_vec_clear(vec, ctx);

        /* gr_ctx_println(ore_ctx); */
        gr_test_ring(ore_ctx, reps, 0 * GR_TEST_VERBOSE);

        gr_ore_poly_ctx_clear(ore_ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
