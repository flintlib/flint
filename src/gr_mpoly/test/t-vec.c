/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_mpoly.h"

TEST_FUNCTION_START(gr_mpoly_vec, state)
{
    slong iter;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t cctx;
        gr_mpoly_ctx_t ctx;
        gr_mpoly_vec_t vec, vec2;
        gr_mpoly_t f;
        slong i, len;
        int status = GR_SUCCESS;

        gr_ctx_init_random(cctx, state);
        gr_mpoly_ctx_init_rand(ctx, state, cctx, 3);

        len = n_randint(state, 6);

        gr_mpoly_vec_init(vec, len, ctx);
        gr_mpoly_vec_init(vec2, 0, ctx);
        gr_mpoly_init(f, ctx);

        if (gr_mpoly_vec_length(vec, ctx) != len)
        {
            flint_printf("FAIL: length after init\n");
            fflush(stdout);
            flint_abort();
        }

        /* fill entries and append copies into vec2 */
        for (i = 0; i < len; i++)
        {
            status |= gr_mpoly_randtest_bits(gr_mpoly_vec_entry_ptr(vec, i, ctx),
                                             state, 1 + n_randint(state, 6),
                                             2 + n_randint(state, 6), ctx);
            status |= gr_mpoly_vec_append(vec2, gr_mpoly_vec_entry_srcptr(vec, i, ctx), ctx);
        }

        if (status == GR_SUCCESS && gr_mpoly_vec_length(vec2, ctx) != len)
        {
            flint_printf("FAIL: length after append\n");
            fflush(stdout);
            flint_abort();
        }

        /* set then compare entrywise */
        {
            gr_mpoly_vec_t vec3;
            gr_mpoly_vec_init(vec3, 0, ctx);
            status |= gr_mpoly_vec_set(vec3, vec, ctx);

            if (status == GR_SUCCESS)
            {
                for (i = 0; i < len; i++)
                {
                    if (gr_mpoly_equal(gr_mpoly_vec_entry_srcptr(vec3, i, ctx),
                                       gr_mpoly_vec_entry_srcptr(vec, i, ctx), ctx) == T_FALSE)
                    {
                        flint_printf("FAIL: set mismatch at %wd\n", i);
                        fflush(stdout);
                        flint_abort();
                    }
                }
            }

            gr_mpoly_vec_clear(vec3, ctx);
        }

        /* set_length grows with zeros and shrinks */
        gr_mpoly_vec_set_length(vec2, len + 3, ctx);
        for (i = len; i < len + 3; i++)
        {
            if (gr_mpoly_is_zero(gr_mpoly_vec_entry_ptr(vec2, i, ctx), ctx) == T_FALSE)
            {
                flint_printf("FAIL: grown entries not zero\n");
                fflush(stdout);
                flint_abort();
            }
        }
        gr_mpoly_vec_set_length(vec2, 1, ctx);
        if (gr_mpoly_vec_length(vec2, ctx) != 1)
        {
            flint_printf("FAIL: shrink\n");
            fflush(stdout);
            flint_abort();
        }

        gr_mpoly_clear(f, ctx);
        gr_mpoly_vec_clear(vec, ctx);
        gr_mpoly_vec_clear(vec2, ctx);
        gr_mpoly_ctx_clear(ctx);
        gr_ctx_clear(cctx);
    }

    TEST_FUNCTION_END(state);
}
