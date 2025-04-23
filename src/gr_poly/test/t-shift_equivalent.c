/*
    Copyright (C) 2025 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "gr_poly.h"

TEST_FUNCTION_START(gr_poly_shift_equivalent, state)
{
    for (slong i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        gr_ctx_t ctx;
        gr_poly_t f, g;
        gr_ptr a;
        fmpz_t s0, s1;
        truth_t equiv = T_UNKNOWN;

        gr_ctx_init_random(ctx, state);

        gr_poly_init(f, ctx);
        gr_poly_init(g, ctx);
        a = gr_heap_init(ctx);
        fmpz_init(s0);
        fmpz_init(s1);

        int status = GR_SUCCESS;

        status |= gr_poly_randtest(f, state, n_randint(state, 10), ctx);

        if (n_randint(state, 2))
        {
            fmpz_randtest(s0, state, 16);
            status |= gr_set_fmpz(a, s0, ctx);
            status |= gr_poly_taylor_shift(g, f, a, ctx);
            equiv = T_TRUE;
        }
        else
        {
            status |= gr_poly_randtest(g, state, n_randint(state, 10), ctx);
        }

        if (status != GR_SUCCESS)
            goto cleanup;

        fmpz * s = n_randint(state, 4) ? s1 : NULL;

        int res = gr_poly_shift_equivalent(s, f, g, ctx);

        if (res == T_FALSE && equiv == T_TRUE)
        {
            flint_printf("FAIL (false negative)\n");
            fflush(stdout);
            flint_abort();
        }

        if (res == T_UNKNOWN && ctx->which_ring == GR_CTX_FMPZ)
        {
            flint_printf("FAIL (unexpected failure)\n");
            fflush(stdout);
            flint_abort();
        }

        if (res == T_TRUE && s != NULL)
        {
            status |= gr_set_fmpz(a, s, ctx);
            status |= gr_poly_taylor_shift(f, f, a, ctx);

            if (status == GR_SUCCESS && gr_poly_equal(f, g, ctx) == T_FALSE)
            {
                flint_printf("FAIL (false positive or incorrect shift)\n");
                fflush(stdout);
                flint_abort();
            }
        }

cleanup:
        fmpz_clear(s0);
        fmpz_clear(s1);
        gr_heap_clear(a, ctx);
        gr_poly_clear(f, ctx);
        gr_poly_clear(g, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
