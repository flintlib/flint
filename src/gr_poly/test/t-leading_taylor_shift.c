/*
    Copyright (C) 2025 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_poly.h"

TEST_FUNCTION_START(gr_poly_leading_taylor_shift, state)
{
    for (slong i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        gr_ctx_t ctx;
        gr_poly_t f, g;
        gr_ptr s;

        gr_ctx_init_random(ctx, state);

        s = gr_heap_init(ctx);
        gr_poly_init(f, ctx);
        gr_poly_init(g, ctx);

        int status = GR_SUCCESS;

        status |= gr_poly_randtest(f, state, n_randint(state, 10), ctx);
        status |= gr_poly_randtest(g, state, n_randint(state, 10), ctx);
        slong len = f->length;
        status |= gr_poly_truncate(g, g, len, ctx);

        if (len >= 1 && g->length == len && n_randint(state, 8))
            status |= gr_set(GR_ENTRY(g->coeffs, len - 1, ctx->sizeof_elem),
                             GR_ENTRY(f->coeffs, len - 1, ctx->sizeof_elem),
                             ctx);

        status |= gr_poly_leading_taylor_shift(s, f, g, ctx);

        if (status == GR_UNABLE && ctx->which_ring == GR_CTX_FMPZ)
        {
            flint_printf("FAIL (unexpected failure)\n");
            fflush(stdout);
            flint_abort();
        }

        status |= gr_poly_taylor_shift(f, f, s, ctx);
        if (len >= 2)
        {
            status |= gr_poly_shift_right(f, f, len - 2, ctx);
            status |= gr_poly_shift_right(g, g, len - 2, ctx);
        }

        if (status == GR_SUCCESS && gr_poly_equal(f, g, ctx) == T_FALSE)
        {
            flint_printf("FAIL (incorrect result)\n");
            fflush(stdout);
            flint_abort();
        }

        gr_poly_clear(f, ctx);
        gr_poly_clear(g, ctx);
        gr_heap_clear(s, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}

