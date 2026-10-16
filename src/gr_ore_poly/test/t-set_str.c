/*
    Copyright (C) 2025 Ricardo Buring

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_ore_poly.h"

/* This is mostly a copy of src/gr_poly/test/t-set_str.c */

TEST_FUNCTION_START(gr_ore_poly_set_str, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        int status = GR_SUCCESS;
        gr_ctx_t ctx;
        gr_ore_poly_ctx_t ore_ctx;
        gr_ore_poly_t f, g;
        slong glen;
        char * s;

        gr_ore_poly_ctx_init_randtest2(ctx, ore_ctx, state);
        gr_ore_poly_init(f, ore_ctx);
        gr_ore_poly_init(g, ore_ctx);

        status |= gr_ore_poly_randtest(f, state, 1 + n_randint(state, 10), ore_ctx);
        status |= gr_ore_poly_randtest(g, state, 1 + n_randint(state, 10), ore_ctx);

        if (n_randint(state, 2))
        {
            status |= gr_ore_poly_get_str(&s, f, ore_ctx);
            status |= gr_ore_poly_set_str(g, s, ore_ctx);
            glen = -1;

            if ((ctx->which_ring == GR_CTX_FMPZ || ctx->which_ring == GR_CTX_FMPQ) &&
                status != GR_SUCCESS)
            {
                flint_printf("FAIL: unable over ZZ or QQ\n\n");
                flint_printf("f = %{gr_ore_poly}\n\n", f, ctx);
                flint_printf("s = %s\n\n", s);
                flint_abort();
            }
        }
        else
        {
            status |= _gr_ore_poly_get_str(&s, f->coeffs, f->length, ore_ctx);
            glen = n_randint(state, 10);
            gr_ore_poly_fit_length(g, glen, ore_ctx);
            status |= _gr_ore_poly_set_str(g->coeffs, s, glen, ore_ctx);
            _gr_ore_poly_set_length(g, glen, ore_ctx);
            _gr_ore_poly_normalise(g, ore_ctx);
        }

        if (status == GR_SUCCESS && gr_ore_poly_equal(f, g, ore_ctx) == T_FALSE)
        {
            flint_printf("FAIL: get_str, set_str roundtrip\n\n");
            flint_printf("f = %{gr_ore_poly}\n\n", f, ore_ctx);
            flint_printf("s = %s\n\n", s);
            flint_printf("glen = %wd\n\n", glen);
            flint_printf("g = %{gr_ore_poly}\n\n", g, ore_ctx);
            flint_abort();
        }

        flint_free(s);

        gr_ore_poly_clear(f, ore_ctx);
        gr_ore_poly_clear(g, ore_ctx);
        gr_ctx_clear(ore_ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
