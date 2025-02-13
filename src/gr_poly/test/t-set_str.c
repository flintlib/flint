/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "gr_vec.h"
#include "gr_poly.h"

TEST_FUNCTION_START(gr_poly_set_str, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        int status = GR_SUCCESS;
        gr_ctx_t ctx;
        gr_poly_t f, g;
        slong glen;
        char * s;

        gr_ctx_init_random(ctx, state);
        gr_poly_init(f, ctx);
        gr_poly_init(g, ctx);

        status |= gr_poly_randtest(f, state, 1 + n_randint(state, 10), ctx);
        status |= gr_poly_randtest(g, state, 1 + n_randint(state, 10), ctx);

        if (n_randint(state, 2))
        {
            status |= gr_poly_get_str(&s, f, "XX", ctx);
            status |= gr_poly_set_str(g, s, "XX", ctx);
            glen = -1;

            if ((ctx->which_ring == GR_CTX_FMPZ || ctx->which_ring == GR_CTX_FMPQ) &&
                status != GR_SUCCESS)
            {
                flint_printf("FAIL: unable over ZZ or QQ\n\n");
                flint_printf("f = %{gr_poly}\n\n", f, ctx);
                flint_printf("s = %s\n\n", s);
                flint_abort();
            }
        }
        else
        {
            status |= _gr_poly_get_str(&s, f->coeffs, f->length, "XX", ctx);
            glen = n_randint(state, 10);
            gr_poly_fit_length(g, glen, ctx);
            status |= _gr_poly_set_str(g->coeffs, s, "XX", glen, ctx);
            _gr_poly_set_length(g, glen, ctx);
            _gr_poly_normalise(g, ctx);
        }

        if (status == GR_SUCCESS && gr_poly_equal(f, g, ctx) == T_FALSE)
        {
            flint_printf("FAIL: get_str, set_str roundtrip\n\n");
            flint_printf("f = %{gr_poly}\n\n", f, ctx);
            flint_printf("s = %s\n\n", s);
            flint_printf("glen = %wd\n\n", glen);
            flint_printf("g = %{gr_poly}\n\n", g, ctx);
            flint_abort();
        }

        flint_free(s);

        gr_poly_clear(f, ctx);
        gr_poly_clear(g, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
