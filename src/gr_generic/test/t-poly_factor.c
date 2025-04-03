/*
    Copyright (C) 2025 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr.h"
#include "gr_generic.h"
#include "gr_poly.h"
#include "gr_vec.h"

TEST_FUNCTION_START(gr_generic_poly_factor, state)
{
    slong iter;
    slong count_success = 0, count_unable = 0, count_domain = 0;

    for (iter = 0; iter < 1000; iter++)
    {
        gr_ctx_t ctx, pctx, ZZ;
        gr_poly_t pol, f;
        gr_ptr c;
        gr_vec_t fac;
        gr_vec_t mult;

        int status = GR_SUCCESS;

        gr_ctx_init_random(ctx, state);
        gr_ctx_init_gr_poly(pctx, ctx);
        gr_ctx_init_fmpz(ZZ);

        gr_poly_init(pol, ctx);
        gr_poly_init(f, ctx);
        c = gr_heap_init(pctx);
        gr_vec_init(fac, 0, pctx);
        gr_vec_init(mult, 0, ZZ);

        slong nfac = n_randint(state, 3);
        int split = n_randint(state, 2);
        status |= gr_poly_one(pol, ctx);
        for (slong i = 0; i < nfac; i++)
        {
            status |= gr_poly_randtest(f, state,
                                       2 + (split ? 0 : n_randint(state, 3)),
                                       ctx);
            status |= gr_poly_mul(pol, pol, f, ctx);
        }

        status |= gr_generic_poly_factor_roots(c, fac, mult, pol, 0, ctx);

        if (status != GR_SUCCESS)
            goto cleanup;

        for (slong i = 0; i < mult->length; i++)
        {
            gr_ptr g = gr_vec_entry_ptr(fac, i, pctx);
            status |= gr_poly_pow_fmpz(g, g, gr_vec_entry_ptr(mult, i, ZZ), ctx);
            status |= gr_poly_mul(c, c, g, ctx);
        }

        if (status == GR_SUCCESS && gr_equal(pol, c, pctx) == T_FALSE)
        {
            flint_printf("FAIL\npol = %{gr}", pol, pctx);
            flint_abort();
        }

cleanup:

        gr_heap_clear(c, pctx);
        gr_vec_clear(fac, pctx);
        gr_vec_clear(mult, ZZ);
        gr_poly_clear(f, ctx);
        gr_poly_clear(pol, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_unable, count_domain);
}

