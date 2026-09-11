/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_poly.h"

/* Products in the cyclic transformed representation agree with
   multiplication modulo x^L - 1. */
TEST_FUNCTION_START(gr_transformed_nmod_poly_cyclic, state)
{
    slong iter;

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx, tctx;
        gr_poly_t A, B, C, D;
        gr_ptr x, y;
        gr_transformed_poly_workload_t wl = {{ 2, 1, 1, 4, 0, 1 }};
        slong L, i, lenw;
        int status = GR_SUCCESS;

        if (n_randint(state, 2))
        {
            gr_ctx_init_nmod(ctx, n_randtest_not_zero(state));
        }
        else
        {
            fmpz_t m;
            fmpz_init(m);
            do {
                fmpz_randbits(m, state, FLINT_BITS + 1 + n_randint(state, 3 * FLINT_BITS));
                fmpz_abs(m, m);
            } while (fmpz_cmp_ui(m, 1) <= 0 || gr_ctx_init_mpn_mod(ctx, m) != GR_SUCCESS);
            fmpz_clear(m);
        }
        gr_poly_init(A, ctx);
        gr_poly_init(B, ctx);
        gr_poly_init(C, ctx);
        gr_poly_init(D, ctx);

        L = 1 + n_randint(state, n_randint(state, 4) ? 100 : 2000);

        /* the cyclic length may be rounded up: bound the accumulation generously */
        if (gr_ctx_init_gr_poly_transformed_cyclic_repr(tctx, ctx, &L, 4096, wl) != GR_SUCCESS)
        {
            gr_poly_clear(A, ctx); gr_poly_clear(B, ctx); gr_poly_clear(C, ctx); gr_poly_clear(D, ctx);
            gr_ctx_clear(ctx);
            continue;
        }

        status |= gr_poly_randtest(A, state, 1 + n_randint(state, L), ctx);
        status |= gr_poly_randtest(B, state, 1 + n_randint(state, L), ctx);

        /* C = A B mod (x^L - 1) by folding the plain product */
        status |= gr_poly_mul(C, A, B, ctx);
        for (i = L; i < C->length; i++)
            status |= gr_add(GR_ENTRY(C->coeffs, i - L, ctx->sizeof_elem),
                             GR_ENTRY(C->coeffs, i - L, ctx->sizeof_elem),
                             GR_ENTRY(C->coeffs, i, ctx->sizeof_elem), ctx);
        status |= gr_poly_truncate(C, C, L, ctx);

        x = gr_heap_init(tctx);
        y = gr_heap_init(tctx);
        status |= _gr_set_gr_poly(x, A->coeffs, A->length, ctx, tctx);
        status |= _gr_set_gr_poly(y, B->coeffs, B->length, ctx, tctx);
        status |= gr_mul(x, x, y, tctx);

        /* windowed and full conversion out */
        lenw = n_randint(state, L + 1);
        gr_poly_fit_length(D, L, ctx);
        status |= _gr_get_gr_poly_window(D->coeffs, x, 0, lenw, ctx, tctx);
        status |= _gr_get_gr_poly_window_destructive(GR_ENTRY(D->coeffs, lenw, ctx->sizeof_elem), x, lenw, L, ctx, tctx);
        _gr_poly_set_length_normalise(D, L, ctx);

        if (status != GR_SUCCESS || gr_poly_equal(C, D, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n\n");
            gr_ctx_println(ctx);
            flint_printf("L = %wd, lenw = %wd, status = %d\n", L, lenw, status);
            flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
            flint_printf("B = "); gr_poly_print(B, ctx); flint_printf("\n");
            flint_printf("C = "); gr_poly_print(C, ctx); flint_printf("\n");
            flint_printf("D = "); gr_poly_print(D, ctx); flint_printf("\n");
            flint_abort();
        }

        gr_heap_clear(x, tctx);
        gr_heap_clear(y, tctx);
        gr_ctx_clear(tctx);
        gr_poly_clear(A, ctx);
        gr_poly_clear(B, ctx);
        gr_poly_clear(C, ctx);
        gr_poly_clear(D, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
