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
#include "gr_vec.h"
#include "gr_poly.h"

TEST_FUNCTION_START(gr_poly_compose_mod_brent_kung_vec_preinv, state)
{
    slong iter;

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        gr_poly_t B, C, Cinv, D;
        gr_poly_vec_t A, R;
        slong i, l, maxlen;
        int status = GR_SUCCESS;

        if (n_randint(state, 2))
        {
            gr_ctx_init_random_finite_field(ctx, state);
        }
        else
        {
            while (1)
            {
                gr_ctx_init_random_commutative_ring(ctx, state);
                if (gr_ctx_is_finite(ctx) == T_TRUE || gr_ctx_has_real_prec(ctx) == T_TRUE)
                    break;
                gr_ctx_clear(ctx);
            }
        }

        gr_poly_init(B, ctx);
        gr_poly_init(C, ctx);
        gr_poly_init(Cinv, ctx);
        gr_poly_init(D, ctx);

        l = 1 + n_randint(state, 6);
        gr_poly_vec_init(A, l, ctx);
        gr_poly_vec_init(R, l, ctx);

        /* occasionally use long polynomials to exercise rectangular splitting */
        maxlen = (n_randint(state, 10) == 0 && gr_ctx_is_finite(ctx) == T_TRUE) ? 100 : 12;

        GR_MUST_SUCCEED(gr_poly_randtest(C, state, 3 + n_randint(state, maxlen), ctx));

        for (i = 0; i < l; i++)
            GR_MUST_SUCCEED(gr_poly_randtest(A->entries + i, state, n_randint(state, FLINT_MAX(C->length, 1)), ctx));

        GR_MUST_SUCCEED(gr_poly_randtest(B, state, 1 + n_randint(state, maxlen), ctx));

        status |= gr_poly_reverse(Cinv, C, C->length, ctx);
        status |= gr_poly_inv_series(Cinv, Cinv, C->length, ctx);

        status |= gr_poly_compose_mod_brent_kung_vec_preinv(R->entries, A->entries, l, l, B, C, Cinv, ctx);

        if (status == GR_SUCCESS)
        {
            for (i = 0; i < l; i++)
            {
                status |= gr_poly_compose_mod(D, A->entries + i, B, C, ctx);

                if (status == GR_SUCCESS && gr_poly_equal(D, R->entries + i, ctx) == T_FALSE)
                {
                    flint_printf("FAIL\n\n");
                    gr_ctx_println(ctx);
                    flint_printf("A[%wd] = ", i); gr_poly_print(A->entries + i, ctx); flint_printf("\n");
                    flint_printf("B = "); gr_poly_print(B, ctx); flint_printf("\n");
                    flint_printf("C = "); gr_poly_print(C, ctx); flint_printf("\n");
                    flint_printf("D = "); gr_poly_print(D, ctx); flint_printf("\n");
                    flint_printf("R = "); gr_poly_print(R->entries + i, ctx); flint_printf("\n");
                    flint_abort();
                }
            }
        }

        gr_poly_clear(B, ctx);
        gr_poly_clear(C, ctx);
        gr_poly_clear(Cinv, ctx);
        gr_poly_clear(D, ctx);
        gr_poly_vec_clear(A, ctx);
        gr_poly_vec_clear(R, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
