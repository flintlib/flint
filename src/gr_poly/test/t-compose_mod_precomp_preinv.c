/*
    Copyright (C) 2023, 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr_mat.h"
#include "gr_poly.h"

int
test_compose_mod_precomp_preinv(flint_rand_t state, int which)
{
    gr_ctx_t ctx;
    gr_poly_t A, B, C, D, E, Cinv;
    gr_mat_t M;
    int status = GR_SUCCESS;

    while (1)
    {
        gr_ctx_init_random(ctx, state);
        if (gr_ctx_is_finite(ctx) == T_TRUE || gr_ctx_has_real_prec(ctx) == T_TRUE)
            break;
        gr_ctx_clear(ctx);
    }

    gr_poly_init(A, ctx);
    gr_poly_init(B, ctx);
    gr_poly_init(C, ctx);
    gr_poly_init(D, ctx);
    gr_poly_init(E, ctx);
    gr_poly_init(Cinv, ctx);

    GR_MUST_SUCCEED(gr_poly_randtest(A, state, 1 + n_randint(state, 15), ctx));
    GR_MUST_SUCCEED(gr_poly_randtest(B, state, 1 + n_randint(state, 15), ctx));
    GR_MUST_SUCCEED(gr_poly_randtest(C, state, 1 + n_randint(state, 15), ctx));
    GR_MUST_SUCCEED(gr_poly_randtest(D, state, 1 + n_randint(state, 15), ctx));
    GR_MUST_SUCCEED(gr_poly_randtest(Cinv, state, 1 + n_randint(state, 15), ctx));

    status |= gr_poly_reverse(Cinv, C, C->length, ctx);
    status |= gr_poly_inv_series(Cinv, Cinv, C->length, ctx);

    if (C->length != 0)
    {
        gr_mat_init(M, n_sqrt(C->length-1) + 1, C->length - 1, ctx);
        status |= gr_poly_precompute_matrix(M, B, C, Cinv, ctx);
        status |= gr_poly_rem(A, A, C, ctx);

        switch (which)
        {
            case 0:
                status |= gr_poly_compose_mod_brent_kung_precomp_preinv(D, A, M, C, Cinv, ctx);
                break;
            case 1:
                status |= gr_poly_set(D, A, ctx);
                status |= gr_poly_compose_mod_brent_kung_precomp_preinv(D, D, M, C, Cinv, ctx);
                break;
            case 2:
                status |= gr_poly_set(D, C, ctx);
                status |= gr_poly_compose_mod_brent_kung_precomp_preinv(D, A, M, D, Cinv, ctx);
                break;
            case 3:
                status |= gr_poly_set(D, Cinv, ctx);
                status |= gr_poly_compose_mod_brent_kung_precomp_preinv(D, A, M, C, D, ctx);
                break;

            default:
                flint_abort();
        }

            status |= gr_poly_compose_mod(E, A, B, C, ctx);

        if (status == GR_SUCCESS && gr_poly_equal(D, E, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n\n");
            flint_printf("which = %d\n\n", which);
            gr_ctx_println(ctx);
            flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
            flint_printf("B = "); gr_poly_print(B, ctx); flint_printf("\n");
            flint_printf("C = "); gr_poly_print(C, ctx); flint_printf("\n");
            flint_printf("D = "); gr_poly_print(D, ctx); flint_printf("\n");
            flint_printf("E = "); gr_poly_print(E, ctx); flint_printf("\n");
            flint_abort();
        }

        gr_mat_clear(M, ctx);
    }

    gr_poly_clear(A, ctx);
    gr_poly_clear(B, ctx);
    gr_poly_clear(C, ctx);
    gr_poly_clear(D, ctx);
    gr_poly_clear(E, ctx);
    gr_poly_clear(Cinv, ctx);

    gr_ctx_clear(ctx);

    return status;
}

TEST_FUNCTION_START(gr_poly_compose_mod_precomp_preinv, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        test_compose_mod_precomp_preinv(state, n_randint(state, 4));
    }

    TEST_FUNCTION_END(state);
}
