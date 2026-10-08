/*
    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* generated using Claude Opus 4.8 */

#include "test_helpers.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_ore_poly.h"

TEST_GR_FUNCTION_START(gr_ore_poly_ddx_to_euler, state, count_success, count_domain, count_unable)
{
    for (slong iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t cctx;
        gr_ore_poly_ctx_t ctx_d, ctx_e;
        gr_vec_t gens;
        slong len, i, j, sz, ngens, var;
        int status = GR_SUCCESS;

        switch (n_randint(state, 8))
        {
            case 0:
                gr_ctx_init_random(cctx, state);
                break;
            case 1:
                gr_ctx_init_random_mpoly(cctx, state);
                break;
            default:
                gr_ctx_init_random_poly(cctx, state);
                break;
        }

        sz = cctx->sizeof_elem;

        ngens = 0;
        status |= gr_ctx_ngens(&ngens, cctx);
        var = (ngens >= 1) ? (slong) n_randint(state, ngens) : 0;

        gr_ore_poly_ctx_init(ctx_d, cctx, var, ORE_ALGEBRA_DERIVATIVE);
        gr_ore_poly_ctx_init(ctx_e, cctx, var, ORE_ALGEBRA_EULER_DERIVATIVE);

        gr_ore_poly_t P_d, P_e, P_d2, scaled_P;
        gr_ore_poly_init(P_d, ctx_d);
        gr_ore_poly_init(P_e, ctx_e);
        gr_ore_poly_init(P_d2, ctx_d);
        gr_ore_poly_init(scaled_P, ctx_d);

        gr_ptr xpow = gr_heap_init(cctx);
        gr_ptr f = gr_heap_init(cctx);
        gr_ptr g_d = gr_heap_init(cctx);
        gr_ptr g_e = gr_heap_init(cctx);
        gr_ptr scaled = gr_heap_init(cctx);

        /* main run */

        status |= gr_ore_poly_randtest(P_d, state, 1 + n_randint(state, 5), ctx_d);
        len = P_d->length;

        gr_ore_poly_fit_length(P_e, len, ctx_e);
        status |= _gr_ore_poly_ddx_to_euler(P_e->coeffs, P_d->coeffs, len, var, cctx);
        _gr_ore_poly_set_length(P_e, len, ctx_e);
        _gr_ore_poly_normalise(P_e, ctx_e);

        int expect_success = (var == 0 && cctx->which_ring == GR_CTX_GR_POLY
            && (POLYNOMIAL_ELEM_CTX(cctx)->which_ring == GR_CTX_FMPZ
                || POLYNOMIAL_ELEM_CTX(cctx)->which_ring == GR_CTX_CC_ACB
                || POLYNOMIAL_ELEM_CTX(cctx)->which_ring == GR_CTX_NMOD));
        if (expect_success && status != GR_SUCCESS)
        {
            flint_printf("FAIL: unexpected failure\n");
            flint_abort();
        }

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        /* correcting factor */

        gr_vec_init(gens, 0, cctx);
        status |= gr_gens(gens, cctx);
        if (var < gens->length)
            status |= gr_set(xpow, gr_vec_entry_srcptr(gens, var, cctx), cctx);
        else
            status |= GR_UNABLE;
        gr_vec_clear(gens, cctx);
        status |= gr_pow_ui(xpow, xpow, (len >= 1) ? (ulong) (len - 1) : 0, cctx);

        /* test round trip */

        gr_ore_poly_fit_length(P_d2, len, ctx_d);
        status |= _gr_ore_poly_euler_to_ddx(P_d2->coeffs, P_e->coeffs, len, var, cctx);
        _gr_ore_poly_set_length(P_d2, len, ctx_d);
        _gr_ore_poly_normalise(P_d2, ctx_d);

        gr_ore_poly_fit_length(scaled_P, len, ctx_d);
        for (i = 0; i < len; i++)
            status |= gr_mul(GR_ENTRY(scaled_P->coeffs, i, sz), xpow, GR_ENTRY(P_d->coeffs, i, sz), cctx);
        _gr_ore_poly_set_length(scaled_P, len, ctx_d);
        _gr_ore_poly_normalise(scaled_P, ctx_d);

        if (status == GR_SUCCESS && gr_ore_poly_equal(P_d2, scaled_P, ctx_d) == T_FALSE)
        {
            flint_printf("FAIL: round trip\n");
            flint_abort();
        }

        /* test action */

        for (j = 0; j < 4; j++)
        {
            status |= gr_randtest_not_zero(f, state, cctx);
            status |= gr_ore_poly_apply(g_d, P_d, f, ctx_d);
            status |= gr_ore_poly_apply(g_e, P_e, f, ctx_e);

            status |= gr_mul(scaled, xpow, g_d, cctx);
            if (status == GR_SUCCESS && gr_equal(g_e, scaled, cctx) == T_FALSE)
            {
                flint_printf("FAIL: application identity\n");
                flint_abort();
            }
        }

        gr_heap_clear(xpow, cctx);
        gr_heap_clear(f, cctx);
        gr_heap_clear(g_d, cctx);
        gr_heap_clear(g_e, cctx);
        gr_heap_clear(scaled, cctx);
        gr_ore_poly_clear(P_d, ctx_d);
        gr_ore_poly_clear(P_e, ctx_e);
        gr_ore_poly_clear(P_d2, ctx_d);
        gr_ore_poly_clear(scaled_P, ctx_d);
        gr_ore_poly_ctx_clear(ctx_d);
        gr_ore_poly_ctx_clear(ctx_e);
        gr_ctx_clear(cctx);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_domain, count_unable);
}
