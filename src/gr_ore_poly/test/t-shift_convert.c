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
#include "gr_poly.h"
#include "gr_ore_poly.h"

static const ore_algebra_t shift_algs[4] = {
    ORE_ALGEBRA_FORWARD_SHIFT, ORE_ALGEBRA_BACKWARD_SHIFT,
    ORE_ALGEBRA_FORWARD_DIFFERENCE, ORE_ALGEBRA_BACKWARD_DIFFERENCE
};

TEST_GR_FUNCTION_START(gr_ore_poly_shift_convert, state, count_success, count_domain, count_unable)
{
    for (slong iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t cctx;
        gr_ore_poly_ctx_t ctx_src, ctx_dst;
        ore_algebra_t src_alg, dst_alg;
        slong len, j, epow, epow2, ngens, var;
        int status = GR_SUCCESS;

        /* occasionally use a multivariate base ring so that the conversion is
           exercised with a generator other than the default one (which is not
           yet implemented and returns GR_UNABLE) */
        switch (n_randint(state, 8))
        {
            case 0:
                gr_ctx_init_random(cctx, state);
                break;
            case 1:
            case 2:
                gr_ctx_init_random_mpoly(cctx, state);
                break;
            default:
                gr_ctx_init_random_poly(cctx, state);
                break;
        }

        ngens = 0;
        status |= gr_ctx_ngens(&ngens, cctx);
        var = (ngens >= 1) ? (slong) n_randint(state, ngens) : 0;

        src_alg = shift_algs[n_randint(state, 4)];
        dst_alg = shift_algs[n_randint(state, 4)];

        gr_ore_poly_ctx_init(ctx_src, cctx, var, src_alg);
        gr_ore_poly_ctx_init(ctx_dst, cctx, var, dst_alg);

        gr_ore_poly_t P, Pd, P2;
        gr_ore_poly_init(P, ctx_src);
        gr_ore_poly_init(Pd, ctx_dst);
        gr_ore_poly_init(P2, ctx_src);

        status |= gr_ore_poly_randtest(P, state, 1 + n_randint(state, 5), ctx_src);
        len = P->length;

        gr_ore_poly_fit_length(Pd, len, ctx_dst);
        status |= _gr_ore_poly_shift_convert(Pd->coeffs, &epow, P->coeffs, len, src_alg, dst_alg, var, cctx);
        _gr_ore_poly_set_length(Pd, len, ctx_dst);
        _gr_ore_poly_normalise(Pd, ctx_dst);

        gr_ore_poly_fit_length(P2, len, ctx_src);
        status |= _gr_ore_poly_shift_convert(P2->coeffs, &epow2, Pd->coeffs, len, dst_alg, src_alg, var, cctx);
        _gr_ore_poly_set_length(P2, len, ctx_src);
        _gr_ore_poly_normalise(P2, ctx_src);

        /* these representative base rings never fail for the implemented var 0 */
        int expect_success = (var == 0 && cctx->which_ring == GR_CTX_GR_POLY
            && (POLYNOMIAL_ELEM_CTX(cctx)->which_ring == GR_CTX_FMPZ
                || POLYNOMIAL_ELEM_CTX(cctx)->which_ring == GR_CTX_CC_ACB
                || POLYNOMIAL_ELEM_CTX(cctx)->which_ring == GR_CTX_NMOD));
        if (expect_success && status != GR_SUCCESS)
        {
            flint_printf("FAIL: shift_convert unexpected failure\n");
            flint_abort();
        }

        /* the round trip recovers the operator: the powers of S cancel */
        if (status == GR_SUCCESS && epow + epow2 != 0)
        {
            flint_printf("FAIL: shift_convert round-trip power\n");
            flint_abort();
        }
        if (status == GR_SUCCESS && gr_ore_poly_equal(P2, P, ctx_src) == T_FALSE)
        {
            flint_printf("FAIL: shift_convert round trip\n");
            flint_abort();
        }

        /* application identity dst(f) = S^{-epow}(src(f)), over polynomial base
           rings (S^epow * dst = src as operators, so the correction is S^{-epow}
           applied on the left of the value) */
        if (cctx->which_ring == GR_CTX_GR_POLY)
        {
            gr_ctx_struct * sctx = POLYNOMIAL_ELEM_CTX(cctx);
            gr_ptr f = gr_heap_init(cctx);
            gr_ptr shifted = gr_heap_init(cctx);
            gr_ptr g1 = gr_heap_init(cctx);
            gr_ptr g2 = gr_heap_init(cctx);
            gr_ptr c = gr_heap_init(sctx);

            status |= gr_set_si(c, -epow, sctx);

            for (j = 0; j < 4; j++)
            {
                status |= gr_randtest(f, state, cctx);
                status |= gr_ore_poly_apply(g1, Pd, f, ctx_dst);
                status |= gr_ore_poly_apply(g2, P, f, ctx_src);
                status |= gr_poly_taylor_shift((gr_poly_struct *) shifted, (gr_poly_struct *) g2, c, sctx);
                if (status == GR_SUCCESS && gr_equal(g1, shifted, cctx) == T_FALSE)
                {
                    flint_printf("FAIL: shift_convert application identity\n");
                    flint_abort();
                }
            }

            gr_heap_clear(f, cctx);
            gr_heap_clear(shifted, cctx);
            gr_heap_clear(g1, cctx);
            gr_heap_clear(g2, cctx);
            gr_heap_clear(c, sctx);
        }

        /* a non-shift algebra type is rejected cleanly */
        if (_gr_ore_poly_shift_convert(P2->coeffs, &epow2, P->coeffs, len, ORE_ALGEBRA_DERIVATIVE, dst_alg, var, cctx) != GR_DOMAIN)
        {
            flint_printf("FAIL: shift_convert expected GR_DOMAIN\n");
            flint_abort();
        }

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        gr_ore_poly_clear(P, ctx_src);
        gr_ore_poly_clear(Pd, ctx_dst);
        gr_ore_poly_clear(P2, ctx_src);
        gr_ore_poly_ctx_clear(ctx_src);
        gr_ore_poly_ctx_clear(ctx_dst);
        gr_ctx_clear(cctx);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_domain, count_unable);
}
