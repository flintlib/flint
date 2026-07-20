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
#include "gr_vec.h"
#include "gr_ore_poly.h"

TEST_GR_FUNCTION_START(gr_ore_poly_convert, state, count_success, count_domain, count_unable)
{
    for (slong iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t cctx, cctx2;
        gr_ore_poly_ctx_t ctx_src, ctx_dst;
        ore_algebra_t sa, da;
        slong power = 0, j;
        int independent, status = GR_SUCCESS;

        independent = (n_randint(state, 4) == 0);

        gr_ore_poly_ctx_init_randtest2(cctx, ctx_src, state);
        if (independent)
            gr_ore_poly_ctx_init_randtest2(cctx2, ctx_dst, state);
        else
            gr_ore_poly_ctx_init_randtest(ctx_dst, state, cctx);

        sa = GR_ORE_POLY_CTX(ctx_src)->which_algebra;
        da = GR_ORE_POLY_CTX(ctx_dst)->which_algebra;

        gr_ore_poly_t op, res;
        gr_ore_poly_init(op, ctx_src);
        gr_ore_poly_init(res, ctx_dst);

        status |= gr_ore_poly_randtest(op, state, 1 + n_randint(state, 5), ctx_src);

        status = gr_ore_poly_convert(res, &power, op, ctx_dst, ctx_src);

        int sa_diff = (sa == ORE_ALGEBRA_DERIVATIVE || sa == ORE_ALGEBRA_EULER_DERIVATIVE);
        int da_diff = (da == ORE_ALGEBRA_DERIVATIVE || da == ORE_ALGEBRA_EULER_DERIVATIVE);
        int sa_shift = (sa == ORE_ALGEBRA_FORWARD_SHIFT || sa == ORE_ALGEBRA_BACKWARD_SHIFT
                     || sa == ORE_ALGEBRA_FORWARD_DIFFERENCE || sa == ORE_ALGEBRA_BACKWARD_DIFFERENCE);
        int da_shift = (da == ORE_ALGEBRA_FORWARD_SHIFT || da == ORE_ALGEBRA_BACKWARD_SHIFT
                     || da == ORE_ALGEBRA_FORWARD_DIFFERENCE || da == ORE_ALGEBRA_BACKWARD_DIFFERENCE);
        int expect_success = (!independent && ((sa_diff && da_diff) || (sa_shift && da_shift))
            && cctx->which_ring == GR_CTX_GR_POLY
            && (POLYNOMIAL_ELEM_CTX(cctx)->which_ring == GR_CTX_FMPZ
                || POLYNOMIAL_ELEM_CTX(cctx)->which_ring == GR_CTX_CC_ACB
                || POLYNOMIAL_ELEM_CTX(cctx)->which_ring == GR_CTX_NMOD));
        if (expect_success && status != GR_SUCCESS)
        {
            flint_printf("FAIL: unexpected failure\n");
            flint_printf("sa = %d, da = %d\n", sa, da);
            flint_abort();
        }

        /* the converted operator must act on the base ring consistently with
           the original */

        if (status == GR_SUCCESS)
        {
            gr_ptr f = gr_heap_init(cctx);
            gr_ptr g_src = gr_heap_init(cctx);
            gr_ptr g_dst = gr_heap_init(cctx);
            gr_ptr corrected = gr_heap_init(cctx);

            for (j = 0; j < 4; j++)
            {
                status |= gr_randtest(f, state, cctx);
                status |= gr_ore_poly_apply(g_src, op, f, ctx_src);
                status |= gr_ore_poly_apply(g_dst, res, f, ctx_dst);

                if (status != GR_SUCCESS)
                    continue;

                if (power == 0)
                {
                    status |= gr_set(corrected, g_src, cctx);
                }
                else if (sa_diff && da_diff)
                {
                    /* for currently implemented conversions */
                    FLINT_ASSERT(power <= 0);

                    slong var = GR_ORE_POLY_ORE_DATA(ctx_src)->base_var;
                    gr_vec_t gens;

                    gr_vec_init(gens, 0, cctx);
                    status |= gr_gens(gens, cctx);
                    if (var < gens->length)
                        status |= gr_pow_si(corrected, gr_vec_entry_srcptr(gens, var, cctx), -power, cctx);
                    else
                        status |= GR_UNABLE;
                    gr_vec_clear(gens, cctx);

                    status |= gr_mul(corrected, corrected, g_src, cctx);
                }
                else if (sa_shift && da_shift)
                {
                    gr_ctx_struct * sctx = POLYNOMIAL_ELEM_CTX(cctx);
                    gr_ptr c = gr_heap_init(sctx);

                    status |= gr_set_si(c, -power, sctx);
                    status |= gr_poly_taylor_shift((gr_poly_struct *) corrected,
                                                   (gr_poly_struct *) g_src, c, sctx);
                    gr_heap_clear(c, sctx);
                }

                if (status == GR_SUCCESS && gr_equal(g_dst, corrected, cctx) == T_FALSE)
                {
                    flint_printf("FAIL: application\n");
                    flint_printf("sa = %d, da = %d, power = %wd\n", sa, da, power);
                    flint_abort();
                }
            }

            gr_heap_clear(f, cctx);
            gr_heap_clear(g_src, cctx);
            gr_heap_clear(g_dst, cctx);
            gr_heap_clear(corrected, cctx);
        }

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        gr_ore_poly_clear(op, ctx_src);
        gr_ore_poly_clear(res, ctx_dst);
        gr_ore_poly_ctx_clear(ctx_src);
        gr_ore_poly_ctx_clear(ctx_dst);
        gr_ctx_clear(cctx);
        if (independent)
            gr_ctx_clear(cctx2);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_domain, count_unable);
}
