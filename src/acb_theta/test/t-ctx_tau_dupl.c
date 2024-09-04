/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_mat.h"
#include "acb_theta.h"

static int
acb_theta_ctx_tau_overlaps(const acb_theta_ctx_tau_t ctx1, const acb_theta_ctx_tau_t ctx2)
{
    slong g = ctx2->g;
    slong n = 1 << g;
    int res;

    FLINT_ASSERT(ctx1->g == g);
    FLINT_ASSERT(ctx1->allow_shift == ctx2->allow_shift);

    res = arb_mat_overlaps(&ctx1->yinv, &ctx2->yinv)
        && acb_mat_overlaps(ctx1->exp_tau_div_4, ctx2->exp_tau_div_4)
        && acb_mat_overlaps(ctx1->exp_tau_div_2, ctx2->exp_tau_div_2)
        && acb_mat_overlaps(ctx1->exp_tau, ctx2->exp_tau);

    if (g > 1 && res)
    {
        res = arb_mat_overlaps(&ctx1->cho, &ctx2->cho)
            && acb_mat_overlaps(ctx1->exp_tau_div_4_inv, ctx2->exp_tau_div_4_inv)
            && acb_mat_overlaps(ctx1->exp_tau_div_2_inv, ctx2->exp_tau_div_2_inv)
            && acb_mat_overlaps(ctx1->exp_tau_inv, ctx2->exp_tau_inv);
    }
    if (ctx1->allow_shift && res)
    {
        res = _acb_vec_overlaps(ctx1->exp_tau_a_div_2, ctx2->exp_tau_a_div_2, n * g)
            && _acb_vec_overlaps(ctx1->exp_tau_a, ctx2->exp_tau_a, n * g)
            && _acb_vec_overlaps(ctx1->exp_tau_a_div_2_inv, ctx2->exp_tau_a_div_2_inv, n * g)
            && _acb_vec_overlaps(ctx1->exp_tau_a_inv, ctx2->exp_tau_a_inv, n * g)
            && _acb_vec_overlaps(ctx1->exp_a_tau_a_div_4, ctx2->exp_a_tau_a_div_4, n);
    }

    return res;
}

TEST_FUNCTION_START(acb_theta_ctx_tau_dupl, state)
{
    slong iter;

    /* Test: matches with acb_theta_ctx_tau_set with doubled input */
    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong prec = 100 + n_randint(state, 100);
        slong mag_bits = n_randint(state, 5);
        acb_mat_t tau;
        acb_theta_ctx_tau_t ctx1, ctx2;

        acb_mat_init(tau, g, g);
        acb_theta_ctx_tau_init(ctx1, 1, g);
        acb_theta_ctx_tau_init(ctx2, 1, g);

        acb_siegel_randtest_reduced(tau, state, prec, mag_bits);
        acb_theta_ctx_tau_set(ctx1, tau, prec);
        acb_theta_ctx_tau_dupl(ctx1, prec);

        acb_mat_scalar_mul_2exp_si(tau, tau, 1);
        acb_theta_ctx_tau_set(ctx2, tau, prec);

        if (!acb_theta_ctx_tau_overlaps(ctx1, ctx2))
        {
            flint_printf("FAIL\n");
            flint_abort();
        }

        acb_mat_clear(tau);
        acb_theta_ctx_tau_clear(ctx1);
        acb_theta_ctx_tau_clear(ctx2);
    }

    TEST_FUNCTION_END(state);
}
