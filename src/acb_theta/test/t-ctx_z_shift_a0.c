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

TEST_FUNCTION_START(acb_theta_ctx_z_shift_a0, state)
{
    slong iter;

    /* Test: matches with acb_theta_ctx_z_set with shifted input */
    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        slong g = 2 + n_randint(state, 3);
        slong n = 1 << g;
        slong prec = 100 + n_randint(state, 200);
        slong mag_bits = n_randint(state, 5);
        ulong a = n_randint(state, n);
        acb_mat_t tau;
        acb_ptr z, z_shift;
        acb_theta_ctx_tau_t ctx_tau;
        acb_theta_ctx_z_t ctx1, ctx2;
        acb_t c;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        z_shift = _acb_vec_init(g);
        acb_theta_ctx_tau_init(ctx_tau, 1, g);
        acb_theta_ctx_z_init(ctx1, g);
        acb_theta_ctx_z_init(ctx2, g);
        acb_init(c);

        acb_siegel_randtest_reduced(tau, state, prec, mag_bits);
        acb_siegel_randtest_vec(z, state, g, prec);
        acb_theta_char_get_acb(z_shift, a, g);
        acb_mat_vector_mul_col(z_shift, tau, z_shift, prec);

        acb_theta_ctx_tau_set(ctx_tau, tau, prec);
        acb_theta_ctx_z_set(ctx1, z, ctx_tau, prec);
        acb_theta_ctx_z_shift_a0(ctx1, c, ctx1, ctx_tau, a, prec);

        _acb_vec_add(z, z, z_shift, g, prec);
        acb_theta_ctx_z_set(ctx2, z, ctx_tau, prec);

        /* Check: contexts match except u, uinv */
        if (!_acb_vec_overlaps(ctx1->exp_z, ctx2->exp_z, g))
        {
            flint_printf("FAIL (exp_z)\n");
            flint_abort();
        }
        if (!_acb_vec_overlaps(ctx1->exp_z_inv, ctx2->exp_z_inv, g))
        {
            flint_printf("FAIL (exp_z_inv)\n");
            flint_abort();
        }
        if (!_acb_vec_overlaps(ctx1->exp_2z, ctx2->exp_2z, g))
        {
            flint_printf("FAIL (exp_2z)\n");
            flint_abort();
        }
        if (!_acb_vec_overlaps(ctx1->exp_2z_inv, ctx2->exp_2z_inv, g))
        {
            flint_printf("FAIL (exp_2z_inv)\n");
            flint_abort();
        }
        if (!_arb_vec_overlaps(ctx1->v, ctx2->v, g))
        {
            flint_printf("FAIL (v)\n");
            _arb_vec_printd(ctx1->v, g, 5);
            _arb_vec_printd(ctx2->v, g, 5);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _acb_vec_clear(z_shift, g);
        acb_theta_ctx_tau_clear(ctx_tau);
        acb_theta_ctx_z_clear(ctx1);
        acb_theta_ctx_z_clear(ctx2);
        acb_clear(c);
    }

    TEST_FUNCTION_END(state);
}
