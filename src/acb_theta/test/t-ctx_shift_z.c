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

TEST_FUNCTION_START(acb_theta_ctx_shift_z, state)
{
    slong iter;

    /* Test: matches with acb_theta_ctx_set with shifted input */
    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong n = 1 << g;
        slong nb_full = 1 + n_randint(state, 10);
        slong nb = 1 + n_randint(state, nb_full);
        slong start = n_randint(state, nb_full - nb + 1);
        slong prec = 100 + n_randint(state, 200);
        slong mag_bits = n_randint(state, 5);
        ulong a = n_randint(state, n);
        acb_mat_t tau;
        acb_ptr z, z_shift, new_z;
        acb_theta_ctx_t ctx, ctx1, ctx2;
        int same_as;
        slong k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(nb_full * g);
        z_shift = _acb_vec_init(g);
        new_z = _acb_vec_init(nb * g);
        acb_theta_ctx_init(ctx, nb_full, g);
        acb_theta_ctx_init(ctx1, nb, g);
        acb_theta_ctx_init(ctx2, nb, g);

        acb_siegel_randtest_reduced(tau, state, prec, mag_bits);
        acb_siegel_randtest_vec(z, state, nb_full * g, prec);
        acb_theta_char_get_acb(z_shift, a, g);
        acb_mat_vector_mul_col(z_shift, tau, z_shift, prec);
        for (k = 0; k < nb; k++)
        {
            _acb_vec_add(new_z + k * g, z + (start + k) * g, z_shift, g, prec);
        }

        acb_theta_ctx_set_tau(ctx, tau, prec);
        for (k = 0; k < nb_full; k++)
        {
            acb_theta_ctx_set_z(ctx, z + k * g, k, prec);
        }

        acb_theta_ctx_set_tau(ctx1, tau, prec);
        for (k = 0; k < nb; k++)
        {
            acb_theta_ctx_set_z(ctx1, new_z + k * g, k, prec);
        }

        acb_theta_ctx_copy_tau(ctx2, ctx);
        acb_theta_ctx_shift_z(ctx2, ctx, start, nb, a, prec);

        /* Check:
           - information on tau matches
           - if the a vectors are identical, then information on z matches except c and u */
        same_as = _arb_vec_equal(acb_theta_ctx_as(ctx1), acb_theta_ctx_as(ctx2), nb * g);

        if (!acb_mat_overlaps(acb_theta_ctx_tau(ctx1), acb_theta_ctx_tau(ctx2)))
        {
            flint_printf("FAIL (tau)\n");
            acb_mat_printd(acb_theta_ctx_tau(ctx1), 5);
            acb_mat_printd(acb_theta_ctx_tau(ctx2), 5);
            flint_abort();
        }
        if (!arb_mat_overlaps(acb_theta_ctx_y(ctx1), acb_theta_ctx_y(ctx2)))
        {
            flint_printf("FAIL (Y)\n");
            flint_abort();
        }
        if (!arb_mat_overlaps(acb_theta_ctx_yinv(ctx1), acb_theta_ctx_yinv(ctx2)))
        {
            flint_printf("FAIL (Yinv)\n");
            flint_abort();
        }
        if (!acb_mat_overlaps(acb_theta_ctx_exp_tau_div_4(ctx1), acb_theta_ctx_exp_tau_div_4(ctx2)))
        {
            flint_printf("FAIL (exp_tau_div_4)\n");
            flint_abort();
        }
        if (!acb_mat_overlaps(acb_theta_ctx_exp_tau_div_2(ctx1), acb_theta_ctx_exp_tau_div_2(ctx2)))
        {
            flint_printf("FAIL (exp_tau_div_2)\n");
            flint_abort();
        }
        if (!acb_mat_overlaps(acb_theta_ctx_exp_tau(ctx1), acb_theta_ctx_exp_tau(ctx2)))
        {
            flint_printf("FAIL (exp_tau)\n");
            flint_abort();
        }
        if (same_as && !_acb_vec_overlaps(acb_theta_ctx_exp_zs(ctx1), acb_theta_ctx_exp_zs(ctx2), nb * g))
        {
            flint_printf("FAIL (exp_zs)\n");
            _acb_vec_printd(acb_theta_ctx_exp_zs(ctx1), nb * g, 5);
            _acb_vec_printd(acb_theta_ctx_exp_zs(ctx2), nb * g, 5);
            _acb_vec_printd(acb_theta_ctx_exp_zs(ctx) + start * g, nb * g, 5);
            flint_abort();
        }
        if (same_as && !_acb_vec_overlaps(acb_theta_ctx_exp_zs_inv(ctx1), acb_theta_ctx_exp_zs_inv(ctx2), nb * g))
        {
            flint_printf("FAIL (exp_zs_inv)\n");
            flint_abort();
        }
        if (same_as && !_acb_vec_overlaps(acb_theta_ctx_exp_2zs(ctx1), acb_theta_ctx_exp_2zs(ctx2), nb * g))
        {
            flint_printf("FAIL (exp_2zs)\n");
            flint_abort();
        }
        if (same_as && !_acb_vec_overlaps(acb_theta_ctx_exp_2zs_inv(ctx1), acb_theta_ctx_exp_2zs_inv(ctx2), nb * g))
        {
            flint_printf("FAIL (exp_2zs_inv)\n");
            _acb_vec_printd(acb_theta_ctx_exp_2zs(ctx1), nb * g, 5);
            _acb_vec_printd(acb_theta_ctx_exp_2zs_inv(ctx1), nb * g, 5);
            _acb_vec_printd(acb_theta_ctx_exp_2zs_inv(ctx2), nb * g, 5);
            flint_abort();
        }
        if (g > 1 && !arb_mat_overlaps(acb_theta_ctx_cho(ctx1), acb_theta_ctx_cho(ctx2)))
        {
            flint_printf("FAIL (cho)\n");
            flint_abort();
        }
        if (g > 1 && !arb_mat_overlaps(acb_theta_ctx_choinv(ctx1), acb_theta_ctx_choinv(ctx2)))
        {
            flint_printf("FAIL (choinv)\n");
            flint_abort();
        }
        if (g > 1 && !acb_mat_overlaps(acb_theta_ctx_exp_tau_inv(ctx1), acb_theta_ctx_exp_tau_inv(ctx2)))
        {
            flint_printf("FAIL (exp_tau_inv)\n");
            flint_abort();
        }
        if (g > 1 && same_as && !_arb_vec_overlaps(acb_theta_ctx_vs(ctx1), acb_theta_ctx_vs(ctx2), nb * g))
        {
            flint_printf("FAIL (vs)\n");
            acb_mat_printd(tau, 5);
            _acb_vec_printd(z, g, 5);
            _arb_vec_printd(acb_theta_ctx_vs(ctx1), nb * g, 5);
            _arb_vec_printd(acb_theta_ctx_vs(ctx2), nb * g, 5);
            flint_abort();
        }
        /* d0 and d are not copied. */

        acb_mat_clear(tau);
        _acb_vec_clear(z, nb_full * g);
        _acb_vec_clear(z_shift, g);
        _acb_vec_clear(new_z, nb * g);
        acb_theta_ctx_clear(ctx);
        acb_theta_ctx_clear(ctx1);
        acb_theta_ctx_clear(ctx2);
    }

    TEST_FUNCTION_END(state);
}
