/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "arb_mat.h"
#include "acb_mat.h"
#include "acb_theta.h"

void
acb_theta_ctx_set_tau(acb_theta_ctx_t ctx, const acb_mat_t tau, slong prec)
{
    slong g = acb_theta_ctx_g(ctx);
    arb_t pi;
    acb_t x;
    slong j, k;
    int res;

    FLINT_ASSERT(g == acb_mat_nrows(tau));
    arb_init(pi);
    acb_init(x);

    /* Set tau, Y, exp_tau */
    acb_mat_set(acb_theta_ctx_tau(ctx), tau);
    acb_mat_get_imag(acb_theta_ctx_y(ctx), tau);

    for (j = 0; j < g; j++)
    {
        for (k = j; k < g; k++)
        {
            acb_set_round(x, acb_mat_entry(tau, j, k), prec);
            acb_mul_2exp_si(x, x, (k == j ? -2 : -1));
            acb_exp_pi_i(acb_mat_entry(acb_theta_ctx_exp_tau_div_4(ctx), j, k), x, prec);

	    acb_sqr(acb_mat_entry(acb_theta_ctx_exp_tau_div_2(ctx), j, k),
		acb_mat_entry(acb_theta_ctx_exp_tau_div_4(ctx), j, k), prec);
	    acb_sqr(acb_mat_entry(acb_theta_ctx_exp_tau(ctx), j, k),
		acb_mat_entry(acb_theta_ctx_exp_tau_div_2(ctx), j, k), prec);
	}
    }

    /* Set C, Cinv, Yinv, exp_tau_inv */
    if (g == 1)
    {
        arb_inv(arb_mat_entry(acb_theta_ctx_yinv(ctx), 0, 0),
            arb_mat_entry(acb_theta_ctx_y(ctx), 0, 0), prec);
    }
    else
    {
	arb_const_pi(pi, prec);
        acb_siegel_cho(acb_theta_ctx_cho(ctx), tau, prec); /* has a factor sqrt(pi) */
        res = arb_mat_inv(acb_theta_ctx_choinv(ctx), acb_theta_ctx_cho(ctx), prec);
        if (!res)
        {
            arb_mat_indeterminate(acb_theta_ctx_choinv(ctx));
        }
        arb_mat_transpose(acb_theta_ctx_yinv(ctx), acb_theta_ctx_choinv(ctx));
        arb_mat_mul(acb_theta_ctx_yinv(ctx), acb_theta_ctx_choinv(ctx), acb_theta_ctx_yinv(ctx), prec);
	arb_mat_scalar_mul_arb(acb_theta_ctx_yinv(ctx), acb_theta_ctx_yinv(ctx), pi, prec);

        for (j = 0; j < g; j++)
        {
            for (k = j + 1; k < g; k++)
            {
                if (acb_is_real(acb_mat_entry(tau, j, k)))
                {
                    acb_conj(acb_mat_entry(acb_theta_ctx_exp_tau_inv(ctx), j, k),
                        acb_mat_entry(acb_theta_ctx_exp_tau(ctx), j, k));
                }
                else
                {
                    acb_inv(acb_mat_entry(acb_theta_ctx_exp_tau_inv(ctx), j, k),
                        acb_mat_entry(acb_theta_ctx_exp_tau(ctx), j, k), prec);
                }
            }
        }
    }

    arb_clear(pi);
    acb_clear(x);
}
