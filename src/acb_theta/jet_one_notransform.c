/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_mat.h"
#include "acb_modular.h"
#include "acb_theta.h"

void acb_theta_jet_one_notransform(acb_ptr th, acb_srcptr zs, slong nb,
    const acb_mat_t tau, slong ord, ulong ab, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong nbth = acb_theta_jet_nb(ord, g);
    slong j;
    ulong b = ab % (1 << g);
    ulong a = ab >> g;

    FLINT_ASSERT(nb >= 0);
    FLINT_ASSERT(ab >= 0 && ab < (1 << (2 * g)));

    if (g == 1)
    {
        /* call acb_modular_theta_sum directly */
        acb_theta_ctx_tau_t ctx_tau;
        acb_theta_ctx_z_t ctx;
        acb_ptr res;
        arb_ptr r;

        acb_theta_ctx_tau_init(ctx_tau, g);
        acb_theta_ctx_z_init(ctx, g);
        res = _acb_vec_init(4 * nbth);
        r = _arb_vec_init(g);

        acb_theta_ctx_tau_set(ctx_tau, tau, prec);
        for (j = 0; j < nb; j++)
        {
            acb_theta_ctx_z_set(ctx, zs + j * g, ctx_tau, prec);
            /* acb_modular_theta_sum recomputes the inverse of exp_z */
            acb_modular_theta_sum(res, res + nbth, res + 2 * nbth, res + 3 * nbth,
                acb_theta_ctx_exp_z(ctx), acb_theta_ctx_is_real(ctx),
                acb_mat_entry(acb_theta_ctx_exp_tau(ctx_tau), 0, 0), ord + 1, prec);
            if (ab == 0)
            {
                _acb_vec_set(th + j * nbth, res + 2 * nbth, nbth);
            }
            else if (ab == 1)
            {
                _acb_vec_set(th + j * nbth, res + 3 * nbth, nbth);
            }
            else if (ab == 2)
            {
                _acb_vec_set(th + j * nbth, res + nbth, nbth);
            }
            else
            {
                _acb_vec_neg(th + j * nbth, res, nbth);
            }
            if (ab >= 2)
            {
                _acb_vec_scalar_mul(th + j * nbth, th + j * nbth, nbth,
                    acb_mat_entry(acb_theta_ctx_exp_tau_div_4(ctx_tau), 0, 0), prec);
            }
            _acb_vec_scalar_mul(th + j * nbth, th + j * nbth, nbth,
                acb_theta_ctx_c(ctx), prec);

            /* Adjust derivatives due to reduction */
            _arb_vec_neg(r, acb_theta_ctx_r(ctx), g);

            /* flint_printf("(jet_one_notransform) r: ");
               _arb_vec_printd(r, g, 5); */

            _arb_vec_scalar_mul_2exp_si(r, r, g, 1);
            acb_theta_jet_exp_pi_i(res, r, ord, g, prec);
            acb_theta_jet_mul(th + j * nbth, th + j * nbth, res, ord, g, prec);
        }

        acb_theta_ctx_tau_clear(ctx_tau);
        acb_theta_ctx_z_clear(ctx);
        _acb_vec_clear(res, 4 * nbth);
        _arb_vec_clear(r, g);
    }
    else if (ab == 0)
    {
        /* Call jet_00_notransform directly */
        acb_theta_jet_00_notransform(th, zs, nb, tau, ord, prec);
    }
    else
    {
        /* theta_ab(z, tau) = exp(pi i a^T tau a/4) exp(2 pi i a^T (z + b/2))
           theta_00(z + tau a/2 + b/2, tau) */
        acb_ptr new_zs, v, w, aux;
        arb_ptr u;
        acb_t c, x;

        new_zs = _acb_vec_init(nb * g);
        v = _acb_vec_init(g);
        w = _acb_vec_init(g);
        aux = _acb_vec_init(nbth);
        u = _arb_vec_init(g);
        acb_init(c);
        acb_init(x);

        acb_theta_char_get_acb(v, a, g);
        acb_mat_vector_mul_col(v, tau, v, prec); /* tau.a/2 */
        acb_theta_char_get_acb(w, b, g);
        _acb_vec_add(w, v, w, g, prec);
        for (j = 0; j < nb; j++)
        {
            _acb_vec_add(new_zs + j * g, zs + j * g, w, g, prec);
        }

        acb_theta_jet_00_notransform(th, new_zs, nb, tau, ord, prec);

        acb_theta_char_dot_acb(c, a, v, g, prec);
        for (j = 0; j < nb; j++)
        {
            acb_theta_char_get_acb(w, b, g);
            _acb_vec_add(w, w, zs + j * g, g, prec);
            acb_theta_char_dot_acb(x, a, w, g, prec);
            acb_mul_2exp_si(x, x, 1);
            acb_add(x, x, c, prec);
            acb_exp_pi_i(x, x, prec);
            _acb_vec_scalar_mul(th + j * nbth, th + j * nbth, nbth, x, prec);
        }

        acb_theta_char_get_arb(u, a, g);
        _arb_vec_scalar_mul_2exp_si(u, u, g, 1);
        acb_theta_jet_exp_pi_i(aux, u, ord, g, prec);
        for (j = 0; j < nb; j++)
        {
            acb_theta_jet_mul(th + j * nbth, th + j * nbth, aux, ord, g, prec);
        }

        _acb_vec_clear(new_zs, nb * g);
        _acb_vec_clear(aux, nbth);
        _acb_vec_clear(v, g);
        _acb_vec_clear(w, g);
        _arb_vec_clear(u, g);
        acb_clear(c);
        acb_clear(x);
    }
}
