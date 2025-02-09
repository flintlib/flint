/*
    Copyright (C) 2025 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_mat.h"
#include "acb_theta.h"

/* We make a choice between direct summation or jet_notransform_ql */
/* See p-jet_notransform_ql */

static void
acb_theta_jet_notransform_ax(acb_ptr th, acb_srcptr zs, slong nb,
    const acb_mat_t tau, slong ord, int all, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong guard = ACB_THETA_LOW_PREC;
    slong * pattern;
    acb_theta_ctx_tau_t ctx_tau;
    acb_theta_ctx_z_struct * vec;
    slong j;
    int use_sum = 1;

    pattern = flint_malloc(g * sizeof(slong));

    acb_theta_ql_nb_steps(pattern, tau, 0, prec);
    if (pattern[0] >= 8 + ord)
    {
        use_sum = 0;
    }
    if (g >= 2 && pattern[1] >= 6 + ord)
    {
        use_sum = 0;
    }
    if (g >= 3 && pattern[2] >= 7)
    {
        use_sum = 0;
    }

    if (use_sum)
    {
        acb_theta_ctx_tau_init(ctx_tau, 0, g);
        vec = acb_theta_ctx_z_vec_init(nb, g);

        acb_theta_ctx_tau_set(ctx_tau, tau, prec + guard);
        for (j = 0; j < nb; j++)
        {
            acb_theta_ctx_z_set(&vec[j], zs + j * g, ctx_tau, prec + guard);
        }

        acb_theta_sum_jet(th, vec, nb, ctx_tau, ord, 1, all, prec);

        acb_theta_ctx_tau_clear(ctx_tau);
        acb_theta_ctx_z_vec_clear(vec, nb);
    }
    else
    {
        acb_theta_jet_notransform_ql(th, zs, nb, tau, ord, all, prec);
    }

    flint_free(pattern);
}

/* We use the formula:
   theta_00(z, tau) = sum_a theta_{a,0}(2z, 4tau) */

static void
acb_theta_jet_notransform_00(acb_ptr th, acb_srcptr zs, slong nb,
    const acb_mat_t tau, slong ord, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    slong nbjet = acb_theta_jet_nb(ord, g);
    slong * tups;
    acb_ptr new_zs, aux;
    acb_mat_t new_tau;
    slong j, k, a;

    tups = flint_malloc(g * nbjet * sizeof(slong));
    new_zs = _acb_vec_init(nb * g);
    aux = _acb_vec_init(nb * n * nbjet);
    acb_mat_init(new_tau, g, g);

    _acb_vec_scalar_mul_2exp_si(new_zs, zs, nb * g, 1);
    acb_mat_scalar_mul_2exp_si(new_tau, tau, 2);
    acb_theta_jet_notransform_ax(aux, new_zs, nb, new_tau, ord, 0, prec);

    _acb_vec_zero(th, nb * nbjet);
    acb_theta_jet_tuples(tups, ord, g);
    for (j = 0; j < nb; j++)
    {
        for (k = 0; k < nbjet; k++)
        {
            for (a = 0; a < n; a++)
            {
                acb_add(&th[j * nbjet + k], &th[j * nbjet + k],
                    &aux[j * n * nbjet + a * nbjet + k], prec);
            }
            acb_mul_2exp_si(&th[j * nbjet + k], &th[j * nbjet + k],
                acb_theta_jet_total_order(tups + k * g, g));
        }
    }

    flint_free(tups);
    _acb_vec_clear(new_zs, nb * g);
    _acb_vec_clear(aux, nb * n * nbjet);
    acb_mat_clear(new_tau);
}

/* We use the formula:
   theta_ab(z, tau) = exp(pi i a^T tau a/4) exp(2 pi i a^T (z + b/2))
           theta_00(z + tau a/2 + b/2, tau) */

static void
acb_theta_jet_notransform_one(acb_ptr th, acb_srcptr zs, slong nb,
    const acb_mat_t tau, slong ord, ulong ab, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong nbjet = acb_theta_jet_nb(ord, g);
    ulong b = ab % (1 << g);
    ulong a = ab >> g;
    acb_ptr new_zs, v, w, aux;
    arb_ptr u;
    acb_t c, x;
    slong j;

    new_zs = _acb_vec_init(nb * g);
    v = _acb_vec_init(g);
    w = _acb_vec_init(g);
    aux = _acb_vec_init(nbjet);
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

    acb_theta_jet_notransform_00(th, new_zs, nb, tau, ord, prec);

    acb_theta_char_dot_acb(c, a, v, g, prec);
    for (j = 0; j < nb; j++)
    {
        acb_theta_char_get_acb(w, b, g);
        _acb_vec_add(w, w, zs + j * g, g, prec);
        acb_theta_char_dot_acb(x, a, w, g, prec);
        acb_mul_2exp_si(x, x, 1);
        acb_add(x, x, c, prec);
        acb_exp_pi_i(x, x, prec);
        _acb_vec_scalar_mul(th + j * nbjet, th + j * nbjet, nbjet, x, prec);
    }

    if (ord > 0)
    {
        acb_theta_char_get_arb(u, a, g);
        _arb_vec_scalar_mul_2exp_si(u, u, g, 1);
        acb_theta_jet_exp_pi_i(aux, u, ord, g, prec);
        for (j = 0; j < nb; j++)
        {
            acb_theta_jet_mul(th + j * nbjet, th + j * nbjet, aux, ord, g, prec);
        }
    }

    _acb_vec_clear(new_zs, nb * g);
    _acb_vec_clear(v, g);
    _acb_vec_clear(w, g);
    _acb_vec_clear(aux, nbjet);
    _arb_vec_clear(u, g);
    acb_clear(c);
    acb_clear(x);
}

void
acb_theta_jet_notransform(acb_ptr th, acb_srcptr zs, slong nb,
    const acb_mat_t tau, slong ord, ulong ab, int all, slong prec)
{
    if (all)
    {
        acb_theta_jet_notransform_ax(th, zs, nb, tau, ord, 1, prec);
    }
    else if (ab == 0)
    {
        acb_theta_jet_notransform_00(th, zs, nb, tau, ord, prec);
    }
    else
    {
        acb_theta_jet_notransform_one(th, zs, nb, tau, ord, ab, prec);
    }
}
