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
    acb_theta_ql_jet(aux, new_zs, nb, new_tau, ord, 0, prec);

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
acb_theta_char_dot_acb(acb_t x, ulong a, acb_srcptr z, slong g, slong prec)
{
    slong * v;
    slong j;

    v = flint_malloc(g * sizeof(slong));

    for (j = 0; j < g; j++)
    {
        v[j] = acb_theta_char_bit(a, j, g);
    }
    acb_dot_si(x, NULL, 0, z, 1, v, 1, g, prec);
    acb_mul_2exp_si(x, x, -1);

    flint_free(v);
}

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
    const acb_mat_t tau, slong ord, ulong ab, int all, int sqr, slong prec)
{
    if (nb <= 0)
    {
        return;
    }
    else if (all && ord == 0 && sqr)
    {
        /* Use duplication formula */
        slong g = acb_mat_nrows(tau);
        slong n = 1 << g;
        acb_ptr new_zs, new_th;
        acb_mat_t new_tau;
        slong j;
        int add_zero;

        add_zero = !_acb_vec_is_zero(zs, g);
        new_zs = _acb_vec_init((nb + add_zero) * g);
        new_th = _acb_vec_init((nb + add_zero) * n);
        acb_mat_init(new_tau, g, g);

        _acb_vec_scalar_mul_2exp_si(new_zs + add_zero * g, zs, nb * g, 1);
        acb_mat_scalar_mul_2exp_si(new_tau, tau, 1);
        acb_theta_ql_jet(new_th, new_zs, nb + add_zero, new_tau,
            0, 0, prec);

        for (j = 0; j < nb; j++)
        {
            acb_theta_agm_mul(th + j * n * n, new_th,
                new_th + (j + add_zero) * n, g, 1, prec);
        }

        _acb_vec_clear(new_zs, (nb + add_zero) * g);
        _acb_vec_clear(new_th, (nb + add_zero) * n);
        acb_mat_clear(new_tau);
    }
    else if (all)
    {
        acb_theta_ql_jet(th, zs, nb, tau, ord, 1, prec);
    }
    else /* just one theta value */
    {
        if (ab == 0)
        {
            acb_theta_jet_notransform_00(th, zs, nb, tau, ord, prec);
        }
        else
        {
            acb_theta_jet_notransform_one(th, zs, nb, tau, ord, ab, prec);
        }
        if (ord == 0 && sqr)
        {
            _acb_vec_sqr(th, th, nb, prec);
        }
    }
}
