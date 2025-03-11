/*
    Copyright (C) 2025 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "ulong_extras.h"
#include "acb.h"
#include "arb_mat.h"
#include "acb_mat.h"
#include "acb_dft.h"
#include "acb_theta.h"

void
acb_theta_ql_inexact(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau,
    acb_srcptr dth, int all, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    slong nbth = (all ? n * n : n);
    slong nbjet = acb_theta_jet_nb(2, g);
    slong lp = ACB_THETA_LOW_PREC;
    acb_mat_t new_tau;
    acb_ptr new_zs, new_th;
    slong * pattern;
    arb_ptr err;
    arb_mat_t cho, yinv;
    arb_ptr y, w;
    arb_t u, pi;
    slong j, k;
    int add_zero;

    FLINT_ASSERT(nb > 0);

    add_zero = !_acb_vec_is_zero(zs, g);
    acb_mat_init(new_tau, g, g);
    new_zs = _acb_vec_init((nb + add_zero) * g);
    new_th = _acb_vec_init((nb + add_zero) * nbth);
    pattern = flint_malloc(g * sizeof(slong));
    err = _arb_vec_init(nb * nbth);
    arb_mat_init(cho, g, g);
    arb_mat_init(yinv, g, g);
    y = _arb_vec_init(g);
    w = _arb_vec_init(g);
    arb_init(u);
    arb_init(pi);

    /* Strip tau, z of error bounds */
    for (j = 0; j < g; j++)
    {
        for (k = j; k < g; k++)
        {
            acb_get_mid(acb_mat_entry(new_tau, j, k), acb_mat_entry(tau, j, k));
            acb_set(acb_mat_entry(new_tau, k, j), acb_mat_entry(new_tau, j, k));
        }
    }
    for (j = 0; j < nb * g; j++)
    {
        acb_get_mid(&new_zs[g * add_zero + j], &zs[j]);
    }

    /* Get error bounds */
    for (j = 0; j < nb; j++)
    {
        for (k = 0; k < nbth; k++)
        {
            acb_theta_ql_jet_error(err + j * nbth + k, zs + j * g, tau,
                dth + j * (nbth * nbjet) + k * nbjet, 0, lp);
        }
    }

    /* Call ql_exact */
    /* Todo: adjust current working precision so that 2^(-prec) is roughly err ? */
    acb_theta_ql_nb_steps(pattern, new_tau, 0, prec);
    acb_theta_ql_exact(new_th, new_zs, nb + add_zero, new_tau, pattern, all, 0, prec);
    _acb_vec_set(th, new_th + add_zero * nbth, nb * nbth);

    /* Multiply by exp(pi y^T Yinv y) */
    acb_siegel_cho_yinv(cho, yinv, tau, prec);
    arb_const_pi(pi, prec);
    for (j = 0; j < nb; j++)
    {
        _acb_vec_get_imag(y, new_zs + (add_zero + j) * g, g);
        arb_mat_vector_mul_col(w, yinv, y, prec);
        arb_dot(u, NULL, 0, y, 1, w, 1, g, prec);
        arb_mul(u, u, pi, prec);
        arb_exp(u, u, prec);
        _acb_vec_scalar_mul_arb(th + j * nbth, th + j * nbth, nbth, u, prec);
    }

    /* Add error bounds */
    for (j = 0; j < nb * nbth; j++)
    {
        acb_add_error_arb(&th[j], &err[j]);
    }

    acb_mat_clear(new_tau);
    _acb_vec_clear(new_zs, (nb + add_zero) * g);
    _acb_vec_clear(new_th, (nb + add_zero) * nbth);
    flint_free(pattern);
    _arb_vec_clear(err, nb * nbth);
    arb_mat_clear(cho);
    arb_mat_clear(yinv);
    _arb_vec_clear(y, g);
    _arb_vec_clear(w, g);
    arb_clear(u);
    arb_clear(pi);
}

/* Given values of f at (x_1 + eps zeta^{n_1}, ..., x_g + eps zeta^{n_g}), make
   Fourier transforms to get Taylor coefficients and add error bound */

static void
acb_theta_jet_finite_diff(acb_ptr dth, const arf_t eps, const arf_t err,
    const arb_t rho, acb_srcptr val, slong ord, slong g, slong prec)
{
    slong nbjet = acb_theta_jet_nb(ord, g);
    slong nbaux = n_pow(ord + 1, g);
    slong lp = ACB_THETA_LOW_PREC;
    acb_ptr aux;
    arb_t t, e;
    slong * tups;
    slong * cyc;
    slong j, i, l;
    slong k;

    aux = _acb_vec_init(nbaux);
    arb_init(t);
    arb_init(e);
    tups = flint_malloc(g * nbjet * sizeof(slong));
    cyc = flint_malloc(g * sizeof(slong));

    for (j = 0; j < g; j++)
    {
        cyc[j] = ord + 1;
    }
    acb_dft_prod(aux, val, cyc, g, prec);
    arb_set_si(t, nbaux);
    _acb_vec_scalar_div_arb(aux, aux, nbaux, t, prec);

    /* Get Taylor coefficients, divide by eps^k, add error */
    acb_theta_jet_tuples(tups, ord, g);
    k = 0;
    arb_one(t);
    arb_pow_ui(e, rho, ord, lp);
    arb_mul_arf(e, e, err, lp);
    for (j = 0; j < nbjet; j++)
    {
        l = 0;
        for (i = 0; i < g; i++)
        {
            l *= ord + 1;
            l += tups[j * g + i];
        }
        acb_set(&dth[j], &aux[l]);

        if (acb_theta_jet_total_order(tups + j * g, g) > k)
        {
            k++;
            arb_mul_arf(t, t, eps, prec);
            arb_pow_ui(e, rho, ord - k, lp);
            arb_mul_arf(e, e, err, lp);
        }
        acb_div_arb(&dth[j], &dth[j], t, prec);
        acb_add_error_arb(&dth[j], e);
    }

    _acb_vec_clear(aux, nbaux);
    arb_clear(t);
    arb_clear(e);
    flint_free(tups);
    flint_free(cyc);
}

static void
acb_theta_jet_finite_diff_radius(arf_t eps, arf_t err, const arb_t c, const arb_t rho,
    slong ord, slong g, slong prec)
{
    slong lp = ACB_THETA_LOW_PREC;
    arb_t x, y;

    arb_init(x);
    arb_init(y);

    /* Set x to min of (1/2g)^(1/b)*rho, (2^(-prec)/2cg)^(1/b)*rho^(2b-1)/b
       where b = ord + 1 */
    arb_set_si(x, 2 * g);
    arb_inv(x, x, lp);
    arb_root_ui(x, x, ord + 1, lp);
    arb_mul(x, x, rho, lp);

    arb_pow_ui(y, rho, 2 * ord + 1, prec);
    arb_mul_2exp_si(y, y, -prec);
    arb_div(y, y, c, lp);
    arb_div_si(y, y, 2 * g, lp);
    arb_root_ui(y, y, ord + 1, lp);

    arb_min(x, x, y, lp);
    arb_get_lbound_arf(eps, x, lp);

    arf_one(err);
    arf_mul_2exp_si(err, err, -prec);

    arb_clear(x);
    arb_clear(y);
}

static void
acb_theta_ql_jet_exact(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau,
    acb_srcptr dth, slong ord, int all, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    slong hprec;
    slong nbth = (all ? n * n : n);
    slong nbaux = n_pow(ord + 1, g);
    slong nbjet = acb_theta_jet_nb(ord, g);
    slong nbjet_2 = acb_theta_jet_nb(ord + 2, g);
    slong nbjet_0 = acb_theta_jet_nb(2, g);
    arb_ptr c, rho, eps, err;
    arb_t t;
    arf_t e;
    acb_ptr zetas, new_zs, new_th, new_dth, dft;
    slong j, k, kmod, l;
    int res = 1;

    c = _arb_vec_init(nb);
    rho = _arb_vec_init(nb);
    eps = _arb_vec_init(nb);
    err = _arb_vec_init(nb);
    arb_init(t);
    arf_init(e);
    zetas = _acb_vec_init(ord + 1);
    new_zs = _acb_vec_init(g * nb * nbaux);
    new_th = _acb_vec_init(nbth * nb * nbaux);
    new_dth = _acb_vec_init(nb * nbth * nbaux * nbjet_0);
    dft = _acb_vec_init(nbaux);

    /* Get bounds and high precision, fail if too large */
    hprec = prec;
    for (j = 0; (j < nb) && res; j++)
    {
        acb_theta_ql_local_bound(&c[j], &rho[j], zs + j * g, tau, ord);
        acb_theta_jet_finite_diff_radius(arb_midref(&eps[j]), arb_midref(&err[j]),
            &c[j], &rho[j], ord, g, prec);

        arb_log_base_ui(t, &eps[j], 2, ACB_THETA_LOW_PREC);
        arb_neg(t, t);

        /* we expect that the second bound holds if finite, but still check */
        res = arb_is_finite(t) && arf_cmpabs_2exp_si(arb_midref(t), 20) <= 0;
        if (res)
        {
            hprec = FLINT_MAX(hprec, prec + ord * (arf_get_si(arb_midref(t), ARF_RND_CEIL) + g));
        }
    }

    /* Fill in new_zs */
    _acb_vec_unit_roots(zetas, ord + 1, ord + 1, hprec);
    for (j = 0; (j < nb) && res; j++)
    {
        for (k = 0; k < nbaux; k++)
        {
            kmod = k;
            for (l = g - 1; l >= 0; l--)
            {
                acb_set(&new_zs[j * g * nbaux + k * g + l], &zetas[kmod % (ord + 1)]);
                kmod = kmod / (ord + 1);
            }
            _acb_vec_scalar_mul_arb(new_zs + j * g * nbaux + k * g,
                new_zs + j * g * nbaux + k * g, g, &eps[j], hprec);
            /* eps[j] should be much smaller than 2^(-prec / (ord + 1)), but still check it. */
            for (l = 0; (l < g) && res; l++)
            {
                acb_get_abs_ubound_arf(e, &new_zs[j * g * nbaux + k * g + l], ACB_THETA_LOW_PREC);
                res = arf_cmp_2exp_si(e, -floor(prec / (ord + 1))) <= 0;
            }
            _acb_vec_add(new_zs + j * g * nbaux + k * g,
                new_zs + j * g * nbaux + k * g, zs + j * g, g, hprec);
        }
    }
    if (res)
    {
        /* Call ql_inexact (as zetas is inexact). We can reuse entries from dth
           here, because all the auxiliary points are contained in the vector
           zs where dth was initially computed. */
        /* If ord = 1 or ord = 3, then zetas is exact. We don't make a special
           case because ql_inexact doesn't have a lot of overhead */
        for (j = 0; j < nb; j++)
        {
            for (k = 0; k < nbaux; k++)
            {
                for (l = 0; l < nbth; l++)
                {
                    _acb_vec_set(new_dth + j * nbth * nbaux * nbjet_0
                        + k * nbth * nbjet_0 + l * nbjet_0,
                        dth + j * nbth * nbjet_2 + l * nbjet_2,
                        nbjet_0);
                }
            }
        }
        acb_theta_ql_inexact(new_th, new_zs, nb * nbaux, tau, new_dth, all, hprec);

        /* Make finite differences */
        for (j = 0; j < nb; j++)
        {
            for (k = 0; k < nbth; k++)
            {
                for (l = 0; l < nbaux; l++)
                {
                    acb_set(&dft[l], &new_th[j * nbth * nbaux + l * nbth + k]);
                }
                acb_theta_jet_finite_diff(th + j * nbth * nbjet + k * nbjet,
                    arb_midref(&eps[j]), arb_midref(&err[j]), &rho[j], dft, ord, g, hprec);
            }
        }
    }
    else
    {
        /* Should not happen in tests */
        _acb_vec_indeterminate(th, nb * nbth * nbjet);
    }

    _arb_vec_clear(c, nb);
    _arb_vec_clear(rho, nb);
    _arb_vec_clear(eps, nb);
    _arb_vec_clear(err, nb);
    arb_clear(t);
    arf_clear(e);
    _acb_vec_clear(zetas, ord + 1);
    _acb_vec_clear(new_zs, g * nb * nbaux);
    _acb_vec_clear(new_th, nbth * nb * nbaux);
    _acb_vec_clear(new_dth, nb * nbth * nbaux * nbjet_0);
    _acb_vec_clear(dft, nbaux);
}

void
acb_theta_ql_jet_fd(acb_ptr th, acb_srcptr zs, slong nb,
    const acb_mat_t tau, slong ord, int all, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    slong nbth = (all ? n * n : n);
    slong nbjet = acb_theta_jet_nb(ord, g);
    slong nbjet_2 = acb_theta_jet_nb(ord + 2, g);
    acb_ptr new_zs;
    acb_theta_ctx_tau_t ctx_tau;
    acb_theta_ctx_z_struct * vec;
    acb_ptr dth;
    arb_ptr err;
    acb_mat_t tau_mid;
    arf_t e;
    slong j, k, lp;
    int res = 0;

    if (nb <= 0)
    {
        return;
    }

    new_zs = _acb_vec_init(nb * g);
    acb_theta_ctx_tau_init(ctx_tau, 0, g);
    vec = acb_theta_ctx_z_vec_init(nb, g);
    dth = _acb_vec_init(nb * nbjet_2 * nbth);
    err = _arb_vec_init(nb * nbjet * nbth);
    acb_mat_init(tau_mid, g, g);
    arf_init(e);

    /* Compute low-precision derivatives of theta */
    /* Add an error of 2^(-prec/(ord + 1)) around z, so that the auxiliary
       evaluation points used in finite differences are contained in new_zs */
    arf_one(e);
    arf_mul_2exp_si(e, e, -floor(prec / (ord + 1)));
    _acb_vec_set(new_zs, zs, nb * g);
    for (j = 0; j < nb * g; j++)
    {
        acb_add_error_arf(&new_zs[j], e);
    }
    for (lp = 8; lp <= 64; lp *= 2)
    {
        acb_theta_ctx_tau_set(ctx_tau, tau, lp + ACB_THETA_LOW_PREC);
        for (j = 0; j < nb; j++)
        {
            acb_theta_ctx_z_set(&vec[j], zs + j * g, ctx_tau, lp + ACB_THETA_LOW_PREC);
        }
        acb_theta_sum_jet(dth, vec, nb, ctx_tau, ord + 2, 1, all, lp);
        if (_acb_vec_is_finite(dth, nb * nbth * nbjet_2))
        {
            res = 1;
            break;
        }
    }

    /* Get error bounds */
    lp = ACB_THETA_LOW_PREC;
    for (j = 0; j < nb; j++)
    {
        for (k = 0; k < nbth; k++)
        {
            acb_theta_ql_jet_error(err + j * nbth * nbjet + k * nbjet, zs + j * g,
                tau, dth + j * nbth * nbjet_2 + k * nbjet_2, ord, lp);
        }
    }
    /* Todo: adjust current working precision so that 2^(-prec) is roughly err ? */

    if (res && ord == 0)
    {
        /* Just call ql_inexact */
        acb_theta_ql_inexact(th, zs, nb, tau, dth, all, prec);
    }
    else if (res)
    {
        /* We want to call ql_jet_exact. Strip tau, z of error bounds */
        for (j = 0; j < g; j++)
        {
            for (k = j; k < g; k++)
            {
                acb_get_mid(acb_mat_entry(tau_mid, j, k), acb_mat_entry(tau, j, k));
                acb_set(acb_mat_entry(tau_mid, k, j), acb_mat_entry(tau_mid, j, k));
            }
        }
        for (j = 0; j < nb * g; j++)
        {
            acb_get_mid(&new_zs[j], &zs[j]);
        }

        /* Call ql_jet_exact and add error */
        acb_theta_ql_jet_exact(th, new_zs, nb, tau_mid, dth, ord, all, prec);
        for (j = 0; j < nb * nbth * nbjet; j++)
        {
            acb_add_error_arb(&th[j], &err[j]);
        }
    }
    else
    {
        /* Should not happen in tests */
        _acb_vec_indeterminate(th, nb * nbth * nbjet);
    }

    _acb_vec_clear(new_zs, nb * g);
    acb_theta_ctx_tau_clear(ctx_tau);
    acb_theta_ctx_z_vec_clear(vec, nb);
    _acb_vec_clear(dth, nb * nbjet_2 * nbth);
    _arb_vec_clear(err, nb * nbjet * nbth);
    acb_mat_clear(tau_mid);
    arf_clear(e);
}
