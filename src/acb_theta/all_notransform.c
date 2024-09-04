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

static void
acb_theta_all_dupl(acb_ptr th2, acb_srcptr th0, acb_srcptr th, slong g, slong prec)
{
    slong n = 1 << g;
    acb_ptr v;
    ulong a, b;

    v = _acb_vec_init(n);

    for (a = 0; a < n; a++)
    {
        _acb_vec_set(v, th, n);
        for (b = 0; b < n; b++)
        {
            if (acb_theta_char_dot(a, b, g) % 2 == 1)
            {
                acb_neg(&v[b], &v[b]);
            }
        }
        acb_theta_agm_mul(v, th0, v, g, prec);
        for (b = 0; b < n; b++)
        {
            acb_set(&th2[b * n + a], &v[b]);
        }
    }
    _acb_vec_scalar_mul_2exp_si(th2, th2, n * n, g);

    _acb_vec_clear(v, n);
}

static void
acb_theta_all_add_err(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau,
    int sqr, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    slong lp = 8;
    slong nb_der = acb_theta_jet_nb(2, g);
    acb_ptr dth;
    acb_theta_ctx_tau_t ctx_tau;
    acb_theta_ctx_z_struct * vec;
    arb_t err, x;
    slong j, k;

    dth = _acb_vec_init(n * n * nb * nb_der);
    acb_theta_ctx_tau_init(ctx_tau, 0, g);
    vec = acb_theta_ctx_z_vec_init(nb, g);
    arb_init(err);
    arb_init(x);

    /* Attempt to get finite derivatives of theta */
    k = 0;
    do
    {
        k++;
        lp *= 2;
        acb_theta_ctx_tau_set(ctx_tau, tau, lp);
        for (j = 0; j < nb; j++)
        {
            acb_theta_ctx_z_set(&vec[j], zs + j * g, ctx_tau, lp);
        }
        acb_theta_sum_jet_all(dth, vec, nb, ctx_tau, 2, lp);
    }
    while (!_acb_vec_is_finite(dth, n * n * nb * nb_der) && k < 4);

    /* Get error bounds */
    for (j = 0; j < nb; j++)
    {
        for (k = 0; k < n * n; k++)
        {
            acb_theta_jet_error(err, zs + j * g, tau,
                dth + j * (n * n * nb_der) + k * nb_der, 0, lp);
            if (sqr)
            {
                acb_abs(x, &th[j * n * n + k], lp);
                arb_mul_2exp_si(x, x, 1);
                arb_add(x, x, err, lp);
                arb_mul(err, x, err, lp);
            }
            acb_add_error_arb(&th[j * n * n + k], err);
        }
    }

    _acb_vec_clear(dth, n * n * nb * nb_der);
    acb_theta_ctx_tau_clear(ctx_tau);
    acb_theta_ctx_z_vec_clear(vec, nb);
    arb_clear(err);
    arb_clear(x);
}

static void
acb_theta_all_sum(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau,
    int sqr, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_theta_ctx_tau_t ctx_tau;
    acb_theta_ctx_z_struct * vec;
    acb_ptr zero, new_th;
    arb_ptr distances;
    int add_zero;
    slong j;

    FLINT_ASSERT(nb > 0);

    add_zero = sqr && !_acb_vec_is_zero(zs, g);
    acb_theta_ctx_tau_init(ctx_tau, 1, g);
    vec = acb_theta_ctx_z_vec_init(nb + add_zero, g);
    distances = _arb_vec_init(n); /* set to zero */
    zero = _acb_vec_init(g);

    acb_theta_ctx_tau_set(ctx_tau, tau, prec);
    if (add_zero)
    {
        acb_theta_ctx_z_set(&vec[0], zero, ctx_tau, prec);
    }
    for (j = 0; j < nb; j++)
    {
        acb_theta_ctx_z_set(&vec[add_zero + j], zs + j * g, ctx_tau, prec);
    }

    if (!sqr)
    {
        acb_theta_sum_all_tilde(th, vec, nb, ctx_tau, distances, prec);
    }
    else
    {
        new_th = _acb_vec_init((nb + add_zero) * n);

        acb_theta_ctx_tau_dupl(ctx_tau, prec);
        for (j = 0; j < nb + add_zero; j++)
        {
            acb_theta_ctx_z_dupl(&vec[j], prec);
        }

        acb_theta_sum_a0_tilde(new_th, vec, nb + add_zero, ctx_tau, distances, prec);
        for (j = 0; j < nb; j++)
        {
            acb_theta_all_dupl(th + j * n * n, new_th,
                new_th + (add_zero + j) * n, g, prec);
        }

        _acb_vec_clear(new_th, (nb + add_zero) * n);
    }

    for (j = 0; j < nb; j++)
    {
        /* if sqr, then u was squared during duplication. */
        _acb_vec_scalar_mul_arb(th + j * n * n, th + j * n * n, n * n,
            &(&vec[add_zero + j])->u, prec);
    }

    acb_theta_ctx_tau_clear(ctx_tau);
    acb_theta_ctx_z_vec_clear(vec, nb);
    _arb_vec_clear(distances, n);
    _acb_vec_clear(zero, g);
}

static int
acb_theta_all_mid_err(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau,
    const slong * pattern, int sqr, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_mat_t new_tau;
    acb_ptr new_z, new_th;
    arb_mat_t cho, yinv;
    arb_ptr y, w;
    arb_t u, pi;
    int add_zero;
    slong j, k;
    int res;

    FLINT_ASSERT(nb > 0);

    add_zero = !_acb_vec_is_zero(zs, g);
    acb_mat_init(new_tau, g, g);
    new_z = _acb_vec_init((nb + add_zero) * g);
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
        acb_get_mid(&new_z[g * add_zero + j], &zs[j]);
    }

    /* Call ql_exact, with an extra duplication step if sqr is set */
    if (sqr)
    {
        new_th = _acb_vec_init((nb + add_zero) * n);

        /* Duplication */
        acb_mat_scalar_mul_2exp_si(new_tau, new_tau, 1);
        _acb_vec_scalar_mul_2exp_si(new_z, new_z, (nb + add_zero) * g, 1);
        res = acb_theta_ql_exact(new_th, new_z, nb + add_zero, new_tau, pattern, 0, 0, prec);
        if (res)
        {
            for (j = 0; j < nb; j++)
            {
                acb_theta_all_dupl(th + j * n * n, new_th,
                    new_th + (j + add_zero) * n, g, prec);
            }
        }
        acb_mat_scalar_mul_2exp_si(new_tau, new_tau, -1);
        _acb_vec_scalar_mul_2exp_si(new_z, new_z, (nb + add_zero) * g, -1);

        _acb_vec_clear(new_th, (nb + add_zero) * n);
    }
    else
    {
        new_th = _acb_vec_init((nb + add_zero) * n * n);

        res = acb_theta_ql_exact(new_th, new_z, nb + add_zero, new_tau, pattern, 1, 0, prec);

        /*flint_printf("(all_notransform) result of ql_exact: %wd, add_zero = %wd, got values\n",res, add_zero);
        _acb_vec_printd(new_th, (nb + add_zero) * n * n, 5);
        flint_printf("input:\n");
        acb_mat_printd(new_tau, 5);
        _acb_vec_printd(new_z, (nb + add_zero) * g, 5);*/

        if (res)
        {
            _acb_vec_set(th, new_th + add_zero * n * n, nb * n * n);
        }

        _acb_vec_clear(new_th, (nb + add_zero) * n * n);
    }

    if (res)
    {
        /* Multiply by exp(pi y^T Yinv y) */
        acb_siegel_cho_yinv(cho, yinv, tau, prec);
        arb_const_pi(pi, prec);
        for (j = 0; j < nb; j++)
        {
            _acb_vec_get_imag(y, zs + j * g, g);
            arb_mat_vector_mul_col(w, yinv, y, prec);
            arb_dot(u, NULL, 0, y, 1, w, 1, g, prec);
            arb_mul(u, u, pi, prec);
            if (sqr)
            {
                arb_mul_2exp_si(u, u, 1);
            }
            arb_exp(u, u, prec);
            _acb_vec_scalar_mul_arb(th + j * n * n, th + j * n * n, n * n, u, prec);
        }
        /* Add error bounds */
        acb_theta_all_add_err(th, zs, nb, tau, sqr, prec);
    }

    acb_mat_clear(new_tau);
    _acb_vec_clear(new_z, (nb + add_zero) * g);
    arb_mat_clear(cho);
    arb_mat_clear(yinv);
    _arb_vec_clear(y, g);
    _arb_vec_clear(w, g);
    arb_clear(u);
    arb_clear(pi);
    return res;
}

void
acb_theta_all_notransform(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau,
    int sqr, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    slong * pattern;
    slong j;
    int use_sum = 1;
    int res;

    if (nb <= 0)
    {
        return;
    }

    pattern = flint_malloc(g * sizeof(slong));

    res = acb_theta_ql_nb_steps(pattern, tau, prec);

    flint_printf("(all_notransform) pattern:");
    for (j = 0; j < g; j++)
    {
        flint_printf(" %wd", pattern[j]);
    }
    flint_printf("\n");

    if (res && sqr) /* duplication formula means one step less */
    {
        for (j = 0; j < g; j++)
        {
            pattern[j] = FLINT_MAX(0, pattern[j] - 1);
        }
    }

    for (j = 0; j < g; j++)
    {
        if (pattern[j] > 0)
        {
            use_sum = 0;
        }
    }

    if (res && use_sum)
    {
        acb_theta_all_sum(th, zs, nb, tau, sqr, prec);
    }
    else if (res)
    {
        res = acb_theta_all_mid_err(th, zs, nb, tau, pattern, sqr, prec);
    }

    if (!res)
    {
        /* Use sum_bound to avoid returning NaN */
        arb_t c, rho;
        ulong ab;

        arb_init(c);
        arb_init(rho);

        for (j = 0; j < nb; j++)
        {
            acb_theta_sum_bound(c, rho, zs + j * g, tau, 0);
            if (sqr)
            {
                arb_sqr(c, c, prec);
            }
            for (ab = 0; ab < n * n; ab++)
            {
                arb_zero_pm_one(acb_realref(&th[j * n * n + ab]));
                arb_zero_pm_one(acb_imagref(&th[j * n * n + ab]));
            }
            _acb_vec_scalar_mul_arb(th + j * n * n, th + j * n * n, n * n, c, prec);
        }

        arb_clear(c);
        arb_clear(rho);
    }

    flint_free(pattern);
    return;
}
