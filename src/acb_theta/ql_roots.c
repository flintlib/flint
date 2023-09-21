/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"
#include "profiler.h"

/* todo: move out? */
static int
_acb_vec_contains_zero(acb_srcptr v, slong n)
{
    slong k;

    for (k = 0; k < n; k++)
    {
        if (acb_contains_zero(&v[k]))
        {
            return 1;
        }
    }

    return 0;
}

static void
acb_theta_ql_roots_step(acb_ptr r, acb_srcptr exp_z, const acb_mat_t exp_tau,
    const arb_mat_t cho, const arb_ptr offsets, const arb_ptr mul_err, slong prec,
    slong* hprecs)
{
    slong g = acb_mat_nrows(exp_tau);
    slong n = 1 << g;
    slong* nb_pts;
    slong** pts;
    acb_theta_eld_t E;
    slong* a_vec;
    arf_t R2;
    arf_t eps;
    arb_ptr err;
    acb_t t;
    slong* exps;
    acb_ptr* powers;
    slong a, k, i, j, e;

    arf_init(R2);
    arf_init(eps);
    err = _arb_vec_init(n);
    acb_init(t);
    a_vec = flint_malloc(g * sizeof(slong));
    pts = flint_malloc(n * sizeof(slong*));
    nb_pts = flint_malloc(n * sizeof(slong));
    exps = flint_malloc((2 * g * g + 2 * g) * sizeof(slong));
    powers = flint_malloc((2 * g * g + 2 * g) * sizeof(acb_ptr));

    /* Collect points */
    for (a = 0; a < n; a++)
    {
        acb_theta_eld_init(E, g, g);
        acb_theta_naive_radius(R2, eps, cho, 0, hprecs[a]);
        acb_theta_eld_fill(E, cho, R2, &offsets[a * g], prec);
        nb_pts[a] = acb_theta_eld_nb_pts(E);
        arb_mul_arf(&err[a], &mul_err[a], eps, prec);

        flint_printf("(ql_roots_step) a = %wd: nb_pts = %wd\n", a, nb_pts[a]);
        pts[a] = flint_malloc(g * nb_pts[a] * sizeof(slong));
        acb_theta_eld_points(pts[a], E);
        acb_theta_eld_clear(E);
    }

    /* For each entry (i,j) in matrix, find out the largest power needed */
    /* Use addition sequences here? */
    for (i = 0; i < g; i++)
    {
        exps[2 * g * g + 2 * i] = 0;
        exps[2 * g * g + 2 * i + 1] = 0;
        for (j = i; j < g; j++)
        {
            exps[2 * (i * g + j)] = 0;
            exps[2 * (i * g + j) + 1] = 0;
        }

        for (a = 0; a < n; a++)
        {
            acb_theta_char_get_slong(a_vec, a, g);
            for (k = 0; k < nb_pts[a]; k++)
            {
                e = (2 * pts[a][k * g + i] + a_vec[i]);
                if (e >= 0)
                {
                    exps[2 * (g * g + i)] = FLINT_MAX(exps[2 * (g * g + i)], e);
                }
                else
                {
                    exps[2 * (g * g + i) + 1] = FLINT_MAX(exps[2 * (g * g + i) + 1], -e);
                }
                for (j = i; j < g; j++)
                {
                    e = (2 * pts[a][k * g + i] + a_vec[i])
                        * (2 * pts[a][k * g + j] + a_vec[j]); /* todo: overflow? */
                    if (e >= 0)
                    {
                        exps[2 * (i * g + j)] = FLINT_MAX(exps[2 * (i * g + j)], e);
                    }
                    else
                    {
                        exps[2 * (i * g + j) + 1] = FLINT_MAX(exps[2 * (i * g + j) + 1], -e);
                    }
                }
            }
        }
        flint_printf("(ql_roots_step) i = %wd, exponents -%wd to %wd\n",
            i, j, exps[2 * (g * g + i) + 1], exps[2 * (g * g + i)]);
    }

    /* Compute powers */
    for (i = 0; i < g; i++)
    {
        for (j = i; j < g; j++)
        {
            powers[2 * (i * g + j)] = _acb_vec_init(exps[2 * (i * g + j)] + 1);
            _acb_vec_set_powers(powers[2 * (i * g + j)], acb_mat_entry(exp_tau, i, j),
                exps[2 * (i * g + j)], prec);
            acb_inv(t, acb_mat_entry(exp_tau, i, j), prec);
            powers[2 * (i * g + j) + 1] = _acb_vec_init(exps[2 * (i * g + j) + 1] + 1);
            _acb_vec_set_powers(powers[2 * (i * g + j) + 1], t,
                exps[2 * (i * g + j) + 1], prec);
        }
        powers[2 * g * g + 2 * i] = _acb_vec_init(exps[2 * g * g + 2 * i] + 1);
        _acb_vec_set_powers(powers[2 * g * g + 2 * i], &exp_z[i], exps[2 * g * g + 2 * i], prec);
        acb_inv(t, &exp_z[i], prec);
        powers[2 * g * g + 2 * i + 1] = _acb_vec_init(exps[2 * g * g + 2 * i + 1] + 1);
        _acb_vec_set_powers(powers[2 * g * g + 2 * i + 1], t, exps[2 * g * g + 2 * i + 1], prec);
    }

    /* For each a, sum at precision hprecs[a] and add error */
    for (a = 0; a < n; a++)
    {
        acb_theta_char_get_slong(a_vec, a, g);
        acb_zero(&r[a]);
        for (k = 0; k < nb_pts[a]; k++)
        {
            acb_one(t);
            for (i = 0; i < g; i++)
            {
                for (j = i; j < g; j++)
                {
                    e = (2 * pts[a][k * g + i] + a_vec[i])
                        * (2 * pts[a][k * g + j] + a_vec[j]);
                    if (e >= 0)
                    {
                        acb_mul(t, t, powers[2 * (i * g + j)] + e, prec);
                    }
                    else
                    {
                        acb_mul(t, t, powers[2 * (i * g + j) + 1] - e, prec);
                    }
                }
                e = (2 * pts[a][k * g + i] + a_vec[i]);
                if (e >= 0)
                {
                    acb_mul(t, t, powers[2 * g * g + 2 * i] + e, prec);
                }
                else
                {
                    acb_mul(t, t, powers[2 * g * g + 2 * i + 1] - e, prec);
                }
            }
            acb_add(&r[a], &r[a], t, hprecs[a]);
        }
        acb_add_error_arb(&r[a], &err[a]);
    }

    arf_clear(R2);
    arf_clear(eps);
    _arb_vec_clear(err, n);
    acb_clear(t);
    flint_free(a_vec);
    for (a = 0; a < n; a++)
    {
        flint_free(pts[a]);
    }
    flint_free(pts);
    flint_free(nb_pts);
    for (i = 0; i < g; i++)
    {
        for (j = i; j < g; j++)
        {
            _acb_vec_clear(powers[2 * (i * g + j)], exps[2 * (i * g + j)] + 1);
            _acb_vec_clear(powers[2 * (i * g + j) + 1], exps[2 * (i * g + j) + 1] + 1);
        }
        _acb_vec_clear(powers[2 * g * g + 2 * i], exps[2 * g * g + 2 * i] + 1);
        _acb_vec_clear(powers[2 * g * g + 2 * i + 1], exps[2 * g * g + 2 * i + 1] + 1);
    }
    flint_free(exps);
    flint_free(powers);
}

static int
acb_theta_ql_new_roots_1(acb_ptr r, acb_srcptr z, arb_srcptr dist,
    const acb_t f, const acb_mat_t tau, slong nb_steps, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_mat_t exp_tau;
    acb_ptr exp_z;
    arb_mat_t cho, Y, Yinv;
    arb_ptr y, t1, t2;
    arb_ptr offsets;
    arb_ptr mul_err;
    acb_t c;
    arb_t u, d;
    slong* hprecs;
    slong i, j, a, k;
    int res = 1;

    acb_mat_init(exp_tau, g, g);
    exp_z = _acb_vec_init(g);
    arb_mat_init(cho, g, g);
    arb_mat_init(Y, g, g);
    arb_mat_init(Yinv, g, g);
    y = _arb_vec_init(g);
    t1 = _arb_vec_init(g);
    t2 = _arb_vec_init(g);
    offsets = _arb_vec_init(n * g);
    mul_err = _arb_vec_init(n);
    acb_init(c);
    arb_init(u);
    arb_init(d);
    hprecs = flint_malloc(n * sizeof(slong));

    /* Initialize things at step 0 */
    for (i = 0; i < g; i++)
    {
        for (j = i; j < g; j++)
        {
            if (i == j)
            {
                acb_mul_2exp_si(c, acb_mat_entry(tau, i, j), -2);
            }
            else
            {
                acb_mul_2exp_si(c, acb_mat_entry(tau, i, j), -1);
            }
            acb_exp_pi_i(acb_mat_entry(exp_tau, i, j), c, prec);
        }
        acb_exp_pi_i(&exp_z[i], &z[i], prec);
    }
    /* Offsets are cho.(Y^{-1} y + a/2) */
    acb_theta_eld_cho(cho, tau, prec);
    acb_mat_get_imag(Y, tau);
    arb_mat_inv(Yinv, Y, prec);
    _acb_vec_get_imag(y, z, g);
    arb_mat_vector_mul_col(y, Yinv, y, prec);

    arb_const_pi(u, prec);
    for (a = 0; a < n; a++)
    {
        acb_theta_char_get_arb(t1, a, g);
        _arb_vec_add(t1, t1, y, g, prec);
        arb_mat_vector_mul_col(&offsets[a * g], cho, t1, prec);

        arb_mat_vector_mul_col(t2, Y, t1, prec);
        arb_dot(&mul_err[a], NULL, 0, t1, 1, t2, 1, g, prec);
        arb_mul(&mul_err[a], &mul_err[a], u, prec);
        arb_exp(&mul_err[a], &mul_err[a], prec);
    }

    arb_sqrt_ui(u, 2, prec);
    for (k = 0; (k < nb_steps) && res; k++)
    {
        /* Update data, set hprecs */
        if (k > 0)
        {
            for (i = 0; i < g; i++)
            {
                for (j = i; j < g; j++)
                {
                    acb_sqr(acb_mat_entry(exp_tau, i, j), acb_mat_entry(exp_tau, i, j), prec);
                }
                acb_sqr(&exp_z[i], &exp_z[i], prec);
            }
            arb_mat_scalar_mul_arb(cho, cho, u, prec);
            _arb_vec_scalar_mul(offsets, offsets, n * g, u, prec);
            for (a = 0; a < n; a++)
            {
                arb_sqr(&mul_err[a], &mul_err[a], prec);
            }
        }
        for (a = 0; a < n; a++)
        {
            arb_mul_2exp_si(d, &dist[a], k);
            hprecs[a] = prec + acb_theta_dist_addprec(d);
        }
        acb_mul_2exp_si(c, f, k);
        acb_exp_pi_i(c, c, prec);

        acb_theta_ql_roots_step(r + k * n, exp_z, exp_tau, cho, offsets, mul_err, prec, hprecs);

        if (_acb_vec_contains_zero(r + k * n, n))
        {
            res = 0;
        }
        _acb_vec_scalar_mul(r + k * n, r + k * n, n, c, prec);
    }

    acb_mat_clear(exp_tau);
    _acb_vec_clear(exp_z, g);
    arb_mat_clear(cho);
    arb_mat_clear(Y);
    arb_mat_clear(Yinv);
    _arb_vec_clear(y, g);
    _arb_vec_clear(t1, g);
    _arb_vec_clear(t2, g);
    _arb_vec_clear(offsets, n * g);
    _arb_vec_clear(mul_err, n);
    acb_clear(c);
    arb_clear(u);
    arb_clear(d);
    flint_free(hprecs);
    return res;
}

static int
acb_theta_ql_roots_1(acb_ptr r, acb_srcptr z, arb_srcptr dist,
    const acb_t f, const acb_mat_t tau, slong nb_steps, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_mat_t w;
    acb_ptr x;
    acb_t c;
    arb_t d;
    slong hprec;
    slong k, a;
    int res = 1;

    acb_mat_init(w, g, g);
    x = _acb_vec_init(g);
    acb_init(c);
    arb_init(d);

    for (k = 0; (k < nb_steps) && res; k++)
    {
        acb_mat_scalar_mul_2exp_si(w, tau, k);
        _acb_vec_scalar_mul_2exp_si(x, z, g, k);
        acb_mul_2exp_si(c, f, k);
        acb_exp_pi_i(c, c, prec);

        for (a = 0; a < n; a++)
        {
            arb_mul_2exp_si(d, &dist[a], k);
            hprec = prec + acb_theta_dist_addprec(d);

            flint_printf("(ql_roots) k = %wd, a = %wd, hprec = %wd:\n",
                k, a, hprec);
            acb_theta_naive_ind(&r[k * n + a], a << g, x, 1, w, hprec);
            if (acb_contains_zero(&r[k * n + a]))
            {
                res = 0;
                break;
            }
        }

        _acb_vec_scalar_mul(r + k * n, r + k * n, n, c, prec);
    }

    acb_mat_clear(w);
    _acb_vec_clear(x, g);
    acb_clear(c);
    arb_clear(d);
    return res;
}

static int
acb_theta_ql_roots_3(acb_ptr r, acb_srcptr t, acb_srcptr z, arb_srcptr dist,
    const acb_mat_t tau, slong nb_steps, slong guard, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    int has_t = !_acb_vec_is_zero(t, g);
    acb_ptr x;
    acb_t f;
    slong k;
    int res = 1;

    x = _acb_vec_init(g);
    acb_init(f);

    acb_theta_ql_log_rescale(f, z, tau, prec);

    if (!has_t)
    {
        res = acb_theta_ql_roots_1(r, z, dist, f, tau, nb_steps, guard);
    }
    else
    {
        for (k = 1; (k < 3) && res; k++)
        {
            _acb_vec_scalar_mul_ui(x, t, g, k, prec);
            _acb_vec_add(x, x, z, g, prec);
            res = acb_theta_ql_roots_1(r + (k - 1) * nb_steps * n, x, dist,
                f, tau, nb_steps, guard);
        }
    }

    _acb_vec_clear(x, g);
    acb_clear(f);
    return res;
}

int
acb_theta_ql_roots(acb_ptr r, acb_srcptr t, acb_srcptr z, arb_srcptr dist0,
    arb_srcptr dist, const acb_mat_t tau, slong nb_steps, slong guard, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    int has_z = !_acb_vec_is_zero(z, g);
    int has_t = !_acb_vec_is_zero(t, g);
    slong nb_r = (has_t ? 2 : 1);
    acb_ptr x;
    int res;

    x = _acb_vec_init(g);

    res = acb_theta_ql_roots_3(r, t, x, dist0, tau, nb_steps, guard, prec);
    if (res && has_z)
    {
        res = acb_theta_ql_roots_3(r + nb_r * n * nb_steps, t, z, dist, tau,
            nb_steps, guard, prec);
    }

    _acb_vec_clear(x, g);
    return res;
}
