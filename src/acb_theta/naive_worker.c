/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "ulong_extras.h"
#include "acb_mat.h"
#include "acb_theta.h"

static slong
acb_theta_naive_fullprec(const acb_theta_eld_t E, slong prec)
{
    return prec + FLINT_MAX(prec + ceil(n_flog(1 + acb_theta_eld_nb_pts(E), 2)),
        ACB_THETA_LOW_PREC);
}


ACB_INLINE slong
acb_theta_naive_newprec(slong prec, slong coord, slong dist, slong max_dist, slong ord)
{
    double r = ((double) FLINT_MAX(0, dist - 1)) / (max_dist + 2);
    double neg = r * r * prec;
    double pos = ord * n_clog(1 + FLINT_ABS(coord), 2);

    return FLINT_MAX(ACB_THETA_LOW_PREC, ceil((double) prec - neg + pos));
}

/* Call worker in dimension 1: make vectors to use in acb_dot  */

static void
acb_theta_naive_call_dim1(acb_ptr th, acb_ptr v1, acb_ptr v2, slong * precs,
    const acb_t lin, const acb_t lin_inv, const acb_t cf, acb_srcptr sqr_pow,
    const acb_theta_eld_t E, slong ord, slong prec, slong fullprec,
    acb_theta_naive_worker_t worker)
{
    slong *coords;
    slong g = acb_theta_eld_ambient_dim(E);
    slong min = acb_theta_eld_min(E);
    slong mid = acb_theta_eld_mid(E);
    slong max = acb_theta_eld_max(E);
    slong len = acb_theta_eld_nb_pts(E);
    slong k;

    coords = flint_malloc(g * sizeof(slong));
    coords[0] = min;
    for (k = 1; k < g; k++)
    {
        coords[k] = acb_theta_eld_coord(E, k);
    }

    /* Store lin^k in v1 and square powers in v2; then call worker */
    for (k = mid; k <= max; k++)
    {
        precs[k - min] = acb_theta_naive_newprec(prec, k, k - mid, max - mid, ord);
        if ((k == mid) && (k >= 0))
        {
            acb_pow_ui(&v1[mid - min], lin, mid, prec);
        }
        else if ((k == mid) && (k <= 0))
        {
            acb_pow_ui(&v1[mid - min], lin_inv, -mid, prec);
        }
        else if ((k > FLINT_MAX(2 * mid, 0)) && (k % 2 == 0))
        {
            acb_sqr(&v1[k - min], &v1[(k / 2) - min], precs[k - min]);
        }
        else
        {
            acb_mul(&v1[k - min], &v1[k - 1 - min], lin, precs[k - min]);
        }
        acb_set_round(&v2[k - min], &sqr_pow[FLINT_ABS(k)], precs[k - min]);
    }
    for (k = mid - 1; k >= min; k--)
    {
        precs[k - min] = acb_theta_naive_newprec(prec, k, k - mid, max - mid, ord);
        if ((k < FLINT_MIN(2 * mid, 0)) && (k % 2 == 0))
        {
            acb_sqr(&v1[k - min], &v1[(k / 2) - min], precs[k - min]);
        }
        else
        {
            acb_mul(&v1[k - min], &v1[k + 1 - min], lin_inv, precs[k - min]);
        }
        acb_set_round(&v2[k - min], &sqr_pow[FLINT_ABS(k)], precs[k - min]);
    }

    worker(th, v1, v2, precs, len, cf, coords, ord, g, prec, fullprec);

    flint_free(coords);
}

/* Recursive call to smaller dimension; fall back to dim1 when appropriate */

static void
acb_theta_naive_worker_rec(acb_ptr th, acb_ptr v1, acb_ptr v2, slong * precs,
    acb_mat_t lin_pow, acb_mat_t lin_pow_inv, const acb_t cf, acb_srcptr exp_z,
    acb_srcptr exp_z_inv, const acb_mat_t exp_tau, const acb_mat_t exp_tau_inv,
    const acb_ptr * sqr_pow, const acb_theta_eld_t E, slong ord, slong prec,
    slong fullprec, acb_theta_naive_worker_t worker)
{
    slong d = acb_theta_eld_dim(E);
    slong g = acb_theta_eld_ambient_dim(E);
    slong nr = acb_theta_eld_nr(E);
    slong nl = acb_theta_eld_nl(E);
    slong min = acb_theta_eld_min(E);
    slong mid = acb_theta_eld_mid(E);
    slong max = acb_theta_eld_max(E);
    acb_t start_cf, diff_cf, diff_cf_inv, lin_cf, lin_cf_inv, full_cf;
    acb_ptr start_lin_pow, start_lin_pow_inv, diff_lin_pow, diff_lin_pow_inv;
    slong newprec;
    slong k, j, c;

    /* Catch cases: no points in ellipsoid; d=1 */
    if (acb_theta_eld_nb_pts(E) == 0)
    {
        return;
    }
    else if (d == 1)
    {
        acb_init(lin_cf);
        acb_init(lin_cf_inv);

        acb_set(lin_cf, &exp_z[0]);
        acb_set(lin_cf_inv, &exp_z_inv[0]);
        for (k = 1; k < g; k++)
        {
            acb_mul(lin_cf, lin_cf, acb_mat_entry(lin_pow, 0, k), prec);
            acb_mul(lin_cf_inv, lin_cf_inv, acb_mat_entry(lin_pow_inv, 0, k), prec);
        }
        acb_theta_naive_call_dim1(th, v1, v2, precs, lin_cf, lin_cf_inv, cf,
            sqr_pow[0], E, ord, prec, fullprec, worker);

        acb_clear(lin_cf);
        acb_clear(lin_cf_inv);
        return;
    }

    acb_init(start_cf);
    acb_init(diff_cf);
    acb_init(diff_cf_inv);
    acb_init(lin_cf);
    acb_init(full_cf);
    start_lin_pow = _acb_vec_init(d - 1);
    start_lin_pow_inv = _acb_vec_init(d - 1);
    diff_lin_pow = _acb_vec_init(d - 1);
    diff_lin_pow_inv = _acb_vec_init(d - 1);

    /* Set up things for new cofactor */
    acb_set(diff_cf, &exp_z[d - 1]);
    acb_set(diff_cf_inv, &exp_z_inv[d - 1]);
    for (k = d; k < g; k++)
    {
        acb_mul(diff_cf, diff_cf, acb_mat_entry(lin_pow, d - 1, k), prec);
        acb_mul(diff_cf_inv, diff_cf_inv, acb_mat_entry(lin_pow_inv, d - 1, k), prec);
    }
    if (mid >= 0)
    {
        acb_pow_ui(start_cf, diff_cf, mid, prec);
    }
    else
    {
        acb_pow_ui(start_cf, diff_cf_inv, -mid, prec);
    }
    acb_mul(start_cf, start_cf, cf, prec);

    /* Set up things to update entries (k,d) of lin_pow, k < d */
    for (k = 0; k < d - 1; k++)
    {
        acb_set(&diff_lin_pow[k], acb_mat_entry(exp_tau, k, d - 1));
        acb_set(&diff_lin_pow_inv[k], acb_mat_entry(exp_tau_inv, k, d - 1));
        if (mid >= 0)
        {
            acb_pow_ui(&start_lin_pow[k], &diff_lin_pow[k], mid, prec);
            acb_pow_ui(&start_lin_pow_inv[k], &diff_lin_pow_inv[k], mid, prec);
        }
        else
        {
            acb_pow_ui(&start_lin_pow[k], &diff_lin_pow_inv[k], -mid, prec);
            acb_pow_ui(&start_lin_pow_inv[k], &diff_lin_pow[k], -mid, prec);
        }
    }

    /* Right loop */
    acb_set(lin_cf, start_cf);
    for (k = 0; k < d - 1; k++)
    {
        acb_set(acb_mat_entry(lin_pow, k, d - 1), &start_lin_pow[k]);
        acb_set(acb_mat_entry(lin_pow_inv, k, d - 1), &start_lin_pow_inv[k]);
    }
    for (k = 0; k < nr; k++)
    {
        c = mid + k;
        newprec = acb_theta_naive_newprec(prec, c, c - mid, max - mid, ord);
        if (k > 0) /* Update lin_cf, lin_pow using diff */
        {
            for (j = 0; j < d - 1; j++)
            {
                acb_mul(acb_mat_entry(lin_pow, j, d - 1),
                    acb_mat_entry(lin_pow, j, d - 1), &diff_lin_pow[j], newprec);
                acb_mul(acb_mat_entry(lin_pow_inv, j, d - 1),
                    acb_mat_entry(lin_pow_inv, j, d - 1), &diff_lin_pow_inv[j], newprec);
            }
            acb_mul(lin_cf, lin_cf, diff_cf, newprec);
        }

        acb_mul(full_cf, lin_cf, &sqr_pow[d - 1][FLINT_ABS(c)], newprec);
        acb_theta_naive_worker_rec(th, v1, v2, precs, lin_pow, lin_pow_inv, full_cf,
            exp_z, exp_z_inv, exp_tau, exp_tau_inv, sqr_pow, acb_theta_eld_rchild(E, k),
            ord, newprec, fullprec, worker);
    }

    /* Left loop */
    acb_set(lin_cf, start_cf);
    for (k = 0; k < d - 1; k++)
    {
        acb_set(acb_mat_entry(lin_pow, k, d - 1), &start_lin_pow[k]);
        acb_set(acb_mat_entry(lin_pow_inv, k, d - 1), &start_lin_pow_inv[k]);
    }
    for (k = 0; k < nl; k++)
    {
        c = mid - (k + 1);
        newprec = acb_theta_naive_newprec(prec, c, mid - c, mid - min, ord);
        for (j = 0; j < d - 1; j++)
        {
            acb_mul(acb_mat_entry(lin_pow, j, d - 1),
                acb_mat_entry(lin_pow, j, d - 1), &diff_lin_pow_inv[j], newprec);
            acb_mul(acb_mat_entry(lin_pow_inv, j, d - 1),
                acb_mat_entry(lin_pow_inv, j, d - 1), &diff_lin_pow[j], newprec);
        }
        acb_mul(lin_cf, lin_cf, diff_cf_inv, newprec);

        acb_mul(full_cf, lin_cf, &sqr_pow[d - 1][FLINT_ABS(c)], newprec);
        acb_theta_naive_worker_rec(th, v1, v2, precs, lin_pow, lin_pow_inv, full_cf,
            exp_z, exp_z_inv, exp_tau, exp_tau_inv, sqr_pow, acb_theta_eld_lchild(E, k),
            ord, newprec, fullprec, worker);
    }

    acb_clear(start_cf);
    acb_clear(diff_cf);
    acb_clear(diff_cf_inv);
    acb_clear(lin_cf);
    acb_clear(full_cf);
    _acb_vec_clear(start_lin_pow, d - 1);
    _acb_vec_clear(start_lin_pow_inv, d - 1);
    _acb_vec_clear(diff_lin_pow, d - 1);
    _acb_vec_clear(diff_lin_pow_inv, d - 1);
}

static void
acb_theta_naive_precompute(acb_mat_t exp_tau, acb_mat_t exp_tau_inv,
    acb_ptr * sqr_pow, const acb_mat_t tau, const acb_theta_eld_t E, slong prec)
{
    slong g = acb_mat_nrows(tau);
    acb_t c, dc, ddc;
    slong k, j;

    acb_init(c);
    acb_init(dc);
    acb_init(ddc);

    for (k = 0; k < g; k++)
    {
        for (j = k; j < g; j++)
        {
            acb_set(c, acb_mat_entry(tau, k, j));
            if (k != j)
            {
                acb_mul_2exp_si(c, c, 1);
            }
            acb_exp_pi_i(acb_mat_entry(exp_tau, k, j), c, prec);
            acb_inv(acb_mat_entry(exp_tau_inv, k, j), acb_mat_entry(exp_tau, k, j), prec);
        }
    }

    /* Addition chains do not make a huge difference here. */
    for (k = 0; k < g; k++)
    {
        acb_one(c);
        acb_set(dc, acb_mat_entry(exp_tau, k, k));
        acb_sqr(ddc, dc, prec);
        for (j = 0; j <= acb_theta_eld_box(E, k); j++)
        {
            acb_set(&sqr_pow[k][j], c);
            acb_mul(c, c, dc, prec);
            acb_mul(dc, dc, ddc, prec);
        }
    }

    acb_clear(c);
    acb_clear(dc);
    acb_clear(ddc);
}

/* User function */

void
acb_theta_naive_worker(acb_ptr th, slong len, acb_srcptr zs, slong nb,
    const acb_mat_t tau, const acb_theta_eld_t E, slong ord, slong prec,
    acb_theta_naive_worker_t worker)
{
    slong g = acb_theta_eld_ambient_dim(E);
    slong fullprec = acb_theta_naive_fullprec(E, prec);
    slong width = 0;
    acb_mat_t exp_tau, exp_tau_inv, lin_pow, lin_pow_inv;
    acb_ptr * sqr_pow;
    acb_ptr v1, v2, exp_z, exp_z_inv, res;
    slong * precs;
    acb_t cf;
    slong j, k;

    for (j = 0; j < g; j++)
    {
        width = FLINT_MAX(width, 2 * acb_theta_eld_box(E, j) + 1);
    }

    acb_mat_init(exp_tau, g, g);
    acb_mat_init(exp_tau_inv, g, g);
    acb_mat_init(lin_pow, g, g);
    acb_mat_init(lin_pow_inv, g, g);
    sqr_pow = flint_malloc(g * sizeof(acb_ptr));
    for (j = 0; j < g; j++)
    {
        sqr_pow[j] = _acb_vec_init(acb_theta_eld_box(E, j) + 1);
    }
    v1 = _acb_vec_init(width);
    v2 = _acb_vec_init(width);
    exp_z = _acb_vec_init(g);
    exp_z_inv = _acb_vec_init(g);
    res = _acb_vec_init(len * nb);
    acb_init(cf);
    precs = flint_malloc(width * sizeof(slong));

    acb_theta_naive_precompute(exp_tau, exp_tau_inv, sqr_pow, tau, E, prec);
    acb_one(cf);

    for (j = 0; j < nb; j++)
    {
        for (k = 0; k < g; k++)
        {
            acb_mul_2exp_si(&exp_z[k], &zs[j * g + k], 1);
            acb_exp_pi_i(&exp_z[k], &exp_z[k], prec);
            acb_inv(&exp_z_inv[k], &exp_z[k], prec);
        }
        acb_mat_set(lin_pow, exp_tau);
        acb_mat_set(lin_pow_inv, exp_tau_inv);

        acb_theta_naive_worker_rec(res + j * len, v1, v2, precs, lin_pow, lin_pow_inv,
            cf, exp_z, exp_z_inv, exp_tau, exp_tau_inv, sqr_pow, E, ord,
            fullprec, fullprec, worker);
    }
    _acb_vec_set(th, res, len * nb);

    acb_mat_clear(exp_tau);
    acb_mat_clear(exp_tau_inv);
    acb_mat_clear(lin_pow);
    acb_mat_clear(lin_pow_inv);
    for (j = 0; j < g; j++)
    {
        _acb_vec_clear(sqr_pow[j], acb_theta_eld_box(E, j) + 1);
    }
    flint_free(sqr_pow);
    _acb_vec_clear(v1, width);
    _acb_vec_clear(v2, width);
    _acb_vec_clear(exp_z, g);
    _acb_vec_clear(exp_z_inv, g);
    _acb_vec_clear(res, len * nb);
    acb_clear(cf);
    flint_free(precs);
}
