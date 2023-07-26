/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

static slong
acb_theta_naive_newprec(slong prec, slong coord, slong dist, slong max_dist, slong ord)
{
    double r = ((double) dist) / (max_dist + 2);
    double neg = r * r * prec;
    double pos = ord * n_clog(1 + FLINT_ABS(coord), 2);

    return FLINT_MAX(ACB_THETA_LOW_PREC, ceil((double) prec - neg + pos));
}

/* Work in dimension 1: compute exponentiel terms with two
   multiplications per term only, at just the necessary precision.
   Each term is: cofactor * lin^k * x^(k^2), and square
   powers of x are precomputed. */

static void
acb_theta_naive_worker_dim1(acb_ptr th, const acb_theta_eld_t E,
    const acb_theta_precomp_t D, const acb_t lin, const acb_t cofactor,
    slong ord, slong prec, slong fullprec, acb_theta_naive_worker_t worker_dim0)
{
    acb_t start, diff, aff, term;
    slong *coords;
    slong g = acb_theta_eld_ambient_dim(E);
    slong min = acb_theta_eld_min(E);
    slong mid = acb_theta_eld_mid(E);
    slong max = acb_theta_eld_max(E);
    slong newprec;
    slong k;

    if (acb_theta_eld_nb_pts(E) == 0)
    {
        return;
    }

    acb_init(start);
    acb_init(diff);
    acb_init(aff);
    acb_init(term);
    coords = flint_malloc(g * sizeof(slong));

    for (k = 1; k < g; k++)
    {
        coords[k] = acb_theta_eld_coord(E, k);
    }

    acb_pow_si(start, lin, mid, prec);
    acb_mul(start, start, cofactor, prec);
    acb_set(diff, lin);

    /* Right loop */
    acb_set(aff, start);
    for (k = mid; k <= max; k++)
    {
        coords[0] = k;
        newprec = acb_theta_naive_newprec(prec, k, k - mid, max - mid, ord);
        if (k > mid)
        {
            acb_mul(aff, aff, diff, newprec);
        }

        acb_mul(term, aff, acb_theta_precomp_sqr_pow(D, 0, FLINT_ABS(k)), newprec);
        worker_dim0(th, term, coords, g, ord, newprec, fullprec);
    }

    /* Left loop */
    acb_set(aff, start);
    acb_inv(diff, diff, prec);
    for (k = mid - 1; k >= min; k--)
    {
        coords[0] = k;
        newprec = acb_theta_naive_newprec(prec, k, mid - k, mid - min, ord);
        acb_mul(aff, aff, diff, newprec);

        acb_mul(term, aff, acb_theta_precomp_sqr_pow(D, 0, FLINT_ABS(k)), newprec);
        worker_dim0(th, term, coords, g, ord, newprec, fullprec);
    }

    acb_clear(start);
    acb_clear(diff);
    acb_clear(aff);
    acb_clear(term);
    flint_free(coords);
}

/* Recursive call to smaller dimension; fall back to dim1 when appropriate */

static void
acb_theta_naive_worker_rec(acb_ptr th, acb_mat_t lin_powers,
    const acb_theta_eld_t E, const acb_theta_precomp_t D, acb_srcptr exp_z,
    const acb_t cofactor, slong ord, slong prec, slong fullprec,
    acb_theta_naive_worker_t worker_dim0)
{
    slong d = acb_theta_eld_dim(E);
    slong g = acb_theta_eld_ambient_dim(E);
    slong nr = acb_theta_eld_nr(E);
    slong nl = acb_theta_eld_nl(E);
    slong min = acb_theta_eld_min(E);
    slong mid = acb_theta_eld_mid(E);
    slong max = acb_theta_eld_max(E);
    acb_t start_cf, diff_cf, lin_cf, full_cf;   /* Set up next cofactor */
    acb_ptr start_lin_powers, diff_lin_powers;  /* Set up next lin_powers */
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
        acb_set(lin_cf, &exp_z[0]);
        for (k = 1; k < g; k++)
        {
            acb_mul(lin_cf, lin_cf, acb_mat_entry(lin_powers, 0, k), prec);
        }
        acb_theta_naive_worker_dim1(th, E, D, lin_cf, cofactor, ord, prec,
            fullprec, worker_dim0);
        acb_clear(lin_cf);
        return;
    }

    acb_init(start_cf);
    acb_init(diff_cf);
    acb_init(lin_cf);
    acb_init(full_cf);
    start_lin_powers = _acb_vec_init(d - 1);
    diff_lin_powers = _acb_vec_init(d - 1);

    /* Set up things for new cofactor */
    acb_set(diff_cf, &exp_z[d - 1]);
    for (k = d; k < g; k++)
    {
        acb_mul(diff_cf, diff_cf, acb_mat_entry(lin_powers, d - 1, k), prec);
    }
    acb_pow_si(start_cf, diff_cf, mid, prec);
    acb_mul(start_cf, start_cf, cofactor, prec);

    /* Set up things to update entries (k,d) of lin_powers, k < d */
    for (k = 0; k < d - 1; k++)
    {
        acb_set(&diff_lin_powers[k], acb_mat_entry(acb_theta_precomp_exp_mat(D), k, d - 1));
        acb_pow_si(&start_lin_powers[k], &diff_lin_powers[k], mid, prec);
    }

    /* Right loop */
    acb_set(lin_cf, start_cf);
    for (k = 0; k < d - 1; k++)
    {
        acb_set(acb_mat_entry(lin_powers, k, d - 1), &start_lin_powers[k]);
    }
    for (k = 0; k < nr; k++)
    {
        c = mid + k;
        newprec = acb_theta_naive_newprec(prec, c, c - mid, max - mid, ord);
        if (k > 0) /* Update lin_cf, lin_powers using diff */
        {
            for (j = 0; j < d - 1; j++)
            {
                acb_mul(acb_mat_entry(lin_powers, j, d - 1),
                    acb_mat_entry(lin_powers, j, d - 1), &diff_lin_powers[j], newprec);
            }
            acb_mul(lin_cf, lin_cf, diff_cf, newprec);
        }

        acb_mul(full_cf, lin_cf,
            acb_theta_precomp_sqr_pow(D, d - 1, FLINT_ABS(c)), newprec);
        acb_theta_naive_worker_rec(th, lin_powers, acb_theta_eld_rchild(E, k),
            D, exp_z, full_cf, ord, newprec, fullprec, worker_dim0);
    }

    /* Left loop */
    acb_set(lin_cf, start_cf);
    for (k = 0; k < d - 1; k++)
    {
        acb_set(acb_mat_entry(lin_powers, k, d - 1), &start_lin_powers[k]);
    }
    acb_inv(diff_cf, diff_cf, prec);
    for (k = 0; k < d - 1; k++)
    {
        acb_inv(&diff_lin_powers[k], &diff_lin_powers[k], prec);
    }
    for (k = 0; k < nl; k++)
    {
        c = mid - (k + 1);
        newprec = acb_theta_naive_newprec(prec, c, mid - c, mid - min, ord);
        for (j = 0; j < d - 1; j++)
        {
            acb_mul(acb_mat_entry(lin_powers, j, d - 1),
                acb_mat_entry(lin_powers, j, d - 1), &diff_lin_powers[j], newprec);
        }
        acb_mul(lin_cf, lin_cf, diff_cf, newprec);

        acb_mul(full_cf, lin_cf,
            acb_theta_precomp_sqr_pow(D, d - 1, FLINT_ABS(c)), newprec);
        acb_theta_naive_worker_rec(th, lin_powers, acb_theta_eld_lchild(E, k),
            D, exp_z, full_cf, ord, newprec, fullprec, worker_dim0);
    }

    acb_clear(start_cf);
    acb_clear(diff_cf);
    acb_clear(lin_cf);
    acb_clear(full_cf);
    _acb_vec_clear(start_lin_powers, d - 1);
    _acb_vec_clear(diff_lin_powers, d - 1);
}

/* User function */

void
acb_theta_naive_worker(acb_ptr th, slong nb, const acb_t c, const arb_t u,
    const acb_theta_eld_t E, const acb_theta_precomp_t D, slong k,
    slong ord, slong prec, acb_theta_naive_worker_t worker_dim0)
{
    slong g = acb_theta_eld_ambient_dim(E);
    acb_mat_t lin_powers;
    acb_t cofactor;
    slong j;

    acb_mat_init(lin_powers, g, g);
    acb_init(cofactor);

    acb_mat_set(lin_powers, acb_theta_precomp_exp_mat(D));
    acb_one(cofactor);

    for (j = 0; j < nb; j++)
    {
        acb_zero(&th[j]);
    }

    acb_theta_naive_worker_rec(th, lin_powers, E, D, acb_theta_precomp_exp_z(D, k, 0),
        cofactor, ord, prec, prec, worker_dim0);

    for (j = 0; j < nb; j++)
    {
        acb_mul(&th[j], &th[j], c, prec);
        acb_add_error_arb(&th[j], u);
    }

    acb_mat_clear(lin_powers);
    acb_clear(cofactor);
}
