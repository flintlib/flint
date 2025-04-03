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

/* Same as ql_exact when pattern suggests that summation must be used */

static void
acb_theta_ql_exact_sum(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau,
    arb_srcptr distances, int all, int shifted_prec, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    slong nba = (all ? n : 1);
    acb_theta_ctx_tau_t ctx_tau;
    acb_theta_ctx_z_struct * vec;
    slong j;

    FLINT_ASSERT(nb >= 1);

    acb_theta_ctx_tau_init(ctx_tau, 1, g);
    vec = acb_theta_ctx_z_vec_init(nb, g);

    acb_theta_ctx_tau_set(ctx_tau, tau, prec);
    for (j = 0; j < nb; j++)
    {
        acb_theta_ctx_z_set(&vec[j], zs + j * g, ctx_tau, prec);
    }

    if (shifted_prec)
    {
        for (j = 0; j < nb; j++)
        {
            acb_theta_sum(th + j * nba * n, &vec[j], 1, ctx_tau,
                distances + j * n, 1, all, 1, prec);
        }
    }
    else
    {
        /* distances are all set to zero; use one ellipsoid */
        acb_theta_sum(th, vec, nb, ctx_tau, distances, 1, all, 1, prec);
    }

    acb_theta_ctx_tau_clear(ctx_tau);
    acb_theta_ctx_z_vec_clear(vec, nb);
}

/* Same as ql_exact when pattern suggests that duplication must be used */

static void
acb_theta_ql_exact_lower_dim(acb_ptr th, acb_srcptr zs, slong nb,
    const acb_mat_t tau, arb_srcptr distances, slong s,
    const slong * pattern, int all, int shifted_prec, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong nba = 1 << (g - s);
    slong n0 = 1 << s;
    slong n = 1 << g;
    slong nbth = (all ? n * n : n);
    acb_ptr * new_zs;
    acb_ptr * cofactors;
    slong * nb_pts;
    slong ** pts;
    slong * fullprec;
    arf_struct * err;
    acb_ptr th0, z0s;
    slong nb0, nbth0;
    acb_mat_t tau0;
    acb_t x;
    slong j, a, index;
    int res = 1;

    new_zs = flint_malloc(nb * nba * sizeof(acb_ptr));
    cofactors = flint_malloc(nb * nba * sizeof(acb_ptr));
    nb_pts = flint_malloc(nb * nba * sizeof(slong));
    err = flint_malloc(nb * nba * sizeof(arf_struct));
    pts = flint_malloc(nb * nba * sizeof(slong *));
    for (j = 0; j < nb * nba; j++)
    {
        arf_init(&err[j]);
    }
    fullprec = flint_malloc(nb * nba * sizeof(arf_struct));
    acb_mat_window_init(tau0, tau, 0, 0, s, s);
    acb_init(x);

    nb0 = 1; /* always start with zero vector */
    for (j = 0; j < nb; j++)
    {
        for (a = 0; a < nba; a++)
        {
            if (res)
            {
                res = acb_theta_ql_lower_dim(&new_zs[j * nba + a], &cofactors[j * nba + a],
                    &pts[j * nba + a], &nb_pts[j * nba + a], &err[j * nba + a], &fullprec[j * nba + a],
                    zs + j * g, tau, distances + j * n, s, a, prec);
                nb0 += nb_pts[j * nba + a];
            }
            else
            {
                /* Initialize with length 0 to be able to free later. */
                /* Should not happen in tests */
                new_zs[j * nba + a] = _acb_vec_init(0);
                cofactors[j * nba + a] = _acb_vec_init(0);
                nb_pts[j * nba + a] = 0;
                fullprec[j * nba + a] = 0;
                pts[j * nba + a] = flint_malloc(0);
            }
        }
    }

    if (!res)
    {
        /* Should not happen in tests */
        nb0 = 0;
    }
    z0s = _acb_vec_init(nb0 * s);
    nbth0 = (all ? n0 * n0 : n0);
    th0 = _acb_vec_init(nb0 * nbth0);

    if (res)
    {
        /* Put everything together */
        index = 1; /* always start with zero vector */
        for (j = 0; j < nb; j++)
        {
            for (a = 0; a < nba; a++)
            {
                _acb_vec_set(z0s + index * s, new_zs[j * nba + a], nb_pts[j * nba + a] * s);
                index += nb_pts[j * nba + a];
            }
        }

        /* Call acb_theta_ql_exact in dimension s */
        acb_theta_ql_exact(th0, z0s, nb0, tau0, pattern, all, shifted_prec, prec);

        /* Recombine into th */
        index = 1; /* do not use result from zero vector */
        _acb_vec_zero(th, nb * nbth);
        for (j = 0; j < nb; j++)
        {
            for (a = 0; a < nba; a++)
            {
                acb_theta_ql_recombine(th + j * nbth, th0 + index * nbth0, cofactors[j * nba + a],
                    pts[j * nba + a], nb_pts[j * nba + a], &err[j * nba + a], fullprec[j * nba + a],
                    s, a, all, g, prec);
                index += nb_pts[j * nba + a];
            }
        }
    }
    else
    {
        /* Should not happen in tests */
        acb_theta_ql_exact_sum(th, zs, nb, tau, distances, all, shifted_prec, prec);
    }

    /* Clear */
    for (j = 0; j < nb; j++)
    {
        for (a = 0; a < nba; a++)
        {
            _acb_vec_clear(new_zs[j * nba + a], nb_pts[j * nba + a] * s);
            _acb_vec_clear(cofactors[j * nba + a], nb_pts[j * nba + a]);
        }
    }
    flint_free(new_zs);
    flint_free(cofactors);
    flint_free(nb_pts);
    for (j = 0; j < nb * nba; j++)
    {
        arf_clear(&err[j]);
        flint_free(pts[j]);
    }
    flint_free(pts);
    flint_free(err);
    flint_free(fullprec);
    acb_mat_window_clear(tau0);
    _acb_vec_clear(z0s, nb0 * s);
    _acb_vec_clear(th0, nb0 * nbth0);
    acb_clear(x);
}

/* Now assume that pattern suggests to perform duplication steps. Set up
   input vector to acb_theta_ql_perform_steps using sum */
/* We take advantage of the fact that some zs differ by a real value; this is
   why we don't use acb_theta_ql_exact_sum here */

static void
acb_theta_ql_steps_input_sum(acb_ptr th_init, acb_srcptr zs, slong nb,
    acb_srcptr t, const acb_mat_t tau, arb_srcptr distances, slong nb_steps,
    const slong * easy_steps, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_theta_ctx_tau_t ctx_tau;
    acb_theta_ctx_z_struct * aux;
    acb_theta_ctx_z_t ctxt;
    acb_mat_t new_tau;
    acb_ptr new_z;
    arb_ptr d;
    slong j;

    acb_theta_ctx_tau_init(ctx_tau, 1, g);
    aux = acb_theta_ctx_z_vec_init(3, g);
    acb_mat_init(new_tau, g, g);
    new_z = _acb_vec_init(g);
    d = _arb_vec_init(n);
    if (easy_steps[0] < nb_steps)
    {
        acb_theta_ctx_z_init(ctxt, g);
    }

    acb_mat_scalar_mul_2exp_si(new_tau, tau, nb_steps);
    acb_theta_ctx_tau_set(ctx_tau, new_tau, prec);
    if (easy_steps[0] < nb_steps)
    {
        _acb_vec_scalar_mul_2exp_si(new_z, t, g, nb_steps);
        acb_theta_ctx_z_set(ctxt, new_z, ctx_tau, prec);
    }

    for (j = 0; j < nb; j++)
    {
        _acb_vec_scalar_mul_2exp_si(new_z, zs + j * g, g, nb_steps);
        acb_theta_ctx_z_set(&aux[0], new_z, ctx_tau, prec);
        _arb_vec_scalar_mul_2exp_si(d, distances + j * n, n, nb_steps);
        if (easy_steps[j] == nb_steps)
        {
            acb_theta_sum(th_init + 3 * n * j, aux, 1, ctx_tau, d, 1, 0, 1, prec);
        }
        else
        {
            acb_theta_ctx_z_add_real(&aux[1], &aux[0], ctxt, prec);
            acb_theta_ctx_z_add_real(&aux[2], &aux[1], ctxt, prec);
            if (j == 0)
            {
                acb_theta_sum(th_init + 3 * n * j, aux, 3, ctx_tau, d, 1, 0, 1, prec);
            }
            else
            {
                acb_theta_sum(th_init + 3 * n * j + n, aux + 1, 2, ctx_tau, d, 1, 0, 1, prec);
            }
        }
    }

    acb_theta_ctx_tau_clear(ctx_tau);
    acb_theta_ctx_z_vec_clear(aux, 3);
    acb_mat_clear(new_tau);
    _acb_vec_clear(new_z, g);
    _arb_vec_clear(d, n);
    if (easy_steps[0] < nb_steps)
    {
        acb_theta_ctx_z_clear(ctxt);
    }
}

/* Still assume that pattern suggests to perform duplication steps. Set up
   input vector to acb_theta_ql_perform_steps using dimension-lowering */
/* We don't take advantage of the fact that some zs differ by a real value,
   in contrast with acb_theta_ql_steps_input_sum */

static void
acb_theta_ql_steps_input_lower_dim(acb_ptr th_init, acb_srcptr zs, slong nb,
    acb_srcptr t, const acb_mat_t tau, arb_srcptr distances, slong split,
    slong nb_steps, const slong * pattern, const slong * easy_steps, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_ptr new_th, new_z;
    acb_mat_t new_tau;
    arb_ptr new_distances;
    slong new_nb, add;
    slong * new_pattern;
    slong j, k;

    /* Count number of vectors */
    new_nb = 0;
    for (j = 0; j < nb; j++)
    {
        if (j == 0 && easy_steps[j] < nb_steps)
        {
            new_nb += 3;
        }
        else if (easy_steps[j] < nb_steps)
        {
            new_nb += 2;
        }
        else
        {
            new_nb += 1;
        }
    }

    acb_mat_init(new_tau, g, g);
    new_th = _acb_vec_init(n * new_nb);
    new_z = _acb_vec_init(g * new_nb);
    new_distances = _arb_vec_init(n * new_nb);
    new_pattern = flint_malloc(split * sizeof(slong));

    /* Set up input of ql_exact_lower_dim */
    for (j = 0; j < split; j++)
    {
        new_pattern[j] = FLINT_MAX(0, pattern[j] - nb_steps);
    }
    new_nb = 0;
    for (j = 0; j < nb; j++)
    {
        if (j == 0 && easy_steps[j] < nb_steps)
        {
            _acb_vec_zero(new_z + new_nb * g, g);
            _acb_vec_set(new_z + (new_nb + 1) * g, t, g);
            _acb_vec_scalar_mul_2exp_si(new_z + (new_nb + 2) * g, t, g, 1);
            add = 3;
        }
        else if (easy_steps[j] < nb_steps)
        {
            _acb_vec_add(new_z + new_nb * g, zs + j * g, t, g, prec);
            _acb_vec_add(new_z + (new_nb + 1) * g, new_z + new_nb * g, t, g, prec);
            add = 2;
        }
        else
        {
            _acb_vec_set(new_z + new_nb * g, zs + j * g, g);
            add = 1;
        }
        for (k = 0; k < add; k++)
        {
            _arb_vec_set(new_distances + (new_nb + k) * n, distances + j * n, n);
        }
        new_nb += add;
    }

    acb_mat_scalar_mul_2exp_si(new_tau, tau, nb_steps);
    _acb_vec_scalar_mul_2exp_si(new_z, new_z, new_nb * g, nb_steps);
    _arb_vec_scalar_mul_2exp_si(new_distances, new_distances, new_nb * n, nb_steps);

    /* Call ql_exact_lower_dim */
    acb_theta_ql_exact_lower_dim(new_th, new_z, new_nb, new_tau, new_distances,
        split, new_pattern, 0, 1, prec);

    /* Set up th_init from computed data */
    new_nb = 0;
    for (j = 0; j < nb; j++)
    {
        if (j == 0 && easy_steps[j] < nb_steps)
        {
            _acb_vec_set(th_init, new_th, 3 * n);
            new_nb += 3;
        }
        else if (easy_steps[j] < nb_steps)
        {
            _acb_vec_set(th_init + 3 * j * n + n, new_th + new_nb * n, 2 * n);
            new_nb += 2;
        }
        else
        {
            _acb_vec_set(th_init + 3 * j * n, new_th + new_nb * n, n);
            new_nb += 1;
        }
    }

    acb_mat_clear(new_tau);
    _acb_vec_clear(new_th, new_nb * n);
    _acb_vec_clear(new_z, new_nb * g);
    _arb_vec_clear(new_distances, new_nb * n);
    flint_free(new_pattern);
}

static void
acb_theta_ql_step_1(acb_ptr res, acb_srcptr th0, acb_srcptr th, acb_srcptr rts,
    arb_srcptr d0, arb_srcptr d, slong g, slong prec)
{
    slong n = 1 << g;

    acb_theta_agm_mul_tight(res, th0, th, d0, d, g, 0, prec);
    acb_theta_agm_sqrt(res, res, rts, n, prec);
}

static void
acb_theta_ql_perform_steps(acb_ptr th, acb_ptr th_init, acb_srcptr rts,
    acb_srcptr rts_all, slong nb, slong nb_steps, arb_srcptr distances,
    const slong * easy_steps, int all, slong g, slong prec)
{
    slong n = 1 << g;
    acb_ptr th_next, th_temp, aux;
    arb_ptr d0, d;
    slong j, k, a;

    FLINT_ASSERT(nb_steps >= 1);

    th_next = _acb_vec_init(3 * n * nb);
    d = _arb_vec_init(n);
    d0 = _arb_vec_init(n);
    if (all)
    {
        aux = _acb_vec_init(n * n);
    }

    for (k = nb_steps - 1; k >= 0; k--)
    {
        _arb_vec_scalar_mul_2exp_si(d0, distances, n, k + 1);
        for (j = 0; j < nb; j++)
        {
            _arb_vec_scalar_mul_2exp_si(d, distances + j * n, n, k + 1);
            if (k == 0 && all && easy_steps[j] > k)
            {
                /* Compute all theta_ab with an easy step */
                acb_theta_agm_mul_tight(th + j * n * n, th_init,
                    th_init + 3 * j * n, d0, d, g, 1, prec);
                acb_theta_agm_sqrt(th + j * n * n, th + j * n * n,
                    rts_all + j * n * n, n * n, prec);
            }
            else if (k == 0 && all)
            {
                /* Compute all theta_ab with a harder step */
                /* Store theta_ab(2t, tau) in aux */
                acb_theta_agm_mul_tight(aux, th_init,
                    th_init + 3 * j * n + 2 * n, d0, d, g, 1, prec);
                acb_theta_agm_sqrt(aux, aux, rts_all + j * n * n, n * n, prec);
                acb_theta_agm_mul_tight(th + j * n * n, th_init + n,
                    th_init + 3 * j * n + n, d0, d, g, 1, prec);
                for (a = 0; a < n * n; a++)
                {
                    acb_div(&th[j * n * n + a], &th[j * n * n + a],
                        &aux[a], prec);
                }
            }
            else if (easy_steps[j] > k)
            {
                /* theta(z, tau)^2 = sum of theta(0, 2tau) theta(2z, 2tau) */
                acb_theta_ql_step_1(th_next + 3 * j * n, th_init,
                    th_init + 3 * j * n, rts + j * (3 * n * nb_steps) + k * (3 * n),
                    d0, d, g, prec);
            }
            else
            {
                if (k > 0)
                {
                    /* theta(z + t, tau)^2 = sum of theta(0, 2tau) theta(2z + 2t, 2tau) */
                    acb_theta_ql_step_1(th_next + 3 * j * n + n, th_init,
                        th_init + 3 * j * n + n,
                        rts + j * (3 * n * nb_steps) + k * (3 * n) + n, d0, d, g, prec);
                }
                /* theta(z + 2t, tau)^2 = sum of theta(0, 2tau) theta(2z + 4t, 2tau) */
                acb_theta_ql_step_1(th_next + 3 * j * n + 2 * n, th_init,
                    th_init + 3 * j * n + 2 * n,
                    rts + j * (3 * n * nb_steps) + k * (3 * n) + 2 * n, d0, d, g, prec);
                if ((easy_steps[j] == k) || (j == 0))
                {
                    /* theta(z, tau) theta(z + 2t, tau)
                       = sum of theta(2t, 2tau) theta(2z + 2t, 2tau) */
                    acb_theta_agm_mul_tight(th_next + 3 * j * n, th_init + n,
                        th_init + 3 * j * n + n, d0, d, g, 0, prec);
                    for (a = 0; a < n; a++)
                    {
                        acb_div(&th_next[3 * j * n + a], &th_next[3 * j * n + a],
                            &th_next[3 * j * n + 2 * n + a], prec);
                    }
                }
            }
        }
        /* Swap th_next and th_init */
        th_temp = th_init;
        th_init = th_next;
        th_next = th_temp;
    }

    if (!all)
    {
        for (j = 0; j < nb; j++)
        {
            _acb_vec_set(th + j * n, th_init + 3 * j * n, n);
        }
    }

    if (nb_steps % 2 == 0)
    {
        _acb_vec_clear(th_next, 3 * n * nb);
    }
    else
    {
        _acb_vec_clear(th_init, 3 * n * nb);
    }
    _arb_vec_clear(d, n);
    _arb_vec_clear(d0, n);
    if (all)
    {
        _acb_vec_clear(aux, n * n);
    }
}

static void
acb_theta_ql_exact_steps(acb_ptr th, acb_srcptr zs, slong nb,
    const acb_mat_t tau, arb_srcptr distances, slong split,
    const slong * pattern, int all, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    slong nb_steps = pattern[g - 1];
    acb_ptr rts, th_init, t, rts_all;
    slong * easy_steps;
    acb_mat_t new_tau;
    slong guard;
    int res;

    FLINT_ASSERT(nb_steps > 0);

    rts = _acb_vec_init(3 * n * nb * nb_steps);
    t = _acb_vec_init(g);
    th_init = _acb_vec_init(3 * n * nb);
    easy_steps = flint_malloc(nb * sizeof(slong));
    acb_mat_init(new_tau, g, g);
    if (all)
    {
        rts_all = _acb_vec_init(nb * n * n);
    }

    res = acb_theta_ql_setup(rts, rts_all, t, &guard, easy_steps, zs, nb, tau, distances,
        nb_steps, all, prec);

    if (res)
    {
        if (split > 0)
        {
            acb_theta_ql_steps_input_lower_dim(th_init, zs, nb, t, tau,
                distances, split, nb_steps, pattern, easy_steps, prec + guard);
        }
        else
        {
            acb_theta_ql_steps_input_sum(th_init, zs, nb, t, tau, distances,
                nb_steps, easy_steps, prec + guard);
        }

        acb_theta_ql_perform_steps(th, th_init, rts, rts_all, nb, nb_steps,
            distances, easy_steps, all, g, prec + guard);
    }
    else
    {
        /* Setup or intput did not succeed, fall back to summation */
        /* Should not happen in tests */
        acb_theta_ql_exact_sum(th, zs, nb, tau, distances, all, 1, prec);
    }

    _acb_vec_clear(rts, 3 * n * nb * nb_steps);
    _acb_vec_clear(t, g);
    _acb_vec_clear(th_init, 3 * n * nb);
    flint_free(easy_steps);
    acb_mat_clear(new_tau);
    if (all)
    {
        _acb_vec_clear(rts_all, nb * n * n);
    }
}

void
acb_theta_ql_exact(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau,
    const slong * pattern, int all, int shifted_prec, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    arb_ptr distances;
    slong nb_steps, split;
    slong lp = ACB_THETA_LOW_PREC;

    FLINT_ASSERT(nb >= 1);
    FLINT_ASSERT(_acb_vec_is_zero(zs, g));

    distances = _arb_vec_init(n * nb);

    nb_steps = pattern[g - 1];
    split = g - 1;
    while ((split > 0) && (pattern[split - 1] <= nb_steps))
    {
        split --;
    }

    if (nb_steps > 0 || shifted_prec)
    {
        acb_theta_eld_distances(distances, zs, nb, tau, lp);
    }

    if (nb_steps == 0 && split == 0)
    {
        acb_theta_ql_exact_sum(th, zs, nb, tau, distances, all, shifted_prec, prec);
    }
    else if (nb_steps == 0 && split > 0)
    {
        acb_theta_ql_exact_lower_dim(th, zs, nb, tau, distances, split,
            pattern, all, shifted_prec, prec);
    }
    else
    {
        acb_theta_ql_exact_steps(th, zs, nb, tau, distances, split, pattern, all, prec);
    }

    _arb_vec_clear(distances, n * nb);
}
