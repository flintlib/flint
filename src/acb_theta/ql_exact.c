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
acb_theta_ql_exact_sum(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau,
    arb_srcptr distances, int all, int shifted_prec, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
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
            if (all)
            {
                acb_theta_sum_all_tilde(th + j * n * n, &vec[j], 1, ctx_tau, distances + j * n, prec);
            }
            else
            {
                acb_theta_sum_a0_tilde(th + j * n, &vec[j], 1, ctx_tau, distances + j * n, prec);
            }
        }
    }
    else
    {
        /* distances are all set to zero; use one ellipsoid */
        if (all)
        {
            acb_theta_sum_all_tilde(th, vec, nb, ctx_tau, distances, prec);
        }
        else
        {
            acb_theta_sum_a0_tilde(th, vec, nb, ctx_tau, distances, prec);
        }
    }

    acb_theta_ctx_tau_clear(ctx_tau);
    acb_theta_ctx_z_vec_clear(vec, nb);
}

static int
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
        /* flint_printf("(ql_exact_lower_dim) calling ql_exact on %wd vectors in dimension %wd\n", nb0, s); */
        res = acb_theta_ql_exact(th0, z0s, nb0, tau0, pattern, all, shifted_prec, prec);
    }

    if (res)
    {
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
        /* flint_printf("WARNING: ql_lower_dim failed, falling back to summation\n");
           flint_printf("g = %wd, nb = %wd, tau:\n", g, nb);
           acb_mat_printd(tau, 5);
        flint_printf("zs:\n");
        for (j = 0; j < nb; j++)
        {
            _acb_vec_printd(zs + j * g, g, 5);
            } */

        acb_theta_ql_exact_sum(th, zs, nb, tau, distances, all, shifted_prec, prec);
        res = 1;
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
    return res;
}

static int
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
    slong guard, hp;
    slong j, k;
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
    hp = prec + nb_steps * guard;
    acb_mat_scalar_mul_2exp_si(new_tau, tau, nb_steps);

    flint_printf("(ql_exact_steps) g = %wd, split = %wd, setup: %wd, nb_steps = %wd, guard = %wd, prec = %wd, hp = %wd\n",
        g, split, res, nb_steps, guard, prec, hp);
    /* if (easy_steps[0] < nb_steps)
    {
        for (j = 0; j < nb; j++)
        {
            flint_printf("%wd -> %wd\n", j, easy_steps[j]);
        }
    }
    flint_printf("distances:\n");
    _arb_vec_printd(distances, nb * n, 5); */

    if (res && (split > 0))
    {
        /* Set list of zs for which we want theta to be computed. */
        /* We don't take advantage of the fact that some thetas differ by a real value */
        acb_ptr new_th, new_z;
        arb_ptr new_distances;
        slong new_nb, add;
        slong * new_pattern;

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
                _acb_vec_add(new_z + new_nb * g, zs + j * g, t, g, hp);
                _acb_vec_add(new_z + (new_nb + 1) * g, new_z + new_nb * g, t, g, hp);
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

        _acb_vec_scalar_mul_2exp_si(new_z, new_z, new_nb * g, nb_steps);
        _arb_vec_scalar_mul_2exp_si(new_distances, new_distances, new_nb * n, nb_steps);

        /* Recursive call */
        res = acb_theta_ql_exact_lower_dim(new_th, new_z, new_nb, new_tau,
            new_distances, split, new_pattern, 0, 1, hp);

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

        _acb_vec_clear(new_th, new_nb * n);
        _acb_vec_clear(new_z, new_nb * g);
        _arb_vec_clear(new_distances, new_nb * n);
        flint_free(new_pattern);
    }
    else if (res) /* split = 0; set up th_init efficiently using sum_a0_tilde */
    {
        acb_theta_ctx_tau_t ctx_tau;
        acb_theta_ctx_z_struct * aux;
        acb_theta_ctx_z_t ctxt;
        acb_ptr new_z;
        arb_ptr d;

        acb_theta_ctx_tau_init(ctx_tau, 1, g);
        aux = acb_theta_ctx_z_vec_init(3, g);
        new_z = _acb_vec_init(g);
        d = _arb_vec_init(n);
        if (easy_steps[0] < nb_steps)
        {
            acb_theta_ctx_z_init(ctxt, g);
        }

        /* flint_printf("(ql_exact_steps) guard = %wd, nb_steps = %wd, hp = %wd, new_tau, new_zs:\n", guard, nb_steps, hp);
           acb_mat_printd(new_tau, 5);*/

        acb_theta_ctx_tau_set(ctx_tau, new_tau, hp);
        if (easy_steps[0] < nb_steps)
        {
            _acb_vec_scalar_mul_2exp_si(new_z, t, g, nb_steps);
            acb_theta_ctx_z_set(ctxt, new_z, ctx_tau, hp);

            /* _acb_vec_printd(new_z, g, 5);*/
        }


        for (j = 0; j < nb; j++)
        {
            _acb_vec_scalar_mul_2exp_si(new_z, zs + j * g, g, nb_steps);
            acb_theta_ctx_z_set(&aux[0], new_z, ctx_tau, hp);
            _arb_vec_scalar_mul_2exp_si(d, distances + j * n, n, nb_steps);
            if (easy_steps[j] == nb_steps)
            {
                acb_theta_sum_a0_tilde(th_init + 3 * n * j, aux, 1, ctx_tau, d, hp);
            }
            else
            {
                acb_theta_ctx_z_add_real(&aux[1], &aux[0], ctxt, hp);
                acb_theta_ctx_z_add_real(&aux[2], &aux[1], ctxt, hp);
                if (j == 0)
                {
                    acb_theta_sum_a0_tilde(th_init + 3 * n * j, aux, 3, ctx_tau, d, hp);
                }
                else
                {
                    acb_theta_sum_a0_tilde(th_init + 3 * n * j + n, aux + 1, 2, ctx_tau, d, hp);
                }
            }
        }

        acb_theta_ctx_tau_clear(ctx_tau);
        acb_theta_ctx_z_vec_clear(aux, 3);
        _acb_vec_clear(new_z, g);
        _arb_vec_clear(d, n);
        if (easy_steps[0] < nb_steps)
        {
            acb_theta_ctx_z_clear(ctxt);
        }
    }

    if (res) /* th_init is set; now perform steps */
    {
        /* flint_printf("(ql_exact_steps) got th_init:\n");
           _acb_vec_printd(th_init, 3 * n * nb, 5);*/

        acb_theta_ql_steps(th, th_init, rts, rts_all, nb, nb_steps, distances,
            easy_steps, all, g, hp);
    }

    else /* setup did not succeed: fall back to summation */
    {
        /* flint_printf("WARNING: ql_setup failed, falling back to summation\n");
        flint_printf("g = %wd, nb = %wd, tau:\n", g, nb);
        acb_mat_printd(tau, 5);
        flint_printf("zs:\n");
        for (j = 0; j < nb; j++)
        {
            _acb_vec_printd(zs + j * g, g, 5);
            } */

        acb_theta_ql_exact_sum(th, zs, nb, tau, distances, all, 1, prec);
        res = 1;
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
    return res;
}

int acb_theta_ql_exact(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau,
    const slong * pattern, int all, int shifted_prec, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    arb_ptr distances;
    slong nb_steps, split;
    slong lp = ACB_THETA_LOW_PREC;
    int res = 1;

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
        acb_theta_agm_distances(distances, zs, nb, tau, lp);
    }

    if (nb_steps == 0 && split == 0)
    {
        acb_theta_ql_exact_sum(th, zs, nb, tau, distances, all, shifted_prec, prec);
    }
    else if (nb_steps == 0 && split > 0)
    {
        res = acb_theta_ql_exact_lower_dim(th, zs, nb, tau, distances, split,
            pattern, all, shifted_prec, prec);
    }
    else
    {
        res = acb_theta_ql_exact_steps(th, zs, nb, tau, distances, split, pattern, all, prec);
    }

    _arb_vec_clear(distances, n * nb);
    return res;
}
