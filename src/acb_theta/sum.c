/*
    Copyright (C) 2025 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "profiler.h"
#include "fmpz.h"
#include "acb.h"
#include "arb_mat.h"
#include "acb_mat.h"
#include "acb_modular.h"
#include "acb_theta.h"

static void
acb_theta_sum_00_worker(acb_ptr th, acb_srcptr v1, acb_srcptr v2,
    const slong * precs, slong len, const acb_t cofactor, const slong * coords,
    slong ord, slong g, slong prec, slong fullprec)
{
    acb_t sum;

    acb_init(sum);

    acb_dot(sum, NULL, 0, v1, 1, v2, 1, len, prec);
    acb_addmul(th, sum, cofactor, fullprec);

    acb_clear(sum);
}

static void
acb_theta_sum_0b_worker(acb_ptr th, acb_srcptr v1, acb_srcptr v2,
    const slong * precs, slong len, const acb_t cofactor, const slong * coords,
    slong ord, slong g, slong prec, slong fullprec)
{
    slong n = 1 << g;
    acb_t s0, s1, add, sub;
    ulong b;
    slong dot;

    acb_init(s0);
    acb_init(s1);
    acb_init(add);
    acb_init(sub);

    /* Compute alternate sums to adjust signs */
    acb_dot(s0, NULL, 0, v1, 2, v2, 2, (len + 1) / 2, prec);
    acb_dot(s1, NULL, 0, v1 + 1, 2, v2 + 1, 2, len / 2, prec);
    acb_add(add, s0, s1, prec);
    acb_sub(sub, s0, s1, prec);
    acb_mul(add, add, cofactor, prec);
    acb_mul(sub, sub, cofactor, prec);

    for (b = 0; b < n; b++)
    {
        dot = acb_theta_char_dot_slong(b, coords, g) % 2;
        if ((b >> (g - 1)) && dot)
        {
            acb_sub(&th[b], &th[b], sub, fullprec);
        }
        else if ((b >> (g - 1)))
        {
            acb_add(&th[b], &th[b], sub, fullprec);
        }
        else if (dot)
        {
            acb_sub(&th[b], &th[b], add, fullprec);
        }
        else
        {
            acb_add(&th[b], &th[b], add, fullprec);
        }
    }

    acb_clear(s0);
    acb_clear(s1);
    acb_clear(add);
    acb_clear(sub);
}

void
acb_theta_sum_0x(acb_ptr th, const acb_theta_ctx_z_struct * vec, slong nb,
    const acb_theta_ctx_tau_t ctx_tau, arb_srcptr distance, int all_b, slong prec)
{
    /* Call sum_work at the right precision */
    slong g = ctx_tau->g;
    slong n = (all_b ? 1 << g : 1);
    slong guard = ACB_THETA_LOW_PREC;
    acb_theta_eld_t E;
    arf_t R2, eps;
    arb_t err;
    arb_ptr v;
    int b;
    slong j, k;

    acb_theta_eld_init(E, g, g);
    arf_init(R2);
    arf_init(eps);
    arb_init(err);
    v = _arb_vec_init(g);

    acb_theta_ctx_z_common_v(v, vec, nb, prec + guard);
    acb_theta_sum_radius(R2, eps, &ctx_tau->cho, 0,
        prec + FLINT_MAX(0, acb_theta_sum_addprec(distance)));
    b = acb_theta_eld_set(E, &ctx_tau->cho, R2, v);

    /* flint_printf("(sum) ellipsoid radius: ");
       arf_printd(R2, 5);
       flint_printf("\n(sum) number of points: %wd\n", E->nb_pts); */

    if (b)
    {
        for (j = 0; j < nb; j++)
        {
            acb_theta_sum_work(th + j * n, n, (&vec[j])->exp_2z,
                (&vec[j])->exp_2z_inv, 1, ctx_tau->exp_tau,
                ctx_tau->exp_tau_inv, E, 0, prec,
                (all_b ? acb_theta_sum_0b_worker : acb_theta_sum_00_worker));
            arb_mul_arf(err, &(&vec[j])->u, eps, guard);
            /* flint_printf("(sum) error: ");
               arb_printd(err, 5);
               flint_printf("\n"); */
            for (k = 0; k < n; k++)
            {
                acb_add_error_arb(&th[j * n + k], err);
            }
            /* flint_printf("(sum) 1st value: ");
               acb_printd(&th[j * n], 5);
               flint_printf("\n"); */
        }
    }
    else
    {
        for (j = 0; j < nb * n; j++)
        {
            acb_indeterminate(&th[j]);
        }
    }

    acb_theta_eld_clear(E);
    arf_clear(R2);
    arf_clear(eps);
    arb_clear(err);
    _arb_vec_clear(v, g);
}

void
acb_theta_sum(acb_ptr th, const acb_theta_ctx_z_struct * vec, slong nb,
    const acb_theta_ctx_tau_t ctx_tau, arb_srcptr distances, int all_a,
    int all_b, int tilde, slong prec)
{
    slong g = ctx_tau->g;
    slong n = 1 << g;
    slong nba = (all_a ? n : 1);
    slong nbb = (all_b ? n : 1);
    slong nbab = nba * nbb;
    slong guard = ACB_THETA_LOW_PREC;
    acb_theta_ctx_z_struct * new_vec;
    acb_ptr cs, res;
    slong new_prec, dot;
    slong j, a, b;

    FLINT_ASSERT(nb >= 0);
    if (nb == 0)
    {
        return;
    }
    if (g == 1)
    {
        /* Find out if we want to use acb_modular_theta or acb_theta_sum_0x */
        new_prec = FLINT_MAX(prec, prec + acb_theta_sum_addprec(&distances[0]));
        if (all_a)
        {
            new_prec = FLINT_MAX(new_prec, prec + acb_theta_sum_addprec(&distances[1]));
        }
    }

    if (g == 1 && new_prec <= 4 * prec)
    {
        /* Call acb_modular_theta_sum directly: we accept to run computations
           at a slightly higher precision than necessary. */
        res = _acb_vec_init(4);
        for (j = 0; j < nb; j++)
        {
            acb_modular_theta_sum(&res[0], &res[1], &res[2], &res[3],
                (&vec[j])->exp_z, (&vec[j])->is_real,
                acb_mat_entry(ctx_tau->exp_tau, 0, 0), 1, new_prec);

            acb_set(&th[nbab * j], &res[2]);
            if (all_a && all_b)
            {
                acb_set(&th[nbab * j + 1], &res[3]);
                acb_set(&th[nbab * j + 2], &res[1]);
                acb_neg(&th[nbab * j + 3], &res[0]);
                _acb_vec_scalar_mul(th + nbab * j + 2, th + nbab * j + 2, 2,
                    acb_mat_entry(ctx_tau->exp_tau_div_4, 0, 0), prec + guard);
            }
            else if (all_a)
            {
                acb_set(&th[nbab * j + 1], &res[1]);
                acb_mul(&th[nbab * j + 1], &th[nbab * j + 1],
                    acb_mat_entry(ctx_tau->exp_tau_div_4, 0, 0), prec + guard);
            }
            else if (all_b)
            {
                acb_set(&th[nbab * j + 1], &res[3]);
            }
        }
        _acb_vec_clear(res, 4);
    }
    else if (all_a)
    {
        new_vec = acb_theta_ctx_z_vec_init(nb, g);
        res = _acb_vec_init(nbb * nb);
        cs = _acb_vec_init(nb);

        for (a = 0; a < n; a++)
        {
            for (j = 0; j < nb; j++)
            {
                acb_theta_ctx_z_shift_a0(&new_vec[j], &cs[j], &vec[j], ctx_tau, a, prec + guard);
            }
            acb_theta_sum_0x(res, new_vec, nb, ctx_tau, &distances[a], all_b, prec);
            for (j = 0; j < nb; j++)
            {
                _acb_vec_scalar_mul(th + nbab * j + nbb * a, res + nbb * j,
                    nbb, &cs[j], prec + guard);
            }
            for (b = 1; b < nbb; b++) /* No sign changes for b=0 */
            {
                dot = acb_theta_char_dot(a, b, g);
                for (j = 0; j < nb; j++)
                {
                    acb_mul_i_pow_si(&th[nbab * j + nbb * a + b], &th[nbab * j + nbb * a + b], dot);
                }
            }
        }

        acb_theta_ctx_z_vec_clear(new_vec, nb);
        _acb_vec_clear(res, nbb * nb);
        _acb_vec_clear(cs, nb);
    }
    else
    {
        acb_theta_sum_0x(th, vec, nb, ctx_tau, &distances[0], all_b, prec);
    }

    if (tilde)
    {
        for (j = 0; j < nb; j++)
        {
            _acb_vec_scalar_mul_arb(th + nbab * j, th + nbab * j, nbab,
                &(&vec[j])->uinv, prec + guard);
        }
    }
}
