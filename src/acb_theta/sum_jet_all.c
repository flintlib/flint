/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_mat.h"
#include "acb_modular.h"
#include "acb_theta.h"

void
acb_theta_sum_jet_all(acb_ptr th, const acb_theta_ctx_t ctx, slong ord, slong prec)
{
    slong g = acb_theta_ctx_g(ctx);
    slong n2 = 1 << (2 * g);
    slong nbth = acb_theta_jet_nb(ord, g);
    slong nbz = acb_theta_ctx_nb(ctx);
    arb_ptr a;
    acb_ptr aux;
    arb_t u;
    acb_t c;
    slong j, k, l;

    if (nbz == 0)
    {
        return;
    }

    a = _arb_vec_init(g);
    aux = _acb_vec_init(nbth);
    arb_init(u);
    acb_init(c);

    if (g == 1)
    {
        acb_ptr res;

        res = _acb_vec_init(4 * nbth);
        for (j = 0; j < nbz; j++)
        {
            /* acb_modular_theta_sum recomputes the inverse of exp_z */
            /* todo: store w_is_unit as part of context */
            acb_modular_theta_sum(res, res + nbth, res + 2 * nbth, res + 3 * nbth,
                &acb_theta_ctx_exp_zs(ctx)[j], 0,
                acb_mat_entry(acb_theta_ctx_exp_tau(ctx), 0, 0), ord + 1, prec);
            _acb_vec_set(th + 4 * j * nbth, res + 2 * nbth, nbth);
            _acb_vec_set(th + (4 * j + 1) * nbth, res + 3 * nbth, nbth);
            _acb_vec_set(th + (4 * j + 2) * nbth, res + nbth, nbth);
            _acb_vec_neg(th + (4 * j + 3) * nbth, res, nbth);
            _acb_vec_scalar_mul(th + (4 * j + 2) * nbth, th + (4 * j + 2) * nbth, 2 * nbth,
                acb_mat_entry(acb_theta_ctx_exp_tau_div_4(ctx), 0, 0), prec);
            _acb_vec_scalar_mul(th + 4 * j * nbth, th + 4 * j * nbth, 4 * nbth,
                &acb_theta_ctx_cs(ctx)[j], prec);
        }
        _acb_vec_clear(res, 4 * nbth);
    }
    else
    {
        acb_theta_eld_t E;
        arf_t R2, eps;
        arb_ptr v;
        acb_mat_t exp_tau_div_4_inv;
        slong * tups;
        fmpz_t m, t;
        int b;

        acb_theta_eld_init(E, g, g);
        v = _arb_vec_init(g);
        acb_mat_init(exp_tau_div_4_inv, g, g);
        arf_init(R2);
        arf_init(eps);
        tups = flint_malloc(g * nbth * sizeof(slong));
        fmpz_init(m);
        fmpz_init(t);

        /* Take into account that everything is duplicated in worker */
        acb_theta_ctx_common_v(v, ctx, prec);
        acb_theta_jet_naive_radius(R2, eps, acb_theta_ctx_cho(ctx), v, ord, prec);
        _arb_vec_scalar_mul_2exp_si(v, v, g, 1);
        arf_mul_2exp_si(R2, R2, 2);
        b = acb_theta_eld_set(E, acb_theta_ctx_cho(ctx), R2, v);

        if (b)
        {
            for (k = 0; k < g; k++)
            {
                for (j = k + 1; j < g; j++)
                {
                    acb_inv(acb_mat_entry(exp_tau_div_4_inv, k, j),
                        acb_mat_entry(acb_theta_ctx_exp_tau_div_4(ctx), k, j), prec);
                }
            }

            /* Duplication in worker: use exp_z instead of exp_2z, exp_tau_div_4 instead of exp_tau */
            acb_theta_sum_work(th, n2 * nbth, acb_theta_ctx_exp_zs(ctx), acb_theta_ctx_exp_zs_inv(ctx),
                nbz, acb_theta_ctx_exp_tau_div_4(ctx), exp_tau_div_4_inv, E, ord,
                prec, acb_theta_sum_jet_all_worker);
            for (j = 0; j < nbz; j++)
            {
                _acb_vec_scalar_mul(th + j * n2 * nbth, th + j * n2 * nbth, n2 * nbth,
                    &acb_theta_ctx_cs(ctx)[j], prec);
                arb_mul_arf(u, &acb_theta_ctx_us(ctx)[j], eps, prec);
                for (k = 0; k < n2 * nbth; k++)
                {
                    acb_add_error_arb(&th[j * n2 * nbth + k], u);
                }
            }

            acb_theta_jet_tuples(tups, ord, g);
            for (k = 0; k < nbth; k++)
            {
                acb_const_pi(c, prec); /* not 2 pi because of rescaling */
                acb_mul_onei(c, c);
                acb_pow_ui(c, c, acb_theta_jet_total_order(tups + k * g, g), prec);
                fmpz_one(m);
                for (j = 0; j < g; j++)
                {
                    fmpz_fac_ui(t, tups[k * g + j]);
                    fmpz_mul(m, m, t);
                }
                acb_div_fmpz(c, c, m, prec);
                for (j = 0; j < nbz * n2; j++)
                {
                    acb_mul(&th[j * nbth + k], &th[j * nbth + k], c, prec);
                }
            }
        }
        else
        {
            for (k = 0; k < nbz * n2 * nbth; k++)
            {
                acb_indeterminate(&th[k]);
            }
        }

        acb_theta_eld_clear(E);
        arf_clear(R2);
        arf_clear(eps);
        _arb_vec_clear(v, g);
        acb_mat_clear(exp_tau_div_4_inv);
        flint_free(tups);
        fmpz_clear(m);
        fmpz_clear(t);
    }

    /* Both for genus 1 and higher, adjust derivatives after reduction */
    for (j = 0; j < nbz; j++)
    {
        _arb_vec_neg(a, acb_theta_ctx_as(ctx) + j * g, g);
        _arb_vec_scalar_mul_2exp_si(a, a, g, 1);
        acb_theta_jet_exp_pi_i(aux, a, ord, g, prec);
        for (k = 0; k < n2; k++)
        {
            acb_theta_jet_mul(th + j * n2 * nbth + k * nbth, th + j * n2 * nbth + k * nbth,
                aux, ord, g, prec);
            arb_zero(u);
            for (l = 0; l < g; l++)
            {
                if ((k >> (g - 1 - l)) % 2 == 1)
                {
                    arb_add(u, u, &a[l], prec);
                }
            }
            acb_onei(c);
            acb_pow_arb(c, c, u, prec);
            _acb_vec_scalar_mul(th + j * n2 * nbth + k * nbth,
                th + j * n2 * nbth + k * nbth, nbth, c, prec);
        }
    }

    _acb_vec_clear(aux, nbth);
    _arb_vec_clear(a, g);
    arb_clear(u);
    acb_clear(c);
}
