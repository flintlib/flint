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
acb_theta_sum_jet_00(acb_ptr th, const acb_theta_ctx_t ctx, slong ord, slong prec)
{
    slong g = acb_theta_ctx_g(ctx);
    slong nbth = acb_theta_jet_nb(ord, g);
    slong nbz = acb_theta_ctx_nb(ctx);
    arb_ptr a;
    acb_ptr aux;
    slong j, k;

    if (nbz == 0)
    {
        return;
    }

    a = _arb_vec_init(g);
    aux = _acb_vec_init(nbth);

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
            _acb_vec_scalar_mul(th + j * nbth, res + 2 * nbth, nbth, &acb_theta_ctx_cs(ctx)[j], prec);
        }
        _acb_vec_clear(res, 4 * nbth);
    }
    else
    {
        acb_theta_eld_t E;
        arf_t R2, eps;
        arb_t err;
        arb_ptr v;
        slong * tups;
        fmpz_t m, t;
        acb_t c;
        int b;

        acb_theta_eld_init(E, g, g);
        arf_init(R2);
        arf_init(eps);
        arb_init(err);
        v = _arb_vec_init(g);
        tups = flint_malloc(g * nbth * sizeof(slong));
        acb_init(c);
        fmpz_init(m);
        fmpz_init(t);

        acb_theta_ctx_common_v(v, ctx, prec);
        acb_theta_jet_naive_radius(R2, eps, acb_theta_ctx_cho(ctx), v, ord, prec);
        b = acb_theta_eld_set(E, acb_theta_ctx_cho(ctx), R2, v);

        if (b)
        {
            /* Sum series, rescale by c factor */
            acb_theta_sum_work(th, nbth, acb_theta_ctx_exp_2zs(ctx), acb_theta_ctx_exp_2zs_inv(ctx), nbz,
                acb_theta_ctx_exp_tau(ctx), acb_theta_ctx_exp_tau_inv(ctx), E, ord,
                prec, acb_theta_sum_jet_00_worker);
            for (j = 0; j < nbz; j++)
            {
                _acb_vec_scalar_mul(th + j * nbth, th + j * nbth, nbth, &acb_theta_ctx_cs(ctx)[j], prec);
                arb_mul_arf(err, &acb_theta_ctx_us(ctx)[j], eps, prec);
                for (k = 0; k < nbth; k++)
                {
                    acb_add_error_arb(&th[j * nbth + k], err);
                }
            }

            /* Rescale by factorials and powers of 2pi*i */
            acb_theta_jet_tuples(tups, ord, g);
            for (k = 0; k < nbth; k++)
            {
                acb_const_pi(c, prec);
                acb_mul_2exp_si(c, c, 1);
                acb_mul_onei(c, c);
                acb_pow_ui(c, c, acb_theta_jet_total_order(tups + k * g, g), prec);
                fmpz_one(m);
                for (j = 0; j < g; j++)
                {
                    fmpz_fac_ui(t, tups[k * g + j]);
                    fmpz_mul(m, m, t);
                }
                acb_div_fmpz(c, c, m, prec);
                for (j = 0; j < nbz; j++)
                {
                    acb_mul(&th[j * nbth + k], &th[j * nbth + k], c, prec);
                }
            }
        }
        else
        {
            for (j = 0; j < nbz * nbth; j++)
            {
                acb_indeterminate(&th[j]);
            }
        }

        acb_theta_eld_clear(E);
        arf_clear(R2);
        arf_clear(eps);
        arb_clear(err);
        _arb_vec_clear(v, g);
        flint_free(tups);
        acb_clear(c);
        fmpz_clear(m);
        fmpz_clear(t);
    }

    /* Both for genus 1 and higher, adjust derivatives due to reduction */
    for (j = 0; j < nbz; j++)
    {
        _arb_vec_neg(a, acb_theta_ctx_as(ctx) + j * g, g);
        _arb_vec_scalar_mul_2exp_si(a, a, g, 1);
        acb_theta_jet_exp_pi_i(aux, a, ord, g, prec);
        acb_theta_jet_mul(th + j * nbth, th + j * nbth, aux, ord, g, prec);
    }

    _arb_vec_clear(a, g);
    _acb_vec_clear(aux, nbth);
}
