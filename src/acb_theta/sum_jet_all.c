/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "acb.h"
#include "acb_mat.h"
#include "acb_modular.h"
#include "acb_theta.h"

/* Use a big ellipsoid to avoid complicated formulas for derivatives; this
   introduces powers of i in worker */
static ulong
acb_theta_char_get_a(const slong * n, slong g)
{
    slong k;
    ulong a = 0;

    for (k = 0; k < g; k++)
    {
        a *= 2;
        a += ((n[k] % 2) + 2) % 2;
    }

    return a;
}

static void
acb_theta_sum_jet_all_worker(acb_ptr th, acb_srcptr v1, acb_srcptr v2,
    const slong * precs, slong len, const acb_t cofactor, const slong * coords,
    slong ord, slong g, slong prec, slong fullprec)
{
    slong n = 1 << g;
    slong nb = acb_theta_jet_nb(ord, g);
    slong * tups;
    slong a0, a1;
    slong * dots;
    acb_ptr v3, aux;
    acb_t x, y;
    fmpz_t num, t;
    slong j, i;
    ulong b;

    tups = flint_malloc(g * nb * sizeof(slong));
    dots = flint_malloc(n * sizeof(slong));
    v3 = _acb_vec_init(len);
    aux = _acb_vec_init(nb * n * n);
    acb_init(x);
    acb_init(y);
    fmpz_init(num);
    fmpz_init(t);

    /* Precompute a0, a1, dots */
    a0 = acb_theta_char_get_a(coords, g);
    a1 = a0 ^ (1 << (g - 1));
    for (b = 0; b < n; b++)
    {
        dots[b] = acb_theta_char_dot_slong(b, coords, g);
    }

    /* Compute products in v3 */
    for (i = 0; i < len; i++)
    {
        acb_mul(&v3[i], &v1[i], &v2[i], precs[i]);
    }

    acb_theta_jet_tuples(tups, ord, g);
    for (j = 0; j < nb; j++)
    {
        fmpz_one(num);
        for (i = 1; i < g; i++)
        {
            fmpz_set_si(t, coords[i]);
            fmpz_pow_ui(t, t, tups[j * g + i]);
            fmpz_mul(num, num, t);
        }

        /* Loop over lattice points */
        for (i = 0; i < len; i++)
        {
            fmpz_set_si(t, coords[0] + i);
            fmpz_pow_ui(t, t, tups[j * g]);
            acb_mul_fmpz(x, &v3[i], t, precs[i]);
            /* Loop over b, adding coefficients in both a0b and a1b */
            for (b = 0; b < n; b++)
            {
                acb_mul_i_pow_si(y, x, (dots[b] + i * (b >> (g - 1))) % 4);
                if (i % 2 == 0)
                {
                    acb_add(&aux[(n * a0 + b) * nb + j],
                        &aux[(n * a0 + b) * nb + j], y, prec);
                }
                else
                {
                    acb_add(&aux[(n * a1 + b) * nb + j],
                        &aux[(n * a1 + b) * nb + j], y, prec);
                }
            }
        }

        /* Multiply by cofactor * num */
        acb_mul_fmpz(x, cofactor, num, prec);
        for (b = 0; b < n; b++)
        {
            acb_mul(&aux[(n * a0 + b) * nb + j], &aux[(n * a0 + b) * nb + j], x, prec);
            acb_mul(&aux[(n * a1 + b) * nb + j], &aux[(n * a1 + b) * nb + j], x, prec);
        }
    }

    _acb_vec_add(th, th, aux, nb * n * n, fullprec);

    flint_free(tups);
    flint_free(dots);
    _acb_vec_clear(v3, len);
    _acb_vec_clear(aux, nb * n * n);
    acb_clear(x);
    acb_clear(y);
    fmpz_clear(num);
    fmpz_clear(t);
}

void
acb_theta_sum_jet_all(acb_ptr th, const acb_theta_ctx_z_struct * vec, slong nb,
    const acb_theta_ctx_tau_t ctx_tau, slong ord, slong prec)
{
    slong g = ctx_tau->g;
    slong n2 = 1 << (2 * g);
    slong nbth = acb_theta_jet_nb(ord, g);
    slong j, k;

    FLINT_ASSERT(nb >= 0);
    if (nb == 0)
    {
        return;
    }

    if (g == 1)
    {
        acb_ptr res;

        res = _acb_vec_init(4 * nbth);
        for (j = 0; j < nb; j++)
        {
            /* acb_modular_theta_sum recomputes the inverse of exp_z */
            acb_modular_theta_sum(res, res + nbth, res + 2 * nbth, res + 3 * nbth,
                (&vec[j])->exp_z, (&vec[j])->is_real,
                acb_mat_entry(ctx_tau->exp_tau, 0, 0), ord + 1, prec);
            _acb_vec_set(th + 4 * j * nbth, res + 2 * nbth, nbth);
            _acb_vec_set(th + (4 * j + 1) * nbth, res + 3 * nbth, nbth);
            _acb_vec_set(th + (4 * j + 2) * nbth, res + nbth, nbth);
            _acb_vec_neg(th + (4 * j + 3) * nbth, res, nbth);
            _acb_vec_scalar_mul(th + (4 * j + 2) * nbth, th + (4 * j + 2) * nbth, 2 * nbth,
                acb_mat_entry(ctx_tau->exp_tau_div_4, 0, 0), prec);
        }
        _acb_vec_clear(res, 4 * nbth);
    }
    else
    {
        acb_theta_eld_t E;
        arf_t R2, eps;
        arb_ptr v;
        slong * tups;
        fmpz_t m, t;
        arb_t u;
        acb_t c;
        int b;

        acb_theta_eld_init(E, g, g);
        v = _arb_vec_init(g);
        arf_init(R2);
        arf_init(eps);
        tups = flint_malloc(g * nbth * sizeof(slong));
        fmpz_init(m);
        fmpz_init(t);
        arb_init(u);
        acb_init(c);

        /* Take into account that everything is duplicated in worker */
        acb_theta_ctx_z_common_v(v, vec, nb, prec);
        acb_theta_sum_jet_radius(R2, eps, &ctx_tau->cho, v, ord, prec);
        _arb_vec_scalar_mul_2exp_si(v, v, g, 1);
        arf_mul_2exp_si(R2, R2, 2);
        b = acb_theta_eld_set(E, &ctx_tau->cho, R2, v);

        if (b)
        {
            for (j = 0; j < nb; j++)
            {
                /* Duplication in worker: use exp_z instead of exp_2z,
                   exp_tau_div_4 instead of exp_tau */
                acb_theta_sum_work(th + j * n2 * nbth, n2 * nbth,
                    (&vec[j])->exp_z, (&vec[j])->exp_z_inv,
                    1, ctx_tau->exp_tau_div_4, ctx_tau->exp_tau_div_4_inv, E, ord,
                    prec, acb_theta_sum_jet_all_worker);
                arb_mul_arf(u, &(&vec[j])->u, eps, prec);
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
                for (j = 0; j < nb * n2; j++)
                {
                    acb_mul(&th[j * nbth + k], &th[j * nbth + k], c, prec);
                }
            }
        }
        else
        {
            for (k = 0; k < nb * n2 * nbth; k++)
            {
                acb_indeterminate(&th[k]);
            }
        }

        acb_theta_eld_clear(E);
        arf_clear(R2);
        arf_clear(eps);
        _arb_vec_clear(v, g);
        flint_free(tups);
        fmpz_clear(m);
        fmpz_clear(t);
        arb_clear(u);
        acb_clear(c);
    }
}
