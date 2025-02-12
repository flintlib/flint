/*
    Copyright (C) 2025 Jean Kieffer

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

static void
acb_theta_sum_jet_0x_worker(acb_ptr th, acb_srcptr v1, acb_srcptr v2,
    const slong * precs, slong len, const acb_t cofactor, const slong * coords,
    slong ord, int all_b, slong g, slong prec, slong fullprec)
{
    slong nbth = (all_b ? (1 << g) : 1);
    slong nbjet = acb_theta_jet_nb(ord, g);
    slong * tups;
    slong * dots;
    acb_ptr v3, aux;
    acb_t x, y;
    fmpz_t num, t;
    slong j, i, b;

    tups = flint_malloc(g * nbjet * sizeof(slong));
    dots = flint_malloc(nbth * sizeof(slong));
    v3 = _acb_vec_init(len);
    aux = _acb_vec_init(nbth * nbjet);
    acb_init(x);
    acb_init(y);
    fmpz_init(num);
    fmpz_init(t);

    for (b = 0; b < nbth; b++)
    {
        dots[b] = acb_theta_char_dot_slong(b, coords, g);
    }

    /* Compute products in v3 */
    for (i = 0; i < len; i++)
    {
        acb_mul(&v3[i], &v1[i], &v2[i], precs[i]);
    }

    acb_theta_jet_tuples(tups, ord, g);
    for (j = 0; j < nbjet; j++)
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
            for (b = 0; b < nbth; b++)
            {
                acb_mul_i_pow_si(y, x, 2 * ((dots[b] + i * acb_theta_char_bit(b, 0, g)) % 4));
                acb_add(&aux[b * nbjet + j], &aux[b * nbjet + j], y, prec);
            }
        }

        /* Multiply by cofactor * num */
        acb_mul_fmpz(x, cofactor, num, prec);
        for (b = 0; b < nbth; b++)
        {
            acb_mul(&aux[b * nbjet + j], &aux[b * nbjet + j], x, prec);
        }
    }
    _acb_vec_add(th, th, aux, nbth * nbjet, fullprec);

    flint_free(tups);
    flint_free(dots);
    _acb_vec_clear(v3, len);
    _acb_vec_clear(aux, nbth * nbjet);
    acb_clear(x);
    acb_clear(y);
    fmpz_clear(num);
    fmpz_clear(t);

}

static void
acb_theta_sum_jet_00_worker(acb_ptr th, acb_srcptr v1, acb_srcptr v2,
    const slong * precs, slong len, const acb_t cofactor, const slong * coords,
    slong ord, slong g, slong prec, slong fullprec)
{
    acb_theta_sum_jet_0x_worker(th, v1, v2, precs, len, cofactor, coords, ord,
        0, g, prec, fullprec);
}

static void
acb_theta_sum_jet_0b_worker(acb_ptr th, acb_srcptr v1, acb_srcptr v2,
    const slong * precs, slong len, const acb_t cofactor, const slong * coords,
    slong ord, slong g, slong prec, slong fullprec)
{
    acb_theta_sum_jet_0x_worker(th, v1, v2, precs, len, cofactor, coords, ord,
        1, g, prec, fullprec);
}

/* To compute derivatives of all theta values, we use a big ellipsoid to avoid
   complicated formulas for composing derivatives; this introduces powers of i
   in the worker */

static ulong
acb_theta_char_get_a(const slong * n, slong g)
{
    slong k;
    ulong a = 0;

    for (k = 0; k < g; k++)
    {
        a *= 2;
        a += n[k] & 1;
    }

    return a;
}

static void
acb_theta_sum_jet_ax_worker(acb_ptr th, acb_srcptr v1, acb_srcptr v2,
    const slong * precs, slong len, const acb_t cofactor, const slong * coords,
    slong ord, int all_b, slong g, slong prec, slong fullprec)
{
    slong n = 1 << g;
    slong nbb = (all_b ? n : 1);
    slong nbjet = acb_theta_jet_nb(ord, g);
    slong * tups;
    slong a0, a1;
    slong * dots;
    acb_ptr v3, aux;
    acb_t x, y;
    fmpz_t num, t;
    slong j, i;
    ulong b;

    tups = flint_malloc(g * nbjet * sizeof(slong));
    dots = flint_malloc(nbb * sizeof(slong));
    v3 = _acb_vec_init(len);
    aux = _acb_vec_init(nbjet * nbb * n);
    acb_init(x);
    acb_init(y);
    fmpz_init(num);
    fmpz_init(t);

    /* Precompute a0, a1, dots */
    a0 = acb_theta_char_get_a(coords, g);
    a1 = a0 ^ (1 << (g - 1));
    for (b = 0; b < nbb; b++)
    {
        dots[b] = acb_theta_char_dot_slong(b, coords, g);
    }

    /* Compute products in v3 */
    for (i = 0; i < len; i++)
    {
        acb_mul(&v3[i], &v1[i], &v2[i], precs[i]);
    }

    acb_theta_jet_tuples(tups, ord, g);
    for (j = 0; j < nbjet; j++)
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
            for (b = 0; b < nbb; b++)
            {
                acb_mul_i_pow_si(y, x, (dots[b] + i * acb_theta_char_bit(b, 0, g)) % 4);
                if (i % 2 == 0)
                {
                    acb_add(&aux[(nbb * a0 + b) * nbjet + j],
                        &aux[(nbb * a0 + b) * nbjet + j], y, prec);
                }
                else
                {
                    acb_add(&aux[(nbb * a1 + b) * nbjet + j],
                        &aux[(nbb * a1 + b) * nbjet + j], y, prec);
                }
            }
        }

        /* Multiply by cofactor * num */
        acb_mul_fmpz(x, cofactor, num, prec);
        for (b = 0; b < nbb; b++)
        {
            acb_mul(&aux[(nbb * a0 + b) * nbjet + j],
                &aux[(nbb * a0 + b) * nbjet + j], x, prec);
            acb_mul(&aux[(nbb * a1 + b) * nbjet + j],
                &aux[(nbb * a1 + b) * nbjet + j], x, prec);
        }
    }

    _acb_vec_add(th, th, aux, nbjet * nbb * n, fullprec);

    flint_free(tups);
    flint_free(dots);
    _acb_vec_clear(v3, len);
    _acb_vec_clear(aux, nbjet * nbb * n);
    acb_clear(x);
    acb_clear(y);
    fmpz_clear(num);
    fmpz_clear(t);
}


static void
acb_theta_sum_jet_a0_worker(acb_ptr th, acb_srcptr v1, acb_srcptr v2,
    const slong * precs, slong len, const acb_t cofactor, const slong * coords,
    slong ord, slong g, slong prec, slong fullprec)
{
    acb_theta_sum_jet_ax_worker(th, v1, v2, precs, len, cofactor, coords, ord,
        0, g, prec, fullprec);
}

static void
acb_theta_sum_jet_all_worker(acb_ptr th, acb_srcptr v1, acb_srcptr v2,
    const slong * precs, slong len, const acb_t cofactor, const slong * coords,
    slong ord, slong g, slong prec, slong fullprec)
{
    acb_theta_sum_jet_ax_worker(th, v1, v2, precs, len, cofactor, coords, ord,
        1, g, prec, fullprec);
}

void
acb_theta_sum_jet(acb_ptr th, const acb_theta_ctx_z_struct * vec, slong nb,
    const acb_theta_ctx_tau_t ctx_tau, slong ord, int all_a, int all_b, slong prec)
{
    slong g = ctx_tau->g;
    slong n = 1 << g;
    slong nbjet = acb_theta_jet_nb(ord, g);
    slong nbth = (all_a ? n : 1) * (all_b ? n : 1);
    slong guard = ACB_THETA_LOW_PREC;
    slong j, k;

    FLINT_ASSERT(nb >= 0);
    if (nb == 0)
    {
        return;
    }

    if (g == 1)
    {
        acb_ptr res;

        res = _acb_vec_init(4 * nbjet);
        for (j = 0; j < nb; j++)
        {
            /* acb_modular_theta_sum recomputes the inverse of exp_z */
            acb_modular_theta_sum(res, res + nbjet, res + 2 * nbjet, res + 3 * nbjet,
                (&vec[j])->exp_z, (&vec[j])->is_real,
                acb_mat_entry(ctx_tau->exp_tau, 0, 0), ord + 1, prec);
            _acb_vec_set(th + j * nbth * nbjet, res + 2 * nbjet, nbjet);
            if (all_a && all_b)
            {
                _acb_vec_set(th + (4 * j + 1) * nbjet, res + 3 * nbjet, nbjet);
                _acb_vec_set(th + (4 * j + 2) * nbjet, res + nbjet, nbjet);
                _acb_vec_neg(th + (4 * j + 3) * nbjet, res, nbjet);
                _acb_vec_scalar_mul(th + (4 * j + 2) * nbjet, th + (4 * j + 2) * nbjet, 2 * nbjet,
                    acb_mat_entry(ctx_tau->exp_tau_div_4, 0, 0), prec);
            }
            else if (all_a)
            {
                _acb_vec_set(th + (2 * j + 1) * nbjet, res + nbjet, nbjet);
                _acb_vec_scalar_mul(th + (2 * j + 1) * nbjet, th + (2 * j + 1) * nbjet, nbjet,
                    acb_mat_entry(ctx_tau->exp_tau_div_4, 0, 0), prec);
            }
            else if (all_b)
            {
                _acb_vec_set(th + (2 * j + 1) * nbjet, res + 3 * nbjet, nbjet);
            }
        }
        _acb_vec_clear(res, 4 * nbjet);
    }
    else
    {
        acb_theta_eld_t E;
        acb_ptr * sqr_pow;
        arf_t R2, eps;
        arb_t err;
        arb_ptr v;
        slong * tups;
        fmpz_t m, t;
        acb_t c;
        int b;

        acb_theta_eld_init(E, g, g);
        sqr_pow = flint_malloc(g * sizeof(acb_ptr));
        arf_init(R2);
        arf_init(eps);
        arb_init(err);
        v = _arb_vec_init(g);
        tups = flint_malloc(g * nbjet * sizeof(slong));
        acb_init(c);
        fmpz_init(m);
        fmpz_init(t);

        acb_theta_ctx_z_common_v(v, vec, nb, prec + guard);
        acb_theta_sum_jet_radius(R2, eps, &ctx_tau->cho, v, ord, prec);
        if (all_a)
        {
            /* Take into account that everything is duplicated in worker */
            _arb_vec_scalar_mul_2exp_si(v, v, g, 1);
            arf_mul_2exp_si(R2, R2, 2);
        }
        b = acb_theta_eld_set(E, &ctx_tau->cho, R2, v);

        if (b)
        {
            for (j = 0; j < g; j++)
            {
                sqr_pow[j] = _acb_vec_init(acb_theta_eld_box(E, j) + 1);
            }
            /* If all_a, use exp_z, exp_tau_div_4 instead of exp_2z, exp_tau */
            acb_theta_sum_sqr_pow(sqr_pow,
                (all_a ? ctx_tau->exp_tau_div_4 : ctx_tau->exp_tau),
                E, prec + guard);

            /* Sum series, rescale by c factor */
            for (j = 0; j < nb; j++)
            {
                if (all_a)
                {
                    acb_theta_sum_work(th + j * nbth * nbjet, nbth * nbjet,
                        (&vec[j])->exp_z, (&vec[j])->exp_z_inv,
                        ctx_tau->exp_tau_div_4, ctx_tau->exp_tau_div_4_inv,
                        sqr_pow, E, ord, prec,
                        (all_b ? acb_theta_sum_jet_all_worker : acb_theta_sum_jet_a0_worker));
                }
                else
                {
                    acb_theta_sum_work(th + j * nbth * nbjet, nbth * nbjet,
                        (&vec[j])->exp_2z,(&vec[j])->exp_2z_inv,
                        ctx_tau->exp_tau, ctx_tau->exp_tau_inv,
                        sqr_pow, E, ord, prec,
                        (all_b ? acb_theta_sum_jet_0b_worker : acb_theta_sum_jet_00_worker));
                }
                arb_mul_arf(err, &(&vec[j])->u, eps, guard);
                for (k = 0; k < nbth * nbjet; k++)
                {
                    acb_add_error_arb(&th[j * nbth * nbjet + k], err);
                }
            }

            /* Rescale by factorials and powers of 2pi*i */
            acb_theta_jet_tuples(tups, ord, g);
            for (k = 0; k < nbjet; k++)
            {
                acb_const_pi(c, prec);
                if (!all_a)
                {
                    acb_mul_2exp_si(c, c, 1);
                }
                acb_mul_onei(c, c);
                acb_pow_ui(c, c, acb_theta_jet_total_order(tups + k * g, g), prec);
                fmpz_one(m);
                for (j = 0; j < g; j++)
                {
                    fmpz_fac_ui(t, tups[k * g + j]);
                    fmpz_mul(m, m, t);
                }
                acb_div_fmpz(c, c, m, prec);
                for (j = 0; j < nb * nbth; j++)
                {
                    acb_mul(&th[j * nbjet + k], &th[j * nbjet + k], c, prec);
                }
            }

            for (j = 0; j < g; j++)
            {
                _acb_vec_clear(sqr_pow[j], acb_theta_eld_box(E, j) + 1);
            }
        }
        else
        {
            /* Should not happen in tests */
            _acb_vec_indeterminate(th, nb * nbth * nbjet);
        }

        acb_theta_eld_clear(E);
        flint_free(sqr_pow);
        arf_clear(R2);
        arf_clear(eps);
        arb_clear(err);
        _arb_vec_clear(v, g);
        flint_free(tups);
        acb_clear(c);
        fmpz_clear(m);
        fmpz_clear(t);
    }
}
