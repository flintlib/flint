/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

/* Use a big ellipsoid to avoid complicated formulas for derivatives; this
   introduces powers of i in worker_dim1 */

static void
worker_dim0(acb_ptr dth, slong len, const acb_t term, slong* coords, slong g,
    slong ord, slong prec, slong fullprec)
{
    slong n = 1 << g;
    slong nb_max = acb_theta_jet_nb(ord, g);
    slong nb_tot = acb_theta_jet_nb(ord, g + 1);
    acb_t x;
    fmpz_t m;
    acb_ptr f;
    ulong a, b;
    slong k, j, i, nb, ind;
    slong* orders;

    acb_init(x);
    fmpz_init(m);
    orders = flint_malloc(g * nb_max * sizeof(slong));
    f = _acb_vec_init(nb_max);

    a = acb_theta_char_get_a(coords, g);
    ind = 0;
    for (k = 0; k <= ord; k++)
    {
        /* Get list of orders */
        nb = acb_theta_jet_nb(k, g);
        acb_theta_jet_orders(orders, k, g);

        /* Compute factor for each tuple */
        for (j = 0; j < nb; j++)
        {
            acb_one(&f[j]);
            for (i = 0; i < g; i++)
            {
                fmpz_set_si(m, coords[i]);
                fmpz_pow_ui(m, m, orders[j * g + i]);
                acb_mul_fmpz(&f[j], &f[j], m, prec);
                fmpz_fac_ui(m, orders[j * g + i]);
                acb_div_fmpz(&f[j], &f[j], m, prec);
            }
        }
        acb_const_pi(x, prec);
        acb_mul_onei(x, x);
        acb_pow_ui(x, x, k, prec);
        _acb_vec_scalar_mul(f, f, nb, x, prec);

        /* Loop over b */
        for (b = 0; b < n; b++)
        {
            acb_mul_powi(x, term, acb_theta_char_dot_slong(b, coords, g) % 4);
            for (j = 0; j < nb; j++)
            {
                acb_addmul(&dth[(n * a + b) * nb_tot + ind + j], x, &f[j], fullprec);
            }
        }

        ind += nb;
    }

    acb_clear(x);
    fmpz_clear(m);
    flint_free(orders);
    _acb_vec_clear(f, nb_max);
}


static void
worker_dim1(acb_ptr dth, acb_srcptr v1, acb_srcptr v2, const slong* precs, slong len,
    const acb_t cofactor, const slong* coords, slong ord, slong g, slong prec, slong fullprec)
{
    slong n = 1 << g;
    slong nb_max = acb_theta_jet_nb(ord, g);
    slong nb_tot = acb_theta_jet_nb(ord, g + 1);
    slong* orders;
    slong a0, a1;
    slong* dots;
    acb_ptr v3, aux;
    acb_t x, y;
    fmpz_t num, den, t;
    slong k, j, i, nb, ind;
    ulong b;

    orders = flint_malloc(g * nb_max * sizeof(slong));
    dots = flint_malloc(n * sizeof(slong));
    v3 = _acb_vec_init(len);
    aux = _acb_vec_init(nb_tot * n * n);
    acb_init(x);
    acb_init(y);
    fmpz_init(num);
    fmpz_init(den);
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

    ind = 0;
    for (k = 0; k <= ord; k++)
    {
        /* Get list of orders */
        nb = acb_theta_jet_nb(k, g);
        acb_theta_jet_orders(orders, k, g);
        _acb_vec_zero(aux, nb_tot * n * n);

        /* Process each tuple of orders */
        for (j = 0; j < nb; j++)
        {
            /* Set cofactors as fmpz's */
            fmpz_one(num);
            fmpz_one(den);
            for (i = 1; i < g; i++)
            {
                fmpz_set_si(t, coords[i]);
                fmpz_pow_ui(t, t, orders[j * g + i]);
                fmpz_mul(num, num, t);
            }
            for (i = 0; i < g; i++)
            {
                fmpz_fac_ui(t, orders[j * g + i]);
                fmpz_mul(den, den, t);
            }

            /* Now loop over coordinates; compute right cofactor */
            for (i = 0; i < len; i++)
            {
                fmpz_set_si(t, coords[0] + i);
                fmpz_pow_ui(t, t, orders[j * g]);
                acb_mul_fmpz(x, &v3[i], t, precs[i]);
                /* Now loop over b, adding coefficients in both a0b and a1b */
                for (b = 0; b < n; b++)
                {
                    acb_mul_powi(y, x, (dots[b] + i * (b >> (g - 1))) % 4);
                    if (i % 2 == 0)
                    {
                        acb_add(&aux[(n * a0 + b) * nb_tot + ind + j],
                            &aux[(n * a0 + b) * nb_tot + ind + j], y, prec);
                    }
                    else
                    {
                        acb_add(&aux[(n * a1 + b) * nb_tot + ind + j],
                            &aux[(n * a1 + b) * nb_tot + ind + j], y, prec);
                    }
                }
            }

            /* Finally, loop over b to multiply by num/den */
            for (b = 0; b < n; b++)
            {
                acb_mul_fmpz(&aux[(n * a0 + b) * nb_tot + ind + j],
                    &aux[(n * a0 + b) * nb_tot + ind + j], num, prec);
                acb_mul_fmpz(&aux[(n * a1 + b) * nb_tot + ind + j],
                    &aux[(n * a1 + b) * nb_tot + ind + j], num, prec);
                acb_div_fmpz(&aux[(n * a0 + b) * nb_tot + ind + j],
                    &aux[(n * a0 + b) * nb_tot + ind + j], den, prec);
                acb_div_fmpz(&aux[(n * a1 + b) * nb_tot + ind + j],
                    &aux[(n * a1 + b) * nb_tot + ind + j], den, prec);
            }
        }

        /* Multiply the whole thing by cofactor * (i pi)^k and add to dth */
        acb_const_pi(x, prec);
        acb_mul_onei(x, x);
        acb_pow_ui(x, x, k, prec);
        acb_mul(x, x, cofactor, prec);
        _acb_vec_scalar_mul(aux, aux, nb_tot * n * n, x, prec);
        _acb_vec_add(dth, dth, aux, nb_tot * n * n, fullprec);
        ind += nb;
    }

    flint_free(orders);
    flint_free(dots);
    _acb_vec_clear(v3, len);
    _acb_vec_clear(aux, nb_tot * n * n);
    acb_clear(x);
    acb_clear(y);
    fmpz_clear(num);
    fmpz_clear(den);
    fmpz_clear(t);
}

static void
acb_theta_jet_naive_all_gen(acb_ptr dth, acb_srcptr z, const acb_mat_t tau,
    slong ord, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_theta_eld_t E;
    acb_theta_precomp_t D;
    acb_t c;
    arb_t u;
    acb_mat_t new_tau;
    acb_ptr new_z;
    slong nb = n * n * acb_theta_jet_nb(ord, g + 1);

    acb_theta_eld_init(E, g, g);
    acb_theta_precomp_init(D, 1, g);
    acb_init(c);
    arb_init(u);
    acb_mat_init(new_tau, g, g);
    new_z = _acb_vec_init(g);

    _acb_vec_scalar_mul_2exp_si(new_z, z, g, -1);
    acb_mat_scalar_mul_2exp_si(new_tau, tau, -2);

    acb_theta_jet_ellipsoid(E, u, new_z, new_tau, ord, prec);
    prec = acb_theta_naive_fullprec(E, prec);
    acb_theta_precomp_set(D, new_z, new_tau, E, prec);
    acb_one(c);

    acb_theta_naive_worker_new(dth, nb, c, u, E, D, 0, ord, prec, worker_dim1);

    acb_theta_eld_clear(E);
    acb_theta_precomp_clear(D);
    acb_clear(c);
    arb_clear(u);
    acb_mat_clear(new_tau);
    _acb_vec_clear(new_z, g);
}

void
acb_theta_jet_naive_all(acb_ptr dth, acb_srcptr z, const acb_mat_t tau,
    slong ord, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong nb = acb_theta_jet_nb(ord, g + 1);
    acb_ptr res;

    if (g == 1)
    {
        res = _acb_vec_init(4 * nb);

        acb_modular_theta_jet(res, res + nb, res + 2 * nb, res + 3 * nb,
            z, acb_mat_entry(tau, 0, 0), nb, prec);
        _acb_vec_set(dth, res + 2 * nb, nb);
        _acb_vec_set(dth + nb, res + 3 * nb, nb);
        _acb_vec_set(dth + 2 * nb, res + nb, nb);
        _acb_vec_neg(dth + 3 * nb, res, nb);

        _acb_vec_clear(res, 4 * nb);
    }
    else
    {
        acb_theta_jet_naive_all_gen(dth, z, tau, ord, prec);
    }
}
