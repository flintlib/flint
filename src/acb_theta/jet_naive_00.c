/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

static void
worker_dim1(acb_ptr dth, acb_srcptr v1, acb_srcptr v2, const slong* precs, slong len,
    const acb_t cofactor, const slong* coords, slong ord, slong g, slong prec, slong fullprec)
{
    slong nb = acb_theta_jet_nb(ord, g);
    slong* orders;
    acb_ptr v3, aux;
    acb_t x;
    fmpz_t num, den, t;
    slong j, i;

    orders = flint_malloc(g * nb * sizeof(slong));
    v3 = _acb_vec_init(len);
    aux = _acb_vec_init(nb);
    acb_init(x);
    fmpz_init(num);
    fmpz_init(den);
    fmpz_init(t);

    /* Compute products in v3 */
    for (i = 0; i < len; i++)
    {
        acb_mul(&v3[i], &v1[i], &v2[i], precs[i]);
    }

    acb_theta_jet_orders(orders, ord, g);
    for (j = 0; j < nb; j++)
    {
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

        /* Loop over lattice points */
        for (i = 0; i < len; i++)
        {
            fmpz_set_si(t, coords[0] + i);
            fmpz_pow_ui(t, t, orders[j * g]);
            acb_mul_fmpz(x, &v3[i], t, precs[i]);
            acb_add(&aux[j], &aux[j], x, prec);
        }

        /* Get cofactor * (2 i pi)^k * num/den */
        acb_const_pi(x, prec);
        acb_mul_onei(x, x);
        acb_mul_2exp_si(x, x, 1);
        acb_pow_ui(x, x, acb_theta_jet_total_order(orders + j * g, g), prec);
        acb_mul(x, x, cofactor, prec);
        acb_mul_fmpz(x, x, num, prec);
        acb_div_fmpz(x, x, den, prec);
        acb_mul(&aux[j], &aux[j], x, prec);
    }
    _acb_vec_add(dth, dth, aux, nb, fullprec);

    flint_free(orders);
    _acb_vec_clear(v3, len);
    _acb_vec_clear(aux, nb);
    acb_clear(x);
    fmpz_clear(num);
    fmpz_clear(den);
    fmpz_clear(t);
}

static void
acb_theta_jet_naive_00_gen(acb_ptr dth, acb_srcptr z, const acb_mat_t tau,
    slong ord, slong prec)
{
    slong g = acb_mat_nrows(tau);
    acb_theta_eld_t E;
    acb_theta_precomp_t D;
    acb_t c;
    arb_t u;
    slong nb = acb_theta_jet_nb(ord, g);

    acb_theta_eld_init(E, g, g);
    acb_theta_precomp_init(D, 1, g);
    acb_init(c);
    arb_init(u);

    acb_theta_jet_ellipsoid(E, u, z, tau, ord, prec);
    prec = acb_theta_naive_fullprec(E, prec);
    acb_theta_precomp_set(D, z, tau, E, prec);
    acb_one(c);

    acb_theta_naive_worker(dth, nb, c, u, E, D, 0, ord, prec, worker_dim1);

    acb_theta_eld_clear(E);
    acb_theta_precomp_clear(D);
    acb_clear(c);
    arb_clear(u);
}

void
acb_theta_jet_naive_00(acb_ptr dth, acb_srcptr z, const acb_mat_t tau,
    slong ord, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong nb = acb_theta_jet_nb(ord, g);
    acb_ptr res;

    if (g == 1)
    {
        res = _acb_vec_init(4 * nb);

        acb_modular_theta_jet(res, res + nb, res + 2 * nb, res + 3 * nb,
            z, acb_mat_entry(tau, 0, 0), nb, prec);
        _acb_vec_set(dth, res + 2 * nb, nb);

        _acb_vec_clear(res, 4 * nb);
    }
    else
    {
        acb_theta_jet_naive_00_gen(dth, z, tau, ord, prec);
    }
}
