/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"
#include "acb_modular.h"
#include "acb_theta.h"

static void
worker(acb_ptr dth, acb_srcptr v1, acb_srcptr v2, const slong * precs, slong len,
    const acb_t cofactor, const slong * coords, slong ord, slong g, slong prec, slong fullprec)
{
    slong nb = acb_theta_jet_nb(ord, g);
    slong * tups;
    acb_ptr v3, aux;
    acb_t x;
    fmpz_t num, t;
    slong j, i;

    tups = flint_malloc(g * nb * sizeof(slong));
    v3 = _acb_vec_init(len);
    aux = _acb_vec_init(nb);
    acb_init(x);
    fmpz_init(num);
    fmpz_init(t);

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
            acb_add(&aux[j], &aux[j], x, prec);
        }

        /* Multiply by cofactor * num */
        acb_mul_fmpz(x, cofactor, num, prec);
        acb_mul(&aux[j], &aux[j], x, prec);
    }
    _acb_vec_add(dth, dth, aux, nb, fullprec);

    flint_free(tups);
    _acb_vec_clear(v3, len);
    _acb_vec_clear(aux, nb);
    acb_clear(x);
    fmpz_clear(num);
    fmpz_clear(t);
}

static void
acb_theta_jet_naive_00_gen(acb_ptr dth, acb_srcptr z, const acb_mat_t tau,
    slong ord, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong nb = acb_theta_jet_nb(ord, g);
    slong * tups;
    acb_theta_eld_t E;
    arb_mat_t C;
    arf_t R2, eps;
    acb_ptr new_z, aux;
    arb_ptr v, a;
    acb_t c;
    arb_t u;
    fmpz_t m, t;
    slong j, k;
    int b;

    tups = flint_malloc(g * nb * sizeof(slong));
    acb_theta_eld_init(E, g, g);
    arb_mat_init(C, g, g);
    arf_init(R2);
    arf_init(eps);
    new_z = _acb_vec_init(g);
    aux = _acb_vec_init(nb);
    a = _arb_vec_init(g);
    v = _arb_vec_init(g);
    acb_init(c);
    arb_init(u);
    fmpz_init(m);
    fmpz_init(t);

    acb_siegel_cho(C, tau, prec);
    acb_theta_naive_reduce(v, new_z, a, c, u, z, 1, tau, prec);
    acb_theta_jet_naive_radius(R2, eps, C, v, ord, prec);
    b = acb_theta_eld_set(E, C, R2, v);

    if (b)
    {
        acb_theta_naive_worker(dth, nb, new_z, 1, tau, E, ord, prec, worker);

        arb_mul_arf(u, u, eps, prec);
        for (k = 0; k < nb; k++)
        {
            acb_mul(&dth[k], &dth[k], c, prec);
            acb_add_error_arb(&dth[k], u);
        }

        acb_theta_jet_tuples(tups, ord, g);
        for (k = 0; k < nb; k++)
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
            acb_mul(&dth[k], &dth[k], c, prec);
        }

        _arb_vec_neg(a, a, g);
        _arb_vec_scalar_mul_2exp_si(a, a, g, 1);
        acb_theta_jet_exp_pi_i(aux, a, ord, g, prec);
        acb_theta_jet_mul(dth, dth, aux, ord, g, prec);
    }
    else
    {
        for (k = 0; k < nb; k++)
        {
            acb_indeterminate(&dth[k]);
        }
    }

    flint_free(tups);
    acb_theta_eld_clear(E);
    arb_mat_clear(C);
    arf_clear(R2);
    arf_clear(eps);
    _acb_vec_clear(new_z, g);
    _acb_vec_clear(aux, nb);
    _arb_vec_clear(v, g);
    _arb_vec_clear(a, g);
    acb_clear(c);
    arb_clear(u);
    fmpz_clear(m);
    fmpz_clear(t);
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
