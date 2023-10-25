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
acb_theta_jet_exp(acb_ptr res, acb_srcptr z, const acb_mat_t tau, ulong a, ulong b,
    slong ord, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong nb = acb_theta_jet_nb(ord, g);
    slong* tups;
    acb_ptr v, w;
    acb_t c, x;
    fmpz_t num, den, t;
    slong k, l;

    tups = flint_malloc(g * nb * sizeof(slong));
    v = _acb_vec_init(g);
    w = _acb_vec_init(g);
    acb_init(c);
    acb_init(x);
    fmpz_init(num);
    fmpz_init(den);
    fmpz_init(t);

    /* Get exponential factor */
    acb_theta_char_get_acb(v, a, g);
    acb_mat_vector_mul_col(v, tau, v, prec);
    acb_theta_char_dot_acb(c, a, v, g, prec);

    acb_theta_char_get_acb(w, b, g);
    _acb_vec_add(w, w, z, g, prec);
    acb_theta_char_dot_acb(x, a, w, g, prec);
    acb_mul_2exp_si(x, x, 1);
    acb_add(x, x, c, prec);
    acb_exp_pi_i(&res[0], x, prec);

    /* Get other coefficients */
    acb_theta_jet_tuples(tups, ord, g);
    for (k = 1; k < nb; k++)
    {
        fmpz_one(num);
        fmpz_one(den);
        for (l = 0; l < g; l++)
        {
            fmpz_ui_pow_ui(t, (a >> (g - 1 - l)) % 2, tups[k * g + l]);
            fmpz_mul(num, num, t);
            fmpz_fac_ui(t, tups[k * g + l]);
            fmpz_mul(den, den, t);
        }

        acb_const_pi(c, prec);
        acb_mul_onei(c, c);
        acb_pow_ui(c, c, acb_theta_jet_total_order(tups + k * g, g), prec);
        acb_mul(&res[k], &res[0], c, prec);

        acb_set_fmpz(c, num);
        acb_div_fmpz(c, c, den, prec);
        acb_mul(&res[k], &res[k], c, prec);
    }

    flint_free(tups);
    _acb_vec_clear(v, g);
    _acb_vec_clear(w, g);
    acb_clear(c);
    acb_clear(x);
    fmpz_clear(num);
    fmpz_clear(den);
    fmpz_clear(t);
}

void
acb_theta_jet_naive_fixed_ab(acb_ptr dth, ulong ab, acb_srcptr z, const acb_mat_t tau,
    slong ord, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong nb = acb_theta_jet_nb(ord, g);
    ulong a = ab >> g;
    ulong b = ab;
    acb_ptr v, new_z, aux;

    v = _acb_vec_init(g);
    new_z = _acb_vec_init(g);
    aux = _acb_vec_init(nb);

    acb_theta_char_get_acb(v, a, g);
    acb_mat_vector_mul_col(new_z, tau, v, prec);
    acb_theta_char_get_acb(v, b, g);
    _acb_vec_add(new_z, new_z, v, g, prec);
    _acb_vec_add(new_z, new_z, z, g, prec);

    acb_theta_jet_exp(aux, z, tau, a, b, ord, prec);
    acb_theta_jet_naive_00(dth, new_z, tau, ord, prec);
    acb_theta_jet_mul(dth, dth, aux, ord, g, prec);

    _acb_vec_clear(new_z, g);
    _acb_vec_clear(v, g);
    _acb_vec_clear(aux, nb);
}
