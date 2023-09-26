/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void
acb_theta_jet_naive_ind(acb_ptr dth, ulong ab, acb_srcptr z, const acb_mat_t tau,
    slong ord, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong nb = acb_theta_jet_nb(ord, g);
    ulong a = ab >> g;
    ulong b = ab;
    slong* orders;
    acb_ptr new_z, v, w, aux;
    acb_t c, x;
    fmpz_t m;
    slong j, k, l;
    int le;

    orders = flint_malloc(nb * g * sizeof(slong));
    new_z = _acb_vec_init(g);
    v = _acb_vec_init(g);
    w = _acb_vec_init(g);
    aux = _acb_vec_init(nb);
    acb_init(c);
    acb_init(x);
    fmpz_init(m);

    acb_theta_char_get_acb(v, a, g);
    acb_mat_vector_mul_col(v, tau, v, prec); /* v = tau.a/2 */
    acb_theta_char_get_acb(new_z, b, g);
    _acb_vec_add(new_z, new_z, v, g, prec);
    _acb_vec_add(new_z, new_z, z, g, prec);

    acb_theta_jet_naive_00(aux, new_z, tau, ord, prec);

    /* Multiply everything by exponential factor */
    acb_theta_char_dot_acb(c, a, v, g, prec);
    acb_theta_char_get_acb(w, b, g);
    _acb_vec_add(w, w, z, g, prec);
    acb_theta_char_dot_acb(x, a, w, g, prec);
    acb_mul_2exp_si(x, x, 1);
    acb_add(x, x, c, prec);
    acb_exp_pi_i(x, x, prec);
    _acb_vec_scalar_mul(aux, aux, nb, x, prec);

    /* Make linear combinations */
    acb_theta_jet_orders(orders, ord, g);
    _acb_vec_zero(dth, nb);
    for (j = 0; j < nb; j++)
    {
        for (k = 0; k < nb; k++)
        {
            le = 1;
            for (l = 0; (l < g && le); l++)
            {
                if (orders[k * g + l] > orders[j * g + l])
                {
                    le = 0;
                }
            }
            if (!le)
            {
                continue;
            }

            acb_one(x);
            for (l = 0; l < g; l++)
            {
                acb_set_si(c, (a >> (g - 1 - l)) % 2);
                acb_mul_2exp_si(c, c, -1);
                acb_pow_ui(c, c, orders[j * g + l] - orders[k * g + l], prec);
                fmpz_fac_ui(m, orders[j * g + l] - orders[k * g + l]);
                acb_div_fmpz(c, c, m, prec);
                acb_mul(x, x, c, prec);
            }
            acb_const_pi(c, prec);
            acb_mul_onei(c, c);
            acb_mul_2exp_si(c, c, 1);
            acb_pow_ui(c, c, acb_theta_jet_total_order(orders + j * g, g)
                - acb_theta_jet_total_order(orders + k * g, g), prec);
            acb_mul(x, x, c, prec);
            acb_addmul(&dth[j], &aux[k], x, prec);
        }
    }

    flint_free(orders);
    _acb_vec_clear(new_z, g);
    _acb_vec_clear(v, g);
    _acb_vec_clear(w, g);
    _acb_vec_clear(aux, nb);
    acb_clear(c);
    acb_clear(x);
    fmpz_clear(m);
}
