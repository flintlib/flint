/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

/* Given values of f at (x_1 + eps zeta^{n_1}, ..., x_g + eps zeta^{n_g}), make
   Fourier transforms to get Taylor coefficients and add error bound */

void
acb_theta_jet_fd(acb_ptr dth, const arf_t eps, const arb_t c,
    const arb_t rho, acb_srcptr val, slong ord, slong g, slong prec)
{
    acb_ptr aux;
    arb_t t;
    slong nb_max = acb_theta_jet_nb(ord, g);
    slong b = ord + 1;
    slong* orders;
    slong k, j, i, l, nb, ind;

    aux = _acb_vec_init(n_pow(b, g));
    arb_init(t);
    orders = flint_malloc(g * nb_max * sizeof(slong));

    acb_theta_jet_fourier(aux, val, ord, g, prec);

    ind = 0;
    for (k = 0; k <= ord; k++)
    {
        /* Get list of orders */
        nb = acb_theta_jet_nb(k, g);
        acb_theta_jet_orders(orders, k, g);

        /* Get Taylor coefficients, divide by eps^k */
        for (j = 0; j < nb; j++)
        {
            l = 0;
            for (i = 0; i < g; i++)
            {
                l *= b;
                l += orders[j * g + i];
            }
            acb_set(&dth[ind + j], &aux[l]);
        }
        arb_set_arf(t, eps);
        arb_pow_ui(t, t, k, prec);
        _acb_vec_scalar_div_arb(dth + ind, dth + ind, nb, t, prec);

        /* Add error bound */
        arb_one(t);
        arb_mul_2exp_si(t, t, -prec);
        for (j = 0; j < nb; j++)
        {
            acb_add_error_arb(&dth[ind + j], t);
        }

        ind += nb;
    }

    _acb_vec_clear(aux, n_pow(b, g));
    arb_clear(t);
    flint_free(orders);
}
