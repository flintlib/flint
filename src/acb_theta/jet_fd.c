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
acb_theta_jet_fd(acb_ptr dth, const arf_t eps, const arf_t err, acb_srcptr val,
    slong ord, slong g, slong prec)
{
    acb_ptr aux;
    arb_t t;
    slong nb = acb_theta_jet_nb(ord, g);
    slong b = ord + 1;
    slong* tups;
    slong j, i, l;
    slong k = 0;

    aux = _acb_vec_init(n_pow(b, g));
    arb_init(t);
    tups = flint_malloc(g * nb * sizeof(slong));

    acb_theta_jet_fourier(aux, val, ord, g, prec);
    arb_set_si(t, n_pow(b, g));
    _acb_vec_scalar_div_arb(aux, aux, n_pow(b, g), t, prec);

    acb_theta_jet_tuples(tups, ord, g);

    /* Get Taylor coefficients, divide by eps^k */
    k = 0;
    arb_one(t);
    for (j = 0; j < nb; j++)
    {
        l = 0;
        for (i = 0; i < g; i++)
        {
            l *= b;
            l += tups[j * g + i];
        }
        acb_set(&dth[j], &aux[l]);

        if (acb_theta_jet_total_order(tups + j * g, g) > k)
        {
            k++;
            arb_mul_arf(t, t, eps, prec);
        }
        acb_div_arb(&dth[j], &dth[j], t, prec);
        acb_add_error_arf(&dth[j], err);
    }

    _acb_vec_clear(aux, n_pow(b, g));
    arb_clear(t);
    flint_free(tups);
}
