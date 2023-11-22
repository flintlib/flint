/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dft.h"
#include "acb_theta.h"

/* Given values of f at (x_1 + eps zeta^{n_1}, ..., x_g + eps zeta^{n_g}), make
   Fourier transforms to get Taylor coefficients and add error bound */

void
acb_theta_jet_ql_finite_diff(acb_ptr dth, const arf_t eps, const arf_t err,
    const arb_t rho, acb_srcptr val, slong ord, slong g, slong prec)
{
    slong nb = acb_theta_jet_nb(ord, g);
    slong lp = ACB_THETA_LOW_PREC;
    acb_ptr aux;
    arb_t t, e;
    slong b = ord + 1;
    slong * tups;
    slong * cyc;
    slong j, i, l;
    slong k;

    aux = _acb_vec_init(n_pow(b, g));
    arb_init(t);
    arb_init(e);
    tups = flint_malloc(g * nb * sizeof(slong));
    cyc = flint_malloc(g * sizeof(slong));

    for (j = 0; j < g; j++)
    {
        cyc[j] = b;
    }
    acb_dft_prod(aux, val, cyc, g, prec);
    arb_set_si(t, n_pow(b, g));
    _acb_vec_scalar_div_arb(aux, aux, n_pow(b, g), t, prec);

    /* Get Taylor coefficients, divide by eps^k, add error */
    acb_theta_jet_tuples(tups, ord, g);
    k = 0;
    arb_one(t);
    arb_pow_ui(e, rho, ord, lp);
    arb_mul_arf(e, e, err, lp);
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
            arb_pow_ui(e, rho, ord - k, lp);
            arb_mul_arf(e, e, err, lp);
        }
        acb_div_arb(&dth[j], &dth[j], t, prec);
        acb_add_error_arb(&dth[j], e);
    }

    _acb_vec_clear(aux, n_pow(b, g));
    arb_clear(t);
    arb_clear(e);
    flint_free(tups);
    flint_free(cyc);
}
