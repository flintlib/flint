/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

/* val is organized as follows: writing j = a_{g-1}...a_0 in basis ord + 1, val[j]
   corresponds to coefficient (a_0,...,a_{g-1}) */

void
acb_theta_jet_fourier(acb_ptr res, acb_srcptr val, slong ord, slong g, slong prec)
{
    acb_ptr aux;
    acb_ptr zetas, y;
    acb_poly_t pol;
    slong b = ord + 1;
    slong k, j, i, l;

    aux = _acb_vec_init(n_pow(b, g));
    zetas = _acb_vec_init(b);
    y = _acb_vec_init(b);
    acb_poly_init(pol);

    _acb_vec_unit_roots(zetas, -b, b, prec);
    _acb_vec_set(aux, val, n_pow(b, g));

    for (k = 0; k < g; k++)
    {
        /* Fourier transform on variable number k */
        for (j = 0; j < n_pow(b, g - 1); j++)
        {
            /* Make polynomial in X_k of degree ord */
            acb_poly_zero(pol);
            for (i = 0; i < b; i++)
            {
                l = (j % n_pow(b, k)) + i * n_pow(b, k)
                    + n_pow(b, k + 1) * (j / n_pow(b, k));
                acb_poly_set_coeff_acb(pol, i, &aux[l]);
            }
            /* Evaluate and update aux */
            acb_poly_evaluate_vec_fast(y, pol, zetas, b, prec);
            for (i = 0; i < b; i++)
            {
                l = (j % n_pow(b, k)) + i * n_pow(b, k)
                    + n_pow(b, k + 1) * (j / n_pow(b, k));
                acb_set(&aux[l], &y[i]);
            }
        }
    }
    _acb_vec_set(res, aux, n_pow(b, g));

    _acb_vec_clear(aux, n_pow(b, g));
    _acb_vec_clear(zetas, b);
    _acb_vec_clear(y, b);
    acb_poly_clear(pol);
}
