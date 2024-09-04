/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_mat.h"
#include "acb_theta.h"

/* Compute jet of exp (z^T N z) */
void
acb_theta_jet_exp_qf(acb_ptr res, acb_srcptr z, const acb_mat_t N, slong ord, slong prec)
{
    slong g = acb_mat_nrows(N);
    slong nb = acb_theta_jet_nb(ord, g);
    acb_mat_t tp;
    acb_poly_t pol;
    acb_ptr aux;
    acb_ptr y;
    acb_t c;
    slong * tup;
    slong j, k, l, i;

    acb_mat_init(tp, g, g);
    acb_poly_init(pol);
    aux = _acb_vec_init(nb);
    y = _acb_vec_init(g);
    acb_init(c);
    tup = flint_malloc(g * sizeof(slong));

    /* exp((z+h)^T N (z+h)) = exp(z^T N z) exp(z^T (N+N^T) h) exp(h^T N h) */
    _acb_vec_zero(res, nb);
    acb_mat_vector_mul_col(y, N, z, prec);
    acb_dot(&res[0], NULL, 0, z, 1, y, 1, g, prec);
    acb_exp(&res[0], &res[0], prec);

    acb_mat_transpose(tp, N);
    acb_mat_add(tp, tp, N, prec);
    acb_mat_vector_mul_row(y, z, tp, prec);
    for (j = 0; j < g; j++)
    {
        _acb_vec_zero(aux, nb);
        acb_poly_zero(pol);
        acb_poly_set_coeff_acb(pol, 1, &y[j]);
        acb_poly_exp_series(pol, pol, ord + 1, prec);
        for (l = 0; l <= ord; l++)
        {
            for (i = 0; i < g; i++)
            {
                tup[i] = 0;
            }
            tup[j] = l;
            acb_poly_get_coeff_acb(&aux[acb_theta_jet_index(tup, g)], pol, l);
        }
        acb_theta_jet_mul(res, res, aux, ord, g, prec);
    }

    for (j = 0; j < g; j++)
    {
        for (k = j; k < g; k++)
        {
            _acb_vec_zero(aux, nb);
            acb_poly_zero(pol);
            acb_add(c, acb_mat_entry(N, k, j), acb_mat_entry(N, j, k), prec);
            if (j == k)
            {
                acb_mul_2exp_si(c, c, -1);
            }
            acb_poly_set_coeff_acb(pol, 1, c);
            acb_poly_exp_series(pol, pol, (ord / 2) + 1, prec);
            for (l = 0; l <= (ord / 2); l++)
            {
                for (i = 0; i < g; i++)
                {
                    tup[i] = 0;
                }
                tup[j] += l;
                tup[k] += l;
                acb_poly_get_coeff_acb(&aux[acb_theta_jet_index(tup, g)], pol, l);
            }
            acb_theta_jet_mul(res, res, aux, ord, g, prec);
        }
    }

    acb_mat_clear(tp);
    acb_poly_clear(pol);
    _acb_vec_clear(aux, nb);
    _acb_vec_clear(y, g);
    acb_clear(c);
    flint_free(tup);
}
