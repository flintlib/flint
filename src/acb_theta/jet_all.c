/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_mat.h"
#include "acb_theta.h"

/* Compute jet of exp (z^T N z) */
static void
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

void
acb_theta_jet_all(acb_ptr dth, acb_srcptr z, const acb_mat_t tau, slong ord, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n2 = 1 << (2 * g);
    slong nb = acb_theta_jet_nb(ord, g);
    fmpz_mat_t mat, gamma;
    acb_mat_t w, c, cinv, N;
    acb_ptr aux, x, units;
    acb_t s, t;
    ulong ab, image_ab;
    slong kappa, e;

    fmpz_mat_init(mat, 2 * g, 2 * g);
    acb_mat_init(w, g, g);
    acb_mat_init(c, g, g);
    acb_mat_init(cinv, g, g);
    acb_mat_init(N, g, g);
    x = _acb_vec_init(g);
    aux = _acb_vec_init(n2 * nb);
    units = _acb_vec_init(8);
    acb_init(s);
    acb_init(t);

    acb_siegel_reduce(mat, tau, prec);
    acb_siegel_transform_cocycle_inv(w, c, cinv, mat, tau, prec);
    _acb_vec_unit_roots(units, 8, 8, prec);

    if (acb_siegel_is_reduced(w, -10, prec))
    {
        sp2gz_inv(mat, mat);
        acb_mat_transpose(cinv, cinv);
        acb_mat_vector_mul_col(x, cinv, z, prec);

        acb_theta_jet_ql_all(aux, x, w, ord, prec);

        kappa = acb_theta_transform_kappa(s, mat, w, prec);
        for (ab = 0; ab < n2; ab++)
        {
            image_ab = acb_theta_transform_char(&e, mat, ab);
            acb_mul(t, s, &units[(kappa + e) % 8], prec);
            _acb_vec_scalar_mul(dth + ab * nb, aux + image_ab * nb, nb, t, prec);
            acb_theta_jet_compose(dth + ab * nb, dth + ab * nb, cinv, ord, prec);
        }

        fmpz_mat_window_init(gamma, mat, g, 0, 2 * g, g);
        acb_mat_set_fmpz_mat(N, gamma);
        acb_mat_mul(N, N, cinv, prec);
        acb_const_pi(t, prec);
        acb_mul_onei(t, t);
        acb_mat_scalar_mul_acb(N, N, t, prec);
        fmpz_mat_window_clear(gamma);

        acb_theta_jet_exp_qf(aux, z, N, ord, prec);

        for (ab = 0; ab < n2; ab++)
        {
            acb_theta_jet_mul(dth + ab * nb, dth + ab * nb, aux, ord, g, prec);
        }
    }
    else
    {
        _acb_vec_indeterminate(dth, n2 * nb);
    }


    fmpz_mat_clear(mat);
    acb_mat_clear(w);
    acb_mat_clear(c);
    acb_mat_clear(cinv);
    acb_mat_clear(N);
    _acb_vec_clear(x, g);
    _acb_vec_clear(aux, n2 * nb);
    _acb_vec_clear(units, 8);
    acb_clear(s);
    acb_clear(t);
}
