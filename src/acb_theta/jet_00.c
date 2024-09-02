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
acb_theta_jet_00(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau,
    slong ord, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong nbth = acb_theta_jet_nb(ord, g);
    fmpz_mat_t mat, gamma;
    acb_mat_t new_tau, c, cinv, N;
    acb_ptr new_zs, aux, units;
    acb_t s, t;
    ulong image_ab;
    slong kappa, e, j;

    if (nb <= 0)
    {
        return;
    }

    fmpz_mat_init(mat, 2 * g, 2 * g);
    acb_mat_init(new_tau, g, g);
    acb_mat_init(c, g, g);
    acb_mat_init(cinv, g, g);
    acb_mat_init(N, g, g);
    new_zs = _acb_vec_init(nb * g);
    aux = _acb_vec_init(nb * nbth);
    units = _acb_vec_init(8);
    acb_init(s);
    acb_init(t);

    acb_siegel_reduce(mat, tau, prec);
    acb_siegel_transform_cocycle_inv(new_tau, c, cinv, mat, tau, prec);
    _acb_vec_unit_roots(units, 8, 8, prec);

    acb_mat_transpose(cinv, cinv);
    for (j = 0; j < nb; j++)
    {
        acb_mat_vector_mul_col(new_zs + j * g, cinv, zs + j * g, prec);
    }

    if (acb_siegel_is_reduced(new_tau, -10, prec))
    {
        /* todo: reduce z here */

        sp2gz_inv(mat, mat);
        kappa = acb_theta_transform_kappa(s, mat, new_tau, prec);
        image_ab = acb_theta_transform_char(&e, mat, 0);

        fmpz_mat_window_init(gamma, mat, g, 0, 2 * g, g);
        acb_mat_set_fmpz_mat(N, gamma);
        acb_mat_mul(N, N, cinv, prec);
        acb_const_pi(t, prec);
        acb_mul_onei(t, t);
        acb_mat_scalar_mul_acb(N, N, t, prec);
        fmpz_mat_window_clear(gamma);

        acb_theta_jet_one_notransform(aux, new_zs, nb, new_tau, ord, image_ab, prec);

        for (j = 0; j < nb; j++)
        {
            acb_mul(t, s, &units[(kappa + e) % 8], prec);
            _acb_vec_scalar_mul(th + j * nbth, aux + j * nbth, nbth, t, prec);
            acb_theta_jet_compose(th + j * nbth, th + j * nbth, cinv, ord, prec);
            acb_theta_jet_exp_qf(aux + j * nbth, zs + j * g, N, ord, prec);
            acb_theta_jet_mul(th + j * nbth, th + j * nbth, aux + j * nbth, ord, g, prec);
        }
    }
    else
    {
        _acb_vec_indeterminate(th, nb * nbth);
    }

    fmpz_mat_clear(mat);
    acb_mat_clear(new_tau);
    acb_mat_clear(c);
    acb_mat_clear(cinv);
    acb_mat_clear(N);
    _acb_vec_clear(new_zs, nb * g);
    _acb_vec_clear(aux, nb * nbth);
    _acb_vec_clear(units, 8);
    acb_clear(s);
    acb_clear(t);
}
