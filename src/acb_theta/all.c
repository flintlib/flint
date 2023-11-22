/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_mat.h"
#include "acb_theta.h"

void
acb_theta_all(acb_ptr th, acb_srcptr z, const acb_mat_t tau, int sqr, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n2 = 1 << (2 * g);
    fmpz_mat_t mat, gamma;
    acb_mat_t w, c, N;
    acb_ptr x, y, aux, units;
    acb_t s, t;
    ulong ab, image_ab;
    slong kappa, e;

    fmpz_mat_init(mat, 2 * g, 2 * g);
    acb_mat_init(w, g, g);
    acb_mat_init(c, g, g);
    acb_mat_init(N, g, g);
    x = _acb_vec_init(g);
    y = _acb_vec_init(g);
    aux = _acb_vec_init(n2);
    units = _acb_vec_init(8);
    acb_init(s);
    acb_init(t);

    acb_siegel_reduce(mat, tau, prec);
    acb_siegel_transform_z(x, w, mat, z, tau, prec);
    acb_siegel_cocycle(c, mat, tau, prec);
    _acb_vec_unit_roots(units, 8, 8, prec);

    if (acb_siegel_is_reduced(w, -10, prec))
    {
        sp2gz_inv(mat, mat);

        fmpz_mat_window_init(gamma, mat, g, 0, 2 * g, g);
        acb_mat_set_fmpz_mat(N, gamma);
        fmpz_mat_window_clear(gamma);
        acb_mat_mul(N, c, N, prec);
        acb_mat_vector_mul_col(y, N, x, prec);
        acb_dot(t, NULL, 0, x, 1, y, 1, g, prec);

        acb_theta_ql_all(aux, x, w, sqr, prec);

        if (sqr)
        {
            kappa = acb_theta_transform_kappa2(mat);
            acb_siegel_cocycle(c, mat, w, prec);
            acb_mat_det(s, c, prec);
            acb_mul_2exp_si(t, t, 1);
        }
        else
        {
            kappa = acb_theta_transform_kappa(s, mat, w, prec);
        }

        acb_exp_pi_i(t, t, prec);
        acb_mul(s, s, t, prec);

        for (ab = 0; ab < n2; ab++)
        {
            image_ab = acb_theta_transform_char(&e, mat, ab);
            acb_mul(t, s, &units[((sqr ? 2 : 1) * (kappa + e)) % 8], prec);
            acb_mul(&th[ab], &aux[image_ab], t, prec);
        }
    }
    else
    {
        _acb_vec_indeterminate(th, n2);
    }


    fmpz_mat_clear(mat);
    acb_mat_clear(w);
    acb_mat_clear(c);
    acb_mat_clear(N);
    _acb_vec_clear(x, g);
    _acb_vec_clear(y, g);
    _acb_vec_clear(aux, n2);
    _acb_vec_clear(units, 8);
    acb_clear(s);
    acb_clear(t);
}
