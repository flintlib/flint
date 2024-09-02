/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_mat.h"
#include "acb_theta.h"

void
acb_theta_00(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    fmpz_mat_t mat, gamma;
    acb_mat_t new_tau, c, cinv, N;
    acb_ptr new_zs, y;
    acb_ptr aux, units;
    acb_t s, t;
    slong kappa, e, ab;
    slong j;

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
    y = _acb_vec_init(g);
    aux = _acb_vec_init(nb);
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
        /* todo: reduce z here. */

        sp2gz_inv(mat, mat);
        ab = acb_theta_transform_char(&e, mat, 0);

        acb_theta_one_notransform(aux, new_zs, nb, new_tau, ab, prec);

        kappa = acb_theta_transform_kappa(s, mat, new_tau, prec);

        fmpz_mat_window_init(gamma, mat, g, 0, 2 * g, g);
        acb_mat_set_fmpz_mat(N, gamma);
        acb_mat_mul(N, c, N, prec);
        fmpz_mat_window_clear(gamma);

        for (j = 0; j < nb; j++)
        {
            acb_mat_vector_mul_col(y, N, new_zs + j * g, prec);
            acb_dot(t, NULL, 0, new_zs + j * g, 1, y, 1, g, prec);
            acb_exp_pi_i(t, t, prec);
            acb_mul(t, t, s, prec);
            acb_mul(t, t, &units[(kappa + e) % 8], prec);
            acb_mul(&th[j], &aux[j], t, prec);
        }
    }
    else
    {
        /* todo: replace with upper bound */
        _acb_vec_indeterminate(th, nb);
    }

    fmpz_mat_clear(mat);
    acb_mat_clear(new_tau);
    acb_mat_clear(c);
    acb_mat_clear(cinv);
    acb_mat_clear(N);
    _acb_vec_clear(new_zs, nb * g);
    _acb_vec_clear(y, g);
    _acb_vec_clear(aux, nb);
    _acb_vec_clear(units, 8);
    acb_clear(s);
    acb_clear(t);
}
