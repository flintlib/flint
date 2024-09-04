/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "acb.h"
#include "acb_mat.h"
#include "acb_theta.h"

int
acb_theta_reduce_tau(acb_ptr new_zs, acb_mat_t new_tau, fmpz_mat_t mat, acb_mat_t N,
    acb_mat_t ct, acb_ptr exps, acb_srcptr zs, slong nb, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    fmpz_mat_t gamma;
    acb_mat_t c;
    acb_ptr y;
    acb_t ipi;
    slong j;
    int res;

    FLINT_ASSERT(nb >= 0);

    acb_mat_init(c, g, g);
    y = _acb_vec_init(g);
    acb_init(ipi);

    /* Set new_tau, mat, ct, res */
    acb_siegel_reduce(mat, tau, prec);
    acb_siegel_transform_cocycle_inv(new_tau, c, ct, mat, tau, prec);
    sp2gz_inv(mat, mat);
    acb_mat_transpose(ct, ct);
    res = acb_siegel_is_reduced(new_tau, -10, prec);

    if (res)
    {
        /* Set N */
        fmpz_mat_window_init(gamma, mat, g, 0, 2 * g, g);
        acb_mat_set_fmpz_mat(N, gamma);
        acb_mat_mul(N, N, ct, prec);
        acb_const_pi(ipi, prec);
        acb_mul_onei(ipi, ipi);
        acb_mat_scalar_mul_acb(N, N, ipi, prec);

        /* Set exps and nex_zs */
        for (j = 0; j < nb; j++)
        {
            acb_mat_vector_mul_col(y, N, zs + j * g, prec);
            acb_dot(&exps[j], NULL, 0, zs + j * g, 1, y, 1, g, prec);
            acb_exp(&exps[j], &exps[j], prec);
            acb_mat_vector_mul_col(new_zs + j * g, ct, zs + j * g, prec);
        }
    }

    fmpz_mat_window_clear(gamma);
    acb_mat_clear(c);
    _acb_vec_clear(y, g);
    acb_clear(ipi);
    return res;
}
