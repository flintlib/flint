/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "arb_mat.h"
#include "acb_mat.h"
#include "acb_theta.h"

/* Use duplication formula in case of squared theta values */

void
acb_theta_all_notransform(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau,
    int sqr, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;

    if (nb <= 0)
    {
        return;
    }

    if (sqr)
    {
        acb_ptr new_zs, new_th;
        acb_mat_t new_tau;
        slong j;
        int add_zero;

        add_zero = !_acb_vec_is_zero(zs, g);
        new_zs = _acb_vec_init((nb + add_zero) * g);
        new_th = _acb_vec_init((nb + add_zero) * n);
        acb_mat_init(new_tau, g, g);

        _acb_vec_scalar_mul_2exp_si(new_zs + add_zero * g, zs, nb * g, 1);
        acb_mat_scalar_mul_2exp_si(new_tau, tau, 1);
        acb_theta_ql_jet(new_th, new_zs, nb + add_zero, new_tau,
            0, 0, prec);

        for (j = 0; j < nb; j++)
        {
            acb_theta_agm_mul_all(th + j * n * n, new_th,
                new_th + (j + add_zero) * n, g, prec);
        }
    }
    else
    {
        acb_theta_jet_notransform(th, zs, nb, tau, 0, 0, 1, prec);
    }
}
