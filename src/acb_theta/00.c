/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "fmpz_mat.h"
#include "acb_mat.h"
#include "acb_theta.h"

void
acb_theta_00(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    fmpz_mat_t mat;
    acb_mat_t new_tau, N, ct;
    acb_ptr new_zs, exps, cs, units;
    arb_ptr rs;
    acb_t s;
    slong kappa, e;
    ulong ab;
    slong j;
    int res;

    if (nb <= 0)
    {
        return;
    }

    fmpz_mat_init(mat, 2 * g, 2 * g);
    acb_mat_init(new_tau, g, g);
    acb_mat_init(N, g, g);
    acb_mat_init(ct, g, g);
    new_zs = _acb_vec_init(nb * g);
    exps = _acb_vec_init(nb);
    cs = _acb_vec_init(nb);
    units = _acb_vec_init(8);
    rs = _arb_vec_init(nb * g);
    acb_init(s);

    /* Reduce tau then z */
    res = acb_theta_reduce_tau(new_zs, new_tau, mat, N, ct, exps, zs, nb, tau, prec);
    if (res)
    {
        res = acb_theta_reduce_z(new_zs, rs, cs, new_zs, nb, new_tau, prec);
    }

    if (res)
    {
        /* Setup */
        acb_theta_char_table(&ab, &e, mat, 0);
        _acb_vec_unit_roots(units, 8, 8, prec);
        kappa = acb_siegel_kappa(s, mat, new_tau, prec);

        acb_theta_one_notransform(th, new_zs, nb, new_tau, ab, prec);

        /* Account for reduce_z */
        for (j = 0; j < nb; j++)
        {
            acb_mul(&th[j], &th[j], &cs[j], prec);
        }

        /* Account for reduce_tau */
        for (j = 0; j < nb; j++)
        {
            acb_mul(&th[j], &th[j], &exps[j], prec);
            acb_mul(&th[j], &th[j], s, prec);
            acb_mul(&th[j], &th[j], &units[(kappa + e) % 8], prec);
        }
    }
    else
    {
        /* Use sum_bound to avoid returning NaN */
        arb_t c, rho;
        arb_init(c);
        arb_init(rho);

        for (j = 0; j < nb; j++)
        {
            acb_theta_sum_bound(c, rho, zs + j * g, tau, 0);
            arb_zero_pm_one(acb_realref(&th[j]));
            arb_zero_pm_one(acb_imagref(&th[j]));
            acb_mul_arb(&th[j], &th[j], c, prec);
        }

        arb_clear(c);
        arb_clear(rho);
    }

    fmpz_mat_clear(mat);
    acb_mat_clear(new_tau);
    acb_mat_clear(N);
    acb_mat_clear(ct);
    _acb_vec_clear(new_zs, nb * g);
    _acb_vec_clear(exps, nb);
    _acb_vec_clear(cs, nb);
    _acb_vec_clear(units, 8);
    _arb_vec_clear(rs, nb * g);
    acb_clear(s);
}
