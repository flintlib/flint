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
acb_theta_all(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau,
    int sqr, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n2 = 1 << (2 * g);
    fmpz_mat_t mat;
    acb_mat_t new_tau, N, ct;
    acb_ptr new_zs, exps, cs, aux, units;
    arb_ptr rs;
    acb_t s;
    ulong * ch;
    slong * e;
    slong kappa;
    slong j, ab;
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
    aux = _acb_vec_init(n2 * nb);
    units = _acb_vec_init(8);
    rs = _arb_vec_init(nb * g);
    ch = flint_malloc(n2 * sizeof(ulong));
    e = flint_malloc(n2 * sizeof(slong));
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
        _acb_vec_unit_roots(units, 8, 8, prec);
        acb_theta_char_table(ch, e, mat, -1);
        if (sqr)
        {
            kappa = acb_siegel_kappa2(mat);
            acb_mat_det(s, ct, prec);
        }
        else
        {
            kappa = acb_siegel_kappa(s, mat, new_tau, prec);
        }

        acb_theta_all_notransform(aux, new_zs, nb, new_tau, sqr, prec);

        /* Account for reduce_z */
        for (j = 0; j < nb; j++)
        {
            if (sqr)
            {
                acb_sqr(&cs[j], &cs[j], prec);
            }
            _acb_vec_scalar_mul(aux + j * n2, aux + j * n2, n2, &cs[j], prec);
        }

        /* Account for reduce_tau */
        for (j = 0; j < nb; j++)
        {
            if (sqr)
            {
                acb_sqr(&exps[j], &exps[j], prec);
            }
            for (ab = 0; ab < n2; ab++)
            {
                acb_mul(&th[j * n2 + ab], &aux[j * n2 + ch[ab]], &exps[j], prec);
                acb_mul(&th[j * n2 + ab], &th[j * n2 + ab], s, prec);
                acb_mul(&th[j * n2 + ab], &th[j * n2 + ab],
                    &units[((sqr ? 2 : 1) * (kappa + e[ab])) % 8], prec);
            }
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
            for (ab = 0; ab < n2; ab++)
            {
                arb_zero_pm_one(acb_realref(&th[j * n2 + ab]));
                arb_zero_pm_one(acb_imagref(&th[j * n2 + ab]));
            }
            _acb_vec_scalar_mul_arb(th + j * n2, th + j * n2, n2, c, prec);
        }

        arb_clear(c);
        arb_clear(rho);
    }


    fmpz_mat_clear(mat);
    acb_mat_clear(new_tau);
    acb_mat_clear(ct);
    acb_mat_clear(N);
    _acb_vec_clear(new_zs, nb * g);
    _acb_vec_clear(exps, nb);
    _acb_vec_clear(cs, nb);
    _acb_vec_clear(aux, n2 * nb);
    _acb_vec_clear(units, 8);
    _arb_vec_clear(rs, nb * g);
    acb_clear(s);
    flint_free(e);
    flint_free(ch);
}
