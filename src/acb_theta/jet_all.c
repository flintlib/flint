/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "acb_poly.h"
#include "acb_mat.h"
#include "acb_theta.h"

void
acb_theta_jet_all(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau,
    slong ord, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n2 = 1 << (2 * g);
    slong nbth = acb_theta_jet_nb(ord, g);
    fmpz_mat_t mat;
    acb_mat_t new_tau, N, ct;
    acb_ptr new_zs, exps, cs, aux, units, jet;
    arb_ptr rs, r;
    acb_t s, t;
    ulong ab;
    ulong * ch;
    slong * e;
    slong kappa, j;
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
    aux = _acb_vec_init(n2 * nb * nbth);
    units = _acb_vec_init(8);
    jet = _acb_vec_init(nbth);
    rs = _arb_vec_init(nb * g);
    r = _arb_vec_init(g);
    acb_init(s);
    acb_init(t);
    ch = flint_malloc(n2 * sizeof(ulong));
    e = flint_malloc(n2 * sizeof(slong));

    res = acb_theta_reduce_tau(new_zs, new_tau, mat, N, ct, exps, zs, nb, tau, prec);
    if (res)
    {
        res = acb_theta_reduce_z(new_zs, rs, cs, new_zs, nb, new_tau, prec);
    }

    if (res)
    {
        /* Setup */
        _acb_vec_unit_roots(units, 8, 8, prec);
        kappa = acb_siegel_kappa(s, mat, new_tau, prec);
        acb_theta_char_table(ch, e, mat, -1);

        acb_theta_jet_all_notransform(aux, new_zs, nb, new_tau, ord, prec);

        /* Account for reduce_z */
        for (j = 0; j < nb; j++)
        {
            _acb_vec_scalar_mul(aux + j * n2 * nbth, aux + j * n2 * nbth,
                n2 * nbth, &cs[j], prec);
            _arb_vec_neg(r, rs + j * g, g);
            _arb_vec_scalar_mul_2exp_si(r, r, g, 1);
            acb_theta_jet_exp_pi_i(jet, r, ord, g, prec);
            for (ab = 0; ab < n2; ab++)
            {
                acb_theta_jet_mul(aux + j * n2 * nbth + ab * nbth,
                    aux + j * n2 * nbth + ab * nbth, jet, ord, g, prec);
                /* No signs because 2r is divisible by 4 */
            }
        }

        /* Account for reduce_tau */
        for (j = 0; j < nb; j++)
        {
            acb_theta_jet_exp_qf(jet, zs + j * g, N, ord, prec);

            for (ab = 0; ab < n2; ab++)
            {
                acb_mul(t, s, &units[(kappa + e[ab]) % 8], prec);
                _acb_vec_scalar_mul(th + j * n2 * nbth + ab * nbth,
                    aux + j * n2 * nbth + ch[ab] * nbth, nbth, t, prec);
                acb_theta_jet_compose(th + j * n2 * nbth + ab * nbth,
                    th + j * n2 * nbth + ab * nbth, ct, ord, prec);
                acb_theta_jet_mul(th + j * n2 * nbth + ab * nbth,
                    th + j * n2 * nbth + ab * nbth, jet, ord, g, prec);
            }
        }
    }
    else
    {
        _acb_vec_indeterminate(th, nb * n2 * nbth);
    }

    fmpz_mat_clear(mat);
    acb_mat_clear(new_tau);
    acb_mat_clear(N);
    acb_mat_clear(ct);
    _acb_vec_clear(new_zs, nb * g);
    _acb_vec_clear(exps, nb);
    _acb_vec_clear(cs, nb);
    _acb_vec_clear(aux, nb * n2 * nbth);
    _acb_vec_clear(units, 8);
    _acb_vec_clear(jet, nbth);
    _arb_vec_clear(rs, nb * g);
    _arb_vec_clear(r, g);
    acb_clear(s);
    acb_clear(t);
    flint_free(ch);
    flint_free(e);
}
