/*
    Copyright (C) 2025 Jean Kieffer

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
acb_theta_jet(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau,
    slong ord, ulong ab, int all, int sqr, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    slong nbth = (all ? n * n : 1);
    slong nbjet = acb_theta_jet_nb(ord, g);
    slong nbu;
    fmpz_mat_t mat;
    acb_mat_t new_tau, N, ct;
    acb_ptr new_zs, exps, cs, aux, units, jet;
    arb_ptr rs, r;
    acb_t s, t;
    ulong * ch;
    slong * e;
    slong kappa, j, k;
    int res;

    if (nb <= 0)
    {
        return;
    }
    sqr = sqr && (ord == 0);
    nbu = (sqr ? 4 : 8);

    fmpz_mat_init(mat, 2 * g, 2 * g);
    acb_mat_init(new_tau, g, g);
    acb_mat_init(N, g, g);
    acb_mat_init(ct, g, g);
    new_zs = _acb_vec_init(nb * g);
    exps = _acb_vec_init(nb);
    cs = _acb_vec_init(nb);
    aux = _acb_vec_init(nbth * nb * nbjet);
    units = _acb_vec_init(nbu);
    jet = _acb_vec_init(nbjet);
    rs = _arb_vec_init(nb * g);
    r = _arb_vec_init(g);
    acb_init(s);
    acb_init(t);
    ch = flint_malloc(nbth * sizeof(ulong));
    e = flint_malloc(nbth * sizeof(slong));

    res = acb_theta_reduce_tau(new_zs, new_tau, mat, N, ct, exps, zs, nb, tau, prec);
    if (res)
    {
        res = acb_theta_reduce_z(new_zs, rs, cs, new_zs, nb, new_tau, prec);
    }

    if (res)
    {
        /* Setup */
        _acb_vec_unit_roots(units, nbu, nbu, prec);
        kappa = acb_siegel_kappa(s, mat, new_tau, sqr, prec);
        acb_theta_char_table(ch, e, mat, ab, all);

        acb_theta_jet_notransform(aux, new_zs, nb, new_tau, ord, *ch, all, sqr, prec);

        /* Account for reduce_z */
        for (j = 0; j < nb; j++)
        {
            if (sqr)
            {
                acb_sqr(&cs[j], &cs[j], prec);
            }
            _acb_vec_scalar_mul(aux + j * nbth * nbjet, aux + j * nbth * nbjet,
                nbth * nbjet, &cs[j], prec);
            _arb_vec_neg(r, rs + j * g, g);
            _arb_vec_scalar_mul_2exp_si(r, r, g, 1);
            acb_theta_jet_exp_pi_i(jet, r, ord, g, prec);
            for (k = 0; k < nbth; k++)
            {
                acb_theta_jet_mul(aux + j * nbth * nbjet + k * nbjet,
                    aux + j * nbth * nbjet + k * nbjet, jet, ord, g, prec);
                /* No signs because 2r is divisible by 4 */
            }
        }

        /* Account for reduce_tau */
        for (j = 0; j < nb; j++)
        {
            acb_theta_jet_exp_qf(jet, zs + j * g, N, ord, prec);
            if (sqr)
            {
                acb_sqr(&jet[0], &jet[0], prec);
            }

            for (k = 0; k < nbth; k++)
            {
                acb_mul(t, s, &units[(kappa + e[k]) % (sqr ? 4 : 8)], prec);
                _acb_vec_scalar_mul(th + j * nbth * nbjet + k * nbjet,
                    aux + j * nbth * nbjet + (all ? ch[k] : 0) * nbjet,
                    nbjet, t, prec);
                acb_theta_jet_compose(th + j * nbth * nbjet + k * nbjet,
                    th + j * nbth * nbjet + k * nbjet, ct, ord, prec);
                acb_theta_jet_mul(th + j * nbth * nbjet + k * nbjet,
                    th + j * nbth * nbjet + k * nbjet, jet, ord, g, prec);
            }
        }
    }
    else
    {
        /* Should not happen in tests */
        _acb_vec_indeterminate(th, nb * nbth * nbjet);
    }

    fmpz_mat_clear(mat);
    acb_mat_clear(new_tau);
    acb_mat_clear(N);
    acb_mat_clear(ct);
    _acb_vec_clear(new_zs, nb * g);
    _acb_vec_clear(exps, nb);
    _acb_vec_clear(cs, nb);
    _acb_vec_clear(aux, nb * nbth * nbjet);
    _acb_vec_clear(units, nbu);
    _acb_vec_clear(jet, nbjet);
    _arb_vec_clear(rs, nb * g);
    _arb_vec_clear(r, g);
    acb_clear(s);
    acb_clear(t);
    flint_free(ch);
    flint_free(e);
}
