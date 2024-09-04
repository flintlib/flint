/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"
#include "acb.h"
#include "acb_mat.h"
#include "acb_theta.h"

void
acb_theta_jet_error(arb_ptr err, acb_srcptr z, const acb_mat_t tau,
    acb_srcptr dth, slong ord, slong prec)
{
    slong g = acb_mat_nrows(tau);
    arb_ptr abs_der;
    arb_mat_t tau_err;
    arb_ptr z_err;
    arb_t e, f;
    arf_t half;
    slong nb = acb_theta_jet_nb(ord, g);
    slong nb_dth = acb_theta_jet_nb(ord + 2, g);
    slong * tups;
    slong * new_tups;
    slong j, l, m, i;

    abs_der = _arb_vec_init(nb_dth);
    arb_mat_init(tau_err, g, g);
    z_err = _arb_vec_init(g);
    arb_init(e);
    arb_init(f);
    arf_init(half);
    tups = flint_malloc(nb * g * sizeof(slong));
    new_tups = flint_malloc(g * sizeof(slong));

    /* Get input errors on z, tau */
    arf_set_si(half, 1);
    arf_mul_2exp_si(half, half, -1);
    for (l = 0; l < g; l++)
    {
        for (m = l; m < g; m++)
        {
            arb_zero(e);
            acb_get_rad_ubound_arf(arb_midref(e), acb_mat_entry(tau, l, m), prec);
            arb_set(arb_mat_entry(tau_err, l, m), e);
        }
        arb_zero(e);
        arb_zero(f);
        /* Use that theta is 1-periodic in re(z) */
        arf_set_mag(arb_midref(e), arb_radref(acb_imagref(&z[l])));
        arf_set_mag(arb_midref(f), arb_radref(acb_realref(&z[l])));
        arf_min(arb_midref(f), arb_midref(f), half);
        arb_sqr(e, e, prec);
        arb_sqr(f, f, prec);
        arb_add(e, e, f, prec);
        arb_sqrt(&z_err[l], e, prec);
    }

    /* We need order ord + 2 to use the heat equation. */
    for (j = 0; j < nb_dth; j++)
    {
        acb_get_abs_ubound_arf(arb_midref(&abs_der[j]), &dth[j], prec);
    }

    /* Loop over tuples to compute the correct bounds */
    acb_theta_jet_tuples(tups, ord, g);
    for (j = 0; j < nb; j++)
    {
        arb_zero(&err[j]);
        /* Add error corresponding to entries of tau */
        for (l = 0; l < g; l++)
        {
            for (m = l; m < g; m++)
            {
                /* Heat equation: d/dzl d/dzm = 2pi i (1 + delta) d/dtaulm */
                for (i = 0; i < g; i++)
                {
                    new_tups[i] = tups[j * g + i];
                }
                new_tups[l] += 1;
                new_tups[m] += 1;
                i = acb_theta_jet_index(new_tups, g);

                arb_mul(e, arb_mat_entry(tau_err, l, m), &abs_der[i], prec);
                arb_const_pi(f, prec);
                if (l == m)
                {
                    arb_mul_2exp_si(f, f, 2);
                    arb_mul_si(e, e, new_tups[l] * (new_tups[l] - 1), prec);
                }
                else
                {
                    arb_mul_2exp_si(f, f, 1);
                    arb_mul_si(e, e, new_tups[l] * new_tups[m], prec);
                }
                arb_div(e, e, f, prec);
                arb_add(&err[j], &err[j], e, prec);
            }
        }
        /* Add error corresponding to entries of z */
        for (l = 0; l < g; l++)
        {
            for (i = 0; i < g; i++)
            {
                new_tups[i] = tups[j * g + i];
            }
            new_tups[l] += 1;
            i = acb_theta_jet_index(new_tups, g);

            arb_mul(e, &z_err[l], &abs_der[i], prec);
            arb_mul_si(e, e, new_tups[l], prec);
            arb_add(&err[j], &err[j], e, prec);
        }
    }

    _arb_vec_clear(abs_der, nb_dth);
    arb_mat_clear(tau_err);
    _arb_vec_clear(z_err, g);
    arb_clear(e);
    arb_clear(f);
    arf_clear(half);
    flint_free(tups);
    flint_free(new_tups);
}
