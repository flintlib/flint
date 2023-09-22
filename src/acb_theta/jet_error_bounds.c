/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void acb_theta_jet_error_bounds(arb_ptr err, acb_srcptr z, const acb_mat_t tau,
    slong ord, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << (2 * g);
    acb_ptr der;
    arb_ptr abs_der;
    arb_mat_t tau_err;
    arb_ptr z_err;
    arb_t e, f;
    slong nb_der = acb_theta_jet_nb(ord + 2, g + 1);
    slong nb;
    slong nb_max = acb_theta_jet_nb(ord, g);
    slong nb_tot = acb_theta_jet_nb(ord, g + 1);
    slong* orders;
    slong* new_orders;
    slong k, j, l, m, a, i, ind, ind1, ind2;

    der = _acb_vec_init(n * nb_der);
    abs_der = _arb_vec_init(n * nb_der);
    arb_mat_init(tau_err, g, g);
    z_err = _arb_vec_init(g);
    arb_init(e);
    arb_init(f);
    orders = flint_malloc(nb_max * g * sizeof(slong));
    new_orders = flint_malloc(g * sizeof(slong));

    /* Get input errors on z, tau */
    for (l = 0; l < g; l++)
    {
        for (m = l; m < g; m++)
        {
            acb_get_rad_ubound_arf(arb_midref(e), acb_mat_entry(tau, l, m), prec);
            arb_set(arb_mat_entry(tau_err, l, m), e);
        }
        acb_get_rad_ubound_arf(arb_midref(e), &z[l], prec);
        arb_set(&z_err[l], e);
    }

    /* We need order ord + 2 to use the heat equation. */
    acb_theta_jet_naive_all(der, z, tau, ord + 2, prec);
    for (k = 0; k < n * nb_der; k++)
    {
        acb_get_abs_ubound_arf(arb_midref(&abs_der[k]), &der[k], prec);
    }

    /* Loop over orders to compute the correct bounds */
    ind = 0;
    ind1 = acb_theta_jet_nb(0, g);
    ind2 = ind1 + acb_theta_jet_nb(1, g);
    for (k = 0; k <= ord; k++)
    {
        nb = acb_theta_jet_nb(k, g);
        acb_theta_jet_orders(orders, k, g);

        for (a = 0; a < n; a++)
        {
            for (j = 0; j < nb; j++)
            {
                arb_zero(&err[a * nb_tot + ind + j]);
               /* Add error corresponding to entries of tau */
                for (l = 0; l < g; l++)
                {
                    for (m = l; m < g; m++)
                    {
                        /* Heat equation: d/dzl d/dzm = 2pi i (1 + delta) d/dtaulm */
                        for (i = 0; i < g; i++)
                        {
                            new_orders[i] = orders[j * g + i];
                        }
                        new_orders[l] += 1;
                        new_orders[m] += 1;
                        i = ind2 + acb_theta_jet_index(new_orders, g);

                        arb_mul(e, arb_mat_entry(tau_err, l, m), &abs_der[a * nb_der + i], prec);
                        arb_const_pi(f, prec);
                        if (l == m)
                        {
                            arb_mul_2exp_si(f, f, 2);
                            arb_mul_si(e, e, new_orders[l] * (new_orders[l] - 1), prec);
                        }
                        else
                        {
                            arb_mul_2exp_si(f, f, 1);
                            arb_mul_si(e, e, new_orders[l] * new_orders[m], prec);
                        }
                        arb_div(e, e, f, prec);
                        arb_add(&err[a * nb_tot + ind + j], &err[a * nb_tot + ind + j], e, prec);
                    }
                }
                /* Add error corresponding to entries of z */
                for (l = 0; l < g; l++)
                {
                    for (i = 0; i < g; i++)
                    {
                        new_orders[i] = orders[j * g + i];
                    }
                    new_orders[l] += 1;
                    i = ind1 + acb_theta_jet_index(new_orders, g);

                    arb_mul(e, &z_err[l], &abs_der[a * nb_der + i], prec);
                    arb_mul_si(e, e, new_orders[l], prec);
                    arb_add(&err[a * nb_tot + ind + j], &err[a * nb_tot + ind + j], e, prec);
                }
            }
        }

        ind += nb;
        ind1 += acb_theta_jet_nb(k + 1, g);
        ind2 += acb_theta_jet_nb(k + 2, g);
    }

    _acb_vec_clear(der, n * nb_der);
    _arb_vec_clear(abs_der, n * nb_der);
    arb_mat_clear(tau_err);
    _arb_vec_clear(z_err, g);
    arb_clear(e);
    flint_free(orders);
    flint_free(new_orders);
}
