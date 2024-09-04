/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "arb_mat.h"
#include "acb_modular.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_sum_jet_all, state)
{
    slong iter;

    /* Test: matches acb_modular_theta_jet in genus 1;
       agrees with product of dim 1 case on diagonal matrices */
    for (iter = 0; iter < 25 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong n2 = 1 << (2 * g);
        slong mprec = 100 + n_randint(state, 100);
        slong prec = mprec + 50;
        slong ord = n_randint(state, 3);
        slong bits = n_randint(state, 4);
        slong nbz = n_randint(state, 3);
        slong nbth = acb_theta_jet_nb(ord, g);
        acb_mat_t tau, tau11;
        acb_ptr zs;
        acb_theta_ctx_tau_t ctx_tau;
        acb_theta_ctx_z_struct * vec;
        acb_ptr th, aux, test;
        slong j, k, l;

        acb_mat_init(tau, g, g);
        acb_mat_init(tau11, 1, 1);
        zs = _acb_vec_init(nbz * g);
        acb_theta_ctx_tau_init(ctx_tau, g);
        vec = acb_theta_ctx_z_vec_init(nbz, g);
        th = _acb_vec_init(nbz * nbth * n2);
        aux = _acb_vec_init(nbz * 4 * (ord + 1));
        test = _acb_vec_init(nbz * nbth * n2);

        /* Sample diagonal matrix tau */
        for (j = 0; j < g; j++)
        {
            acb_siegel_randtest_reduced(tau11, state, prec, bits);
            acb_set(acb_mat_entry(tau, j, j), acb_mat_entry(tau11, 0, 0));
        }
        for (j = 0; j < nbz; j++)
        {
            acb_siegel_randtest_vec_reduced(zs + j * g, state, tau, 0, prec);
        }

        /* Call sum_jet_all at precision mprec */
        acb_theta_ctx_tau_set(ctx_tau, tau, prec);
        for (j = 0; j < nbz; j++)
        {
            acb_theta_ctx_z_set(&vec[j], zs + j * g, ctx_tau, prec);
        }
        acb_theta_sum_jet_all(th, vec, nbz, ctx_tau, ord, mprec);

        if (g == 1)
        {
            /* Call acb_modular_theta_jet */
            for (j = 0; j < nbz; j++)
            {
                acb_modular_theta_jet(aux + 3 * nbth, aux + 2 * nbth, aux,
                    aux + nbth, &zs[j], acb_mat_entry(tau, 0, 0), nbth, prec);
                _acb_vec_neg(aux + 3 * nbth, aux + 3 * nbth, nbth);
                _acb_vec_set(test + j * 4 * nbth, aux, 4 * nbth);
            }
            if (!_acb_vec_overlaps(th, test, n2 * nbz * nbth))
            {
                flint_printf("FAIL (g = 1, ord = %wd, nbz = %wd, mprec = %wd, prec = %wd)\n",
                    ord, nbz, mprec, prec);
                acb_mat_printd(tau, 5);
                _acb_vec_printd(zs, nbz * g, 5);
                flint_printf("th: ");
                _acb_vec_printd(th, nbz * n2 * nbth, 5);
                flint_printf("test: ");
                _acb_vec_printd(test, nbz * n2 * nbth, 5);
                flint_abort();
            }
        }
        else
        {
            /* Make test vector using products of derivatives wrt each variable */
            slong * tups;
            acb_theta_ctx_tau_t ctx_tau11;
            acb_theta_ctx_z_struct * vec_g1;
            ulong ab, a1b1;

            tups = flint_malloc(nbth * g * sizeof(slong));
            acb_theta_ctx_tau_init(ctx_tau11, 1);
            vec_g1 = acb_theta_ctx_z_vec_init(nbz, 1);

            acb_theta_jet_tuples(tups, ord, g);
            for (j = 0; j < n2 * nbth * nbz; j++)
            {
                acb_one(&test[j]);
            }

            for (k = 0; k < g; k++)
            {
                acb_set(acb_mat_entry(tau11, 0, 0), acb_mat_entry(tau, k, k));
                acb_theta_ctx_tau_set(ctx_tau11, tau11, prec);
                for (j = 0; j < nbz; j++)
                {
                    acb_theta_ctx_z_set(&vec_g1[j], &zs[j * g + k], ctx_tau11, prec);
                }

                acb_theta_sum_jet_all(aux, vec_g1, nbz, ctx_tau11, ord, prec);

                for (j = 0; j < nbz; j++)
                {
                    for (ab = 0; ab < n2; ab++)
                    {
                        a1b1 = 2 * ((ab >> (2 * g - k - 1)) % 2) + ((ab >> (g - k - 1)) % 2);
                        for (l = 0; l < nbth; l++)
                        {
                            acb_mul(&test[j * n2 * nbth + ab * nbth + l],
                                &test[j * n2 * nbth + ab * nbth + l],
                                &aux[j * 4 * (ord + 1) + a1b1 * (ord + 1) + tups[l * g + k]], prec);
                        }
                    }
                }
            }
            if (!_acb_vec_overlaps(th, test, n2 * nbz * nbth))
            {
                flint_printf("FAIL (g = %wd, ord = %wd, nbz = %wd, mprec = %wd, prec = %wd)\n",
                    g, ord, nbz, mprec, prec);
                acb_mat_printd(tau, 5);
                _acb_vec_printd(zs, nbz * g, 5);
                flint_printf("th: ");
                _acb_vec_printd(th, nbz * n2 * nbth, 5);
                flint_printf("test: ");
                _acb_vec_printd(test, nbz * n2 * nbth, 5);
                flint_abort();
            }

            flint_free(tups);
            acb_theta_ctx_tau_clear(ctx_tau11);
            acb_theta_ctx_z_vec_clear(vec_g1, nbz);
        }

        acb_mat_clear(tau);
        acb_mat_clear(tau11);
        _acb_vec_clear(zs, nbz * g);
        acb_theta_ctx_tau_clear(ctx_tau);
        acb_theta_ctx_z_vec_clear(vec, nbz);
        _acb_vec_clear(th, nbz * n2 * nbth);
        _acb_vec_clear(aux, nbz * 4 * (ord + 1));
        _acb_vec_clear(test, nbz * n2 * nbth);
    }

    TEST_FUNCTION_END(state);
}
