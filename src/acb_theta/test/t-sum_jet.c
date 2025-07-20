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
#include "acb_mat.h"
#include "acb_modular.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_sum_jet, state)
{
    slong iter;

    /* Test: matches acb_modular_theta_jet in genus 1;
       agrees with product of dim 1 case on diagonal matrices */
    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong n = 1 << g;
        slong mprec = 100 + n_randint(state, 100);
        slong prec = mprec + 50;
        slong ord = n_randint(state, 3);
        slong bits = n_randint(state, 4);
        slong nbz = n_randint(state, 3);
        slong nbjet = acb_theta_jet_nb(ord, g);
        int all_a = n_randint(state, 2);
        int all_b = n_randint(state, 2);
        slong nbth = (all_a ? n : 1) * (all_b ? n : 1);
        acb_mat_t tau, tau11;
        acb_ptr zs;
        acb_theta_ctx_tau_t ctx_tau;
        acb_theta_ctx_z_struct * vec;
        acb_ptr th, aux, test;
        slong j, k, l;

        acb_mat_init(tau, g, g);
        acb_mat_init(tau11, 1, 1);
        zs = _acb_vec_init(nbz * g);
        acb_theta_ctx_tau_init(ctx_tau, 0, g);
        vec = acb_theta_ctx_z_vec_init(nbz, g);
        th = _acb_vec_init(nbz * nbjet * nbth);
        aux = _acb_vec_init(nbz * 4 * (ord + 1));
        test = _acb_vec_init(nbz * nbjet * nbth);

        /* Sample diagonal matrix tau */
        for (j = 0; j < g; j++)
        {
            acb_siegel_randtest_reduced(tau11, state, prec, bits);
            acb_set(acb_mat_entry(tau, j, j), acb_mat_entry(tau11, 0, 0));
        }
        acb_siegel_randtest_vec_reduced(zs, state, nbz, tau, 0, prec);

        /* Call sum_jet at precision mprec */
        acb_theta_ctx_tau_set(ctx_tau, tau, prec);
        for (j = 0; j < nbz; j++)
        {
            acb_theta_ctx_z_set(&vec[j], zs + j * g, ctx_tau, prec);
        }
        acb_theta_sum_jet(th, vec, nbz, ctx_tau, ord, all_a, all_b, mprec);

        if (g == 1)
        {
            /* Make test vector using acb_modular_theta_jet */
            for (j = 0; j < nbz; j++)
            {
                acb_modular_theta_jet(aux + 3 * nbjet, aux + 2 * nbjet, aux,
                    aux + nbjet, &zs[j], acb_mat_entry(tau, 0, 0), nbjet, prec);
                _acb_vec_neg(aux + 3 * nbjet, aux + 3 * nbjet, nbjet);
                if (all_a && !all_b)
                {
                    _acb_vec_set(test + 2 * j * nbjet, aux, nbjet);
                    _acb_vec_set(test + (2 * j + 1) * nbjet, aux + 2 * nbjet, nbjet);
                }
                else
                {
                    _acb_vec_set(test + j * nbth * nbjet, aux, nbth * nbjet);
                }
            }
        }
        else
        {
            /* Make test vector using products of derivatives wrt each variable */
            slong * tups;
            acb_theta_ctx_tau_t ctx_tau11;
            acb_theta_ctx_z_struct * vec_g1;
            ulong ab, a1b1;

            tups = flint_malloc(nbjet * g * sizeof(slong));
            acb_theta_ctx_tau_init(ctx_tau11, 0, 1);
            vec_g1 = acb_theta_ctx_z_vec_init(nbz, 1);

            acb_theta_jet_tuples(tups, ord, g);
            for (j = 0; j < nbth * nbjet * nbz; j++)
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

                acb_theta_sum_jet(aux, vec_g1, nbz, ctx_tau11, ord, 1, 1, prec);

                for (j = 0; j < nbz; j++)
                {
                    for (ab = 0; ab < nbth; ab++)
                    {
                        if (all_a && !all_b)
                        {
                            /* ab actually encodes the characteristic n * ab */
                            a1b1 = 2 * acb_theta_char_bit(ab, k, g);
                        }
                        else
                        {
                            a1b1 = 2 * acb_theta_char_bit(ab, k, 2 * g)
                                + acb_theta_char_bit(ab, g + k, 2 * g);
                        }
                        for (l = 0; l < nbjet; l++)
                        {
                            acb_mul(&test[j * nbth * nbjet + ab * nbjet + l],
                                &test[j * nbth * nbjet + ab * nbjet + l],
                                &aux[j * 4 * (ord + 1) + a1b1 * (ord + 1) + tups[l * g + k]], prec);
                        }
                    }
                }
            }

            flint_free(tups);
            acb_theta_ctx_tau_clear(ctx_tau11);
            acb_theta_ctx_z_vec_clear(vec_g1, nbz);
        }

        if (!_acb_vec_overlaps(th, test, nbth * nbz * nbjet)
            || !_acb_vec_is_finite(th, nbth * nbz * nbjet)
            || !_acb_vec_is_finite(test, nbth * nbz * nbjet))
        {
            flint_printf("FAIL (g = %wd, ord = %wd, nbz = %wd, all_a = %wd, all_b = %wd, mprec = %wd, prec = %wd)\n",
                g, ord, nbz, all_a, all_b, mprec, prec);
            acb_mat_printd(tau, 5);
            _acb_vec_printd(zs, nbz * g, 5);
            flint_printf("th: ");
            _acb_vec_printd(th, nbz * nbth * nbjet, 5);
            flint_printf("test: ");
            _acb_vec_printd(test, nbz * nbth * nbjet, 5);
            flint_abort();
        }

        acb_mat_clear(tau);
        acb_mat_clear(tau11);
        _acb_vec_clear(zs, nbz * g);
        acb_theta_ctx_tau_clear(ctx_tau);
        acb_theta_ctx_z_vec_clear(vec, nbz);
        _acb_vec_clear(th, nbz * nbth * nbjet);
        _acb_vec_clear(aux, nbz * 4 * (ord + 1));
        _acb_vec_clear(test, nbz * nbth * nbjet);
    }

    TEST_FUNCTION_END(state);
}
