/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int main(void)
{
    slong iter;
    flint_rand_t state;

    flint_printf("jet_naive_all....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: values match acb_modular_theta_jet on diagonal matrices */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong mprec = ACB_THETA_LOW_PREC + n_randint(state, 100);
        slong prec = mprec + 50;
        slong bits = n_randint(state, 4);
        slong ord = n_randint(state, 4);
        slong n2 = 1 << (2 * g);
        slong nb = acb_theta_jet_nb(ord, g);
        acb_mat_t tau, tau11;
        acb_ptr z, dth, dth_g1, test;
        acb_t prod, t;
        slong* tups;
        slong k, j, l;

        acb_mat_init(tau, g, g);
        acb_mat_init(tau11, 1, 1);
        z = _acb_vec_init(g);
        dth = _acb_vec_init(nb * n2);
        dth_g1 = _acb_vec_init((ord + 1) * g * n2);
        test = _acb_vec_init(nb * n2);
        acb_init(prod);
        acb_init(t);
        tups = flint_malloc(nb * g * sizeof(slong));

        for (k = 0; k < g; k++)
        {
            acb_siegel_randtest(tau11, state, prec, bits);
            acb_set(acb_mat_entry(tau, k, k), acb_mat_entry(tau11, 0, 0));
            acb_urandom(&z[k], state, prec);
        }

        acb_theta_jet_naive_all(dth, z, tau, ord, mprec);
        for (k = 0; k < g; k++)
        {
            acb_set(acb_mat_entry(tau11, 0, 0), acb_mat_entry(tau, k, k));
            acb_theta_jet_naive_all(dth_g1 + k * (ord + 1) * n2, &z[k], tau11, ord, prec);
        }

        /* Make test vector using products of derivatives wrt each variable */
        for (j = 0; j < nb; j++)
        {
            for (k = 0; k < n2; k++)
            {
                acb_one(prod);
                for (l = 0; l < g; l++)
                {
                    acb_mul(prod, prod, &dth_g1[l * (ord + 1) * n2 + k * (ord + 1) + j], prec);
                }
                acb_set(&test[k * nb + j], prod);
            }
        }

        if (!_acb_vec_overlaps(dth, test, n2 * nb))
        {
            flint_printf("FAIL (overlap)\n");
            flint_printf("g = %wd, prec = %wd, ord = %wd\n", g, prec, ord);
            acb_mat_printd(tau, 5);
            _acb_vec_printd(z, g, 5);
            flint_printf("jet_naive_all:\n");
            _acb_vec_printd(dth, nb * n2, 5);
            flint_printf("test:\n");
            _acb_vec_printd(test, nb * n2, 5);
            flint_abort();
        }

        acb_mat_clear(tau);
        acb_mat_clear(tau11);
        _acb_vec_clear(z, g);
        _acb_vec_clear(dth, nb * n2);
        _acb_vec_clear(dth_g1, (ord + 1) * g * n2);
        _acb_vec_clear(test, nb * n2);
        acb_clear(prod);
        acb_clear(t);
        flint_free(tups);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
