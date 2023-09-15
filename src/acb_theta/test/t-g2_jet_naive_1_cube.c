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

    flint_printf("g2_jet_naive_1_cube....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with jet_naive_1 */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 2;
        slong n = 1 << (2 * g);
        slong p = 1 + n_randint(state, 10);
        slong prec = 100 + n_randint(state, 100);
        slong bits = n_randint(state, 4);
        acb_mat_t tau;
        acb_mat_t w;
        acb_ptr dth, test;
        slong k;

        acb_mat_init(tau, g, g);
        acb_mat_init(w, g, g);
        dth = _acb_vec_init(n * 3 * n_pow(p, 3));
        test = _acb_vec_init(n * 3);

        acb_siegel_randtest_reduced(tau, state, prec, bits);

        acb_theta_g2_jet_naive_1_cube(dth, tau, p, prec);
        
                acb_mat_printd(tau, 5);

        for (k = 0; k < n_pow(p, 3); k++)
        {
            acb_set_si(acb_mat_entry(w, 0, 0), k % p);
            acb_set_si(acb_mat_entry(w, 0, 1), (k / p) % p);
            acb_set_si(acb_mat_entry(w, 1, 0), (k / p) % p);
            acb_set_si(acb_mat_entry(w, 1, 1), (k / n_pow(p, 2)) % p);
            acb_mat_scalar_div_si(w, w, p, prec);
            acb_mat_add(w, w, tau, prec);

            acb_theta_g2_jet_naive_1(test, w, prec);
            if (!_acb_vec_overlaps(test, dth + 3 * n * k, 3 * n))
            {
                flint_printf("FAIL (p = %wd, k = %wd)\n", p, k);
                acb_mat_printd(tau, 5);
                flint_printf("values:\n");
                _acb_vec_printd(test, 3 * n, 5);
                _acb_vec_printd(dth + 3 * n * k, 3 * n, 5);
                flint_abort();
            }
        }

        acb_mat_clear(tau);
        acb_mat_clear(w);
        _acb_vec_clear(dth, n * 3 * n_pow(p, 3));
        _acb_vec_clear(test, n * 3);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
