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

    flint_printf("g2_jet_naive_1....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with usual jet_naive at the right indices, up to pi*i */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 2;
        slong n = 1 << (2 * g);
        slong nb = acb_theta_jet_nb(1, g + 1);
        slong prec = 100 + n_randint(state, 100);
        slong bits = n_randint(state, 4);
        acb_mat_t tau;
        acb_ptr z, dth, test;
        acb_t c;
        slong k;
        int res;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        dth = _acb_vec_init(n * nb);
        test = _acb_vec_init(n * nb);
        acb_init(c);

        acb_siegel_randtest_reduced(tau, state, prec, bits);
        acb_const_pi(c, prec);
        acb_mul_onei(c, c);

        acb_theta_g2_jet_naive_1(dth, tau, prec);
        acb_theta_jet_naive_all(test, z, tau, 1, prec);

        for (k = 0; k < n; k++)
        {
            if (acb_theta_char_is_even(k, 2))
            {
                res = acb_overlaps(&dth[3 * k], &test[3 * k]);
            }
            else
            {
                _acb_vec_scalar_mul(&dth[3 * k + 1], &dth[3 * k + 1], 2, c, prec);
                res = _acb_vec_overlaps(&dth[3 * k + 1], &test[3 * k + 1], 2);
            }

            if (!res)
            {
                flint_printf("FAIL (k = %wd)\n", k);
                acb_mat_printd(tau, 5);
                flint_printf("values:\n");
                _acb_vec_printd(&dth[3 * k], 3, 5);
                _acb_vec_printd(&test[3 * k], 3, 5);
                flint_abort();
            }
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _acb_vec_clear(dth, n * nb);
        _acb_vec_clear(test, n * nb);
        acb_clear(c);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
