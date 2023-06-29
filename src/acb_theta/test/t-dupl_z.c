/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int
main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("dupl_z....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with naive algorithm */
    for (iter = 0; iter < 5 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong n = 1 << (2 * g);
        slong prec = 200 + n_randint(state, 1000);
        slong mag_bits = n_randint(state, 2);

        acb_mat_t tau;
        acb_ptr z;
        acb_ptr th;
        acb_ptr dupl;
        acb_ptr test;
        arf_t rad;
        slong k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(2 * g);
        th = _acb_vec_init(2 * n);
        dupl = _acb_vec_init(2 * n);
        test = _acb_vec_init(2 * n);
        arf_init(rad);

        acb_siegel_randtest(tau, state, prec, mag_bits);
        arf_one(rad);
        for (k = 0; k < g; k++)
        {
            acb_urandom(&z[k], state, prec);
        }

        acb_theta_naive_all(th, z, 2, tau, prec);
        acb_theta_dupl_z(dupl, th, g, prec);
        _acb_vec_scalar_mul_2exp_si(z, z, g, 1);
        acb_theta_naive_all(test, z, 2, tau, prec);

        if (!_acb_vec_overlaps(dupl, test, 2 * n))
        {
            flint_printf("FAIL (overlap)\n");
            flint_printf("tau:\n");
            acb_mat_printd(tau, 10);
            flint_printf("z:\n");
            _acb_vec_printd(z, g, 10);
            flint_printf("Before dupl:\n");
            _acb_vec_printd(th, 2 * n, 10);
            flint_printf("Comparison:\n");
            _acb_vec_printd(test, 2 * n, 10);
            _acb_vec_printd(dupl, 2 * n, 10);
            fflush(stdout);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, 2 * g);
        _acb_vec_clear(th, 2 * n);
        _acb_vec_clear(dupl, 2 * n);
        _acb_vec_clear(test, 2 * n);
        arf_clear(rad);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
