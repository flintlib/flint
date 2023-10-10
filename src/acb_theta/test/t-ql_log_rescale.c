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

    flint_printf("ql_log_rescale....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: if z = C^t x, should find i|z|^2 */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong prec = 100;
        slong bits = 1 + n_randint(state, 4);
        acb_mat_t tau;
        arb_mat_t C;
        acb_ptr x, z;
        acb_t r, s, t;

        acb_mat_init(tau, g, g);
        x = _acb_vec_init(g);
        acb_init(r);
        acb_init(s);
        acb_init(t);

        acb_siegel_randtest_reduced(tau, state, prec, bits);
        for (k = 0; k < g; k++)
        {
            acb_urandom(&x[k], state, prec);
        }
        acb_theta_eld_cho(C, tau, prec);
        acb_mat_transpose(C, C);
        acb_mat_vector_mul_col(z, C, x, prec);

        acb_theta_ql_log_rescale(r, z, tau, prec);
        acb_dot(t, NULL, 0, x, 1, x, 1, prec);
        acb_const_pi(s, prec);
        acb_mul_onei(s, s);
        acb_mul(t, t, s, prec);

        if (!acb_overlaps(r, t))
        {
            flint_printf("FAIL\n");
            acb_printd(r, 5);
            flint_printf("\n");
            acb_printd(t, 5);
            flint_printf("\n");
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(x, g);
        acb_clear(r);
        acb_clear(s);
        acb_clear(t);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}

