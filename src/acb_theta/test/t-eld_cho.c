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

    flint_printf("eld_cho....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: C^T C = pi Im(tau) on good input, not finite on phony input */
    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 10);
        slong prec = 200;
        slong mag_bits = n_randint(state, 4);
        acb_mat_t tau;
        arb_mat_t C;
        arb_mat_t im, test;
        arb_t pi;

        acb_mat_init(tau, g, g);
        arb_mat_init(C, g, g);
        arb_mat_init(im, g, g);
        arb_mat_init(test, g, g);
        arb_init(pi);

        acb_siegel_randtest(tau, state, prec, mag_bits);
        acb_theta_eld_cho(C, tau, prec);
        arb_mat_transpose(test, C);
        arb_mat_mul(test, test, C, prec);
        acb_mat_get_imag(im, tau);
        arb_const_pi(pi, prec);
        arb_mat_scalar_mul_arb(im, im, pi, prec);

        if (!arb_mat_overlaps(im, test))
        {
            flint_printf("FAIL\n");
            acb_mat_printd(tau, 5);
            arb_mat_printd(C, 5);
            flint_abort();
        }

        acb_zero(acb_mat_entry(tau, 0, 0));
        acb_theta_eld_cho(C, tau, prec);

        if (arb_mat_is_finite(C))
        {
            flint_printf("FAIL (not infinite)\n");
            flint_abort();
        }

        acb_mat_clear(tau);
        arb_mat_clear(C);
        arb_mat_clear(im);
        arb_mat_clear(test);
        arb_clear(pi);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
