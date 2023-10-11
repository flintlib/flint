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

    flint_printf("jet_fd_radius....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: inequalities are satisfied */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong prec = 100 + n_randint(state, 500);
        slong mag_bits = n_randint(state, 10);
        slong ord = n_randint(state, 10);
        slong g = n_randint(state, 10);
        arb_t c, rho, t;
        arf_t eps, err;

        arb_init(c);
        arb_init(rho);
        arb_init(t);
        arf_init(eps);
        arf_init(err);

        arb_randtest_positive(c, state, prec, mag_bits);
        arb_randtest_positive(c, state, prec, mag_bits);

        acb_theta_jet_fd_radius(eps, err, c, rho, ord, g, prec);

        arb_set_si(t, 2 * g);
        arb_root_ui(t, t, ord + 1, prec);
        arb_mul_arf(t, t, eps, prec);

        if (arb_gt(t, rho))
        {
            flint_printf("FAIL (1st bound)\n");
            flint_printf("c, rho, eps, err:\n");
            arb_printd(c, 10);
            flint_printf("\n");
            arb_printd(rho, 10);
            flint_printf("\n");
            arf_printd(eps, 10);
            flint_printf("\n");
            flint_abort();
        }

        arb_set_arf(t, eps);
        arb_div(t, t, rho, prec);
        arb_pow_ui(t, t, ord + 1, prec);
        arb_mul(t, t, c, prec);
        arb_mul_si(t, t, 2 * g, prec);
        arb_sub_arf(t, t, err, prec);

        if (arb_is_positive(t))
        {
            flint_printf("FAIL (2nd bound)\n");
            flint_printf("c, rho, eps, err:\n");
            arb_printd(c, 10);
            flint_printf("\n");
            arb_printd(rho, 10);
            flint_printf("\n");
            arf_printd(eps, 10);
            flint_printf("\n");
            arf_printd(err, 10);
            flint_abort();
        }

        arb_clear(c);
        arb_clear(rho);
        arb_clear(t);
        arf_clear(eps);
        arf_clear(err);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
