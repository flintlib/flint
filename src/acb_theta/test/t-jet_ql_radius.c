/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_jet_ql_radius, state)
{
    slong iter;

    /* Test: inequalities are satisfied */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong prec = 100 + n_randint(state, 500);
        slong mag_bits = n_randint(state, 10);
        slong ord = n_randint(state, 10);
        slong g = n_randint(state, 10);
        arb_t c, rho, t, u;
        arf_t eps, err;

        arb_init(c);
        arb_init(rho);
        arb_init(t);
        arb_init(u);
        arf_init(eps);
        arf_init(err);

        arb_randtest_positive(c, state, prec, mag_bits);
        arb_randtest_positive(rho, state, prec, mag_bits);

        acb_theta_jet_ql_radius(eps, err, c, rho, ord, g, prec);

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
        arb_pow_ui(t, t, ord + 1, prec);
        arb_pow_ui(u, rho, 2 * ord + 1, prec);
        arb_div(t, t, u, prec);
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
        arb_clear(u);
        arf_clear(eps);
        arf_clear(err);
    }

    TEST_FUNCTION_END(state);
}
