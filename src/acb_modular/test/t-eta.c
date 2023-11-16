/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "acb.h"
#include "acb_modular.h"

TEST_FUNCTION_START(acb_modular_eta, state)
{
    slong iter;

    /* Test functional equation */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t tau1, tau2, z1, z2, z3, t;
        fmpq_t arg;
        slong e0, prec0, prec1, prec2;
        psl2z_t g;

        psl2z_init(g);
        fmpq_init(arg);
        acb_init(tau1);
        acb_init(tau2);
        acb_init(z1);
        acb_init(z2);
        acb_init(z3);
        acb_init(t);

        e0 = 1 + n_randint(state, 200);
        prec0 = 2 + n_randint(state, 2000);
        prec1 = 2 + n_randint(state, 2000);
        prec2 = 2 + n_randint(state, 2000);

        acb_randtest(tau1, state, prec0, e0);
        acb_randtest(tau2, state, prec0, e0);
        acb_randtest(z1, state, prec0, e0);
        acb_randtest(z2, state, prec0, e0);

        psl2z_randtest(g, state, 1 + n_randint(state, 200));
        acb_modular_transform(tau2, g, tau1, prec0);

        acb_modular_eta(z1, tau1, prec1);
        acb_modular_eta(z2, tau2, prec2);

        /* apply transformation */
        fmpq_set_si(arg, acb_modular_epsilon_arg(g), 12);
        arb_sin_cos_pi_fmpq(acb_imagref(t), acb_realref(t), arg, prec1);
        acb_mul(z3, z1, t, prec1);

        acb_mul_fmpz(t, tau1, &g->c, prec1);
        acb_add_fmpz(t, t, &g->d, prec1);
        acb_sqrt(t, t, prec1);
        acb_mul(z3, z3, t, prec1);

        if (!acb_overlaps(z3, z2))
        {
            flint_printf("FAIL (overlap)\n");
            flint_printf("tau1 = "); acb_printd(tau1, 15); flint_printf("\n\n");
            flint_printf("tau2 = "); acb_printd(tau2, 15); flint_printf("\n\n");
            flint_printf("g = "); psl2z_print(g); flint_printf("\n\n");
            flint_printf("z1 = "); acb_printd(z1, 15); flint_printf("\n\n");
            flint_printf("z2 = "); acb_printd(z2, 15); flint_printf("\n\n");
            flint_printf("z3 = "); acb_printd(z3, 15); flint_printf("\n\n");
            flint_abort();
        }

        acb_modular_eta(tau1, tau1, prec2);

        if (!acb_overlaps(z1, tau1))
        {
            flint_printf("FAIL (aliasing)\n");
            flint_printf("tau1 = "); acb_print(tau1); flint_printf("\n\n");
            flint_printf("tau2 = "); acb_print(tau2); flint_printf("\n\n");
            flint_printf("z1 = "); acb_print(z1); flint_printf("\n\n");
            flint_printf("z2 = "); acb_print(z2); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(tau1);
        acb_clear(tau2);
        acb_clear(z1);
        acb_clear(z2);
        acb_clear(z3);
        acb_clear(t);
        psl2z_clear(g);
        fmpq_clear(arg);
    }

    /* Test special values */
    for (iter = 0; iter < 100 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t tau, z;
        arb_t t, u;
        slong prec;

        acb_init(tau);
        acb_init(z);
        arb_init(t);
        arb_init(u);

        prec = 2 + n_randint(state, 10000);

        acb_randtest(z, state, prec, 10);

        acb_onei(tau);
        acb_modular_eta(z, tau, prec);

        arb_one(t);
        arb_mul_2exp_si(t, t, -2);
        arb_gamma(t, t, prec);
        arb_const_pi(u, prec);
        arb_root_ui(u, u, 4, prec);
        arb_pow_ui(u, u, 3, prec);
        arb_div(t, t, u, prec);
        arb_mul_2exp_si(t, t, -1);

        if (!arb_overlaps(acb_realref(z), t) ||
            !arb_contains_zero(acb_imagref(z)))
        {
            flint_printf("FAIL (value 1)\n");
            flint_printf("tau = "); acb_print(tau); flint_printf("\n\n");
            flint_printf("z = "); acb_print(z); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(tau);
        acb_clear(z);
        arb_clear(t);
        arb_clear(u);
    }

    TEST_FUNCTION_END(state);
}
