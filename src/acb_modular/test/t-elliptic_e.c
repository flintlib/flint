/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_modular.h"

TEST_FUNCTION_START(acb_modular_elliptic_e, state)
{
    slong iter;

    for (iter = 0; iter < 500 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t m, w, K, Kp, E, Ep, r, pi2;
        slong prec;

        acb_init(m);
        acb_init(w);
        acb_init(K);
        acb_init(Kp);
        acb_init(E);
        acb_init(Ep);
        acb_init(r);
        acb_init(pi2);

        prec = 2 + n_randint(state, 1000);
        acb_randtest(m, state, prec, 1 + n_randint(state, 100));

        acb_sub_ui(w, m, 1, prec);
        acb_neg(w, w);
        acb_const_pi(pi2, prec);
        acb_mul_2exp_si(pi2, pi2, -1);

        acb_modular_elliptic_k(K, m, prec);
        acb_modular_elliptic_k(Kp, w, prec);
        acb_modular_elliptic_e(E, m, prec);
        acb_modular_elliptic_e(Ep, w, prec);

        acb_mul(r, K, Ep, prec);
        acb_addmul(r, E, Kp, prec);
        acb_submul(r, K, Kp, prec);

        if (!acb_overlaps(r, pi2))
        {
            flint_printf("FAIL (overlap)\n\n");

            flint_printf("m = "); acb_printd(m, 30); flint_printf("\n\n");
            flint_printf("w = "); acb_printd(w, 30); flint_printf("\n\n");
            flint_printf("K = "); acb_printd(K, 30); flint_printf("\n\n");
            flint_printf("Kp = "); acb_printd(Kp, 30); flint_printf("\n\n");
            flint_printf("E = "); acb_printd(E, 30); flint_printf("\n\n");
            flint_printf("Ep = "); acb_printd(Ep, 30); flint_printf("\n\n");
            flint_printf("r = "); acb_printd(r, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(m);
        acb_clear(w);
        acb_clear(K);
        acb_clear(Kp);
        acb_clear(E);
        acb_clear(Ep);
        acb_clear(r);
        acb_clear(pi2);
    }

    TEST_FUNCTION_END(state);
}
