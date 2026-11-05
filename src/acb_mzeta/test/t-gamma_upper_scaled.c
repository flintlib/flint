/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_hypgeom.h"
#include "acb_mzeta.h"

TEST_FUNCTION_START(acb_mzeta_gamma_upper_scaled, state)
{
    slong iter;

    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        acb_t a, r1, r2, tt;
        arb_t t;
        slong prec = 2 + n_randint(state, 400);

        acb_init(a);
        acb_init(r1);
        acb_init(r2);
        acb_init(tt);
        arb_init(t);

        /* a: real or complex, integers (incl. nonpositive) sometimes */
        if (n_randint(state, 4) == 0)
            acb_set_si(a, (slong) n_randint(state, 12) - 8);
        else
        {
            arb_set_si(acb_realref(a), (slong) n_randint(state, 4000) - 2000);
            arb_mul_2exp_si(acb_realref(a), acb_realref(a), -8);
            if (n_randint(state, 3) == 0)
            {
                arb_set_si(acb_imagref(a), (slong) n_randint(state, 400) - 200);
                arb_mul_2exp_si(acb_imagref(a), acb_imagref(a), -6);
            }
        }

        /* t in (0, ~300], sometimes an inexact ball */
        arb_set_ui(t, 1 + n_randint(state, 1 << 16));
        arb_mul_2exp_si(t, t, -(slong) n_randint(state, 16) - 4 + (slong) n_randint(state, 9));
        if (n_randint(state, 2))
        {
            arb_const_pi(acb_realref(tt), prec + 20);
            arb_mul(t, t, acb_realref(tt), prec + 20);
        }

        _acb_mzeta_gamma_upper_scaled(r1, a, t, NULL, prec);

        acb_set_arb(tt, t);
        acb_hypgeom_gamma_upper(r2, a, tt, 2, prec + 100 + 2 * arf_abs_bound_lt_2exp_si(arb_midref(t)));

        if (!acb_overlaps(r1, r2) || (acb_is_finite(r2) && acb_rel_accuracy_bits(r2) > prec + 20
            && arb_is_exact(t) && acb_rel_accuracy_bits(r1) < prec - 12))
        {
            flint_printf("FAIL\n");
            flint_printf("prec = %wd\n", prec);
            flint_printf("a = "); acb_printd(a, 20); flint_printf("\n");
            flint_printf("t = "); arb_printd(t, 20); flint_printf("\n");
            flint_printf("r1 = "); acb_printd(r1, 20); flint_printf("\n");
            flint_printf("r2 = "); acb_printd(r2, 20); flint_printf("\n");
            flint_abort();
        }

        acb_clear(a);
        acb_clear(r1);
        acb_clear(r2);
        acb_clear(tt);
        arb_clear(t);
    }

    TEST_FUNCTION_END(state);
}
