/*
    Copyright (C) 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* generated using Claude Code (claude-opus-5-5) */

#include <math.h>
#include "test_helpers.h"
#include "acb.h"
#include "acb_ode.h"
#include "acb_poly.h"
#include "arb.h"
#include "mag.h"

/* Operators a·θ^r + b·x^k·θ^i + (terms below the Newton polygon), whose series
   solutions have coefficients ≈ (b/a)^(n/k) / n!^((r-i)/k) */

TEST_FUNCTION_START(acb_ode_solution_growth, state)
{
    for (slong iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, c;
        arb_t ref;
        mag_t order, base;

        acb_init(a);
        acb_init(b);
        acb_init(c);
        arb_init(ref);
        mag_init(order);
        mag_init(base);

        slong r = 1 + n_randint(state, 4);
        slong i = n_randint(state, r);
        slong k = 1 + n_randint(state, 4);
        int has_edge = n_randint(state, 8);

        /* todo: generalize generator and expose as acb_ode_randtest_something */

        acb_poly_struct * dop = _acb_poly_vec_init(r + 1);

        acb_randtest_precise(a, state, 100, 4);
        acb_randtest_precise(b, state, 100, 4);
        if (acb_is_zero(a))
            acb_one(a);
        if (acb_is_zero(b))
            acb_one(b);

        acb_poly_set_acb(dop + r, a);
        for (slong l = 0; l < r; l++)
        {
            acb_randtest_precise(c, state, 100, 4);
            acb_poly_set_coeff_acb(dop + l, 0, c);
            if (has_edge && l <= i)
                for (slong j = 1; j <= k; j++)
                {
                    acb_randtest_precise(c, state, 100, 4);
                    acb_poly_set_coeff_acb(dop + l, j, c);
                }
        }
        if (has_edge)
            acb_poly_set_coeff_acb(dop + i, k, b);

        _acb_ode_solution_growth(order, base, dop, r + 1);

        if (!has_edge)
        {
            if (!mag_is_zero(base))
            {
                flint_printf("FAIL (no edge)\n\n");
                goto fail;
            }
        }
        else
        {
            acb_div(c, b, a, 100);
            acb_abs(ref, c, 100);
            arb_root_ui(ref, ref, k, 100);
            double ref_base = arf_get_d(arb_midref(ref), ARF_RND_NEAR);

            if (fabs(mag_get_d(order) * (r - i) / k - 1) > 1e-6
                || fabs(mag_get_d(base) / ref_base - 1) > 1e-6)
            {
                flint_printf("FAIL\n\n");
                flint_printf("expected order = %wd/%wd, base = %{arb}\n",
                             k, r - i, ref);
                goto fail;
            }
        }

        _acb_poly_vec_clear(dop, r + 1);
        mag_clear(base);
        mag_clear(order);
        arb_clear(ref);
        acb_clear(c);
        acb_clear(b);
        acb_clear(a);

        continue;

fail:
        for (slong l = 0; l <= r; l++)
            flint_printf("dop[%wd] = %{acb_poly}\n", l, dop + l);
        flint_printf("order = %{mag}, base = %{mag}\n", order, base);
        flint_abort();
    }

    TEST_FUNCTION_END(state);
}
