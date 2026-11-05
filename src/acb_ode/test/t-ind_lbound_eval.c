/*
    Copyright (C) 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* generated using Claude Code (claude-opus-5-5) */

#include "test_helpers.h"
#include "acb.h"
#include "acb_ode.h"
#include "acb_poly.h"
#include "arb.h"
#include "mag.h"

TEST_FUNCTION_START(acb_ode_ind_lbound_eval, state)
{
    for (slong iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        acb_ode_exponents_t expos;
        acb_ode_ind_lbound_t ind_lbound;
        acb_poly_t ind;
        acb_t e, val;
        arb_t ref, lb, kpow;
        mag_t res;

        acb_ode_exponents_init(expos);
        acb_ode_ind_lbound_init(ind_lbound);
        acb_poly_init(ind);
        acb_init(e);
        acb_init(val);
        arb_init(ref);
        arb_init(lb);
        arb_init(kpow);
        mag_init(res);

        slong prec = 2 + n_randint(state, 100);
        slong len = 1 + n_randint(state, 8);
        slong disp = n_randint(state, 50);
        acb_ode_exponents_randtest(expos, state, len, disp, prec,
                                   1 + n_randint(state, 6));
        acb_ode_indicial_polynomial_from_exponents(ind, expos, prec);

        slong grp = n_randint(state, expos->ngroups);
        acb_ode_group_struct * group = expos->groups + grp;

        if (n_randint(state, 2))  /* overwrite an existing object */
            acb_ode_ind_lbound_precompute(ind_lbound, expos,
                                          n_randint(state, expos->ngroups),
                                          prec);
        acb_ode_ind_lbound_precompute(ind_lbound, expos, grp, prec);

        slong n = n_randint(state, 4) ? n_randint(state, 2 * disp + 10)
                                      : n_randint(state, 1000);
        acb_ode_ind_lbound_eval(res, ind_lbound, n, prec);
        arf_set_mag(arb_midref(lb), res);

        for (slong i = 0; i < 40; i++)
        {
            slong k = (i < 30) ? n + i : n + n_randint(state, 1000000);

            if (acb_ode_group_multiplicity(group, k) > 0)
                continue;

            acb_add_si(e, group->leader, k, prec);
            acb_poly_evaluate(val, ind, e, prec);
            acb_abs(ref, val, prec);
            arb_ui_pow_ui(kpow, k, len, prec);
            arb_div(ref, ref, kpow, prec);

            if (arb_lt(ref, lb))
            {
                flint_printf("FAIL\n\n");
                flint_printf("expos = ");
                acb_ode_exponents_println(expos);
                flint_printf("grp = %wd, n = %wd, k = %wd\n", grp, n, k);
                flint_printf("res = %{mag}\n", res);
                flint_printf("|q(k)|/k^r = %{arb}\n", ref);
                flint_abort();
            }
        }

        mag_clear(res);
        arb_clear(kpow);
        arb_clear(lb);
        arb_clear(ref);
        acb_clear(val);
        acb_clear(e);
        acb_poly_clear(ind);
        acb_ode_ind_lbound_clear(ind_lbound);
        acb_ode_exponents_clear(expos);
    }

    TEST_FUNCTION_END(state);
}
