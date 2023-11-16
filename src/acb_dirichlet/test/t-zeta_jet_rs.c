/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_poly.h"
#include "acb_dirichlet.h"

TEST_FUNCTION_START(acb_dirichlet_zeta_jet_rs, state)
{
    slong iter;

    for (iter = 0; iter < 200 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t s1, s2, a;
        acb_ptr r1, r2;
        slong prec1, prec2, len;

        acb_init(a);
        acb_init(s1);
        acb_init(s2);
        r1 = _acb_vec_init(2);
        r2 = _acb_vec_init(2);

        len = 1 + n_randint(state, 2);

        if (n_randint(state, 2))
            arb_set_d(acb_realref(s1), 0.5);
        else
            arb_randtest(acb_realref(s1), state, 2 + n_randint(state, 100), 2);

        arb_randtest(acb_imagref(s1), state, 2 + n_randint(state, 100), 2);
        arb_add_ui(acb_imagref(s1), acb_imagref(s1), n_randtest(state) % 1000, 100);
        acb_set(s2, s1);

        acb_randtest(r1, state, 2 + n_randint(state, 100), 3);
        acb_randtest(r2, state, 2 + n_randint(state, 100), 3);

        if (len > 1)
        {
            acb_randtest(r1 + 1, state, 2 + n_randint(state, 100), 3);
            acb_randtest(r2 + 1, state, 2 + n_randint(state, 100), 3);
        }

        prec1 = 2 + n_randint(state, 150);
        prec2 = 2 + n_randint(state, 150);

        if (n_randint(state, 4) == 0)
        {
            mag_zero(arb_radref(acb_realref(s2)));
            mag_zero(arb_radref(acb_imagref(s2)));
        }

        if (n_randint(state, 4) == 0)
        {
            if (n_randint(state, 2))
                arb_get_ubound_arf(arb_midref(acb_realref(s2)), acb_realref(s2), ARF_PREC_EXACT);
            else
                arb_get_lbound_arf(arb_midref(acb_realref(s2)), acb_realref(s2), ARF_PREC_EXACT);
            mag_zero(arb_radref(acb_realref(s2)));
        }

        if (n_randint(state, 4) == 0)
        {
            if (n_randint(state, 2))
                arb_get_ubound_arf(arb_midref(acb_imagref(s2)), acb_imagref(s2), ARF_PREC_EXACT);
            else
                arb_get_lbound_arf(arb_midref(acb_imagref(s2)), acb_imagref(s2), ARF_PREC_EXACT);
            mag_zero(arb_radref(acb_imagref(s2)));
        }

        acb_dirichlet_zeta_jet_rs(r1, s1, len, prec1);
        acb_one(a);

        _acb_poly_zeta_cpx_series(r2, s2, a, 0, len, prec2);

        if (!acb_overlaps(r1, r2) || !acb_overlaps(r1 + 1, r2 + 1))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("iter = %wd, len = %wd\n", iter, len);
            flint_printf("s1 = "); acb_printn(s1, 50, 0); flint_printf("\n\n");
            flint_printf("s2 = "); acb_printn(s2, 50, 0); flint_printf("\n\n");
            flint_printf("r1[0] = "); acb_printn(r1, 50, 0); flint_printf("\n\n");
            flint_printf("r2[0] = "); acb_printn(r2, 50, 0); flint_printf("\n\n");
            flint_printf("r1[1] = "); acb_printn(r1 + 1, 50, 0); flint_printf("\n\n");
            flint_printf("r2[1] = "); acb_printn(r2 + 1, 50, 0); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(a);
        acb_clear(s1);
        acb_clear(s2);
        _acb_vec_clear(r1, 2);
        _acb_vec_clear(r2, 2);
    }

    TEST_FUNCTION_END(state);
}
