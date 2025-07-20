/*
    Copyright (C) 2025 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"
#include "acb.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_jet_exp_pi_i, state)
{
    slong iter;

    /* Test: compatible with exponential of a sum */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong ord = n_randint(state, 4);
        slong nb = acb_theta_jet_nb(ord, g);
        slong prec = 200;
        slong mag_bits = 2;
        arb_ptr a1, a2, a3;
        acb_ptr v1, v2, v3, test;
        slong k;

        a1 = _arb_vec_init(g);
        a2 = _arb_vec_init(g);
        a3 = _arb_vec_init(g);
        v1 = _acb_vec_init(nb);
        v2 = _acb_vec_init(nb);
        v3 = _acb_vec_init(nb);
        test = _acb_vec_init(nb);

        for (k = 0; k < g; k++)
        {
            arb_randtest_precise(&a1[k], state, prec, mag_bits);
            arb_randtest_precise(&a2[k], state, prec, mag_bits);
            arb_add(&a3[k], &a1[k], &a2[k], prec);
        }

        acb_theta_jet_exp_pi_i(v1, a1, ord, g, prec);
        acb_theta_jet_exp_pi_i(v2, a2, ord, g, prec);
        acb_theta_jet_exp_pi_i(v3, a3, ord, g, prec);
        acb_theta_jet_mul(test, v1, v2, ord, g, prec);

        if (!_acb_vec_overlaps(test, v3, nb)
            || !_acb_vec_is_finite(test, nb)
            || !_acb_vec_is_finite(v3, nb))
        {
            flint_printf("FAIL (g = %wd, ord = %wd)\n", g, ord);
            _acb_vec_printd(v3, nb, 5);
            _acb_vec_printd(test, nb, 5);
            flint_abort();
        }

        _arb_vec_clear(a1, g);
        _arb_vec_clear(a2, g);
        _arb_vec_clear(a3, g);
        _acb_vec_clear(v1, nb);
        _acb_vec_clear(v2, nb);
        _acb_vec_clear(v3, nb);
        _acb_vec_clear(test, nb);
    }

    TEST_FUNCTION_END(state);
}
