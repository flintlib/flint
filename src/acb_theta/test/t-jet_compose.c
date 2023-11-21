/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_mat.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_jet_compose, state)
{
    slong iter;

    /* Test: chain rule */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong ord = n_randint(state, 4);
        slong nb = acb_theta_jet_nb(ord, g);
        slong prec = 200;
        slong mag_bits = 2;
        acb_mat_t N1, N2, N3;
        acb_ptr v1, v2, v3, test;
        slong k;

        acb_mat_init(N1, g, g);
        acb_mat_init(N2, g, g);
        acb_mat_init(N3, g, g);
        v1 = _acb_vec_init(nb);
        v2 = _acb_vec_init(nb);
        v3 = _acb_vec_init(nb);
        test = _acb_vec_init(nb);

        for (k = 0; k < nb; k++)
        {
            acb_randtest_precise(&v3[k], state, prec, mag_bits);
        }
        acb_mat_randtest(N1, state, prec, mag_bits);
        acb_mat_randtest(N2, state, prec, mag_bits);
        acb_mat_mul(N3, N2, N1, prec);

        acb_theta_jet_compose(v2, v3, N2, ord, prec);
        acb_theta_jet_compose(v1, v2, N1, ord, prec);
        acb_theta_jet_compose(test, v3, N3, ord, prec);

        if (!_acb_vec_overlaps(test, v1, nb))
        {
            flint_printf("FAIL (g = %wd, ord = %wd)\n", g, ord);
            _acb_vec_printd(v3, nb, 5);
            _acb_vec_printd(test, nb, 5);
            flint_abort();
        }

        acb_mat_clear(N1);
        acb_mat_clear(N2);
        acb_mat_clear(N3);
        _acb_vec_clear(v1, nb);
        _acb_vec_clear(v2, nb);
        _acb_vec_clear(v3, nb);
        _acb_vec_clear(test, nb);
    }

    TEST_FUNCTION_END(state);
}
