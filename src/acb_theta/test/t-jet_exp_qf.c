/*
    Copyright (C) 2025 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_mat.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_jet_exp_qf, state)
{
    slong iter;

    /* Test: compatible with exponential of a sum, and with jet_compose */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong ord = n_randint(state, 4);
        slong nb = acb_theta_jet_nb(ord, g);
        slong prec = 200;
        slong mag_bits = 2;
        acb_ptr z, u, v1, v2, v3, test;
        acb_mat_t N1, N2, N3;
        slong k;

        z = _acb_vec_init(g);
        u = _acb_vec_init(g);
        v1 = _acb_vec_init(nb);
        v2 = _acb_vec_init(nb);
        v3 = _acb_vec_init(nb);
        test = _acb_vec_init(nb);
        acb_mat_init(N1, g, g);
        acb_mat_init(N2, g, g);
        acb_mat_init(N3, g, g);

        for (k = 0; k < g; k++)
        {
            acb_randtest_precise(&z[k], state, prec, mag_bits);
        }
        acb_mat_randtest(N1, state, prec, mag_bits);
        acb_mat_randtest(N2, state, prec, mag_bits);

        acb_mat_add(N3, N1, N2, prec);
        acb_theta_jet_exp_qf(v1, z, N1, ord, prec);
        acb_theta_jet_exp_qf(v2, z, N2, ord, prec);
        acb_theta_jet_exp_qf(v3, z, N3, ord, prec);
        acb_theta_jet_mul(test, v1, v2, ord, g, prec);

        if (!_acb_vec_overlaps(test, v3, nb)
            || !_acb_vec_is_finite(test, nb)
            || !_acb_vec_is_finite(v3, nb))
        {
            flint_printf("FAIL (addition)\n");
            flint_printf("g = %wd, ord = %wd\n", g, ord);
            _acb_vec_printd(v3, nb, 5);
            _acb_vec_printd(test, nb, 5);
            flint_abort();
        }

        acb_mat_vector_mul_col(u, N2, z, prec);
        acb_theta_jet_exp_qf(v1, u, N1, ord, prec);
        acb_theta_jet_compose(v3, v1, N2, ord, prec);

        acb_mat_transpose(N3, N2);
        acb_mat_mul(N3, N3, N1, prec);
        acb_mat_mul(N3, N3, N2, prec);
        acb_theta_jet_exp_qf(test, z, N3, ord, prec);

        if (!_acb_vec_overlaps(test, v3, nb)
            || !_acb_vec_is_finite(test, nb)
            || !_acb_vec_is_finite(v3, nb))
        {
            flint_printf("FAIL (composition)\n");
            flint_printf("g = %wd, ord = %wd\n", g, ord);
            _acb_vec_printd(v3, nb, 5);
            _acb_vec_printd(test, nb, 5);
            flint_abort();
        }

        _acb_vec_clear(z, g);
        _acb_vec_clear(u, g);
        _acb_vec_clear(v1, nb);
        _acb_vec_clear(v2, nb);
        _acb_vec_clear(v3, nb);
        _acb_vec_clear(test, nb);
        acb_mat_clear(N1);
        acb_mat_clear(N2);
        acb_mat_clear(N3);
    }

    TEST_FUNCTION_END(state);
}
