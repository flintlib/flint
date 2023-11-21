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

TEST_FUNCTION_START(acb_theta_g2_chi35, state)
{
    slong iter;

    /* Test: transforms like a modular form */
    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        slong g = 2;
        slong n2 = 1 << (2 * g);
        slong prec = 100 + n_randint(state, 500);
        slong mag_bits = n_randint(state, 2);
        fmpz_mat_t mat;
        acb_mat_t tau;
        acb_ptr th, z;
        acb_t r, s;

        fmpz_mat_init(mat, 2 * g, 2 * g);
        acb_mat_init(tau, g, g);
        th = _acb_vec_init(n2);
        z = _acb_vec_init(g);
        acb_init(r);
        acb_init(s);

        sp2gz_randtest(mat, state, mag_bits);
        acb_siegel_randtest_reduced(tau, state, prec, mag_bits);
        acb_theta_all(th, z, tau, 0, prec);

        acb_theta_g2_chi35(r, th, prec);
        acb_theta_transform_proj(th, mat, th, 0, prec);
        acb_theta_g2_chi35(s, th, prec);
        acb_mul_i_pow_si(s, s, -acb_theta_transform_kappa2(mat));

        if (!acb_overlaps(r, s))
        {
            flint_printf("FAIL\n");
            acb_printd(r, 10);
            flint_printf("\n");
            acb_printd(s, 10);
            flint_printf("\n");
            flint_abort();
        }

        fmpz_mat_clear(mat);
        acb_mat_clear(tau);
        _acb_vec_clear(th, n2);
        _acb_vec_clear(z, g);
        acb_clear(r);
        acb_clear(s);
    }

    TEST_FUNCTION_END(state);
}
