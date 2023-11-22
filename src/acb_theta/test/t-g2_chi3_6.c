/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_poly.h"
#include "acb_mat.h"
#include "acb_theta.h"

static void
acb_theta_g2_chi8_6(acb_poly_t res, const acb_mat_t tau, slong prec)
{
    acb_ptr z, dth;
    acb_t c;
    slong k;

    dth = _acb_vec_init(3 * 16);
    z = _acb_vec_init(2);
    acb_init(c);

    acb_theta_jet_all(dth, z, tau, 1, prec);
    acb_theta_g2_chi3_6(res, dth, prec);
    for (k = 0; k < 16; k++)
    {
        acb_set(&dth[k], &dth[3 * k]);
    }
    acb_theta_g2_chi5(c, dth, prec);
    acb_poly_scalar_mul(res, res, c, prec);

    _acb_vec_clear(dth, 3 * 16);
    _acb_vec_clear(z, 2);
    acb_clear(c);
}

TEST_FUNCTION_START(acb_theta_g2_chi3_6, state)
{
    slong iter;

    /* Test: chi5 * chi3_6 transforms like a modular form */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        slong g = 2;
        slong prec = 100 + n_randint(state, 200);
        slong mag_bits = n_randint(state, 2);
        fmpz_mat_t mat;
        acb_mat_t tau, w, c, cinv;
        acb_poly_t r, s;

        fmpz_mat_init(mat, 2 * g, 2 * g);
        acb_mat_init(tau, g, g);
        acb_mat_init(w, g, g);
        acb_mat_init(c, g, g);
        acb_mat_init(cinv, g, g);
        acb_poly_init(r);
        acb_poly_init(s);

        sp2gz_randtest(mat, state, mag_bits);
        acb_siegel_randtest_reduced(tau, state, prec, mag_bits);

        acb_theta_g2_chi8_6(r, tau, prec);
        acb_siegel_transform_cocycle_inv(w, c, cinv, mat, tau, prec);
        acb_theta_g2_chi8_6(s, w, prec);
        acb_theta_g2_detk_symj(s, cinv, s, 8, 6, prec);

        if (!acb_poly_overlaps(r, s))
        {
            flint_printf("FAIL\n");
            acb_mat_printd(tau, 5);
            fmpz_mat_print_pretty(mat);
            flint_printf("values at tau, m*tau:\n");
            acb_poly_printd(r, 10);
            flint_printf("\n");
            acb_poly_printd(s, 10);
            flint_printf("\n");
            flint_abort();
        }

        fmpz_mat_clear(mat);
        acb_mat_clear(tau);
        acb_mat_clear(w);
        acb_mat_clear(c);
        acb_mat_clear(cinv);
        acb_poly_clear(r);
        acb_poly_clear(s);
    }

    TEST_FUNCTION_END(state);
}

