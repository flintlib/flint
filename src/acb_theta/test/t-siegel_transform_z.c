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

TEST_FUNCTION_START(acb_theta_siegel_transform_z, state)
{
    slong iter;

    /* Test: matches siegel_transform, inverse matrix gives inverse transformation */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 6);
        slong prec = 100 + n_randint(state, 200);
        slong bits = n_randint(state, 10);
        acb_mat_t tau1, w, tau2;
        acb_ptr z1, r, z2;
        fmpz_mat_t m;

        acb_mat_init(tau1, g, g);
        acb_mat_init(w, g, g);
        acb_mat_init(tau2, g, g);
        z1 = _acb_vec_init(g);
        r = _acb_vec_init(g);
        z2 = _acb_vec_init(g);
        fmpz_mat_init(m, 2 * g, 2 * g);

        acb_siegel_randtest(tau1, state, prec, bits);
        acb_siegel_randtest_vec(z1, state, g, prec);

        sp2gz_randtest(m, state, bits);
        acb_siegel_transform_z(r, w, m, z1, tau1, prec);

        /* Test: agrees with transform */
        acb_siegel_transform(tau2, m, tau1, prec);
        if (!acb_mat_overlaps(tau2, w))
        {
            flint_printf("FAIL (transform)\n\n");
            acb_mat_printd(w, 10);
            flint_printf("\n");
            acb_mat_printd(tau2, 10);
            flint_printf("\n");
            flint_abort();
        }

        /* Test: inverse transformation */
        sp2gz_inv(m, m);
        acb_siegel_transform_z(z2, tau2, m, r, w, prec);
        if (!acb_mat_contains(tau2, tau1) || !_acb_vec_contains(z2, z1, g))
        {
            flint_printf("FAIL (inverse)\n\n");
            acb_mat_printd(tau1, 10);
            flint_printf("\n");
            acb_mat_printd(tau2, 10);
            flint_printf("\n\n");
            _acb_vec_printd(z1, g, 10);
            flint_printf("\n\n");
            _acb_vec_printd(z2, g, 10);
            flint_printf("\n");
            flint_abort();
        }

        acb_mat_clear(tau1);
        acb_mat_clear(w);
        acb_mat_clear(tau2);
        _acb_vec_clear(z1, g);
        _acb_vec_clear(r, g);
        _acb_vec_clear(z2, g);
        fmpz_mat_clear(m);
    }

    TEST_FUNCTION_END(state);
}
