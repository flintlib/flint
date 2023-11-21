/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_g2_psi4, state)
{
    slong iter;

    /* Test: is a modular form */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong g = 2;
        slong n2 = 16;
        slong prec = 100 + n_randint(state, 500);
        slong mag_bits = n_randint(state, 4);
        fmpz_mat_t mat;
        acb_ptr th2;
        acb_t r, s;
        slong k;

        fmpz_mat_init(mat, 2 * g, 2 * g);
        th2 = _acb_vec_init(n2);
        acb_init(r);
        acb_init(s);

        sp2gz_randtest(mat, state, mag_bits);
        for (k = 0; k < n2; k++)
        {
            acb_randtest_precise(&th2[k], state, prec, mag_bits);
        }

        acb_theta_g2_psi4(r, th2, prec);
        acb_theta_transform_proj(th2, mat, th2, 1, prec);
        acb_theta_g2_psi4(s, th2, prec);

        if (!acb_overlaps(r, s))
        {
            flint_printf("FAIL\n");
            fmpz_mat_print_pretty(mat);
            acb_printd(r, 10);
            flint_printf("\n");
            acb_printd(s, 10);
            flint_printf("\n");
            flint_abort();
        }

        fmpz_mat_clear(mat);
        _acb_vec_clear(th2, n2);
        acb_clear(r);
        acb_clear(s);
    }

    TEST_FUNCTION_END(state);
}
