/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mat.h"
#include "acb.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_g2_even_weight, state)
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
        acb_ptr th2, mf, test;
        slong k;

        fmpz_mat_init(mat, 2 * g, 2 * g);
        th2 = _acb_vec_init(n2);
        mf = _acb_vec_init(4);
        test = _acb_vec_init(4);

        sp2gz_randtest(mat, state, mag_bits);
        for (k = 0; k < n2; k++)
        {
            acb_randtest_precise(&th2[k], state, prec, mag_bits);
        }

        acb_theta_g2_even_weight(&mf[0], &mf[1], &mf[2], &mf[3], th2, prec);
        acb_theta_char_shuffle(th2, mat, th2, 1, prec);
        acb_theta_g2_even_weight(&test[0], &test[1], &test[2], &test[3], th2, prec);
        if (acb_siegel_kappa2(mat) % 2 == 1)
        {
            /* Negate weights 2 mod 4 */
            acb_neg(&test[1], &test[1]);
            acb_neg(&test[3], &test[3]);
        }

        if (!_acb_vec_overlaps(test, mf, 4)
            || !_acb_vec_is_finite(test, 4)
            || !_acb_vec_is_finite(mf, 4))
        {
            flint_printf("FAIL\n");
            fmpz_mat_print_pretty(mat);
            _acb_vec_printd(test, 4, 5);
            _acb_vec_printd(mf, 4, 5);
            flint_abort();
        }

        fmpz_mat_clear(mat);
        _acb_vec_clear(th2, n2);
        _acb_vec_clear(mf, 4);
        _acb_vec_clear(test, 4);
    }

    TEST_FUNCTION_END(state);
}
