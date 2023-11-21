/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_transform_proj, state)
{
    slong iter;

    /* Test: inverse matrix gives back the same projective point */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong n2 = 1 << (2 * g);
        slong prec = 100;
        slong bits = n_randint(state, 5);
        int sqr = iter % 2;
        fmpz_mat_t mat, inv;
        acb_ptr th, aux, test;
        acb_t scal;
        slong k;

        fmpz_mat_init(mat, 2 * g, 2 * g);
        fmpz_mat_init(inv, 2 * g, 2 * g);
        th = _acb_vec_init(n2);
        aux = _acb_vec_init(n2);
        test = _acb_vec_init(n2);
        acb_init(scal);

        sp2gz_randtest(mat, state, bits);
        sp2gz_inv(inv, mat);
        for (k = 0; k < n2; k++)
        {
            acb_urandom(&test[k], state, prec);
        }

        acb_theta_transform_proj(aux, mat, test, sqr, prec);
        acb_theta_transform_proj(th, inv, aux, sqr, prec);
        acb_div(scal, &test[0], &th[0], prec);
        _acb_vec_scalar_mul(th, th, n2, scal, prec);

        if (!_acb_vec_overlaps(th, test, n2))
        {
            flint_printf("FAIL (sqr = %wd)\n", sqr);
            flint_printf("test, th:\n");
            _acb_vec_printd(test, n2, 5);
            _acb_vec_printd(th, n2, 5);
            flint_abort();
        }

        fmpz_mat_clear(mat);
        fmpz_mat_clear(inv);
        _acb_vec_clear(th, n2);
        _acb_vec_clear(aux, n2);
        _acb_vec_clear(test, n2);
        acb_clear(scal);
    }

    TEST_FUNCTION_END(state);
}
