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

TEST_FUNCTION_START(acb_mat_vector_mul, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        slong nrow = n_randint(state, 10);
        slong ncol = n_randint(state, 10);
        slong bits = n_randint(state, 10);
        slong prec = 100 + n_randint(state, 200);
        acb_mat_t A, B;
        acb_ptr v, res, t;
        slong k;

        acb_mat_init(A, nrow, ncol);
        acb_mat_init(B, ncol, nrow);
        v = _acb_vec_init(ncol);
        res = _acb_vec_init(nrow);
        t = _acb_vec_init(nrow);

        acb_mat_randtest(A, state, prec, bits);
        for (k = 0; k < ncol; k++)
        {
            acb_randtest_precise(&v[k], state, prec, bits);
        }

        /* Test: should be equal for transpose */
        acb_mat_vector_mul_col(res, A, v, prec);
        acb_mat_transpose(B, A);
        acb_mat_vector_mul_row(t, v, B, prec);

        if (!_acb_vec_overlaps(res, t, nrow))
        {
            flint_printf("FAIL\n");
            flint_abort();
        }

        acb_mat_clear(A);
        acb_mat_clear(B);
        _acb_vec_clear(v, ncol);
        _acb_vec_clear(res, nrow);
        _acb_vec_clear(t, nrow);
    }

    TEST_FUNCTION_END(state);
}

