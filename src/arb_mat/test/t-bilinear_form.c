/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

int main(void)
{
    slong iter;
    flint_rand_t state;

    flint_printf("bilinear_form....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        slong nrow = n_randint(state, 10);
        slong ncol = n_randint(state, 10);
        slong bits = n_randint(state, 10);
        slong prec = 100 + n_randint(state, 200);
        arb_mat_t A, B;
        arb_ptr v1, v2;
        arb_t x, t;
        slong k;

        arb_mat_init(A, nrow, ncol);
        arb_mat_init(B, ncol, nrow);
        v1 = _arb_vec_init(nrow);
        v2 = _arb_vec_init(ncol);
        arb_init(x);
        arb_init(t);

        arb_mat_randtest(A, state, prec, bits);
        for (k = 0; k < nrow; k++)
        {
            arb_randtest_precise(&v1[k], state, prec, bits);
        }
        for (k = 0; k < ncol; k++)
        {
            arb_randtest_precise(&v2[k], state, prec, bits);
        }

        /* Test: should be equal for transpose */
        arb_mat_bilinear_form(x, A, v1, v2, prec);
        arb_mat_transpose(B, A);
        arb_mat_bilinear_form(t, B, v2, v1, prec);

        if (!arb_overlaps(x, t))
        {
            flint_printf("FAIL\n");
            flint_abort();
        }

        arb_mat_clear(A);
        arb_mat_clear(B);
        _arb_vec_clear(v1, nrow);
        _arb_vec_clear(v2, ncol);
        arb_clear(x);
        arb_clear(t);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
