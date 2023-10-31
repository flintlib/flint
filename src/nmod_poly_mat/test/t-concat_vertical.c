/*
    Copyright (C) 2015 Anubhav Srivastava
    Copyright (C) 2015 Elena Sergeicheva

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_poly_mat.h"

TEST_FUNCTION_START(nmod_poly_mat_concat_vertical, state)
{
    nmod_poly_mat_t A, B, C;
    nmod_poly_mat_t window1, window2;
    slong i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        slong r1, r2, c1;
        mp_limb_t mod;

        r1 = n_randint(state, 10);
        r2 = n_randint(state, 10);
        c1 = n_randint(state, 10);
        mod = n_randtest_prime(state, 0);

        nmod_poly_mat_init(A, r1, c1, mod);
        nmod_poly_mat_init(B, r2, c1, mod);
        nmod_poly_mat_init(C, r1 + r2, c1, mod);

        nmod_poly_mat_randtest(A, state, n_randint(state, 10) + 1);
        nmod_poly_mat_randtest(B, state, n_randint(state, 10) + 1);

        nmod_poly_mat_randtest(C, state, n_randint(state, 10) + 1);

        nmod_poly_mat_concat_vertical(C, A, B);

        nmod_poly_mat_window_init(window1, C, 0, 0, r1, c1);
        nmod_poly_mat_window_init(window2, C, r1, 0, r1 + r2, c1);

        if (!(nmod_poly_mat_equal(window1, A) && nmod_poly_mat_equal(window2, B)))
        {
            flint_printf("FAIL: results not equal\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(B);
        nmod_poly_mat_clear(C);

        nmod_poly_mat_window_clear(window1);
        nmod_poly_mat_window_clear(window2);
    }

    TEST_FUNCTION_END(state);
}
