/*
    Copyright (C) 2015 Elena Sergeicheva

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_mat.h"

TEST_FUNCTION_START(nmod_mat_concat_horizontal, state)
{
    nmod_mat_t A, B, C;
    nmod_mat_t window1, window2;
    slong i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        slong c1, c2, r1, n;

        c1 = n_randint(state, 50);
        c2 = n_randint(state, 50);
        r1 = n_randint(state, 50);
        n = n_randint(state, 50) + 1;

        nmod_mat_init(A, r1, c1, n);
        nmod_mat_init(B, r1, c2, n);
        nmod_mat_init(C, r1, c1 + c2, n);

        nmod_mat_randtest(A, state);
        nmod_mat_randtest(B, state);

        nmod_mat_randtest(C, state);

        nmod_mat_concat_horizontal(C, A, B);

        nmod_mat_window_init(window1, C, 0, 0, r1, c1);
        nmod_mat_window_init(window2, C, 0, c1, r1, c1 + c2);

        if (!(nmod_mat_equal(window1, A) && nmod_mat_equal(window2, B)))
        {
            flint_printf("A = \n");
            nmod_mat_print_pretty(A);
            flint_printf("B = \n");
            nmod_mat_print_pretty(B);
            flint_printf("A concat_horizontal B = \n");
            nmod_mat_print_pretty(C);
            flint_printf("FAIL: results not equal\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);

        nmod_mat_window_clear(window1);
        nmod_mat_window_clear(window2);
    }

    TEST_FUNCTION_END(state);
}
