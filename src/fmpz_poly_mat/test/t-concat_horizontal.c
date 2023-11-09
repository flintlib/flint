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
#include "fmpz_poly_mat.h"

TEST_FUNCTION_START(fmpz_poly_mat_concat_horizontal, state)
{
    fmpz_poly_mat_t A, B, C;
    fmpz_poly_mat_t window1, window2;
    slong i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        slong c1, c2, r1, bits;

        c1 = n_randint(state, 10);
        c2 = n_randint(state, 10);
        r1 = n_randint(state, 10);
        bits = 1 + n_randint(state, 20);

        fmpz_poly_mat_init(A, r1, c1);
        fmpz_poly_mat_init(B, r1, c2);
        fmpz_poly_mat_init(C, r1, (c1 + c2));

        fmpz_poly_mat_randtest(A, state, n_randint(state, 10) + 1, bits);
        fmpz_poly_mat_randtest(B, state, n_randint(state, 10) + 1, bits);

        fmpz_poly_mat_randtest(C, state, n_randint(state, 10) + 1, bits);

        fmpz_poly_mat_concat_horizontal(C, A, B);

        fmpz_poly_mat_window_init(window1, C, 0, 0, r1, c1);
        fmpz_poly_mat_window_init(window2, C, 0, c1, r1, (c1 + c2));

        if (!(fmpz_poly_mat_equal(window1, A) && fmpz_poly_mat_equal(window2, B)))
        {
            flint_printf("FAIL: results not equal\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_mat_clear(A);
        fmpz_poly_mat_clear(B);
        fmpz_poly_mat_clear(C);

        fmpz_poly_mat_window_clear(window1);
        fmpz_poly_mat_window_clear(window2);
    }

    TEST_FUNCTION_END(state);
}
