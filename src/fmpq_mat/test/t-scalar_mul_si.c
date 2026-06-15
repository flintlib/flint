/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2026 Lars Göttgens

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "long_extras.h"
#include "fmpq_mat.h"

TEST_FUNCTION_START(fmpq_mat_scalar_mul_si, state)
{
    int i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t C, B, A;
        slong rows, cols;
        slong c;

        rows = n_randint(state, 10);
        cols = n_randint(state, 10);

        fmpq_mat_init(A, rows, cols);
        fmpq_mat_init(B, rows, cols);
        fmpq_mat_init(C, rows, cols);

        c = z_randtest(state);
        fmpq_mat_randtest(A, state, 100);

        fmpq_mat_scalar_mul_si(B, A, c);

        if (c == 0)
        {
            if (!fmpq_mat_is_zero(B))
            {
                flint_printf("FAIL!\n");
                fflush(stdout);
                flint_abort();
            }
        }
        else
        {
            fmpq_mat_scalar_div_si(C, B, c);

            if (!fmpq_mat_equal(C, A))
            {
                flint_printf("FAIL!\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
        fmpq_mat_clear(C);
    }

    TEST_FUNCTION_END(state);
}
