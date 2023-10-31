/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpfr.h"
#include "gmpcompat.h"
#include "mpf_mat.h"

TEST_FUNCTION_START(mpf_mat_init_clear, state)
{
    int i;

    for (i = 0; i < 1000; i++)
    {
        mpf_mat_t a;
        slong j, k;
        slong rows = n_randint(state, 10);
        slong cols = n_randint(state, 10);
        mp_prec_t prec = n_randint(state, 200) + 2;

        mpf_mat_init(a, rows, cols, prec);

        for (j = 0; j < rows; j++)
            for (k = 0; k < cols; k++)
                flint_mpf_set_ui(a->rows[j] + k, 0);

        mpf_mat_clear(a);
    }

    TEST_FUNCTION_END(state);
}
