/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include <mpfr.h>
#include "mpfr_mat.h"

TEST_FUNCTION_START(mpfr_mat_init_clear, state)
{
    int i;

    for (i = 0; i < 1000; i++)
    {
        mpfr_mat_t a;
        slong j, k;
        slong rows = n_randint(state, 10);
        slong cols = n_randint(state, 10);
        mp_prec_t prec = n_randint(state, 200) + MPFR_PREC_MIN;

        mpfr_mat_init(a, rows, cols, prec);

        for (j = 0; j < rows; j++)
            for (k = 0; k < cols; k++)
                mpfr_set_ui(a->rows[j] + k, 0, GMP_RNDN);

        mpfr_mat_clear(a);
    }

    TEST_FUNCTION_END(state);
}
