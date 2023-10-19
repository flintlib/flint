/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef __GNUC__
# define fabs __builtin_fabs
#else
# include <math.h>
#endif

#include "test_helpers.h"
#include "d_vec.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "fmpz_lll.h"

#define FMPZ_LLL_HD_EPS (1.0E-9)

TEST_FUNCTION_START(fmpz_lll_heuristic_dot, state)
{
    int i;
    fmpz_mat_t B;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        int result, expo1, expo2;
        double *v1, *v2;
        slong rows, cols, r1, r2;
        double d1, d2, d3, d4;

        rows = n_randint(state, 50) + 1;
        cols = n_randint(state, 50) + 2;
        r1 = n_randint(state, rows);
        r2 = n_randint(state, rows);

        fmpz_mat_init(B, rows, cols);
        fmpz_mat_randtest(B, state, 200);

        v1 = _d_vec_init(cols);
        v2 = _d_vec_init(cols);

        expo1 = _fmpz_vec_get_d_vec_2exp(v1, B->rows[r1], cols);
        expo2 = _fmpz_vec_get_d_vec_2exp(v2, B->rows[r2], cols);

        d1 = fmpz_lll_heuristic_dot(v1, v2, cols, B, r1, r2, expo1 + expo2);
        d2 = fmpz_lll_heuristic_dot(v1, v1, cols, B, r1, r1, expo1 + expo1);
        d3 = fmpz_lll_heuristic_dot(v2, v2, cols, B, r2, r2, expo2 + expo2);

        _d_vec_add(v2, v1, v2, cols);
        _fmpz_vec_add(B->rows[r2], B->rows[r1], B->rows[r2], cols);

        d4 = fmpz_lll_heuristic_dot(v2, v2, cols, B, r2, r2, expo2 + expo2);

        result = (fabs(d4 - d3 - d2 - 2 * d1) < FMPZ_LLL_HD_EPS);

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf
                ("expo1 = %d, expo2 = %d\nd1 = %g, d2 = %g, d3 = %g, d4 = %g\n",
                 expo1, expo2, d1, d2, d3, d4);
            flint_printf("%g\n", fabs(d4 - d3 - d2 - 2 * d1));
            fflush(stdout);
            flint_abort();
        }

        _d_vec_clear(v1);
        _d_vec_clear(v2);
        fmpz_mat_clear(B);
    }

    TEST_FUNCTION_END(state);
}
