/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "d_vec.h"
#include "ulong_extras.h"

#ifdef __GNUC__
# define fabs __builtin_fabs
#else
# include <math.h>
#endif

#define D_VEC_NORM_EPS (1e-14)

TEST_FUNCTION_START(d_vec_norm, state)
{
    int i, result;

    /* check sum of norms of parts of a vector is equal to the
       norm of the vector */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        double *a;
        double res1, res2, res3;
        slong len = n_randint(state, 100);
        if (!len)
            continue;

        a = _d_vec_init(len);
        _d_vec_randtest(a, state, len, 0, 0);

        res1 = _d_vec_norm(a, len - 1);
        res2 = _d_vec_norm(a + len - 1, 1);
        res3 = _d_vec_norm(a, len);

        result = fabs(res1 + res2 - res3) < D_VEC_NORM_EPS;
        if (!result)
        {
            flint_printf("FAIL:\n");
            printf("%g\n", fabs(res1 + res2 - res3));
            fflush(stdout);
            flint_abort();
        }

        _d_vec_clear(a);
    }

    TEST_FUNCTION_END(state);
}
