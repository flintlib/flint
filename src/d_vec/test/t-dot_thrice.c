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

#define D_VEC_SP_EPS (1e-14)

TEST_FUNCTION_START(d_vec_dot_thrice, state)
{
    int i, result;

    /* check sum of scalar products of parts of vectors is equal to the
       scalar product of vectors */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        double *a, *b;
        double res1, res2, res3, err1, err2, err3;
        slong len = n_randint(state, 100);
        if (!len)
            continue;

        a = _d_vec_init(len);
        b = _d_vec_init(len);
        _d_vec_randtest(a, state, len, 0, 0);
        _d_vec_randtest(b, state, len, 0, 0);

        res1 = _d_vec_dot_thrice(a, b, len - 1, &err1);
        res2 = _d_vec_dot_thrice(a + len - 1, b + len - 1, 1, &err2);
        res3 = _d_vec_dot_thrice(a, b, len, &err3);

        result = fabs(res1 + res2 - res3) < D_VEC_SP_EPS;

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("%g\n", fabs(res1 + res2 - res3));
            flint_printf("%g\n", res1);
            flint_printf("%g\n", res2);
            flint_printf("%g\n", res3);
            flint_printf("%g\n", err1);
            flint_printf("%g\n", err2);
            flint_printf("%g\n", err3);
            fflush(stdout);
            flint_abort();
        }

        _d_vec_clear(a);
        _d_vec_clear(b);
    }

    TEST_FUNCTION_END(state);
}
