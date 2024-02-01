/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include <float.h>
#include "ulong_extras.h"
#include "double_extras.h"
#include "d_vec.h"

static int
_check_results(double res1, double res2)
{
    if (d_is_nan(res2))
        return d_is_nan(res1);
    else
#if !defined(_MSC_VER)
        return res1 == res2;
#else
        /* MSVC's ldexp doesn't give the same result as multiplying by 2^e
           for subnormals, so relax the test */
        return (res1 == res2) || fabs(res1 - res2) < 1e-300;
#endif
}

TEST_FUNCTION_START(d_mul_2exp, state)
{
    double x, res1, res2;
    int e;
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        x = d_randtest_special(state, -1024, 1024);
        e = n_randint(state, 3000) - 1500;

        res1 = d_mul_2exp(x, e);
        res2 = ldexp(x, e);

        if (!_check_results(res1, res2))
            TEST_FUNCTION_FAIL("x = %.20g\n res1 = %.20g\n res2 = %.20g\n", x, res1, res2);
    }

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        do {
            x = d_randtest_signed(state, -256, 256);
        } while (x == 0.0);
        e = n_randint(state, 512) - 256;

        res1 = d_mul_2exp_inrange2(x, e);
        res2 = ldexp(x, e);

        if (!(res1 == res2))
            TEST_FUNCTION_FAIL("x = %.20g\n res1 = %.20g\n res2 = %.20g\n", x, res1, res2);
    }

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        x = d_randtest_special(state, -1024, 1024);
        e = n_randint(state, D_MAX_NORMAL_EXPONENT - D_MIN_NORMAL_EXPONENT + 1) + D_MIN_NORMAL_EXPONENT;

        res1 = d_mul_2exp_inrange(x, e);
        res2 = ldexp(x, e);

        if (!_check_results(res1, res2))
            TEST_FUNCTION_FAIL("x = %.20g\n res1 = %.20g\n res2 = %.20g\n", x, res1, res2);
    }

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        double v[5] = { 0, 0, 0, 0, 0 };
        double r[10] = { 0, 0, 0, 0, 0 };
        slong i, n = n_randint(state, 5);

        for (i = 0; i < n; i++)
            x = d_randtest_special(state, -1024, 1024);

        e = n_randint(state, 3000) - 1500;

        _d_vec_mul_2exp(r, v, n, e);

        for (i = 0; i < n; i++)
        {
            res1 = r[i];
            res2 = ldexp(v[i], e);

            if (!_check_results(res1, res2))
                TEST_FUNCTION_FAIL("x = %.20g\n res1 = %.20g\n res2 = %.20g\n", x, res1, res2);
        }
    }

    TEST_FUNCTION_END(state);
}
