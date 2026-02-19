/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "test_helpers.h"
#include "radix.h"

TEST_FUNCTION_START(radix_rsqrt_1_approx, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        nn_ptr A, B, err;
        ulong a;
        slong n, errn, err_limbs = 10;
        double dbl_err, dbl_powB;
        slong j;

        radix_t radix;
        radix_init_randtest(radix, state);

        if (n_randint(state, 100))
            n = 1 + n_randint(state, 50);
        else
            n = 1 + n_randint(state, 500);

        a = n_randint(state, LIMB_RADIX(radix));

        if (a > 1)
        {
            A = flint_malloc(sizeof(ulong) * (n + err_limbs));
            B = flint_malloc(sizeof(ulong) * (n + err_limbs));
            err = flint_malloc(sizeof(ulong) * (n + err_limbs));

            radix_randtest_limbs(A, state, n, radix);
            radix_randtest_limbs(B, state, n, radix);

            radix_rsqrt_1_approx_basecase(A, a, n + err_limbs, radix);
            flint_mpn_zero(B, n + err_limbs);
            flint_mpn_zero(B, err_limbs);
            radix_rsqrt_1_approx(B + err_limbs, a, n, radix);

            if (mpn_cmp(A, B, n + err_limbs) >= 0)
                radix_sub(err, A, n + err_limbs, B, n + err_limbs, radix);
            else
                radix_sub(err, B, n + err_limbs, B, n + err_limbs, radix);

            errn = n + err_limbs;
            MPN_NORM(err, errn);

            dbl_err = 0.0;
            dbl_powB = pow(LIMB_RADIX(radix), -err_limbs);
            for (j = 0; j < errn; j++)
            {
                dbl_err += dbl_powB * err[j];
                dbl_powB *= LIMB_RADIX(radix);
            }

            if (dbl_err > 1.7)
            {
                flint_printf("FAIL: too large error: radix %wu^%wd, n = %wd\n", radix->B.n, radix->exp, n);
                flint_printf("a = %wu\n", a);
                flint_printf("A = %{ulong*}\n", A, n + err_limbs);
                flint_printf("B = %{ulong*}\n", B, n + err_limbs);
                flint_printf("err = %{ulong*}\n", err, n + err_limbs);
                flint_printf("dbl_err: %f\n", dbl_err);
                flint_abort();
            }

            flint_free(A);
            flint_free(B);
            flint_free(err);
        }

        radix_clear(radix);
    }

    TEST_FUNCTION_END(state);
}
