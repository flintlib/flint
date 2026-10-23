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

TEST_FUNCTION_START(radix_sqrt_approx, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        nn_ptr A, B, C, T, err;
        slong errn, err_limbs = 10;
        double dbl_err, dbl_a, dbl_powB;
        slong j, N;

        radix_t radix;
        radix_init_randtest(radix, state);

        slong n, an;

        if (n_randint(state, 100))
        {
            n = 1 + n_randint(state, 100);
            an = 1 + n_randint(state, 100);
        }
        else
        {
            n = 1 + n_randint(state, 700);
            an = 1 + n_randint(state, 700);
        }

        N = n + err_limbs;

        A = flint_malloc(sizeof(ulong) * an);
        B = flint_malloc(sizeof(ulong) * (N + 2));
        C = flint_malloc(sizeof(ulong) * (N + 2));
        T = flint_malloc(sizeof(ulong) * (2 * N + 6));
        err = flint_malloc(sizeof(ulong) * (N + 2));

        radix_randtest_limbs(A, state, an, radix);

        if (an >= 2 && n_randint(state, 2))
        {
            A[an - 1] = 0;
            if (A[an - 2] == 0)
                A[an - 2] = 1 + n_randint(state, LIMB_RADIX(radix) - 1);
        }
        else
        {
            if (A[an - 1] == 0)
                A[an - 1] = 1 + n_randint(state, LIMB_RADIX(radix) - 1);
        }

        /* Reference: multiply a basecase reciprocal square root at higher
           precision by A, mirroring radix_sqrt_approx_rsqrtmul. */
        {
            slong A2n = FLINT_MIN(an, N + 2);
            nn_srcptr A2 = A + an - A2n;
            radix_rsqrt_approx_basecase(B, A, an, N, radix);
            radix_mul(T, B, N + 2, A2, A2n, radix);
            flint_mpn_copyi(B, T + A2n, N + 2);
        }

        flint_mpn_zero(C, err_limbs);
        radix_sqrt_approx(C + err_limbs, A, an, n, radix);

        if (mpn_cmp(B, C, N + 2) >= 0)
            radix_sub(err, B, N + 2, C, N + 2, radix);
        else
            radix_sub(err, C, N + 2, B, N + 2, radix);

        errn = N + 2;
        MPN_NORM(err, errn);

        dbl_err = 0.0;
        dbl_powB = pow(LIMB_RADIX(radix), -err_limbs);
        for (j = 0; j < errn; j++)
        {
            dbl_err += dbl_powB * err[j];
            dbl_powB *= LIMB_RADIX(radix);
        }

        dbl_a = 0.0;
        dbl_powB = 1.0 / LIMB_RADIX(radix);
        for (j = 0; j < an; j++)
        {
            dbl_a += dbl_powB * A[an - 1 - j];
            dbl_powB /= LIMB_RADIX(radix);
        }

        /* |computed - sqrt(A)| * sqrt(A) should be a few ulps */
        dbl_err *= sqrt(dbl_a);

        if (dbl_err > 3.0)
        {
            flint_printf("FAIL: too large error: radix %wu^%wd, n = %wd, an = %wd\n", DIGIT_RADIX(radix), radix->exp, n, an);
            flint_printf("A = %{ulong*}\n", A, an);
            flint_printf("B = %{ulong*}\n", B, N + 2);
            flint_printf("C = %{ulong*}\n", C, N + 2);
            flint_printf("err = %{ulong*}\n", err, N + 2);
            flint_printf("dbl_err: %f\n", dbl_err);
            flint_abort();
        }

        flint_free(A);
        flint_free(B);
        flint_free(C);
        flint_free(T);
        flint_free(err);

        radix_clear(radix);
    }

    TEST_FUNCTION_END(state);
}
