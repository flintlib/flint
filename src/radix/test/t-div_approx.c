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

TEST_FUNCTION_START(radix_div_approx, state)
{
    slong iter;

    //double max_err = 0.0;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        nn_ptr A, D, B, C, err;
        slong errn, err_limbs = 10;
        double dbl_err, dbl_a, dbl_d, dbl_powB;
        slong j;

        radix_t radix;
        radix_init_randtest(radix, state);

        slong n, dn, an;

        if (n_randint(state, 100))
        {
            n = 1 + n_randint(state, 10);
            an = 1 + n_randint(state, 10);
            dn = 1 + n_randint(state, 10);
        }
        else
        {
            n = 1 + n_randint(state, 100);
            an = 1 + n_randint(state, 100);
            dn = 1 + n_randint(state, 100);
        }

        A = flint_malloc(sizeof(ulong) * an);
        D = flint_malloc(sizeof(ulong) * dn);
        B = flint_malloc(sizeof(ulong) * (n + 2 + err_limbs));
        C = flint_malloc(sizeof(ulong) * (n + 2 + err_limbs));
        err = flint_malloc(sizeof(ulong) * (n + 2 + err_limbs));

        radix_randtest_limbs(A, state, an, radix);
        if (A[an - 1] == 0)
            A[an - 1] = 1;
        radix_randtest_limbs(D, state, dn, radix);

        radix_div_approx_invmul(B, D, dn, A, an, n + err_limbs, radix);

        flint_mpn_zero(C, err_limbs);
        radix_div_approx(C + err_limbs, D, dn, A, an, n, radix);

        if (mpn_cmp(B, C, n + 2 + err_limbs) >= 0)
            radix_sub(err, B, n + 2 + err_limbs, C, n + 2 + err_limbs, radix);
        else
            radix_sub(err, C, n + 2 + err_limbs, B, n + 2 + err_limbs, radix);

        errn = n + 2 + err_limbs;
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

        dbl_d = 0.0;
        dbl_powB = 1.0 / LIMB_RADIX(radix);
        for (j = 0; j < dn; j++)
        {
            dbl_d += dbl_powB * D[dn - 1 - j];
            dbl_powB /= LIMB_RADIX(radix);
        }

        dbl_err *= dbl_a;

        //if (dbl_err > max_err)
        //{
        //    max_err = dbl_err;
        //    flint_printf("MAX %.6f    radix %wu^%wd, n = %wd\n", max_err, radix->b.n, radix->exp, n);
        //}

        if (dbl_err > 3.5)
        {
            flint_printf("FAIL: too large error: radix %wu^%wd, n = %wd\n", radix->b.n, radix->exp, n);
            flint_printf("A = %{ulong*}\n", A, an);
            flint_printf("D = %{ulong*}\n", D, dn);
            flint_printf("B = %{ulong*}\n", B, n + err_limbs + 2);
            flint_printf("C = %{ulong*}\n", C, n + err_limbs + 2);
            flint_printf("err = %{ulong*}\n", err, n + err_limbs + 2);
            flint_printf("dbl_err: %f\n", dbl_err);
            flint_printf("dbl_a: %f\n", dbl_a);
            flint_printf("dbl_d: %f\n", dbl_d);
            flint_abort();
        }

        flint_free(A);
        flint_free(B);
        flint_free(C);
        flint_free(D);
        flint_free(err);

        radix_clear(radix);
    }

    TEST_FUNCTION_END(state);
}
