/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "mpn_extras.h"
#include "ulong_extras.h"

int main(void)
{
    slong iter;
    FLINT_TEST_INIT(state);

    flint_printf("sqrtrem2....");
    fflush(stdout);

    for (iter = 0; iter < 100000 * flint_test_multiplier(); iter++)
    {
        mp_limb_t X[2], S1[1], S2[1], R1[2], R2[2];
        mp_size_t L1, L2;

        X[0] = n_randtest(state);
        X[1] = n_randtest_not_zero(state);

        /* Generate perfect squares, or +/- 1 */
        if (n_randint(state, 2))
        {
            mpn_sqrtrem(S1, NULL, X, 2);
            umul_ppmm(X[1], X[0], S1[0], S1[0]);
            if (n_randint(state, 2))
                add_ssaaaa(X[1], X[0], X[1], X[0], 0, 1);
            if (n_randint(state, 2))
                sub_ddmmss(X[1], X[0], X[1], X[0], 0, 1);
            if (X[1] == 0)
                X[1] = n_randtest_not_zero(state);
        }

        L1 = mpn_sqrtrem(S1, R1, X, 2);
        L2 = flint_mpn_sqrtrem2(S2, R2, X);

        if (S1[0] != S2[0] || L1 != L2 || (L1 > 0 && R1[0] != R2[0]) || (L1 == 2 && R1[1] != R2[1]))
        {
            flint_printf("FAIL\n");
            flint_printf("X = %wu, %wu\n", X[1], X[0]);
            flint_printf("L1 = %wd, L2 = %wd\n", L1, L2);
            flint_printf("S1 = %wu, S2 = %wu\n", S1[0], S2[0]);
            flint_printf("R1[0] = %wd, R2[0] = %wu\n", R1[0], R2[0]);
            flint_printf("R1[1] = %wd, R2[1] = %wu\n", R1[1], R2[1]);
            flint_abort();
        }

        L1 = (mpn_sqrtrem(S1, NULL, X, 2) != 0);
        L2 = (flint_mpn_sqrtrem2(S2, NULL, X) != 0);

        if (S1[0] != S2[0] || L1 != L2)
        {
            flint_printf("FAIL\n");
            flint_printf("X = %wu, %wu\n", X[1], X[0]);
            flint_printf("L1 = %wd, L2 = %wd\n", L1, L2);
            flint_printf("S1 = %wu, S2 = %wu\n", S1[0], S2[0]);
            flint_abort();
        }
    }
    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
