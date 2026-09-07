/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.1

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "test_helpers.h"
#include "fmpz.h"
#include "ecpp.h"

TEST_FUNCTION_START(ecpp_cornacchia, state)
{
    slong iter;

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        fmpz_t n, D, sqrtD, t, v, u;
        slong Dv;
        int r;

        fmpz_init(n); fmpz_init(D); fmpz_init(sqrtD);
        fmpz_init(t); fmpz_init(v); fmpz_init(u);

        /* random prime and small discriminant with (D/n) = 1 */
        fmpz_randprime(n, state, 10 + n_randint(state, 300), 0);
        do
        {
            Dv = -(slong) (3 + n_randint(state, 2000));
            if ((-Dv) % 4 != 3 && (-Dv) % 4 != 0)
                continue;
            fmpz_set_si(D, Dv);
            fmpz_mod(D, D, n);
        } while (fmpz_jacobi(D, n) != 1);

        if (!fmpz_sqrtmod(sqrtD, D, n))
            flint_abort();

        r = ecpp_cornacchia(t, v, n, Dv, sqrtD);

        if (r)
        {
            /* t^2 + |D| v^2 = 4n */
            fmpz_mul(u, t, t);
            fmpz_mul(v, v, v);
            fmpz_mul_si(v, v, -Dv);
            fmpz_add(u, u, v);
            fmpz_mul_2exp(t, n, 2);
            if (!fmpz_equal(u, t))
            {
                flint_printf("FAIL: t^2 + |D| v^2 != 4n\n");
                flint_printf("n = "); fmpz_print(n); flint_printf(" D = %wd\n", Dv);
                fflush(stdout);
                flint_abort();
            }
        }
        else
        {
            /* verify there is really no solution by brute force for small n */
            if (fmpz_bits(n) <= 20)
            {
                slong nn = fmpz_get_si(n), tt, found = 0;
                for (tt = 0; tt * tt <= 4 * nn; tt++)
                {
                    slong rem = 4 * nn - tt * tt;
                    if (rem % (-Dv) == 0)
                    {
                        slong vv = rem / (-Dv), sq = (slong) sqrt((double) vv);
                        while (sq * sq < vv) sq++;
                        while (sq * sq > vv) sq--;
                        if (sq * sq == vv)
                            found = 1;
                    }
                }
                if (found)
                {
                    flint_printf("FAIL: missed a solution, n = %wd D = %wd\n", nn, Dv);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        fmpz_clear(n); fmpz_clear(D); fmpz_clear(sqrtD);
        fmpz_clear(t); fmpz_clear(v); fmpz_clear(u);
    }

    TEST_FUNCTION_END(state);
}
