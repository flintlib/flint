/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "gmpcompat.h"
#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_ll_small_powmod, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        ulong exp[2], m[2], minv[2], y[3][2];
        mpz_t zb, ze, zm, zr, zy;
        ulong b[3], bmax;
        int i;

        exp[0] = n_randtest(state);
        exp[1] = n_randtest(state);

        while (1)
        {
            if (n_randint(state, 2))
                bmax = n_randint(state, (FLINT_BITS == 64) ? 565 : 23);
            else
                bmax = n_randtest(state);
            bmax = FLINT_MAX(bmax, 2);

            m[0] = n_randtest(state);
            m[1] = n_randtest_not_zero(state);

            if (m[1] == 1 && m[0] == 0)
                continue;

            double bound;

            bound = 6.0 * bmax * (m[0] + ldexp(1.0, FLINT_BITS) * m[1]);
            bound = bound * bound;

            if (bound < ldexp(1.0, 3 * FLINT_BITS))
                break;
        }

        b[2] = n_randint(state, bmax + 1);
        b[1] = n_randint(state, b[2] + 1);
        b[0] = n_randint(state, b[1] + 1);

        n_ll_small_preinv(minv, m);

        mpz_init(zb);
        mpz_init(ze);
        mpz_init(zm);
        mpz_init(zr);
        mpz_init(zy);

        flint_mpz_set_ui(zm, m[1]);
        mpz_mul_2exp(zm, zm, FLINT_BITS);
        flint_mpz_add_ui(zm, zm, m[0]);

        flint_mpz_set_ui(ze, exp[1]);
        mpz_mul_2exp(ze, ze, FLINT_BITS);
        flint_mpz_add_ui(ze, ze, exp[0]);

        n_ll_small_2_powmod(y[0], exp, m, minv);

        mpz_set_ui(zb, 2);
        mpz_powm(zr, zb, ze, zm);
        flint_mpz_set_ui(zy, y[0][1]);
        mpz_mul_2exp(zy, zy, FLINT_BITS);
        flint_mpz_add_ui(zy, zy, y[0][0]);

        if (mpz_cmp(zr, zy))
        {
            flint_printf("FAIL: base 2\n");
            gmp_printf("exp = %Zd\n", ze);
            gmp_printf("m = %Zd\n", zm);
            gmp_printf("zr = %Zd\n", zr);
            gmp_printf("zy = %Zd\n", zy);
            flint_abort();
        }

        n_ll_small_powmod_triple(y[0], y[1], y[2], b[0], b[1], b[2], exp, m, minv);

        for (i = 0; i < 3; i++)
        {
            mpz_set_ui(zb, b[i]);
            mpz_powm(zr, zb, ze, zm);
            flint_mpz_set_ui(zy, y[i][1]);
            mpz_mul_2exp(zy, zy, FLINT_BITS);
            flint_mpz_add_ui(zy, zy, y[i][0]);

            if (mpz_cmp(zr, zy))
            {
                flint_printf("FAIL\n");
                gmp_printf("exp = %Zd\n", ze);
                flint_printf("b = %wu\n", b[i]);
                gmp_printf("m = %Zd\n", zm);
                gmp_printf("zr = %Zd\n", zr);
                gmp_printf("zy = %Zd\n", zy);
                flint_abort();
            }
        }

        mpz_clear(zb);
        mpz_clear(ze);
        mpz_clear(zm);
        mpz_clear(zr);
        mpz_clear(zy);
    }

    TEST_FUNCTION_END(state);
}
