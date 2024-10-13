/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_sqrtmod_primepow, state)
{
    slong ix;
    int result;

    for (ix = 0; ix < 1000 * flint_test_multiplier(); ix++) /* Test random squares mod a power of 2 */
    {
        ulong a, b, p, pow, pow2, pinv;
        slong exp, num, i;
        ulong * sqrt;
        int btest;

        p = 2;
        exp = n_randint(state, FLINT_BITS - 1) + 1;
        pow = n_pow(p, exp);
        b = n_randtest(state) % pow;

        pow2 = p;
        while (FLINT_BIT_COUNT(p*pow2) <= 12)
            pow2 *= p;

        if ((b % (p*pow2)) == 0)
        {
            b += pow2;
            b %= pow;
        }

        pinv = n_preinvert_limb(pow);
        a = n_mulmod2_preinv(b, b, pow, pinv);

        num = n_sqrtmod_primepow(&sqrt, a, p, exp);

        btest = 0;
        for (i = 0; i < num; i++)
        {
            if (a != n_mulmod2_preinv(sqrt[i], sqrt[i], pow, pinv))
                break;
            if (sqrt[i] == b)
                btest = 1;
        }

        result = btest & (i == num);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "p = %wu\n"
                    "exp = %wd\n"
                    "a = %wu\n"
                    "b = %wu\n"
                    "num = %wd\n"
                    "btest = %d\n"
                    "i != num = %d\n",
                    p, exp, a, b, num, btest, i != num);

        flint_free(sqrt);
    }

    for (ix = 0; ix < 1000 * flint_test_multiplier(); ix++) /* Test random squares mod other prime powers */
    {
        ulong a, b, p, pow, pow2, pinv;
        slong exp, maxexp, num, i;
        flint_bitcnt_t bits;
        ulong * sqrt;
        int btest;

        bits = n_randint(state, 18) + 2;
        p = n_randprime(state, bits, 0);
        maxexp = FLINT_BITS/bits;
        exp = n_randint(state, maxexp) + 1;
        pow = n_pow(p, exp);
        b = n_randtest(state) % pow;

        if (bits <= FLINT_BITS/2)
        {
            pow2 = p;
            while (FLINT_BIT_COUNT(p*pow2) <= 12)
                pow2 *= p;

            if ((b % (p*pow2)) == 0)
                b += pow2;

            b %= pow;
        }

        pinv = n_preinvert_limb(pow);
        a = n_mulmod2_preinv(b, b, pow, pinv);

        num = n_sqrtmod_primepow(&sqrt, a, p, exp);

        btest = 0;
        for (i = 0; i < num; i++)
        {
            if (a != n_mulmod2_preinv(sqrt[i], sqrt[i], pow, pinv))
                break;
            if (sqrt[i] == b)
                btest = 1;
        }

        result = btest & (i == num);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "p = %wu\n"
                    "exp = %wd\n"
                    "a = %wu\n"
                    "b = %wu\n"
                    "num = %wd\n"
                    "btest = %d\n"
                    "i != num = %d\n",
                    p, exp, a, b, num, btest, i != num);

        flint_free(sqrt);
    }

    for (ix = 0; ix < 500 * flint_test_multiplier(); ix++) /* Test random nonsquares */
    {
        ulong a, b, p, pow, pinv;
        slong exp, maxexp;
        flint_bitcnt_t bits;
        ulong * sqrt;

        bits = n_randint(state, 18) + 2;
        p = n_randprime(state, bits, 0);
        maxexp = 20/bits;
        exp = n_randint(state, maxexp) + 1 + (p == 2);
        pow = n_pow(p, exp);

        pinv = n_preinvert_limb(pow);

        a = n_randtest(state) % pow;
        while (n_sqrtmod_primepow(&sqrt, a, p, exp))
        {
            if (n_mulmod2_preinv(sqrt[0], sqrt[0], pow, pinv) != a)
                TEST_FUNCTION_FAIL("%wu^2 is not %wu mod %wu\n", sqrt[0], a, pow);

            flint_free(sqrt);
            a = n_randtest(state) % pow;
        }

        for (b = 0; b < pow; b++)
        {
            if (n_mulmod2_preinv(b, b, pow, pinv) == a)
                break;
        }

        result = (b == pow);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "p = %wu\n"
                    "exp = %wd\n"
                    "a = %wu\n"
                    "b = %wu\n",
                    p, exp, a, b);

        flint_free(sqrt);
    }

    TEST_FUNCTION_END(state);
}
