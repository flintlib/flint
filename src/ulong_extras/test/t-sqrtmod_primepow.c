/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_sqrtmod_primepow, state)
{
    int i, result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test random squares mod a power of 2 */
    {
        mp_limb_t a, b, p, pow, pow2, pinv;
        slong exp, num, i;
        mp_limb_t * sqrt;
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
        {
            flint_printf("FAIL:\n");
            flint_printf("p = %wu\n", p);
            flint_printf("exp = %wd\n", exp);
            flint_printf("a = %wu\n", a);
            flint_printf("b = %wu\n", b);
            flint_printf("num = %wd\n", num);

            if (!btest)
                flint_printf("Square root not found.\n");
            if (i != num)
                flint_printf("%wu not a square root of %wu mod %wu\n", sqrt[i], a, pow);

            fflush(stdout);
            flint_abort();
        }

        flint_free(sqrt);
    }

    for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test random squares mod other prime powers */
    {
        mp_limb_t a, b, p, pow, pow2, pinv;
        slong exp, maxexp, num, i;
        flint_bitcnt_t bits;
        mp_limb_t * sqrt;
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
        {
            flint_printf("FAIL:\n");
            flint_printf("p = %wu\n", p);
            flint_printf("exp = %wd\n", exp);
            flint_printf("a = %wu\n", a);
            flint_printf("b = %wu\n", b);
            flint_printf("num = %wd\n", num);

            if (!btest)
                flint_printf("Square root not found.\n");
            if (i != num)
                flint_printf("%wu not a square root of %wu mod %wu\n", sqrt[i], a, pow);

            fflush(stdout);
            flint_abort();
        }

        flint_free(sqrt);
    }

    for (i = 0; i < 500 * flint_test_multiplier(); i++) /* Test random nonsquares */
    {
        mp_limb_t a, b, p, pow, pinv;
        slong exp, maxexp;
        flint_bitcnt_t bits;
        mp_limb_t * sqrt;

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
            {
                flint_printf("FAIL:\n");
                flint_printf("%wu^2 is not %wu mod %wu\n", sqrt[0], a, pow);
                fflush(stdout);
                flint_abort();
            }

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
        {
            flint_printf("FAIL:\n");
            flint_printf("p = %wu\n", p);
            flint_printf("exp = %wd\n", exp);
            flint_printf("a = %wu\n", a);
            flint_printf("b = %wu\n", b);

            fflush(stdout);
            flint_abort();
        }

        flint_free(sqrt);
    }

    TEST_FUNCTION_END(state);
}
