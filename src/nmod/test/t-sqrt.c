/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod.h"

/* A random prime; with some probability, one with a large 2-adic valuation
   of n - 1 or one near a boundary between the internal algorithm variants. */
static ulong
_randtest_sqrt_prime(flint_rand_t state)
{
    ulong n, c;
    slong s, bits;

    switch (n_randint(state, 6))
    {
        case 0:
        case 1:
        case 2:
            return n_randtest_prime(state, 0);
        case 3:
            /* n = c 2^s + 1 */
            bits = 2 + n_randint(state, FLINT_BITS - 1);
            s = 1 + n_randint(state, bits - 1);
            for (;;)
            {
                c = n_randlimb(state) | 1;
                if (bits - s < FLINT_BITS)
                    c &= (UWORD(1) << (bits - s)) - 1;
                if (c == 0)
                    continue;
                n = (c << s) + 1;
                if (n > 2 && n_is_prime(n))
                    return n;
                if (n_randint(state, 50) == 0)
                    return n_randtest_prime(state, 0);
            }
        case 4:
        {
            /* just below 2^b for b at a variant boundary */
            static const int b[] = { 30, 31, 32, 33, 62, 63, FLINT_BITS };
            bits = b[n_randint(state, 7)];
            if (bits > FLINT_BITS)
                bits = FLINT_BITS;
            n = (bits == FLINT_BITS) ? UWORD_MAX : (UWORD(1) << bits) - 1;
            n -= n_randint(state, 1000);
            n |= 1;
            while (!n_is_prime(n))
                n -= 2;
            return n;
        }
        default:
        {
            /* small primes, including the Fermat primes */
            static const ulong f[] = { 3, 5, 17, 257, 65537 };
            if (n_randint(state, 2))
                return f[n_randint(state, 5)];
            return n_nth_prime(1 + n_randint(state, 200));
        }
    }
}

TEST_FUNCTION_START(nmod_sqrt, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * flint_test_multiplier(); iter++)
    {
        nmod_t mod;
        ulong n, a, b, r;
        int success, square;

        n = _randtest_sqrt_prime(state);
        nmod_init(&mod, n);

        if (n_randint(state, 2))
        {
            /* a random square */
            b = n_randtest(state) % n;
            a = nmod_mul(b, b, mod);
            square = 1;
        }
        else
        {
            a = n_randtest(state) % n;
            square = (a == 0) || (n == 2) || n_jacobi_unsigned(a, n) == 1;
        }

        success = nmod_sqrt(&r, a, mod);

        if (success != square || (success && nmod_mul(r, r, mod) != a)
            || (success && r > n / 2) || (!success && r != 0))
        {
            TEST_FUNCTION_FAIL("n = %wu, a = %wu, square = %d, success = %d, r = %wu\n",
                n, a, square, success, r);
        }

        /* n_sqrtmod agrees */
        if (n_sqrtmod(a, n) != r)
            TEST_FUNCTION_FAIL("n_sqrtmod: n = %wu, a = %wu, r = %wu, %wu\n",
                n, a, r, n_sqrtmod(a, n));
    }

    /* primes with large 2-adic valuation of n - 1, exercising all
       variants of Tonelli-Shanks as well as Mueller's algorithm */
    {
        static const ulong special[] = {
            17, 257, 65537, UWORD(998244353), UWORD(2013265921),
            UWORD(3221225473),
#if FLINT_BITS == 64
            UWORD(0xFFFFFFFF00000001), UWORD(4179340454199820289),
            UWORD(0x3fffc00000000001), UWORD(2394789101854261249),
            UWORD(4619848792751996929), UWORD(9224603489877884929),
#endif
            0 };
        slong j;

        for (j = 0; special[j] != 0; j++)
        {
            nmod_t mod;
            ulong n, a, b, r;
            int success, square;

            n = special[j];
            if (!n_is_prime(n))
                continue;
            nmod_init(&mod, n);

            for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
            {
                if (iter % 2)
                {
                    b = n_randtest(state) % n;
                    a = nmod_mul(b, b, mod);
                }
                else
                {
                    a = n_randtest(state) % n;
                }

                square = (a == 0) || n_jacobi_unsigned(a, n) == 1;
                success = nmod_sqrt(&r, a, mod);

                if (success != square || (success && nmod_mul(r, r, mod) != a)
                    || (success && r > n / 2))
                    TEST_FUNCTION_FAIL("n = %wu, a = %wu, square = %d, success = %d, r = %wu\n",
                        n, a, square, success, r);
            }
        }
    }

    /* exhaustive for small primes */
    for (iter = 1; iter <= 60; iter++)
    {
        nmod_t mod;
        ulong n, a, r, x;
        int success, square;

        n = n_nth_prime(iter);
        nmod_init(&mod, n);

        for (a = 0; a < n; a++)
        {
            square = 0;
            for (x = 0; x < n && !square; x++)
                square = (nmod_mul(x, x, mod) == a);

            success = nmod_sqrt(&r, a, mod);

            if (success != square || (success && nmod_mul(r, r, mod) != a) || (success && r > n / 2))
                TEST_FUNCTION_FAIL("n = %wu, a = %wu, square = %d, success = %d, r = %wu\n",
                    n, a, square, success, r);
        }
    }

    /* composite moduli: any root returned must be genuine */
    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        nmod_t mod;
        ulong n, a, r;
        int success;

        n = n_randtest_not_zero(state);

        if (n_randint(state, 4) == 0)
        {
            /* odd squares, for which every Jacobi symbol is 0 or 1 */
            n = n_randtest(state) % (UWORD(1) << (FLINT_BITS / 2 - 1));
            n = n * n;
            n |= 1;
        }
        else if (n_randint(state, 4) == 0)
        {
            /* products of two primes = 1 mod 8 */
            ulong p, q;
            do { p = n_randprime(state, 2 + n_randint(state, FLINT_BITS / 2 - 3), 0); } while (p % 8 != 1);
            do { q = n_randprime(state, 2 + n_randint(state, FLINT_BITS / 2 - 3), 0); } while (q % 8 != 1);
            n = p * q;
        }

        nmod_init(&mod, n);

        a = n_randtest(state) % n;
        if (n_randint(state, 2))
            a = nmod_mul(a, a, mod);

        success = nmod_sqrt(&r, a, mod);

        if (success && nmod_mul(r, r, mod) != a)
            TEST_FUNCTION_FAIL("n = %wu, a = %wu, success = %d, r = %wu\n",
                n, a, success, r);
    }

    TEST_FUNCTION_END(state);
}
