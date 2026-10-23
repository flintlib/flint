/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpz_sqrtmod, state)
{
    int i, result;

    for (i = 0; i < 100 * flint_test_multiplier(); i++) /* Test random integers */
    {
        int ans;
        fmpz_t a, b, c, p;
        ulong prime;

        prime = n_randint(state, UWORD(1) << (FLINT_BITS - 1));
        prime = n_nextprime(prime, 1);

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(p);

        fmpz_set_ui(p, prime);
        fmpz_randm(a, state, p);

        ans = fmpz_sqrtmod(b, a, p);

        fmpz_mul(c, b, b);
        fmpz_mod(c, c, p);

        result = (ans == 0 || fmpz_equal(a, c)) && _fmpz_is_canonical(b);
        if (!result)
        {
            flint_printf("FAIL (random):\n");
            flint_printf("p = "), fmpz_print(p), flint_printf("\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("c = "), fmpz_print(c), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(p);
    }

    for (i = 0; i < 100 * flint_test_multiplier(); i++) /* Test random squares */
    {
        int ans;
        fmpz_t a, b, c, d, p;
        ulong prime;

        prime = n_randint(state, UWORD(1) << (FLINT_BITS - 1));
        prime = n_nextprime(prime, 1);

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(d);
        fmpz_init(p);

        fmpz_set_ui(p, prime);
        do
            fmpz_randm(b, state, p);
        while (fmpz_is_zero(b));

        fmpz_mul(a, b, b);
        fmpz_mod(a, a, p);

        /* check a special case */
        if (i == 0)
        {
            fmpz_set_str(p, "15951355998396157", 10);
            fmpz_set_str(a, "7009303413761286", 10);
        }

        ans = fmpz_sqrtmod(c, a, p);

        fmpz_mul(d, c, c);
        fmpz_mod(d, d, p);

        result = (ans && fmpz_equal(a, d)) && _fmpz_is_canonical(c);
        if (!result)
        {
            flint_printf("FAIL (squares):\n");
            flint_printf("p            = "), fmpz_print(p), flint_printf("\n");
            flint_printf("a (= b^2)    = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b            = "), fmpz_print(b), flint_printf("\n");
            flint_printf("c (= sqrt(a) = "), fmpz_print(c), flint_printf("\n");
            flint_printf("d (= c^2)    = "), fmpz_print(d), flint_printf("\n");
            flint_printf("ans          = %d\n", ans);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(d);
        fmpz_clear(p);
    }

    /*
        Large primes, with the 2-adic valuation of p - 1 controlled so that
        every branch of flint_mpn_sqrtmod is reached: val2 = 1 is the
        p = 3 mod 4 shortcut, val2 = 2 is Atkin's formula, and larger values
        are Tonelli and Shanks.
    */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        int ans;
        fmpz_t a, b, c, p, q;
        slong val2;
        flint_bitcnt_t bits;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(p);
        fmpz_init(q);

        bits = 64 + n_randint(state, 320);
        val2 = 1 + (slong) n_randint(state, 10);

        while (1)
        {
            fmpz_randbits(q, state, bits - val2);
            fmpz_abs(q, q);
            fmpz_setbit(q, 0);
            fmpz_setbit(q, bits - val2 - 1);
            fmpz_mul_2exp(p, q, val2);
            fmpz_add_ui(p, p, 1);

            if (fmpz_is_probabprime(p))
                break;
        }

        fmpz_randm(a, state, p);

        if (i % 2 == 0)     /* force a square */
        {
            fmpz_mul(a, a, a);
            fmpz_mod(a, a, p);
        }

        ans = fmpz_sqrtmod(b, a, p);

        /* the Jacobi symbol decides, the modulus being prime */
        result = (ans == (fmpz_jacobi(a, p) != -1)) && _fmpz_is_canonical(b);

        if (result && ans)
        {
            fmpz_mul(c, b, b);
            fmpz_mod(c, c, p);
            result = fmpz_equal(c, a);
        }

        if (result && !ans)
            result = fmpz_is_zero(b);

        if (!result)
        {
            flint_printf("FAIL (large primes):\n");
            flint_printf("val2 = %wd\n", val2);
            flint_printf("p = "), fmpz_print(p), flint_printf("\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("ans = %d\n", ans);
            fflush(stdout);
            flint_abort();
        }

        /* aliasing must not change the answer */
        fmpz_set(c, a);
        FLINT_TEST(fmpz_sqrtmod(c, c, p) == ans);
        FLINT_TEST(fmpz_equal(c, b));

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(p);
        fmpz_clear(q);
    }

    TEST_FUNCTION_END(state);
}
