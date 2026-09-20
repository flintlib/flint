/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"
#include "fmpz.h"
#include "ulong_extras.h"

/* A nonnegative fmpz written to exactly n limbs. */
static void
_set_fmpz(nn_ptr r, mp_size_t n, const fmpz_t x)
{
    mpz_t t;

    mpz_init(t);
    fmpz_get_mpz(t, x);
    flint_mpn_copyi(r, t->_mp_d, t->_mp_size);
    flint_mpn_zero(r + t->_mp_size, n - t->_mp_size);
    mpz_clear(t);
}

/*
    A prime of the given number of limbs. With val2 > 0 the prime is
    congruent to 1 modulo 2^val2 and not modulo 2^(val2+1), which is what
    selects the branch taken: val2 = 1 is the p = 3 mod 4 shortcut, val2 = 2
    is Atkin's formula, and larger values are Tonelli and Shanks.
*/
static void
_random_prime(fmpz_t p, flint_rand_t state, mp_size_t n, slong val2)
{
    flint_bitcnt_t bits = n * FLINT_BITS;
    fmpz_t q;

    fmpz_init(q);

    while (1)
    {
        if (val2 <= 0)
        {
            fmpz_randprime(p, state, bits, 0);

            if (fmpz_size(p) == n && fmpz_is_odd(p))
                break;

            continue;
        }

        fmpz_randbits(q, state, bits - val2);
        fmpz_abs(q, q);
        fmpz_setbit(q, 0);
        fmpz_setbit(q, bits - val2 - 1);
        fmpz_mul_2exp(p, q, val2);
        fmpz_add_ui(p, p, 1);

        if (fmpz_size(p) == n && fmpz_is_probabprime(p))
            break;
    }

    fmpz_clear(q);
}

TEST_FUNCTION_START(flint_mpn_sqrtmod, state)
{
    slong iter;

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        fmpz_t p, x, y;
        nn_ptr d, dnormed, dinv, a, r;
        mp_size_t n;
        slong val2, j;
        flint_bitcnt_t norm;

        n = 1 + n_randint(state, 6);
        val2 = (slong) n_randint(state, FLINT_MIN(12, n * FLINT_BITS - 8));

        fmpz_init(p);
        fmpz_init(x);
        fmpz_init(y);

        _random_prime(p, state, n, val2);

        d = flint_malloc(5 * n * sizeof(ulong));
        dnormed = d + n;
        dinv = dnormed + n;
        a = dinv + n;
        r = a + n;

        _set_fmpz(d, n, p);

        norm = flint_clz(d[n - 1]);

        if (norm)
            mpn_lshift(dnormed, d, n, norm);
        else
            flint_mpn_copyi(dnormed, d, n);

        flint_mpn_preinvn(dinv, dnormed, n);

        for (j = 0; j < 8; j++)
        {
            int square, got, got_preinv;

            if (j == 0)
                fmpz_zero(x);
            else if (j == 1)
                fmpz_one(x);
            else
            {
                fmpz_randm(x, state, p);

                if (j % 2 == 0)     /* force a square */
                {
                    fmpz_mul(x, x, x);
                    fmpz_mod(x, x, p);
                }
            }

            _set_fmpz(a, n, x);

            /* the Jacobi symbol decides, the modulus being prime */
            square = (fmpz_jacobi(x, p) != -1);

            if (flint_mpn_is_square_mod(a, d, n) != square)
            {
                flint_printf("FAIL: is_square_mod\n");
                flint_printf("p = "); fmpz_print(p); flint_printf("\n");
                flint_printf("x = "); fmpz_print(x); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            /* the modulus is prime, so the algorithm cannot fail */
            got = flint_mpn_sqrtmod(r, a, d, n);

            if (got != square)
            {
                flint_printf("FAIL: sqrtmod reports %d, expected %d\n", got, square);
                flint_printf("p = "); fmpz_print(p); flint_printf("\n");
                flint_printf("x = "); fmpz_print(x); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            if (got)
            {
                fmpz_set_ui_array(y, r, n);
                fmpz_mul(y, y, y);
                fmpz_mod(y, y, p);

                if (!fmpz_equal(y, x))
                {
                    flint_printf("FAIL: sqrtmod(x)^2 != x\n");
                    flint_printf("p = "); fmpz_print(p); flint_printf("\n");
                    flint_printf("x = "); fmpz_print(x); flint_printf("\n");
                    fflush(stdout);
                    flint_abort();
                }
            }
            else
                FLINT_TEST(flint_mpn_zero_p(r, n));

            /* the preinverse must change nothing */
            flint_mpn_zero(r, n);
            got_preinv = flint_mpn_sqrtmod_preinv(r, a, d, n, dinv, norm);

            if (got_preinv != got)
            {
                flint_printf("FAIL: sqrtmod_preinv disagrees with sqrtmod\n");
                flint_printf("p = "); fmpz_print(p); flint_printf("\n");
                flint_printf("x = "); fmpz_print(x); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            if (got_preinv)
            {
                fmpz_set_ui_array(y, r, n);
                fmpz_mul(y, y, y);
                fmpz_mod(y, y, p);
                FLINT_TEST(fmpz_equal(y, x));
            }
        }

        flint_free(d);
        fmpz_clear(p);
        fmpz_clear(x);
        fmpz_clear(y);
    }

    /*
        An odd modulus that is not prime. Nothing is promised except that a
        Jacobi symbol of -1, which is what a zero return reports, really
        does rule out a square root.
    */
    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        fmpz_t nn, p1, p2, x;
        nn_ptr d, a, r;
        mp_size_t n;

        fmpz_init(nn);
        fmpz_init(p1);
        fmpz_init(p2);
        fmpz_init(x);

        fmpz_randprime(p1, state, 2 + n_randint(state, 78), 0);
        fmpz_randprime(p2, state, 2 + n_randint(state, 78), 0);
        fmpz_mul(nn, p1, p2);

        if (fmpz_is_odd(nn) && fmpz_cmp_ui(nn, 8) > 0)
        {
            n = fmpz_size(nn);

            d = flint_malloc(3 * n * sizeof(ulong));
            a = d + n;
            r = a + n;

            _set_fmpz(d, n, nn);
            fmpz_randm(x, state, nn);
            _set_fmpz(a, n, x);

            if (flint_mpn_is_square_mod(a, d, n) == 0)
                FLINT_TEST(fmpz_jacobi(x, nn) == -1);

            /* a modulus that is not prime may only fail, never lie */
            if (n > 1 && flint_mpn_sqrtmod(r, a, d, n) == 1)
            {
                fmpz_t y;
                fmpz_init(y);
                fmpz_set_ui_array(y, r, n);
                fmpz_mul(y, y, y);
                fmpz_mod(y, y, nn);
                FLINT_TEST(fmpz_equal(y, x) || fmpz_jacobi(x, nn) != -1);
                fmpz_clear(y);
            }

            flint_free(d);
        }

        fmpz_clear(nn);
        fmpz_clear(p1);
        fmpz_clear(p2);
        fmpz_clear(x);
    }

    TEST_FUNCTION_END(state);
}
