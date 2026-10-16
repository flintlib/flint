/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(flint_mpn_divisible, state)
{
    slong iter;

    for (iter = 0; iter < 5000 * flint_test_multiplier(); iter++)
    {
        fmpz_t a, b, t;
        mpz_t ma, mb;
        int res, res2;
        slong abits, bbits;

        fmpz_init(a); fmpz_init(b); fmpz_init(t); mpz_init(ma); mpz_init(mb);

        switch (n_randint(state, 4))
        {
            case 0: bbits = 1 + n_randint(state, 130); abits = 1 + n_randint(state, 260); break;
            case 1: bbits = 1 + n_randint(state, 2000); abits = 1 + n_randint(state, 4000); break;
            case 2: bbits = 1 + n_randint(state, 20000); abits = 1 + n_randint(state, 60000); break;
            /* rarely, huge operands (the Newton-based code) */
            case 3:
                if (n_randint(state, 100) == 0)
                {
                    if (n_randint(state, 2))
                    {
                        /* long quotients */
                        bbits = FLINT_BITS * (FLINT_MPN_DIVEXACT_UNBALANCED_CUTOFF + n_randint(state, 100));
                        abits = bbits * (8 + n_randint(state, 3));
                    }
                    else
                    {
                        /* balanced, including quotients shorter than b */
                        bbits = FLINT_BITS * (FLINT_MPN_DIVEXACT_NEWTON_CUTOFF + n_randint(state, 100));
                        abits = FLINT_BITS * FLINT_MPN_DIVEXACT_NEWTON_CUTOFF
                            + n_randint(state, bbits - FLINT_BITS * FLINT_MPN_DIVEXACT_NEWTON_CUTOFF + FLINT_BITS);
                    }
                    break;
                }
                FLINT_FALLTHROUGH;
            default: bbits = 1 + n_randint(state, 300); abits = 1 + n_randint(state, 20000); break;
        }

        fmpz_randbits(b, state, bbits);
        fmpz_abs(b, b);
        if (fmpz_is_zero(b))
            fmpz_one(b);
        if (n_randint(state, 3) == 0)
            fmpz_mul_2exp(b, b, n_randint(state, 200));   /* even divisors, zero limbs */
        if (n_randint(state, 5) == 0)
            fmpz_mul_ui(b, b, 3 * 7 * 11);

        /* same number of limbs with the top limb of a below 2^v, so that
           the shifted a' is one limb shorter than b' */
        if (n_randint(state, 8) == 0)
        {
            slong k = 1 + n_randint(state, 40);
            unsigned int v = 1 + n_randint(state, FLINT_BITS - 1);

            fmpz_randbits(b, state, FLINT_BITS * k);
            fmpz_abs(b, b);
            fmpz_setbit(b, FLINT_BITS * k - 1);
            fmpz_tdiv_q_2exp(b, b, v);
            fmpz_mul_2exp(b, b, v);
            fmpz_setbit(b, v);   /* exactly v trailing zero bits */

            fmpz_randbits(a, state, FLINT_BITS * (k - 1) + v);
            fmpz_abs(a, a);
            fmpz_setbit(a, FLINT_BITS * (k - 1) + v - 1);   /* k limbs */
            fmpz_tdiv_q_2exp(a, a, v);
            fmpz_mul_2exp(a, a, v);
            if (fmpz_is_zero(a))
                fmpz_mul_2exp(a, b, 1);

            goto compare;
        }

        fmpz_randbits(a, state, abits);
        fmpz_abs(a, a);
        if (n_randint(state, 2))
        {
            /* divisible or nearly */
            fmpz_mul(a, a, b);
            if (n_randint(state, 3) == 0)
            {
                fmpz_randbits(t, state, n_randint(state, bbits + 1));
                fmpz_add(a, a, t);
            }
            else if (n_randint(state, 3) == 0)
            {
                /* one bit off */
                fmpz_combit(a, n_randint(state, fmpz_bits(a) + 1));
            }
        }

compare:
        fmpz_get_mpz(ma, a);
        fmpz_get_mpz(mb, b);

        res2 = mpz_divisible_p(ma, mb);
        res = flint_mpn_divisible(ma->_mp_d, ma->_mp_size, mb->_mp_d, mb->_mp_size);

        if (res != res2)
            TEST_FUNCTION_FAIL("res = %d, res2 = %d\na = %{fmpz}\nb = %{fmpz}\n", res, res2, a, b);

        /* a with zero top limbs (and a = 0) */
        if (n_randint(state, 4) == 0)
        {
            mp_size_t an = ma->_mp_size, pad = n_randint(state, 3);
            mp_ptr ap = flint_calloc(an + pad + 1, sizeof(mp_limb_t));
            if (an > 0)
                flint_mpn_copyi(ap, ma->_mp_d, an);
            res = _flint_mpn_divisible(ap, an + pad, mb->_mp_d, mb->_mp_size);
            flint_free(ap);
            if (res != res2)
                TEST_FUNCTION_FAIL("zero-padded: res = %d, res2 = %d\na = %{fmpz}\nb = %{fmpz}\npad = %wd\n",
                    res, res2, a, b, pad);
        }

        /* fmpz level */
        if (fmpz_divisible(a, b) != res2)
            TEST_FUNCTION_FAIL("fmpz_divisible: a = %{fmpz}\nb = %{fmpz}\n", a, b);

        fmpz_clear(a); fmpz_clear(b); fmpz_clear(t); mpz_clear(ma); mpz_clear(mb);
    }

    /* the exported kernels with a two-limb divisor whose odd part has one
       limb (flint_mpn_divisible handles that shape itself) */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        mp_limb_t b[2], c;
        mp_ptr a, x;
        mp_size_t an, xn;
        unsigned int v;
        mpz_t ma, mb;
        int res, res2;

        /* b = c 2^v with c odd of one full limb, so that b[1] < 2^v */
        v = 1 + n_randint(state, FLINT_BITS - 1);
        c = n_randtest(state) | 1 | (UWORD(1) << (FLINT_BITS - 1));
        b[0] = c << v;
        b[1] = c >> (FLINT_BITS - v);

        xn = 1 + n_randint(state, 12);
        a = flint_malloc((xn + 2) * sizeof(mp_limb_t));
        x = flint_malloc(xn * sizeof(mp_limb_t));
        flint_mpn_rrandom(x, state, xn);
        if (x[xn - 1] == 0)
            x[xn - 1] = 1;

        if (n_randint(state, 2))
        {
            /* a multiple of b, possibly off by a multiple of 2^v */
            if (xn >= 2)
                flint_mpn_mul(a, x, xn, b, 2);
            else
                flint_mpn_mul(a, b, 2, x, xn);
            an = xn + 2;
            if (n_randint(state, 3) == 0)
                mpn_add_1(a, a, an, UWORD(1) << v);
        }
        else
        {
            flint_mpn_copyi(a, x, xn);
            a[xn] = n_randtest(state);
            an = xn + 1;
            a[0] &= ~((UWORD(1) << v) - 1);   /* 2^v divides a */
        }
        while (an > 2 && a[an - 1] == 0)
            an--;
        if (a[an - 1] == 0)
            a[an - 1] = 1;

        mpz_init(ma); mpz_init(mb);
        mpz_import(ma, an, -1, sizeof(mp_limb_t), 0, 0, a);
        mpz_import(mb, 2, -1, sizeof(mp_limb_t), 0, 0, b);
        res2 = mpz_divisible_p(ma, mb);

        res = _flint_mpn_divisible_small(a, an, b, 2, v);
        if (res != res2)
            TEST_FUNCTION_FAIL("_flint_mpn_divisible_small: res = %d, res2 = %d\na = %{ulong*}\nb = %{ulong*}\n",
                res, res2, a, an, b, 2);

        res = _flint_mpn_divisible_bdiv(a, an, b, 2, v);
        if (res != res2)
            TEST_FUNCTION_FAIL("_flint_mpn_divisible_bdiv: res = %d, res2 = %d\na = %{ulong*}\nb = %{ulong*}\n",
                res, res2, a, an, b, 2);

        mpz_clear(ma); mpz_clear(mb);
        flint_free(a);
        flint_free(x);
    }

    TEST_FUNCTION_END(state);
}
