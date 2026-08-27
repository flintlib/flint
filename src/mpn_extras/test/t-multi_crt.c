/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "mpn_extras.h"

/* tests both flint_mpn_multi_mod and flint_mpn_multi_crt */
TEST_FUNCTION_START(flint_mpn_multi_crt, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        slong n, i, bits, xn;
        ulong * primes, * res, * x, * out;
        flint_mpn_crt_t C;
        fmpz_t X, Y, P, T;
        int sign, neg, flags;

        flint_set_num_threads(1 + n_randint(state, 3));

        n = 1 + n_randint(state, 100);
        if (n_randint(state, 4) == 0)
            n = 1 + n_randint(state, 20);
        if (n_randint(state, 100) == 0)
            n = 2000 + n_randint(state, 6000);   /* exercise the threaded paths */
        bits = 2 + n_randint(state, FLINT_BITS - 1);
        if (n_randint(state, 5) == 0)
            bits = FLINT_BITS;

        primes = flint_malloc(n * sizeof(ulong));

        if (n_randint(state, 4) != 0)
        {
            /* random pairwise coprime moduli of varying sizes */
            n = FLINT_MIN(n, 30);
            for (i = 0; i < n; i++)
            {
                slong k, tries = 0;
                do {
                    slong b = 1 + n_randint(state, bits);
                    primes[i] = n_randtest_bits(state, b);
                    if (primes[i] == 0)
                        primes[i] = 1;
                    for (k = 0; k < i; k++)
                        if (n_gcd(primes[i], primes[k]) != 1)
                            break;
                    tries++;
                } while (k < i && tries < 100);
                if (k < i)
                {
                    n = i;
                    break;
                }
            }
        }
        else
        {
            ulong p = n_nextprime(UWORD(1) << (bits - 1), 1);
            if (bits == FLINT_BITS && n_randint(state, 2))
                p = n_nextprime(UWORD(1) << (FLINT_BITS - 2), 1) ;
            for (i = 0; i < n; i++)
            {
                primes[i] = p;
                p = n_nextprime(p, 1);
                if (p < primes[i] || (bits < FLINT_BITS && p >= (UWORD(1) << bits)))
                {
                    n = i + 1;
                    break;
                }
            }
        }

        flags = 1 + n_randint(state, 3);

        if (n_randint(state, 2))
            flint_mpn_crt_init2(C, primes, n, flags);
        else
            flint_mpn_crt_init_tuned(C, primes, n, flags,
                64 * (1 + n_randint(state, 40)),
                64 * (1 + n_randint(state, 80)),
                1 + n_randint(state, 30));

        fmpz_init(X);
        fmpz_init(Y);
        fmpz_init(P);
        fmpz_init(T);
        fmpz_set_ui_array(P, C->prod, C->prod_len);

        fmpz_one(T);
        for (i = 0; i < n; i++)
            fmpz_mul_ui(T, T, primes[i]);
        if (!fmpz_equal(P, T))
            TEST_FUNCTION_FAIL("product\n");

        res = flint_malloc(n * sizeof(ulong));
        out = flint_malloc(C->prod_len * sizeof(ulong));

        /* multi_mod: random x, possibly larger than P */
        if (flags & FLINT_MPN_CRT_MOD)
        {
            slong xbits = n_randint(state, 2 * fmpz_bits(P) + 200);
            fmpz_randtest_unsigned(X, state, xbits);
            if (n_randint(state, 4) == 0)
                fmpz_randm(X, state, P);
            xn = fmpz_size(X);
            x = flint_malloc((xn + 1) * sizeof(ulong));
            fmpz_get_ui_array(x, xn + 1, X);

            if (n_randint(state, 4) == 0)
                flint_mpn_multi_mod_once(res, x, xn, primes, n);
            else
                flint_mpn_multi_mod(res, x, xn, C, NULL);
            for (i = 0; i < n; i++)
            {
                ulong r = fmpz_fdiv_ui(X, primes[i]);
                if (r != res[i])
                    TEST_FUNCTION_FAIL("mod: n = %wd, bits = %wd, i = %wd, got %wu, want %wu\n", n, bits, i, res[i], r);
            }
            flint_free(x);
        }

        /* multi_crt */
        if (flags & FLINT_MPN_CRT_CRT)
        {
            sign = n_randint(state, 2);
            fmpz_randm(X, state, P);
            /* small values (exercise the shortcut), including negative ones */
            if (n_randint(state, 3) == 0)
            {
                fmpz_randtest(X, state, n_randint(state, 600));
                fmpz_mod(X, X, P);
            }
            for (i = 0; i < n; i++)
                res[i] = fmpz_fdiv_ui(X, primes[i]);
            if (n_randint(state, 4) == 0)
            {
                nn_ptr out2 = flint_malloc(2 * n * sizeof(ulong));
                slong outn;
                neg = flint_mpn_multi_crt_once(out2, &outn, out2 + n, res, primes, n, sign);
                if (outn != C->prod_len || mpn_cmp(out2 + n, C->prod, C->prod_len) != 0)
                    TEST_FUNCTION_FAIL("crt_once: product\n");
                flint_mpn_copyi(out, out2, C->prod_len);
                flint_free(out2);
            }
            else
                neg = flint_mpn_multi_crt(out, res, C, sign, NULL);
            fmpz_set_ui_array(Y, out, C->prod_len);
            if (neg)
                fmpz_neg(Y, Y);
            if (sign)
                _fmpz_smod(T, X, P, 1, Y);   /* Y clobbered as temp, recompute */
            else
                fmpz_set(T, X);
            fmpz_set_ui_array(Y, out, C->prod_len);
            if (neg)
                fmpz_neg(Y, Y);
            if (!fmpz_equal(T, Y))
                TEST_FUNCTION_FAIL("crt: n = %wd, bits = %wd, sign = %d\nX = %{fmpz}\nY = %{fmpz}\n", n, bits, sign, X, Y);
        }

        /* vector versions against the scalar ones */
        if ((flags & FLINT_MPN_CRT_MOD) && (flags & FLINT_MPN_CRT_CRT) && n_randint(state, 3) == 0)
        {
            slong len = 1 + n_randint(state, 10), vxn = n_randint(state, C->prod_len + 2), k;
            slong ostride = len + n_randint(state, 3), cstride = C->prod_len + n_randint(state, 3);
            nn_ptr xv = flint_malloc(len * FLINT_MAX(vxn, 1) * sizeof(ulong));
            nn_ptr rv = flint_malloc(ostride * n * sizeof(ulong));
            nn_ptr rv2 = flint_malloc(n * sizeof(ulong));
            nn_ptr ov = flint_malloc(len * cstride * sizeof(ulong));
            nn_ptr ov2 = flint_malloc(C->prod_len * sizeof(ulong));
            int * negv = flint_malloc(len * sizeof(int));

            for (k = 0; k < len * vxn; k++)
                xv[k] = n_randtest(state);
            /* make some entries small */
            for (k = 0; k < len; k++)
            {
                slong keep = FLINT_MIN(vxn, n_randint(state, 4));
                if (n_randint(state, 2))
                    flint_mpn_zero(xv + k * vxn + keep, vxn - keep);
            }
            flint_mpn_multi_mod_vec(rv, ostride, xv, vxn, len, C, NULL);
            for (k = 0; k < len; k++)
            {
                flint_mpn_multi_mod(rv2, xv + k * vxn, vxn, C, NULL);
                for (i = 0; i < n; i++)
                    if (rv2[i] != rv[i * ostride + k])
                        TEST_FUNCTION_FAIL("mod_vec: n = %wd, vxn = %wd, len = %wd\n", n, vxn, len);
            }

            sign = n_randint(state, 2);
            flint_mpn_multi_crt_vec(ov, cstride, negv, rv, ostride, len, C, sign, NULL);
            for (k = 0; k < len; k++)
            {
                for (i = 0; i < n; i++)
                    rv2[i] = rv[i * ostride + k];
                neg = flint_mpn_multi_crt(ov2, rv2, C, sign, NULL);
                if (neg != negv[k] || mpn_cmp(ov2, ov + k * cstride, C->prod_len) != 0)
                    TEST_FUNCTION_FAIL("crt_vec: n = %wd, len = %wd\n", n, len);
            }

            flint_free(xv);
            flint_free(rv);
            flint_free(rv2);
            flint_free(ov);
            flint_free(ov2);
            flint_free(negv);
        }

        fmpz_clear(X);
        fmpz_clear(Y);
        fmpz_clear(P);
        fmpz_clear(T);
        flint_free(res);
        flint_free(out);
        flint_free(primes);
        flint_mpn_crt_clear(C);
    }

    TEST_FUNCTION_END(state);
}
