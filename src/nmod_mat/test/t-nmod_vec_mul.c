/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "nmod_mat/impl.h"

/* moduli at the limits of the vectorized paths of the products with few
   rows (64-bit: u52 and fp50 up to 2^52 and 2^50, u64 beyond, in its
   split mode up to 2^61) */
static ulong
_vec_mul_modulus(flint_rand_t state)
{
    ulong n;

    switch (n_randint(state, 20))
    {
        case 0: n = UWORD(1) << 26; break;
        case 1: n = (UWORD(1) << 26) + 1; break;
#if FLINT_BITS == 64
        case 2: n = UWORD(1) << 32; break;
        case 3: n = (UWORD(1) << 50) - 1; break;
        case 4: n = UWORD(1) << 50; break;
        case 5: n = (UWORD(1) << 50) + 1; break;
        case 6: n = UWORD(1) << 52; break;
        case 7: n = (UWORD(1) << 52) + 1; break;
        case 8: n = UWORD(1) << 61; break;
        case 9: n = (UWORD(1) << 61) + 1; break;
        case 10: n = (UWORD(1) << 62) - 1; break;
        case 11: n = (UWORD(1) << 62) + 1; break;
        case 12: n = (UWORD(1) << 63) - 1; break;
        case 13: n = (UWORD(1) << 63) + 1; break;
        case 14: n = UWORD_MAX - 58; break;  /* 2^64 - 59 */
        case 15: n = UWORD_MAX; break;
        case 16: n = n_randbits(state, 50 + n_randint(state, 15)); break;
#endif
        default:
            n = n_randbits(state, 2 + n_randint(state, FLINT_BITS - 1));
            break;
    }

    return FLINT_MAX(n, UWORD(2));
}

TEST_FUNCTION_START(nmod_mat_nmod_vec_mul, state)
{
    slong i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        ulong p;
        nmod_mat_t A, B, C;
        ulong * a, * c;
        ulong ** aa, ** cc;
        slong j, m, n, alen;

        p = n_randtest_not_zero(state);
        m = n_randint(state, 50);
        n = n_randint(state, 50);
        alen = n_randint(state, 50);

        nmod_mat_init(C, 1, n, p);
        nmod_mat_init(A, 1, m, p);
        nmod_mat_init(B, m, n, p);
        c = _nmod_vec_init(n);
        a = _nmod_vec_init(alen);

        nmod_mat_randtest(B, state);
        _nmod_vec_randtest(c, state, n, B->mod);
        _nmod_vec_randtest(a, state, alen, B->mod);

        cc = FLINT_ARRAY_ALLOC(n, ulong*);
        for (j = 0; j < n; j++)
        {
            cc[j] = FLINT_ARRAY_ALLOC(1, ulong);
            cc[j][0] = c[j];
        }

        aa = FLINT_ARRAY_ALLOC(alen, ulong*);
        for (j = 0; j < alen; j++)
        {
            aa[j] = FLINT_ARRAY_ALLOC(1, ulong);
            aa[j][0] = a[j];
        }

        nmod_mat_nmod_vec_mul(c, a, alen, B);
        nmod_mat_nmod_vec_mul_ptr(cc, (const ulong * const *)aa, alen, B);

        /* supposed to match mul of the chopped or zero-extended a */
        for (j = 0; j < m && j < alen; j++)
            nmod_mat_entry(A, 0, j) = a[j];

        nmod_mat_mul(C, A, B);

        for (j = 0; j < n; j++)
            if (nmod_mat_entry(C, 0, j) != c[j] || nmod_mat_entry(C, 0, j) != cc[j][0])
                TEST_FUNCTION_FAIL("");

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
        _nmod_vec_clear(c);
        _nmod_vec_clear(a);

        for (j = 0; j < n; j++)
        {
            flint_free(cc[j]);
        }
        flint_free(cc);

        for (j = 0; j < alen; j++)
        {
            flint_free(aa[j]);
        }
        flint_free(aa);
    }

    /* against one modular multiplication per entry, on the dimensions
       and moduli of the vectorized paths: rows beyond the chunks of the
       accumulators, columns beyond the blocks, moduli at the limits of
       each method, entries n - 1 (largest products) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_mat_t B;
        ulong * a, * c, * d;
        ulong n;
        slong j, k, m;
        int mode;

        n = _vec_mul_modulus(state);

        k = 1 + n_randint(state, n_randint(state, 4) ? 80 : 1200);
        m = 1 + n_randint(state, n_randint(state, 4) ? 80 : 4500);
        mode = n_randint(state, 3);

        nmod_mat_init(B, k, m, n);
        a = _nmod_vec_init(k);
        c = _nmod_vec_init(m + 1);
        d = _nmod_vec_init(m);

        for (j = 0; j < k * m; j++)
            B->entries[j] = mode ? n - 1 - (mode == 2 ? n_randint(state, 2) : 0)
                                 : n_randint(state, n);
        for (j = 0; j < k; j++)
            a[j] = mode ? n - 1 - (mode == 2 ? n_randint(state, 2) : 0)
                        : n_randint(state, n);
        c[m] = 12345;

        nmod_mat_nmod_vec_mul(c, a, k, B);

        _nmod_vec_scalar_mul_nmod(d, nmod_mat_entry_ptr(B, 0, 0), m, a[0], B->mod);
        for (j = 1; j < k; j++)
            _nmod_vec_scalar_addmul_nmod(d, nmod_mat_entry_ptr(B, j, 0), m, a[j], B->mod);

        if (!_nmod_vec_equal(c, d, m) || c[m] != 12345)
            TEST_FUNCTION_FAIL("k: %wd, m: %wd, n: %wu, mode: %d\n", k, m, n, mode);

        nmod_mat_clear(B);
        _nmod_vec_clear(a);
        _nmod_vec_clear(c);
        _nmod_vec_clear(d);
    }

    /* _nmod_mat_mul_rows_simd, for 1 to NMOD_MAT_MUL_ROWS_MAX rows, on a
       window of a wider matrix, against one modular multiplication per
       entry: lengths around multiples of the folds of the u64 path (every
       2^(63-b) rows for b-bit moduli up to 61 bits) and of its chunks of
       16 rows (B larger than L2), entries n - 1 (largest products) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_mat_t PB, B;
        ulong * a[NMOD_MAT_MUL_ROWS_MAX], * c[NMOD_MAT_MUL_ROWS_MAX];
        ulong * d;
        ulong n;
        slong j, t, k, m, r, extra, off;
        int mode;

        n = _vec_mul_modulus(state);
        r = 1 + n_randint(state, NMOD_MAT_MUL_ROWS_MAX);

        switch (n_randint(state, 4))
        {
            case 0:
                k = (WORD(1) << n_randint(state, 7)) * (1 + n_randint(state, 4))
                    + n_randint(state, 3) - 1;
                m = 1 + n_randint(state, 40);
                break;
            case 1:
#if FLINT_BITS == 64
                /* folds every 2048 rows for 52 bits */
                n = (UWORD(1) << 51) + n_randbits(state, 51);
#endif
                k = 2047 + n_randint(state, 3);
                m = 1 + n_randint(state, 20);
                break;
            case 2:
                k = 250 + n_randint(state, 150);
                m = 600 + n_randint(state, 1500);
                break;
            default:
                k = 1 + n_randint(state, 100);
                m = 1 + n_randint(state, 200);
                break;
        }
        k = FLINT_MAX(k, 1);
        extra = n_randint(state, 2) ? 0 : 1 + n_randint(state, 9);
        off = n_randint(state, extra + 1);
        mode = n_randint(state, 3);

        nmod_mat_init(PB, k, m + extra, n);
        for (j = 0; j < k * (m + extra); j++)
            PB->entries[j] = mode ? n - 1 - (mode == 2 ? n_randint(state, 2) : 0)
                                  : n_randint(state, n);
        nmod_mat_window_init(B, PB, 0, off, k, off + m);

        d = _nmod_vec_init(m);
        for (t = 0; t < r; t++)
        {
            a[t] = _nmod_vec_init(k);
            c[t] = _nmod_vec_init(m + 1);
            for (j = 0; j < k; j++)
                a[t][j] = mode ? n - 1 - (mode == 2 ? n_randint(state, 2) : 0)
                               : n_randint(state, n);
            c[t][m] = 12345;
        }

        if (_nmod_mat_mul_rows_simd(c, (const ulong * const *) a, r, k, B))
        {
            for (t = 0; t < r; t++)
            {
                _nmod_vec_scalar_mul_nmod(d, nmod_mat_entry_ptr(B, 0, 0), m,
                                          a[t][0], B->mod);
                for (j = 1; j < k; j++)
                    _nmod_vec_scalar_addmul_nmod(d, nmod_mat_entry_ptr(B, j, 0), m,
                                                 a[t][j], B->mod);

                if (!_nmod_vec_equal(c[t], d, m) || c[t][m] != 12345)
                    TEST_FUNCTION_FAIL("rows: r: %wd, row: %wd, k: %wd, m: %wd, "
                                       "extra: %wd, n: %wu, mode: %d\n",
                                       r, t, k, m, extra, n, mode);
            }
        }
        else if (NMOD_MAT_NMOD_VEC_MUL_IS_SIMD(n))
            TEST_FUNCTION_FAIL("rows: declined, n: %wu\n", n);

        /* nmod_mat_nmod_vec_mul on the same window (its other paths) */
        nmod_mat_nmod_vec_mul(c[0], a[0], k, B);
        _nmod_vec_scalar_mul_nmod(d, nmod_mat_entry_ptr(B, 0, 0), m,
                                  a[0][0], B->mod);
        for (j = 1; j < k; j++)
            _nmod_vec_scalar_addmul_nmod(d, nmod_mat_entry_ptr(B, j, 0), m,
                                         a[0][j], B->mod);
        if (!_nmod_vec_equal(c[0], d, m) || c[0][m] != 12345)
            TEST_FUNCTION_FAIL("window vec_mul: k: %wd, m: %wd, extra: %wd, "
                               "n: %wu, mode: %d\n", k, m, extra, n, mode);

        for (t = 0; t < r; t++)
        {
            _nmod_vec_clear(a[t]);
            _nmod_vec_clear(c[t]);
        }
        _nmod_vec_clear(d);
        nmod_mat_window_clear(B);
        nmod_mat_clear(PB);
    }

#if FLINT_BITS == 64
    /* the folds of the u64 path every 2^(63-b) rows for b-bit moduli are
       tight: moduli 2^b - 1, whose entries n - 1 have nearly full limbs,
       with all entries n - 1 and one row more than the interval; since
       (n - 1)^2 = 1 mod n, every entry of the product is k */
    {
        slong b;

        for (b = 51; b <= 61; b++)
        {
            ulong n = (UWORD(1) << b) - 1;
            slong k = (WORD(1) << (63 - b)) + 1;
            slong m = 8 + n_randint(state, 40), j, t, r;
            ulong * a[NMOD_MAT_MUL_ROWS_MAX], * c[NMOD_MAT_MUL_ROWS_MAX];
            nmod_mat_t B;

            nmod_mat_init(B, k, m, n);
            for (j = 0; j < k * m; j++)
                B->entries[j] = n - 1;

            for (t = 0; t < NMOD_MAT_MUL_ROWS_MAX; t++)
            {
                a[t] = _nmod_vec_init(k);
                c[t] = _nmod_vec_init(m);
                for (j = 0; j < k; j++)
                    a[t][j] = n - 1;
            }

            nmod_mat_nmod_vec_mul(c[0], a[0], k, B);
            for (j = 0; j < m; j++)
                if (c[0][j] != (ulong) k)
                    TEST_FUNCTION_FAIL("2^b - 1: vec_mul, b: %wd, m: %wd\n", b, m);

            for (r = 1; r <= NMOD_MAT_MUL_ROWS_MAX; r++)
            {
                if (!_nmod_mat_mul_rows_simd(c, (const ulong * const *) a, r, k, B))
                    continue;

                for (t = 0; t < r; t++)
                    for (j = 0; j < m; j++)
                        if (c[t][j] != (ulong) k)
                            TEST_FUNCTION_FAIL("2^b - 1: rows, b: %wd, r: %wd, "
                                               "m: %wd\n", b, r, m);
            }

            for (t = 0; t < NMOD_MAT_MUL_ROWS_MAX; t++)
            {
                _nmod_vec_clear(a[t]);
                _nmod_vec_clear(c[t]);
            }
            nmod_mat_clear(B);
        }
    }
#endif

    TEST_FUNCTION_END(state);
}
