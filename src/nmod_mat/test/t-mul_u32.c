/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_mat.h"

#if FLINT_BITS == 64
# define TWO32 (UWORD(1) << 32)
#else
/* unused in effect: the 32-bit build declines every multiplication */
# define TWO32 UWORD(0)
#endif

/*
    Entries of maximal absolute value in the symmetric lift: n/2 lifts to
    +floor(n/2) and n - floor(n/2) lifts to -floor(n/2) (for n odd) or
    +n/2 (for n even). With all entries like this a dot product grows by
    floor(n/2)^2 at every step, which is exactly the bound the delayed
    reduction is designed against.
*/
static void
nmod_mat_randfull_extreme(nmod_mat_t mat, flint_rand_t state, int mode)
{
    slong i, j;
    ulong n = mat->mod.n;
    ulong pos = n / 2;
    ulong neg = n - n / 2;

    for (i = 0; i < mat->r; i++)
        for (j = 0; j < mat->c; j++)
        {
            ulong v;

            if (mode == 0)
                v = pos;
            else if (mode == 1)
                v = neg;
            else
                v = n_randint(state, 2) ? pos : neg;

            nmod_mat_entry(mat, i, j) = v % n;
        }
}

static ulong
random_modulus(flint_rand_t state)
{
    ulong bits, n;

    switch (n_randint(state, 8))
    {
        case 0:
            /* full range of bit sizes */
            bits = 1 + n_randint(state, 32);
            n = n_randbits(state, bits);
            break;
        case 1:
            /* just below a power of two: small 2^32 mod n, single fold
               round */
            bits = 2 + n_randint(state, 31);
            n = (UWORD(1) << bits) - 1 - n_randint(state, 16);
            break;
        case 2:
            /* just above a power of two: large 2^32 mod n, two rounds */
            bits = 1 + n_randint(state, 31);
            n = (UWORD(1) << bits) + 1 + n_randint(state, 16);
            break;
        case 3:
            /* the largest moduli: reduction after every 1-2 products */
            n = TWO32 - 1 - n_randint(state, 1000);
            break;
        case 4:
            /* 31-32 bits, arbitrary: exercises the n >= 2^31 correction
               in the final reduction */
            n = (UWORD(1) << 31) + n_randbits(state, 31);
            break;
        case 5:
            /* the range where this method is meant to be used */
            bits = 24 + n_randint(state, 9);
            n = n_randbits(state, bits);
            break;
        case 6:
            /* powers of two and neighbours */
            bits = 1 + n_randint(state, 32);
            n = (bits == 32) ? TWO32 - 1
                             : (UWORD(1) << bits) + n_randint(state, 3) - 1;
            break;
        default:
            /* tiny */
            n = 1 + n_randint(state, 20);
            break;
    }

    if (n == 0)
        n = 1;
    if (FLINT_BITS == 64 && n >= TWO32)
        n = TWO32 - 1;

    return n;
}

TEST_FUNCTION_START(nmod_mat_mul_u32, state)
{
    slong i, max_threads = 5;

    for (i = 0; i < 300 * flint_test_multiplier(); i++)
    {
        nmod_mat_t A, B, C, D;
        ulong modulus;
        slong m, k, n;
        int ret;

        /* mostly small shapes, sometimes a long inner dimension so that
           several k-blocks and many in-kernel folds are exercised */
        switch (n_randint(state, 4))
        {
            case 0:
                m = 1 + n_randint(state, 20);
                k = 1 + n_randint(state, 20);
                n = 1 + n_randint(state, 20);
                break;
            case 1:
                m = 1 + n_randint(state, 40);
                k = 1 + n_randint(state, 700);
                n = 1 + n_randint(state, 40);
                break;
            default:
                m = 1 + n_randint(state, 100);
                k = 1 + n_randint(state, 100);
                n = 1 + n_randint(state, 100);
                break;
        }

        if (n_randint(state, 50) == 0)
        {
            /* occasionally a dimension zero */
            switch (n_randint(state, 3))
            {
                case 0: m = 0; break;
                case 1: k = 0; break;
                default: n = 0; break;
            }
        }

        modulus = random_modulus(state);

        nmod_mat_init(A, m, k, modulus);
        nmod_mat_init(B, k, n, modulus);
        nmod_mat_init(C, m, n, modulus);
        nmod_mat_init(D, m, n, modulus);

        switch (n_randint(state, 4))
        {
            case 0:
                nmod_mat_randfull(A, state);
                nmod_mat_randfull(B, state);
                break;
            case 1:
                nmod_mat_randfull_extreme(A, state, n_randint(state, 3));
                nmod_mat_randfull_extreme(B, state, n_randint(state, 3));
                break;
            case 2:
                nmod_mat_randfull_extreme(A, state, n_randint(state, 3));
                nmod_mat_randfull(B, state);
                break;
            default:
                nmod_mat_randtest(A, state);
                nmod_mat_randtest(B, state);
                break;
        }

        /* garbage in C must be ignored */
        nmod_mat_randtest(C, state);

        flint_set_num_threads(n_randint(state, max_threads) + 1);

        ret = nmod_mat_mul_u32(C, A, B);

#if FLINT_BITS == 64
        if (!ret)
            TEST_FUNCTION_FAIL("mul_u32 should have worked\n"
                               "m: %wd, k: %wd, n: %wd, mod: %wu\n",
                               m, k, n, modulus);
#endif

        if (ret)
        {
            nmod_mat_mul_classical(D, A, B);

            if (!nmod_mat_equal(C, D))
                TEST_FUNCTION_FAIL("m: %wd, k: %wd, n: %wd, mod: %wu, "
                                   "threads: %d\n", m, k, n, modulus,
                                   flint_get_num_threads());
        }

        /* aliasing */
        if (ret && m == k && k == n && m > 0)
        {
            nmod_mat_set(D, A);
            nmod_mat_mul_u32(D, D, B);
            if (!nmod_mat_equal(C, D))
                TEST_FUNCTION_FAIL("aliasing C = A: m: %wd, mod: %wu\n",
                                   m, modulus);

            nmod_mat_set(D, B);
            nmod_mat_mul_u32(D, A, D);
            if (!nmod_mat_equal(C, D))
                TEST_FUNCTION_FAIL("aliasing C = B: m: %wd, mod: %wu\n",
                                   m, modulus);
        }

        /* the uint32 entry point, on the same data, with row strides
           larger than the row lengths */
        if (ret)
        {
            slong lda = k + n_randint(state, 3), ldb = n + n_randint(state, 3);
            slong ldc = n + n_randint(state, 3);
            uint32_t * a = flint_malloc((m * lda + 1) * sizeof(uint32_t));
            uint32_t * b = flint_malloc((k * ldb + 1) * sizeof(uint32_t));
            uint32_t * c = flint_malloc((m * ldc + 1) * sizeof(uint32_t));
            slong i, j;
            int ok = 1;

            for (i = 0; i < m; i++)
                for (j = 0; j < k; j++)
                    a[i * lda + j] = (uint32_t) nmod_mat_entry(A, i, j);
            for (i = 0; i < k; i++)
                for (j = 0; j < n; j++)
                    b[i * ldb + j] = (uint32_t) nmod_mat_entry(B, i, j);
            for (i = 0; i < m * ldc + 1; i++)
                c[i] = (uint32_t) n_randtest(state);

            if (!_nmod_mat_mul_u32(c, ldc, a, lda, b, ldb, m, k, n, C->mod))
                TEST_FUNCTION_FAIL("_nmod_mat_mul_u32 should have worked\n"
                                   "m: %wd, k: %wd, n: %wd, mod: %wu\n",
                                   m, k, n, modulus);

            for (i = 0; i < m; i++)
                for (j = 0; j < n; j++)
                    if (c[i * ldc + j] != nmod_mat_entry(D, i, j))
                        ok = 0;

            /* aliasing, square only */
            if (ok && m == k && k == n && m > 0 && lda == ldb)
            {
                for (i = 0; i < m * lda; i++)
                    c[i] = a[i];
                _nmod_mat_mul_u32(c, lda, c, lda, b, ldb, m, k, n, C->mod);
                for (i = 0; i < m; i++)
                    for (j = 0; j < n; j++)
                        if (c[i * lda + j] != nmod_mat_entry(D, i, j))
                            ok = 0;
            }

            if (!ok)
                TEST_FUNCTION_FAIL("uint32 entries: m: %wd, k: %wd, n: %wd, "
                                   "mod: %wu, threads: %d\n", m, k, n,
                                   modulus, flint_get_num_threads());

            flint_free(a);
            flint_free(b);
            flint_free(c);
        }

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
        nmod_mat_clear(D);
    }

#if FLINT_BITS == 64
    /* a modulus of 32 bits or more must be declined */
    {
        nmod_mat_t A, B, C;
        ulong modulus = TWO32;

        nmod_mat_init(A, 5, 5, modulus);
        nmod_mat_init(B, 5, 5, modulus);
        nmod_mat_init(C, 5, 5, modulus);

        if (nmod_mat_mul_u32(C, A, B) != 0)
            TEST_FUNCTION_FAIL("mul_u32 accepted a modulus of 2^32\n");

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
    }
#endif

    TEST_FUNCTION_END(state);
}
