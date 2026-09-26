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
# define TWO52 (UWORD(1) << 52)
# define TWO50 (UWORD(1) << 50)
#else
# define TWO52 UWORD(0)
# define TWO50 UWORD(0)
#endif

/* all entries n - 1: every product is the largest possible one */
static void
randfull_max_k52(nmod_mat_t mat, flint_rand_t state, int mixed)
{
    slong i, j;
    ulong n = mat->mod.n;

    for (i = 0; i < mat->r; i++)
        for (j = 0; j < mat->c; j++)
            nmod_mat_entry(mat, i, j) =
                (mixed && n_randint(state, 2)) ? n_randint(state, n) : n - 1;
}

static ulong
random_modulus_k52(flint_rand_t state)
{
    ulong bits, n;

    switch (n_randint(state, 8))
    {
        case 0:
            /* full range of bit sizes */
            bits = 1 + n_randint(state, 52);
            n = n_randbits(state, bits);
            break;
        case 1:
            /* just below a power of two: small 2^52 mod n */
            bits = 2 + n_randint(state, 51);
            n = (UWORD(1) << bits) - 1 - n_randint(state, 16);
            break;
        case 2:
            /* just above a power of two */
            bits = 1 + n_randint(state, 51);
            n = (UWORD(1) << bits) + 1 + n_randint(state, 16);
            break;
        case 3:
            /* the largest moduli, including 2^52 itself */
            n = TWO52 - n_randint(state, 1000);
            break;
        case 4:
            /* around 2^27, where the high limb appears */
            n = (UWORD(1) << 27) - 8 + n_randint(state, 17);
            break;
        case 5:
            /* the range where this method is meant to be used */
            bits = 31 + n_randint(state, 22);
            n = n_randbits(state, bits);
            break;
        case 6:
            /* powers of two */
            bits = 1 + n_randint(state, 52);
            n = UWORD(1) << bits;
            break;
        default:
            /* tiny */
            n = 1 + n_randint(state, 20);
            break;
    }

    if (n == 0)
        n = 1;
    if (FLINT_BITS == 64 && n > TWO52)
        n = TWO52;

    return n;
}

/*
    The product A*B = D again, through windows at random offsets inside
    larger matrices (row strides larger than the row lengths, most of the
    time), with random entries around them: returns 0 if the result is
    wrong or if an entry of the parent of C outside its window was
    touched.
*/
static int
mul_window_check_k52(const nmod_mat_t A, const nmod_mat_t B,
                     const nmod_mat_t D, flint_rand_t state)
{
    nmod_mat_t PA, PB, PC, PC0, WA, WB, WC;
    slong m = A->r, k = A->c, n = B->c, i, j;
    slong ra = n_randint(state, 4), ca = n_randint(state, 4);
    slong rb = n_randint(state, 4), cb = n_randint(state, 4);
    slong rc = n_randint(state, 4), cc = n_randint(state, 4);
    ulong modulus = A->mod.n;
    int ok;

    nmod_mat_init(PA, ra + m + n_randint(state, 4),
                      ca + k + n_randint(state, 4), modulus);
    nmod_mat_init(PB, rb + k + n_randint(state, 4),
                      cb + n + n_randint(state, 4), modulus);
    nmod_mat_init(PC, rc + m + n_randint(state, 4),
                      cc + n + n_randint(state, 4), modulus);
    nmod_mat_randtest(PA, state);
    nmod_mat_randtest(PB, state);
    nmod_mat_randtest(PC, state);

    nmod_mat_window_init(WA, PA, ra, ca, ra + m, ca + k);
    nmod_mat_window_init(WB, PB, rb, cb, rb + k, cb + n);
    nmod_mat_window_init(WC, PC, rc, cc, rc + m, cc + n);
    nmod_mat_set(WA, A);
    nmod_mat_set(WB, B);
    nmod_mat_init_set(PC0, PC);

    ok = nmod_mat_mul_k52(WC, WA, WB) && nmod_mat_equal(WC, D);

    for (i = 0; i < PC->r; i++)
        for (j = 0; j < PC->c; j++)
            if ((i < rc || i >= rc + m || j < cc || j >= cc + n)
                    && nmod_mat_entry(PC, i, j) != nmod_mat_entry(PC0, i, j))
                ok = 0;

    nmod_mat_window_clear(WA);
    nmod_mat_window_clear(WB);
    nmod_mat_window_clear(WC);
    nmod_mat_clear(PA);
    nmod_mat_clear(PB);
    nmod_mat_clear(PC);
    nmod_mat_clear(PC0);

    return ok;
}

TEST_FUNCTION_START(nmod_mat_mul_k52, state)
{
    slong i, max_threads = 5;
    int available = 0;

    for (i = 0; i < 300 * flint_test_multiplier(); i++)
    {
        nmod_mat_t A, B, C, D;
        ulong modulus;
        slong m, k, n;
        int ret;

        switch (n_randint(state, 4))
        {
            case 0:
                m = 1 + n_randint(state, 20);
                k = 1 + n_randint(state, 20);
                n = 1 + n_randint(state, 20);
                break;
            case 1:
                /* several k-blocks */
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
            switch (n_randint(state, 3))
            {
                case 0: m = 0; break;
                case 1: k = 0; break;
                default: n = 0; break;
            }
        }

        modulus = random_modulus_k52(state);

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
                randfull_max_k52(A, state, 0);
                randfull_max_k52(B, state, 0);
                break;
            case 2:
                randfull_max_k52(A, state, 1);
                randfull_max_k52(B, state, 1);
                break;
            default:
                nmod_mat_randtest(A, state);
                nmod_mat_randtest(B, state);
                break;
        }

        /* garbage in C must be ignored */
        nmod_mat_randtest(C, state);

        flint_set_num_threads(n_randint(state, max_threads) + 1);

        ret = nmod_mat_mul_k52(C, A, B);

        if (i == 0)
            available = ret;
        else if (ret != available)
            TEST_FUNCTION_FAIL("availability changed: %d then %d\n"
                               "m: %wd, k: %wd, n: %wd, mod: %wu\n",
                               available, ret, m, k, n, modulus);

        if (ret)
        {
            nmod_mat_mul_classical(D, A, B);

            if (!nmod_mat_equal(C, D))
                TEST_FUNCTION_FAIL("m: %wd, k: %wd, n: %wd, mod: %wu, "
                                   "threads: %d\n", m, k, n, modulus,
                                   flint_get_num_threads());

            /* the same product on windows of larger matrices */
            if (n_randint(state, 2)
                    && !mul_window_check_k52(A, B, D, state))
                TEST_FUNCTION_FAIL("windows: m: %wd, k: %wd, n: %wd, "
                                   "mod: %wu, threads: %d\n", m, k, n,
                                   modulus, flint_get_num_threads());

            /* aliasing */
            if (m == k && k == n && m > 0)
            {
                nmod_mat_set(D, A);
                nmod_mat_mul_k52(D, D, B);
                if (!nmod_mat_equal(C, D))
                    TEST_FUNCTION_FAIL("aliasing C = A: m: %wd, mod: %wu\n",
                                       m, modulus);

                nmod_mat_set(D, B);
                nmod_mat_mul_k52(D, A, D);
                if (!nmod_mat_equal(C, D))
                    TEST_FUNCTION_FAIL("aliasing C = B: m: %wd, mod: %wu\n",
                                       m, modulus);
            }
        }

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
        nmod_mat_clear(D);
    }

#if FLINT_BITS == 64
    /* a modulus out of range must be declined */
    {
        nmod_mat_t A, B, C;
        ulong modulus = TWO52 + 1;

        nmod_mat_init(A, 5, 5, modulus);
        nmod_mat_init(B, 5, 5, modulus);
        nmod_mat_init(C, 5, 5, modulus);

        if (nmod_mat_mul_k52(C, A, B) != 0)
            TEST_FUNCTION_FAIL("mul_k52 accepted a modulus of 2^52 + 1\n");

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
    }
#endif

    TEST_FUNCTION_END(state);
}
