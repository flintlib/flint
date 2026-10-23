/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5 and Opus 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod.h"
#include "nmod_mat.h"

/* runtime Strassen cutoffs of the packed engines, deliberately not in
   the header: at the shipped defaults the recursion needs matrices in
   the thousands, so full coverage lowers them here */
FLINT_DLL extern slong nmod_mat_gf2_strassen_cutoff;
FLINT_DLL extern slong nmod_mat_gf3_strassen_cutoff;

/*
    The engines dispatch on modulus and shape, so the checks below aim
    every regime explicitly rather than trusting random sizes to reach
    them: the shuffle basecase and its tiny-size takeover of moduli 2
    and 3; the bit-packed serial paths above that takeover, including
    the Strassen recursion via the runtime cutoff hook, since the
    shipped cutoffs need matrices in the thousands; the threaded packed
    driver, whose tile grid needs a 2048 side; the sgemm basecase for
    moduli above 15, in its single-block, accumulated multi-block and
    (very large k) accumulator-fold forms; the modular add/sub pair
    above modulus 128, reached through the byte Strassen; the squaring
    shortcuts, which alias A and B; and the ulong <-> byte conversions
    of the public wrapper, against the raw byte entry.
*/

static void
u8_image(uint8_t * U, const nmod_mat_t M, slong stride)
{
    slong i, j;

    for (i = 0; i < M->r; i++)
        for (j = 0; j < M->c; j++)
            U[i * stride + j] = (uint8_t) nmod_mat_entry(M, i, j);
}

/* check both entries against classical on given operands */
static void
check_one(nmod_mat_t A, nmod_mat_t B, flint_rand_t state,
          const char * ctx)
{
    slong m = A->r, k = A->c, n = B->c, i, j;
    slong As = k + n_randint(state, 3);
    slong Bs = n + n_randint(state, 3);
    slong Cs = n + n_randint(state, 3);
    int same = (A == B);
    nmod_mat_t C, D;
    uint8_t * Au;
    uint8_t * Bu;
    uint8_t * Cu;

    nmod_mat_init(C, m, n, A->mod.n);
    nmod_mat_init(D, m, n, A->mod.n);
    nmod_mat_randtest(C, state);

    nmod_mat_mul_classical(D, A, B);

    nmod_mat_mul_u8(C, A, B);

    if (!nmod_mat_equal(C, D))
        do { flint_printf("FAIL: " "wrapper wrong (%s): "
                           "m: %wd, k: %wd, n: %wd, mod: %wu\n",
                           ctx, m, k, n, A->mod.n); flint_abort(); } while (0);

    /* raw byte entry, slack strides; squaring keeps A and B aliased */
    if (same)
        Bs = As;

    Au = flint_malloc(m * As + (same ? 0 : k * Bs) + m * Cs);
    Bu = same ? Au : Au + m * As;
    Cu = Bu + k * Bs;

    u8_image(Au, A, As);
    if (!same)
        u8_image(Bu, B, Bs);
    memset(Cu, 0xee, m * Cs);

    _nmod_mat_mul_u8(Cu, Cs, Au, As, Bu, Bs, m, k, n, A->mod);

    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            if (Cu[i * Cs + j] != nmod_mat_entry(D, i, j))
                do { flint_printf("FAIL: " "raw entry wrong (%s): "
                                   "m: %wd, k: %wd, n: %wd, mod: %wu\n",
                                   ctx, m, k, n, A->mod.n); flint_abort(); } while (0);

    flint_free(Au);
    nmod_mat_clear(C);
    nmod_mat_clear(D);
}

static void
check_square(slong d, ulong p, flint_rand_t state, const char * ctx)
{
    nmod_mat_t A;

    nmod_mat_init(A, d, d, p);
    nmod_mat_randfull(A, state);
    check_one(A, A, state, ctx);
    nmod_mat_clear(A);
}

static void
check_dims(slong m, slong k, slong n, ulong p, flint_rand_t state,
           const char * ctx)
{
    nmod_mat_t A, B;

    nmod_mat_init(A, m, k, p);
    nmod_mat_init(B, k, n, p);
    nmod_mat_randfull(A, state);
    nmod_mat_randfull(B, state);
    check_one(A, B, state, ctx);
    nmod_mat_clear(A);
    nmod_mat_clear(B);
}

/*
    Aliased raw-entry calls (C == A, C == B, and the aliased square
    C == A == B) against the unaliased product; d x d so both
    aliasings are dimensionally possible.
*/
static void
check_aliased(slong d, ulong p, flint_rand_t state, const char * ctx)
{
    nmod_t mod;
    uint8_t * M;
    uint8_t * N;
    uint8_t * ref;
    slong i, mode;

    nmod_init(&mod, p);
    M = flint_malloc(3 * d * d);
    N = M + d * d;
    ref = N + d * d;

    for (mode = 0; mode < 3; mode++)
    {
        for (i = 0; i < d * d; i++)
        {
            M[i] = (uint8_t) n_randint(state, p);
            N[i] = (uint8_t) n_randint(state, p);
        }

        if (mode == 0)      /* C == A */
        {
            _nmod_mat_mul_u8(ref, d, M, d, N, d, d, d, d, mod);
            _nmod_mat_mul_u8(M, d, M, d, N, d, d, d, d, mod);
        }
        else if (mode == 1) /* C == B */
        {
            _nmod_mat_mul_u8(ref, d, M, d, N, d, d, d, d, mod);
            _nmod_mat_mul_u8(N, d, M, d, N, d, d, d, d, mod);
        }
        else                /* C == A == B */
        {
            _nmod_mat_mul_u8(ref, d, M, d, M, d, d, d, d, mod);
            _nmod_mat_mul_u8(M, d, M, d, M, d, d, d, d, mod);
        }

        if (memcmp(mode == 1 ? N : M, ref, d * d) != 0)
        { flint_printf("FAIL: aliased (%s): d: %wd, mod: %wu, "
                       "mode: %wd\n", ctx, d, p, mode); flint_abort(); }
    }

    flint_free(M);
}

TEST_FUNCTION_START(nmod_mat_mul_u8, state)
{
    slong i;

    /* random shapes across the full modulus range, dispatch's choice */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        slong m = n_randint(state, 90) + 1;
        slong k = n_randint(state, 90) + 1;
        slong n = n_randint(state, 90) + 1;
        ulong p = 2 + n_randint(state, 254);

        flint_set_num_threads(n_randint(state, 4) + 1);
        check_dims(m, k, n, p, state, "random small");
    }

    flint_set_num_threads(1);

    /* the tiny-size takeover boundaries for the packed moduli, and
       degenerate inner dimensions on both sides of them */
    {
        slong d2[] = { 1, 63, 64, 65 };
        slong d3[] = { 1, 191, 192, 193 };

        for (i = 0; i < 4; i++)
        {
            check_dims(d2[i], 1 + n_randint(state, 200), 40, 2, state,
                       "gf2 boundary");
            check_dims(d3[i], 1 + n_randint(state, 200), 40, 3, state,
                       "gf3 boundary");
            check_square(d2[i], 2, state, "gf2 boundary squaring");
            check_square(d3[i], 3, state, "gf3 boundary squaring");
        }

        check_dims(300, 1, 300, 2, state, "gf2 k = 1");
        check_dims(300, 1, 300, 3, state, "gf3 k = 1");
    }

    /* serial packed leaves above the takeover, ragged and squared */
    for (i = 0; i < 2 * flint_test_multiplier(); i++)
    {
        check_dims(65 + n_randint(state, 300), 1 + n_randint(state, 300),
                   1 + n_randint(state, 300), 2, state, "gf2 serial");
        check_dims(193 + n_randint(state, 300), 1 + n_randint(state, 300),
                   1 + n_randint(state, 300), 3, state, "gf3 serial");
        check_square(65 + n_randint(state, 300), 2, state, "gf2 squaring");
        check_square(193 + n_randint(state, 300), 3, state, "gf3 squaring");
    }

    /* packed Strassen: unreachable at the shipped cutoffs, so lower
       them; ragged sizes force the padding, squares share planes only
       at the leaves */
    {
        slong save2 = nmod_mat_gf2_strassen_cutoff;
        slong save3 = nmod_mat_gf3_strassen_cutoff;

        nmod_mat_gf2_strassen_cutoff = 32;
        nmod_mat_gf3_strassen_cutoff = 64;

    for (i = 0; i < 2 * flint_test_multiplier(); i++)
    {
        check_dims(130 + n_randint(state, 130), 130 + n_randint(state, 130),
                   130 + n_randint(state, 130), 2, state, "gf2 strassen");
        check_dims(260 + n_randint(state, 130), 260 + n_randint(state, 130),
                   260 + n_randint(state, 130), 3, state, "gf3 strassen");
        check_square(131 + n_randint(state, 130), 2, state,
                     "gf2 strassen square");
    }
        nmod_mat_gf2_strassen_cutoff = save2;
        nmod_mat_gf3_strassen_cutoff = save3;
    }

    /* threaded packed driver: the grid wants a 2048 side, k stays
       small to keep this cheap; compare threaded against one thread */
    {
        slong ks[] = { 3, 70 };
        ulong ps[] = { 2, 3 };
        slong j;

        for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
        {
            nmod_mat_t A, B, C1, C2;

            nmod_mat_init(A, 2048, ks[i], ps[j]);
            nmod_mat_init(B, ks[i], 2048, ps[j]);
            nmod_mat_init(C1, 2048, 2048, ps[j]);
            nmod_mat_init(C2, 2048, 2048, ps[j]);
            nmod_mat_randfull(A, state);
            nmod_mat_randfull(B, state);

            flint_set_num_threads(1);
            nmod_mat_mul_u8(C1, A, B);

            flint_set_num_threads(2 + (int) n_randint(state, 3));
            nmod_mat_mul_u8(C2, A, B);

            if (!nmod_mat_equal(C1, C2))
                TEST_FUNCTION_FAIL("threaded disagrees with serial: "
                                   "k: %wd, mod: %wu\n", ks[i], ps[j]);

            nmod_mat_clear(A);
            nmod_mat_clear(B);
            nmod_mat_clear(C1);
            nmod_mat_clear(C2);
        }

        /* threaded squaring shares the packed planes; at this size the
           serial u8 result (verified against classical at small sizes
           and via the cutoff hook) is the reference, a classical
           product being far too slow here */
        for (j = 0; j < 2; j++)
        {
            nmod_mat_t A, C1, C2;

            nmod_mat_init(A, 2048, 2048, ps[j]);
            nmod_mat_init(C1, 2048, 2048, ps[j]);
            nmod_mat_init(C2, 2048, 2048, ps[j]);
            nmod_mat_randfull(A, state);

            flint_set_num_threads(1);
            nmod_mat_mul_u8(C1, A, A);
            flint_set_num_threads(3);
            nmod_mat_mul_u8(C2, A, A);

            if (!nmod_mat_equal(C1, C2))
                TEST_FUNCTION_FAIL("threaded squaring disagrees: "
                                   "mod: %wu\n", ps[j]);

            nmod_mat_clear(A);
            nmod_mat_clear(C1);
            nmod_mat_clear(C2);
        }
        flint_set_num_threads(1);
    }

    /* byte Strassen and the modular add/sub above modulus 128 */
    for (i = 0; i < 2 * flint_test_multiplier(); i++)
    {
        ulong p = 129 + n_randint(state, 127);

        check_dims(200 + n_randint(state, 120), 200 + n_randint(state, 120),
                   200 + n_randint(state, 120), p, state, "high add/sub");
        check_square(200 + n_randint(state, 120), p, state,
                     "high squaring");
    }

    /* sgemm basecase: single block, accumulated blocks, and -- with k
       past 128 blocks at modulus 251 -- the accumulator fold */
    for (i = 0; i < 2 * flint_test_multiplier(); i++)
    {
        check_dims(30 + n_randint(state, 40), 40 + n_randint(state, 40),
                   30 + n_randint(state, 40), 17 + n_randint(state, 100),
                   state, "sgemm single block");
        check_dims(20 + n_randint(state, 20), 900 + n_randint(state, 300),
                   20 + n_randint(state, 20), 200 + n_randint(state, 56),
                   state, "sgemm multi block");
        check_square(60 + n_randint(state, 60), 17 + n_randint(state, 239),
                     state, "sgemm squaring");
    }
    check_dims(4, 35000, 4, 251, state, "sgemm accumulator fold");

    /* aliased raw-entry calls across the engines: shuffle basecase
       and byte Strassen (buffered), sgemm, and -- exempt from the
       buffer -- the packed moduli, serial and threaded */
    for (i = 0; i < flint_test_multiplier(); i++)
    {
        check_aliased(40 + n_randint(state, 60),
                      5 + n_randint(state, 11), state, "shuffle");
        check_aliased(200 + n_randint(state, 120),
                      5 + n_randint(state, 11), state, "byte strassen");
        check_aliased(200 + n_randint(state, 120),
                      17 + n_randint(state, 239), state, "sgemm");
        check_aliased(100 + n_randint(state, 200), 2, state, "gf2");
        check_aliased(250 + n_randint(state, 200), 3, state, "gf3");
    }
    flint_set_num_threads(3);
    check_aliased(2048, 2, state, "gf2 threaded aliased");
    flint_set_num_threads(1);

    /* the flat shared-operand gemm engages above 256 x 256 with
       several threads when one exact sgemm covers the product; the
       larger-modulus shapes stay on the tiled path and are checked
       threaded here too */
    flint_set_num_threads(2 + (int) n_randint(state, 3));
    for (i = 0; i < flint_test_multiplier(); i++)
    {
        check_dims(260 + n_randint(state, 60), 100 + n_randint(state, 400),
                   260 + n_randint(state, 60), 17 + n_randint(state, 40),
                   state, "flat gemm");
        check_square(260 + n_randint(state, 60),
                     17 + n_randint(state, 40), state,
                     "flat gemm squaring");
        check_dims(260 + n_randint(state, 60), 600 + n_randint(state, 400),
                   260 + n_randint(state, 60), 200 + n_randint(state, 56),
                   state, "tiled large modulus threaded");
    }
    flint_set_num_threads(1);

    TEST_FUNCTION_END(state);
}
