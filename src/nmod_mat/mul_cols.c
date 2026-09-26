/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "nmod_mat/impl.h"

/*
    C = A * B for B with few columns (2 <= r <= NMOD_MAT_MUL_COLS_MAX): the
    columns of B are transposed once into contiguous, zero-padded vectors, and
    each row of A is then multiplied with all of them at once, as r dot
    products sharing the loads of the row, with the accumulation strategies of
    the SIMD dot products of nmod_vec:

    - AVX512-IFMA, moduli up to 2^52 (dot_u52.c): low and high 52-bit
      halves of the products, u accumulator pairs per column to hide the
      latency (u = 4, 2, 1 for r <= 2, <= 4, > 4), summed into two-limb
      totals per column every 4032 terms;
    - AVX512-IFMA, moduli above 2^52 (dot_u64.c): the four products of
      the 32-bit halves of the entries, eight accumulators per column,
      two columns per pass over A, three-limb totals every 2000 terms;
    - otherwise a vector backend of mul_fp_vec.h, moduli below 2^50
      (dot_fp50.c): products reduced in double precision as they are
      formed, u accumulators per column, reduced every four terms.

    Compared with the kernels, nothing is packed and no tile is padded to
    NR columns; compared with the classical code, the columns of B are
    read contiguously and A is streamed once instead of once per column.
*/

#if FLINT_BITS == 64 && NMOD_MAT_HAVE_MUL_U52

#include <immintrin.h>

#define COLS_U52_MAXACC 8

/*
    One row a of A (k entries) against the r transposed columns bt
    (stride kp, zero-padded to kp entries): c[j] = sum_i a[i] bt[j][i].
    r and u (accumulator pairs per column, r * u <= 8) are constants at
    each call site. The lanes are summed every chunk of iterations, when
    the u pairs of a column hold at most 508 halves per lane (below 2^61).
*/
FLINT_FORCE_INLINE void
cols_u52_row(ulong * c, const ulong * a, const ulong * bt, slong kp,
             slong k, int r, int u, nmod_t mod)
{
    const slong chunk = (504 / u) * (8 * u);   /* terms per chunk */
    __m512i lo[COLS_U52_MAXACC], hi[COLS_U52_MAXACC], x, y;
    ulong t1[8], t0[8], sl, sh;
    slong i = 0, stop;
    int j, v;

    for (j = 0; j < r; j++)
        t1[j] = t0[j] = 0;
    for (j = 0; j < r * u; j++)
        lo[j] = hi[j] = _mm512_setzero_si512();

    while (i < k)
    {
        stop = FLINT_MIN(k, i + chunk);

        for ( ; i + 8 * u <= stop; i += 8 * u)
        {
            for (v = 0; v < u; v++)
            {
                x = _mm512_loadu_si512((const void *) (a + i + 8 * v));
                for (j = 0; j < r; j++)
                {
                    y = _mm512_loadu_si512((const void *) (bt + j * kp + i + 8 * v));
                    lo[j * u + v] = _mm512_madd52lo_epu64(lo[j * u + v], x, y);
                    hi[j * u + v] = _mm512_madd52hi_epu64(hi[j * u + v], x, y);
                }
            }
        }

        if (i < stop)   /* fewer than u vectors remain: the last chunk */
        {
            for (v = 0; i < stop; i += 8, v++)
            {
                const __mmask8 m = (stop - i >= 8) ? 0xFF
                                   : (__mmask8) ((1 << (stop - i)) - 1);
                x = _mm512_maskz_loadu_epi64(m, (const void *) (a + i));
                for (j = 0; j < r; j++)
                {
                    y = _mm512_loadu_si512((const void *) (bt + j * kp + i));
                    lo[j * u + v] = _mm512_madd52lo_epu64(lo[j * u + v], x, y);
                    hi[j * u + v] = _mm512_madd52hi_epu64(hi[j * u + v], x, y);
                }
            }
        }

        for (j = 0; j < r; j++)
        {
            x = lo[j * u];
            y = hi[j * u];
            for (v = 1; v < u; v++)
            {
                x = _mm512_add_epi64(x, lo[j * u + v]);
                y = _mm512_add_epi64(y, hi[j * u + v]);
            }
            sl = _mm512_reduce_add_epi64(x);
            sh = _mm512_reduce_add_epi64(y);
            add_ssaaaa(t1[j], t0[j], t1[j], t0[j], sh >> 12, sh << 52);
            add_ssaaaa(t1[j], t0[j], t1[j], t0[j], UWORD(0), sl);
        }
        if (i < k)
            for (j = 0; j < r * u; j++)
                lo[j] = hi[j] = _mm512_setzero_si512();
    }

    for (j = 0; j < r; j++)
        NMOD2_RED2(c[j], t1[j], t0[j], mod);
}

static void
cols_u52(nmod_mat_t C, const nmod_mat_t A, const ulong * bt, slong kp, int r)
{
    const slong m = A->r, k = A->c;
    const nmod_t mod = A->mod;
    slong i;

    for (i = 0; i < m; i++)
    {
        ulong * c = nmod_mat_entry_ptr(C, i, 0);
        const ulong * a = nmod_mat_entry_ptr(A, i, 0);

        switch (r)
        {
            case 1: cols_u52_row(c, a, bt, kp, k, 1, 4, mod); break;
            case 2: cols_u52_row(c, a, bt, kp, k, 2, 4, mod); break;
            case 3: cols_u52_row(c, a, bt, kp, k, 3, 2, mod); break;
            case 4: cols_u52_row(c, a, bt, kp, k, 4, 2, mod); break;
            case 5: cols_u52_row(c, a, bt, kp, k, 5, 1, mod); break;
            case 6: cols_u52_row(c, a, bt, kp, k, 6, 1, mod); break;
            case 7: cols_u52_row(c, a, bt, kp, k, 7, 1, mod); break;
            case 8: cols_u52_row(c, a, bt, kp, k, 8, 1, mod); break;
            default: FLINT_UNREACHABLE;
        }
    }
}

/*
    Moduli above 2^52: as in dot_u64.c, x y = x0 y0 + 2^32 (x0 y1 + x1 y0)
    + 2^64 x1 y1 on the 32-bit halves, eight accumulators per column, two
    columns per pass over A (16 accumulators; three spill); every 250
    iterations the lanes (below 2^60) go to three-limb totals.
*/
#define COLS_U64_CHUNK (250 * 8)

#define COLS_U64_ADD_SHIFTED(t2, t1, t0, s, shift)                      \
    add_sssaaaaaa(t2, t1, t0, t2, t1, t0,                               \
                  UWORD(0), (s) >> (64 - (shift)), (s) << (shift))
#define COLS_U64_ADD_SHIFTED64(t2, t1, t0, s, shift)                    \
    add_sssaaaaaa(t2, t1, t0, t2, t1, t0,                               \
                  (s) >> (64 - (shift)), (s) << (shift), UWORD(0))

FLINT_FORCE_INLINE void
cols_u64_row(ulong * c, const ulong * a, const ulong * bt, slong kp,
             slong k, int r, nmod_t mod)
{
    const __m512i m32 = _mm512_set1_epi64(UWORD(0xFFFFFFFF));
    __m512i l00[2], h00[2], l01[2], h01[2], l10[2], h10[2], l11[2], h11[2];
    __m512i x, x0, x1, y, y0, y1;
    ulong t2[2], t1[2], t0[2], s;
    slong i = 0, stop;
    int j;

    for (j = 0; j < r; j++)
    {
        t2[j] = t1[j] = t0[j] = 0;
        l00[j] = h00[j] = l01[j] = h01[j] = _mm512_setzero_si512();
        l10[j] = h10[j] = l11[j] = h11[j] = _mm512_setzero_si512();
    }

    while (i < k)
    {
        stop = FLINT_MIN(k, i + COLS_U64_CHUNK);

        for ( ; i < stop; i += 8)
        {
            const __mmask8 m = (stop - i >= 8) ? 0xFF
                               : (__mmask8) ((1 << (stop - i)) - 1);
            x = _mm512_maskz_loadu_epi64(m, (const void *) (a + i));
            x0 = _mm512_and_si512(x, m32);
            x1 = _mm512_srli_epi64(x, 32);
            for (j = 0; j < r; j++)
            {
                y = _mm512_loadu_si512((const void *) (bt + j * kp + i));
                y0 = _mm512_and_si512(y, m32);
                y1 = _mm512_srli_epi64(y, 32);
                l00[j] = _mm512_madd52lo_epu64(l00[j], x0, y0);
                h00[j] = _mm512_madd52hi_epu64(h00[j], x0, y0);
                l01[j] = _mm512_madd52lo_epu64(l01[j], x0, y1);
                h01[j] = _mm512_madd52hi_epu64(h01[j], x0, y1);
                l10[j] = _mm512_madd52lo_epu64(l10[j], x1, y0);
                h10[j] = _mm512_madd52hi_epu64(h10[j], x1, y0);
                l11[j] = _mm512_madd52lo_epu64(l11[j], x1, y1);
                h11[j] = _mm512_madd52hi_epu64(h11[j], x1, y1);
            }
        }

        for (j = 0; j < r; j++)
        {
            s = _mm512_reduce_add_epi64(l00[j]);
            add_sssaaaaaa(t2[j], t1[j], t0[j], t2[j], t1[j], t0[j], UWORD(0), UWORD(0), s);
            s = _mm512_reduce_add_epi64(h00[j]);
            COLS_U64_ADD_SHIFTED(t2[j], t1[j], t0[j], s, 52);
            s = _mm512_reduce_add_epi64(_mm512_add_epi64(l01[j], l10[j]));
            COLS_U64_ADD_SHIFTED(t2[j], t1[j], t0[j], s, 32);
            s = _mm512_reduce_add_epi64(_mm512_add_epi64(h01[j], h10[j]));
            COLS_U64_ADD_SHIFTED64(t2[j], t1[j], t0[j], s, 20);
            s = _mm512_reduce_add_epi64(l11[j]);
            add_sssaaaaaa(t2[j], t1[j], t0[j], t2[j], t1[j], t0[j], UWORD(0), s, UWORD(0));
            s = _mm512_reduce_add_epi64(h11[j]);
            COLS_U64_ADD_SHIFTED64(t2[j], t1[j], t0[j], s, 52);
            if (i < k)
            {
                l00[j] = h00[j] = l01[j] = h01[j] = _mm512_setzero_si512();
                l10[j] = h10[j] = l11[j] = h11[j] = _mm512_setzero_si512();
            }
        }
    }

    for (j = 0; j < r; j++)
    {
        NMOD_RED(t2[j], t2[j], mod);
        NMOD_RED3(c[j], t2[j], t1[j], t0[j], mod);
    }
}

static void
cols_u64(nmod_mat_t C, const nmod_mat_t A, const ulong * bt, slong kp, int r)
{
    const slong m = A->r, k = A->c;
    const nmod_t mod = A->mod;
    slong i;
    int j0;

    for (j0 = 0; j0 < r; j0 += 2)
    {
        const int rr = FLINT_MIN(r - j0, 2);
        const ulong * b = bt + j0 * kp;

        for (i = 0; i < m; i++)
        {
            ulong * c = nmod_mat_entry_ptr(C, i, j0);
            const ulong * a = nmod_mat_entry_ptr(A, i, 0);

            switch (rr)
            {
                case 1: cols_u64_row(c, a, b, kp, k, 1, mod); break;
                case 2: cols_u64_row(c, a, b, kp, k, 2, mod); break;
                default: FLINT_UNREACHABLE;
            }
        }
    }
}

#endif  /* IFMA */

#if FLINT_BITS == 64 && !NMOD_MAT_HAVE_MUL_U52

#include "nmod_mat/mul_fp_vec.h"

#if NMOD_MAT_HAVE_FPV
#define COLS_HAVE_FP50 1

/*
    Moduli below 2^50 in double precision (dot_fp50.c): every product is
    in (-9/8 n, 9/8 n), u accumulators per column (r * u <= 8) take four
    terms each before being brought back to (-0.51 n, 0.51 n). The tail of
    fewer than FPV_VL terms is scalar. The horizontal sum of a column is an
    integer offset by 8 n (nonnegative, below 2^54) plus the tail.
*/
FLINT_FORCE_INLINE void
cols_fp50_row(ulong * c, const ulong * a, const ulong * bt, slong kp,
              slong k, int r, int u, nmod_t mod)
{
    const ulong n = mod.n;
    const fpv nv = fpv_set1((double) n);
    const fpv ninv = fpv_set1(1.0 / (double) n);
    fpv acc[8], x, s;
    slong i = 0;
    int j, v, cnt = 0;

    for (j = 0; j < r * u; j++)
        acc[j] = fpv_zero();

    for ( ; i + FPV_VL * u <= k; i += FPV_VL * u)
    {
        for (v = 0; v < u; v++)
        {
            x = fpv_load_u64(a + i + FPV_VL * v);
            for (j = 0; j < r; j++)
                acc[j * u + v] = fpv_add(acc[j * u + v], fpv_mulmod(x,
                        fpv_load_u64(bt + j * kp + i + FPV_VL * v), nv, ninv));
        }
        if (++cnt == 4)
        {
            for (j = 0; j < r * u; j++)
                acc[j] = fpv_reduce_pm1n(acc[j], nv, ninv);
            cnt = 0;
        }
    }
    if (cnt != 0)
        for (j = 0; j < r * u; j++)
            acc[j] = fpv_reduce_pm1n(acc[j], nv, ninv);

    /* remaining full vectors (fewer than u), reduced every time */
    for ( ; i + FPV_VL <= k; i += FPV_VL)
    {
        x = fpv_load_u64(a + i);
        for (j = 0; j < r; j++)
            acc[j * u] = fpv_reduce_pm1n(fpv_add(acc[j * u], fpv_mulmod(x,
                    fpv_load_u64(bt + j * kp + i), nv, ninv)), nv, ninv);
    }

    for (j = 0; j < r; j++)
    {
        double buf[FPV_VL], sum = 0.0;
        ulong t1 = 0, t0;
        slong p;
        int l;

        s = acc[j * u];
        for (v = 1; v < u; v++)
            s = fpv_reduce_pm1n(fpv_add(s, acc[j * u + v]), nv, ninv);
        fpv_storeu(buf, s);
        for (l = 0; l < FPV_VL; l++)
            sum += buf[l];
        t0 = (ulong) ((slong) sum + 8 * (slong) n);

        for (p = i; p < k; p++)
        {
            ulong s1, s0;
            umul_ppmm(s1, s0, a[p], bt[j * kp + p]);
            add_ssaaaa(t1, t0, t1, t0, s1, s0);
        }
        NMOD2_RED2(c[j], t1, t0, mod);
    }
}

static void
cols_fp50(nmod_mat_t C, const nmod_mat_t A, const ulong * bt, slong kp, int r)
{
    const slong m = A->r, k = A->c;
    const nmod_t mod = A->mod;
    slong i;

    for (i = 0; i < m; i++)
    {
        ulong * c = nmod_mat_entry_ptr(C, i, 0);
        const ulong * a = nmod_mat_entry_ptr(A, i, 0);

        switch (r)
        {
            case 1: cols_fp50_row(c, a, bt, kp, k, 1, 4, mod); break;
            case 2: cols_fp50_row(c, a, bt, kp, k, 2, 4, mod); break;
            case 3: cols_fp50_row(c, a, bt, kp, k, 3, 2, mod); break;
            case 4: cols_fp50_row(c, a, bt, kp, k, 4, 2, mod); break;
            case 5: cols_fp50_row(c, a, bt, kp, k, 5, 1, mod); break;
            case 6: cols_fp50_row(c, a, bt, kp, k, 6, 1, mod); break;
            case 7: cols_fp50_row(c, a, bt, kp, k, 7, 1, mod); break;
            case 8: cols_fp50_row(c, a, bt, kp, k, 8, 1, mod); break;
            default: FLINT_UNREACHABLE;
        }
    }
}

#else
#define COLS_HAVE_FP50 0
#endif

#else
#define COLS_HAVE_FP50 0
#endif  /* FLINT_BITS == 64 && !NMOD_MAT_HAVE_MUL_U52 */

int
_nmod_mat_mul_cols_simd(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    const slong k = A->c, r = B->c;

    FLINT_ASSERT(1 <= r && r <= NMOD_MAT_MUL_COLS_MAX && k >= 1);

#if FLINT_BITS == 64 && (NMOD_MAT_HAVE_MUL_U52 || COLS_HAVE_FP50)
    if (NMOD_MAT_MUL_COLS_IS_SIMD(A->mod.n))
    {
        /* the columns of B, contiguous and zero-padded to a multiple of 8 */
        const slong kp = (k + 7) & ~(slong) 7;
        ulong * bt = flint_malloc(r * kp * sizeof(ulong));
        slong i, j;

        for (j = 0; j < r; j++)
        {
            for (i = 0; i < k; i++)
                bt[j * kp + i] = nmod_mat_entry(B, i, j);
            for ( ; i < kp; i++)
                bt[j * kp + i] = 0;
        }

#if NMOD_MAT_HAVE_MUL_U52
        /* u52 needs the unreduced products to fit two limbs */
        if (A->mod.n <= (UWORD(1) << 52)
                && FLINT_BIT_COUNT(k) + 2 * FLINT_BIT_COUNT(A->mod.n - 1) <= 128)
            cols_u52(C, A, bt, kp, (int) r);
        else
            cols_u64(C, A, bt, kp, (int) r);
#else
        cols_fp50(C, A, bt, kp, (int) r);
#endif

        flint_free(bt);
        return 1;
    }
#endif

    (void) C; (void) A; (void) B; (void) k; (void) r;
    return 0;
}
