/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5 and Opus 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Matrix multiplication over Z/nZ for small n with the entries packed
    one per byte.

    The underscore method takes strided uint8 arrays, so it can be called
    directly by gr_mat over nmod8 (whose entries are already uint8 with a
    row stride); nmod_mat and gr_mat over nmod reach it through a
    ulong <-> uint8 roundtrip, which is O(mn + mk + kn) against O(mkn) of
    work.

    Kernels. With every entry canonical and less than 16, a byte times a
    scalar is one in-lane table lookup: prod[i] = lut_c[b[i]] where
    lut_c[j] = c*j mod n. One shuffle plus one add per vector replaces the
    multiply entirely, and because the table is already reduced, an
    accumulator byte grows by at most n-1 per step rather than (n-1)^2,
    so the reduction cadence is 255/(n-1) rather than 255/(n-1)^2. This
    is what makes small moduli much faster here than in a floating-point
    gemm: at n = 7 the cadence is 41 steps.

    Without a byte shuffle instruction the portable kernel accumulates in
    16-bit lanes with raw products, giving cadence 65535/(n-1)^2, which
    is still 291 steps at n = 16.

    Above n = 16 the entries no longer index a 16-byte table and the
    advantage disappears (a double-shuffle variant was measured at 0.65x
    of a float gemm at n = 31 and worse beyond), so this method reports
    that it does not apply and callers fall back.

    Strassen is applied on top with ceil splits: quadrant views carry
    their true dimensions, temporaries are allocated at nominal sizes,
    additions zero-extend short operands, and the recursion normalizes
    the operand dimensions to those of the destination, so no peeling or
    padding is required anywhere.
*/

#include <string.h>
#include "nmod_mat.h"
#include "nmod.h"
#include "machine_vectors.h"
#include "thread_pool.h"
#include "thread_support.h"

#if defined(__AVX2__) && FLINT_BITS == 64
# define NMOD_MAT_U8_SHUFFLE 1
# define NMOD_MAT_U8_X86 1
#elif (defined(__ARM_NEON) || defined(_M_ARM64)) && FLINT_BITS == 64
/* vqtbl1q_u8 is the 16-byte table lookup, the same primitive vpshufb
   provides; entries below 16 index it directly */
# define NMOD_MAT_U8_SHUFFLE 1
# define NMOD_MAT_U8_NEON 1
#endif

/*
    AVX512VBMI's vpermi2b looks up a 128-entry byte table across two
    registers in a single instruction, so with it the multiply costs one
    instruction for any modulus up to 128, rather than the two shuffles
    plus add plus canonicalisation that synthesising a 128-entry lookup
    from 16-byte shuffles requires. Note that the useful range is
    narrower than the legal one: the accumulator grows by n-1 per step,
    so the reduction cadence is 255/(n-1), which is 41 at n = 7 but only
    2 at n = 128 -- past roughly n = 32 the reductions dominate and a
    floating-point gemm is the better tool.
*/
#if defined(__AVX512VBMI__) && defined(__AVX512BW__) && FLINT_BITS == 64
# define NMOD_MAT_U8_VBMI 1
#endif

/*
    Largest modulus a byte kernel is worth using for. The lookup width
    permits more (16 with vpshufb, 128 with vpermi2b), but the
    accumulator gains n-1 per step and so must be reduced every
    255/(n-1) - 1 steps: 41 at n = 7, 20 at n = 13, 7 at n = 31, 1 at
    n = 127, and measurement puts the crossover against a
    floating-point multiplication between 13 and 17.
*/
#define NMOD_MAT_U8_MAX_SHUFFLE_MOD 16

/*
    Above that, the leaves are computed in single precision instead:
    convert the block to float, multiply, reduce back to bytes. The
    Strassen layer above stays byte-wide, which is the point -- its
    additions are byte operations and its temporaries are a quarter the
    size of float ones.
*/

/* microkernel tile */
#ifndef U8_MR
# define U8_MR 5
#endif
#ifndef U8_NR
# define U8_NR 64
#endif

/* blocking; KC is bounded by the reduction cadence at run time */
#define U8_KC 256
#define U8_MC 64
#define U8_NC 1024

/* Strassen below this many rows/cols is not worth the temporaries */
#ifndef U8_STRASSEN_CUTOFF
# define U8_STRASSEN_CUTOFF 128
#endif

/*
    Smallest tile the parallel driver will create. This defaults to a
    size that still recurses, but it is a separate knob from the
    Strassen cutoff on purpose: one is about having enough work per
    thread, the other about when recursion pays. Tying them together
    means that raising the cutoff to disable Strassen also disables
    threading, which makes the two effects impossible to measure apart.
*/
/*
    Separate Strassen cutoff for the single-precision leaf. The leaf
    converts its blocks from bytes to float, which is O(leaf^2) against
    O(leaf^3) of multiplication, so it is amortized only when leaves are
    large -- with the byte path's cutoff of 128 the conversion dominates
    and the whole thing loses to a plain sgemm. Large leaves are also
    what makes this worth doing at all: the point is to keep the
    recursion's operands and temporaries in bytes, a quarter the size of
    float ones, which is what allows huge nmod8 matrices to be
    multiplied at all.
*/
#ifndef U8_PACKED_CUTOFF
# define U8_PACKED_CUTOFF 2048
#endif

#ifndef U8_SGEMM_CUTOFF
# define U8_SGEMM_CUTOFF 1024
#endif

#ifndef U8_PAR_MIN_TILE
# define U8_PAR_MIN_TILE (2 * U8_STRASSEN_CUTOFF)
#endif

/*
    Below these row counts the packed path loses to the byte shuffle
    basecase: the M4R table builds are a fixed cost per k-group
    (2^GF2_T resp. GF3_NT * 3^GF3_T entries) that the absorb amortizes
    over the rows, so few-row problems -- including degenerate shapes
    like k = 1 or small m with large k, n -- pay for tables they barely
    read. Crossovers measured on an 8-core Zen 3; the shuffle basecase
    handles every modulus up to 15, so 2 and 3 are always in range.
*/
#ifndef NMOD_MAT_U8_GF2_BASECASE_M
# define NMOD_MAT_U8_GF2_BASECASE_M 64
#endif
#ifndef NMOD_MAT_U8_GF3_BASECASE_M
# define NMOD_MAT_U8_GF3_BASECASE_M 192
#endif

/*
    Minimum tile dimension for the packed (modulus 2 and 3) path. Row
    splits build the M4R tables once per row strip and column splits
    narrow the table rows the absorb streams, so tiles below ~1024 pay
    for their parallelism twice; above it both costs are a few percent.
*/
#ifndef U8_PACKED_PAR_MIN_TILE
# define U8_PACKED_PAR_MIN_TILE 1024
#endif

/*
    Cached tables for one modulus: mul[c][j] = c*j mod n for c, j < 16,
    and red[v] = v mod n for v < 256 as two nibble tables. Building them
    is 256 divisions, which is significant for a small multiplication, so
    the last modulus is cached per thread.
*/
typedef struct
{
    ulong n;
#if defined(NMOD_MAT_U8_VBMI)
    /* 128-entry rows for vpermi2b; 16 KiB, built once per modulus */
    uint8_t mul128[255][128];
#endif
    uint8_t mul[16][16];
    uint8_t red_lo[16];
    uint8_t red_hi[16];
}
nmod_mat_u8_tables;

static FLINT_TLS_PREFIX nmod_mat_u8_tables * u8_tables = NULL;

static void
u8_tables_cleanup(void)
{
    flint_free(u8_tables);
    u8_tables = NULL;
}

static const nmod_mat_u8_tables *
u8_get_tables(ulong n)
{
    nmod_mat_u8_tables * t = u8_tables;
    ulong c, j;

    if (t != NULL && t->n == n)
        return t;

    if (t == NULL)
    {
        t = flint_malloc(sizeof(nmod_mat_u8_tables));
        u8_tables = t;
        flint_register_cleanup_function(u8_tables_cleanup);
    }

    t->n = n;

    /*
        Every entry walks in arithmetic progressions mod n, so the
        tables build with a running value and a conditional subtract
        in place of two remainders per entry; it only runs once per
        modulus, but the VBMI table is 32K entries.
    */
    {
        ulong cr = 0;                     /* c mod n */

        for (c = 0; c < 16; c++)
        {
            ulong v = 0;                  /* cr * j mod n */

            for (j = 0; j < 16; j++)
            {
                t->mul[c][j] = (uint8_t) v;
                v += cr;
                if (v >= n)
                    v -= n;
            }

            cr++;
            if (cr == n)
                cr = 0;
        }
    }

#if defined(NMOD_MAT_U8_VBMI)
    if (n > 16)
        for (c = 0; c < n; c++)
        {
            ulong jr = 0, v = 0;          /* j mod n, c * (j mod n) mod n */

            for (j = 0; j < 128; j++)
            {
                t->mul128[c][j] = (uint8_t) v;
                jr++;
                if (jr == n)
                    jr = 0, v = 0;
                else
                {
                    v += c;
                    if (v >= n)
                        v -= n;
                }
            }
        }
#endif

    {
        ulong lo = 0, hi = 0, h16 = 16 % n;

        for (j = 0; j < 16; j++)
        {
            t->red_lo[j] = (uint8_t) lo;
            t->red_hi[j] = (uint8_t) hi;
            lo++;
            if (lo == n)
                lo = 0;
            hi += h16;
            if (hi >= n)
                hi -= n;
        }
    }

    return t;
}

/* ------------------------------------------------------------------ */
/* basecase                                                           */
/* ------------------------------------------------------------------ */

/*
    C = A * B on packed panels. ap holds U8_MR rows of A as bytes, bp a
    kc x U8_NR panel of B, and the U8_MR x U8_NR accumulator tile is
    reduced every cad steps.
*/
static void
u8_kernel(uint8_t * C, slong Cstride, const uint8_t * ap, slong apstride,
          const uint8_t * bp, slong kc, slong cad, slong rows, slong cols,
          const nmod_mat_u8_tables * tab, ulong n, int first)
{
#if defined(NMOD_MAT_U8_VBMI)
    if (n > 16)
    {
        __m512i acc[U8_MR][U8_NR / 64];
        __m512i vn = _mm512_set1_epi8((char) n);
        __m512i mask = _mm512_set1_epi8(0x0f);
        __m512i tlo = _mm512_broadcast_i32x4(
            _mm_loadu_si128((const __m128i *) tab->red_lo));
        __m512i thi = _mm512_broadcast_i32x4(
            _mm_loadu_si128((const __m128i *) tab->red_hi));
        slong l0, r;
        uint8_t stage[U8_MR][U8_NR];

        for (r = 0; r < U8_MR; r++)
        {
            int h;

            for (h = 0; h < U8_NR / 64; h++)
            {
                if (first || r >= rows)
                    acc[r][h] = _mm512_setzero_si512();
                else if (cols == U8_NR)
                    acc[r][h] = _mm512_loadu_si512(C + r * Cstride + 64 * h);
                else
                {
                    memcpy(stage[r], C + r * Cstride, cols);
                    memset(stage[r] + cols, 0, U8_NR - cols);
                    acc[r][h] = _mm512_loadu_si512(stage[r] + 64 * h);
                }
            }
        }

        for (l0 = 0; l0 < kc; l0 += cad)
        {
            slong le = FLINT_MIN(kc, l0 + cad), l;

            for (l = l0; l < le; l++)
            {
                int h;
                __m512i b[U8_NR / 64];

                for (h = 0; h < U8_NR / 64; h++)
                    b[h] = _mm512_loadu_si512(bp + l * U8_NR + 64 * h);

                for (r = 0; r < U8_MR; r++)
                {
                    const uint8_t * row = tab->mul128[ap[r * apstride + l]];
                    __m512i t0 = _mm512_loadu_si512(row);
                    __m512i t1 = _mm512_loadu_si512(row + 64);

                    for (h = 0; h < U8_NR / 64; h++)
                        acc[r][h] = _mm512_add_epi8(acc[r][h],
                            _mm512_permutex2var_epi8(t0, b[h], t1));
                }
            }

            for (r = 0; r < U8_MR; r++)
            {
                int h;

                for (h = 0; h < U8_NR / 64; h++)
                {
                    __m512i v = acc[r][h];
                    __m512i lo = _mm512_and_si512(v, mask);
                    __m512i hi = _mm512_and_si512(_mm512_srli_epi16(v, 4),
                                                  mask);
                    __m512i sm = _mm512_add_epi8(
                        _mm512_shuffle_epi8(tlo, lo),
                        _mm512_shuffle_epi8(thi, hi));

                    acc[r][h] = _mm512_min_epu8(sm,
                        _mm512_sub_epi8(sm, vn));
                }
            }
        }

        for (r = 0; r < rows; r++)
        {
            int h;

            for (h = 0; h < U8_NR / 64; h++)
                _mm512_storeu_si512(stage[r] + 64 * h, acc[r][h]);

            memcpy(C + r * Cstride, stage[r], cols);
        }

        return;
    }
#endif

#if defined(NMOD_MAT_U8_NEON)
    {
        uint8x16_t acc[U8_MR][U8_NR / 16];
        uint8x16_t vn = vdupq_n_u8((uint8_t) n);
        uint8x16_t mask = vdupq_n_u8(0x0f);
        uint8x16_t tlo = vld1q_u8(tab->red_lo);
        uint8x16_t thi = vld1q_u8(tab->red_hi);
        slong l0, r;
        uint8_t stage[U8_MR][U8_NR];

        for (r = 0; r < U8_MR; r++)
        {
            int h;

            for (h = 0; h < U8_NR / 16; h++)
            {
                if (first || r >= rows)
                    acc[r][h] = vdupq_n_u8(0);
                else if (cols == U8_NR)
                    acc[r][h] = vld1q_u8(C + r * Cstride + 16 * h);
                else
                {
                    memcpy(stage[r], C + r * Cstride, cols);
                    memset(stage[r] + cols, 0, U8_NR - cols);
                    acc[r][h] = vld1q_u8(stage[r] + 16 * h);
                }
            }
        }

        for (l0 = 0; l0 < kc; l0 += cad)
        {
            slong le = FLINT_MIN(kc, l0 + cad), l;

            for (l = l0; l < le; l++)
            {
                int h;
                uint8x16_t b[U8_NR / 16];

                for (h = 0; h < U8_NR / 16; h++)
                    b[h] = vld1q_u8(bp + l * U8_NR + 16 * h);

                for (r = 0; r < U8_MR; r++)
                {
                    uint8x16_t t = vld1q_u8(tab->mul[ap[r * apstride + l]]);

                    for (h = 0; h < U8_NR / 16; h++)
                        acc[r][h] = vaddq_u8(acc[r][h], vqtbl1q_u8(t, b[h]));
                }
            }

            for (r = 0; r < U8_MR; r++)
            {
                int h;

                for (h = 0; h < U8_NR / 16; h++)
                {
                    uint8x16_t v = acc[r][h];
                    uint8x16_t lo = vandq_u8(v, mask);
                    uint8x16_t hi = vshrq_n_u8(v, 4);
                    uint8x16_t sm = vaddq_u8(vqtbl1q_u8(tlo, lo),
                                             vqtbl1q_u8(thi, hi));

                    acc[r][h] = vminq_u8(sm, vsubq_u8(sm, vn));
                }
            }
        }

        for (r = 0; r < rows; r++)
        {
            int h;

            for (h = 0; h < U8_NR / 16; h++)
                vst1q_u8(stage[r] + 16 * h, acc[r][h]);

            memcpy(C + r * Cstride, stage[r], cols);
        }

        return;
    }
#elif defined(NMOD_MAT_U8_X86)
    __m256i acc[U8_MR][U8_NR / 32];
    __m256i vn = _mm256_set1_epi8((char) n);
    __m256i mask = _mm256_set1_epi8(0x0f);
    __m256i tlo = _mm256_broadcastsi128_si256(
        _mm_loadu_si128((const __m128i *) tab->red_lo));
    __m256i thi = _mm256_broadcastsi128_si256(
        _mm_loadu_si128((const __m128i *) tab->red_hi));
    slong l0, r;
    uint8_t stage[U8_MR][U8_NR];

    /*
        Continuing a k block: start from what C already holds, which is
        canonical and so counts as one step's worth of growth. A partial
        strip or short row block must be read through the padded staging
        tile, since a vector load would run past the end of C.
    */
    for (r = 0; r < U8_MR; r++)
    {
        int h;
        const uint8_t * src = NULL;

        if (!first && r < rows)
        {
            if (cols == U8_NR)
                src = C + r * Cstride;
            else
            {
                memcpy(stage[r], C + r * Cstride, cols);
                memset(stage[r] + cols, 0, U8_NR - cols);
                src = stage[r];
            }
        }

        for (h = 0; h < U8_NR / 32; h++)
            acc[r][h] = (src == NULL) ? _mm256_setzero_si256()
                : _mm256_loadu_si256((const __m256i *) (src + 32 * h));
    }

    for (l0 = 0; l0 < kc; l0 += cad)
    {
        slong le = FLINT_MIN(kc, l0 + cad), l;

        for (l = l0; l < le; l++)
        {
            __m256i b[U8_NR / 32];
            int h;

            for (h = 0; h < U8_NR / 32; h++)
                b[h] = _mm256_load_si256((const __m256i *)
                                         (bp + l * U8_NR + 32 * h));

            /* several 32-byte strips per table load: the shuffle is the
               work and the broadcast is overhead, so the strip is as
               wide as the register budget allows */
            for (r = 0; r < U8_MR; r++)
            {
                __m256i t = _mm256_broadcastsi128_si256(
                    _mm_loadu_si128((const __m128i *)
                        tab->mul[ap[r * apstride + l]]));

                for (h = 0; h < U8_NR / 32; h++)
                    acc[r][h] = _mm256_add_epi8(acc[r][h],
                                                _mm256_shuffle_epi8(t, b[h]));
            }
        }

        for (r = 0; r < U8_MR; r++)
        {
            int h;

            for (h = 0; h < U8_NR / 32; h++)
            {
                __m256i v = acc[r][h];
                __m256i lo = _mm256_and_si256(v, mask);
                __m256i hi = _mm256_and_si256(_mm256_srli_epi16(v, 4), mask);
                __m256i s = _mm256_add_epi8(_mm256_shuffle_epi8(tlo, lo),
                                            _mm256_shuffle_epi8(thi, hi));
                acc[r][h] = _mm256_min_epu8(s, _mm256_sub_epi8(s, vn));
            }
        }
    }

    if (cols == U8_NR)
    {
        /* full strip: store straight into C */
        for (r = 0; r < rows; r++)
        {
            int h;

            for (h = 0; h < U8_NR / 32; h++)
                _mm256_storeu_si256((__m256i *) (C + r * Cstride + 32 * h),
                                    acc[r][h]);
        }
    }
    else
    {
        for (r = 0; r < rows; r++)
        {
            int h;

            for (h = 0; h < U8_NR / 32; h++)
                _mm256_storeu_si256((__m256i *) (stage[r] + 32 * h), acc[r][h]);

            memcpy(C + r * Cstride, stage[r], cols);
        }
    }
#else
    uint16_t acc[U8_MR][U8_NR];
    slong l0, r, j;

    memset(acc, 0, sizeof(acc));

    if (!first)
        for (r = 0; r < rows; r++)
            for (j = 0; j < cols; j++)
                acc[r][j] = C[r * Cstride + j];

    for (l0 = 0; l0 < kc; l0 += cad)
    {
        slong le = FLINT_MIN(kc, l0 + cad), l;

        for (l = l0; l < le; l++)
            for (r = 0; r < U8_MR; r++)
            {
                uint16_t c = ap[r * apstride + l];

                for (j = 0; j < U8_NR; j++)
                    acc[r][j] += c * (uint16_t) bp[l * U8_NR + j];
            }

        for (r = 0; r < U8_MR; r++)
            for (j = 0; j < U8_NR; j++)
                acc[r][j] %= (uint16_t) n;
    }

    for (r = 0; r < rows; r++)
        for (j = 0; j < cols; j++)
            C[r * Cstride + j] = (uint8_t) acc[r][j];
#endif
}

/*
    Column strips outermost: a strip of B is packed once for the whole
    inner dimension and stays resident while every row block of A sweeps
    it, so the accumulators live in registers across the entire product
    and each entry of C is written exactly once. A needs no packing at
    all, its rows already being contiguous bytes.
*/
static void
u8_mul_basecase(uint8_t * C, slong Cstride,
                const uint8_t * A, slong Astride,
                const uint8_t * B, slong Bstride,
                slong m, slong k, slong n, ulong nn)
{
    const nmod_mat_u8_tables * tab = u8_get_tables(nn);
    slong jc, i, l, cad, tail;
    uint8_t * bp;
    uint8_t * apad;
    TMP_INIT;

    /*
        A chunk begins with whatever the previous reduction left, which
        is below nn, and then adds cad products. The byte kernels have
        products below nn (the tables are reduced), the 16-bit one has
        products below (nn-1)^2, so the bound is
            (nn-1) + cad * bound(product) <= limit
        and dropping the leading term -- as an earlier version did --
        overflows for near-maximal entries. That went unnoticed while
        the modulus was capped at 16, where it needs almost the largest
        possible entries to show up; at nn = 92 it fails at once.
    */
#if defined(NMOD_MAT_U8_SHUFFLE)
    cad = (nn <= 1) ? 1 : (255 - (slong) (nn - 1)) / (slong) (nn - 1);
#else
    cad = (nn <= 1) ? 1 : (65535 - (slong) (nn - 1))
                          / (slong) ((nn - 1) * (nn - 1));
#endif
    cad = FLINT_MAX(cad, 1);

    TMP_START;
    /*
        One packed strip of B, plus a padded copy of the final short row
        block of A. The kernel always reads U8_MR rows, so a block with
        fewer must be given full-height storage rather than a pointer
        into A, which would read past the matrix.
    */
    tail = m % U8_MR;
    bp = (uint8_t *) TMP_ALLOC(k * U8_NR + (tail ? U8_MR * k : 0) + 64);
    bp = (uint8_t *) ((((uintptr_t) bp) + 31) & ~(uintptr_t) 31);
    apad = tail ? bp + k * U8_NR : NULL;

    if (tail)
    {
        for (i = 0; i < U8_MR; i++)
        {
            if (i < tail)
                memcpy(apad + i * k, A + (m - tail + i) * Astride, k);
            else
                memset(apad + i * k, 0, k);
        }
    }

    for (jc = 0; jc < n; jc += U8_NR)
    {
        slong cols = FLINT_MIN(U8_NR, n - jc);

        for (l = 0; l < k; l++)
        {
            memcpy(bp + l * U8_NR, B + l * Bstride + jc, cols);
            memset(bp + l * U8_NR + cols, 0, U8_NR - cols);
        }

        for (i = 0; i + U8_MR <= m; i += U8_MR)
            u8_kernel(C + i * Cstride + jc, Cstride,
                      A + i * Astride, Astride, bp, k, cad,
                      U8_MR, cols, tab, nn, 1);

        if (tail)
            u8_kernel(C + (m - tail) * Cstride + jc, Cstride,
                      apad, k, bp, k, cad, tail, cols, tab, nn, 1);
    }

    TMP_END;
}



/* byte -> float, the widening the compiler tends not to vectorize */
static void
u8_to_float(float * dst, const uint8_t * src, slong len)
{
    slong i = 0;

#if defined(NMOD_MAT_U8_X86)
    for ( ; i + 8 <= len; i += 8)
        _mm256_storeu_ps(dst + i,
            _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(
                _mm_loadl_epi64((const __m128i *) (src + i)))));
#elif defined(NMOD_MAT_U8_NEON)
    for ( ; i + 8 <= len; i += 8)
    {
        uint16x8_t w = vmovl_u8(vld1_u8(src + i));

        vst1q_f32(dst + i, vcvtq_f32_u32(vmovl_u16(vget_low_u16(w))));
        vst1q_f32(dst + i + 4, vcvtq_f32_u32(vmovl_u16(vget_high_u16(w))));
    }
#endif

    for ( ; i < len; i++)
        dst[i] = (float) src[i];
}

/*
    Single precision basecase, for moduli too large for a byte
    accumulator. A float32 mantissa holds every integer below 2^24
    exactly, so the inner dimension is chunked at
        kmax = (2^24 - 1) / (n-1)^2
    which is 466033 at n = 7 but still 260 at n = 255, i.e. never a
    real constraint. Only the leaf is float; A, B and C stay bytes, so
    the recursion above is unaffected.
*/
/*
    Reduce a row of nonnegative integer values below 2^31, held as
    floats (sgemm output, exact below 2^24) or as uint32 accumulations
    of such blocks, to bytes mod nn. The vector path runs four doubles
    at a time through machine_vectors' reduce_to_0n -- one fnmadd of a
    rounded quotient and a blend fixup, exact in this range -- and
    packs sixteen results to bytes per iteration; per entry that
    replaces a scalar division whose latency dominated the multi-block
    range. The uint32 inputs stay below 2^31 (the accumulator folds
    earlier than that), so the signed epi32 -> double conversion is
    safe.
*/
static void
u8_reduce_row_f32(uint8_t * dst, const float * src, slong n, nmod_t mod)
{
    slong j = 0;

#if defined(NMOD_MAT_U8_X86)
    {
        vec4d vn = vec4d_set_d((double) mod.n);
        vec4d vninv = vec4d_set_d(1.0 / (double) mod.n);

        for (; j + 16 <= n; j += 16)
        {
            __m128i w0, w1;
            slong t;
            __m128i r[4];

            for (t = 0; t < 4; t++)
            {
                vec4d x = _mm256_cvtps_pd(
                              _mm_loadu_ps(src + j + 4 * t));
                x = vec4d_reduce_to_0n(x, vn, vninv);
                r[t] = _mm256_cvttpd_epi32(x);
            }

            w0 = _mm_packus_epi32(r[0], r[1]);
            w1 = _mm_packus_epi32(r[2], r[3]);
            _mm_storeu_si128((__m128i *) (dst + j),
                             _mm_packus_epi16(w0, w1));
        }
    }
#endif

    for (; j < n; j++)
    {
        ulong r;
        NMOD_RED(r, (ulong) src[j], mod);
        dst[j] = (uint8_t) r;
    }
}

static void
u8_reduce_row_u32(uint8_t * dst, const uint32_t * src, slong n, nmod_t mod)
{
    slong j = 0;

#if defined(NMOD_MAT_U8_X86)
    {
        vec4d vn = vec4d_set_d((double) mod.n);
        vec4d vninv = vec4d_set_d(1.0 / (double) mod.n);

        for (; j + 16 <= n; j += 16)
        {
            __m128i w0, w1;
            slong t;
            __m128i r[4];

            for (t = 0; t < 4; t++)
            {
                vec4d x = _mm256_cvtepi32_pd(
                              _mm_loadu_si128((const __m128i *)
                                              (src + j + 4 * t)));
                x = vec4d_reduce_to_0n(x, vn, vninv);
                r[t] = _mm256_cvttpd_epi32(x);
            }

            w0 = _mm_packus_epi32(r[0], r[1]);
            w1 = _mm_packus_epi32(r[2], r[3]);
            _mm_storeu_si128((__m128i *) (dst + j),
                             _mm_packus_epi16(w0, w1));
        }
    }
#endif

    for (; j < n; j++)
    {
        ulong r;
        NMOD_RED(r, (ulong) src[j], mod);
        dst[j] = (uint8_t) r;
    }
}

static void
u8_mul_basecase_sgemm(uint8_t * C, slong Cstride,
                      const uint8_t * A, slong Astride,
                      const uint8_t * B, slong Bstride,
                      slong m, slong k, slong n, ulong nn)
{
    slong kmax, pc, i, j;
    float * af;
    float * bf;
    float * cf;
    uint32_t * cacc;
    ulong bound;
    nmod_t mod;
    int same;
    TMP_INIT;

    nmod_init(&mod, nn);

    kmax = (slong) ((16777215.0) / ((double) (nn - 1) * (double) (nn - 1)));
    kmax = FLINT_MAX(kmax, 1);
    kmax = FLINT_MIN(kmax, k);

    /*
        Squaring: convert the matrix to floats once and hand sgemm
        strided views of it -- the A block is columns pc..pc+kc of all
        rows and the B block is rows pc..pc+kc, so with the whole
        matrix resident both are (F + pc, lda = k) and
        (F + pc * k, ldb = k). Halves the conversion work; the packed
        per-block slices are only needed when the operands differ.
    */
    same = (A == B && Astride == Bstride && m == k && k == n);

    TMP_START;
    if (same)
    {
        af = (float *) TMP_ALLOC((m * k + m * n) * sizeof(float)
                             + ((kmax < k) ? m * n * sizeof(uint32_t) : 0));
        bf = af;
        cf = af + m * k;

        for (i = 0; i < m; i++)
            u8_to_float(af + i * k, A + i * Astride, k);
    }
    else
    {
        af = (float *) TMP_ALLOC((m * kmax + kmax * n + m * n)
                                 * sizeof(float)
                             + ((kmax < k) ? m * n * sizeof(uint32_t) : 0));
        bf = af + m * kmax;
        cf = bf + kmax * n;
    }
    cacc = (uint32_t *) (cf + m * n);

    /*
        Partial sums per sgemm block stay below 2^24 by the choice of
        kmax; blocks accumulate exactly in a uint32 image of C and the
        modular reduction runs once at the end. The earlier version
        reduced C after every block, a scalar division per entry per
        block that dominated the multi-block range (measured at about
        half the total at modulus 251, k = 1000) and let
        nmod_mat_mul_blas overtake what is otherwise the same gemm.
        The accumulator is folded in the rare case the running bound
        would overflow 32 bits.
    */
    bound = 0;

    for (pc = 0; pc < k; pc += kmax)
    {
        slong kc = FLINT_MIN(kmax, k - pc);

        if (same)
        {
            flint_sgemm(m, kc, n, af + pc, k, bf + pc * k, k, cf, n);
        }
        else
        {
            for (i = 0; i < m; i++)
                u8_to_float(af + i * kc, A + i * Astride + pc, kc);

            for (i = 0; i < kc; i++)
                u8_to_float(bf + i * n, B + (pc + i) * Bstride, n);

            flint_sgemm(m, kc, n, af, kc, bf, n, cf, n);
        }

        if (kmax >= k)
        {
            /* single block: reduce straight into C */
            for (i = 0; i < m; i++)
                u8_reduce_row_f32(C + i * Cstride, cf + i * n, n, mod);
            break;
        }

        if (bound + UWORD(16777216) > (UWORD(1) << 31))
        {
            for (i = 0; i < m * n; i++)
            {
                ulong r;
                NMOD_RED(r, (ulong) cacc[i], mod);
                cacc[i] = (uint32_t) r;
            }
            bound = nn - 1;
        }

        if (pc == 0)
        {
            for (i = 0; i < m; i++)
                for (j = 0; j < n; j++)
                    cacc[i * n + j] = (uint32_t) cf[i * n + j];
        }
        else
        {
            for (i = 0; i < m; i++)
                for (j = 0; j < n; j++)
                    cacc[i * n + j] += (uint32_t) cf[i * n + j];
        }

        bound += UWORD(16777216);
    }

    if (kmax < k)
    {
        for (i = 0; i < m; i++)
            u8_reduce_row_u32(C + i * Cstride, cacc + i * n, n, mod);
    }

    TMP_END;
}

/* ------------------------------------------------------------------ */
/* Flat threaded gemm for moduli above the shuffle range               */
/* ------------------------------------------------------------------ */

/*
    The tiled path parallelizes moduli above 15 by running the
    Strassen + sgemm-leaf recursion per tile of C, which re-converts
    slices per tile and scaled to ~2.3x on 8 Zen 3 cores where
    nmod_mat_mul_blas -- one shared lift, one pool-threaded gemm --
    reached 4.6x. This driver mirrors the blas shape in the byte
    domain (bytes are an eighth of the ulong traffic): convert both
    operands to a float image once, in parallel over rows, run one
    pool-threaded flint_sgemm, reduce in parallel. Squaring shares the
    image.

    It only engages when the whole product fits a single exact sgemm,
    (n - 1)^2 k < 2^24: a multi-block variant that accumulated exactly
    in uint32 was measured neutral at 8 threads and 30% slower at 4 on
    the same Zen 3 -- at modulus 251 the exactness bound forces either
    the accumulate traffic or double precision at half the rate, the
    same wall mul_blas hits, so the shared lift buys nothing there and
    the tiled path keeps that range. (The structural fix for large
    moduli would be an integer gemm on VNNI-class hardware, where
    vpdpbusd accumulates byte products in uint32 natively; future
    work.) It also declines with one thread, where the tiled path's
    Strassen saves real work, and below 256 x 256, where the
    conversions do not amortize.
*/

typedef struct
{
    const uint8_t * A; slong Astride;
    const uint8_t * B; slong Bstride;
    uint8_t * C; slong Cstride;
    slong m, k, n;
    float * af;
    float * bf;
    float * cf;
    nmod_t mod;
    int share;
}
u8_flat_ctx;

typedef struct
{
    u8_flat_ctx * ctx;
    int phase;            /* 0 convert, 1 reduce */
    slong start, end;     /* row ranges; phase 0 counts A then B rows */
}
u8_flat_job;

static void
u8_flat_worker(void * varg)
{
    u8_flat_job * job = (u8_flat_job *) varg;
    u8_flat_ctx * ctx = job->ctx;
    slong r;

    if (job->phase == 0)
    {
        for (r = job->start; r < job->end; r++)
        {
            if (r < ctx->m)
                u8_to_float(ctx->af + r * ctx->k,
                            ctx->A + r * ctx->Astride, ctx->k);
            else
                u8_to_float(ctx->bf + (r - ctx->m) * ctx->n,
                            ctx->B + (r - ctx->m) * ctx->Bstride,
                            ctx->n);
        }
    }
    else
    {
        for (r = job->start; r < job->end; r++)
            u8_reduce_row_f32(ctx->C + r * ctx->Cstride,
                              ctx->cf + r * ctx->n, ctx->n, ctx->mod);
    }
}

static void
u8_flat_run_phase(u8_flat_ctx * ctx, int phase, slong total,
                  thread_pool_handle * handles, slong nw,
                  u8_flat_job * jobs)
{
    slong tid, nt = nw + 1;

    for (tid = 0; tid < nt; tid++)
    {
        jobs[tid].ctx = ctx;
        jobs[tid].phase = phase;
        jobs[tid].start = (total * tid) / nt;
        jobs[tid].end = (total * (tid + 1)) / nt;
    }

    for (tid = 0; tid < nw; tid++)
        thread_pool_wake(global_thread_pool, handles[tid], 0,
                         u8_flat_worker, jobs + tid);

    u8_flat_worker(jobs + nw);

    for (tid = 0; tid < nw; tid++)
        thread_pool_wait(global_thread_pool, handles[tid]);
}

static int
u8_mul_flat_sgemm(uint8_t * C, slong Cstride,
                  const uint8_t * A, slong Astride,
                  const uint8_t * B, slong Bstride,
                  slong m, slong k, slong n, ulong nn)
{
    slong nthreads = flint_get_num_threads();
    slong nw, max_workers;
    u8_flat_ctx ctx;
    u8_flat_job * jobs;
    thread_pool_handle * handles;
    TMP_INIT;

    if (nthreads <= 1 || m < 256 || n < 256
        || (double) (nn - 1) * (double) (nn - 1) * (double) k
           >= 16777216.0)
        return 0;

    ctx.A = A; ctx.Astride = Astride;
    ctx.B = B; ctx.Bstride = Bstride;
    ctx.C = C; ctx.Cstride = Cstride;
    ctx.m = m; ctx.k = k; ctx.n = n;
    nmod_init(&ctx.mod, nn);
    ctx.share = (A == B && Astride == Bstride && m == k && k == n);

    TMP_START;
    ctx.af = (float *) TMP_ALLOC(
        (m * k + (ctx.share ? 0 : k * n) + m * n) * sizeof(float));
    ctx.bf = ctx.share ? ctx.af : ctx.af + m * k;
    ctx.cf = ctx.bf + k * n;

    nw = flint_request_threads(&handles, nthreads);
    max_workers = nw;
    jobs = (u8_flat_job *) TMP_ALLOC((nw + 1) * sizeof(u8_flat_job));

    u8_flat_run_phase(&ctx, 0, m + (ctx.share ? 0 : k), handles, nw,
                      jobs);

    /* the fallback gemm threads through the pool itself: hand the
       workers back around the call (see mul_blas.c) */
    flint_give_back_threads(handles, nw);
    flint_sgemm(m, k, n, ctx.af, k, ctx.bf, n, ctx.cf, n);
    nw = flint_request_threads(&handles, max_workers + 1);

    u8_flat_run_phase(&ctx, 1, m, handles, nw, jobs);

    flint_give_back_threads(handles, nw);

    TMP_END;

    return 1;
}

/*
    Runtime copies of the packed Strassen cutoffs (the rationale for
    the defaults sits with the engines below), deliberately not in the
    header: tests and tuning runs read and write them through a local
    "FLINT_DLL extern slong ..." declaration. Defined unconditionally:
    the test references them on every platform, so a 32-bit build,
    where the packed engines are compiled out and the values are
    inert, must still export the symbols.
*/
#ifndef GF2_STRASSEN_CUTOFF
# define GF2_STRASSEN_CUTOFF 4096
#endif
#ifndef GF3_STRASSEN_CUTOFF
# define GF3_STRASSEN_CUTOFF 8192
#endif

FLINT_DLL slong nmod_mat_gf2_strassen_cutoff = GF2_STRASSEN_CUTOFF;
FLINT_DLL slong nmod_mat_gf3_strassen_cutoff = GF3_STRASSEN_CUTOFF;

/*
    The bit-packed engines below store 64 entries per ulong and assume
    64-bit limbs throughout (word counts, shifts, the SWAR conversion
    constants); a 32-bit build routes moduli 2 and 3 through the byte
    shuffle machinery instead, which handles them correctly if without
    the packed speed.
*/
#if FLINT_BITS == 64

/* ------------------------------------------------------------------ */
/* GF(2): one bit per entry, method of four Russians                   */
/* ------------------------------------------------------------------ */

/*
    At modulus 2 the byte representation wastes eight bits per entry and
    the multiplication is an XOR, so a bit-packed kernel is a different
    algorithm rather than a tuning of the same one.

    The method of four Russians processes t columns of A at a time: for
    each group of t rows of B it precomputes all 2^t linear combinations
    of them -- one XOR per table entry, walking the indices in Gray code
    order so consecutive entries differ in one row -- after which a row
    of C absorbs t columns of A with a single table lookup and XOR.
    That replaces t row-XORs by one, and the table cost is amortized
    over all m rows.
*/

#ifndef GF2_T
# define GF2_T 8
#endif

/* Gray-code tables built and consumed together */
#ifndef GF2_NT
# define GF2_NT 8
#endif

/*
    Leaf size for the bit-packed Strassen at modulus 2. Independent of
    U8_PACKED_CUTOFF, which governs the byte-domain recursion above the
    packed GF(3) basecase: the mod-2 leaf improves through 4096 (the
    multi-table absorb pass is C-row-traffic bound, and doubling the
    leaf halves the number of times each row of C is revisited faster
    than the table set outgrows L2).
    (Default 4096; the runtime copy is defined above the packed
    region.)
*/

#if defined(__AVX512BW__) && FLINT_BITS == 64
# define NMOD_MAT_U8_GF2_AVX512 1
#endif

/*
    The same conversions straight from reduced ulong entries, for the
    public wrapper: an nmod_mat holds 8-byte entries, and converting
    them to a byte image only to re-pack the bytes to bits doubles the
    dominant memory traffic of small-modulus products. Entries are
    assumed reduced (below the modulus), as everywhere in this file.
*/
static void
gf2_pack_row_ul(ulong * dst, const ulong * src, slong n, slong w)
{
    slong j = 0, words = (n + 63) / 64;
    ulong acc = 0;

#if defined(NMOD_MAT_U8_GF2_AVX512)
    for (; j + 64 <= n; j += 64)
    {
        ulong word = 0;
        slong t;

        for (t = 0; t < 64; t += 8)
        {
            __m512i v = _mm512_loadu_si512((const void *) (src + j + t));
            word |= (ulong) _mm512_test_epi64_mask(v,
                                _mm512_set1_epi64(1)) << t;
        }

        dst[j / 64] = word;
    }
#elif defined(NMOD_MAT_U8_X86)
    for (; j + 64 <= n; j += 64)
    {
        ulong word = 0;
        slong t;

        for (t = 0; t < 64; t += 4)
        {
            __m256i v = _mm256_loadu_si256((const __m256i *)
                                           (src + j + t));
            v = _mm256_slli_epi64(v, 63);
            word |= (ulong) _mm256_movemask_pd(_mm256_castsi256_pd(v))
                    << t;
        }

        dst[j / 64] = word;
    }
#endif

    acc = 0;
    for (; j < n; j++)
    {
        acc |= (src[j] & 1) << (j & 63);
        if ((j & 63) == 63)
        {
            dst[j / 64] = acc;
            acc = 0;
        }
    }
    if (n & 63)
        dst[n / 64] = acc;

    for (j = words; j < w; j++)
        dst[j] = 0;
}

static void
gf2_unpack_row_ul(ulong * dst, const ulong * src, slong n)
{
    slong j = 0;

#if defined(NMOD_MAT_U8_GF2_AVX512)
    for (; j + 8 <= n; j += 8)
    {
        __mmask8 msk = (__mmask8) (src[j / 64] >> (j & 63));
        _mm512_storeu_si512((void *) (dst + j),
                            _mm512_maskz_set1_epi64(msk, 1));
    }
#endif

    for (; j < n; j++)
        dst[j] = (src[j / 64] >> (j & 63)) & 1;
}

static void
gf3_pack_row_ul(ulong * d1, ulong * d2, const ulong * src, slong n,
                slong w)
{
    slong j = 0, words = (n + 63) / 64;
    ulong a1 = 0, a2 = 0;

#if defined(NMOD_MAT_U8_GF2_AVX512)
    for (; j + 64 <= n; j += 64)
    {
        ulong w1 = 0, w2 = 0;
        slong t;

        for (t = 0; t < 64; t += 8)
        {
            __m512i v = _mm512_loadu_si512((const void *) (src + j + t));
            w1 |= (ulong) _mm512_test_epi64_mask(v,
                              _mm512_set1_epi64(1)) << t;
            w2 |= (ulong) _mm512_test_epi64_mask(v,
                              _mm512_set1_epi64(2)) << t;
        }

        d1[j / 64] = w1;
        d2[j / 64] = w2;
    }
#elif defined(NMOD_MAT_U8_X86)
    for (; j + 64 <= n; j += 64)
    {
        ulong w1 = 0, w2 = 0;
        slong t;

        for (t = 0; t < 64; t += 4)
        {
            __m256i v = _mm256_loadu_si256((const __m256i *)
                                           (src + j + t));
            __m256i b1 = _mm256_slli_epi64(v, 63);
            __m256i b2 = _mm256_slli_epi64(v, 62);

            w1 |= (ulong) _mm256_movemask_pd(_mm256_castsi256_pd(b1))
                  << t;
            w2 |= (ulong) _mm256_movemask_pd(_mm256_castsi256_pd(b2))
                  << t;
        }

        d1[j / 64] = w1;
        d2[j / 64] = w2;
    }
#endif

    a1 = a2 = 0;
    for (; j < n; j++)
    {
        a1 |= (src[j] & 1) << (j & 63);
        a2 |= ((src[j] >> 1) & 1) << (j & 63);
        if ((j & 63) == 63)
        {
            d1[j / 64] = a1;
            d2[j / 64] = a2;
            a1 = a2 = 0;
        }
    }
    if (n & 63)
    {
        d1[n / 64] = a1;
        d2[n / 64] = a2;
    }

    for (j = words; j < w; j++)
        d1[j] = d2[j] = 0;
}

static void
gf3_unpack_row_ul(ulong * dst, const ulong * s1, const ulong * s2,
                  slong n)
{
    slong j = 0;

#if defined(NMOD_MAT_U8_GF2_AVX512)
    for (; j + 8 <= n; j += 8)
    {
        __mmask8 m1 = (__mmask8) (s1[j / 64] >> (j & 63));
        __mmask8 m2 = (__mmask8) (s2[j / 64] >> (j & 63));
        __m512i v = _mm512_or_si512(_mm512_maskz_set1_epi64(m1, 1),
                                    _mm512_maskz_set1_epi64(m2, 2));

        _mm512_storeu_si512((void *) (dst + j), v);
    }
#endif

    for (; j < n; j++)
        dst[j] = ((s1[j / 64] >> (j & 63)) & 1)
               + 2 * ((s2[j / 64] >> (j & 63)) & 1);
}

/*
    Row conversions between canonical mod-2 bytes and packed bits.
    These are the O(n^2) boundary of the bit-packed kernel, and done a
    bit at a time they cost more than the entire multiplication below
    n ~ 4096 -- pack A + pack B + unpack C measured 10x the M4R kernel
    at 1024 and was the whole of a 15x gap against M4RI, whose
    interface is bit-packed and pays no conversion. Byte tests become
    one instruction per 64 entries with AVX-512BW (vptestmb) and one
    per 32 with AVX2 (vpcmpgtb + vpmovmskb); the portable form packs 8
    bytes at a time with the usual multiply gather, which places the
    LSB of byte i at bit 56 + i (the partial products land at
    56 + 8i - 7j, each position hit at most once, so no carries).

    The word tail beyond n bits is zeroed; the padded Strassen layout
    relies on that.
*/
FLINT_FORCE_INLINE void
gf2_pack_row(ulong * dst, const uint8_t * src, slong n, slong words)
{
    slong j = 0, w = 0;

#if defined(NMOD_MAT_U8_GF2_AVX512)
    for (; j + 64 <= n; j += 64)
    {
        __m512i v = _mm512_loadu_si512((const void *) (src + j));
        dst[w++] = (ulong) _mm512_test_epi8_mask(v, v);
    }
#elif defined(NMOD_MAT_U8_X86)
    for (; j + 64 <= n; j += 64)
    {
        __m256i lo = _mm256_loadu_si256((const __m256i *) (src + j));
        __m256i hi = _mm256_loadu_si256((const __m256i *) (src + j + 32));
        __m256i z = _mm256_setzero_si256();
        ulong mlo, mhi;

        mlo = (unsigned int) _mm256_movemask_epi8(_mm256_cmpgt_epi8(lo, z));
        mhi = (unsigned int) _mm256_movemask_epi8(_mm256_cmpgt_epi8(hi, z));
        dst[w++] = mlo | (mhi << 32);
    }
#else
    for (; j + 64 <= n; j += 64)
    {
        ulong bits = 0;
        slong b;

        for (b = 0; b < 64; b += 8)
        {
            ulong x;

            memcpy(&x, src + j + b, 8);
            bits |= ((x * UWORD(0x0102040810204080)) >> 56) << b;
        }

        dst[w++] = bits;
    }
#endif

    if (j < n)
    {
        ulong bits = 0;
        slong b;

        for (b = 0; j < n; j++, b++)
            bits |= (ulong) (src[j] & 1) << b;

        dst[w++] = bits;
    }

    for (; w < words; w++)
        dst[w] = 0;
}

FLINT_FORCE_INLINE void
gf2_unpack_row(uint8_t * dst, const ulong * src, slong n)
{
    slong j = 0;

#if defined(NMOD_MAT_U8_GF2_AVX512)
    for (; j + 64 <= n; j += 64)
        _mm512_storeu_si512((void *) (dst + j),
                _mm512_maskz_set1_epi8((__mmask64) src[j / 64], 1));
#elif defined(NMOD_MAT_U8_X86)
    {
        const __m256i sel = _mm256_setr_epi8(
            0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,
            2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3);
        const __m256i bit = _mm256_setr_epi8(
            1, 2, 4, 8, 16, 32, 64, (char) 128,
            1, 2, 4, 8, 16, 32, 64, (char) 128,
            1, 2, 4, 8, 16, 32, 64, (char) 128,
            1, 2, 4, 8, 16, 32, 64, (char) 128);
        const __m256i one = _mm256_set1_epi8(1);

        for (; j + 64 <= n; j += 64)
        {
            ulong x = src[j / 64];
            __m256i lo = _mm256_set1_epi32((int) (uint32_t) x);
            __m256i hi = _mm256_set1_epi32((int) (uint32_t) (x >> 32));

            lo = _mm256_shuffle_epi8(lo, sel);
            hi = _mm256_shuffle_epi8(hi, sel);
            lo = _mm256_min_epu8(_mm256_and_si256(lo, bit), one);
            hi = _mm256_min_epu8(_mm256_and_si256(hi, bit), one);
            _mm256_storeu_si256((__m256i *) (dst + j), lo);
            _mm256_storeu_si256((__m256i *) (dst + j + 32), hi);
        }
    }
#else
    for (; j + 64 <= n; j += 64)
    {
        ulong x = src[j / 64];
        slong b;

        for (b = 0; b < 64; b += 8)
        {
            /* spread bits to byte LSB positions, then to 0/1 bytes */
            ulong t = ((x >> b) & 0xff) * UWORD(0x0101010101010101);

            t &= UWORD(0x8040201008040201);
            t = ((t + UWORD(0x7f7f7f7f7f7f7f7f)) >> 7)
                & UWORD(0x0101010101010101);
            memcpy(dst + j + b, &t, 8);
        }
    }
#endif

    for (; j < n; j++)
        dst[j] = (uint8_t) ((src[j / 64] >> (j % 64)) & 1);
}

/* ------------------------------------------------------------------ */
/* GF(3): bit-packed Strassen over two planes                          */
/* ------------------------------------------------------------------ */

/*
    Elements are held as two bit-planes, 0 = (0,0), 1 = (1,0),
    2 = (0,1), so one pair of words carries 64 elements. Addition is
    six bitwise operations on the pair and negation is a swap of the
    planes, so subtraction is the same six operations reading the
    subtrahend's planes in the other order -- at the view level it is
    free. The encoding and its adder are the complemented KAT form,
    verified exhaustively over all 9 input pairs in the prototype.

    The multiplication is the same shape as the modulus-2 one: pack
    once, run Strassen in the packed domain where a quadrant addition
    costs six ops per 64 elements, and use a multi-table method of four
    Russians leaf. The differences are the plane pair everywhere a
    word was, tables of 3^t entries walked one trit at a time, and a
    Winograd schedule that keeps its signs.
*/

#ifndef GF3_T
# define GF3_T 5
#endif
#if GF3_T == 4
# define GF3_ENTRIES 81
#elif GF3_T == 5
# define GF3_ENTRIES 243
#elif GF3_T == 6
# define GF3_ENTRIES 729
#else
# error "GF3_T must be 4, 5 or 6"
#endif

/*
    Trit-decomposition tables built and consumed together. Unlike the
    modulus-2 absorb, whose XOR chain the compiler reassociates freely,
    the six-op adder chain is serial, and the absorb's logic-op count
    is invariant in GF3_NT; measurements put 2 ahead of 4 and 8.
*/
#ifndef GF3_NT
# define GF3_NT 2
#endif

/*
    Leaf size for the bit-packed Strassen at modulus 3. Unlike the
    modulus-2 case the recursion is slow to pay here: the M4R table
    build is O(dim^2) work that grows 7/4 per level while the absorb
    only shrinks 7/8, and the leaf loses a little efficiency each time
    its row width halves. On the large-L3 machine this was tuned on the
    leaf alone won through 8192; hosts with small last-level caches may
    prefer 4096.
    (Default 8192; the runtime copy is defined above the packed
    region.)
*/

typedef struct { ulong p1, p2; } gf3_pair;

FLINT_FORCE_INLINE gf3_pair
gf3_add(gf3_pair a, gf3_pair b)
{
    ulong t0 = a.p1 ^ b.p1;
    ulong t1 = a.p2 ^ b.p2;
    gf3_pair s;

    s.p1 = (t0 ^ a.p2) & ~t1;
    s.p2 = (t1 ^ a.p1) & ~t0;

    return s;
}

/*
    Row conversions between canonical mod-3 bytes and the two planes.
    Canonical values make the planes pure bit selections -- bit 0 of the
    byte is the plane-1 indicator and bit 1 the plane-2 indicator -- so
    packing is two mask tests per 64 bytes on AVX-512BW, two
    compare-and-movemask pairs per 32 bytes on AVX2, and two of the
    same multiply gathers as the modulus-2 case portably. Unpacking
    rebuilds the byte as p1 + 2 p2 from disjoint planes.

    Word tails beyond n bits are zeroed; the padded layout relies on it.
*/
FLINT_FORCE_INLINE void
gf3_pack_row(ulong * d1, ulong * d2, const uint8_t * src, slong n,
             slong words)
{
    slong j = 0, w = 0;

#if defined(NMOD_MAT_U8_GF2_AVX512)
    for (; j + 64 <= n; j += 64)
    {
        __m512i v = _mm512_loadu_si512((const void *) (src + j));

        d1[w] = (ulong) _mm512_test_epi8_mask(v, _mm512_set1_epi8(1));
        d2[w] = (ulong) _mm512_test_epi8_mask(v, _mm512_set1_epi8(2));
        w++;
    }
#elif defined(NMOD_MAT_U8_X86)
    for (; j + 64 <= n; j += 64)
    {
        __m256i lo = _mm256_loadu_si256((const __m256i *) (src + j));
        __m256i hi = _mm256_loadu_si256((const __m256i *) (src + j + 32));
        __m256i c1 = _mm256_set1_epi8(1);
        __m256i c2 = _mm256_set1_epi8(2);
        ulong m1lo, m1hi, m2lo, m2hi;

        m1lo = (unsigned int) _mm256_movemask_epi8(
                                  _mm256_cmpeq_epi8(lo, c1));
        m1hi = (unsigned int) _mm256_movemask_epi8(
                                  _mm256_cmpeq_epi8(hi, c1));
        m2lo = (unsigned int) _mm256_movemask_epi8(
                                  _mm256_cmpeq_epi8(lo, c2));
        m2hi = (unsigned int) _mm256_movemask_epi8(
                                  _mm256_cmpeq_epi8(hi, c2));
        d1[w] = m1lo | (m1hi << 32);
        d2[w] = m2lo | (m2hi << 32);
        w++;
    }
#else
    for (; j + 64 <= n; j += 64)
    {
        ulong b1 = 0, b2 = 0;
        slong b;

        for (b = 0; b < 64; b += 8)
        {
            ulong x;

            memcpy(&x, src + j + b, 8);
            b1 |= (((x & UWORD(0x0101010101010101))
                    * UWORD(0x0102040810204080)) >> 56) << b;
            b2 |= ((((x >> 1) & UWORD(0x0101010101010101))
                    * UWORD(0x0102040810204080)) >> 56) << b;
        }

        d1[w] = b1;
        d2[w] = b2;
        w++;
    }
#endif

    if (j < n)
    {
        ulong b1 = 0, b2 = 0;
        slong b;

        for (b = 0; j < n; j++, b++)
        {
            b1 |= (ulong) (src[j] & 1) << b;
            b2 |= (ulong) ((src[j] >> 1) & 1) << b;
        }

        d1[w] = b1;
        d2[w] = b2;
        w++;
    }

    for (; w < words; w++)
        d1[w] = d2[w] = 0;
}

FLINT_FORCE_INLINE void
gf3_unpack_row(uint8_t * dst, const ulong * s1, const ulong * s2, slong n)
{
    slong j = 0;

#if defined(NMOD_MAT_U8_GF2_AVX512)
    for (; j + 64 <= n; j += 64)
    {
        __m512i v1 = _mm512_maskz_set1_epi8((__mmask64) s1[j / 64], 1);
        __m512i v2 = _mm512_maskz_set1_epi8((__mmask64) s2[j / 64], 2);

        _mm512_storeu_si512((void *) (dst + j), _mm512_or_si512(v1, v2));
    }
#elif defined(NMOD_MAT_U8_X86)
    {
        const __m256i sel = _mm256_setr_epi8(
            0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,
            2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3);
        const __m256i bit = _mm256_setr_epi8(
            1, 2, 4, 8, 16, 32, 64, (char) 128,
            1, 2, 4, 8, 16, 32, 64, (char) 128,
            1, 2, 4, 8, 16, 32, 64, (char) 128,
            1, 2, 4, 8, 16, 32, 64, (char) 128);
        const __m256i one = _mm256_set1_epi8(1);

        for (; j + 64 <= n; j += 64)
        {
            ulong x1 = s1[j / 64], x2 = s2[j / 64];
            __m256i lo1 = _mm256_set1_epi32((int) (uint32_t) x1);
            __m256i hi1 = _mm256_set1_epi32((int) (uint32_t) (x1 >> 32));
            __m256i lo2 = _mm256_set1_epi32((int) (uint32_t) x2);
            __m256i hi2 = _mm256_set1_epi32((int) (uint32_t) (x2 >> 32));
            __m256i lo, hi;

            lo1 = _mm256_min_epu8(_mm256_and_si256(
                      _mm256_shuffle_epi8(lo1, sel), bit), one);
            hi1 = _mm256_min_epu8(_mm256_and_si256(
                      _mm256_shuffle_epi8(hi1, sel), bit), one);
            lo2 = _mm256_min_epu8(_mm256_and_si256(
                      _mm256_shuffle_epi8(lo2, sel), bit), one);
            hi2 = _mm256_min_epu8(_mm256_and_si256(
                      _mm256_shuffle_epi8(hi2, sel), bit), one);
            lo = _mm256_or_si256(lo1, _mm256_add_epi8(lo2, lo2));
            hi = _mm256_or_si256(hi1, _mm256_add_epi8(hi2, hi2));
            _mm256_storeu_si256((__m256i *) (dst + j), lo);
            _mm256_storeu_si256((__m256i *) (dst + j + 32), hi);
        }
    }
#else
    for (; j + 64 <= n; j += 64)
    {
        ulong x1 = s1[j / 64], x2 = s2[j / 64];
        slong b;

        for (b = 0; b < 64; b += 8)
        {
            ulong t1 = ((x1 >> b) & 0xff) * UWORD(0x0101010101010101);
            ulong t2 = ((x2 >> b) & 0xff) * UWORD(0x0101010101010101);
            ulong t;

            t1 &= UWORD(0x8040201008040201);
            t2 &= UWORD(0x8040201008040201);
            t1 = ((t1 + UWORD(0x7f7f7f7f7f7f7f7f)) >> 7)
                 & UWORD(0x0101010101010101);
            t2 = ((t2 + UWORD(0x7f7f7f7f7f7f7f7f)) >> 7)
                 & UWORD(0x0101010101010101);
            t = t1 | (t2 << 1);
            memcpy(dst + j + b, &t, 8);
        }
    }
#endif

    for (; j < n; j++)
        dst[j] = (uint8_t) (((s1[j / 64] >> (j % 64)) & 1)
                            + 2 * ((s2[j / 64] >> (j % 64)) & 1));
}

/*
    A matrix of plane pairs: d1 and d2 address the two planes with a
    common row stride. Negation is gf3_view_neg, which swaps the plane
    pointers and costs nothing, so a subtraction is an addition of the
    negated view.
*/
typedef struct
{
    ulong * d1;
    ulong * d2;
    slong r;      /* rows */
    slong w;      /* words per row per plane */
    slong stride; /* words between rows within a plane */
}
gf3_view;

FLINT_FORCE_INLINE gf3_view
gf3_sub_view(gf3_view v, slong r0, slong w0, slong r, slong w)
{
    gf3_view s;

    s.d1 = v.d1 + r0 * v.stride + w0;
    s.d2 = v.d2 + r0 * v.stride + w0;
    s.r = r;
    s.w = w;
    s.stride = v.stride;

    return s;
}

FLINT_FORCE_INLINE gf3_view
gf3_view_neg(gf3_view v)
{
    gf3_view s = v;

    s.d1 = v.d2;
    s.d2 = v.d1;

    return s;
}

/*
    C = A + B elementwise. C may alias A exactly (same base and
    stride), which the schedule uses; the loop reads both operands of
    an iteration before storing, and carries no restrict for that
    reason.
*/
static void
gf3_view_add(gf3_view C, gf3_view A, gf3_view B)
{
    slong i, j;

    for (i = 0; i < C.r; i++)
    {
        ulong * c1 = C.d1 + i * C.stride;
        ulong * c2 = C.d2 + i * C.stride;
        const ulong * a1 = A.d1 + i * A.stride;
        const ulong * a2 = A.d2 + i * A.stride;
        const ulong * b1 = B.d1 + i * B.stride;
        const ulong * b2 = B.d2 + i * B.stride;

        for (j = 0; j < C.w; j++)
        {
            ulong t0 = a1[j] ^ b1[j];
            ulong t1 = a2[j] ^ b2[j];
            ulong r1 = (t0 ^ a2[j]) & ~t1;
            ulong r2 = (t1 ^ a1[j]) & ~t0;

            c1[j] = r1;
            c2[j] = r2;
        }
    }
}

#define gf3_view_sub(C, A, B) gf3_view_add((C), (A), gf3_view_neg(B))

/*
    Leaf: multi-table M4R over the planes. The structure is the
    modulus-2 leaf's: GF3_NT tables of all 3^GF3_T combinations of
    GF3_T rows of B are built per group -- one plane-pair addition per
    entry, walking indices so e and e - 3^d differ in one trit -- and a
    row of C absorbs GF3_NT * GF3_T columns of A in one fixed-shape
    pass, every lookup XORed in unconditionally (entry 0 is the zero
    row). Only the k tail, fewer than GF3_T digits, runs the one-table
    variable path.

    A's digits are unpacked once per leaf into a byte scratch: GF3_T
    does not divide 64, so packed trit groups straddle word boundaries,
    and one vectorized unpack pass is cheaper than per-group bit
    surgery. B and C never leave the packed domain.
*/
/*
    Table indices for every GF3_T-digit group of A, extracted once per
    leaf: the digits pass through a row-sized buffer (so the vectorized
    unpack still does the plane surgery) and the base-3 horner step
    runs here rather than once per absorb pass. The absorb then reads
    two uint16 per lookup, which also avoids an m x kbits byte scratch
    -- at 8192 that scratch was 64 MB of page faults and re-reads.
    Trailing partial groups read zero digits from the padded words.

    Split from the multiply so the threaded path can run it once over
    shared packed A and hand every tile the same array.
*/
static void
gf3_leaf_indices(gf3_view A, uint16_t * aidx, uint8_t * rowbuf)
{
    slong i, g, kbits = A.w * 64;
    slong ng = (kbits + GF3_T - 1) / GF3_T;

    for (i = 0; i < A.r; i++)
    {
        uint16_t * xi = aidx + i * ng;

        gf3_unpack_row(rowbuf, A.d1 + i * A.stride, A.d2 + i * A.stride,
                       kbits);
        memset(rowbuf + kbits, 0, GF3_T);

        for (g = 0; g < ng; g++)
        {
            const uint8_t * ad = rowbuf + g * GF3_T;
            slong idx = 0, l;

            for (l = GF3_T - 1; l >= 0; l--)
                idx = 3 * idx + ad[l];

            xi[g] = (uint16_t) idx;
        }
    }
}

/*
    Leaf: multi-table M4R over the planes, digits of A precomputed as
    table indices. The structure is the modulus-2 leaf's: GF3_NT tables
    of all 3^GF3_T combinations of GF3_T rows of B are built per group
    -- one plane-pair addition per entry, walking indices so e and
    e - 3^d differ in one trit -- and a row of C absorbs
    GF3_NT * GF3_T columns of A in one fixed-shape pass, every lookup
    added in unconditionally (entry 0 is the zero row). Only the k
    tail, fewer than GF3_T digits, runs the one-table variable path.
*/
static void
gf3_leaf_mul(gf3_view C, gf3_view A, gf3_view B, ulong * table,
             const uint16_t * aidx)
{
    slong kg, i, j, t, e, g, kbits = A.w * 64;
    slong ng = (kbits + GF3_T - 1) / GF3_T;
    slong w = C.w;
    slong src3[GF3_ENTRIES], row3[GF3_ENTRIES], pow3[GF3_T + 1];

    pow3[0] = 1;
    for (e = 1; e <= GF3_T; e++)
        pow3[e] = 3 * pow3[e - 1];

    for (e = 1; e < GF3_ENTRIES; e++)
    {
        slong d = 0;

        while ((e / pow3[d]) % 3 == 0)
            d++;

        row3[e] = d;
        src3[e] = e - pow3[d];
    }

    for (i = 0; i < C.r; i++)
        for (j = 0; j < w; j++)
            C.d1[i * C.stride + j] = C.d2[i * C.stride + j] = 0;

    for (kg = 0, g = 0; kg + GF3_T * GF3_NT <= kbits;
         kg += GF3_T * GF3_NT, g += GF3_NT)
    {
        for (t = 0; t < GF3_NT; t++)
        {
            ulong * tb = table + t * GF3_ENTRIES * 2 * w;
            slong kb = kg + GF3_T * t;

            for (j = 0; j < 2 * w; j++)
                tb[j] = 0;

            for (e = 1; e < GF3_ENTRIES; e++)
            {
                ulong * dst = tb + e * 2 * w;
                const ulong * sp = tb + src3[e] * 2 * w;
                const ulong * q1 = B.d1 + (kb + row3[e]) * B.stride;
                const ulong * q2 = B.d2 + (kb + row3[e]) * B.stride;

                for (j = 0; j < w; j++)
                {
                    gf3_pair a, b, r;

                    a.p1 = sp[j]; a.p2 = sp[w + j];
                    b.p1 = q1[j]; b.p2 = q2[j];
                    r = gf3_add(a, b);
                    dst[j] = r.p1; dst[w + j] = r.p2;
                }
            }
        }

        for (i = 0; i < C.r; i++)
        {
            const uint16_t * xi = aidx + i * ng + g;
            ulong * c1 = C.d1 + i * C.stride;
            ulong * c2 = C.d2 + i * C.stride;
            const ulong * s[GF3_NT];

            for (t = 0; t < GF3_NT; t++)
                s[t] = table + t * GF3_ENTRIES * 2 * w
                     + (slong) xi[t] * 2 * w;

            for (j = 0; j < w; j++)
            {
                gf3_pair a, b;

                a.p1 = c1[j];
                a.p2 = c2[j];

                for (t = 0; t < GF3_NT; t++)
                {
                    b.p1 = s[t][j];
                    b.p2 = s[t][w + j];
                    a = gf3_add(a, b);
                }

                c1[j] = a.p1;
                c2[j] = a.p2;
            }
        }
    }

    /* k tail: one table of up to GF3_T digits at a time */
    for (; kg < kbits; kg += GF3_T, g++)
    {
        slong td = FLINT_MIN((slong) GF3_T, kbits - kg);
        slong entries = pow3[td];

        for (j = 0; j < 2 * w; j++)
            table[j] = 0;

        for (e = 1; e < entries; e++)
        {
            ulong * dst = table + e * 2 * w;
            const ulong * sp = table + src3[e] * 2 * w;
            const ulong * q1 = B.d1 + (kg + row3[e]) * B.stride;
            const ulong * q2 = B.d2 + (kg + row3[e]) * B.stride;

            for (j = 0; j < w; j++)
            {
                gf3_pair a, b, r;

                a.p1 = sp[j]; a.p2 = sp[w + j];
                b.p1 = q1[j]; b.p2 = q2[j];
                r = gf3_add(a, b);
                dst[j] = r.p1; dst[w + j] = r.p2;
            }
        }

        for (i = 0; i < C.r; i++)
        {
            const ulong * sp;
            ulong * c1 = C.d1 + i * C.stride;
            ulong * c2 = C.d2 + i * C.stride;
            slong idx = aidx[i * ng + g];

            if (idx == 0)
                continue;

            sp = table + idx * 2 * w;

            for (j = 0; j < w; j++)
            {
                gf3_pair a, b, r;

                a.p1 = c1[j]; a.p2 = c2[j];
                b.p1 = sp[j]; b.p2 = sp[w + j];
                r = gf3_add(a, b);
                c1[j] = r.p1; c2[j] = r.p2;
            }
        }
    }
}

static void
gf3_mul_leaf(gf3_view C, gf3_view A, gf3_view B, ulong * table,
             uint16_t * aidx, uint8_t * rowbuf)
{
    gf3_leaf_indices(A, aidx, rowbuf);
    gf3_leaf_mul(C, A, B, table, aidx);
}

/*
    The same Winograd schedule as the byte-domain recursion, signs kept,
    with subtractions realized as additions of plane-swapped views.
    Operands are padded at pack time so every split is even.
*/
static void
gf3_strassen(gf3_view C, gf3_view A, gf3_view B, slong levels,
             ulong * scratch, ulong * table, uint16_t * aidx,
             uint8_t * rowbuf)
{
    slong anr, anw, bnr, bnw, x1w;
    gf3_view A11, A12, A21, A22, B11, B12, B21, B22;
    gf3_view C11, C12, C21, C22, X1a, X1b, X2;
    ulong * child;

    if (levels == 0)
    {
        gf3_mul_leaf(C, A, B, table, aidx, rowbuf);
        return;
    }

    anr = A.r / 2; anw = A.w / 2;
    bnr = B.r / 2; bnw = B.w / 2;

    A11 = gf3_sub_view(A, 0, 0, anr, anw);
    A12 = gf3_sub_view(A, 0, anw, anr, anw);
    A21 = gf3_sub_view(A, anr, 0, anr, anw);
    A22 = gf3_sub_view(A, anr, anw, anr, anw);

    B11 = gf3_sub_view(B, 0, 0, bnr, bnw);
    B12 = gf3_sub_view(B, 0, bnw, bnr, bnw);
    B21 = gf3_sub_view(B, bnr, 0, bnr, bnw);
    B22 = gf3_sub_view(B, bnr, bnw, bnr, bnw);

    C11 = gf3_sub_view(C, 0, 0, anr, bnw);
    C12 = gf3_sub_view(C, 0, bnw, anr, bnw);
    C21 = gf3_sub_view(C, anr, 0, anr, bnw);
    C22 = gf3_sub_view(C, anr, bnw, anr, bnw);

    x1w = FLINT_MAX(anw, bnw);

    X1a.d1 = scratch;
    X1a.d2 = scratch + anr * x1w;
    X1a.r = anr; X1a.w = anw; X1a.stride = x1w;
    X1b = X1a; X1b.w = bnw;
    X2.d1 = scratch + 2 * anr * x1w;
    X2.d2 = X2.d1 + bnr * bnw;
    X2.r = bnr; X2.w = bnw; X2.stride = bnw;
    child = X2.d2 + bnr * bnw;

    gf3_view_add(X1a, A22, A12);
    gf3_view_add(X2, B22, B12);
    gf3_strassen(C21, X1a, X2, levels - 1, child, table, aidx, rowbuf);

    gf3_view_sub(X1a, A22, A21);
    gf3_view_sub(X2, B22, B21);
    gf3_strassen(C22, X1a, X2, levels - 1, child, table, aidx, rowbuf);

    gf3_view_add(X1a, X1a, A12);
    gf3_view_add(X2, X2, B12);
    gf3_strassen(C11, X1a, X2, levels - 1, child, table, aidx, rowbuf);

    gf3_view_sub(X1a, X1a, A11);
    gf3_strassen(C12, X1a, B12, levels - 1, child, table, aidx, rowbuf);

    gf3_strassen(X1b, A12, B21, levels - 1, child, table, aidx, rowbuf);
    gf3_view_add(C11, C11, X1b);
    gf3_view_add(C12, C12, C22);
    gf3_view_sub(C12, C11, C12);
    gf3_view_sub(C11, C21, C11);
    gf3_view_sub(X2, X2, B11);
    gf3_strassen(C21, A21, X2, levels - 1, child, table, aidx, rowbuf);
    gf3_view_sub(C21, C11, C21);
    gf3_view_add(C22, C22, C11);
    gf3_strassen(C11, A11, B11, levels - 1, child, table, aidx, rowbuf);
    gf3_view_add(C11, C11, X1b);
}

/* X1 carries two planes of anr x max(anw, bnw) and X2 two of
   bnr x bnw; bnr is the inner dimension in rows, 64 * anw */
static slong
gf3_scratch_needed(slong ar, slong aw, slong br, slong bw, slong levels)
{
    slong anr, anw, bnr, bnw, x1w;

    if (levels == 0)
        return 0;

    anr = ar / 2; anw = aw / 2;
    bnr = br / 2; bnw = bw / 2;
    x1w = FLINT_MAX(anw, bnw);

    return 2 * (anr * x1w + bnr * bnw)
        + gf3_scratch_needed(anr, anw, bnr, bnw, levels - 1);
}

/* entry point: pack with the padding the recursion needs, recurse,
   unpack */
static void
gf3_mul_strassen(uint8_t * C, slong Cstride,
                 const uint8_t * A, slong Astride,
                 const uint8_t * B, slong Bstride,
                 slong m, slong k, slong n, slong levels, int ulsrc)
{
    slong unit = WORD(1) << levels;
    slong mp = ((m + unit - 1) / unit) * unit;
    slong aw = (((k + 63) / 64 + unit - 1) / unit) * unit;
    slong bw = (((n + 63) / 64 + unit - 1) / unit) * unit;
    slong kp = aw * 64;
    slong i, need, tw, words_total;
    int same;
    gf3_view Av, Bv, Cv;
    ulong * mem;
    ulong * table;
    ulong * scratch;
    uint16_t * aidx;
    uint8_t * rowbuf;
    TMP_INIT;

    need = gf3_scratch_needed(mp, aw, kp, bw, levels);
    tw = bw >> levels;    /* leaf C width: table row width */

    /* squaring at the leaf: as in the modulus-2 entry, pack once with
       the padded B-side row count and share the planes */
    same = (levels == 0 && A == B && Astride == Bstride
            && m == k && k == n);

    words_total = 2 * ((same ? kp * aw : mp * aw + kp * bw) + mp * bw)
                  + (slong) GF3_NT * GF3_ENTRIES * 2 * tw + need;

    TMP_START;
    {
        slong lkbits = (aw >> levels) * 64;
        slong lng = (lkbits + GF3_T - 1) / GF3_T;

        mem = (ulong *) TMP_ALLOC(words_total * sizeof(ulong)
                                  + (mp >> levels) * lng
                                    * (slong) sizeof(uint16_t)
                                  + lkbits + GF3_T);
        aidx = (uint16_t *) (mem + words_total);
        rowbuf = (uint8_t *) (aidx + (mp >> levels) * lng);
    }

    Av.d1 = mem;
    Av.d2 = Av.d1 + (same ? kp : mp) * aw;
    Av.r = mp; Av.w = aw; Av.stride = aw;
    Bv.d1 = same ? Av.d1 : Av.d2 + mp * aw;
    Bv.d2 = same ? Av.d2 : Bv.d1 + kp * bw;
    Bv.r = kp; Bv.w = bw; Bv.stride = bw;
    Cv.d1 = Bv.d2 + kp * bw;
    Cv.d2 = Cv.d1 + mp * bw;
    Cv.r = mp; Cv.w = bw; Cv.stride = bw;
    table = Cv.d2 + mp * bw;
    scratch = table + (slong) GF3_NT * GF3_ENTRIES * 2 * tw;

    for (i = 0; i < m; i++)
    {
        if (ulsrc)
            gf3_pack_row_ul(Av.d1 + i * aw, Av.d2 + i * aw,
                            (const ulong *) A + i * Astride, k, aw);
        else
            gf3_pack_row(Av.d1 + i * aw, Av.d2 + i * aw,
                         A + i * Astride, k, aw);
    }

    if (same)
    {
        memset(Av.d1 + m * aw, 0, (kp - m) * aw * sizeof(ulong));
        memset(Av.d2 + m * aw, 0, (kp - m) * aw * sizeof(ulong));
    }
    else
    {
        memset(Av.d1 + m * aw, 0, (mp - m) * aw * sizeof(ulong));
        memset(Av.d2 + m * aw, 0, (mp - m) * aw * sizeof(ulong));

        for (i = 0; i < k; i++)
        {
            if (ulsrc)
                gf3_pack_row_ul(Bv.d1 + i * bw, Bv.d2 + i * bw,
                                (const ulong *) B + i * Bstride, n, bw);
            else
                gf3_pack_row(Bv.d1 + i * bw, Bv.d2 + i * bw,
                             B + i * Bstride, n, bw);
        }
        memset(Bv.d1 + k * bw, 0, (kp - k) * bw * sizeof(ulong));
        memset(Bv.d2 + k * bw, 0, (kp - k) * bw * sizeof(ulong));
    }

    gf3_strassen(Cv, Av, Bv, levels, scratch, table, aidx, rowbuf);

    for (i = 0; i < m; i++)
    {
        if (ulsrc)
            gf3_unpack_row_ul((ulong *) C + i * Cstride,
                              Cv.d1 + i * bw, Cv.d2 + i * bw, n);
        else
            gf3_unpack_row(C + i * Cstride, Cv.d1 + i * bw,
                           Cv.d2 + i * bw, n);
    }

    TMP_END;
}



/* ------------------------------------------------------------------ */
/* Strassen in the packed representation                               */
/* ------------------------------------------------------------------ */

/*
    Strassen entirely in the bit-packed domain, which is where the
    packing density pays: an addition of two quadrants is one XOR per 64
    elements at modulus 2, so the O(n^2) part of the recursion costs a
    sixty-fourth of what it costs on bytes, while the multiplication
    saved is the same 1/8 per level. Doing the recursion on bytes and
    only the leaves packed -- the obvious shortcut -- throws that away.

    Operands are padded at pack time so that every split lands on a word
    boundary and every quadrant is exactly even: rows to a multiple of
    2^levels, columns to a multiple of 64 * 2^levels. No peeling and no
    ragged views anywhere, at a cost of at most a factor 2^levels/n in
    memory.
*/

typedef struct
{
    ulong * d;
    slong r;      /* rows */
    slong w;      /* words per row */
    slong stride; /* words between rows */
}
gf2_view;

FLINT_FORCE_INLINE gf2_view
gf2_sub_view(gf2_view v, slong r0, slong w0, slong r, slong w)
{
    gf2_view s;

    s.d = v.d + r0 * v.stride + w0;
    s.r = r;
    s.w = w;
    s.stride = v.stride;

    return s;
}

static void
gf2_view_add(gf2_view C, gf2_view A, gf2_view B)
{
    slong i, j;

    /* C may alias A exactly, which the schedule uses; each iteration
       reads before it writes, so no restrict */
    for (i = 0; i < C.r; i++)
    {
        ulong * c = C.d + i * C.stride;
        const ulong * a = A.d + i * A.stride;
        const ulong * b = B.d + i * B.stride;

        for (j = 0; j < C.w; j++)
            c[j] = a[j] ^ b[j];
    }
}

/*
    Leaf: multi-table M4R. Several Gray-code tables are built at once --
    GF2_NT groups of GF2_T columns of A -- and all of their lookups are
    accumulated into a row of C in a single pass, so C is read and
    written once per GF2_NT * GF2_T columns rather than once per
    GF2_T. Per row and per group of tables the traffic falls from
    3 * NT row-widths (a table read plus a read and a write of C, NT
    times) to NT + 2, which is 2.4x less at NT = 8.

    This is what M4RI does, and measuring against it is what identified
    the omission: without it this kernel was 5-12x off M4RI, and the
    table-order sweep behaved as though additions were not the cost,
    which they are not -- the traffic on C is.

    Full groups take a branchless path: every one of the GF2_NT lookups
    is XORed in whether or not its byte of A is zero, since entry 0 of
    each table is the zero row. With dense operands a zero byte has
    probability 2^-8, so the test bought nothing, and routing a
    variable number of sources through a switch meant that any nt
    beyond the largest unrolled case fell to a scalar loop -- which is
    precisely why GF2_NT = 8 measured 2x SLOWER than 4 before this
    form, when the traffic argument says it should win. With the fixed
    shape the compiler unrolls the source XORs and vectorizes, and
    NT = 8 takes its predicted place.

    A is read as bytes out of its packed rows, exact because GF2_T is 8
    and the splits keep everything byte aligned.
*/
static void
gf2_mul_leaf(gf2_view C, gf2_view A, gf2_view B, ulong * table)
{
    slong kg, i, j, t, kbits = A.w * 64;
    slong w = C.w;

    for (i = 0; i < C.r; i++)
        for (j = 0; j < w; j++)
            C.d[i * C.stride + j] = 0;

    for (kg = 0; kg < kbits; kg += 8 * GF2_NT)
    {
        slong nt = FLINT_MIN((slong) GF2_NT, (kbits - kg) / 8);

        for (t = 0; t < nt; t++)
        {
            ulong * tb = table + t * 256 * w;
            slong kb = kg + 8 * t;

            for (j = 0; j < w; j++)
                tb[j] = 0;

            for (i = 1; i < 256; i++)
            {
                slong bit = 0, v = i;
                ulong * restrict dst = tb + i * w;
                const ulong * restrict src;
                const ulong * restrict add;

                while (!(v & 1))
                {
                    v >>= 1;
                    bit++;
                }

                src = tb + (i & (i - 1)) * w;
                add = B.d + (kb + bit) * B.stride;

                for (j = 0; j < w; j++)
                    dst[j] = src[j] ^ add[j];
            }
        }

        if (nt == GF2_NT)
        {
            for (i = 0; i < C.r; i++)
            {
                const unsigned char * ab =
                    (const unsigned char *) (A.d + i * A.stride) + kg / 8;
                ulong * restrict dst = C.d + i * C.stride;
                const ulong * s[GF2_NT];

                for (t = 0; t < GF2_NT; t++)
                    s[t] = table + t * 256 * w + (slong) ab[t] * w;

                for (j = 0; j < w; j++)
                {
                    ulong x = dst[j];

                    for (t = 0; t < GF2_NT; t++)
                        x ^= s[t][j];

                    dst[j] = x;
                }
            }

            continue;
        }

        /* k tail: fewer than GF2_NT tables were built */
        for (i = 0; i < C.r; i++)
        {
            const unsigned char * ab =
                (const unsigned char *) (A.d + i * A.stride) + kg / 8;
            ulong * restrict dst = C.d + i * C.stride;
            const ulong * src[GF2_NT];
            slong used = 0;

            src[0] = NULL;   /* silences a maybe-uninitialized warning
                                at GF2_NT = 1, where the compiler cannot
                                see that used > 0 implies src[0] set */

            for (t = 0; t < nt; t++)
                if (ab[t] != 0)
                    src[used++] = table + t * 256 * w + (slong) ab[t] * w;

            /* one pass over the row of C for the whole group */
            switch (used)
            {
                case 0:
                    break;
                case 1:
                {
                    const ulong * restrict s0 = src[0];

                    for (j = 0; j < w; j++)
                        dst[j] ^= s0[j];
                    break;
                }
#if GF2_NT >= 2
                case 2:
                {
                    const ulong * restrict s0 = src[0];
                    const ulong * restrict s1 = src[1];

                    for (j = 0; j < w; j++)
                        dst[j] ^= s0[j] ^ s1[j];
                    break;
                }
#endif
#if GF2_NT >= 3
                case 3:
                {
                    const ulong * restrict s0 = src[0];
                    const ulong * restrict s1 = src[1];
                    const ulong * restrict s2 = src[2];

                    for (j = 0; j < w; j++)
                        dst[j] ^= s0[j] ^ s1[j] ^ s2[j];
                    break;
                }
#endif
#if GF2_NT >= 4
                case 4:
                {
                    const ulong * restrict s0 = src[0];
                    const ulong * restrict s1 = src[1];
                    const ulong * restrict s2 = src[2];
                    const ulong * restrict s3 = src[3];

                    for (j = 0; j < w; j++)
                        dst[j] ^= s0[j] ^ s1[j] ^ s2[j] ^ s3[j];
                    break;
                }
#endif
                default:
                    for (j = 0; j < w; j++)
                    {
                        ulong x = dst[j];
                        slong u;

                        for (u = 0; u < used; u++)
                            x ^= src[u][j];

                        dst[j] = x;
                    }
            }
        }
    }
}

static void
gf2_strassen(gf2_view C, gf2_view A, gf2_view B, slong levels,
             ulong * scratch, ulong * table)
{
    slong anr, anw, bnr, bnw, x1w;
    gf2_view A11, A12, A21, A22, B11, B12, B21, B22;
    gf2_view C11, C12, C21, C22, X1a, X1b, X2;
    ulong * child;

    if (levels == 0)
    {
        gf2_mul_leaf(C, A, B, table);
        return;
    }

    anr = A.r / 2; anw = A.w / 2;
    bnr = B.r / 2; bnw = B.w / 2;

    A11 = gf2_sub_view(A, 0, 0, anr, anw);
    A12 = gf2_sub_view(A, 0, anw, anr, anw);
    A21 = gf2_sub_view(A, anr, 0, anr, anw);
    A22 = gf2_sub_view(A, anr, anw, anr, anw);

    B11 = gf2_sub_view(B, 0, 0, bnr, bnw);
    B12 = gf2_sub_view(B, 0, bnw, bnr, bnw);
    B21 = gf2_sub_view(B, bnr, 0, bnr, bnw);
    B22 = gf2_sub_view(B, bnr, bnw, bnr, bnw);

    C11 = gf2_sub_view(C, 0, 0, anr, bnw);
    C12 = gf2_sub_view(C, 0, bnw, anr, bnw);
    C21 = gf2_sub_view(C, anr, 0, anr, bnw);
    C22 = gf2_sub_view(C, anr, bnw, anr, bnw);

    x1w = FLINT_MAX(anw, bnw);

    X1a.d = scratch; X1a.r = anr; X1a.w = anw; X1a.stride = x1w;
    X1b.d = scratch; X1b.r = anr; X1b.w = bnw; X1b.stride = x1w;
    X2.d = scratch + anr * x1w; X2.r = bnr; X2.w = bnw; X2.stride = bnw;
    child = X2.d + bnr * bnw;

    /* Bodrato's schedule; over GF(2) addition and subtraction coincide */
    gf2_view_add(X1a, A22, A12);
    gf2_view_add(X2, B22, B12);
    gf2_strassen(C21, X1a, X2, levels - 1, child, table);

    gf2_view_add(X1a, A22, A21);
    gf2_view_add(X2, B22, B21);
    gf2_strassen(C22, X1a, X2, levels - 1, child, table);

    gf2_view_add(X1a, X1a, A12);
    gf2_view_add(X2, X2, B12);
    gf2_strassen(C11, X1a, X2, levels - 1, child, table);

    gf2_view_add(X1a, X1a, A11);
    gf2_strassen(C12, X1a, B12, levels - 1, child, table);

    gf2_strassen(X1b, A12, B21, levels - 1, child, table);
    gf2_view_add(C11, C11, X1b);
    gf2_view_add(C12, C12, C22);
    gf2_view_add(C12, C11, C12);
    gf2_view_add(C11, C21, C11);
    gf2_view_add(X2, X2, B11);
    gf2_strassen(C21, A21, X2, levels - 1, child, table);
    gf2_view_add(C21, C11, C21);
    gf2_view_add(C22, C22, C11);
    gf2_strassen(C11, A11, B11, levels - 1, child, table);
    gf2_view_add(C11, C11, X1b);
}

/* entry point: pack with the padding the recursion needs, recurse,
   unpack */
static void
gf2_mul_strassen(uint8_t * C, slong Cstride,
                 const uint8_t * A, slong Astride,
                 const uint8_t * B, slong Bstride,
                 slong m, slong k, slong n, slong levels, int ulsrc);

/*
    X1 is anr x max(anw, bnw) words and X2 is bnr x bnw, where bnr is
    the inner dimension in rows -- 64 * anw, not a word count. Sizing X2
    with word counts undersizes it by a factor of 64 and corrupts the
    heap as soon as the recursion is entered.
*/
static slong
gf2_scratch_needed(slong ar, slong aw, slong br, slong bw, slong levels)
{
    slong anr, anw, bnr, bnw, x1w;

    if (levels == 0)
        return 0;

    anr = ar / 2; anw = aw / 2;
    bnr = br / 2; bnw = bw / 2;
    x1w = FLINT_MAX(anw, bnw);

    return anr * x1w + bnr * bnw
        + gf2_scratch_needed(anr, anw, bnr, bnw, levels - 1);
}

static void
gf2_mul_strassen(uint8_t * C, slong Cstride,
                 const uint8_t * A, slong Astride,
                 const uint8_t * B, slong Bstride,
                 slong m, slong k, slong n, slong levels, int ulsrc)
{
    slong unit = WORD(1) << levels;
    slong mp = ((m + unit - 1) / unit) * unit;
    slong aw = (((k + 63) / 64 + unit - 1) / unit) * unit;
    slong bw = (((n + 63) / 64 + unit - 1) / unit) * unit;
    slong kp = aw * 64;
    slong i, need;
    gf2_view Av, Bv, Cv;
    ulong * mem;
    ulong * table;
    ulong * scratch;
    int same;
    TMP_INIT;

    need = gf2_scratch_needed(mp, aw, kp, bw, levels);

    /*
        Squaring at the leaf: the same byte matrix packs to the same
        plane (aw == bw when k == n), so pack once with the taller,
        B-side padded row count and share; with levels > 0 the A and B
        paddings differ and the operands stay separate.
    */
    same = (levels == 0 && A == B && Astride == Bstride
            && m == k && k == n);

    TMP_START;
    mem = (ulong *) TMP_ALLOC(((same ? kp * aw : mp * aw + kp * bw)
                               + mp * bw
                               + (slong) GF2_NT * 256 * bw + need)
                              * sizeof(ulong));

    Av.d = mem;                     Av.r = mp; Av.w = aw; Av.stride = aw;
    Bv.d = same ? mem
                : mem + mp * aw;    Bv.r = kp; Bv.w = bw; Bv.stride = bw;
    Cv.d = Bv.d + kp * bw;          Cv.r = mp; Cv.w = bw; Cv.stride = bw;
    table = Cv.d + mp * bw;
    scratch = table + (slong) GF2_NT * 256 * bw;

    for (i = 0; i < m; i++)
    {
        if (ulsrc)
            gf2_pack_row_ul(Av.d + i * aw,
                            (const ulong *) A + i * Astride, k, aw);
        else
            gf2_pack_row(Av.d + i * aw, A + i * Astride, k, aw);
    }

    if (same)
        memset(Av.d + m * aw, 0, (kp - m) * aw * sizeof(ulong));
    else
    {
        memset(Av.d + m * aw, 0, (mp - m) * aw * sizeof(ulong));

        for (i = 0; i < k; i++)
        {
            if (ulsrc)
                gf2_pack_row_ul(Bv.d + i * bw,
                                (const ulong *) B + i * Bstride, n, bw);
            else
                gf2_pack_row(Bv.d + i * bw, B + i * Bstride, n, bw);
        }
        memset(Bv.d + k * bw, 0, (kp - k) * bw * sizeof(ulong));
    }

    gf2_strassen(Cv, Av, Bv, levels, scratch, table);

    for (i = 0; i < m; i++)
    {
        if (ulsrc)
            gf2_unpack_row_ul((ulong *) C + i * Cstride,
                              Cv.d + i * bw, n);
        else
            gf2_unpack_row(C + i * Cstride, Cv.d + i * bw, n);
    }

    TMP_END;
}

/* recursion depth the packed entries are asked for, from the smallest
   dimension and the per-modulus leaf cutoff */
static slong
packed_strassen_levels(slong m, slong k, slong n, slong cutoff)
{
    slong levels = 0, mind = FLINT_MIN(FLINT_MIN(m, k), n);

    while (levels < 6 && (mind >> (levels + 1)) >= cutoff)
        levels++;

    return levels;
}

/* ------------------------------------------------------------------ */
/* Threaded packed multiplication                                      */
/* ------------------------------------------------------------------ */

/*
    Tile grid for the packed path: rows first, columns only for the
    threads left over. The costs are asymmetric. A row strip keeps the
    leaf at full width and reads a contiguous range of the shared index
    array, so its only price is rebuilding the M4R tables, a 2^GF2_T/m
    resp. 3^GF3_T/m fraction of the absorb per extra strip. A column
    strip narrows the table rows the absorb streams, and measured on
    8192 leaves that costs far more than the model of "same total
    words" suggests: each halving of the width roughly adds 10-30% at
    modulus 2 and 30-70% at modulus 3, where the per-row fixed costs --
    for modulus 3 mainly the index array fetch, one cold cache line
    per row pass -- stop being hidden behind the narrowing inner loop.
    Utilization still comes first: an idle thread costs more than
    either split.
*/
static void
packed_tile_grid(slong * rp, slong * cp, slong m, slong n, slong nthreads,
                 ulong nn)
{
    slong rmax = FLINT_MAX(1, m / U8_PACKED_PAR_MIN_TILE);
    slong cmax = FLINT_MAX(1, n / U8_PACKED_PAR_MIN_TILE);
    slong budget, r, c;

    /*
        Each row strip rebuilds the full-width M4R tables, and the
        build is bandwidth-bound, so on a machine whose cores share
        DRAM the duplicated builds start eating the very bandwidth the
        absorb needs: at modulus 2, 4096, eight row strips move more
        build traffic than the absorb itself and 8 threads measured
        barely ahead of 4. The cap keeps the duplicated build a
        fraction of the absorb (the divisors encode build/absorb time
        ratios measured on an 8-core Zen 3; retune with p-mul_u8);
        leftover threads go to column strips, which duplicate nothing.
    */
    budget = 1 + m / ((nn == 2) ? 2048 : 1024);

    r = FLINT_MIN(FLINT_MIN(rmax, nthreads), budget);
    c = FLINT_MIN(cmax, FLINT_MAX(1, nthreads / r));

    *rp = r;
    *cp = c;
}

/*
    Shared-operand tiling. A first version ran the whole serial packed
    multiplication per tile of C, which was simple but repeated the
    conversions: A was packed once per column strip and, worse, the
    modulus-3 index extraction -- a scalar horner over every digit of A
    -- ran once per column strip too, which measured as most of the
    overhead. So the phases share instead:

      1. pack, in parallel over rows: A and B are packed once into
         shared plane storage, and for modulus 3 the A-side table
         indices are extracted once;
      2. multiply, in parallel over a tile grid: each task runs the
         M4R leaf on sub-views of the shared operands, building its
         own tables -- the one duplicated cost, once per row strip,
         bounded by 2^GF2_T/rows resp. 3^GF3_T/rows of the absorb;
      3. unpack, in parallel over rows of C.

    The tiles run the leaf directly, no Strassen: the grid hands out
    tiles at or below the serial cutoffs anyway, and skipping the
    recursion means no power-of-two padding, so the sub-views tile the
    shared operands exactly.

    The k dimension is never split, so no two tasks write the same
    entry of C and the tile phase needs no reduction; the phases are
    separated by the pool's wake/wait.
*/

#ifdef PACKED_PAR_PROF
#include <stdio.h>
#include <time.h>
static double packed_par_prof_now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}
#endif

typedef struct
{
    const uint8_t * A; slong Astride;
    const uint8_t * B; slong Bstride;
    uint8_t * C; slong Cstride;
    slong m, k, n;
    ulong nn;
    ulong * a1; ulong * a2;    /* packed A planes, a2 unused at nn == 2 */
    ulong * b1; ulong * b2;    /* equal to the A planes when squaring */
    int share;                 /* A == B: pack once, share the planes */
    int ulsrc;                 /* operand entries are reduced ulongs
                                  rather than bytes; strides count
                                  elements either way */
    ulong * c1; ulong * c2;
    uint16_t * aidx;           /* nn == 3 */
    slong aw, bw, kp, ng;
    slong rp, cp;
}
packed_par_ctx;

typedef struct
{
    packed_par_ctx * ctx;
    int phase;                 /* 0 pack, 1 multiply, 2 unpack */
    slong start, end;          /* rows in 0 and 2, tiles in 1 */
}
packed_par_job;

/*
    Phase 0 indexes the rows of A and B as one range, A first: a row
    r < m is row r of A, and m <= r < m + kp is row r - m of B, with
    the rows past k being the zero padding the table build reads.
*/
static void
packed_par_pack(packed_par_ctx * ctx, slong start, slong end)
{
    slong r;

    if (ctx->nn == 2)
    {
        const ulong * Aul = (const ulong *) ctx->A;
        const ulong * Bul = (const ulong *) ctx->B;

        for (r = start; r < end; r++)
        {
            if (ctx->share)
            {
                /* one row space 0..kp: rows of A, then zero padding */
                if (r < ctx->k)
                {
                    if (ctx->ulsrc)
                        gf2_pack_row_ul(ctx->a1 + r * ctx->aw,
                                        Aul + r * ctx->Astride, ctx->k,
                                        ctx->aw);
                    else
                        gf2_pack_row(ctx->a1 + r * ctx->aw,
                                     ctx->A + r * ctx->Astride, ctx->k,
                                     ctx->aw);
                }
                else
                    memset(ctx->a1 + r * ctx->aw, 0,
                           ctx->aw * sizeof(ulong));
                continue;
            }

            if (r < ctx->m)
            {
                if (ctx->ulsrc)
                    gf2_pack_row_ul(ctx->a1 + r * ctx->aw,
                                    Aul + r * ctx->Astride, ctx->k,
                                    ctx->aw);
                else
                    gf2_pack_row(ctx->a1 + r * ctx->aw,
                                 ctx->A + r * ctx->Astride, ctx->k,
                                 ctx->aw);
            }
            else if (r - ctx->m < ctx->k)
            {
                if (ctx->ulsrc)
                    gf2_pack_row_ul(ctx->b1 + (r - ctx->m) * ctx->bw,
                                    Bul + (r - ctx->m) * ctx->Bstride,
                                    ctx->n, ctx->bw);
                else
                    gf2_pack_row(ctx->b1 + (r - ctx->m) * ctx->bw,
                                 ctx->B + (r - ctx->m) * ctx->Bstride,
                                 ctx->n, ctx->bw);
            }
            else
                memset(ctx->b1 + (r - ctx->m) * ctx->bw, 0,
                       ctx->bw * sizeof(ulong));
        }
    }
    else
    {
        uint8_t * rowbuf;
        TMP_INIT;

        TMP_START;
        rowbuf = (uint8_t *) TMP_ALLOC(ctx->kp + GF3_T);

        for (r = start; r < end; r++)
        {
            const ulong * Aul = (const ulong *) ctx->A;
            const ulong * Bul = (const ulong *) ctx->B;

            if (ctx->share && r >= ctx->k)
            {
                memset(ctx->a1 + r * ctx->aw, 0, ctx->aw * sizeof(ulong));
                memset(ctx->a2 + r * ctx->aw, 0, ctx->aw * sizeof(ulong));
                continue;
            }

            if (r < ctx->m || ctx->share)
            {
                gf3_view Ar;

                if (ctx->ulsrc)
                    gf3_pack_row_ul(ctx->a1 + r * ctx->aw,
                                    ctx->a2 + r * ctx->aw,
                                    Aul + r * ctx->Astride, ctx->k,
                                    ctx->aw);
                else
                    gf3_pack_row(ctx->a1 + r * ctx->aw,
                                 ctx->a2 + r * ctx->aw,
                                 ctx->A + r * ctx->Astride, ctx->k,
                                 ctx->aw);

                Ar.d1 = ctx->a1 + r * ctx->aw;
                Ar.d2 = ctx->a2 + r * ctx->aw;
                Ar.r = 1; Ar.w = ctx->aw; Ar.stride = ctx->aw;
                gf3_leaf_indices(Ar, ctx->aidx + r * ctx->ng, rowbuf);
            }
            else if (r - ctx->m < ctx->k)
            {
                if (ctx->ulsrc)
                    gf3_pack_row_ul(ctx->b1 + (r - ctx->m) * ctx->bw,
                                    ctx->b2 + (r - ctx->m) * ctx->bw,
                                    Bul + (r - ctx->m) * ctx->Bstride,
                                    ctx->n, ctx->bw);
                else
                    gf3_pack_row(ctx->b1 + (r - ctx->m) * ctx->bw,
                                 ctx->b2 + (r - ctx->m) * ctx->bw,
                                 ctx->B + (r - ctx->m) * ctx->Bstride,
                                 ctx->n, ctx->bw);
            }
            else
            {
                memset(ctx->b1 + (r - ctx->m) * ctx->bw, 0,
                       ctx->bw * sizeof(ulong));
                memset(ctx->b2 + (r - ctx->m) * ctx->bw, 0,
                       ctx->bw * sizeof(ulong));
            }
        }

        TMP_END;
    }
}

static void
packed_par_mul(packed_par_ctx * ctx, slong tile)
{
    slong ri = tile / ctx->cp, ci = tile % ctx->cp;
    slong i0 = (ctx->m * ri) / ctx->rp;
    slong i1 = (ctx->m * (ri + 1)) / ctx->rp;
    slong j0 = (ctx->bw * ci) / ctx->cp;
    slong j1 = (ctx->bw * (ci + 1)) / ctx->cp;
    slong wt = j1 - j0;
    ulong * table;
    TMP_INIT;

    TMP_START;

    if (ctx->nn == 2)
    {
        gf2_view At, Bt, Ct;

        At.d = ctx->a1 + i0 * ctx->aw;
        At.r = i1 - i0; At.w = ctx->aw; At.stride = ctx->aw;
        Bt.d = ctx->b1 + j0;
        Bt.r = ctx->kp; Bt.w = wt; Bt.stride = ctx->bw;
        Ct.d = ctx->c1 + i0 * ctx->bw + j0;
        Ct.r = i1 - i0; Ct.w = wt; Ct.stride = ctx->bw;

        table = (ulong *) TMP_ALLOC((slong) GF2_NT * 256 * wt
                                    * sizeof(ulong));
        gf2_mul_leaf(Ct, At, Bt, table);
    }
    else
    {
        gf3_view At, Bt, Ct;

        At.d1 = ctx->a1 + i0 * ctx->aw;
        At.d2 = ctx->a2 + i0 * ctx->aw;
        At.r = i1 - i0; At.w = ctx->aw; At.stride = ctx->aw;
        Bt.d1 = ctx->b1 + j0;
        Bt.d2 = ctx->b2 + j0;
        Bt.r = ctx->kp; Bt.w = wt; Bt.stride = ctx->bw;
        Ct.d1 = ctx->c1 + i0 * ctx->bw + j0;
        Ct.d2 = ctx->c2 + i0 * ctx->bw + j0;
        Ct.r = i1 - i0; Ct.w = wt; Ct.stride = ctx->bw;

        table = (ulong *) TMP_ALLOC((slong) GF3_NT * GF3_ENTRIES * 2 * wt
                                    * sizeof(ulong));
        gf3_leaf_mul(Ct, At, Bt, table, ctx->aidx + i0 * ctx->ng);
    }

    TMP_END;
}

static void
packed_par_unpack(packed_par_ctx * ctx, slong start, slong end)
{
    slong r;
    ulong * Cul = (ulong *) ctx->C;

    if (ctx->nn == 2)
    {
        for (r = start; r < end; r++)
        {
            if (ctx->ulsrc)
                gf2_unpack_row_ul(Cul + r * ctx->Cstride,
                                  ctx->c1 + r * ctx->bw, ctx->n);
            else
                gf2_unpack_row(ctx->C + r * ctx->Cstride,
                               ctx->c1 + r * ctx->bw, ctx->n);
        }
    }
    else
    {
        for (r = start; r < end; r++)
        {
            if (ctx->ulsrc)
                gf3_unpack_row_ul(Cul + r * ctx->Cstride,
                                  ctx->c1 + r * ctx->bw,
                                  ctx->c2 + r * ctx->bw, ctx->n);
            else
                gf3_unpack_row(ctx->C + r * ctx->Cstride,
                               ctx->c1 + r * ctx->bw,
                               ctx->c2 + r * ctx->bw, ctx->n);
        }
    }
}

static void
packed_par_worker(void * varg)
{
    packed_par_job * job = (packed_par_job *) varg;
    slong i;

    if (job->phase == 0)
        packed_par_pack(job->ctx, job->start, job->end);
    else if (job->phase == 1)
        for (i = job->start; i < job->end; i++)
            packed_par_mul(job->ctx, i);
    else
        packed_par_unpack(job->ctx, job->start, job->end);
}

/*
    Returns 0 without touching C when the grid degenerates or there is
    one thread, in which case the caller runs the serial path.
*/
static int
packed_mul_threaded(uint8_t * C, slong Cstride,
                    const uint8_t * A, slong Astride,
                    const uint8_t * B, slong Bstride,
                    slong m, slong k, slong n, ulong nn, int ulsrc)
{
    slong nthreads = flint_get_num_threads();
    packed_par_ctx ctx;
    packed_par_job * jobs;
    thread_pool_handle * handles;
    slong planes = (nn == 2) ? 1 : 2;
    ulong * mem;
    slong words, nw, nt, tid, ntasks, phase;

    packed_tile_grid(&ctx.rp, &ctx.cp, m, n, nthreads, nn);
    ntasks = ctx.rp * ctx.cp;

    if (ntasks <= 1 || nthreads <= 1)
        return 0;

    ctx.A = A; ctx.Astride = Astride;
    ctx.B = B; ctx.Bstride = Bstride;
    ctx.C = C; ctx.Cstride = Cstride;
    ctx.m = m; ctx.k = k; ctx.n = n;
    ctx.nn = nn;
    ctx.ulsrc = ulsrc;
    ctx.aw = (k + 63) / 64;
    ctx.bw = (n + 63) / 64;
    ctx.kp = ctx.aw * 64;
    ctx.ng = (ctx.kp + GF3_T - 1) / GF3_T;

    /*
        Squaring: the same byte matrix packs to the same planes (rows
        of A and rows of B are the same rows, aw == bw when k == n), so
        allocate the A planes with the padded row count, pack once, and
        point B at them.
    */
    ctx.share = (A == B && Astride == Bstride && m == k && k == n);

    if (ctx.share)
        words = planes * (ctx.kp * ctx.aw + m * ctx.bw);
    else
        words = planes * (m * ctx.aw + ctx.kp * ctx.bw + m * ctx.bw);

    {
        TMP_INIT;

        TMP_START;
        mem = (ulong *) TMP_ALLOC(words * sizeof(ulong)
                  + ((nn == 3) ? m * ctx.ng * (slong) sizeof(uint16_t)
                               : 0));

        if (ctx.share)
        {
            ctx.a1 = mem;
            ctx.a2 = (nn == 2) ? NULL : ctx.a1 + ctx.kp * ctx.aw;
            ctx.b1 = ctx.a1;
            ctx.b2 = ctx.a2;
            ctx.c1 = ctx.a1 + planes * ctx.kp * ctx.aw;
            ctx.c2 = (nn == 2) ? NULL : ctx.c1 + m * ctx.bw;
        }
        else
        {
            ctx.a1 = mem;
            ctx.a2 = (nn == 2) ? NULL : ctx.a1 + m * ctx.aw;
            ctx.b1 = ctx.a1 + planes * m * ctx.aw;
            ctx.b2 = (nn == 2) ? NULL : ctx.b1 + ctx.kp * ctx.bw;
            ctx.c1 = ctx.b1 + planes * ctx.kp * ctx.bw;
            ctx.c2 = (nn == 2) ? NULL : ctx.c1 + m * ctx.bw;
        }
        ctx.aidx = (nn == 3) ? (uint16_t *) (mem + words) : NULL;

        nw = flint_request_threads(&handles, FLINT_MIN(ntasks, nthreads));
        nt = nw + 1;
        jobs = (packed_par_job *) TMP_ALLOC(nt * sizeof(packed_par_job));

        for (phase = 0; phase < 3; phase++)
        {
            slong total = (phase == 0) ? (ctx.share ? ctx.kp : m + ctx.kp)
                        : (phase == 1) ? ntasks : m;
#ifdef PACKED_PAR_PROF
            double ph_t0 = packed_par_prof_now();
#endif

            for (tid = 0; tid < nt; tid++)
            {
                jobs[tid].ctx = &ctx;
                jobs[tid].phase = (int) phase;
                jobs[tid].start = (total * tid) / nt;
                jobs[tid].end = (total * (tid + 1)) / nt;
            }

            for (tid = 0; tid < nw; tid++)
                thread_pool_wake(global_thread_pool, handles[tid], 0,
                                 packed_par_worker, jobs + tid);

            packed_par_worker(jobs + nw);

            for (tid = 0; tid < nw; tid++)
                thread_pool_wait(global_thread_pool, handles[tid]);
#ifdef PACKED_PAR_PROF
            fprintf(stderr, "phase %ld: %.1f ms\n", (long) phase,
                    (packed_par_prof_now() - ph_t0) * 1e3);
#endif
        }

        flint_give_back_threads(handles, nw);

        TMP_END;
    }

    return 1;
}

#endif  /* FLINT_BITS == 64 */

/* ------------------------------------------------------------------ */
/* Strassen                                                           */
/* ------------------------------------------------------------------ */

/*
    C = A +- B with zero extension: rows and columns beyond the operand
    dimensions read as zero. The overlap is a tight loop the compiler
    vectorizes; only the ragged edge, at most one row and one column at
    each recursion level, takes the general path. Keeping the bounds
    tests out of the inner loop matters: these O(n^2) passes are run
    eighteen times per Strassen level.
*/
static void
u8_add(uint8_t * C, slong Cstride, const uint8_t * A, slong Astride,
       const uint8_t * B, slong Bstride, slong ar, slong ac,
       slong br, slong bc, slong m, slong n, ulong nn)
{
    slong i, j, w = FLINT_MIN(FLINT_MIN(ac, bc), n);
    uint8_t nnb = (uint8_t) nn;

    for (i = 0; i < m; i++)
    {
        uint8_t * c = C + i * Cstride;

        if (i < ar && i < br)
        {
            const uint8_t * a = A + i * Astride;
            const uint8_t * b = B + i * Bstride;

            /*
                Branchless when it is safe to be: if 2*(n-1) fits in a
                byte then the sum cannot wrap, so s - n wraps to
                something large exactly when s < n and the unsigned
                minimum picks the right one. Above n = 128 the sum
                itself wraps and the trick silently returns garbage --
                which is what it did when the modulus cap was raised
                from 16 to 255, producing errors a constant n apart.
                These passes run eighteen times per Strassen level, so
                the fast form is worth keeping where it applies.
            */
            if (nn <= 128)
            {
                for (j = 0; j < w; j++)
                {
                    uint8_t t = (uint8_t) (a[j] + b[j]);
                    uint8_t u = (uint8_t) (t - nnb);

                    c[j] = FLINT_MIN(t, u);
                }
            }
            else
            {
                /*
                    Above n = 128 the sum can wrap a byte, but both
                    outcomes of a + b mod n fit one: a - (n - b) when
                    a >= n - b, and the untaken a + b is then < n. The
                    byte-typed select vectorizes at full width
                    (vpcmpeqb/vpblendvb chain); the earlier widening
                    through unsigned int went through dword packs at a
                    quarter of the throughput.
                */
                for (j = 0; j < w; j++)
                {
                    uint8_t d = (uint8_t) (nnb - b[j]);

                    c[j] = (a[j] >= d) ? (uint8_t) (a[j] - d)
                                       : (uint8_t) (a[j] + b[j]);
                }
            }

            for (j = w; j < n; j++)
            {
                ulong t = ((j < ac) ? a[j] : 0) + ((j < bc) ? b[j] : 0);
                c[j] = (uint8_t) (t >= nn ? t - nn : t);
            }
        }
        else if (i < ar)
        {
            for (j = 0; j < n; j++)
                c[j] = (j < ac) ? A[i * Astride + j] : 0;
        }
        else if (i < br)
        {
            for (j = 0; j < n; j++)
                c[j] = (j < bc) ? B[i * Bstride + j] : 0;
        }
        else
            memset(c, 0, n);
    }
}

static void
u8_sub(uint8_t * C, slong Cstride, const uint8_t * A, slong Astride,
       const uint8_t * B, slong Bstride, slong ar, slong ac,
       slong br, slong bc, slong m, slong n, ulong nn)
{
    slong i, j, w = FLINT_MIN(FLINT_MIN(ac, bc), n);
    uint8_t nnb = (uint8_t) nn;

    for (i = 0; i < m; i++)
    {
        uint8_t * c = C + i * Cstride;

        if (i < ar && i < br)
        {
            const uint8_t * a = A + i * Astride;
            const uint8_t * b = B + i * Bstride;

            if (nn <= 128)
            {
                for (j = 0; j < w; j++)
                {
                    uint8_t t = (uint8_t) (a[j] - b[j]);
                    uint8_t u = (uint8_t) (t + nnb);

                    c[j] = FLINT_MIN(t, u);
                }
            }
            else
            {
                /* both outcomes fit a byte: a - b, or a + (n - b) < n */
                for (j = 0; j < w; j++)
                {
                    uint8_t d = (uint8_t) (nnb - b[j]);

                    c[j] = (a[j] >= b[j]) ? (uint8_t) (a[j] - b[j])
                                          : (uint8_t) (a[j] + d);
                }
            }

            for (j = w; j < n; j++)
            {
                ulong x = (j < ac) ? A[i * Astride + j] : 0;
                ulong y = (j < bc) ? B[i * Bstride + j] : 0;
                c[j] = (uint8_t) (x >= y ? x - y : x + nn - y);
            }
        }
        else if (i < ar)
        {
            for (j = 0; j < n; j++)
                c[j] = (j < ac) ? A[i * Astride + j] : 0;
        }
        else if (i < br)
        {
            for (j = 0; j < n; j++)
            {
                ulong y = (j < bc) ? B[i * Bstride + j] : 0;
                c[j] = (uint8_t) (y == 0 ? 0 : nn - y);
            }
        }
        else
            memset(c, 0, n);
    }
}

/*
    Bytes of temporary space the recursion needs. The seven products all
    have dimensions bounded by those of the first quadrant, and the
    recursive calls are sequential, so one buffer of this size serves the
    whole tree: each level takes its two temporaries off the front and
    hands the rest to its children.
*/
static slong
u8_scratch_needed(slong m, slong k, slong n, slong cutoff)
{
    slong anr, anc, bnc, x1w;

    if (m < 2 * cutoff || k < 2 * cutoff || n < 2 * cutoff)
        return 0;

    anr = (m + 1) / 2;
    anc = (k + 1) / 2;
    bnc = (n + 1) / 2;
    x1w = FLINT_MAX(anc, bnc);

    return anr * x1w + anc * bnc
        + u8_scratch_needed(anr, anc, bnc, cutoff);
}

static void
u8_mul_recursive(uint8_t * C, slong Cstride,
                 const uint8_t * A, slong Astride,
                 const uint8_t * B, slong Bstride,
                 slong m, slong k, slong n, ulong nn, slong cutoff,
                 uint8_t * scratch)
{
    slong anr, anc, bnr, bnc, x1w;
    uint8_t * X1;
    uint8_t * X2;
    uint8_t * child;
    const uint8_t *A11, *A12, *A21, *A22, *B11, *B12, *B21, *B22;
    uint8_t *C11, *C12, *C21, *C22;

    if (m < 2 * cutoff || k < 2 * cutoff || n < 2 * cutoff)
    {
        if (nn > NMOD_MAT_U8_MAX_SHUFFLE_MOD)
        {
            u8_mul_basecase_sgemm(C, Cstride, A, Astride, B, Bstride,
                                  m, k, n, nn);
            return;
        }

        u8_mul_basecase(C, Cstride, A, Astride, B, Bstride, m, k, n, nn);
        return;
    }

    anr = (m + 1) / 2; anc = (k + 1) / 2;
    bnr = anc;         bnc = (n + 1) / 2;

    A11 = A;                              A12 = A + anc;
    A21 = A + anr * Astride;              A22 = A21 + anc;
    B11 = B;                              B12 = B + bnc;
    B21 = B + bnr * Bstride;              B22 = B21 + bnc;
    C11 = C;                              C12 = C + bnc;
    C21 = C + anr * Cstride;              C22 = C21 + bnc;

    x1w = FLINT_MAX(anc, bnc);
    X1 = scratch;
    X2 = scratch + anr * x1w;
    child = X2 + bnr * bnc;

    /* dimensions of the quadrants, which may be short by one */
    {
        slong a1r = anr, a2r = m - anr, a1c = anc, a2c = k - anc;
        slong b1r = bnr, b2r = k - bnr, b1c = bnc, b2c = n - bnc;

        (void) b2r; (void) b1c;
        slong c1r = anr, c2r = m - anr, c1c = bnc, c2c = n - bnc;

        u8_add(X1, x1w, A22, Astride, A12, Astride, a2r, a2c, a1r, a1c,
               anr, anc, nn);
        u8_add(X2, bnc, B22, Bstride, B12, Bstride, b2r, b2c, b1r, b1c,
               bnr, bnc, nn);
        u8_mul_recursive(C21, Cstride, X1, x1w, X2, bnc, c2r, anc, c1c,
                         nn, cutoff, child);

        u8_sub(X1, x1w, A22, Astride, A21, Astride, a2r, a2c, a2r, a1c,
               anr, anc, nn);
        u8_sub(X2, bnc, B22, Bstride, B21, Bstride, b2r, b2c, b2r, b1c,
               bnr, bnc, nn);
        u8_mul_recursive(C22, Cstride, X1, x1w, X2, bnc, c2r, anc, c2c,
                         nn, cutoff, child);

        u8_add(X1, x1w, X1, x1w, A12, Astride, anr, anc, a1r, a1c,
               anr, anc, nn);
        u8_add(X2, bnc, X2, bnc, B12, Bstride, bnr, bnc, b1r, b1c,
               bnr, bnc, nn);
        u8_mul_recursive(C11, Cstride, X1, x1w, X2, bnc, c1r, anc, c1c,
                         nn, cutoff, child);

        u8_sub(X1, x1w, X1, x1w, A11, Astride, anr, anc, a1r, a1c,
               anr, anc, nn);
        u8_mul_recursive(C12, Cstride, X1, x1w, B12, Bstride, c1r, anc,
                         c2c, nn, cutoff, child);

        /* A12 has k - anc columns and B21 has k - bnr rows: this is
           the one product whose inner dimension is the *second* half */
        u8_mul_recursive(X1, x1w, A12, Astride, B21, Bstride, c1r, a2c,
                         c1c, nn, cutoff, child);
        u8_add(C11, Cstride, C11, Cstride, X1, x1w, c1r, c1c, c1r, c1c,
               c1r, c1c, nn);
        u8_add(C12, Cstride, C12, Cstride, C22, Cstride, c1r, c2c, c2r, c2c,
               c1r, c2c, nn);
        u8_sub(C12, Cstride, C11, Cstride, C12, Cstride, c1r, c1c, c1r, c2c,
               c1r, c2c, nn);
        u8_sub(C11, Cstride, C21, Cstride, C11, Cstride, c2r, c1c, c1r, c1c,
               c1r, c1c, nn);
        u8_sub(X2, bnc, X2, bnc, B11, Bstride, bnr, bnc, b1r, b1c,
               bnr, bnc, nn);
        u8_mul_recursive(C21, Cstride, A21, Astride, X2, bnc, c2r, a1c,
                         c1c, nn, cutoff, child);
        u8_sub(C21, Cstride, C11, Cstride, C21, Cstride, c1r, c1c, c2r, c1c,
               c2r, c1c, nn);
        u8_add(C22, Cstride, C22, Cstride, C11, Cstride, c2r, c2c, c1r, c1c,
               c2r, c2c, nn);
        u8_mul_recursive(C11, Cstride, A11, Astride, B11, Bstride, c1r, a1c,
                         c1c, nn, cutoff, child);
        u8_add(C11, Cstride, C11, Cstride, X1, x1w, c1r, c1c, c1r, c1c,
               c1r, c1c, nn);
    }

}


/* ------------------------------------------------------------------ */
/* parallel driver                                                     */
/* ------------------------------------------------------------------ */

/*
    The output is tiled and each tile computed by an independent serial
    multiplication: tasks write disjoint parts of C and only read A and
    B, so there is no synchronization beyond the thread pool's wake and
    wait, and no atomics. Tiles are kept at least 2*cutoff on a side
    where possible, so that a task still runs Strassen internally
    instead of degenerating into a plain basecase; this is why the split
    is two-dimensional rather than a row split, which would make the
    blocks too thin to recurse as soon as there are several threads.

    The price, against a Strassen recursion parallelized internally, is
    that the recursion is applied per tile rather than globally, so the
    additions are repeated per tile and the asymptotic saving is that of
    the tile rather than of the whole matrix. In exchange there is no
    shared mutable state at all.
*/
/* the two paths want very different leaf sizes */
static slong
u8_cutoff(ulong nn)
{
    /*
        The bit-packed kernels want their own cutoff: their leaves are
        cheap per element (64 elements per XOR at modulus 2), so the
        recursion's byte-wide additions are relatively dearer than for
        the byte kernel, but the multiplication saved is dearer still.
    */
    if (nn <= 3)
        return U8_PACKED_CUTOFF;

    return (nn <= NMOD_MAT_U8_MAX_SHUFFLE_MOD) ? U8_STRASSEN_CUTOFF
                                               : U8_SGEMM_CUTOFF;
}

typedef struct
{
    uint8_t * C; slong Cstride;
    const uint8_t * A; slong Astride;
    const uint8_t * B; slong Bstride;
    slong k;
    ulong nn;
    slong i0, rows, j0, cols;
}
u8_task;

static void
u8_run_task(const u8_task * t)
{
    slong need = u8_scratch_needed(t->rows, t->k, t->cols,
                                   u8_cutoff(t->nn));
    uint8_t * scratch = NULL;
    TMP_INIT;

    TMP_START;

    if (need > 0)
        scratch = (uint8_t *) TMP_ALLOC(need);

    u8_mul_recursive(t->C + t->i0 * t->Cstride + t->j0, t->Cstride,
                     t->A + t->i0 * t->Astride, t->Astride,
                     t->B + t->j0, t->Bstride,
                     t->rows, t->k, t->cols, t->nn,
                     u8_cutoff(t->nn), scratch);

    TMP_END;
}

/* each participant takes a contiguous range of tiles, so that the work
   is shared even when there are more tiles than threads */
typedef struct
{
    u8_task * tasks;
    slong start;
    slong end;
}
u8_task_range;

static void
u8_worker(void * varg)
{
    u8_task_range * g = (u8_task_range *) varg;
    slong i;

    for (i = g->start; i < g->end; i++)
        u8_run_task(g->tasks + i);
}

/*
    Choose the tile grid.

    Two competing effects. Tiles should be square, because a task
    recurses log2(min(tile)/cutoff) Strassen levels and cutting one
    dimension into T pieces costs log2(T) levels where cutting both into
    sqrt(T) costs half that. But the tile count must divide the work
    evenly among the threads, because every participant takes a whole
    number of tiles: nine tiles on eight threads means one thread does
    two and the whole multiplication waits for it. Balance wins -- an
    unbalanced grid caps the speedup at (tiles per thread) regardless of
    how good the tiles are.

    So: search grids with at most nthreads tiles, take the one whose
    tiles have the largest minimum dimension, and prefer more tiles on a
    tie. With 8 threads this gives 2x4 rather than 3x3 (9 tiles, one
    thread doing two) or 8x1 (deep enough grid but one Strassen level).
*/
static void
u8_tile_grid(slong * rp, slong * cp, slong m, slong n, slong nthreads,
             slong min_tile)
{
    slong rmax = FLINT_MAX(1, m / min_tile);
    slong cmax = FLINT_MAX(1, n / min_tile);
    slong r, bestr = 1, bestc = 1, bestscore = 0, bestn = 1;

    rmax = FLINT_MIN(rmax, nthreads);

    for (r = 1; r <= rmax; r++)
    {
        slong c = FLINT_MIN(cmax, nthreads / r);
        slong score, tiles;

        if (c < 1)
            continue;

        score = FLINT_MIN(m / r, n / c);
        tiles = r * c;

        /* utilization first, shape second: an idle thread costs more
           than a tile with one Strassen level fewer */
        if (tiles > bestn || (tiles == bestn && score > bestscore))
        {
            bestscore = score;
            bestn = tiles;
            bestr = r;
            bestc = c;
        }
    }

    *rp = bestr;
    *cp = bestc;
}

/* ------------------------------------------------------------------ */
/* public interface                                                   */
/* ------------------------------------------------------------------ */


void
_nmod_mat_mul_u8(uint8_t * C, slong Cstride,
                 const uint8_t * A, slong Astride,
                 const uint8_t * B, slong Bstride,
                 slong m, slong k, slong n, nmod_t mod)
{
    slong i;

    FLINT_ASSERT(mod.n <= 255);

    if (m <= 0 || n <= 0)
        return;

    if (k <= 0)
    {
        for (i = 0; i < m; i++)
            memset(C + i * Cstride, 0, n);
        return;
    }

    /* Some of the algorithms below support aliasing, but to keep things
       simple let's just make a copy. */
    if (C == A || C == B)
    {
        uint8_t * Ct;
        TMP_INIT;

        TMP_START;
        Ct = (uint8_t *) TMP_ALLOC(m * n);
        _nmod_mat_mul_u8(Ct, n, A, Astride, B, Bstride, m, k, n, mod);

        for (i = 0; i < m; i++)
            memcpy(C + i * Cstride, Ct + i * n, n);

        TMP_END;
        return;
    }

    /*
        Moduli 2 and 3: the threaded path shares the packed operands
        across a tile grid (see packed_mul_threaded); when the grid
        degenerates or there is one thread it declines and the serial
        packed Strassen runs. With the leaf cutoffs up at 4096 and
        8192 a classical split of C parallelizes better than the
        seven-task Strassen tree would. 32-bit builds skip the packed
        engines and take the generic byte path below.
    */
#if FLINT_BITS == 64
    if (mod.n == 2 || mod.n == 3)
    {
        if (m < ((mod.n == 2) ? NMOD_MAT_U8_GF2_BASECASE_M
                              : NMOD_MAT_U8_GF3_BASECASE_M))
        {
            u8_mul_basecase(C, Cstride, A, Astride, B, Bstride,
                            m, k, n, mod.n);
            return;
        }

        if (packed_mul_threaded(C, Cstride, A, Astride, B, Bstride,
                                m, k, n, mod.n, 0))
            return;

        if (mod.n == 2)
            gf2_mul_strassen(C, Cstride, A, Astride, B, Bstride, m, k, n,
                             packed_strassen_levels(m, k, n,
                                 nmod_mat_gf2_strassen_cutoff), 0);
        else
            gf3_mul_strassen(C, Cstride, A, Astride, B, Bstride, m, k, n,
                             packed_strassen_levels(m, k, n,
                                 nmod_mat_gf3_strassen_cutoff), 0);

        return;
    }
#endif  /* FLINT_BITS == 64 */

    if (mod.n > NMOD_MAT_U8_MAX_SHUFFLE_MOD
        && u8_mul_flat_sgemm(C, Cstride, A, Astride, B, Bstride,
                             m, k, n, mod.n))
        return;

    {
        slong nthreads = flint_get_num_threads();
        slong rp, cp, ntasks;

        u8_tile_grid(&rp, &cp, m, n, nthreads, U8_PAR_MIN_TILE);
        ntasks = rp * cp;

        if (ntasks <= 1 || nthreads <= 1)
        {
            slong need = u8_scratch_needed(m, k, n, u8_cutoff(mod.n));
            uint8_t * scratch = NULL;
            TMP_INIT;

            TMP_START;

            if (need > 0)
                scratch = (uint8_t *) TMP_ALLOC(need);

            u8_mul_recursive(C, Cstride, A, Astride, B, Bstride,
                             m, k, n, mod.n, u8_cutoff(mod.n),
                             scratch);

            TMP_END;
        }
        else
        {
            thread_pool_handle * handles;
            u8_task * tasks;
            u8_task_range * groups;
            slong nw, nt, tid, ri, ci, pos;
            TMP_INIT;

            TMP_START;
            tasks = (u8_task *) TMP_ALLOC(ntasks * sizeof(u8_task));
            groups = (u8_task_range *) TMP_ALLOC((nthreads + 1)
                                                 * sizeof(u8_task_range));

            pos = 0;

            for (ri = 0; ri < rp; ri++)
                for (ci = 0; ci < cp; ci++)
                {
                    u8_task * t = tasks + pos++;

                    t->C = C; t->Cstride = Cstride;
                    t->A = A; t->Astride = Astride;
                    t->B = B; t->Bstride = Bstride;
                    t->k = k; t->nn = mod.n;
                    t->i0 = (m * ri) / rp;
                    t->rows = (m * (ri + 1)) / rp - t->i0;
                    t->j0 = (n * ci) / cp;
                    t->cols = (n * (ci + 1)) / cp - t->j0;
                }

            nw = flint_request_threads(&handles, FLINT_MIN(ntasks, nthreads));
            nt = nw + 1;

            for (tid = 0; tid < nt; tid++)
            {
                groups[tid].tasks = tasks;
                groups[tid].start = (ntasks * tid) / nt;
                groups[tid].end = (ntasks * (tid + 1)) / nt;
            }

            for (tid = 0; tid < nw; tid++)
                thread_pool_wake(global_thread_pool, handles[tid], 0,
                                 u8_worker, groups + tid);

            u8_worker(groups + nw);

            for (tid = 0; tid < nw; tid++)
                thread_pool_wait(global_thread_pool, handles[tid]);

            flint_give_back_threads(handles, nw);

            TMP_END;
        }
    }

    return;
}

/*
    The ulong <-> uint8 roundtrip is O(mk + kn + mn) and was the only
    serial part of nmod_mat_mul_u8: at n = 8192 it is a fifth of the
    time at eight threads, and by Amdahl it caps the speedup of the
    whole routine well below that of the multiplication itself. It is
    trivially data parallel -- disjoint row ranges, no shared state --
    so it is split the same way as the multiplication, with the same
    absence of synchronization.
*/
typedef struct
{
    const nmod_mat_struct * M;
    uint8_t * U;
    slong stride;
    slong r0, r1;
    int to_bytes;
}
u8_conv_task;

/*
    Row conversions between the ulong entries and the byte image. The
    compilers vectorize the plain loops at 128 bits only (the 64 -> 8
    narrowing does not fit their cost models), which left the wrapper's
    conversion overhead around 1 ns per entry; the explicit kernels
    below are within a factor of the pure memory traffic. Entries are
    assumed reduced below 256, as everywhere in this file.
*/
static void
u8_row_to_bytes(uint8_t * dst, const ulong * src, slong c)
{
    slong j = 0;

#if defined(NMOD_MAT_U8_GF2_AVX512)
    for (; j + 8 <= c; j += 8)
    {
        __m512i v = _mm512_loadu_si512((const void *) (src + j));
        _mm_storel_epi64((__m128i *) (dst + j), _mm512_cvtepi64_epi8(v));
    }
#elif defined(NMOD_MAT_U8_X86)
    /*
        16 qwords -> 16 bytes: select the low dword of every qword with
        two cross-lane vpermd, then the usual packusdw/packuswb pair
        with a final vpermq to undo their in-lane interleaving.
    */
    {
        const __m256i sel = _mm256_setr_epi32(0, 2, 4, 6, 1, 3, 5, 7);

        for (; j + 16 <= c; j += 16)
        {
            __m256i a = _mm256_loadu_si256((const __m256i *) (src + j));
            __m256i b = _mm256_loadu_si256((const __m256i *) (src + j + 4));
            __m256i cc = _mm256_loadu_si256((const __m256i *) (src + j + 8));
            __m256i d = _mm256_loadu_si256((const __m256i *) (src + j + 12));
            __m256i ab, cd, w, r;

            a = _mm256_permutevar8x32_epi32(a, sel);
            b = _mm256_permutevar8x32_epi32(b, sel);
            cc = _mm256_permutevar8x32_epi32(cc, sel);
            d = _mm256_permutevar8x32_epi32(d, sel);

            /* low dwords of 8 qwords in each */
            ab = _mm256_inserti128_si256(a, _mm256_castsi256_si128(b), 1);
            cd = _mm256_inserti128_si256(cc, _mm256_castsi256_si128(d), 1);

            w = _mm256_packus_epi32(ab, cd);
            w = _mm256_permute4x64_epi64(w, 0xD8);
            r = _mm256_packus_epi16(w, w);
            r = _mm256_permute4x64_epi64(r, 0xD8);

            _mm_storeu_si128((__m128i *) (dst + j),
                             _mm256_castsi256_si128(r));
        }
    }
#endif

    for (; j < c; j++)
        dst[j] = (uint8_t) src[j];
}

static void
u8_row_from_bytes(ulong * dst, const uint8_t * src, slong c)
{
    slong j = 0;

#if defined(NMOD_MAT_U8_GF2_AVX512)
    for (; j + 8 <= c; j += 8)
    {
        __m128i v = _mm_loadl_epi64((const __m128i *) (src + j));
        _mm512_storeu_si512((void *) (dst + j), _mm512_cvtepu8_epi64(v));
    }
#elif defined(NMOD_MAT_U8_X86)
    for (; j + 4 <= c; j += 4)
    {
        int32_t vi;
        __m128i v;
        memcpy(&vi, src + j, 4);
        v = _mm_cvtsi32_si128(vi);
        _mm256_storeu_si256((__m256i *) (dst + j), _mm256_cvtepu8_epi64(v));
    }
#endif

    for (; j < c; j++)
        dst[j] = src[j];
}

static void
u8_conv_worker(void * varg)
{
    u8_conv_task * t = (u8_conv_task *) varg;
    slong i, c = t->M->c;

    if (t->to_bytes)
    {
        for (i = t->r0; i < t->r1; i++)
            u8_row_to_bytes(t->U + i * t->stride,
                            t->M->entries + i * t->M->stride, c);
    }
    else
    {
        for (i = t->r0; i < t->r1; i++)
            u8_row_from_bytes(t->M->entries + i * t->M->stride,
                              t->U + i * t->stride, c);
    }
}

/* convert one matrix, in parallel when it is large enough to pay for
   the hand-off */
static void
u8_convert(const nmod_mat_t M, uint8_t * U, slong stride, int to_bytes,
           slong nthreads)
{
    slong rows = M->r;
    slong nw = 0, nt, i;
    thread_pool_handle * handles = NULL;
    u8_conv_task * tasks;
    TMP_INIT;

    if (nthreads > 1 && rows * M->c >= 4000000)
        nw = flint_request_threads(&handles, FLINT_MIN(nthreads, rows));

    nt = nw + 1;

    TMP_START;
    tasks = (u8_conv_task *) TMP_ALLOC(nt * sizeof(u8_conv_task));

    for (i = 0; i < nt; i++)
    {
        tasks[i].M = M;
        tasks[i].U = U;
        tasks[i].stride = stride;
        tasks[i].to_bytes = to_bytes;
        tasks[i].r0 = (rows * i) / nt;
        tasks[i].r1 = (rows * (i + 1)) / nt;
    }

    for (i = 0; i < nw; i++)
        thread_pool_wake(global_thread_pool, handles[i], 0,
                         u8_conv_worker, tasks + i);

    u8_conv_worker(tasks + nw);

    for (i = 0; i < nw; i++)
        thread_pool_wait(global_thread_pool, handles[i]);

    if (nw != 0)
        flint_give_back_threads(handles, nw);

    TMP_END;
}

/*
    For nmod_mat the entries are words, so reaching the byte kernels
    costs a ulong -> uint8 conversion of both operands and a conversion
    back. That is worth paying when the byte shuffle kernel applies,
    which is roughly a 1.4x win at moduli up to 16. It is not worth
    paying above that, where the leaf is computed in single precision
    anyway: the route would be ulong -> uint8 -> float, whereas
    nmod_mat_mul_blas converts ulong -> float once. So this reports
    that it does not apply, and callers keep their existing path.

    gr_mat over nmod8 is the opposite case: its entries are already
    bytes, so `_nmod_mat_mul_u8` costs no conversion at all and is used
    over the whole supported range.
*/
void
nmod_mat_mul_u8(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    slong m = A->r, k = A->c, n = B->c;
    uint8_t * Au;
    uint8_t * Bu;
    uint8_t * Cu;
    TMP_INIT;

    /*
        Every modulus up to 255 is handled; the caller is responsible
        for choosing between this and nmod_mat_mul_blas, informed by
        profile/p-mul_u8 (an earlier guard against moduli above 15,
        where blas had once measured faster, was removed after the
        sgemm-backed basecase caught up).
    */
    FLINT_ASSERT(C->mod.n <= 255);

    if (m <= 0 || n <= 0)
        return;

#if FLINT_BITS == 64
    /*
        Moduli 2 and 3, above the tiny-size takeover: pack the bits
        straight from the ulong entries. The multiplication itself is
        nearly free next to the memory traffic here, and the byte
        image the general route builds -- written once and read once
        purely as a relay into the pack -- was measured at nearly
        thrice the cost of the whole raw multiplication at 4096.
        Squaring is detected inside as usual.
    */
    if ((C->mod.n == 2 || C->mod.n == 3)
        && m >= ((C->mod.n == 2) ? NMOD_MAT_U8_GF2_BASECASE_M
                                 : NMOD_MAT_U8_GF3_BASECASE_M))
    {
        const uint8_t * Aul = (const uint8_t *) A->entries;
        const uint8_t * Bul = (const uint8_t *) B->entries;
        uint8_t * Cul = (uint8_t *) C->entries;

        if (!packed_mul_threaded(Cul, C->stride, Aul, A->stride,
                                 Bul, B->stride, m, k, n, C->mod.n, 1))
        {
            if (C->mod.n == 2)
                gf2_mul_strassen(Cul, C->stride, Aul, A->stride,
                                 Bul, B->stride, m, k, n,
                                 packed_strassen_levels(m, k, n,
                                     nmod_mat_gf2_strassen_cutoff), 1);
            else
                gf3_mul_strassen(Cul, C->stride, Aul, A->stride,
                                 Bul, B->stride, m, k, n,
                                 packed_strassen_levels(m, k, n,
                                     nmod_mat_gf3_strassen_cutoff), 1);
        }

        return;
    }
#endif

    /* alloca for small matrices: the roundtrip must not cost an
       allocation when the multiplication itself is a few microseconds */
    TMP_START;
    Au = (uint8_t *) TMP_ALLOC(m * k + ((A == B) ? 0 : k * n) + m * n);
    Bu = (A == B) ? Au : Au + (size_t) m * k;
    Cu = Bu + (size_t) k * n;

    {
        slong nthreads = flint_get_num_threads();

        /* squaring: one operand image, one conversion pass */
        u8_convert(A, Au, k, 1, nthreads);
        if (A != B)
            u8_convert(B, Bu, n, 1, nthreads);

        _nmod_mat_mul_u8(Cu, n, Au, k, Bu, n, m, k, n, C->mod);
        u8_convert(C, Cu, n, 0, nthreads);
    }

    TMP_END;
}
