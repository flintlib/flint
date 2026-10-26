/*
    Copyright (C) 2021 Daniel Schultz
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
    c = a * B for a few vectors a at once (r = 1 for nmod_mat_nmod_vec_mul,
    up to r <= NMOD_MAT_MUL_ROWS_MAX for the products with few rows as
    dispatched by nmod_mat_mul)

    This does a rank-1 accumulation row by row, c[s][j] = sum_i a[s][i]
    B[i][j]. Each row of B is read once and each vector of it loaded once for
    the r rows, in blocks of columns whose accumulators stay in registers, or
    in L1 for very wide B; the reductions are delayed as in the SIMD matrix
    kernels:

    - with AVX512-IFMA and n <= 2^52, the products are accumulated as
      their low and high 52-bit halves (nmod_mat_mul_u52), the pairs of
      accumulators being reduced to a residue in registers every 511 rows
      (mul_u52_vec.h); for n <= 2^26 the products fit in 52 bits and the
      low half suffices (one IFMA per vector instead of two);
    - without IFMA but with AVX2 or AVX-512 and n <= 2^32, the products
      (below 2^64, vpmuludq) are split at bit 56 into two accumulators
      per column as in the split dot products of nmod_vec/dot.c, the low
      one folded into the high one every 255 rows;
    - otherwise with a vector backend of mul_fp_vec.h and n < 2^50, in
      double precision with the mulmod of fft_small (nmod_mat_mul_fp50),
      the accumulators being brought back to (-0.51 n, 0.51 n) every four
      rows;
    - otherwise with AVX2 or AVX-512 (n >= 2^50, or n > 2^52 with IFMA),
      the vpmuludq products of limbs of at most 32 bits of the entries are
      accumulated in three to six lanes per column as in the split limbs
      dot products of nmod_vec/dot_split_limbs.c, carried into each other
      every few rows (n <= 2^61) or never (32-bit limbs), and the totals
      (up to three words) are reduced once at the end; for B larger than
      L2 the rows go by chunks of 16, the accumulators being kept in L2
      in between;
    - otherwise, and for nmod_mat_nmod_vec_mul with B too narrow for
      the SIMD path at hand (vecmat_use_rows_simd), one full modular
      multiplication per entry for short vectors, the classical product
      (nmod_mat_mul_classical) for longer ones.
*/

#define VECMAT_ROWS NMOD_MAT_MUL_ROWS_MAX

/*
    Register blocks read 64-512 bytes per row of B and per block; up to
    VECMAT_REGS_MAX_COLS columns they are used, for chunks of rows of about
    VECMAT_L2_BYTES which the prefetching triggered by the first block leaves
    in L2 for the others (with the narrow blocks of many rows this
    prefetching stops working beyond VECMAT_REGS_WIDE_COLS columns, where the
    rows are then taken 3 at a time); beyond, blocks of columns with
    accumulators in L1 stream the rows of B in runs of a few KB.
    FIXME measured on one machine only: need to tune these knobs for
    other machines?
*/
#ifndef VECMAT_L2_BYTES
# define VECMAT_L2_BYTES (1 << 20)
#endif
#ifndef VECMAT_REGS_MAX_COLS
# define VECMAT_REGS_MAX_COLS 4096
#endif
#ifndef VECMAT_REGS_WIDE_COLS
# define VECMAT_REGS_WIDE_COLS 2048
#endif

#if FLINT_BITS == 64 && NMOD_MAT_HAVE_MUL_U52

#include "nmod_mat/mul_u52_vec.h"

/*
    Rows per chunk: a lane holds a residue plus the halves of up to
    VECMAT_U52_CHUNK products, within the bounds of u52_red2 (mul_u52_vec.h).
*/
#define VECMAT_U52_CHUNK 511

/* the per-chunk reduction of the accumulators of a block; not inlined,
   to keep the many instantiations of vecmat_u52_block small */
FLINT_STATIC_NOINLINE void
vecmat_u52_reduce(__m512i * lo, __m512i * hi, slong cnt, int withhi,
                  const u52_consts * C)
{
    slong v;

    if (withhi)
        for (v = 0; v < cnt; v++)
        {
            lo[v] = u52_red2(lo[v], hi[v], C);
            hi[v] = _mm512_setzero_si512();
        }
    else
        for (v = 0; v < cnt; v++)
            lo[v] = u52_red(lo[v], C);
}

/*
    One block of nv vectors of columns (the last one masked by m) for r
    rows and the rows i0 <= i < i1 of B: the accumulators lo[s * nv + v],
    hi[s * nv + v] start from the residues held in c (zero at first), and
    every chunk of rows the pairs are reduced to a residue kept in lo for
    the next chunk; the last one is stored into c.

    r, withhi and regs are constants at each call site. With regs = 1, nv
    is a constant too (<= 8) and the arrays lo, hi, local to the caller,
    are register allocated: they go through a buffer for the reduction.
    Otherwise they are in L1 and nv is arbitrary.
*/
FLINT_FORCE_INLINE void
vecmat_u52_block(ulong * const * c, __m512i * lo, __m512i * hi,
                 const ulong * const * a, const nmod_mat_t B, slong j0,
                 slong i0, slong i1, int r, slong nv, int withhi, int regs,
                 __mmask8 m, const u52_consts * C)
{
    __m512i av[VECMAT_ROWS], bv;
    slong i, i2, v;
    int s;

    for (s = 0; s < r; s++)
    {
        for (v = 0; v < nv - 1; v++)
            lo[s * nv + v] = _mm512_loadu_si512((const void *) (c[s] + j0 + 8 * v));
        lo[s * nv + v] = _mm512_maskz_loadu_epi64(m, (const void *) (c[s] + j0 + 8 * v));
    }
    for (v = 0; v < r * nv; v++)
        hi[v] = _mm512_setzero_si512();

    for (i = i0; i < i1; i = i2)
    {
        i2 = FLINT_MIN(i1, i + VECMAT_U52_CHUNK);

        for ( ; i < i2; i++)
        {
            const ulong * b = nmod_mat_entry_ptr(B, i, j0);

            for (s = 0; s < r; s++)
                av[s] = _mm512_set1_epi64(a[s][i]);

            for (v = 0; v < nv - 1; v++)
            {
                bv = _mm512_loadu_si512((const void *) (b + 8 * v));
                for (s = 0; s < r; s++)
                {
                    lo[s * nv + v] = u52_madd52lo(lo[s * nv + v], av[s], bv);
                    if (withhi)
                        hi[s * nv + v] = u52_madd52hi(hi[s * nv + v], av[s], bv);
                }
            }
            bv = _mm512_maskz_loadu_epi64(m, (const void *) (b + 8 * v));
            for (s = 0; s < r; s++)
            {
                lo[s * nv + v] = u52_madd52lo(lo[s * nv + v], av[s], bv);
                if (withhi)
                    hi[s * nv + v] = u52_madd52hi(hi[s * nv + v], av[s], bv);
            }
        }

        if (regs)
        {
            __m512i buf[2 * VECMAT_ROWS * 8];

            for (v = 0; v < r * nv; v++)
            {
                buf[v] = lo[v];
                if (withhi)
                    buf[r * nv + v] = hi[v];
            }
            vecmat_u52_reduce(buf, buf + r * nv, r * nv, withhi, C);
            for (v = 0; v < r * nv; v++)
            {
                lo[v] = buf[v];
                if (withhi)
                    hi[v] = _mm512_setzero_si512();
            }
        }
        else
            vecmat_u52_reduce(lo, hi, r * nv, withhi, C);
    }

    for (s = 0; s < r; s++)
    {
        for (v = 0; v < nv - 1; v++)
            _mm512_storeu_si512((void *) (c[s] + j0 + 8 * v), lo[s * nv + v]);
        _mm512_mask_storeu_epi64((void *) (c[s] + j0 + 8 * v), m, lo[s * nv + v]);
    }
}

/*
    Many columns: blocks of columns whose accumulators (16 bytes per
    column and row with the high halves, 8 without) stay in L1, so that
    the rows of B are streamed in runs of a few KB; up to 4 rows.
*/
#ifndef VECMAT_U52_WIDTH
# define VECMAT_U52_WIDTH 512
#endif

#define VECMAT_U52_WIDE_CASE(R) \
    case R: \
        if (withhi) \
            vecmat_u52_block(c, lo, hi, a, B, j0, 0, len, R, nv, 1, 0, m, C); \
        else \
            vecmat_u52_block(c, lo, hi, a, B, j0, 0, len, R, nv, 0, 0, m, C); \
        break;

static void
vecmat_u52_wide(ulong * const * c, const ulong * const * a, slong len,
                const nmod_mat_t B, int r, int withhi, const u52_consts * C)
{
    const slong ncols = B->c;
    /* columns per block: the accumulators of the r rows share the budget */
    const slong wmax = 8 * ((withhi ? VECMAT_U52_WIDTH : 2 * VECMAT_U52_WIDTH) / (8 * r));
    __m512i lo[2 * VECMAT_U52_WIDTH / 8], hi[2 * VECMAT_U52_WIDTH / 8];
    slong j0 = 0;
    int s;

    for (s = 0; s < r; s++)
        _nmod_vec_zero(c[s], ncols);

    while (j0 < ncols)
    {
        const slong w = FLINT_MIN(ncols - j0, wmax);
        const slong nv = (w + 7) / 8;
        const __mmask8 m = (w % 8 == 0) ? 0xFF : (__mmask8) ((1 << (w % 8)) - 1);

        switch (r)
        {
            VECMAT_U52_WIDE_CASE(1)
            VECMAT_U52_WIDE_CASE(2)
            VECMAT_U52_WIDE_CASE(3)
            VECMAT_U52_WIDE_CASE(4)
            default: FLINT_UNREACHABLE;
        }

        j0 += w;
    }
}

/*
    Blocks of 8 nv columns whose accumulators are registers, nv = 8, 6, 4,
    3, 2, 2, 1, 1 for r = 1..8 rows (2 r nv of the 32 zmm), the last,
    narrower block having its own instantiation for each nv (accumulators
    in memory would be bound by the latency of the store-load forwarding
    between rows when there are few of them). Each block reads its 64 nv
    bytes per row of B: when B is not in the caches, this is done for
    chunks of rows of at most VECMAT_L2_BYTES, which the first block brings
    into L2 for the others.
*/
static const int vecmat_u52_regs_nv[VECMAT_ROWS] = {8, 6, 4, 3, 2, 2, 1, 1};

#define VECMAT_U52_REGS_CASE(R, NV) \
    case 8 * (R) + (NV): \
        if (withhi) \
            vecmat_u52_block(c, lo, hi, a, B, j0, i0, i1, R, NV, 1, 1, m, C); \
        else \
            vecmat_u52_block(c, lo, hi, a, B, j0, i0, i1, R, NV, 0, 1, m, C); \
        break;

static void
vecmat_u52_regs(ulong * const * c, const ulong * const * a, slong len,
                const nmod_mat_t B, int r, int withhi, const u52_consts * C)
{
    const slong ncols = B->c;
    const int nvmax = vecmat_u52_regs_nv[r - 1];
    const slong rows = FLINT_MAX(64, VECMAT_L2_BYTES / (8 * ncols));
    slong i0, i1;
    int s;

    for (s = 0; s < r; s++)
        _nmod_vec_zero(c[s], ncols);

    for (i0 = 0; i0 < len; i0 = i1)
    {
        slong j0 = 0;

        i1 = FLINT_MIN(len, i0 + rows);

        while (j0 < ncols)
        {
            __m512i lo[VECMAT_ROWS * 8], hi[VECMAT_ROWS * 8];
            const slong w = FLINT_MIN(ncols - j0, 8 * nvmax);
            const int nv = (int) ((w + 7) / 8);
            const __mmask8 m = (w % 8 == 0) ? 0xFF : (__mmask8) ((1 << (w % 8)) - 1);

            switch (8 * r + nv)
            {
                VECMAT_U52_REGS_CASE(1, 1)
                VECMAT_U52_REGS_CASE(1, 2)
                VECMAT_U52_REGS_CASE(1, 3)
                VECMAT_U52_REGS_CASE(1, 4)
                VECMAT_U52_REGS_CASE(1, 5)
                VECMAT_U52_REGS_CASE(1, 6)
                VECMAT_U52_REGS_CASE(1, 7)
                VECMAT_U52_REGS_CASE(1, 8)
                VECMAT_U52_REGS_CASE(2, 1)
                VECMAT_U52_REGS_CASE(2, 2)
                VECMAT_U52_REGS_CASE(2, 3)
                VECMAT_U52_REGS_CASE(2, 4)
                VECMAT_U52_REGS_CASE(2, 5)
                VECMAT_U52_REGS_CASE(2, 6)
                VECMAT_U52_REGS_CASE(3, 1)
                VECMAT_U52_REGS_CASE(3, 2)
                VECMAT_U52_REGS_CASE(3, 3)
                VECMAT_U52_REGS_CASE(3, 4)
                VECMAT_U52_REGS_CASE(4, 1)
                VECMAT_U52_REGS_CASE(4, 2)
                VECMAT_U52_REGS_CASE(4, 3)
                VECMAT_U52_REGS_CASE(5, 1)
                VECMAT_U52_REGS_CASE(5, 2)
                VECMAT_U52_REGS_CASE(6, 1)
                VECMAT_U52_REGS_CASE(6, 2)
                VECMAT_U52_REGS_CASE(7, 1)
                VECMAT_U52_REGS_CASE(8, 1)
                default: FLINT_UNREACHABLE;
            }

            j0 += w;
        }
    }
}

#endif  /* u52 */

#if FLINT_BITS == 64 && (NMOD_MAT_HAVE_VECMAT_U32 || NMOD_MAT_HAVE_VECMAT_U64)

#include <immintrin.h>

/* 64-bit lanes for the integer paths on x86 (u32 and u64 below); the
   shift counts of vmu_srli / vmu_slli are constants, those of vmu_srlv /
   vmu_sllv vectors (a count of 64 or more gives 0) */
#if defined(__AVX512F__)
# define VMU_VL 8
typedef __m512i vmu_t;
typedef __mmask8 vmu_mask_t;
# define vmu_zero()        _mm512_setzero_si512()
# define vmu_set1(x)       _mm512_set1_epi64(x)
# define vmu_loadu(p)      _mm512_loadu_si512((const void *) (p))
# define vmu_storeu(p, x)  _mm512_storeu_si512((void *) (p), x)
# define vmu_mul(x, y)     _mm512_mul_epu32(x, y)
# define vmu_add(x, y)     _mm512_add_epi64(x, y)
# define vmu_sub(x, y)     _mm512_sub_epi64(x, y)
# define vmu_and(x, y)     _mm512_and_si512(x, y)
# define vmu_or(x, y)      _mm512_or_si512(x, y)
# define vmu_srli(x, c)    _mm512_srli_epi64(x, c)
# define vmu_slli(x, c)    _mm512_slli_epi64(x, c)
# define vmu_srlv(x, c)    _mm512_srlv_epi64(x, c)
# define vmu_sllv(x, c)    _mm512_sllv_epi64(x, c)
/* the first w lanes (0 < w <= VMU_VL), and a load of them (zero above) */
# define vmu_mask(w)       ((__mmask8) ((1 << (w)) - 1))
# define vmu_loadu_mask(p, m) _mm512_maskz_loadu_epi64(m, (const void *) (p))
/* the same lanes, from a load of the whole vector (all of it must be
   readable) */
# define vmu_loadu_clear(p, m) _mm512_maskz_mov_epi64(m, vmu_loadu(p))
#else
# define VMU_VL 4
typedef __m256i vmu_t;
typedef __m256i vmu_mask_t;
# define vmu_zero()        _mm256_setzero_si256()
# define vmu_set1(x)       _mm256_set1_epi64x(x)
# define vmu_loadu(p)      _mm256_loadu_si256((const __m256i *) (p))
# define vmu_storeu(p, x)  _mm256_storeu_si256((__m256i *) (p), x)
# define vmu_mul(x, y)     _mm256_mul_epu32(x, y)
# define vmu_add(x, y)     _mm256_add_epi64(x, y)
# define vmu_sub(x, y)     _mm256_sub_epi64(x, y)
# define vmu_and(x, y)     _mm256_and_si256(x, y)
# define vmu_or(x, y)      _mm256_or_si256(x, y)
# define vmu_srli(x, c)    _mm256_srli_epi64(x, c)
# define vmu_slli(x, c)    _mm256_slli_epi64(x, c)
# define vmu_srlv(x, c)    _mm256_srlv_epi64(x, c)
# define vmu_sllv(x, c)    _mm256_sllv_epi64(x, c)
# define vmu_mask(w)       _mm256_cmpgt_epi64(_mm256_set1_epi64x(w), \
                                              _mm256_set_epi64x(3, 2, 1, 0))
# define vmu_loadu_mask(p, m) _mm256_maskload_epi64((const long long *) (p), m)
# define vmu_loadu_clear(p, m) _mm256_and_si256(vmu_loadu(p), m)
#endif

#endif

#if FLINT_BITS == 64 && NMOD_MAT_HAVE_VECMAT_U32

/*
    Moduli up to 2^32 without IFMA: vpmuludq products below 2^64, split
    at bit 56 (nmod_vec/dot.c): lo accumulates the low 56 bits of the
    products, hi their high 8 bits, and lo is folded into hi every
    VECMAT_U32_FOLD rows. The total lo + 2^56 hi fits two limbs for any
    length. The reduction to the residue is scalar, so the chunks of rows
    between two of them are long (VECMAT_U32_ROWS at least).
*/

#if VMU_VL == 8
/* vectors of columns per register block for r rows (2 r nv + r + 3 of
   32 registers), powers of two */
static const int vecmat_u32_regs_nv[4] = {8, 4, 4, 2};
#else
/* 16 registers */
static const int vecmat_u32_regs_nv[4] = {4, 2, 1, 1};
#endif

#define VECMAT_U32_FOLD 255
#define VECMAT_U32_WIDTH 512
#define VECMAT_U32_ROWS 255
/* the register blocks are narrow (128-512 bytes per row of B): beyond
   this many columns the L1 blocks are used */
#ifndef VECMAT_U32_REGS_MAX_COLS
# define VECMAT_U32_REGS_MAX_COLS (128 * VMU_VL)
#endif

/*
    One block of nv vectors of columns j0 <= j < j0 + w for r rows and the
    rows i0 <= i < i1 of B; lo[s * nv + v] starts from the residues in c
    (zero at first) and the residues of the totals are stored back. r is a
    constant at each call site; with nv also constant and the arrays lo,
    hi local to the caller, they are register allocated (vecmat_u32_regs),
    otherwise they are in L1 and nv is arbitrary (vecmat_u32_wide).

    A last partial vector (w not a multiple of VMU_VL; the block then ends
    at the last column of B) has its lanes above w cleared. It used to be
    copied to a zero-padded buffer in each row, and the wide load of these
    narrow stores could not be forwarded: 1.4-6x slower for ncols = 1-7
    and 9-15 than with a masked load (Cascade Lake). With the masked load
    intrinsic, however, GCC keeps the accumulators lo, hi in memory (a
    store of each per row), so the rows whose whole last vector lies inside
    B (it runs into the next rows) load it whole and clear the lanes (on
    AVX-512, GCC turns this into the masked load instruction, with the
    accumulators in registers); only the last rows, at most VMU_VL - 1 of
    them, use the masked load intrinsic (with GCC and AVX-512; otherwise all
    rows use it). The rows loop vecmat_u32_rows has a copy for each case
    (mode is a constant at the call sites): with a test in the loop, the
    compiler made the load of a full last vector a masked one too, 1.2-2x
    slower. For the same reason, the initial load of the residues from c is
    masked only for a partial vector.
*/

/* the two workarounds below are for GCC with AVX-512 (with clang, or with
   AVX2, they were neutral or slightly slower) */
/* TODO to be re-checked in the future with latest GCC */
#if defined(__GNUC__) && !defined(__clang__) && VMU_VL == 8
# define VECMAT_U32_GCC 1
#else
# define VECMAT_U32_GCC 0
#endif

/* mode 0: full last vector, 1: partial, loaded whole and cleared,
   2: partial, masked load; returns the rows since the last fold */
FLINT_FORCE_INLINE int
vecmat_u32_rows(vmu_t * lo, vmu_t * hi, const ulong * const * a,
                const nmod_mat_t B, slong j0, slong i0, slong i1, int r,
                slong nv, int mode, vmu_mask_t mlast, int cnt)
{
    vmu_t mask = vmu_set1((UWORD(1) << 56) - 1);
    vmu_t av[VECMAT_ROWS], bv, p;
    slong i, v;
    int s;

#if VECMAT_U32_GCC
    /* keep the mask in a register: otherwise GCC may read it from memory
       in each row */
    __asm__ ("" : "+v" (mask));
#endif

    for (i = i0; i < i1; i++)
    {
        const ulong * b = nmod_mat_entry_ptr(B, i, j0);

        for (s = 0; s < r; s++)
            av[s] = vmu_set1(a[s][i]);

        for (v = 0; v < nv - 1; v++)
        {
            bv = vmu_loadu(b + VMU_VL * v);
            for (s = 0; s < r; s++)
            {
                p = vmu_mul(av[s], bv);
                lo[s * nv + v] = vmu_add(lo[s * nv + v], vmu_and(p, mask));
                hi[s * nv + v] = vmu_add(hi[s * nv + v], vmu_srli(p, 56));
            }
        }
        if (mode == 0)
            bv = vmu_loadu(b + VMU_VL * v);
        else if (mode == 1)
            bv = vmu_loadu_clear(b + VMU_VL * v, mlast);
        else
            bv = vmu_loadu_mask(b + VMU_VL * v, mlast);
        for (s = 0; s < r; s++)
        {
            p = vmu_mul(av[s], bv);
            lo[s * nv + v] = vmu_add(lo[s * nv + v], vmu_and(p, mask));
            hi[s * nv + v] = vmu_add(hi[s * nv + v], vmu_srli(p, 56));
        }

        if (++cnt == VECMAT_U32_FOLD)
        {
            for (v = 0; v < r * nv; v++)
            {
                hi[v] = vmu_add(hi[v], vmu_srli(lo[v], 56));
                lo[v] = vmu_and(lo[v], mask);
            }
            cnt = 0;
        }
    }

    return cnt;
}

/* the rows used are those of index < len */
FLINT_FORCE_INLINE void
vecmat_u32_block(ulong * const * c, vmu_t * lo, vmu_t * hi,
                 const ulong * const * a, const nmod_mat_t B, slong j0,
                 slong w, slong i0, slong i1, slong len, int r, slong nv)
{
    const nmod_t mod = B->mod;
    const slong wl = w - VMU_VL * (nv - 1);   /* columns of the last vector */
    const vmu_mask_t mlast = vmu_mask(wl);
    slong v;
    int s;

    for (s = 0; s < r; s++)
    {
        for (v = 0; v < nv - 1; v++)
            lo[s * nv + v] = vmu_loadu(c[s] + j0 + VMU_VL * v);
        lo[s * nv + v] = (wl == VMU_VL)
            ? vmu_loadu(c[s] + j0 + VMU_VL * v)
            : vmu_loadu_mask(c[s] + j0 + VMU_VL * v, mlast);
    }
    for (v = 0; v < r * nv; v++)
        hi[v] = vmu_zero();

    if (wl == VMU_VL)
        vecmat_u32_rows(lo, hi, a, B, j0, i0, i1, r, nv, 0, mlast, 0);
    else
    {
        slong is = i0;
        int cnt = 0;

#if VECMAT_U32_GCC
        /* the whole last vector of row i runs VMU_VL - wl entries past the
           end of the row, which stays inside B (whose rows used end at
           row len - 1) for i < len - ceil((VMU_VL - wl) / stride) */
        const slong stride = B->stride;

        if (stride >= B->c)
        {
            const slong e = VMU_VL - wl;
            is = FLINT_MAX(i0, FLINT_MIN(i1,
                     len - ((stride >= e) ? 1 : (e + stride - 1) / stride)));
            cnt = vecmat_u32_rows(lo, hi, a, B, j0, i0, is, r, nv, 1, mlast, 0);
        }
#else
        (void) len;
#endif
        vecmat_u32_rows(lo, hi, a, B, j0, is, i1, r, nv, 2, mlast, cnt);
    }

    for (s = 0; s < r; s++)
    {
        for (v = 0; v < nv; v++)
        {
            ulong l[VMU_VL], h[VMU_VL];
            int k;

            vmu_storeu(l, lo[s * nv + v]);
            vmu_storeu(h, hi[s * nv + v]);
            for (k = 0; k < VMU_VL && VMU_VL * v + k < w; k++)
            {
                ulong t1, t0, res;
                add_ssaaaa(t1, t0, h[k] >> 8, h[k] << 56, UWORD(0), l[k]);
                NMOD2_RED2(res, t1, t0, mod);
                c[s][j0 + VMU_VL * v + k] = res;
            }
        }
    }
}

#define VECMAT_U32_REGS_CASE(R, NV) \
    case 8 * (R) + (NV): \
        vecmat_u32_block(c, lo, hi, a, B, j0, w, i0, i1, len, R, NV); break;

static void
vecmat_u32_regs(ulong * const * c, const ulong * const * a, slong len,
                const nmod_mat_t B, int r)
{
    const slong ncols = B->c;
    const int nvmax = vecmat_u32_regs_nv[r - 1];
    const slong rows = FLINT_MAX(VECMAT_U32_ROWS,
                                 VECMAT_L2_BYTES / (8 * ncols));
    slong i0, i1;
    int s;

    for (s = 0; s < r; s++)
        _nmod_vec_zero(c[s], ncols);

    for (i0 = 0; i0 < len; i0 = i1)
    {
        slong j0 = 0;

        i1 = FLINT_MIN(len, i0 + rows);

        while (j0 < ncols)
        {
            vmu_t lo[VECMAT_ROWS * 8], hi[VECMAT_ROWS * 8];
            /* nv: nvmax, or for the last columns the largest power of two
               not above their number of vectors */
            const slong t = (ncols - j0 + VMU_VL - 1) / VMU_VL;
            int nv = nvmax;
            slong w;

            while (nv > t)
                nv /= 2;
            w = FLINT_MIN(ncols - j0, VMU_VL * nv);

            switch (8 * r + nv)
            {
                VECMAT_U32_REGS_CASE(1, 1)
                VECMAT_U32_REGS_CASE(1, 2)
                VECMAT_U32_REGS_CASE(1, 4)
                VECMAT_U32_REGS_CASE(2, 1)
                VECMAT_U32_REGS_CASE(2, 2)
                VECMAT_U32_REGS_CASE(3, 1)
                VECMAT_U32_REGS_CASE(4, 1)
#if VMU_VL == 8
                VECMAT_U32_REGS_CASE(1, 8)
                VECMAT_U32_REGS_CASE(2, 4)
                VECMAT_U32_REGS_CASE(3, 2)
                VECMAT_U32_REGS_CASE(3, 4)
                VECMAT_U32_REGS_CASE(4, 2)
#endif
                default: FLINT_UNREACHABLE;
            }

            j0 += w;
        }
    }
}

/* many columns: blocks whose accumulators (16 bytes per column and row)
   stay in L1, the rows of B streamed in runs of a few KB */
static void
vecmat_u32_wide(ulong * const * c, const ulong * const * a, slong len,
                const nmod_mat_t B, int r)
{
    const slong ncols = B->c;
    const slong wmax = VMU_VL * (VECMAT_U32_WIDTH / (VMU_VL * r));
    vmu_t lo[VECMAT_U32_WIDTH / VMU_VL], hi[VECMAT_U32_WIDTH / VMU_VL];
    slong j0 = 0;
    int s;

    for (s = 0; s < r; s++)
        _nmod_vec_zero(c[s], ncols);

    while (j0 < ncols)
    {
        const slong w = FLINT_MIN(ncols - j0, wmax);
        const slong nv = (w + VMU_VL - 1) / VMU_VL;

        switch (r)
        {
            case 1:
                vecmat_u32_block(c, lo, hi, a, B, j0, w, 0, len, len, 1, nv);
                break;
            case 2:
                vecmat_u32_block(c, lo, hi, a, B, j0, w, 0, len, len, 2, nv);
                break;
            case 3:
                vecmat_u32_block(c, lo, hi, a, B, j0, w, 0, len, len, 3, nv);
                break;
            case 4:
                vecmat_u32_block(c, lo, hi, a, B, j0, w, 0, len, len, 4, nv);
                break;
            default: FLINT_UNREACHABLE;
        }

        j0 += w;
    }
}

#endif  /* u32 */

#if FLINT_BITS == 64

#include "nmod_mat/mul_fp_vec.h"

#if NMOD_MAT_HAVE_FPV
#if defined(FPV_GENERIC)
# error "NMOD_MAT_HAVE_FPV and the backend selection of mul_fp_vec.h disagree"
#endif
#define VECMAT_HAVE_FP50 1

/* balanced representative of 0 <= x < n as a double */
FLINT_FORCE_INLINE double
vecmat_sym(ulong x, ulong n)
{
    return (x > n / 2) ? (double) (slong) (x - n) : (double) x;
}

/* columns per block for one row: the accumulators (8 bytes per column
   and row) stay in L1 */
#define VECMAT_FP50_WIDTH 512

/*
    One block of nv vectors of columns j0 <= j < j0 + w for r rows, with
    the accumulators acc[s * nv + v] in L1 (r is a constant at each call
    site, nv is not; register accumulators were tried and bring nothing:
    on AVX2 this loop is bound by the issue width, at about 10 vector
    operations per vector of products). Each term is in (-9/8 n, 9/8 n)
    (fpv_mulmod with |a| <= n/2), and the accumulators are reduced every
    four rows, staying below 5.01 n < 2^53 in absolute value.

    The full vectors go through vecmat_fp50_rows, and a last partial vector
    (w not a multiple of FPV_VL) through a pass of its own over the rows,
    vecmat_fp50_rows_last, which reads it with fpv_load_u64_first (a masked
    load on x86, zero above). It used to be copied to a zero-padded buffer
    in each row, and the wide load of these narrow stores could not be
    forwarded: 1.7-3.5x slower for ncols = 1-7 and 9-15 than now (Cascade
    Lake). With the partial vector in the loop of the full ones, these got
    1.1-1.25x slower (with a test, the compiler made all the loads masked;
    with a copy of the loop per case, it was slower for other reasons).
*/

/* the full vectors v0 <= v < v1 of the rows; acc has nv vectors per row */
FLINT_FORCE_INLINE void
vecmat_fp50_rows(fpv * acc, const ulong * const * a, const nmod_mat_t B,
                 slong j0, slong len, int r, slong nv, slong v0, slong v1,
                 fpv nvec, fpv ninv)
{
    const ulong n = B->mod.n;
    fpv av[VECMAT_ROWS], bv;
    slong i, v;
    int s, cnt = 0;

    for (i = 0; i < len; i++)
    {
        const ulong * b = nmod_mat_entry_ptr(B, i, j0);

        for (s = 0; s < r; s++)
            av[s] = fpv_set1(vecmat_sym(a[s][i], n));

        for (v = v0; v < v1; v++)
        {
            bv = fpv_load_u64(b + FPV_VL * v);
            for (s = 0; s < r; s++)
                acc[s * nv + v] = fpv_add(acc[s * nv + v],
                                          fpv_mulmod(av[s], bv, nvec, ninv));
        }

        if (++cnt == 4)
        {
            for (s = 0; s < r; s++)
                for (v = v0; v < v1; v++)
                    acc[s * nv + v] = fpv_reduce_pm1n(acc[s * nv + v], nvec, ninv);
            cnt = 0;
        }
    }
}

/* the last vector v, with wl < FPV_VL columns, of the rows */
FLINT_FORCE_INLINE void
vecmat_fp50_rows_last(fpv * acc, const ulong * const * a, const nmod_mat_t B,
                      slong j0, slong len, int r, slong nv, slong v, slong wl,
                      fpv nvec, fpv ninv)
{
    const ulong n = B->mod.n;
    fpv bv;
    slong i;
    int s, cnt = 0;

    for (i = 0; i < len; i++)
    {
        bv = fpv_load_u64_first(nmod_mat_entry_ptr(B, i, j0) + FPV_VL * v, wl);
        for (s = 0; s < r; s++)
            acc[s * nv + v] = fpv_add(acc[s * nv + v],
                fpv_mulmod(fpv_set1(vecmat_sym(a[s][i], n)), bv, nvec, ninv));

        if (++cnt == 4)
        {
            for (s = 0; s < r; s++)
                acc[s * nv + v] = fpv_reduce_pm1n(acc[s * nv + v], nvec, ninv);
            cnt = 0;
        }
    }
}

FLINT_FORCE_INLINE void
vecmat_fp50_block(ulong * const * c, fpv * acc, const ulong * const * a,
                  const nmod_mat_t B, slong j0, slong w, slong len, int r,
                  slong nv, fpv nvec, fpv ninv)
{
    const slong wl = w - FPV_VL * (nv - 1);   /* columns of the last vector */
    const slong nfull = (wl == FPV_VL) ? nv : nv - 1;
    ulong buf[FPV_VL];
    slong j, v;
    int s;

    for (v = 0; v < r * nv; v++)
        acc[v] = fpv_zero();

    if (nfull > 0)
        vecmat_fp50_rows(acc, a, B, j0, len, r, nv, 0, nfull, nvec, ninv);
    if (nfull < nv)
        vecmat_fp50_rows_last(acc, a, B, j0, len, r, nv, nv - 1, wl, nvec, ninv);

    for (s = 0; s < r; s++)
    {
        for (v = 0; v < nv - 1; v++)
            fpv_store_u64(c[s] + j0 + FPV_VL * v,
                fpv_reduce_0n(fpv_reduce_pm1n(acc[s * nv + v], nvec, ninv), nvec));
        fpv_store_u64(buf,
            fpv_reduce_0n(fpv_reduce_pm1n(acc[s * nv + v], nvec, ninv), nvec));
        for (j = 0; j < wl; j++)
            c[s][j0 + FPV_VL * v + j] = buf[j];
    }
}

static void
vecmat_fp50(ulong * const * c, const ulong * const * a, slong len,
            const nmod_mat_t B, int r)
{
    const slong ncols = B->c;
    const fpv nvec = fpv_set1((double) B->mod.n);
    const fpv ninv = fpv_set1(1.0 / (double) B->mod.n);
    const slong wmax = FPV_VL * (VECMAT_FP50_WIDTH / (FPV_VL * r));
    fpv acc[VECMAT_FP50_WIDTH / FPV_VL];
    slong j0 = 0;

    while (j0 < ncols)
    {
        const slong w = FLINT_MIN(ncols - j0, wmax);
        const slong nv = (w + FPV_VL - 1) / FPV_VL;

        switch (r)
        {
            case 1: vecmat_fp50_block(c, acc, a, B, j0, w, len, 1, nv, nvec, ninv); break;
            case 2: vecmat_fp50_block(c, acc, a, B, j0, w, len, 2, nv, nvec, ninv); break;
            case 3: vecmat_fp50_block(c, acc, a, B, j0, w, len, 3, nv, nvec, ninv); break;
            case 4: vecmat_fp50_block(c, acc, a, B, j0, w, len, 4, nv, nvec, ninv); break;
            default: FLINT_UNREACHABLE;
        }

        j0 += w;
    }
}

#else
#define VECMAT_HAVE_FP50 0
#endif

#endif  /* FLINT_BITS == 64 */

#if FLINT_BITS == 64 && NMOD_MAT_HAVE_VECMAT_U64

/*
    The other moduli (n >= 2^50, or n > 2^52 with IFMA, and also
    2^32 < n < 2^50 without a vector FMA): the 32 x 32 -> 64-bit products
    (vpmuludq) of limbs of the entries are accumulated in 64-bit lanes
    with no reduction until the end, as in the split limbs dot products of
    nmod_vec/dot_split_limbs.c. Let b be the bit length of n - 1.

    - split, for b <= 61: x = x1 2^s + x0 with s = ceil(b/2), x0 <= l0 =
      2^s - 1 and x1 <= l1 = 2^(b-s) - 1, so that

          a x = a0 x0 + 2^s (a0 x1 + a1 x0) + 2^2s a1 x1;

      L, M, H accumulate the three coefficients, each term being at most
      P = max(l0^2, 2 l0 l1). Every F rows each lane carries its high bits
      into the next accumulator,

          M += L >> s, H += M >> s, Hh += H >> 32,

      keeping the low s, s and 32 bits of L, M, H. This preserves
      L + 2^s M + 2^2s (H + 2^32 Hh) and leaves L, M < 2^s, H < 2^32, so
      that before the next fold the lanes are below 2^32 + F P, plus
      carries below 2^(64-s): F is the largest number of rows with
      2^32 + F P + 2^(64-s) < 2^64 (F = 4 for b = 61, doubling as b
      decreases: 2048 for b = 52). Four products and four additions per
      vector of products, plus the split of the vector of B, shared by
      the rows, and the folds.
    - halves, for b >= 62: 32-bit limbs, whose products fill the lanes;
      each product p is accumulated as X += p (modulo 2^64) and
      Xh += p >> 32, for (X, Xh) = (L, Lh), (M, Mh), (H, Hh). The sum of
      the low halves of the products is then X - 2^32 Xh modulo 2^64,
      exact while below 2^64: no folds, but six accumulators and eight
      more operations per vector of products.

    At the end the totals (below len n^2) are formed on three words in
    the vector registers (in halves, by carrying through their 32-bit
    digits), shifted by the normalization of n and reduced by the scalar
    code: two steps of the division with the precomputed inverse of
    NMOD_RED2, one when the totals are below n 2^64.

    Register blocks of nv vectors of columns share the loads of B between
    up to 4 rows. When the rows of B fit in half of VECMAT_L2_BYTES, or in a
    single register block, each block goes through all of them. Otherwise
    the columns are cut into blocks of at most VECMAT_U64_WIDTH, whose
    accumulators are kept in a buffer, and each of these goes through the
    rows VECMAT_U64_ROWS at a time, its register blocks one after the other:
    B is then read as VECMAT_U64_ROWS streams of consecutive vectors, which
    the hardware prefetchers follow, whereas narrow register blocks going
    through all the rows of a B in L3 or in memory wait for its latency
    (Cascade Lake, 1 to 4 rows: 2.5-7.5x slower for B of 8 to 64 MB).
*/

/* B in at most this many columns at a time when its rows are read by
   chunks: the accumulators (32 or 48 bytes per column and row of a) stay
   in L2 */
#ifndef VECMAT_U64_WIDTH
# define VECMAT_U64_WIDTH 2048
#endif
/* rows per chunk, when the rows are read by chunks */
#ifndef VECMAT_U64_ROWS
# define VECMAT_U64_ROWS 16
#endif
/* most rows of B per call of vecmat_u64: in halves, the 32-bit digits of the
   totals stay below 2^64 up to 2^30 rows */
#ifndef VECMAT_U64_MAX_LEN
# define VECMAT_U64_MAX_LEN (WORD(1) << 28)
#endif

/* rows sharing the loads of B, and vectors of columns per register block
   for r rows (split: 3 r nv accumulators, Hh being only used by the
   folds, halves: 6 r nv, plus 2 r + 4 registers) */
#if VMU_VL == 8
# define VECMAT_U64_GROUP(split) 4
static const int vecmat_u64_nv[2][4] = {{4, 2, 1, 1}, {8, 3, 2, 1}};
#else
# define VECMAT_U64_GROUP(split) ((split) ? 2 : 1)
static const int vecmat_u64_nv[2][4] = {{1, 1, 1, 1}, {2, 1, 1, 1}};
#endif
#define VECMAT_U64_GROUP_MAX 4
#define VECMAT_U64_NV_MAX 8

typedef struct
{
    nmod_t mod;
    ulong nn;       /* n << norm */
    slong F;        /* split: rows between two folds */
    int s;          /* split: bits of the low limb */
    int split;
    int two;        /* the totals are below n 2^64 */
}
vecmat_u64_ctx_struct;

/* for the products of len rows */
static void
vecmat_u64_ctx_init(vecmat_u64_ctx_struct * C, nmod_t mod, slong len)
{
    const int b = FLINT_BIT_COUNT(mod.n - 1);
    ulong t2, t1, t0;

    FLINT_ASSERT(mod.n >= 2 && len <= VECMAT_U64_MAX_LEN);

    C->mod = mod;
    C->nn = mod.n << mod.norm;
    C->split = (b <= 61);
    C->s = (b + 1) / 2;
    C->F = 0;

    if (C->split)
    {
        const ulong l0 = (UWORD(1) << C->s) - 1;
        const ulong l1 = (UWORD(1) << (b - C->s)) - 1;
        const ulong P = FLINT_MAX(l0 * l0, 2 * l0 * l1);
        /* 2^64 - 1 - 2^32 - 2^(64-s) */
        const ulong room = UWORD_MAX - (UWORD(1) << 32)
                                     - (UWORD(1) << (FLINT_BITS - C->s));

        C->F = (slong) FLINT_MIN(room / P, (ulong) VECMAT_U64_MAX_LEN);
    }

    /* the totals are at most len (n - 1)^2 = t2 2^128 + t1 2^64 + ... */
    umul_ppmm(t1, t0, mod.n - 1, mod.n - 1);
    umul_ppmm(t2, t1, t1, (ulong) len);
    add_ssaaaa(t2, t1, t2, t1, UWORD(0), n_mulhi(t0, (ulong) len));
    C->two = (t2 == 0 && t1 < mod.n);
}

/* (u1 2^64 + u0) mod nn for u1 < nn: NMOD_RED2 on normalized words */
FLINT_FORCE_INLINE ulong
vecmat_u64_red2(ulong u1, ulong u0, ulong nn, ulong ninv)
{
    ulong q1, q0, r;

    umul_ppmm(q1, q0, ninv, u1);
    add_ssaaaa(q1, q0, q1, q0, u1, u0);
    r = u0 - (q1 + 1) * nn;
    if (r > q0)
        r += nn;
    if (r >= nn)
        r -= nn;
    return r;
}

/* c[j] = (t2[j] 2^128 + t1[j] 2^64 + t0[j]) mod n for j < cnt, where the
   totals are below n 2^128 (below n 2^64 if C->two) */
FLINT_STATIC_NOINLINE void
vecmat_u64_reduce(ulong * c, slong cnt, const ulong * t2, const ulong * t1,
                  const ulong * t0, const vecmat_u64_ctx_struct * C)
{
    const int norm = C->mod.norm;
    const vmu_t sh = vmu_set1(norm), shc = vmu_set1(FLINT_BITS - norm);
    const vmu_t x2 = vmu_loadu(t2), x1 = vmu_loadu(t1), x0 = vmu_loadu(t0);
    ulong u2[VMU_VL], u1[VMU_VL], u0[VMU_VL];
    slong j;

    /* shifted left by norm (by 64 to the right when norm = 0: zero) */
    vmu_storeu(u2, vmu_or(vmu_sllv(x2, sh), vmu_srlv(x1, shc)));
    vmu_storeu(u1, vmu_or(vmu_sllv(x1, sh), vmu_srlv(x0, shc)));
    vmu_storeu(u0, vmu_sllv(x0, sh));

    if (C->two)
        for (j = 0; j < cnt; j++)
            c[j] = vecmat_u64_red2(u1[j], u0[j], C->nn, C->mod.ninv) >> norm;
    else
        for (j = 0; j < cnt; j++)
            c[j] = vecmat_u64_red2(
                       vecmat_u64_red2(u2[j], u1[j], C->nn, C->mod.ninv),
                       u0[j], C->nn, C->mod.ninv) >> norm;
}

/*
    One block of nv vectors of columns jb <= j < jb + w for r rows and the
    rows i0 <= i < i1 of B, where B and the c[t] point at the first column
    of the block of at most VECMAT_U64_WIDTH columns, whose accumulator state
    is st: acc[k] of row t at st + (t nacc + k) sw + jb. The accumulators
    start from zero (first) or from st; at the end they are stored to st,
    or (last) their totals are reduced into c. The last vector, of
    wl = w - VMU_VL (nv - 1) columns, is loaded with a mask.

    r, nv and split are constants at each call site, so that acc is
    register allocated: split uses acc[0..3] = L, M, H, Hh (Hh, only
    touched by the folds, being spilled if needed), halves acc[0..5] = L,
    Lh, M, Mh, H, Hh.
*/
FLINT_FORCE_INLINE void
vecmat_u64_block(ulong * const * c, ulong * st, slong sw,
                 const ulong * const * a, const ulong * B, slong stride,
                 slong jb, slong w, slong i0, slong i1, int first, int last,
                 int r, int nv, int split, const vecmat_u64_ctx_struct * C)
{
    const int nacc = split ? 4 : 6;
    const slong wl = w - VMU_VL * (nv - 1);
    const vmu_t sv = vmu_set1(C->s), mk = vmu_set1((UWORD(1) << C->s) - 1);
    const vmu_t m32 = vmu_set1(UWORD(0xFFFFFFFF));
    const vmu_mask_t m = vmu_mask(wl);
    vmu_t acc[6][VECMAT_U64_GROUP_MAX * VECMAT_U64_NV_MAX];
    vmu_t a0[VECMAT_U64_GROUP_MAX], a1[VECMAT_U64_GROUP_MAX], bv, b1, p;
    slong i, i2;
    int t, k, v;

    for (t = 0; t < r; t++)
        for (k = 0; k < nacc; k++)
            for (v = 0; v < nv; v++)
                acc[k][t * nv + v] = first ? vmu_zero()
                    : vmu_loadu(st + (t * nacc + k) * sw + jb + VMU_VL * v);

    for (i = i0; i < i1; i = i2)
    {
        i2 = split ? FLINT_MIN(i1, i + C->F) : i1;

        for ( ; i < i2; i++)
        {
            const ulong * b = B + i * stride + jb;

            /* the limbs of a[t][i], split in the vector registers (in
               halves, vpmuludq reads the low half of a0) */
            for (t = 0; t < r; t++)
            {
                a0[t] = vmu_set1(a[t][i]);
                if (split)
                {
                    a1[t] = vmu_srlv(a0[t], sv);
                    a0[t] = vmu_and(a0[t], mk);
                }
                else
                    a1[t] = vmu_srli(a0[t], 32);
            }

            for (v = 0; v < nv; v++)
            {
                if (v < nv - 1 || wl == VMU_VL)
                    bv = vmu_loadu(b + VMU_VL * v);
                else
                    bv = vmu_loadu_mask(b + VMU_VL * v, m);

                if (split)
                {
                    const vmu_t b0 = vmu_and(bv, mk);
                    b1 = vmu_srlv(bv, sv);
                    for (t = 0; t < r; t++)
                    {
                        const int x = t * nv + v;
                        acc[0][x] = vmu_add(acc[0][x], vmu_mul(a0[t], b0));
                        acc[2][x] = vmu_add(acc[2][x], vmu_mul(a1[t], b1));
                        p = vmu_add(vmu_mul(a0[t], b1), vmu_mul(a1[t], b0));
                        acc[1][x] = vmu_add(acc[1][x], p);
                    }
                }
                else
                {
                    b1 = vmu_srli(bv, 32);
                    for (t = 0; t < r; t++)
                    {
                        const int x = t * nv + v;
                        p = vmu_mul(a0[t], bv);
                        acc[0][x] = vmu_add(acc[0][x], p);
                        acc[1][x] = vmu_add(acc[1][x], vmu_srli(p, 32));
                        p = vmu_mul(a1[t], b1);
                        acc[4][x] = vmu_add(acc[4][x], p);
                        acc[5][x] = vmu_add(acc[5][x], vmu_srli(p, 32));
                        p = vmu_mul(a0[t], b1);
                        acc[2][x] = vmu_add(acc[2][x], p);
                        acc[3][x] = vmu_add(acc[3][x], vmu_srli(p, 32));
                        p = vmu_mul(a1[t], bv);
                        acc[2][x] = vmu_add(acc[2][x], p);
                        acc[3][x] = vmu_add(acc[3][x], vmu_srli(p, 32));
                    }
                }
            }
        }

        if (split)
        {
            for (v = 0; v < r * nv; v++)
            {
                acc[1][v] = vmu_add(acc[1][v], vmu_srlv(acc[0][v], sv));
                acc[0][v] = vmu_and(acc[0][v], mk);
                acc[2][v] = vmu_add(acc[2][v], vmu_srlv(acc[1][v], sv));
                acc[1][v] = vmu_and(acc[1][v], mk);
                acc[3][v] = vmu_add(acc[3][v], vmu_srli(acc[2][v], 32));
                acc[2][v] = vmu_and(acc[2][v], m32);
            }
        }
    }

    if (!last)
    {
        for (t = 0; t < r; t++)
            for (k = 0; k < nacc; k++)
                for (v = 0; v < nv; v++)
                    vmu_storeu(st + (t * nacc + k) * sw + jb + VMU_VL * v,
                               acc[k][t * nv + v]);
        return;
    }

    for (t = 0; t < r; t++)
    {
        for (v = 0; v < nv; v++)
        {
            const int x = t * nv + v;
            ulong t2[VMU_VL], t1[VMU_VL], t0[VMU_VL];

            if (split)
            {
                /* L + 2^s M + 2^2s (H + 2^32 Hh) with L, M < 2^s and
                   H < 2^32, and 2 <= 2s <= 62 */
                const vmu_t s2 = vmu_set1(2 * C->s);
                const vmu_t s2c = vmu_set1(FLINT_BITS - 2 * C->s);
                const vmu_t lo = vmu_or(acc[0][x], vmu_sllv(acc[1][x], sv));
                const vmu_t y0 = vmu_or(acc[2][x], vmu_slli(acc[3][x], 32));
                const vmu_t y1 = vmu_srli(acc[3][x], 32);

                vmu_storeu(t0, vmu_or(lo, vmu_sllv(y0, s2)));
                vmu_storeu(t1, vmu_or(vmu_srlv(y0, s2c), vmu_sllv(y1, s2)));
                vmu_storeu(t2, vmu_srlv(y1, s2c));
            }
            else
            {
                /* digits of weight 2^0, 2^32, 2^64, 2^96, with the sums of
                   the low halves X - 2^32 Xh */
                vmu_t d0 = vmu_sub(acc[0][x], vmu_slli(acc[1][x], 32));
                vmu_t d1 = vmu_sub(acc[2][x], vmu_slli(acc[3][x], 32));
                vmu_t d2 = vmu_sub(acc[4][x], vmu_slli(acc[5][x], 32));
                vmu_t d3 = acc[5][x];

                d1 = vmu_add(d1, acc[1][x]);
                d2 = vmu_add(d2, acc[3][x]);

                d1 = vmu_add(d1, vmu_srli(d0, 32));
                d2 = vmu_add(d2, vmu_srli(d1, 32));
                d3 = vmu_add(d3, vmu_srli(d2, 32));
                vmu_storeu(t0, vmu_or(vmu_and(d0, m32), vmu_slli(d1, 32)));
                vmu_storeu(t1, vmu_or(vmu_and(d2, m32), vmu_slli(d3, 32)));
                vmu_storeu(t2, vmu_srli(d3, 32));
            }

            vecmat_u64_reduce(c[t] + jb + VMU_VL * v,
                              (v < nv - 1) ? VMU_VL : wl, t2, t1, t0, C);
        }
    }
}

#define VECMAT_U64_CASE(R, NV, SPLIT) \
    case 64 * (SPLIT) + 8 * (R) + (NV): \
        vecmat_u64_block(cw, st, sw, a, Bw, stride, jb, w, i0, i1, first, \
                         last, R, NV, SPLIT, C); \
        break;

/* c[t] = a[t] * B for t < r rows, where B (len rows of ncols entries and
   stride) has at most VECMAT_U64_MAX_LEN rows */
static void
vecmat_u64(ulong * const * c, const ulong * const * a, const ulong * B,
           slong stride, slong len, slong ncols, int r,
           const vecmat_u64_ctx_struct * C)
{
    const int split = C->split, nacc = split ? 4 : 6;
    const int nvmax = vecmat_u64_nv[split][r - 1];
    const slong W = FLINT_MIN(ncols, VECMAT_U64_WIDTH), sw = W + VMU_VL;
    slong rows, j0, i0, i1;
    ulong * st = NULL;
    TMP_INIT;

    /* B at most half of VECMAT_L2_BYTES: with 1 MB of L2 (Cascade Lake), a
       B of 1 MB was read up to 1.6x slower in one pass than by chunks */
    if (ncols <= VMU_VL * nvmax || len <= VECMAT_L2_BYTES / (16 * ncols))
        rows = len;
    else
        rows = VECMAT_U64_ROWS;

    TMP_START;

    if (rows < len)
        st = TMP_ARRAY_ALLOC(r * nacc * sw, ulong);

    for (j0 = 0; j0 < ncols; j0 += W)
    {
        const slong wj = FLINT_MIN(W, ncols - j0);
        const ulong * Bw = B + j0;
        ulong * cw[VECMAT_U64_GROUP_MAX];
        int t;

        for (t = 0; t < r; t++)
            cw[t] = c[t] + j0;

        for (i0 = 0; i0 < len; i0 = i1)
        {
            const int first = (i0 == 0);
            int last;
            slong jb, w;

            i1 = FLINT_MIN(len, i0 + rows);
            last = (i1 == len);

            for (jb = 0; jb < wj; jb += w)
            {
                /* nv: nvmax, or for the last columns the largest power of
                   two not above their number of vectors */
                const slong tv = (wj - jb + VMU_VL - 1) / VMU_VL;
                int nv = nvmax;

                while (nv > tv)
                    nv = (nv == 3) ? 2 : nv / 2;
                w = FLINT_MIN(wj - jb, VMU_VL * nv);

                switch (64 * split + 8 * r + nv)
                {
                    VECMAT_U64_CASE(1, 1, 1)
                    VECMAT_U64_CASE(1, 2, 1)
                    VECMAT_U64_CASE(2, 1, 1)
                    VECMAT_U64_CASE(1, 1, 0)
#if VMU_VL == 8
                    VECMAT_U64_CASE(1, 4, 1)
                    VECMAT_U64_CASE(1, 8, 1)
                    VECMAT_U64_CASE(2, 2, 1)
                    VECMAT_U64_CASE(2, 3, 1)
                    VECMAT_U64_CASE(3, 1, 1)
                    VECMAT_U64_CASE(3, 2, 1)
                    VECMAT_U64_CASE(4, 1, 1)
                    VECMAT_U64_CASE(1, 2, 0)
                    VECMAT_U64_CASE(1, 4, 0)
                    VECMAT_U64_CASE(2, 1, 0)
                    VECMAT_U64_CASE(2, 2, 0)
                    VECMAT_U64_CASE(3, 1, 0)
                    VECMAT_U64_CASE(4, 1, 0)
#endif
                    default: FLINT_UNREACHABLE;
                }
            }
        }
    }

    TMP_END;
}

#endif  /* u64 */

int
_nmod_mat_mul_rows_simd(ulong * const * c, const ulong * const * a, slong r,
                        slong len, const nmod_mat_t B)
{
    FLINT_ASSERT(1 <= r && r <= NMOD_MAT_MUL_ROWS_MAX);
    FLINT_ASSERT(len >= 1 && len <= B->r && B->c >= 1);

    /* no vectorized path on this target: avoid unused parameter warning */
    (void) c; (void) a; (void) r; (void) len; (void) B;

#if FLINT_BITS == 64
#if NMOD_MAT_HAVE_MUL_U52
    if (B->mod.n <= (UWORD(1) << 52))
    {
        u52_ctx_struct ctx;
        u52_ctx_init(&ctx, B->mod.n);
        const u52_consts C = u52_consts_init(&ctx);

        slong s, g;

        if (B->c <= VECMAT_REGS_MAX_COLS)
        {
            /* with many columns and the high halves, groups of 3 rows keep
               the register blocks 4 vectors wide (see vecmat_u52_regs) */
            g = (!ctx.lo_only && B->c > VECMAT_REGS_WIDE_COLS)
                ? 3 : VECMAT_ROWS;
            for (s = 0; s < r; s += g)
                vecmat_u52_regs(c + s, a + s, len, B, (int) FLINT_MIN(r - s, g),
                                !ctx.lo_only, &C);
        }
        else
        {
            for (s = 0; s < r; s += 4)
                vecmat_u52_wide(c + s, a + s, len, B, (int) FLINT_MIN(r - s, 4),
                                !ctx.lo_only, &C);
        }
        return 1;
    }
#endif
#if NMOD_MAT_HAVE_VECMAT_U32
    if (B->mod.n <= (UWORD(1) << 32)
            && !NMOD_MAT_VECMAT_U32_TO_FP50(len, B->c))
    {
        /* groups of at most 4 rows */
        slong s;
        for (s = 0; s < r; s += 4)
        {
            const int rs = (int) FLINT_MIN(r - s, 4);
            if (B->c <= VECMAT_U32_REGS_MAX_COLS)
                vecmat_u32_regs(c + s, a + s, len, B, rs);
            else
                vecmat_u32_wide(c + s, a + s, len, B, rs);
        }
        return 1;
    }
#endif
#if VECMAT_HAVE_FP50
    if (B->mod.n < (UWORD(1) << 50))
    {
        slong s;
        for (s = 0; s < r; s += 4)
            vecmat_fp50(c + s, a + s, len, B, (int) FLINT_MIN(r - s, 4));
        return 1;
    }
#endif
#if NMOD_MAT_HAVE_VECMAT_U64
    {
        /* all the other moduli; more than VECMAT_U64_MAX_LEN rows by parts,
           the next ones going through a temporary */
        const slong ncols = B->c;
        vecmat_u64_ctx_struct C;
        ulong * cc[NMOD_MAT_MUL_ROWS_MAX];
        const ulong * aa[NMOD_MAT_MUL_ROWS_MAX];
        ulong * t = NULL;
        slong s, g, i0, i1;

        for (i0 = 0; i0 < len; i0 = i1)
        {
            i1 = FLINT_MIN(len, i0 + VECMAT_U64_MAX_LEN);
            vecmat_u64_ctx_init(&C, B->mod, i1 - i0);
            g = VECMAT_U64_GROUP(C.split);

            if (i0 > 0 && t == NULL)
                t = flint_malloc(r * ncols * sizeof(ulong));

            for (s = 0; s < r; s++)
            {
                cc[s] = (i0 == 0) ? c[s] : t + s * ncols;
                aa[s] = a[s] + i0;
            }

            for (s = 0; s < r; s += g)
                vecmat_u64(cc + s, aa + s, nmod_mat_entry_ptr(B, i0, 0),
                           B->stride, i1 - i0, ncols,
                           (int) FLINT_MIN(r - s, g), &C);

            if (i0 > 0)
                for (s = 0; s < r; s++)
                    _nmod_vec_add(c[s], c[s], cc[s], ncols, B->mod);
        }

        flint_free(t);
        return 1;
    }
#endif
#endif

    return 0;
}

/*
    Whether nmod_mat_nmod_vec_mul uses _nmod_mat_mul_rows_simd for a vector
    of length len >= 2 and ncols columns, rather than its fallback (below).
    Every SIMD path pays a per-call setup and a final reduction of each
    entry, and fills whole vectors of columns: for narrow B, the fallback
    is faster, by up to 2-5x at ncols = 1 (all paths, all machines
    measured). The thresholds, per path in the order of
    _nmod_mat_mul_rows_simd, were fitted on Zen 4 and Ice Lake SP (u52,
    u64), Cascade Lake (u32, fp50, u64) and Apple M4 (fp50), for moduli
    of 20 to 64 bits, len = 2..1024 and ncols = 1..2048.

    - u52 (IFMA, n <= 2^52): ncols >= 16 for len = 2, 8 for len = 3, 5 for
      len = 4..5, and from len = 6 on 4 for n <= 2^26 (one IFMA per vector)
      and 3 above;
    - u32 (n <= 2^32): len >= 4 and ncols >= 4, except for short vectors
      (len <= 6, or len <= min(ncols, 32) from 16 columns on), for which
      fp50 replaces it (see NMOD_MAT_VECMAT_U32_TO_FP50);
    - fp50 (n < 2^50): ncols >= 16, or ncols >= 5 and len <= 16 with
      AVX2 or AVX-512; ncols >= 32, or ncols >= 4 and len <= 12 with NEON
      (for longer vectors and fewer columns the classical product, with
      the dot products of nmod_vec, is faster);
    - u64 (the other moduli): ncols >= 5 and len >= 4 for n <= 2^61
      (split limbs), len >= 6 above (halves).
*/
FLINT_FORCE_INLINE int
vecmat_use_rows_simd(slong len, slong ncols, ulong n)
{
#if FLINT_BITS == 64
#if NMOD_MAT_HAVE_MUL_U52
    if (n <= (UWORD(1) << 52))
        return ncols >= ((len == 2) ? 16 : (len == 3) ? 8 : (len <= 5) ? 5
                         : (n <= (UWORD(1) << 26)) ? 4 : 3);
#endif
#if NMOD_MAT_HAVE_VECMAT_U32
    if (n <= (UWORD(1) << 32) && !NMOD_MAT_VECMAT_U32_TO_FP50(len, ncols))
        return len >= 4 && ncols >= 4;
#endif
#if VECMAT_HAVE_FP50
    if (n < (UWORD(1) << 50))
#if defined(FPV_NEON)
        return ncols >= 32 || (ncols >= 4 && len <= 12);
#else
        return ncols >= 16 || (ncols >= 5 && len <= 16);
#endif
#endif
#if NMOD_MAT_HAVE_VECMAT_U64
    if (n <= (UWORD(1) << 61))
        return ncols >= 5 && len >= 4;
    else
        return ncols >= 5 && len >= 6;
#endif
#endif
    (void) len; (void) ncols; (void) n;
    return 0;
}

/*
    The length from which the fallback of nmod_mat_nmod_vec_mul uses the
    classical product rather than one scalar_addmul per entry of a: on x86
    with AVX2, 6 for n <= 2^32, 8 up to 2^62, 10 above (cost of the
    two-word dot products of the classical code); elsewhere 6 up to 2^62
    and 8 above (Apple M4). #2904 used 12 for all moduli, which cost
    1.05-1.8x at len = 6-10 for narrow B.
*/
FLINT_FORCE_INLINE slong
vecmat_classical_min_len(ulong n)
{
#if FLINT_BITS == 64 && NMOD_MAT_HAVE_VECMAT_U64
    if (n <= (UWORD(1) << 32))
        return 6;
    else if (n <= (UWORD(1) << 62))
        return 8;
    else
        return 10;
#elif FLINT_BITS == 64
    return (n <= (UWORD(1) << 62)) ? 6 : 8;
#else
    (void) n;
    return 6;
#endif
}

/* c[j] = sum a[i] B[i, j], one dot product per column (out of line: the
   inlined dot products made the other paths of nmod_mat_nmod_vec_mul
   slower on short vectors) */
FLINT_STATIC_NOINLINE void
vecmat_dots(ulong * c, const ulong * a, slong len, const nmod_mat_t B)
{
    const dot_params_t params = _nmod_vec_dot_params(len, B->mod);
    slong j;

    if (B->c == 1 && B->stride == 1)
        c[0] = _nmod_vec_dot(a, B->entries, len, B->mod, params);
    else
        for (j = 0; j < B->c; j++)
            c[j] = _nmod_vec_dot_strided(a, 1, B->entries + j, B->stride,
                                         len, B->mod, params);
}

#ifndef VECMAT_NARROW_BYTES
# define VECMAT_NARROW_BYTES 16384
#endif

void nmod_mat_nmod_vec_mul(
    ulong * c,
    const ulong * a, slong alen,
    const nmod_mat_t B)
{
    slong i, blk, nblk;
    slong len = FLINT_MIN(B->r, alen);
    slong ncols = B->c;
    nmod_mat_t A1, B1, C1;
    ulong * t;
    TMP_INIT;

    /* scalar_addmul wants non-empty */
    if (ncols < 1)
        return;

    if (len <= 0)
    {
        _nmod_vec_zero(c, ncols);
        return;
    }

    if (len >= 2 && vecmat_use_rows_simd(len, ncols, B->mod.n)
            && _nmod_mat_mul_rows_simd(&c, &a, 1, len, B))
        return;

    /* narrow B: one dot product per column (scalar_addmul and the classical
       product have a setup cost per row or per block of rows), for one
       column from length 3 (2 up to 2^62), for up to 8 columns from length
       4 (8 above 2^62) to 32 */
    if (ncols <= 8)
    {
        const int big = (B->mod.n > (UWORD(1) << (FLINT_BITS - 2)));

        if (ncols == 1 ? (len >= 3 || !big)
                       : (len >= (big ? 8 : 4) && len <= 32))
        {
            vecmat_dots(c, a, len, B);
            return;
        }
    }

    /* below, the classical product does not pay for its overhead; for
       very wide B (2048 columns, 4096 for full-word moduli), scalar_addmul
       streams the rows of B better than the classical product does its
       columns */
    if (len < vecmat_classical_min_len(B->mod.n)
            || ncols >= ((B->mod.norm > 0) ? 2048 : 4096))
    {
        _nmod_vec_scalar_mul_nmod(c, nmod_mat_entry_ptr(B, 0, 0), ncols, a[0], B->mod);

        for (i = 1; i < len; i++)
            _nmod_vec_scalar_addmul_nmod(c, nmod_mat_entry_ptr(B, i, 0), ncols, a[i], B->mod);

        return;
    }

    /*
        nmod_mat_mul_classical reads B by columns; for wide B a column
        spans one page per row, so take few rows at a time to keep them
        in the TLB, in blocks of equal size; for narrow B, larger blocks
        (each one has a setup cost) of up to VECMAT_NARROW_BYTES (in L1)
    */
    blk = (ncols < 16) ? VECMAT_NARROW_BYTES / (8 * ncols)
                       : (ncols < 512) ? 128 : 32;
    if (len <= blk)
        blk = len;
    else
    {
        nblk = (len + blk - 1) / blk;
        blk = (len + nblk - 1) / nblk;
    }

    TMP_START;

    t = TMP_ARRAY_ALLOC(ncols, ulong);

    A1->r = 1;
    A1->mod = B->mod;

    B1->c = ncols;
    B1->stride = B->stride;
    B1->mod = B->mod;

    C1->r = 1;
    C1->c = ncols;
    C1->stride = ncols;
    C1->mod = B->mod;

    for (i = 0; i < len; i += blk)
    {
        A1->entries = (ulong *) a + i;
        A1->c = A1->stride = FLINT_MIN(blk, len - i);

        B1->entries = nmod_mat_entry_ptr(B, i, 0);
        B1->r = A1->c;

        C1->entries = (i == 0) ? c : t;

        nmod_mat_mul_classical(C1, A1, B1);

        if (i != 0)
            _nmod_vec_add(c, c, t, ncols, B->mod);
    }

    TMP_END;
}

void nmod_mat_nmod_vec_mul_ptr(
    ulong * const * c,
    const ulong * const * a, slong alen,
    const nmod_mat_t B)
{
    slong i;
    slong len = FLINT_MIN(B->r, alen);
    slong ncols = B->c;
    ulong * aa, * cc;
    TMP_INIT;

    TMP_START;

    aa = TMP_ARRAY_ALLOC(len, ulong);
    cc = TMP_ARRAY_ALLOC(ncols, ulong);

    for (i = 0; i < len; i++)
        aa[i] = a[i][0];

    nmod_mat_nmod_vec_mul(cc, aa, len, B);

    for (i = 0; i < ncols; i++)
        c[i][0] = cc[i];

    TMP_END;
}
