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
    - otherwise the scalar code: one full modular multiplication per entry
      (nmod_mat_nmod_vec_mul), or the classical product (nmod_mat_mul).
*/

#define VM_ROWS NMOD_MAT_MUL_ROWS_MAX

/*
    Register blocks read 64-512 bytes per row of B and per block; up to
    VM_REGS_MAX_COLS columns they are used, for chunks of rows of about
    VM_L2_BYTES which the prefetching triggered by the first block leaves
    in L2 for the others (with the narrow blocks of many rows this
    prefetching stops working beyond VM_REGS_WIDE_COLS columns, where the
    rows are then taken 3 at a time); beyond, blocks of columns with
    accumulators in L1 stream the rows of B in runs of a few KB.
    FIXME measured on one machine only: need to tune these knobs for
    other machines?
*/
#ifndef VM_L2_BYTES
# define VM_L2_BYTES (1 << 20)
#endif
#ifndef VM_REGS_MAX_COLS
# define VM_REGS_MAX_COLS 4096
#endif
#ifndef VM_REGS_WIDE_COLS
# define VM_REGS_WIDE_COLS 2048
#endif

#if FLINT_BITS == 64 && NMOD_MAT_HAVE_MUL_U52

#include "nmod_mat/mul_u52_vec.h"

/*
    Rows per chunk: a lane holds a residue plus the halves of up to
    VM_U52_CHUNK products, within the bounds of u52_red2 (mul_u52_vec.h).
*/
#define VM_U52_CHUNK 511

/* the per-chunk reduction of the accumulators of a block; not inlined,
   to keep the many instantiations of vm_u52_block small */
FLINT_STATIC_NOINLINE void
vm_u52_reduce(__m512i * lo, __m512i * hi, slong cnt, int withhi,
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
vm_u52_block(ulong * const * c, __m512i * lo, __m512i * hi,
             const ulong * const * a, const nmod_mat_t B, slong j0,
             slong i0, slong i1, int r, slong nv, int withhi, int regs,
             __mmask8 m, const u52_consts * C)
{
    __m512i av[VM_ROWS], bv;
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
        i2 = FLINT_MIN(i1, i + VM_U52_CHUNK);

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
            __m512i buf[2 * VM_ROWS * 8];

            for (v = 0; v < r * nv; v++)
            {
                buf[v] = lo[v];
                if (withhi)
                    buf[r * nv + v] = hi[v];
            }
            vm_u52_reduce(buf, buf + r * nv, r * nv, withhi, C);
            for (v = 0; v < r * nv; v++)
            {
                lo[v] = buf[v];
                if (withhi)
                    hi[v] = _mm512_setzero_si512();
            }
        }
        else
            vm_u52_reduce(lo, hi, r * nv, withhi, C);
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
#ifndef VM_U52_WIDTH
# define VM_U52_WIDTH 512
#endif

#define VM_U52_WIDE_CASE(R) \
    case R: \
        if (withhi) \
            vm_u52_block(c, lo, hi, a, B, j0, 0, len, R, nv, 1, 0, m, C); \
        else \
            vm_u52_block(c, lo, hi, a, B, j0, 0, len, R, nv, 0, 0, m, C); \
        break;

static void
vm_u52_wide(ulong * const * c, const ulong * const * a, slong len,
            const nmod_mat_t B, int r, int withhi, const u52_consts * C)
{
    const slong ncols = B->c;
    /* columns per block: the accumulators of the r rows share the budget */
    const slong wmax = 8 * ((withhi ? VM_U52_WIDTH : 2 * VM_U52_WIDTH) / (8 * r));
    __m512i lo[2 * VM_U52_WIDTH / 8], hi[2 * VM_U52_WIDTH / 8];
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
            VM_U52_WIDE_CASE(1)
            VM_U52_WIDE_CASE(2)
            VM_U52_WIDE_CASE(3)
            VM_U52_WIDE_CASE(4)
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
    chunks of rows of at most VM_L2_BYTES, which the first block brings
    into L2 for the others.
*/
static const int vm_u52_regs_nv[VM_ROWS] = {8, 6, 4, 3, 2, 2, 1, 1};

#define VM_U52_REGS_CASE(R, NV) \
    case 8 * (R) + (NV): \
        if (withhi) \
            vm_u52_block(c, lo, hi, a, B, j0, i0, i1, R, NV, 1, 1, m, C); \
        else \
            vm_u52_block(c, lo, hi, a, B, j0, i0, i1, R, NV, 0, 1, m, C); \
        break;

static void
vm_u52_regs(ulong * const * c, const ulong * const * a, slong len,
            const nmod_mat_t B, int r, int withhi, const u52_consts * C)
{
    const slong ncols = B->c;
    const int nvmax = vm_u52_regs_nv[r - 1];
    const slong rows = FLINT_MAX(64, VM_L2_BYTES / (8 * ncols));
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
            __m512i lo[VM_ROWS * 8], hi[VM_ROWS * 8];
            const slong w = FLINT_MIN(ncols - j0, 8 * nvmax);
            const int nv = (int) ((w + 7) / 8);
            const __mmask8 m = (w % 8 == 0) ? 0xFF : (__mmask8) ((1 << (w % 8)) - 1);

            switch (8 * r + nv)
            {
                VM_U52_REGS_CASE(1, 1)
                VM_U52_REGS_CASE(1, 2)
                VM_U52_REGS_CASE(1, 3)
                VM_U52_REGS_CASE(1, 4)
                VM_U52_REGS_CASE(1, 5)
                VM_U52_REGS_CASE(1, 6)
                VM_U52_REGS_CASE(1, 7)
                VM_U52_REGS_CASE(1, 8)
                VM_U52_REGS_CASE(2, 1)
                VM_U52_REGS_CASE(2, 2)
                VM_U52_REGS_CASE(2, 3)
                VM_U52_REGS_CASE(2, 4)
                VM_U52_REGS_CASE(2, 5)
                VM_U52_REGS_CASE(2, 6)
                VM_U52_REGS_CASE(3, 1)
                VM_U52_REGS_CASE(3, 2)
                VM_U52_REGS_CASE(3, 3)
                VM_U52_REGS_CASE(3, 4)
                VM_U52_REGS_CASE(4, 1)
                VM_U52_REGS_CASE(4, 2)
                VM_U52_REGS_CASE(4, 3)
                VM_U52_REGS_CASE(5, 1)
                VM_U52_REGS_CASE(5, 2)
                VM_U52_REGS_CASE(6, 1)
                VM_U52_REGS_CASE(6, 2)
                VM_U52_REGS_CASE(7, 1)
                VM_U52_REGS_CASE(8, 1)
                default: FLINT_UNREACHABLE;
            }

            j0 += w;
        }
    }
}

#endif  /* u52 */

#if FLINT_BITS == 64 && NMOD_MAT_HAVE_VM_U32

#include <immintrin.h>

/*
    Moduli up to 2^32 without IFMA: vpmuludq products below 2^64, split
    at bit 56 (nmod_vec/dot.c): lo accumulates the low 56 bits of the
    products, hi their high 8 bits, and lo is folded into hi every
    VM_U32_FOLD rows. The total lo + 2^56 hi fits two limbs for any
    length. The reduction to the residue is scalar, so the chunks of rows
    between two of them are long (VM_U32_ROWS at least).
*/

#if defined(__AVX512F__)
# define VMU_VL 8
typedef __m512i vmu_t;
# define vmu_zero()        _mm512_setzero_si512()
# define vmu_set1(x)       _mm512_set1_epi64(x)
# define vmu_loadu(p)      _mm512_loadu_si512((const void *) (p))
# define vmu_storeu(p, x)  _mm512_storeu_si512((void *) (p), x)
# define vmu_mul(x, y)     _mm512_mul_epu32(x, y)
# define vmu_add(x, y)     _mm512_add_epi64(x, y)
# define vmu_and(x, y)     _mm512_and_si512(x, y)
# define vmu_srli56(x)     _mm512_srli_epi64(x, 56)
/* vectors of columns per register block for r rows (2 r nv + r + 3 of
   32 registers), powers of two */
static const int vm_u32_regs_nv[4] = {8, 4, 4, 2};
#else
# define VMU_VL 4
typedef __m256i vmu_t;
# define vmu_zero()        _mm256_setzero_si256()
# define vmu_set1(x)       _mm256_set1_epi64x(x)
# define vmu_loadu(p)      _mm256_loadu_si256((const __m256i *) (p))
# define vmu_storeu(p, x)  _mm256_storeu_si256((__m256i *) (p), x)
# define vmu_mul(x, y)     _mm256_mul_epu32(x, y)
# define vmu_add(x, y)     _mm256_add_epi64(x, y)
# define vmu_and(x, y)     _mm256_and_si256(x, y)
# define vmu_srli56(x)     _mm256_srli_epi64(x, 56)
/* 16 registers */
static const int vm_u32_regs_nv[4] = {4, 2, 1, 1};
#endif

#define VM_U32_FOLD 255
#define VM_U32_WIDTH 512
#define VM_U32_ROWS 255
/* the register blocks are narrow (128-512 bytes per row of B): beyond
   this many columns the L1 blocks are used */
#ifndef VM_U32_REGS_MAX_COLS
# define VM_U32_REGS_MAX_COLS (128 * VMU_VL)
#endif

/*
    One block of nv vectors of columns j0 <= j < j0 + w for r rows and the
    rows i0 <= i < i1 of B; lo[s * nv + v] starts from the residues in c
    (zero at first) and the residues of the totals are stored back. With
    regs = 1, r and nv are constants and the arrays register allocated;
    otherwise they are the caller's (in L1) and nv is arbitrary. A last
    partial vector goes through a zero-padded buffer.
*/
FLINT_FORCE_INLINE void
vm_u32_block(ulong * const * c, vmu_t * lo, vmu_t * hi,
             const ulong * const * a, const nmod_mat_t B, slong j0,
             slong w, slong i0, slong i1, int r, slong nv, int regs)
{
    const nmod_t mod = B->mod;
    const slong wl = w - VMU_VL * (nv - 1);   /* columns of the last vector */
    const vmu_t mask = vmu_set1((UWORD(1) << 56) - 1);
    vmu_t av[VM_ROWS], bv, p;
    ulong buf[VMU_VL];
    slong i, j, v;
    int s, cnt = 0;

    for (s = 0; s < r; s++)
    {
        for (v = 0; v < nv - 1; v++)
            lo[s * nv + v] = vmu_loadu(c[s] + j0 + VMU_VL * v);
        for (j = 0; j < VMU_VL; j++)
            buf[j] = (j < wl) ? c[s][j0 + VMU_VL * v + j] : 0;
        lo[s * nv + v] = vmu_loadu(buf);
    }
    for (v = 0; v < r * nv; v++)
        hi[v] = vmu_zero();

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
                hi[s * nv + v] = vmu_add(hi[s * nv + v], vmu_srli56(p));
            }
        }
        if (wl == VMU_VL)
            bv = vmu_loadu(b + VMU_VL * v);
        else
        {
            for (j = 0; j < VMU_VL; j++)
                buf[j] = (j < wl) ? b[VMU_VL * v + j] : 0;
            bv = vmu_loadu(buf);
        }
        for (s = 0; s < r; s++)
        {
            p = vmu_mul(av[s], bv);
            lo[s * nv + v] = vmu_add(lo[s * nv + v], vmu_and(p, mask));
            hi[s * nv + v] = vmu_add(hi[s * nv + v], vmu_srli56(p));
        }

        if (++cnt == VM_U32_FOLD)
        {
            for (v = 0; v < r * nv; v++)
            {
                hi[v] = vmu_add(hi[v], vmu_srli56(lo[v]));
                lo[v] = vmu_and(lo[v], mask);
            }
            cnt = 0;
        }
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

#define VM_U32_REGS_CASE(R, NV) \
    case 8 * (R) + (NV): \
        vm_u32_block(c, lo, hi, a, B, j0, w, i0, i1, R, NV, 1); break;

static void
vm_u32_regs(ulong * const * c, const ulong * const * a, slong len,
            const nmod_mat_t B, int r)
{
    const slong ncols = B->c;
    const int nvmax = vm_u32_regs_nv[r - 1];
    const slong rows = FLINT_MAX(VM_U32_ROWS, VM_L2_BYTES / (8 * ncols));
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
            vmu_t lo[VM_ROWS * 8], hi[VM_ROWS * 8];
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
                VM_U32_REGS_CASE(1, 1)
                VM_U32_REGS_CASE(1, 2)
                VM_U32_REGS_CASE(1, 4)
                VM_U32_REGS_CASE(2, 1)
                VM_U32_REGS_CASE(2, 2)
                VM_U32_REGS_CASE(3, 1)
                VM_U32_REGS_CASE(4, 1)
#if VMU_VL == 8
                VM_U32_REGS_CASE(1, 8)
                VM_U32_REGS_CASE(2, 4)
                VM_U32_REGS_CASE(3, 2)
                VM_U32_REGS_CASE(3, 4)
                VM_U32_REGS_CASE(4, 2)
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
vm_u32_wide(ulong * const * c, const ulong * const * a, slong len,
            const nmod_mat_t B, int r)
{
    const slong ncols = B->c;
    const slong wmax = VMU_VL * (VM_U32_WIDTH / (VMU_VL * r));
    vmu_t lo[VM_U32_WIDTH / VMU_VL], hi[VM_U32_WIDTH / VMU_VL];
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
            case 1: vm_u32_block(c, lo, hi, a, B, j0, w, 0, len, 1, nv, 0); break;
            case 2: vm_u32_block(c, lo, hi, a, B, j0, w, 0, len, 2, nv, 0); break;
            case 3: vm_u32_block(c, lo, hi, a, B, j0, w, 0, len, 3, nv, 0); break;
            case 4: vm_u32_block(c, lo, hi, a, B, j0, w, 0, len, 4, nv, 0); break;
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
#define VM_HAVE_FP50 1

/* balanced representative of 0 <= x < n as a double */
FLINT_FORCE_INLINE double
vm_sym(ulong x, ulong n)
{
    return (x > n / 2) ? (double) (slong) (x - n) : (double) x;
}

/* columns per block for one row: the accumulators (8 bytes per column
   and row) stay in L1 */
#define VM_FP50_WIDTH 512

/*
    One block of nv vectors of columns j0 <= j < j0 + w for r rows, with
    the accumulators acc[s * nv + v] in L1 (r is a constant at each call
    site, nv is not; register accumulators were tried and bring nothing:
    on AVX2 this loop is bound by the issue width, at about 10 vector
    operations per vector of products). Each term is in (-9/8 n, 9/8 n)
    (fpv_mulmod with |a| <= n/2), and the accumulators are reduced every
    four rows, staying below 5.01 n < 2^53 in absolute value. A last
    partial vector goes through a zero-padded buffer.
*/
FLINT_FORCE_INLINE void
vm_fp50_block(ulong * const * c, fpv * acc, const ulong * const * a,
              const nmod_mat_t B, slong j0, slong w, slong len, int r,
              slong nv, fpv nvec, fpv ninv)
{
    const ulong n = B->mod.n;
    const slong wl = w - FPV_VL * (nv - 1);   /* columns of the last vector */
    fpv av[VM_ROWS], bv;
    ulong buf[FPV_VL];
    slong i, j, v;
    int s, cnt = 0;

    for (v = 0; v < r * nv; v++)
        acc[v] = fpv_zero();

    for (i = 0; i < len; i++)
    {
        const ulong * b = nmod_mat_entry_ptr(B, i, j0);

        for (s = 0; s < r; s++)
            av[s] = fpv_set1(vm_sym(a[s][i], n));

        for (v = 0; v < nv - 1; v++)
        {
            bv = fpv_load_u64(b + FPV_VL * v);
            for (s = 0; s < r; s++)
                acc[s * nv + v] = fpv_add(acc[s * nv + v],
                                          fpv_mulmod(av[s], bv, nvec, ninv));
        }
        if (wl == FPV_VL)
            bv = fpv_load_u64(b + FPV_VL * v);
        else
        {
            for (j = 0; j < FPV_VL; j++)
                buf[j] = (j < wl) ? b[FPV_VL * v + j] : 0;
            bv = fpv_load_u64(buf);
        }
        for (s = 0; s < r; s++)
            acc[s * nv + v] = fpv_add(acc[s * nv + v],
                                      fpv_mulmod(av[s], bv, nvec, ninv));

        if (++cnt == 4)
        {
            for (v = 0; v < r * nv; v++)
                acc[v] = fpv_reduce_pm1n(acc[v], nvec, ninv);
            cnt = 0;
        }
    }

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
vm_fp50(ulong * const * c, const ulong * const * a, slong len,
        const nmod_mat_t B, int r)
{
    const slong ncols = B->c;
    const fpv nvec = fpv_set1((double) B->mod.n);
    const fpv ninv = fpv_set1(1.0 / (double) B->mod.n);
    const slong wmax = FPV_VL * (VM_FP50_WIDTH / (FPV_VL * r));
    fpv acc[VM_FP50_WIDTH / FPV_VL];
    slong j0 = 0;

    while (j0 < ncols)
    {
        const slong w = FLINT_MIN(ncols - j0, wmax);
        const slong nv = (w + FPV_VL - 1) / FPV_VL;

        switch (r)
        {
            case 1: vm_fp50_block(c, acc, a, B, j0, w, len, 1, nv, nvec, ninv); break;
            case 2: vm_fp50_block(c, acc, a, B, j0, w, len, 2, nv, nvec, ninv); break;
            case 3: vm_fp50_block(c, acc, a, B, j0, w, len, 3, nv, nvec, ninv); break;
            case 4: vm_fp50_block(c, acc, a, B, j0, w, len, 4, nv, nvec, ninv); break;
            default: FLINT_UNREACHABLE;
        }

        j0 += w;
    }
}

#else
#define VM_HAVE_FP50 0
#endif

#endif  /* FLINT_BITS == 64 */

int
_nmod_mat_mul_rows_simd(ulong * const * c, const ulong * const * a, slong r,
                        slong len, const nmod_mat_t B)
{
    FLINT_ASSERT(1 <= r && r <= NMOD_MAT_MUL_ROWS_MAX);
    FLINT_ASSERT(len >= 1 && len <= B->r && B->c >= 1);

#if FLINT_BITS == 64
#if NMOD_MAT_HAVE_MUL_U52
    if (B->mod.n <= (UWORD(1) << 52))
    {
        u52_ctx_struct ctx;
        u52_ctx_init(&ctx, B->mod.n);
        const u52_consts C = u52_consts_init(&ctx);

        slong s, g;

        if (B->c <= VM_REGS_MAX_COLS)
        {
            /* with many columns and the high halves, groups of 3 rows keep
               the register blocks 4 vectors wide (see vm_u52_regs) */
            g = (!ctx.lo_only && B->c > VM_REGS_WIDE_COLS) ? 3 : VM_ROWS;
            for (s = 0; s < r; s += g)
                vm_u52_regs(c + s, a + s, len, B, (int) FLINT_MIN(r - s, g),
                            !ctx.lo_only, &C);
        }
        else
        {
            for (s = 0; s < r; s += 4)
                vm_u52_wide(c + s, a + s, len, B, (int) FLINT_MIN(r - s, 4),
                            !ctx.lo_only, &C);
        }
        return 1;
    }
#endif
#if NMOD_MAT_HAVE_VM_U32
    if (B->mod.n <= (UWORD(1) << 32))
    {
        /* groups of at most 4 rows */
        slong s;
        for (s = 0; s < r; s += 4)
        {
            const int rs = (int) FLINT_MIN(r - s, 4);
            if (B->c <= VM_U32_REGS_MAX_COLS)
                vm_u32_regs(c + s, a + s, len, B, rs);
            else
                vm_u32_wide(c + s, a + s, len, B, rs);
        }
        return 1;
    }
#endif
#if VM_HAVE_FP50
    if (B->mod.n < (UWORD(1) << 50))
    {
        slong s;
        for (s = 0; s < r; s += 4)
            vm_fp50(c + s, a + s, len, B, (int) FLINT_MIN(r - s, 4));
        return 1;
    }
#endif
#endif

    return 0;
}

void nmod_mat_nmod_vec_mul(
    ulong * c,
    const ulong * a, slong alen,
    const nmod_mat_t B)
{
    slong i;
    slong len = FLINT_MIN(B->r, alen);
    slong ncols = B->c;

    /* scalar_addmul wants non-empty */
    if (ncols < 1)
        return;

    if (len <= 0)
    {
        _nmod_vec_zero(c, ncols);
        return;
    }

    if (len >= 2 && _nmod_mat_mul_rows_simd(&c, &a, 1, len, B))
        return;

    _nmod_vec_scalar_mul_nmod(c, nmod_mat_entry_ptr(B, 0, 0), ncols, a[0], B->mod);

    for (i = 1; i < len; i++)
        _nmod_vec_scalar_addmul_nmod(c, nmod_mat_entry_ptr(B, i, 0), ncols, a[i], B->mod);
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
