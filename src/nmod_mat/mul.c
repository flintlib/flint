/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2020 William Hart
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "thread_support.h"

#include "longlong.h"
#include "flint-mparam.h"
#include "nmod_mat/impl.h"

#if FLINT_BITS == 64

/*
    Dimension from which nmod_mat_mul_blas is preferred when one dgemm pass
    suffices (0: never), for the current thread count: BLAS_1PASS_CUTOFF
    measured with 1 thread, BLAS_1PASS_CUTOFF_MT with 4, their geometric
    mean for 2 threads, and the latter for 3 threads or more.
*/
static slong
_blas_1pass_cutoff(slong num_threads)
{
    slong c1 = FLINT_NMOD_MAT_MUL_BLAS_1PASS_CUTOFF;
    slong c4 = FLINT_NMOD_MAT_MUL_BLAS_1PASS_CUTOFF_MT;

    if (num_threads <= 1)
        return c1;
    if (num_threads >= 3 || c1 <= 0 || c4 <= 0)
        return c4;
    return (slong) n_sqrt((ulong) c1 * (ulong) c4);
}

#endif

/*
    C = A * B with all dimensions at most tmax (a constant at each call
    site, <= 8): the bound k (n - 1)^2 on the dot products tells whether
    they fit one word, two words with the high word below n (one reduction
    by NMOD_RED2), or three words. The bound with k = tmax decides for most
    moduli; the bound with the actual k only for the others.
*/
FLINT_FORCE_INLINE void
_nmod_mat_mul_tiny_impl(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B,
                        slong tmax)
{
    const slong m = A->r, k = A->c, n = B->c;
    const nmod_t mod = C->mod;
    ulong t[64], bh, bl, kh, kl, s0, s1, s2, h, l;
    slong i, j, q;
    int words;

    FLINT_ASSERT(m <= tmax && k <= tmax && n <= tmax && tmax <= 8);

    /* (bh, bl) = (n - 1)^2 */
    umul_ppmm(bh, bl, mod.n - 1, mod.n - 1);

    if (bh == 0 && bl <= UWORD_MAX / tmax)
        words = 1;
    else if (bh != 0 && bh + 1 <= mod.n / tmax)
        words = 2;
    else
    {
        /* one word if bh = 0 and k bl does not overflow; otherwise two
           words if bh = 0 (then k < n) or k (bh + 1) <= n */
        umul_ppmm(kh, kl, (bh == 0) ? bl : bh + 1, (ulong) k);
        words = (bh == 0) ? ((kh == 0) ? 1 : 2)
                          : ((kh == 0 && kl <= mod.n) ? 2 : 3);
    }

    if (words == 1)
    {
        for (i = 0; i < m; i++)
        {
            const ulong * a = nmod_mat_entry_ptr(A, i, 0);
            for (j = 0; j < n; j++)
            {
                s0 = 0;
                for (q = 0; q < k; q++)
                    s0 += a[q] * nmod_mat_entry(B, q, j);
                NMOD_RED(t[i * n + j], s0, mod);
            }
        }
    }
    else if (words == 2)
    {
        for (i = 0; i < m; i++)
        {
            const ulong * a = nmod_mat_entry_ptr(A, i, 0);
            for (j = 0; j < n; j++)
            {
                s1 = s0 = 0;
                for (q = 0; q < k; q++)
                {
                    umul_ppmm(h, l, a[q], nmod_mat_entry(B, q, j));
                    add_ssaaaa(s1, s0, s1, s0, h, l);
                }
                NMOD_RED2(t[i * n + j], s1, s0, mod);
            }
        }
    }
    else
    {
        for (i = 0; i < m; i++)
        {
            const ulong * a = nmod_mat_entry_ptr(A, i, 0);
            for (j = 0; j < n; j++)
            {
                s2 = s1 = s0 = 0;
                for (q = 0; q < k; q++)
                {
                    umul_ppmm(h, l, a[q], nmod_mat_entry(B, q, j));
                    add_sssaaaaaa(s2, s1, s0, s2, s1, s0, UWORD(0), h, l);
                }
                NMOD_RED3(t[i * n + j], s2, s1, s0, mod);
            }
        }
    }

    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            nmod_mat_entry(C, i, j) = t[i * n + j];
}

/* C = A * B for A of size 1 x 1 (not inlined, to keep the code of
   nmod_mat_mul before it small) */
FLINT_STATIC_NOINLINE void
_nmod_mat_mul_scalar_row(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    _nmod_vec_scalar_mul_nmod(nmod_mat_entry_ptr(C, 0, 0),
                              nmod_mat_entry_ptr(B, 0, 0), B->c,
                              nmod_mat_entry(A, 0, 0), C->mod);
}

void
_nmod_mat_mul_tiny(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    _nmod_mat_mul_tiny_impl(C, A, B, 8);
}

/*
    Whether nmod_mat_mul tries the rows path for m = 2..NMOD_MAT_MUL_ROWS_MAX
    rows, k >= 2 and n columns, when the modulus modn has one (on top of the
    conditions in nmod_mat_mul). Measured on Zen 4 and Ice Lake SP
    (AVX512-IFMA), Cascade Lake (AVX-512 without IFMA) and Apple M4, 20 to
    64 bits, m = 2, 3, 4, 8, k = 2..512, n = 2..32, against what nmod_mat_mul
    does otherwise (the few-columns path, the SIMD kernels, or the classical
    code).

    TODO check with measurements on an AVX2 (non-AVX-512) machine
*/
static int
_nmod_mat_mul_rows_wanted(slong m, slong k, slong n, ulong modn,
                          slong num_threads)
{
    int kernel_shape = (n >= FLINT_NMOD_MAT_MUL_SIMD_MIN_DIM
        && (m >= FLINT_NMOD_MAT_MUL_SIMD_MIN_DIM
            || (2 * m >= FLINT_NMOD_MAT_MUL_SIMD_MIN_DIM
                && n >= 4 * FLINT_NMOD_MAT_MUL_SIMD_MIN_DIM)));
    int x86_no_ifma = NMOD_MAT_HAVE_VECMAT_U64 && !NMOD_MAT_HAVE_MUL_U52;

    (void) kernel_shape;

    if (n >= 2 && n <= NMOD_MAT_MUL_COLS_MAX && num_threads == 1
            && NMOD_MAT_MUL_COLS_IS_SIMD(modn)
            && ((n == 2 && k >= (x86_no_ifma ? 32 : 128))
                || (x86_no_ifma && n <= 4 && k >= 128)))
        return 0;

    if (NMOD_MAT_NMOD_VEC_MUL_IS_U64(modn))
    {
#if FLINT_BITS == 64
        if (modn <= (UWORD(1) << 61))
            return n >= NMOD_MAT_MUL_ROWS_U64_MIN_COLS || (n >= 4 && k <= 128);
        else
            return n >= NMOD_MAT_MUL_ROWS_U64_MIN_COLS
                && !(m >= 4 && k >= 256 && n < 32 && n % 8 != 0);
#else
        return n >= NMOD_MAT_MUL_ROWS_U64_MIN_COLS;
#endif
    }

#if FLINT_BITS == 64
    /* for small k, the rows path takes its fp50 tier (see
       NMOD_MAT_VECMAT_U32_TO_FP50), which also beats the kernels */
    if (NMOD_MAT_HAVE_VECMAT_U32 && modn <= (UWORD(1) << 32))
        return n >= 4 && (NMOD_MAT_VECMAT_U32_TO_FP50(k, n)
                          || !(k <= 4 && kernel_shape));

    /* for k <= 2, the kernels where nmod_mat_mul reaches them (more than
       NMOD_MAT_MUL_COLS_MAX columns, otherwise the few-columns path comes
       first) and where they are fp50 (above 32 bits) */
    if (NMOD_MAT_HAVE_MUL_U52 && modn <= (UWORD(1) << 52))
    {
        flint_bitcnt_t bits = FLINT_BIT_COUNT(modn);
        return !(k <= 2 && kernel_shape && n > NMOD_MAT_MUL_COLS_MAX
                 && bits > 32 && bits <= FLINT_NMOD_MAT_MUL_FP50_MAX_BITS);
    }
#endif

    return 1;
}

void
nmod_mat_mul(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    slong m = A->r;
    slong k = A->c;
    slong n = B->c;
    slong min_dim = FLINT_MIN(FLINT_MIN(m, k), n);
    slong cutoff;

    FLINT_ASSERT(C->r == A->r);
    FLINT_ASSERT(C->c == B->c);
    FLINT_ASSERT(A->c == B->r);

    /* tiny matrices: before the costs of the dispatch below and of the
       setup of the routines it calls */
    if (m <= NMOD_MAT_MUL_TINY_MAX && k <= NMOD_MAT_MUL_TINY_MAX
            && n <= NMOD_MAT_MUL_TINY_MAX)
    {
        _nmod_mat_mul_tiny_impl(C, A, B, NMOD_MAT_MUL_TINY_MAX);
        return;
    }

    /* a scalar times a row of more than NMOD_MAT_MUL_TINY_MAX entries */
    if (m == 1 && k == 1)
    {
        _nmod_mat_mul_scalar_row(C, A, B);
        return;
    }

    /* Todo: optimize nmod_mat_mul_blas. For mod.n >= 17, nmod_mat_mul_blas
       should be faster up to dim about 1000-2000 as nmod_mat_mul_u8 does
       the same sgemm but with an extra uint8 roundtrip. Currently
       mul_blas narrowly loses to mul_u8 due to slow modular reductions. */
    if (C->mod.n <= 255 && min_dim >= 8)
    {
        nmod_mat_mul_u8(C, A, B);
        return;
    }

    slong flint_num_threads = flint_get_num_threads();

    /*
        B is a single column: use nmod_mat_mul_nmod_vec
        TODO handle multithreading in mul_nmod_vec
    */
    if (n == 1 && m >= 8 && k >= 1 && flint_num_threads == 1)
    {
        ulong * t;
        slong i;
        TMP_INIT;

        TMP_START;
        t = TMP_ARRAY_ALLOC(k + m, ulong);
        for (i = 0; i < k; i++)
            t[i] = nmod_mat_entry(B, i, 0);
        nmod_mat_mul_nmod_vec(t + k, A, t, k);
        for (i = 0; i < m; i++)
            nmod_mat_entry(C, i, 0) = t[k + i];
        TMP_END;
        return;
    }

#if FLINT_BITS == 64 && NMOD_MAT_HAVE_FPV && !NMOD_MAT_HAVE_VECMAT_U64 \
        && !NMOD_MAT_HAVE_MUL_U52
    /*
        A has a few rows, NEON (the rows path has only its fp50 tier, two
        lanes per vector): measured on Apple M4 for N < 2^50, m = 2..4,
        k = 2..512, n = 2..32. For N <= 2^32 the u32 kernel where it is
        faster (larger k and n; only from 4 rows for 32-bit N, where its
        folds come every few products), otherwise the classical code where
        it is faster (small n relative to k), otherwise the rows path, also
        for the shapes the generic condition below would send to the
        kernels (m = 4, n >= 32, small k). The few-columns path never won.
    */
    if (m >= 2 && m <= NMOD_MAT_MUL_ROWS_MAX && k >= 2
            && flint_num_threads == 1 && C != A && C != B
            && C->mod.n < (UWORD(1) << 50))
    {
        ulong * c[NMOD_MAT_MUL_ROWS_MAX];
        const ulong * a[NMOD_MAT_MUL_ROWS_MAX];

        if (C->mod.n <= (UWORD(1) << 32)
                && ((m >= 4 && ((k >= 16 && n >= (k >= 512 ? 6 : 4))
                                || (k >= 8 && n >= 16)))
                    || (C->mod.n <= (UWORD(1) << 31)
                        && ((m == 3 && ((k >= 32 && n >= 8)
                                        || (k >= 16 && n >= 20)))
                            || (m == 2 && k >= 128 && n >= 24))))
                && nmod_mat_mul_u32(C, A, B))
            return;

        if ((n <= 12 && n * (m == 2 ? 3 : 5) <= k)
                || (m == 2 && k >= 128 && n <= 20)
                || (k >= 512 && n <= 24))
        {
            nmod_mat_mul_classical(C, A, B);
            return;
        }

        for (slong i = 0; i < m; i++)
        {
            c[i] = nmod_mat_entry_ptr(C, i, 0);
            a[i] = nmod_mat_entry_ptr(A, i, 0);
        }

        if (_nmod_mat_mul_rows_simd(c, a, m, k, B))
            return;
    }
#endif

    /*
        A has a few rows (up to NMOD_MAT_MUL_ROWS_MAX): _nmod_mat_mul_rows_simd
        streams B once for all of them with delayed reductions where it is vectorized
        This avoids calling the kernels which pad these rows to a whole tile
        of MR rows (MR depends on each kernel and on the ISA, but it can be above
        10).

        Where the rows engine only has its floating point tier for moduli
        below 2^32 (no AVX512-IFMA and no x86 integer tier: NEON), the u32
        kernel is faster whenever it takes the shape (Apple M4: 2.0-2.5x
        at 4 rows), see NMOD_MAT_MUL_ROWS_PREFER_KERNEL.

        Where the few-columns path, the kernels or the classical code are
        faster for these shapes, see _nmod_mat_mul_rows_wanted (one row:
        nmod_mat_nmod_vec_mul, which has its own rule).

        NOTE It is single-threaded: with several threads it is still used
        for the shapes the threaded kernels do not take (below), which
        would otherwise go to nmod_mat_mul_classical_threaded (observed 2-7x
        slower than it on Ice Lake, Meteor Lake, 4 threads), and for the
        moduli without SIMD kernel (above 2^52), for which the classical
        code is serial with so few rows.
        TODO thread it (split the columns of B).
    */
    if (m >= 1 && m <= NMOD_MAT_MUL_ROWS_MAX && k >= 2 && n >= 1
            && ((flint_num_threads == 1
                 && !NMOD_MAT_MUL_ROWS_PREFER_KERNEL(C->mod.n))
                || !(NMOD_MAT_MUL_HAVE_SIMD_KERNEL(C->mod.n)
                     && ((m >= FLINT_NMOD_MAT_MUL_SIMD_MIN_DIM
                          && n >= FLINT_NMOD_MAT_MUL_SIMD_MIN_DIM)
                         || (2 * m >= FLINT_NMOD_MAT_MUL_SIMD_MIN_DIM
                             && n >= 4 * FLINT_NMOD_MAT_MUL_SIMD_MIN_DIM))))
            && C != A && C != B && NMOD_MAT_NMOD_VEC_MUL_IS_SIMD(C->mod.n)
            && (m == 1 || _nmod_mat_mul_rows_wanted(m, k, n, C->mod.n,
                                                    flint_num_threads)))
    {
        ulong * c[NMOD_MAT_MUL_ROWS_MAX];
        const ulong * a[NMOD_MAT_MUL_ROWS_MAX];

        for (slong i = 0; i < m; i++)
        {
            c[i] = nmod_mat_entry_ptr(C, i, 0);
            a[i] = nmod_mat_entry_ptr(A, i, 0);
        }

        /* one row: the rule of nmod_mat_nmod_vec_mul for narrow B */
        if (m == 1)
        {
            nmod_mat_nmod_vec_mul(c[0], a[0], k, B);
            return;
        }

        if (_nmod_mat_mul_rows_simd(c, a, m, k, B))
            return;
    }

    /*
        B has a few columns (up to NMOD_MAT_MUL_COLS_MAX):
        _nmod_mat_mul_cols_simd multiplies each row of A with all of them at
        once, as dot products sharing the loads of the row, on the transposed
        columns. This avoids calling the kernels which would pad them to
        NR = 8-24 columns; this also avoids reading B column by column as in
        the classical code.
    */
    if (n >= 2 && n <= NMOD_MAT_MUL_COLS_MAX && k >= 1
            && flint_num_threads == 1 && C != A && C != B
            && NMOD_MAT_MUL_COLS_IS_SIMD(C->mod.n))
    {
        if (_nmod_mat_mul_cols_simd(C, A, B))
            return;
    }

    /*
        Moduli up to 2^52: SIMD kernels with delayed reduction,
        nmod_mat_mul_u32 (any 64-bit target, moduli below 2^32),
        nmod_mat_mul_u52 (AVX512-IFMA, moduli up to 2^52) and, without
        IFMA, nmod_mat_mul_k52 (two-limb integer Karatsuba, moduli up to
        2^52) or nmod_mat_mul_fp50 (all in double precision with the
        mulmod of fft_small, moduli below 2^50). The parameters come from
        flint-mparam.h and were measured with
        src/nmod_mat/profile/p-mul_tune.c; the picture on the machines
        measured so far (Cascade/Ice/Meteor Lake, Zen 4, Apple M4) is:

        - Where nmod_mat_mul_blas needs several dgemm passes and a CRT
          (k*(n/2)^2 >= 2^53, always the case from 25 bits on), the SIMD
          kernels are 2-5x faster than any other method from dimension 8.
        - Where one dgemm pass suffices (up to about 23 bits), u32 wins
          below a dimension of about 100-250 (mul_blas pays O(n^2)
          conversions and thread hand-offs around its gemm) and is
          0.7-1.1x of mul_blas above. The single-IFMA mode of u52 (moduli
          up to 2^26) is 1.1-1.35x faster than u32 and than FLINT's own gemm,
          but an external BLAS can overtake it at large dimensions. Hence
          BLAS_1PASS_CUTOFF(_MT), consulted for u32 always and for u52 only
          with an external BLAS.
        - u52 removes the 31-32 bit cliff of u32 (whose in-kernel folds
          then come every 1-2 products) and extends the single pass to
          52 bits, where the alternative is 4-5 dgemm passes. Its two-IFMA
          mode is slower than u32 between 27 and 30-31 bits (1.5-1.65x on
          Zen 4, up to 1.1x on Ice Lake), hence U52_MIN_BITS.
        - Thin shapes: what the kernels pay for is the padding of the
          last tile of C, of MR rows (4-14) and NR columns (8-24), and
          the packing of A and B, whereas the inner dimension k has no
          such cost. So the kernels are used for any k >= 1 as soon as
          n >= SIMD_MIN_DIM, and m >= SIMD_MIN_DIM or m >= SIMD_MIN_DIM/2
          with n >= 4*SIMD_MIN_DIM (1.1-10x faster than the classical
          code on these shapes on AVX-512, for instance 1000 x 4 times
          4 x 1000 or 4 x 1000 times 1000 x 1000). Below that, notably
          for n < 8 where the padding to NR columns wastes most of the
          tile, the classical code wins.
        - Small inner dimension k on IFMA machines: the per-entry
          reduction of u52 (Barrett through doubles and a Shoup step)
          weighs more than its products, so up to k = 8 u32 (<= 32
          bits) and fp50 (33-50 bits) are 1.15-2.3x faster on Ice Lake
          (1000 x k x 1000); u52 wins from k = 16 on. Hence the
          parameter FLINT_NMOD_MAT_MUL_U52_MIN_K of flint-mparam.h
          TODO measured on one machine so far, could be useful to
          analyze it more thoroughly
        - Tiny inner dimension, moduli up to 2^32: fp50 is 1.1-2x faster
          than u32 at k = 1..3 on AVX2; with AVX-512 and on NEON the range
          depends on the machine (FLINT_NMOD_MAT_MUL_FP50_U32_MIN_K / _MAX_K,
          see impl.h).
        - Single-threaded, Strassen on top (its recursive calls come back
          here) pays from SIMD_STRASSEN_CUTOFF on (e.g. 320-512 on Zen 4).
          TODO with 4 threads it did not pay up to 4096 (which is not that
          large) on several machine, so it is not used with several threads,
          but of course one may handle multi-threading differently within
          Strassen: this requires investigation.
        - Without IFMA, k52 and fp50 are the single-pass options from 33
          bits on, and which of the two wins is a property of the
          instruction set rather than of the modulus (FP50_MAX_BITS).
          On x86 fp50 wins everywhere measured (1.0-2.1x on Cascade Lake
          and Meteor Lake): it keeps one accumulator per tile cell where
          k52 needs three, hence a tile 2.7x wider and fewer operand
          loads per product. On NEON the single-instruction widening
          multiply-add (smlal) reverses this and k52 wins from dimension
          48 on (1.4x on Apple M4). Above 2^50, where fp50 stops, k52 is
          the only single-pass option.
        - These two are 1.4-3.4x faster than blas + CRT up to a
          dimension that grows with the modulus size, and lose beyond it
          (K52_BLAS_CUTOFF, only with an external BLAS): on Apple M4
          that dimension is 320-448 with Accelerate, whose dgemm is far
          out of reach of a NEON kernel; on the x86 machines measured it
          is 768 and beyond, or never. With FLINT's own gemm, blas + CRT
          was 1.3-1.8x slower than k52 on Apple M4 at every dimension.
    */
#if FLINT_BITS == 64
    if (C->mod.n <= (UWORD(1) << 52) && k >= 1
            && n >= FLINT_NMOD_MAT_MUL_SIMD_MIN_DIM
            && (m >= FLINT_NMOD_MAT_MUL_SIMD_MIN_DIM
                || (2 * m >= FLINT_NMOD_MAT_MUL_SIMD_MIN_DIM
                    && n >= 4 * FLINT_NMOD_MAT_MUL_SIMD_MIN_DIM)))
    {
        flint_bitcnt_t bits = FLINT_BIT_COUNT(C->mod.n);
        int (* simd_mul)(nmod_mat_t, const nmod_mat_t, const nmod_mat_t);
        int blas_1pass = 0;

        simd_mul = NULL;

        if (bits <= 32)
        {
            /*
                mul_blas is preferred from BLAS_1PASS_CUTOFF(_MT) on
                (0: never) when it can do it in one dgemm pass,
                k*(n/2)^2 < 2^53. Here half < 2^31, so half^2 does not
                overflow.
            */
            ulong half = C->mod.n / 2;
            slong cut = _blas_1pass_cutoff(flint_num_threads);

            /* min_dim > 100: the BLAS branch below requires it,
               and u32 must not be skipped without BLAS taking over */
            blas_1pass = cut > 0 && min_dim >= cut && min_dim > 100
                && (half == 0
                    || (ulong) k <= ((UWORD(1) << 53) - 1) / (half * half));
        }

        if (NMOD_MAT_HAVE_MUL_U52
                && (bits >= FLINT_NMOD_MAT_MUL_U52_MIN_BITS
                    || bits <= FLINT_NMOD_MAT_MUL_U52_LO_MAX_BITS)
                && (k >= FLINT_NMOD_MAT_MUL_U52_MIN_K
                    || bits > FLINT_NMOD_MAT_MUL_FP50_MAX_BITS
                    || (bits > 32 && !NMOD_MAT_HAVE_FPV)))
        {
            /* FLINT's own gemm was never measured faster than u52 */
#if FLINT_USES_BLAS
            if (!blas_1pass)
#endif
                simd_mul = nmod_mat_mul_u52;
        }
        else if (bits <= 32)
        {
            if (NMOD_MAT_HAVE_FPV
                    && k >= FLINT_NMOD_MAT_MUL_FP50_U32_MIN_K
                    && k <= FLINT_NMOD_MAT_MUL_FP50_U32_MAX_K)
                simd_mul = nmod_mat_mul_fp50;  /* tiny inner dimension */
            else if (!blas_1pass)
                simd_mul = nmod_mat_mul_u32;
        }
        else if (NMOD_MAT_HAVE_MUL_U52 && NMOD_MAT_HAVE_FPV
                 && bits <= FLINT_NMOD_MAT_MUL_FP50_MAX_BITS)
        {
            /* small inner dimension on an IFMA machine, 33-50 bits */
            simd_mul = nmod_mat_mul_fp50;
        }
        else if (NMOD_MAT_HAVE_FPV
                 && FLINT_NMOD_MAT_MUL_K52_MIN_BITS > 0
                 && bits >= FLINT_NMOD_MAT_MUL_K52_MIN_BITS
#if FLINT_USES_BLAS
                 /* FLINT's own gemm never beats k52 / fp50 with its CRT
                    (1.3-1.8x slower on Apple M4 at 512-2048) */
                 && (FLINT_NMOD_MAT_MUL_K52_BLAS_CUTOFF <= 0
                     || min_dim < FLINT_NMOD_MAT_MUL_K52_BLAS_CUTOFF)
#endif
                )
        {
            /*
                33 to 52 bits without IFMA, where the alternative is
                4-5 dgemm passes and a CRT: the floating point kernel
                where the parameters prefer it (it stops below 2^50),
                the integer two-limb one otherwise. Past
                K52_BLAS_CUTOFF the multimodular route wins after all
                and this falls through to the dispatch below.
            */
            if (bits <= FLINT_NMOD_MAT_MUL_FP50_MAX_BITS)
                simd_mul = nmod_mat_mul_fp50;
            else
                simd_mul = nmod_mat_mul_k52;
        }

        if (simd_mul != NULL)
        {
            /* the cutoff must exceed 4: _nmod_mat_mul_strassen_rec sends
               dimensions <= 4 back to nmod_mat_mul (0: never) */
            if (flint_num_threads == 1
                    && FLINT_NMOD_MAT_MUL_SIMD_STRASSEN_CUTOFF > 4
                    && min_dim >= FLINT_NMOD_MAT_MUL_SIMD_STRASSEN_CUTOFF)
            {
                if (C == A || C == B)
                {
                    nmod_mat_t T;
                    nmod_mat_init(T, m, n, A->mod.n);
                    nmod_mat_mul_strassen(T, A, B);
                    nmod_mat_swap_entrywise(C, T);
                    nmod_mat_clear(T);
                }
                else
                    nmod_mat_mul_strassen(C, A, B);
                return;
            }

            if (simd_mul(C, A, B))
                return;
        }
    }
#endif

    /*
        tuning is based on several assumptions:
        (1) the gemm used by nmod_mat_mul_blas is at least as parallel
            as the rest of FLINT (it uses FLINT's thread pool, or an
            external BLAS with at least flint_num_threads threads).
        (2) nmod_mat_mul_blas (with crt) only beats
            nmod_mat_mul_classical on square multiplications
            of large enough dimension
        (3) if nmod_mat_mul_blas beats nmod_mat_mul_classical on
            square multiplications of size d, then it beats it on
            rectangular multiplications as long as all dimensions are >= d
    */
    if (FLINT_BITS == 64 && min_dim > 100)
    {
        flint_bitcnt_t bits = FLINT_BIT_COUNT(A->mod.n);

        if (FLINT_BIT_COUNT(k) + 2*bits < 53 + 5)
        {
            /* mul_blas definitely avoids the slow crt */
            cutoff = 100;
        }
        else if (flint_num_threads > 1)
        {
            /* mul_blas with crt is competing against mul_classical_threaded */
            bits = FLINT_MAX(bits, 32);
            cutoff = 100 + 5*flint_num_threads*bits/2;
        }
        else
        {
            /* mul_blas with crt is competing against mul_strassen */
            cutoff = 450;
        }

        if (min_dim > cutoff && nmod_mat_mul_blas(C, A, B))
            return;
    }

    if (C == A || C == B)
    {
        nmod_mat_t T;
        nmod_mat_init(T, m, n, A->mod.n);
        nmod_mat_mul(T, A, B);
        nmod_mat_swap_entrywise(C, T);
        nmod_mat_clear(T);
        return;
    }

    if (FLINT_BITS == 64 && C->mod.n < 2048)
        cutoff = 400;
    else
        cutoff = 200;

    if (flint_num_threads > 1)
        nmod_mat_mul_classical_threaded(C, A, B);
    else if (min_dim < cutoff)
        nmod_mat_mul_classical(C, A, B);
    else
        nmod_mat_mul_strassen(C, A, B);
}
