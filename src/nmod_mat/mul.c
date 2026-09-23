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
          (K52_BLAS_CUTOFF): on Apple M4 that dimension is 320-448,
          Accelerate's dgemm being far out of reach of a NEON kernel; on
          the x86 machines measured it is 768 and beyond, or never.
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

            blas_1pass = cut > 0 && min_dim >= cut
                && (half == 0
                    || (ulong) k <= ((UWORD(1) << 53) - 1) / (half * half));
        }

        if (NMOD_MAT_HAVE_MUL_U52
                && (bits >= FLINT_NMOD_MAT_MUL_U52_MIN_BITS
                    || bits <= FLINT_NMOD_MAT_MUL_U52_LO_MAX_BITS))
        {
            /* FLINT's own gemm was never measured faster than u52 */
#if FLINT_USES_BLAS
            if (!blas_1pass)
#endif
                simd_mul = nmod_mat_mul_u52;
        }
        else if (bits <= 32)
        {
            if (!blas_1pass)
                simd_mul = nmod_mat_mul_u32;
        }
        else if (FLINT_NMOD_MAT_MUL_K52_MIN_BITS > 0
                 && bits >= FLINT_NMOD_MAT_MUL_K52_MIN_BITS
                 && (FLINT_NMOD_MAT_MUL_K52_BLAS_CUTOFF <= 0
                     || min_dim < FLINT_NMOD_MAT_MUL_K52_BLAS_CUTOFF))
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
            if (flint_num_threads == 1
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
