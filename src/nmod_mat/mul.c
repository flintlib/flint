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
        Moduli up to 2^52: integer SIMD kernels with delayed reduction,
        nmod_mat_mul_u32 (any 64-bit target, moduli below 2^32) and
        nmod_mat_mul_u52 (AVX512-IFMA, moduli up to 2^52). The parameters
        come from flint-mparam.h and were measured with
        src/nmod_mat/profile/p-mul_u32.c; the picture on the machines
        measured so far (Ice Lake, Meteor Lake, Zen 4, Apple M4) is:

        - Where nmod_mat_mul_blas needs several dgemm passes and a CRT
          (k*(n/2)^2 >= 2^53, always the case from 25 bits on), u32 is
          2-5x faster than any other method from dimension 8 on.
        - Where one dgemm pass suffices, u32 wins below a dimension of
          about 100-250 (mul_blas pays O(n^2) conversions and thread
          hand-offs around its gemm) and is 0.7-1.1x of mul_blas above.
        - u52 removes the 31-32 bit cliff of u32 (whose in-kernel folds
          then come every 1-2 products) and extends the single pass to
          52 bits, where the alternative is 4-5 dgemm passes. Below 2^26
          it has a single-IFMA mode, faster than u32 on cores with two
          IFMA ports (Intel) and expected slower on Zen 4.
        - Single-threaded, one Strassen level on top (its recursive
          calls come back here) gains 5-15% from about 768 on. With
          several threads the kernels split C across the pool themselves.
    */
#if FLINT_BITS == 64
    if (min_dim >= FLINT_NMOD_MAT_MUL_U32_MIN_DIM
            && C->mod.n <= (UWORD(1) << 52))
    {
        flint_bitcnt_t bits = FLINT_BIT_COUNT(C->mod.n);
        int (* simd_mul)(nmod_mat_t, const nmod_mat_t, const nmod_mat_t);

        simd_mul = NULL;

        if (NMOD_MAT_HAVE_MUL_U52
                && (bits >= FLINT_NMOD_MAT_MUL_U52_MIN_BITS
                    || bits <= FLINT_NMOD_MAT_MUL_U52_LO_MAX_BITS))
        {
            simd_mul = nmod_mat_mul_u52;
        }
        else if (bits <= 32)
        {
            /* can mul_blas do it in one dgemm pass, k*(n/2)^2 < 2^53 ? */
            ulong half = C->mod.n / 2;
            int one_pass = (half == 0)
                    || ((ulong) k <= ((UWORD(1) << 53) - 1) / (half * half));

            if (!one_pass || min_dim < FLINT_NMOD_MAT_MUL_U32_BLAS_CUTOFF)
                simd_mul = nmod_mat_mul_u32;
        }

        if (simd_mul != NULL)
        {
            if (flint_num_threads == 1
                    && min_dim >= FLINT_NMOD_MAT_MUL_U32_STRASSEN_CUTOFF)
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
