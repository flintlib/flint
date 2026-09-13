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

/* thresholds for nmod_mat_mul_u32: see the "Moduli below 2^32 comment in nmod_mat_mul */
#ifndef NMOD_MAT_MUL_U32_MIN_DIM
# define NMOD_MAT_MUL_U32_MIN_DIM 8
#endif
#ifndef NMOD_MAT_MUL_U32_CRT_BITS
# define NMOD_MAT_MUL_U32_CRT_BITS 23
#endif
#ifndef NMOD_MAT_MUL_U32_SMALL_BITS
# define NMOD_MAT_MUL_U32_SMALL_BITS 20
#endif
#ifndef NMOD_MAT_MUL_U32_MID_DIM
# define NMOD_MAT_MUL_U32_MID_DIM 256
#endif
#ifndef NMOD_MAT_MUL_U32_SMALL_DIM
# define NMOD_MAT_MUL_U32_SMALL_DIM 256
#endif
#ifndef NMOD_MAT_MUL_U32_STRASSEN_DIM
# define NMOD_MAT_MUL_U32_STRASSEN_DIM 512
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
        Moduli below 2^32: integer SIMD kernels (nmod_mat_mul_u32) with
        delayed reduction. Thresholds from src/nmod_mat/profile/p-mul_u32.c
        on a Zen 4 (AVX-512, single thread) and an Intel Sapphire Rapids
        machine; see that file to redo the measurements.

        - Above 2^23, where a single dgemm cannot hold the dot products of
          any useful length and mul_blas goes through CRT, u32 is 3-4x
          faster than mul_blas from min_dim about 8 on.
        - Up to 2^23, u32 beats mul_blas below a dimension of about 256
          (mul_blas pays O(n^2) conversions and thread hand-offs around
          its gemm), and when mul_blas would need the CRT. Above that the
          two are within a few percent of each other on Zen 4, where a
          vpmuldq retires at the same rate as an FMA, while on Intel
          AVX-512 (two FMA units, one vpmuldq port) a single dgemm pass
          is about 1.4x faster, so mul_blas keeps that range.
        - From min_dim 512 on, single threaded, one Strassen level on top
          (its recursive calls come back here and land in u32) gains
          5-10%. With several threads u32 is called directly: it splits
          C across the pool itself, and the threaded case is unmeasured.
    */
#if FLINT_BITS == 64
    if (C->mod.n < (UWORD(1) << 32) && min_dim >= NMOD_MAT_MUL_U32_MIN_DIM)
    {
        flint_bitcnt_t bits = FLINT_BIT_COUNT(C->mod.n);
        int use_u32;

        if (bits > NMOD_MAT_MUL_U32_CRT_BITS)
            use_u32 = 1;
        else if (bits > NMOD_MAT_MUL_U32_SMALL_BITS)
            use_u32 = (FLINT_BIT_COUNT(k) + 2*bits >= 53 + 5)
                      || (min_dim < NMOD_MAT_MUL_U32_MID_DIM);
        else
            use_u32 = (min_dim < NMOD_MAT_MUL_U32_SMALL_DIM);

        if (use_u32)
        {
            if (flint_num_threads == 1
                    && min_dim >= NMOD_MAT_MUL_U32_STRASSEN_DIM)
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

            if (nmod_mat_mul_u32(C, A, B))
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
