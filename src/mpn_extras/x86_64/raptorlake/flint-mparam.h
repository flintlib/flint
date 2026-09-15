/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* parameters found for Intel(R) Core(TM) i5-1335U */

#ifndef FLINT_MPARAM_H
#define FLINT_MPARAM_H

/* TODO these were taken directly from skylake flint-mparam.h  ----> */
#define FLINT_FFT_SMALL_MUL_THRESHOLD           500
#define FLINT_FFT_SMALL_SQR_THRESHOLD           500

#define FLINT_FFT_MUL_THRESHOLD                32000
#define FLINT_FFT_SQR_THRESHOLD                32000
/* <---- these were taken directly from skylake flint-mparam.h  */

#define FFT_TAB \
   { { 4, 4 }, { 4, 3 }, { 3, 2 }, { 2, 2 }, { 1, 1 } }

#define MULMOD_TAB \
   { 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 1, 1 }

#define FFT_N_NUM 19
#define FFT_MULMOD_2EXPP1_CUTOFF 128

#define FLINT_PREINVERT_LIMB_USE_NATIVE 1

#define FLINT_MULMOD_SHOUP_THRESHOLD 0

#define FLINT_MPN_MULHIGH_FFT_SMALL_CUTOFF 500

#define FLINT_MPN_SQRHIGH_FFT_SMALL_CUTOFF 580

/* see x86_64/flint-mparam.h; not tuned on this microarchitecture */
#define FLINT_MPN_TDIV_QR_NEWTON_CUTOFF 1024
#define FLINT_MPN_DIVEXACT_NEWTON_CUTOFF 1024
#define FLINT_MPN_SQRTREM_NEWTON_CUTOFF 2800


/*
    nmod_mat_mul: dispatch of the integer SIMD kernels nmod_mat_mul_u32
    (moduli below 2^32) and nmod_mat_mul_u52 (AVX512-IFMA, moduli up to
    2^52); see src/nmod_mat/mul.c and the profile p-mul_tune.c
    (not measured on this target: defaults).
      U32_MIN_DIM          use the SIMD kernels from this minimal dimension
      U32_BLAS_CUTOFF      when one dgemm pass suffices (k*(n/2)^2 < 2^53),
                           nmod_mat_mul_blas is preferred to u32 from this
                           dimension on (0: never)
      U32_STRASSEN_CUTOFF  single-threaded, one Strassen level is put on
                           top of u32 / u52 from this dimension on
      U52_MIN_BITS         u52 is preferred to u32 from this modulus bit
                           size on (through 52 bits)
      U52_LO_MAX_BITS      u52 in its single-IFMA mode is preferred to u32
                           up to this modulus bit size (0: never)
*/
#define FLINT_NMOD_MAT_MUL_U32_MIN_DIM 8
#define FLINT_NMOD_MAT_MUL_U32_BLAS_CUTOFF 256
#define FLINT_NMOD_MAT_MUL_U32_STRASSEN_CUTOFF 768
#define FLINT_NMOD_MAT_MUL_U52_MIN_BITS 31
#define FLINT_NMOD_MAT_MUL_U52_LO_MAX_BITS 0

#endif
