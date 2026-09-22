/*
    Copyright (C) 2024 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* parameters found for AMD Ryzen 7 PRO 7840U */

#ifndef FLINT_MPARAM_H
#define FLINT_MPARAM_H

/* TODO these were taken directly from zen3 flint-mparam.h  ----> */
#define FLINT_FFT_SMALL_MUL_THRESHOLD            240
#define FLINT_FFT_SMALL_SQR_THRESHOLD            400

#define FLINT_FFT_MUL_THRESHOLD                32000
#define FLINT_FFT_SQR_THRESHOLD                32000
/* <---- these were taken directly from zen3 flint-mparam.h  */

#define FFT_TAB \
   { { 4, 4 }, { 4, 3 }, { 3, 2 }, { 2, 2 }, { 2, 0 } }

#define MULMOD_TAB \
   { 4, 4, 4, 4, 4, 3, 3, 3, 3, 2, 2, 3, 3, 2, 2, 2, 2, 1, 1 }

#define FFT_N_NUM 19
#define FFT_MULMOD_2EXPP1_CUTOFF 128

#define FLINT_PREINVERT_LIMB_USE_NATIVE 1

#define FLINT_MULMOD_SHOUP_THRESHOLD 0

#define FLINT_MPN_MULHIGH_FFT_SMALL_CUTOFF 300

#define FLINT_MPN_SQRHIGH_FFT_SMALL_CUTOFF 470

/* Newton division, Hensel exact division and Newton square root beat
   GMP from these limb counts (divisor / quotient limbs for the first two,
   input limbs for the square root); Zen 3 values; not separately tuned on Zen 4 */
#define FLINT_MPN_TDIV_QR_NEWTON_CUTOFF 500
#define FLINT_MPN_DIVEXACT_NEWTON_CUTOFF 500
#define FLINT_MPN_SQRTREM_NEWTON_CUTOFF 2000


/*
    nmod_mat_mul: dispatch of the SIMD kernels nmod_mat_mul_u32 (moduli
    below 2^32), nmod_mat_mul_u52 (AVX512-IFMA, moduli up to 2^52) and,
    without IFMA, nmod_mat_mul_k52 / nmod_mat_mul_fp50 (moduli up to 2^52
    / below 2^50); see src/nmod_mat/mul.c and the profile p-mul_tune.c
    (measured on AMD Zen 4 (AVX512-IFMA), with an external BLAS (AOCL/BLIS)).
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
      K52_MIN_BITS         without IFMA, from this modulus bit size on
                           (through 52 bits) the two-limb integer kernel
                           nmod_mat_mul_k52 replaces blas + CRT (0: never)
      FP50_MAX_BITS        in that range, the floating point kernel
                           nmod_mat_mul_fp50 is preferred to k52 up to
                           this modulus bit size (through 50; 0: never)
      K52_BLAS_CUTOFF      in that range, nmod_mat_mul_blas and its CRT are
                           preferred to k52 / fp50 from this dimension on
                           (0: never). The crossover grows with the modulus
                           size, since blas needs more primes; the value is
                           the one for the bottom of the range.
*/
#define FLINT_NMOD_MAT_MUL_U32_MIN_DIM 8
#define FLINT_NMOD_MAT_MUL_U32_BLAS_CUTOFF 256
#define FLINT_NMOD_MAT_MUL_U32_STRASSEN_CUTOFF 512
#define FLINT_NMOD_MAT_MUL_U52_MIN_BITS 32
#define FLINT_NMOD_MAT_MUL_U52_LO_MAX_BITS 26
#define FLINT_NMOD_MAT_MUL_K52_MIN_BITS 33
#define FLINT_NMOD_MAT_MUL_FP50_MAX_BITS 50
#define FLINT_NMOD_MAT_MUL_K52_BLAS_CUTOFF 0

#endif
