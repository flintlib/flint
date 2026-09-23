/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FLINT_MPARAM_H
#define FLINT_MPARAM_H

#define FLINT_FFT_SMALL_MUL_THRESHOLD            240
#define FLINT_FFT_SMALL_SQR_THRESHOLD            400

#define FLINT_FFT_MUL_THRESHOLD                32000
#define FLINT_FFT_SQR_THRESHOLD                32000

#define FFT_TAB \
   { {4, 4}, {4, 3}, {3, 2}, {2, 1}, {2, 1} }

#define MULMOD_TAB \
   { 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1 }

#define FFT_N_NUM                                 19
#define FFT_MULMOD_2EXPP1_CUTOFF                 128

#define FLINT_PREINVERT_LIMB_USE_NATIVE 1

/* warning: set by default, likely not optimal         */
/* if you have the relevant architecture, you can help */
/* determine this by running profiling files:          */
/*    nmod_vec's p-scalar_mul, p-scalar_addmul         */
/*    nmod_mat's p-nmod_vec_mul                        */
#define FLINT_MULMOD_SHOUP_THRESHOLD 10

#define FLINT_MPN_MULHIGH_FFT_SMALL_CUTOFF 300

#define FLINT_MPN_SQRHIGH_FFT_SMALL_CUTOFF 470

/* Newton division, Newton Hensel (exact) division and Newton square root
   replace the divide and conquer code from these limb counts (divisor /
   quotient limbs for the first two, divisor limbs for the next two with
   a quotient at least as long, respectively at least four times as long,
   input limbs for the square root) */
#define FLINT_MPN_TDIV_QR_NEWTON_CUTOFF 760
#define FLINT_MPN_TDIV_QR_NEWTON_LONG_CUTOFF 487
#define FLINT_MPN_DIVEXACT_NEWTON_CUTOFF 950
#define FLINT_MPN_DIVEXACT_UNBALANCED_CUTOFF 340
#define FLINT_MPN_SQRTREM_NEWTON_CUTOFF 5155

/* division cutoffs, tuned with src/mpn_extras/tune/tune-div.c */
#define FLINT_MPN_DIV_DC_CUTOFF 20
#define FLINT_MPN_DIVAPPR_DC_CUTOFF 64
#define FLINT_MPN_DIVAPPROX_SHORT_CUTOFF 80
#define FLINT_MPN_TDIV_Q_DC_CUTOFF 112
#define FLINT_MPN_DIVAPPROX_NEWTON_CUTOFF 737
#define FLINT_MPN_INV_NEWTON_CUTOFF 71
#define FLINT_MPN_INV_NEWTON_LONG_CUTOFF 37
#define FLINT_MPN_INV_NEWTON_VERYLONG_CUTOFF 137
#define FLINT_MPN_DC_BDIV_QR_CUTOFF 36
#define FLINT_MPN_DC_BDIV_Q_CUTOFF 120

/* one-limb division by a chain of hardware divisions below these dividend
   lengths (unnormalized / normalized divisor) */
#define FLINT_MPN_DIVREM_1_HW_CUTOFF 28
#define FLINT_MPN_DIVREM_1_NORM_HW_CUTOFF 6

/* two-limb divisors by hardware 2/1 divisions without an inverse below this
   dividend length */
#define FLINT_MPN_DIV_2_HW_CUTOFF 21

/* 3- to 7-limb divisors by hardware 2/1 divisions without an inverse for
   quotients shorter than this */
#define FLINT_MPN_DIV_SMALL_HW_QN_CUTOFF 6

/* two-limb divisors by GMP's assembly mpn_divrem_2 (when available) from
   this dividend length (1000000: never) */
#define FLINT_MPN_DIV_2_GMP_CUTOFF 21


/*
    nmod_mat_mul: dispatch of the SIMD kernels nmod_mat_mul_u32 (moduli
    below 2^32), nmod_mat_mul_u52 (AVX512-IFMA, moduli up to 2^52) and,
    without IFMA, nmod_mat_mul_k52 / nmod_mat_mul_fp50 (moduli up to 2^52
    / below 2^50); see src/nmod_mat/mul.c and the profile p-mul_tune.c
    (not measured on this target: defaults).
      U32_MIN_DIM          use the SIMD kernels when B has at least this
                           many columns and A this many rows (or half as
                           many if B has 4x as many columns), for any inner
                           dimension
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
#define FLINT_NMOD_MAT_MUL_U32_STRASSEN_CUTOFF 768
#define FLINT_NMOD_MAT_MUL_U52_MIN_BITS 31
#define FLINT_NMOD_MAT_MUL_U52_LO_MAX_BITS 0
#define FLINT_NMOD_MAT_MUL_K52_MIN_BITS 33
#define FLINT_NMOD_MAT_MUL_FP50_MAX_BITS 50
#define FLINT_NMOD_MAT_MUL_K52_BLAS_CUTOFF 0

#endif
