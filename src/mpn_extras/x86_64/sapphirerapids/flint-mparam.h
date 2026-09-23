/*
    Copyright (C) 2024 Vincent Neiger
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Sapphire Rapids, Emerald Rapids and Granite Rapids. The multiplication
   parameters are those of the Ice Lake file (Intel(R) Xeon(R) Gold 6354);
   the division and square root cutoffs were measured on an Emerald Rapids
   Xeon (with FLINT_PREINVERT_LIMB_USE_NATIVE, as these cores have a fast
   hardware divider) against GMP 6.3.0 built for the CPU (--host=skylake:
   GMP's config.guess takes these cores for nehalem, and distribution
   builds are generic; FLINT's schoolbook and Newton divisions use GMP's
   multiplication, so the cutoffs depend on it). */

#ifndef FLINT_MPARAM_H
#define FLINT_MPARAM_H

/* TODO these were taken directly from skylake flint-mparam.h  ----> */
#define FLINT_FFT_SMALL_MUL_THRESHOLD           500
#define FLINT_FFT_SMALL_SQR_THRESHOLD           500

#define FLINT_FFT_MUL_THRESHOLD                32000
#define FLINT_FFT_SQR_THRESHOLD                32000
/* <---- these were taken directly from skylake flint-mparam.h  */

#define FFT_TAB \
   { { 4, 4 }, { 4, 3 }, { 3, 2 }, { 2, 2 }, { 2, 1 } }

#define MULMOD_TAB \
   { 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 1, 1 }

#define FFT_N_NUM 19
#define FFT_MULMOD_2EXPP1_CUTOFF 128

#define FLINT_PREINVERT_LIMB_USE_NATIVE 1

#define FLINT_MULMOD_SHOUP_THRESHOLD 0

#define FLINT_MPN_MULHIGH_FFT_SMALL_CUTOFF 500

#define FLINT_MPN_SQRHIGH_FFT_SMALL_CUTOFF 580

/* Newton division, Newton Hensel (exact) division and Newton square root
   replace the divide and conquer code from these limb counts (divisor /
   quotient limbs for the first two, divisor limbs for the next two with
   a quotient at least as long, respectively at least four times as long,
   input limbs for the square root); tuned with
   src/mpn_extras/tune/tune-div.c */
#define FLINT_MPN_TDIV_QR_NEWTON_CUTOFF 950
#define FLINT_MPN_TDIV_QR_NEWTON_LONG_CUTOFF 608
#define FLINT_MPN_DIVEXACT_NEWTON_CUTOFF 1187
#define FLINT_MPN_DIVEXACT_UNBALANCED_CUTOFF 531
#define FLINT_MPN_SQRTREM_NEWTON_CUTOFF 7423

/* division cutoffs, tuned with src/mpn_extras/tune/tune-div.c */
#define FLINT_MPN_DIV_DC_CUTOFF 18
#define FLINT_MPN_DIVAPPR_DC_CUTOFF 34
#define FLINT_MPN_DIVAPPROX_SHORT_CUTOFF 45
#define FLINT_MPN_TDIV_Q_DC_CUTOFF 72
#define FLINT_MPN_DIVAPPROX_NEWTON_CUTOFF 921
#define FLINT_MPN_INV_NEWTON_CUTOFF 37
#define FLINT_MPN_INV_NEWTON_LONG_CUTOFF 24
#define FLINT_MPN_INV_NEWTON_VERYLONG_CUTOFF 171
#define FLINT_MPN_DC_BDIV_QR_CUTOFF 40
#define FLINT_MPN_DC_BDIV_Q_CUTOFF 80

/* one-limb division by a chain of hardware divisions below these dividend
   lengths (unnormalized / normalized divisor) */
#define FLINT_MPN_DIVREM_1_HW_CUTOFF 64
#define FLINT_MPN_DIVREM_1_NORM_HW_CUTOFF 64

/* two-limb divisors by hardware 2/1 divisions without an inverse below this
   dividend length */
#define FLINT_MPN_DIV_2_HW_CUTOFF 44

/* 3- to 7-limb divisors by hardware 2/1 divisions without an inverse for
   quotients shorter than this */
#define FLINT_MPN_DIV_SMALL_HW_QN_CUTOFF 9

/* two-limb divisors by GMP's assembly mpn_divrem_2 (when available) from
   this dividend length (1000000: never) */
#define FLINT_MPN_DIV_2_GMP_CUTOFF 44

/*
    nmod_mat_mul: dispatch of the SIMD kernels nmod_mat_mul_u32 (moduli
    below 2^32), nmod_mat_mul_u52 (AVX512-IFMA, moduli up to 2^52) and,
    without IFMA, nmod_mat_mul_k52 / nmod_mat_mul_fp50 (moduli up to 2^52
    / below 2^50); see src/nmod_mat/mul.c and the profile p-mul_tune.c
    (not measured on Sapphire Rapids yet: the values of x86_64/icelake, the
    closest AVX512-IFMA machine that was).
      SIMD_MIN_DIM         use the SIMD kernels when B has at least this
                           many columns and A this many rows (or half as
                           many if B has 4x as many columns), for any inner
                           dimension
      SIMD_STRASSEN_CUTOFF single-threaded, one Strassen level is put on
                           top of the SIMD kernels from this dimension on
      BLAS_1PASS_CUTOFF    when one dgemm pass suffices (k*(n/2)^2 < 2^53),
                           nmod_mat_mul_blas is preferred to u32 from this
                           dimension on and, with an external BLAS, also to
                           the single-IFMA mode of u52 (0: never)
      BLAS_1PASS_CUTOFF_MT the same with 4 threads or more (with 2 or 3
                           threads: the geometric mean of the two)
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
      K52_BLAS_CUTOFF      in that range and with an external BLAS,
                           nmod_mat_mul_blas and its CRT are preferred to
                           k52 / fp50 from this dimension on (0: never).
                           The crossover grows with the modulus size, since
                           blas needs more primes; the value is the one for
                           the bottom of the range.
*/
#define FLINT_NMOD_MAT_MUL_SIMD_MIN_DIM 8
#define FLINT_NMOD_MAT_MUL_BLAS_1PASS_CUTOFF 0
#define FLINT_NMOD_MAT_MUL_BLAS_1PASS_CUTOFF_MT 0
#define FLINT_NMOD_MAT_MUL_SIMD_STRASSEN_CUTOFF 512
#define FLINT_NMOD_MAT_MUL_U52_MIN_BITS 31
#define FLINT_NMOD_MAT_MUL_U52_LO_MAX_BITS 26
#define FLINT_NMOD_MAT_MUL_K52_MIN_BITS 33
#define FLINT_NMOD_MAT_MUL_FP50_MAX_BITS 50
#define FLINT_NMOD_MAT_MUL_K52_BLAS_CUTOFF 0

#endif
