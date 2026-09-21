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

#define FLINT_FFT_SMALL_MUL_THRESHOLD            810
#define FLINT_FFT_SMALL_SQR_THRESHOLD            865

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

#define FLINT_MPN_MULHIGH_FFT_SMALL_CUTOFF 500

#define FLINT_MPN_SQRHIGH_FFT_SMALL_CUTOFF 580

/* Newton division, Newton Hensel (exact) division and Newton square root
   replace the divide and conquer code from these limb counts (divisor /
   quotient limbs for the first two, divisor limbs for the next two with
   a quotient at least as long, respectively at least four times as long,
   input limbs for the square root). Not tuned on this target: Skylake values
   from src/mpn_extras/tune/tune-div.c, except for
   FLINT_MPN_DIVEXACT_NEWTON_CUTOFF (measured against GMP's mpn_divexact
   before FLINT had its own schoolbook and divide and conquer Hensel
   division). */
#define FLINT_MPN_TDIV_QR_NEWTON_CUTOFF 1024
#define FLINT_MPN_TDIV_QR_NEWTON_LONG_CUTOFF 608
#define FLINT_MPN_DIVEXACT_NEWTON_CUTOFF 700
#define FLINT_MPN_DIVEXACT_UNBALANCED_CUTOFF 425
#define FLINT_MPN_SQRTREM_NEWTON_CUTOFF 5155

/* division cutoffs (Skylake values, not tuned on this target) */
#define FLINT_MPN_DIV_DC_CUTOFF 14
#define FLINT_MPN_DIVAPPR_DC_CUTOFF 20
#define FLINT_MPN_DIVAPPROX_SHORT_CUTOFF 56
#define FLINT_MPN_TDIV_Q_DC_CUTOFF 58
#define FLINT_MPN_DIVAPPROX_NEWTON_CUTOFF 590
#define FLINT_MPN_INV_NEWTON_CUTOFF 57
#define FLINT_MPN_INV_NEWTON_LONG_CUTOFF 57
#define FLINT_MPN_DC_BDIV_QR_CUTOFF 16
#define FLINT_MPN_DC_BDIV_Q_CUTOFF 44

/* one-limb division by a chain of hardware divisions below these dividend
   lengths (unnormalized / normalized divisor); not measured on this target,
   where udiv_qrnnd is the two-step generic C version */
#define FLINT_MPN_DIVREM_1_HW_CUTOFF 4
#define FLINT_MPN_DIVREM_1_NORM_HW_CUTOFF 4

/* two-limb divisors by hardware 2/1 divisions without an inverse below this
   dividend length; disabled on this target, where udiv_qrnnd is the
   two-step generic C version (not measured) */
#define FLINT_MPN_DIV_2_HW_CUTOFF 0

/* 3- to 7-limb divisors by hardware 2/1 divisions without an inverse for
   quotients shorter than this */
#define FLINT_MPN_DIV_SMALL_HW_QN_CUTOFF 0

/* two-limb divisors by GMP's assembly mpn_divrem_2 (when available) from
   this dividend length (1000000: never) */
#define FLINT_MPN_DIV_2_GMP_CUTOFF 1000000

#endif
