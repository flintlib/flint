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

/* Newton division, Newton Hensel (exact) division and Newton square root
   replace the divide and conquer code from these limb counts (divisor /
   quotient limbs for the first two, divisor limbs for the next two with
   a quotient at least as long, respectively at least four times as long,
   input limbs for the square root). Not tuned on this target: Zen 3
   values from src/mpn_extras/tune/tune-div.c. */
#define FLINT_MPN_TDIV_QR_NEWTON_CUTOFF 760
#define FLINT_MPN_TDIV_QR_NEWTON_LONG_CUTOFF 487
#define FLINT_MPN_DIVEXACT_NEWTON_CUTOFF 950
#define FLINT_MPN_DIVEXACT_UNBALANCED_CUTOFF 340
#define FLINT_MPN_SQRTREM_NEWTON_CUTOFF 5155

/* division cutoffs (Zen 3 values, not tuned on this target) */
#define FLINT_MPN_DIV_DC_CUTOFF 20
#define FLINT_MPN_DIVAPPR_DC_CUTOFF 64
#define FLINT_MPN_DIVAPPROX_SHORT_CUTOFF 56
#define FLINT_MPN_TDIV_Q_DC_CUTOFF 112
#define FLINT_MPN_DIVAPPROX_NEWTON_CUTOFF 378
#define FLINT_MPN_INV_NEWTON_CUTOFF 71
#define FLINT_MPN_INV_NEWTON_LONG_CUTOFF 71
#define FLINT_MPN_DC_BDIV_QR_CUTOFF 64
#define FLINT_MPN_DC_BDIV_Q_CUTOFF 120

/* one-limb division by a chain of hardware divisions below these dividend
   lengths (unnormalized / normalized divisor) */
#define FLINT_MPN_DIVREM_1_HW_CUTOFF 20
#define FLINT_MPN_DIVREM_1_NORM_HW_CUTOFF 8

/* two-limb divisors by hardware 2/1 divisions without an inverse below this
   dividend length */
#define FLINT_MPN_DIV_2_HW_CUTOFF 24

/* 3- to 7-limb divisors by hardware 2/1 divisions without an inverse for
   quotients shorter than this */
#define FLINT_MPN_DIV_SMALL_HW_QN_CUTOFF 6

/* two-limb divisors by GMP's assembly mpn_divrem_2 (when available) from
   this dividend length (1000000: never) */
#define FLINT_MPN_DIV_2_GMP_CUTOFF 24

#endif
