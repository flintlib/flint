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

/* fft_small: round quotients with the magic constant instead of vroundpd,
   which is 2 uops with 8 cycle latency on Intel cores since Skylake
   (Haswell/Broadwell: 2 uops, 6 cycles) */
#define FLINT_FFT_SMALL_ROUND_USE_NATIVE 0

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

#endif
