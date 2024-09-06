/*
    Copyright (C) 2024 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* parameters found for Intel(R) Xeon(R) Gold 6354 */

#ifndef FLINT_MPARAM_H
#define FLINT_MPARAM_H

/* TODO these were taken directly from skylake flint-mparam.h  ----> */
#define FLINT_FFT_SMALL_MUL_THRESHOLD           1540
#define FLINT_FFT_SMALL_SQR_THRESHOLD           3080

#define FLINT_FFT_MUL_THRESHOLD                32000
#define FLINT_FFT_SQR_THRESHOLD                32000
/* <---- these were taken directly from skylake flint-mparam.h  */

#define FFT_TAB \
   { { 4, 4 }, { 4, 3 }, { 3, 2 }, { 2, 2 }, { 2, 1 } }

#define MULMOD_TAB \
   { 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 1, 1 }

#define FFT_N_NUM 19
#define FFT_MULMOD_2EXPP1_CUTOFF 128

#define FLINT_MULMOD_SHOUP_THRESHOLD 0

#endif
