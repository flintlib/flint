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

#define FLINT_FFT_SMALL_MUL_THRESHOLD            400
#define FLINT_FFT_SMALL_SQR_THRESHOLD            800

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

#endif
