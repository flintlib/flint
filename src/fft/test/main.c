/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <stdlib.h>

/* Include functions *********************************************************/

#include "t-adjust.c"
#include "t-adjust_sqrt2.c"
#include "t-butterfly.c"
#include "t-butterfly_lshB.c"
#include "t-butterfly_rshB.c"
#include "t-butterfly_sqrt2.c"
#include "t-butterfly_twiddle.c"
#include "t-convolution.c"
#include "t-convolution_precache.c"
#include "t-div_2expmod_2expp1.c"
#include "t-fft_ifft_mfa_truncate_sqrt2.c"
#include "t-fft_ifft_negacyclic.c"
#include "t-fft_ifft_radix2.c"
#include "t-fft_ifft_truncate.c"
#include "t-fft_ifft_truncate_sqrt2.c"
#include "t-mul_2expmod_2expp1.c"
#include "t-mul_fft_main.c"
#include "t-mul_mfa_truncate_sqrt2.c"
#include "t-mulmod_2expp1.c"
#include "t-mul_truncate_sqrt2.c"
#include "t-negmod_2expp1.c"
#include "t-normmod_2expp1.c"
#include "t-split_combine_bits.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(fft_adjust),
    TEST_FUNCTION(fft_adjust_sqrt2),
    TEST_FUNCTION(fft_ifft_butterfly),
    TEST_FUNCTION(butterfly_lshB),
    TEST_FUNCTION(butterfly_rshB),
    TEST_FUNCTION(fft_ifft_butterfly_sqrt2),
    TEST_FUNCTION(fft_ifft_butterfly_twiddle),
    TEST_FUNCTION(fft_convolution),
    TEST_FUNCTION(fft_convolution_precache),
    TEST_FUNCTION(mpn_div_2expmod_2expp1),
    TEST_FUNCTION(fft_ifft_mfa_truncate_sqrt2),
    TEST_FUNCTION(fft_ifft_negacyclic),
    TEST_FUNCTION(fft_ifft_radix2),
    TEST_FUNCTION(fft_ifft_truncate),
    TEST_FUNCTION(fft_ifft_truncate_sqrt2),
    TEST_FUNCTION(mpn_mul_2expmod_2expp1),
    TEST_FUNCTION(flint_mpn_mul_fft_main),
    TEST_FUNCTION(mul_mfa_truncate_sqrt2),
    TEST_FUNCTION(fft_mulmod_2expp1),
    TEST_FUNCTION(mul_truncate_sqrt2),
    TEST_FUNCTION(mpn_negmod_2expp1),
    TEST_FUNCTION(mpn_normmod_2expp1),
    TEST_FUNCTION(fft_split_combine_bits)
};

/* main function *************************************************************/

TEST_MAIN(tests)
