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

#include "t-fmpz_poly_mul.c"
#include "t-mpn_add_inplace_c.c"
#include "t-mul.c"
#include "t-nmod_poly_divrem.c"
#include "t-nmod_poly_mul.c"
#include "t-sd_fft.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(_fmpz_poly_mul_mid_mpn_ctx),
    TEST_FUNCTION(flint_mpn_add_inplace_c),
    TEST_FUNCTION(mpn_ctx_mpn_mul),
    TEST_FUNCTION(_nmod_poly_divrem_mpn_ctx),
    TEST_FUNCTION(_nmod_poly_mul_mid_mpn_ctx),
    TEST_FUNCTION(sd_fft)
};

/* main function *************************************************************/

TEST_MAIN(tests)
