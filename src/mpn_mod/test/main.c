/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <stdlib.h>

/* Include functions *********************************************************/

#include "t-mpn_mod.c"
#include "t-mat.c"
#include "t-mat_lu.c"
#include "t-mat_lu_classical_delayed.c"
#include "t-mat_mul_multi_mod.c"
#include "t-mat_mul_waksman.c"
#include "t-poly_divrem_basecase.c"
#include "t-poly.c"
#include "t-poly_mullow.c"
#include "t-poly_mullow_KS.c"
#include "t-poly_mullow_classical.c"
#include "t-poly_mullow_fft_small.c"
#include "t-poly_mullow_karatsuba.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(mpn_mod),
    TEST_FUNCTION(mpn_mod_mat),
    TEST_FUNCTION(mpn_mod_mat_lu),
    TEST_FUNCTION(mpn_mod_mat_lu_classical_delayed),
    TEST_FUNCTION(mpn_mod_mat_mul_multi_mod),
    TEST_FUNCTION(mpn_mod_mat_mul_waksman),
    TEST_FUNCTION(mpn_mod_poly_divrem_basecase),
    TEST_FUNCTION(mpn_mod_poly),
    TEST_FUNCTION(mpn_mod_poly_mullow),
    TEST_FUNCTION(mpn_mod_poly_mullow_KS),
    TEST_FUNCTION(mpn_mod_poly_mullow_classical),
    TEST_FUNCTION(mpn_mod_poly_mullow_fft_small),
    TEST_FUNCTION(mpn_mod_poly_mullow_karatsuba),
};

/* main function *************************************************************/

TEST_MAIN(tests)
