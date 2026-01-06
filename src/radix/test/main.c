/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Include functions *********************************************************/

#include "t-add.c"
#include "t-divrem_1.c"
#include "t-get_mpn.c"
#include "t-mulmid_classical.c"
#include "t-mulmid_fft_small.c"
#include "t-mulmid_KS.c"
#include "t-neg.c"
#include "t-set_mpn.c"
#include "t-sub.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(radix_add),
    TEST_FUNCTION(radix_divrem_1),
    TEST_FUNCTION(radix_get_mpn),
    TEST_FUNCTION(radix_mulmid_classical),
    TEST_FUNCTION(radix_mulmid_fft_small),
    TEST_FUNCTION(radix_mulmid_KS),
    TEST_FUNCTION(radix_neg),
    TEST_FUNCTION(radix_set_mpn),
    TEST_FUNCTION(radix_sub),
};

/* main function *************************************************************/

TEST_MAIN(tests)
