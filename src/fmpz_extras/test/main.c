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

#include "t-add2_fmpz_si_inline.c"
#include "t-add_inline.c"
#include "t-add_si_inline.c"
#include "t-add_ui_inline.c"
#include "t-lshift_mpn.c"
#include "t-sub_si_inline.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(fmpz_add2_fmpz_si_inline),
    TEST_FUNCTION(fmpz_add_inline),
    TEST_FUNCTION(fmpz_add_si_inline),
    TEST_FUNCTION(fmpz_add_ui_inline),
    TEST_FUNCTION(fmpz_lshift_mpn),
    TEST_FUNCTION(fmpz_sub_si_inline)
};

/* main function *************************************************************/

TEST_MAIN(tests)
