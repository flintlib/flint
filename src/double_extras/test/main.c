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

#include "t-is_nan.c"
#include "t-lambertw.c"
#include "t-log2.c"
#include "t-randtest.c"
#include "t-randtest_signed.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(d_is_nan),
    TEST_FUNCTION(d_lambertw),
    TEST_FUNCTION(d_log2),
    TEST_FUNCTION(d_randtest),
    TEST_FUNCTION(d_randtest_signed)
};

/* main function *************************************************************/

TEST_MAIN(tests)
