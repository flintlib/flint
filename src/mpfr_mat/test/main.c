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

#include "t-entry.c"
#include "t-equal.c"
#include "t-init_clear.c"
#include "t-mul_classical.c"
#include "t-zero.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(mpfr_mat_entry),
    TEST_FUNCTION(mpfr_mat_equal),
    TEST_FUNCTION(mpfr_mat_init_clear),
    TEST_FUNCTION(mpfr_mat_mul_classical),
    TEST_FUNCTION(mpfr_mat_zero)
};

/* main function *************************************************************/

TEST_MAIN(tests)
