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
#include "t-gso.c"
#include "t-init_clear.c"
#include "t-is_empty.c"
#include "t-is_square.c"
#include "t-mul.c"
#include "t-one.c"
#include "t-qr.c"
#include "t-set_fmpz_mat.c"
#include "t-zero.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(mpf_mat_entry),
    TEST_FUNCTION(mpf_mat_equal),
    TEST_FUNCTION(mpf_mat_gso),
    TEST_FUNCTION(mpf_mat_init_clear),
    TEST_FUNCTION(mpf_mat_is_empty),
    TEST_FUNCTION(mpf_mat_is_square),
    TEST_FUNCTION(mpf_mat_mul),
    TEST_FUNCTION(mpf_mat_one),
    TEST_FUNCTION(mpf_mat_qr),
    TEST_FUNCTION(mpf_mat_set_fmpz_mat),
    TEST_FUNCTION(mpf_mat_zero)
};

/* main function *************************************************************/

TEST_MAIN(tests)
