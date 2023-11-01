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

#include "t-add.c"
#include "t-get_set_entry_padic.c"
#include "t-get_set_fmpq_mat.c"
#include "t-mul.c"
#include "t-neg.c"
#include "t-scalar_div_fmpz.c"
#include "t-scalar_mul_fmpz.c"
#include "t-scalar_mul_padic.c"
#include "t-sub.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(padic_mat_add),
    TEST_FUNCTION(padic_mat_get_set_entry_padic),
    TEST_FUNCTION(padic_mat_get_set_fmpq_mat),
    TEST_FUNCTION(padic_mat_mul),
    TEST_FUNCTION(padic_mat_neg),
    TEST_FUNCTION(padic_mat_scalar_div_fmpz),
    TEST_FUNCTION(padic_mat_scalar_mul_fmpz),
    TEST_FUNCTION(padic_mat_scalar_mul_padic),
    TEST_FUNCTION(padic_mat_sub)
};

/* main function *************************************************************/

TEST_MAIN(tests)
