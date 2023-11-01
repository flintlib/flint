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
#include "t-add_fmpq.c"
#include "t-add_fmpz.c"
#include "t-div.c"
#include "t-div_fmpq.c"
#include "t-div_fmpz.c"
#include "t-inv.c"
#include "t-mul.c"
#include "t-mul_fmpq.c"
#include "t-mul_fmpz.c"
#include "t-randtest.c"
#include "t-sub.c"
#include "t-sub_fmpq.c"
#include "t-sub_fmpz.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(fmpz_mpoly_q_add),
    TEST_FUNCTION(fmpz_mpoly_q_add_fmpq),
    TEST_FUNCTION(fmpz_mpoly_q_add_fmpz),
    TEST_FUNCTION(fmpz_mpoly_q_div),
    TEST_FUNCTION(fmpz_mpoly_q_div_fmpq),
    TEST_FUNCTION(fmpz_mpoly_q_div_fmpz),
    TEST_FUNCTION(fmpz_mpoly_q_inv),
    TEST_FUNCTION(fmpz_mpoly_q_mul),
    TEST_FUNCTION(fmpz_mpoly_q_mul_fmpq),
    TEST_FUNCTION(fmpz_mpoly_q_mul_fmpz),
    TEST_FUNCTION(fmpz_mpoly_q_randtest),
    TEST_FUNCTION(fmpz_mpoly_q_sub),
    TEST_FUNCTION(fmpz_mpoly_q_sub_fmpq),
    TEST_FUNCTION(fmpz_mpoly_q_sub_fmpz)
};

/* main function *************************************************************/

TEST_MAIN(tests)
