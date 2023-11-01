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
#include "t-is_zero.c"
#include "t-neg.c"
#include "t-sub.c"
#include "t-swap.c"
#include "t-zero.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(fq_nmod_vec_add),
    TEST_FUNCTION(fq_nmod_vec_is_zero),
    TEST_FUNCTION(fq_nmod_vec_neg),
    TEST_FUNCTION(fq_nmod_vec_sub),
    TEST_FUNCTION(fq_nmod_vec_swap),
    TEST_FUNCTION(fq_nmod_vec_zero)
};

/* main function *************************************************************/

TEST_MAIN(tests)
