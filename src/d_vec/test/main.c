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

#include "t-add.c"
#include "t-dot.c"
#include "t-dot_heuristic.c"
#include "t-dot_thrice.c"
#include "t-init_clear.c"
#include "t-norm.c"
#include "t-set_equal.c"
#include "t-sub.c"
#include "t-zero.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(d_vec_add),
    TEST_FUNCTION(d_vec_dot),
    TEST_FUNCTION(d_vec_dot_heuristic),
    TEST_FUNCTION(d_vec_dot_thrice),
    TEST_FUNCTION(d_vec_init_clear),
    TEST_FUNCTION(d_vec_norm),
    TEST_FUNCTION(d_vec_set_equal),
    TEST_FUNCTION(d_vec_sub),
    TEST_FUNCTION(d_vec_zero)
};

/* main function *************************************************************/

TEST_MAIN(tests)
