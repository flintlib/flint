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

#include "t-all_pairs_longest_walk.c"
#include "t-complement.c"
#include "t-is_diagonal.c"
#include "t-is_nilpotent.c"
#include "t-is_transitive.c"
#include "t-mul.c"
#include "t-nilpotency_degree.c"
#include "t-trace.c"
#include "t-transitive_closure.c"
#include "t-transpose.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(bool_mat_all_pairs_longest_walk),
    TEST_FUNCTION(bool_mat_complement),
    TEST_FUNCTION(bool_mat_is_diagonal),
    TEST_FUNCTION(bool_mat_is_nilpotent),
    TEST_FUNCTION(bool_mat_is_transitive),
    TEST_FUNCTION(bool_mat_mul),
    TEST_FUNCTION(bool_mat_nilpotency_degree),
    TEST_FUNCTION(bool_mat_trace),
    TEST_FUNCTION(bool_mat_transitive_closure),
    TEST_FUNCTION(bool_mat_transpose)
};

/* main function *************************************************************/

TEST_MAIN(tests)
