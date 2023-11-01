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

#include "t-heuristic_dot.c"
#include "t-lll.c"
#include "t-lll_d.c"
#include "t-lll_d_heuristic.c"
#include "t-lll_d_heuristic_with_removal.c"
#include "t-lll_d_with_removal.c"
#include "t-lll_d_with_removal_knapsack.c"
#include "t-lll_mpf.c"
#include "t-lll_mpf_with_removal.c"
#include "t-lll_with_removal.c"
#include "t-wrapper.c"
#include "t-wrapper_with_removal.c"
#include "t-wrapper_with_removal_knapsack.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(fmpz_lll_heuristic_dot),
    TEST_FUNCTION(fmpz_lll),
    TEST_FUNCTION(fmpz_lll_d),
    TEST_FUNCTION(fmpz_lll_d_heuristic),
    TEST_FUNCTION(fmpz_lll_d_heuristic_with_removal),
    TEST_FUNCTION(fmpz_lll_d_with_removal),
    TEST_FUNCTION(fmpz_lll_d_with_removal_knapsack),
    TEST_FUNCTION(fmpz_lll_mpf),
    TEST_FUNCTION(fmpz_lll_mpf_with_removal),
    TEST_FUNCTION(fmpz_lll_with_removal),
    TEST_FUNCTION(fmpz_lll_wrapper),
    TEST_FUNCTION(fmpz_lll_wrapper_with_removal),
    TEST_FUNCTION(fmpz_lll_wrapper_with_removal_knapsack)
};

/* main function *************************************************************/

TEST_MAIN(tests)
