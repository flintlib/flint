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

#include "t-bound_2exp_si.c"
#include "t-fmpq_ui.c"
#include "t-fmpq_ui_multi_mod.c"
#include "t-fmpq_vec.c"
#include "t-mod_p_harvey.c"
#include "t-rev.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(bernoulli_bound_2exp_si),
    TEST_FUNCTION(bernoulli_fmpq_ui),
    TEST_FUNCTION(bernoulli_fmpq_ui_multi_mod),
    TEST_FUNCTION(bernoulli_fmpq_vec),
    TEST_FUNCTION(bernoulli_mod_p_harvey),
    TEST_FUNCTION(bernoulli_rev)
};

/* main function *************************************************************/

TEST_MAIN(tests)
