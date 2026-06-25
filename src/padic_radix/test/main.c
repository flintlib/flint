/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Include functions *********************************************************/

#include "t-dot.c"
#include "t-exp.c"
#include "t-exp_balanced.c"
#include "t-exp_rectangular.c"
#include "t-padic.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(padic_radix_dot),
    TEST_FUNCTION(padic_radix_exp),
    TEST_FUNCTION(padic_radix_exp_balanced),
    TEST_FUNCTION(padic_radix_exp_rectangular),
    TEST_FUNCTION(padic_radix),
};

/* main function *************************************************************/

TEST_MAIN(tests)
