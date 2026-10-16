/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Include functions *********************************************************/

#include "t-decfloat.c"
#include "t-round.c"
#include "t-set_arf.c"
#include "t-arith.c"
#include "t-str.c"
#include "t-conv.c"
#include "t-decmag.c"
#include "t-decball.c"
#include "t-ball_str.c"
#include "t-via_arb.c"
#include "t-deccfloat.c"
#include "t-carith.c"
#include "t-deccball.c"
#include "t-cfunctions.c"
#include "t-api.c"
#include "t-qqbar.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(decfloat),
    TEST_FUNCTION(decfloat_round),
    TEST_FUNCTION(decfloat_set_arf),
    TEST_FUNCTION(decfloat_arith),
    TEST_FUNCTION(decfloat_str),
    TEST_FUNCTION(decfloat_conv),
    TEST_FUNCTION(decmag),
    TEST_FUNCTION(decball),
    TEST_FUNCTION(decball_str),
    TEST_FUNCTION(decimal_via_arb),
    TEST_FUNCTION(deccfloat),
    TEST_FUNCTION(deccfloat_arith),
    TEST_FUNCTION(deccball),
    TEST_FUNCTION(deccfloat_functions),
    TEST_FUNCTION(decimal_api),
    TEST_FUNCTION(decimal_qqbar),
};

/* main function *************************************************************/

TEST_MAIN(tests)
