/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Include functions *********************************************************/

#include "t-ctx.c"
#include "t-point.c"
#include "t-aff_point.c"
#include "t-jac_point.c"
#include "t-convert.c"
#include "t-mul.c"
#include "t-inexact.c"
#include "t-char23.c"
#include "t-generic.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(gr_ec_ctx),
    TEST_FUNCTION(gr_ec_point),
    TEST_FUNCTION(gr_ec_aff_point),
    TEST_FUNCTION(gr_ec_jac_point),
    TEST_FUNCTION(gr_ec_convert),
    TEST_FUNCTION(gr_ec_mul),
    TEST_FUNCTION(gr_ec_inexact),
    TEST_FUNCTION(gr_ec_char23),
    TEST_FUNCTION(gr_ec_generic)
};

/* main function *************************************************************/

TEST_MAIN(tests)
