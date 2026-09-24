/*
    Copyright (C) 2025 Ricardo Buring

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Include functions *********************************************************/

#include "t-ring.c"
#include "t-set_str.c"
#include "t-sigma_delta.c"
#include "t-mul.c"
#include "t-divrem.c"
#include "t-apply.c"
#include "t-ddx_to_euler.c"
#include "t-euler_to_ddx.c"
#include "t-shift_convert.c"
#include "t-differential_shift.c"
#include "t-convert.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(gr_ore_poly_ring),
    TEST_FUNCTION(gr_ore_poly_set_str),
    TEST_FUNCTION(gr_ore_poly_sigma_delta),
    TEST_FUNCTION(gr_ore_poly_mul),
    TEST_FUNCTION(gr_ore_poly_divrem),
    TEST_FUNCTION(gr_ore_poly_apply),
    TEST_FUNCTION(gr_ore_poly_ddx_to_euler),
    TEST_FUNCTION(gr_ore_poly_euler_to_ddx),
    TEST_FUNCTION(gr_ore_poly_shift_convert),
    TEST_FUNCTION(gr_ore_poly_differential_to_shift),
    TEST_FUNCTION(gr_ore_poly_shift_to_differential),
    TEST_FUNCTION(gr_ore_poly_convert)
};

/* main function *************************************************************/

TEST_MAIN(tests)
