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

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(gr_ore_poly_ring),
    TEST_FUNCTION(gr_ore_poly_set_str),
    TEST_FUNCTION(gr_ore_poly_sigma_delta),
};

/* main function *************************************************************/

TEST_MAIN(tests)
