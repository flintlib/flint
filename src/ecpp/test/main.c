/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.1

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Include functions *********************************************************/

#include "t-cert_str.c"
#include "t-class_poly_genus.c"
#include "t-class_poly_tower.c"
#include "t-cornacchia.c"
#include "t-point_mul.c"
#include "t-prove.c"
#include "t-root_radicals.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(ecpp_cert_str),
    TEST_FUNCTION(ecpp_class_poly_genus),
    TEST_FUNCTION(ecpp_class_poly_tower),
    TEST_FUNCTION(ecpp_cornacchia),
    TEST_FUNCTION(ecpp_point_mul),
    TEST_FUNCTION(ecpp_prove),
    TEST_FUNCTION(ecpp_root_radicals)
};

/* main function *************************************************************/

TEST_MAIN(tests)
