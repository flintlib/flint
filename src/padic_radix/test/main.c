/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Include functions *********************************************************/

#include "t-charpoly.c"
#include "t-det.c"
#include "t-minplus.c"
#include "t-solve.c"
#include "t-dot.c"
#include "t-exp.c"
#include "t-exp_balanced.c"
#include "t-exp_rectangular.c"
#include "t-log.c"
#include "t-log_balanced.c"
#include "t-log_rectangular.c"
#include "t-mat_mul.c"
#include "t-ops.c"
#include "t-poly_mul.c"
#include "t-robust.c"
#include "t-padic.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(padic_radix_charpoly),
    TEST_FUNCTION(padic_radix_det),
    TEST_FUNCTION(padic_radix_minplus),
    TEST_FUNCTION(padic_radix_solve),
    TEST_FUNCTION(padic_radix_dot),
    TEST_FUNCTION(padic_radix_exp),
    TEST_FUNCTION(padic_radix_exp_balanced),
    TEST_FUNCTION(padic_radix_exp_rectangular),
    TEST_FUNCTION(padic_radix_log),
    TEST_FUNCTION(padic_radix_log_balanced),
    TEST_FUNCTION(padic_radix_log_rectangular),
    TEST_FUNCTION(padic_radix_mat_mul),
    TEST_FUNCTION(padic_radix_ops),
    TEST_FUNCTION(padic_radix_poly_mul),
    TEST_FUNCTION(padic_radix_robust),
    TEST_FUNCTION(padic_radix),
};

/* main function *************************************************************/

TEST_MAIN(tests)
