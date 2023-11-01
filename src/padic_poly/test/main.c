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

#include "t-add.c"
#include "t-compose.c"
#include "t-compose_pow.c"
#include "t-derivative.c"
#include "t-evaluate_padic.c"
#include "t-get_set_fmpq_poly.c"
#include "t-init_realloc_clear.c"
#include "t-inv_series.c"
#include "t-mul.c"
#include "t-neg.c"
#include "t-one.c"
#include "t-pow.c"
#include "t-shift_left_right.c"
#include "t-sub.c"
#include "t-truncate.c"
#include "t-zero.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(padic_poly_add),
    TEST_FUNCTION(padic_poly_compose),
    TEST_FUNCTION(padic_poly_compose_pow),
    TEST_FUNCTION(padic_poly_derivative),
    TEST_FUNCTION(padic_poly_evaluate_padic),
    TEST_FUNCTION(padic_poly_get_set_fmpq_poly),
    TEST_FUNCTION(padic_poly_init_realloc_clear),
    TEST_FUNCTION(padic_poly_inv_series),
    TEST_FUNCTION(padic_poly_mul),
    TEST_FUNCTION(padic_poly_neg),
    TEST_FUNCTION(padic_poly_one),
    TEST_FUNCTION(padic_poly_pow),
    TEST_FUNCTION(padic_poly_shift_left_right),
    TEST_FUNCTION(padic_poly_sub),
    TEST_FUNCTION(padic_poly_truncate),
    TEST_FUNCTION(padic_poly_zero)
};

/* main function *************************************************************/

TEST_MAIN(tests)
