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
#include "t-dot2.c"
#include "t-dot.c"
#include "t-init_clear.c"
#include "t-norm2.c"
#include "t-norm.c"
#include "t-scalar_mul_2exp.c"
#include "t-scalar_mul_mpf.c"
#include "t-set_equal.c"
#include "t-set_fmpz_vec.c"
#include "t-sub.c"
#include "t-zero.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(mpf_vec_add),
    TEST_FUNCTION(mpf_vec_dot2),
    TEST_FUNCTION(mpf_vec_dot),
    TEST_FUNCTION(mpf_vec_init_clear),
    TEST_FUNCTION(mpf_vec_norm2),
    TEST_FUNCTION(mpf_vec_norm),
    TEST_FUNCTION(mpf_vec_scalar_mul_2exp),
    TEST_FUNCTION(mpf_vec_scalar_mul_mpf),
    TEST_FUNCTION(mpf_vec_set_equal),
    TEST_FUNCTION(mpf_vec_set_fmpz_vec),
    TEST_FUNCTION(mpf_vec_sub),
    TEST_FUNCTION(mpf_vec_zero)
};

/* main function *************************************************************/

TEST_MAIN(tests)
