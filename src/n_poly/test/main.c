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

#include "t-n_fq_poly_add.c"
#include "t-n_fq_poly_divrem.c"
#include "t-n_fq_poly_gcd.c"
#include "t-n_fq_poly_mul.c"
#include "t-n_fq_poly_sub.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(n_fq_poly_add),
    TEST_FUNCTION(n_fq_poly_divrem),
    TEST_FUNCTION(n_fq_poly_gcd),
    TEST_FUNCTION(n_fq_poly_mul),
    TEST_FUNCTION(n_fq_poly_sub),
};

/* main function *************************************************************/

TEST_MAIN(tests)
