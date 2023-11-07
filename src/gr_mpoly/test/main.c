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

#include "t-add_sub.c"
#include "t-gen.c"
#include "t-get_set_coeff.c"
#include "t-mul_johnson.c"
#include "t-mul_monomial.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(gr_mpoly_add_sub),
    TEST_FUNCTION(gr_mpoly_gen),
    TEST_FUNCTION(gr_mpoly_get_set_coeff),
    TEST_FUNCTION(gr_mpoly_mul_johnson),
    TEST_FUNCTION(gr_mpoly_mul_monomial)
};

/* main function *************************************************************/

TEST_MAIN(tests)
