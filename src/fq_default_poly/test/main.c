/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <stdlib.h>

/* Include functions *********************************************************/

#include "t-init.c"
#include "t-inlines.c"
#include "t-set_fmpz_poly.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(fq_default_poly_init),
    TEST_FUNCTION(fq_default_poly_inlines),
    TEST_FUNCTION(fq_default_poly_set_fmpz_poly)
};

/* main function *************************************************************/

TEST_MAIN(tests)
