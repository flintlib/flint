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

#include "t-mpn_mod.c"
#include "t-mat_lu.c"
#include "t-mat_lu_classical_delayed.c"
#include "t-mat_mul_multi_mod.c"
#include "t-mat_mul_waksman.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(mpn_mod),
    TEST_FUNCTION(mpn_mod_mat_lu),
    TEST_FUNCTION(mpn_mod_mat_lu_classical_delayed),
    TEST_FUNCTION(mpn_mod_mat_mul_multi_mod),
    TEST_FUNCTION(mpn_mod_mat_mul_waksman),
};

/* main function *************************************************************/

TEST_MAIN(tests)
