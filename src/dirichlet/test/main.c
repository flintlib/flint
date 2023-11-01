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

#include "t-char.c"
#include "t-properties.c"
#include "t-vec.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(dirichlet_char),
    TEST_FUNCTION(dirichlet_properties),
    TEST_FUNCTION(dirichlet_vec)
};

/* main function *************************************************************/

TEST_MAIN(tests)
