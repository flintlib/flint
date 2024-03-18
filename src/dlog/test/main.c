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

#include "t-dlog.c"
#include "t-modpe.c"
#include "t-vec.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(dlog),
    TEST_FUNCTION(dlog_modpe),
    TEST_FUNCTION(dlog_vec)
};

/* main function *************************************************************/

TEST_MAIN(tests)
