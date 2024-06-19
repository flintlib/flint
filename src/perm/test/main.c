/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Include functions *********************************************************/

#include "t-compose.c"
#include "t-inv.c"
#include "t-parity.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(perm_compose),
    TEST_FUNCTION(perm_inv),
    TEST_FUNCTION(perm_parity)
};

/* main function *************************************************************/

TEST_MAIN(tests)
