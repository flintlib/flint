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

#include "t-builtins.c"
#include "t-call_vec.c"
#include "t-replace.c"
#include "t-set_fmpz.c"
#include "t-write_latex.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(fexpr_builtins),
    TEST_FUNCTION(fexpr_call_vec),
    TEST_FUNCTION(fexpr_replace),
    TEST_FUNCTION(fexpr_set_fmpz),
    TEST_FUNCTION(fexpr_write_latex)
};

/* main function *************************************************************/

TEST_MAIN(tests)
