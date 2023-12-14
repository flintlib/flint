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

#include "t-ecm.c"
#include "t-factor.c"
#include "t-factor_pp1.c"
#include "t-factor_smooth.c"
#include "t-factor_trial.c"
#include "t-io.c"
#include "t-pollard_brent.c"
#include "t-pollard_brent_single.c"
#include "t-refine.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(fmpz_factor),
    TEST_FUNCTION(fmpz_factor_ecm),
    TEST_FUNCTION(fmpz_factor_fprint),
    TEST_FUNCTION(fmpz_factor_pp1),
    TEST_FUNCTION(fmpz_factor_smooth),
    TEST_FUNCTION(fmpz_factor_trial),
    TEST_FUNCTION(fmpz_factor_pollard_brent),
    TEST_FUNCTION(fmpz_factor_pollard_brent_single),
    TEST_FUNCTION(fmpz_factor_refine)
};

/* main function *************************************************************/

TEST_MAIN(tests)
