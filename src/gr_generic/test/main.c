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

#include "fmpz_poly.h"
#include "fmpz_mpoly.h"
#include "gr_generic.h"

/* Include functions *********************************************************/

#include "t-fmpz_mpoly_evaluate.c"
#include "t-fmpz_poly_evaluate.c"
#include "t-pow.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(gr_generic_fmpz_mpoly_evaluate),
    TEST_FUNCTION(gr_generic_fmpz_poly_evaluate),
    TEST_FUNCTION(gr_generic_pow)
};

/* main function *************************************************************/

TEST_MAIN(tests)
