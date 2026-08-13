/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Include functions *********************************************************/

#include "t-factor.c"
#include "t-factor_cubic.c"
#include "t-factor_deflation_capelli.c"
#include "t-factor_squarefree.c"
#include "t-inflation_is_irreducible_capelli.c"
#include "t-factor_zassenhaus.c"
#include "t-zassenhaus_subset.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(fmpz_poly_factor),
    TEST_FUNCTION(fmpz_poly_factor_cubic),
    TEST_FUNCTION(fmpz_poly_factor_deflation_capelli),
    TEST_FUNCTION(fmpz_poly_factor_squarefree),
    TEST_FUNCTION(fmpz_poly_factor_inflation_is_irreducible_capelli),
    TEST_FUNCTION(fmpz_poly_factor_zassenhaus),
    TEST_FUNCTION(fmpz_poly_factor_zassenhaus_subset)
};

/* main function *************************************************************/

TEST_MAIN(tests)
