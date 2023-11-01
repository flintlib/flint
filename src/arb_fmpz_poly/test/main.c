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

#include "t-complex_roots.c"
#include "t-evaluate_acb.c"
#include "t-evaluate_arb.c"
#include "t-gauss_period_minpoly.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(arb_fmpz_poly_complex_roots),
    TEST_FUNCTION(arb_fmpz_poly_evaluate_acb),
    TEST_FUNCTION(arb_fmpz_poly_evaluate_arb),
    TEST_FUNCTION(arb_fmpz_poly_gauss_period_minpoly)
};

/* main function *************************************************************/

TEST_MAIN(tests)
