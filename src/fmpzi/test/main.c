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

#include "t-add_sub.c"
#include "t-divexact.c"
#include "t-divrem_approx.c"
#include "t-divrem.c"
#include "t-gcd_binary.c"
#include "t-gcd.c"
#include "t-gcd_euclidean.c"
#include "t-gcd_euclidean_improved.c"
#include "t-gcd_shortest.c"
#include "t-is_prime.c"
#include "t-is_probabprime.c"
#include "t-mul.c"
#include "t-pow_ui.c"
#include "t-remove_one_plus_i.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(fmpzi_add_sub),
    TEST_FUNCTION(fmpzi_divexact),
    TEST_FUNCTION(fmpzi_divrem_approx),
    TEST_FUNCTION(fmpzi_divrem),
    TEST_FUNCTION(fmpzi_gcd_binary),
    TEST_FUNCTION(fmpzi_gcd),
    TEST_FUNCTION(fmpzi_gcd_euclidean),
    TEST_FUNCTION(fmpzi_gcd_euclidean_improved),
    TEST_FUNCTION(fmpzi_gcd_shortest),
    TEST_FUNCTION(fmpzi_is_prime),
    TEST_FUNCTION(fmpzi_is_probabprime),
    TEST_FUNCTION(fmpzi_mul),
    TEST_FUNCTION(fmpzi_pow_ui),
    TEST_FUNCTION(fmpzi_remove_one_plus_i)
};

/* main function *************************************************************/

TEST_MAIN(tests)
