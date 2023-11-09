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

#include "t-max_degrees_tight.c"
#include "t-max_fields.c"
#include "t-monomial_halves.c"
#include "t-pack_unpack.c"
#include "t-pack_unpack_tight.c"
#include "t-search_monomials.c"
#include "t-test_irreducible.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(mpoly_max_degrees_tight),
    TEST_FUNCTION(mpoly_max_fields),
    TEST_FUNCTION(mpoly_monomial_halves),
    TEST_FUNCTION(mpoly_pack_unpack),
    TEST_FUNCTION(mpoly_pack_unpack_tight),
    TEST_FUNCTION(mpoly_search_monomials),
    TEST_FUNCTION(mpoly_test_irreducible)
};

/* main function *************************************************************/

TEST_MAIN(tests)
