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

#include "t-factor_berlekamp.c"
#include "t-factor.c"
#include "t-factor_cantor_zassenhaus.c"
#include "t-factor_distinct_deg.c"
#include "t-factor_distinct_deg_threaded.c"
#include "t-factor_kaltofen_shoup.c"
#include "t-factor_squarefree.c"
#include "t-interval_threaded.c"
#include "t-is_irreducible.c"
#include "t-is_irreducible_ddf.c"
#include "t-is_irreducible_rabin.c"
#include "t-is_squarefree.c"
#include "t-roots.c"
#include "t-roots_factored.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(nmod_poly_factor_berlekamp),
    TEST_FUNCTION(nmod_poly_factor),
    TEST_FUNCTION(nmod_poly_factor_cantor_zassenhaus),
    TEST_FUNCTION(nmod_poly_factor_distinct_deg),
    TEST_FUNCTION(nmod_poly_factor_distinct_deg_threaded),
    TEST_FUNCTION(nmod_poly_factor_kaltofen_shoup),
    TEST_FUNCTION(nmod_poly_factor_squarefree),
    TEST_FUNCTION(nmod_poly_factor_interval_threaded),
    TEST_FUNCTION(nmod_poly_factor_is_irreducible),
    TEST_FUNCTION(nmod_poly_factor_is_irreducible_ddf),
    TEST_FUNCTION(nmod_poly_factor_is_irreducible_rabin),
    TEST_FUNCTION(nmod_poly_factor_is_squarefree),
    TEST_FUNCTION(nmod_poly_factor_roots),
    TEST_FUNCTION(nmod_poly_factor_roots_factored)
};

/* main function *************************************************************/

TEST_MAIN(tests)
