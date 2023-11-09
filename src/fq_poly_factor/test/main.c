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
#include "t-factor_equal_deg_prob.c"
#include "t-factor_kaltofen_shoup.c"
#include "t-factor_split_single.c"
#include "t-factor_squarefree.c"
#include "t-is_irreducible_ben_or.c"
#include "t-is_irreducible.c"
#include "t-is_irreducible_ddf.c"
#include "t-is_squarefree.c"
#include "t-iterated_frobenius_preinv.c"
#include "t-roots.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(fq_poly_factor_berlekamp),
    TEST_FUNCTION(fq_poly_factor),
    TEST_FUNCTION(fq_poly_factor_cantor_zassenhaus),
    TEST_FUNCTION(fq_poly_factor_distinct_deg),
    TEST_FUNCTION(fq_poly_factor_equal_deg_prob),
    TEST_FUNCTION(fq_poly_factor_kaltofen_shoup),
    TEST_FUNCTION(fq_poly_factor_split_single),
    TEST_FUNCTION(fq_poly_factor_squarefree),
    TEST_FUNCTION(fq_poly_factor_is_irreducible_ben_or),
    TEST_FUNCTION(fq_poly_factor_is_irreducible),
    TEST_FUNCTION(fq_poly_factor_is_irreducible_ddf),
    TEST_FUNCTION(fq_poly_factor_is_squarefree),
    TEST_FUNCTION(fq_poly_factor_iterated_frobenius_preinv),
    TEST_FUNCTION(fq_poly_factor_roots)
};

/* main function *************************************************************/

TEST_MAIN(tests)
