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

#include "t-factor.c"
#include "t-factor_content.c"
#include "t-factor_squarefree.c"
#include "t-factor_wang.c"
#include "t-factor_zassenhaus.c"
#include "t-factor_zippel.c"
#include "t-gcd_subresultant.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(fq_nmod_mpoly_factor),
    TEST_FUNCTION(fq_nmod_mpoly_factor_content),
    TEST_FUNCTION(fq_nmod_mpoly_factor_squarefree),
    TEST_FUNCTION(fq_nmod_mpoly_factor_wang),
    TEST_FUNCTION(fq_nmod_mpoly_factor_zassenhaus),
    TEST_FUNCTION(fq_nmod_mpoly_factor_zippel),
    TEST_FUNCTION(fq_nmod_mpoly_factor_gcd_subresultant)
};

/* main function *************************************************************/

TEST_MAIN(tests)
