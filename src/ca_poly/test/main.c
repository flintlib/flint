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

#include "t-compose.c"
#include "t-divrem.c"
#include "t-div_series.c"
#include "t-evaluate.c"
#include "t-evaluate_horner.c"
#include "t-exp_series.c"
#include "t-factor_squarefree.c"
#include "t-gcd.c"
#include "t-gcd_euclidean.c"
#include "t-inv_series.c"
#include "t-log_series.c"
#include "t-mul.c"
#include "t-pow_ui.c"
#include "t-roots.c"
#include "t-squarefree_part.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(ca_poly_compose),
    TEST_FUNCTION(ca_poly_divrem),
    TEST_FUNCTION(ca_poly_div_series),
    TEST_FUNCTION(ca_poly_evaluate),
    TEST_FUNCTION(ca_poly_evaluate_horner),
    TEST_FUNCTION(ca_poly_exp_series),
    TEST_FUNCTION(ca_poly_factor_squarefree),
    TEST_FUNCTION(ca_poly_gcd),
    TEST_FUNCTION(ca_poly_gcd_euclidean),
    TEST_FUNCTION(ca_poly_inv_series),
    TEST_FUNCTION(ca_poly_log_series),
    TEST_FUNCTION(ca_poly_mul),
    TEST_FUNCTION(ca_poly_pow_ui),
    TEST_FUNCTION(ca_poly_roots),
    TEST_FUNCTION(ca_poly_squarefree_part)
};

/* main function *************************************************************/

TEST_MAIN(tests)
