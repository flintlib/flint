/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Include functions *********************************************************/

#include "t-add_sub_n.c"
#include "t-addmul_submul.c"
#include "t-complex_mat_mul.c"
#include "t-mat_mul.c"
#include "t-nfixed_dot.c"
#include "t-nfixed_mat_mul.c"
#include "t-nfixed_mat_mul_classical.c"
#include "t-nfixed_mat_mul_strassen.c"
#include "t-nfixed_mat_mul_waksman.c"
#include "t-nfloat.c"
#include "t-nfloat_complex.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(add_sub_n),
    TEST_FUNCTION(addmul_submul),
    TEST_FUNCTION(complex_mat_mul),
    TEST_FUNCTION(mat_mul),
    TEST_FUNCTION(nfixed_dot),
    TEST_FUNCTION(nfixed_mat_mul),
    TEST_FUNCTION(nfixed_mat_mul_classical),
    TEST_FUNCTION(nfixed_mat_mul_strassen),
    TEST_FUNCTION(nfixed_mat_mul_waksman),
    TEST_FUNCTION(nfloat),
    TEST_FUNCTION(nfloat_complex),
};

/* main function *************************************************************/

TEST_MAIN(tests)
