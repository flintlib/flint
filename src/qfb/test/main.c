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

#include "t-exponent.c"
#include "t-exponent_element.c"
#include "t-exponent_grh.c"
#include "t-inverse.c"
#include "t-nucomp.c"
#include "t-nudupl.c"
#include "t-pow.c"
#include "t-pow_ui.c"
#include "t-prime_form.c"
#include "t-reduce.c"
#include "t-reduced_forms.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(qfb_exponent),
    TEST_FUNCTION(qfb_exponent_element),
    TEST_FUNCTION(qfb_exponent_grh),
    TEST_FUNCTION(qfb_inverse),
    TEST_FUNCTION(qfb_nucomp),
    TEST_FUNCTION(qfb_nudupl),
    TEST_FUNCTION(qfb_pow),
    TEST_FUNCTION(qfb_pow_ui),
    TEST_FUNCTION(qfb_prime_form),
    TEST_FUNCTION(qfb_reduce),
    TEST_FUNCTION(qfb_reduced_forms)
};

/* main function *************************************************************/

TEST_MAIN(tests)
