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

#include "t-add.c"
#include "t-exp_balanced.c"
#include "t-exp.c"
#include "t-exp_rectangular.c"
#include "t-frobenius.c"
#include "t-inv.c"
#include "t-log_balanced.c"
#include "t-log.c"
#include "t-log_rectangular.c"
#include "t-mul.c"
#include "t-neg.c"
#include "t-norm_analytic.c"
#include "t-norm.c"
#include "t-norm_resultant.c"
#include "t-pow.c"
#include "t-sqrt.c"
#include "t-sub.c"
#include "t-teichmuller.c"
#include "t-trace.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(qadic_add),
    TEST_FUNCTION(qadic_exp_balanced),
    TEST_FUNCTION(qadic_exp),
    TEST_FUNCTION(qadic_exp_rectangular),
    TEST_FUNCTION(qadic_frobenius),
    TEST_FUNCTION(qadic_inv),
    TEST_FUNCTION(qadic_log_balanced),
    TEST_FUNCTION(qadic_log),
    TEST_FUNCTION(qadic_log_rectangular),
    TEST_FUNCTION(qadic_mul),
    TEST_FUNCTION(qadic_neg),
    TEST_FUNCTION(qadic_norm_analytic),
    TEST_FUNCTION(qadic_norm),
    TEST_FUNCTION(qadic_norm_resultant),
    TEST_FUNCTION(qadic_pow),
    TEST_FUNCTION(qadic_sqrt),
    TEST_FUNCTION(qadic_sub),
    TEST_FUNCTION(qadic_teichmuller),
    TEST_FUNCTION(qadic_trace)
};

/* main function *************************************************************/

TEST_MAIN(tests)
