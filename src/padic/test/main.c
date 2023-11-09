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
#include "t-div.c"
#include "t-exp_balanced.c"
#include "t-exp.c"
#include "t-exp_rectangular.c"
#include "t-get_set_fmpz.c"
#include "t-get_set_mpq.c"
#include "t-get_set_mpz.c"
#include "t-get_str.c"
#include "t-inv.c"
#include "t-log_balanced.c"
#include "t-log.c"
#include "t-log_rectangular.c"
#include "t-log_satoh.c"
#include "t-mul.c"
#include "t-neg.c"
#include "t-pow_si.c"
#include "t-randtest.c"
#include "t-shift.c"
#include "t-sqrt.c"
#include "t-sub.c"
#include "t-teichmuller.c"
#include "t-val_fac.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(padic_add),
    TEST_FUNCTION(padic_div),
    TEST_FUNCTION(padic_exp_balanced),
    TEST_FUNCTION(padic_exp),
    TEST_FUNCTION(padic_exp_rectangular),
    TEST_FUNCTION(padic_get_set_fmpz),
    TEST_FUNCTION(padic_get_set_mpq),
    TEST_FUNCTION(padic_get_set_mpz),
    TEST_FUNCTION(padic_get_str),
    TEST_FUNCTION(padic_inv),
    TEST_FUNCTION(padic_log_balanced),
    TEST_FUNCTION(padic_log),
    TEST_FUNCTION(padic_log_rectangular),
    TEST_FUNCTION(padic_log_satoh),
    TEST_FUNCTION(padic_mul),
    TEST_FUNCTION(padic_neg),
    TEST_FUNCTION(padic_pow_si),
    TEST_FUNCTION(padic_randtest),
    TEST_FUNCTION(padic_shift),
    TEST_FUNCTION(padic_sqrt),
    TEST_FUNCTION(padic_sub),
    TEST_FUNCTION(padic_teichmuller),
    TEST_FUNCTION(padic_val_fac)
};

/* main function *************************************************************/

TEST_MAIN(tests)
