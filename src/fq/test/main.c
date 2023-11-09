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
#include "t-assign.c"
#include "t-ctx_init.c"
#include "t-div.c"
#include "t-frobenius.c"
#include "t-get_set_fmpz.c"
#include "t-get_set_fmpz_mod_poly.c"
#include "t-get_set_fmpz_poly.c"
#include "t-inv.c"
#include "t-is_invertible.c"
#include "t-is_invertible_f.c"
#include "t-is_primitive.c"
#include "t-is_square.c"
#include "t-mul.c"
#include "t-mul_fmpz.c"
#include "t-mul_si.c"
#include "t-multiplicative_order.c"
#include "t-mul_ui.c"
#include "t-neg.c"
#include "t-norm.c"
#include "t-pow.c"
#include "t-pth_root.c"
#include "t-sqr.c"
#include "t-sqrt.c"
#include "t-sub.c"
#include "t-trace.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(fq_add),
    TEST_FUNCTION(fq_assign),
    TEST_FUNCTION(fq_ctx_init),
    TEST_FUNCTION(fq_div),
    TEST_FUNCTION(fq_frobenius),
    TEST_FUNCTION(fq_get_set_fmpz),
    TEST_FUNCTION(fq_get_set_fmpz_mod_poly),
    TEST_FUNCTION(fq_get_set_fmpz_poly),
    TEST_FUNCTION(fq_inv),
    TEST_FUNCTION(fq_is_invertible),
    TEST_FUNCTION(fq_is_invertible_f),
    TEST_FUNCTION(fq_is_primitive),
    TEST_FUNCTION(fq_is_square),
    TEST_FUNCTION(fq_mul),
    TEST_FUNCTION(fq_mul_fmpz),
    TEST_FUNCTION(fq_mul_si),
    TEST_FUNCTION(fq_multiplicative_order),
    TEST_FUNCTION(fq_mul_ui),
    TEST_FUNCTION(fq_neg),
    TEST_FUNCTION(fq_norm),
    TEST_FUNCTION(fq_pow),
    TEST_FUNCTION(fq_pth_root),
    TEST_FUNCTION(fq_sqr),
    TEST_FUNCTION(fq_sqrt),
    TEST_FUNCTION(fq_sub),
    TEST_FUNCTION(fq_trace)
};

/* main function *************************************************************/

TEST_MAIN(tests)
