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
#include "t-addmul.c"
#include "t-all.c"
#include "t-derivative.c"
#include "t-div.c"
#include "t-evaluate_fmpq.c"
#include "t-init_clear.c"
#include "t-inv.c"
#include "t-mul.c"
#include "t-neg.c"
#include "t-pow.c"
#include "t-scalar_div_fmpq.c"
#include "t-scalar_div_fmpz.c"
#include "t-scalar_div_si.c"
#include "t-scalar_mul_fmpq.c"
#include "t-scalar_mul_fmpz.c"
#include "t-scalar_mul_si.c"
#include "t-set_equal.c"
#include "t-set_si_equal.c"
#include "t-sub.c"
#include "t-submul.c"
#include "t-swap.c"
#include "t-zero.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(fmpz_poly_q_add),
    TEST_FUNCTION(fmpz_poly_q_addmul),
    TEST_FUNCTION(fmpz_poly_q_all),
    TEST_FUNCTION(fmpz_poly_q_derivative),
    TEST_FUNCTION(fmpz_poly_q_div),
    TEST_FUNCTION(fmpz_poly_q_evaluate_fmpq),
    TEST_FUNCTION(fmpz_poly_q_init_clear),
    TEST_FUNCTION(fmpz_poly_q_inv),
    TEST_FUNCTION(fmpz_poly_q_mul),
    TEST_FUNCTION(fmpz_poly_q_neg),
    TEST_FUNCTION(fmpz_poly_q_pow),
    TEST_FUNCTION(fmpz_poly_q_scalar_div_fmpq),
    TEST_FUNCTION(fmpz_poly_q_scalar_div_fmpz),
    TEST_FUNCTION(fmpz_poly_q_scalar_div_si),
    TEST_FUNCTION(fmpz_poly_q_scalar_mul_fmpq),
    TEST_FUNCTION(fmpz_poly_q_scalar_mul_fmpz),
    TEST_FUNCTION(fmpz_poly_q_scalar_mul_si),
    TEST_FUNCTION(fmpz_poly_q_set_equal),
    TEST_FUNCTION(fmpz_poly_q_set_si_equal),
    TEST_FUNCTION(fmpz_poly_q_sub),
    TEST_FUNCTION(fmpz_poly_q_submul),
    TEST_FUNCTION(fmpz_poly_q_swap),
    TEST_FUNCTION(fmpz_poly_q_zero)
};

/* main function *************************************************************/

TEST_MAIN(tests)
