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
#include "t-concat_horizontal.c"
#include "t-concat_vertical.c"
#include "t-det.c"
#include "t-det_interpolate.c"
#include "t-init_clear.c"
#include "t-inv.c"
#include "t-mul.c"
#include "t-mul_KS.c"
#include "t-mullow.c"
#include "t-neg.c"
#include "t-nullspace.c"
#include "t-one.c"
#include "t-pow.c"
#include "t-pow_trunc.c"
#include "t-prod.c"
#include "t-rank.c"
#include "t-rref.c"
#include "t-solve_fflu.c"
#include "t-sqr.c"
#include "t-sqr_KS.c"
#include "t-sqrlow.c"
#include "t-sub.c"
#include "t-trace.c"
#include "t-window_init_clear.c"
#include "t-zero.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(fmpz_poly_mat_add),
    TEST_FUNCTION(fmpz_poly_mat_concat_horizontal),
    TEST_FUNCTION(fmpz_poly_mat_concat_vertical),
    TEST_FUNCTION(fmpz_poly_mat_det),
    TEST_FUNCTION(fmpz_poly_mat_det_interpolate),
    TEST_FUNCTION(fmpz_poly_mat_init_clear),
    TEST_FUNCTION(fmpz_poly_mat_inv),
    TEST_FUNCTION(fmpz_poly_mat_mul),
    TEST_FUNCTION(fmpz_poly_mat_mul_KS),
    TEST_FUNCTION(fmpz_poly_mat_mullow),
    TEST_FUNCTION(fmpz_poly_mat_neg),
    TEST_FUNCTION(fmpz_poly_mat_nullspace),
    TEST_FUNCTION(fmpz_poly_mat_one),
    TEST_FUNCTION(fmpz_poly_mat_pow),
    TEST_FUNCTION(fmpz_poly_mat_pow_trunc),
    TEST_FUNCTION(fmpz_poly_mat_prod),
    TEST_FUNCTION(fmpz_poly_mat_rank),
    TEST_FUNCTION(fmpz_poly_mat_rref),
    TEST_FUNCTION(fmpz_poly_mat_solve_fflu),
    TEST_FUNCTION(fmpz_poly_mat_sqr),
    TEST_FUNCTION(fmpz_poly_mat_sqr_KS),
    TEST_FUNCTION(fmpz_poly_mat_sqrlow),
    TEST_FUNCTION(fmpz_poly_mat_sub),
    TEST_FUNCTION(fmpz_poly_mat_trace),
    TEST_FUNCTION(fmpz_poly_mat_window_init_clear),
    TEST_FUNCTION(fmpz_poly_mat_zero)
};

/* main function *************************************************************/

TEST_MAIN(tests)
