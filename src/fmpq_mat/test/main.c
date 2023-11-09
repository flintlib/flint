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
#include "t-can_solve.c"
#include "t-can_solve_dixon.c"
#include "t-can_solve_fraction_free.c"
#include "t-can_solve_multi_mod.c"
#include "t-charpoly.c"
#include "t-concat_horizontal.c"
#include "t-concat_vertical.c"
#include "t-det.c"
#include "t-fmpq_vec_mul.c"
#include "t-fmpz_vec_mul.c"
#include "t-gso.c"
#include "t-init_clear.c"
#include "t-inv.c"
#include "t-invert_rows_cols.c"
#include "t-is_integral.c"
#include "t-is_one.c"
#include "t-kronecker_product.c"
#include "t-minpoly.c"
#include "t-mul.c"
#include "t-mul_fmpq_vec.c"
#include "t-mul_fmpz_vec.c"
#include "t-neg.c"
#include "t-one.c"
#include "t-rref.c"
#include "t-scalar_div_fmpz.c"
#include "t-scalar_mul_fmpq.c"
#include "t-scalar_mul_fmpz.c"
#include "t-solve.c"
#include "t-solve_dixon.c"
#include "t-solve_fmpz_mat.c"
#include "t-solve_fmpz_mat_dixon.c"
#include "t-solve_fmpz_mat_fraction_free.c"
#include "t-solve_fmpz_mat_multi_mod.c"
#include "t-solve_fraction_free.c"
#include "t-solve_multi_mod.c"
#include "t-sub.c"
#include "t-trace.c"
#include "t-transpose.c"
#include "t-window_init_clear.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(fmpq_mat_add),
    TEST_FUNCTION(fmpq_mat_can_solve),
    TEST_FUNCTION(fmpq_mat_can_solve_dixon),
    TEST_FUNCTION(fmpq_mat_can_solve_fraction_free),
    TEST_FUNCTION(fmpq_mat_can_solve_multi_mod),
    TEST_FUNCTION(fmpq_mat_charpoly),
    TEST_FUNCTION(fmpq_mat_concat_horizontal),
    TEST_FUNCTION(fmpq_mat_concat_vertical),
    TEST_FUNCTION(fmpq_mat_det),
    TEST_FUNCTION(fmpq_mat_fmpq_vec_mul),
    TEST_FUNCTION(fmpq_mat_fmpz_vec_mul),
    TEST_FUNCTION(fmpq_mat_gso),
    TEST_FUNCTION(fmpq_mat_init_clear),
    TEST_FUNCTION(fmpq_mat_inv),
    TEST_FUNCTION(fmpq_mat_invert_rows_cols),
    TEST_FUNCTION(fmpq_mat_is_integral),
    TEST_FUNCTION(fmpq_mat_is_one),
    TEST_FUNCTION(fmpq_mat_kronecker_product),
    TEST_FUNCTION(fmpq_mat_minpoly),
    TEST_FUNCTION(fmpq_mat_mul),
    TEST_FUNCTION(fmpq_mat_mul_fmpq_vec),
    TEST_FUNCTION(fmpq_mat_mul_fmpz_vec),
    TEST_FUNCTION(fmpq_mat_neg),
    TEST_FUNCTION(fmpq_mat_one),
    TEST_FUNCTION(fmpq_mat_rref),
    TEST_FUNCTION(fmpq_mat_scalar_div_fmpz),
    TEST_FUNCTION(fmpq_mat_scalar_mul_fmpq),
    TEST_FUNCTION(fmpq_mat_scalar_mul_fmpz),
    TEST_FUNCTION(fmpq_mat_solve),
    TEST_FUNCTION(fmpq_mat_solve_dixon),
    TEST_FUNCTION(fmpq_mat_solve_fmpz_mat),
    TEST_FUNCTION(fmpq_mat_solve_fmpz_mat_dixon),
    TEST_FUNCTION(fmpq_mat_solve_fmpz_mat_fraction_free),
    TEST_FUNCTION(fmpq_mat_solve_fmpz_mat_multi_mod),
    TEST_FUNCTION(fmpq_mat_solve_fraction_free),
    TEST_FUNCTION(fmpq_mat_solve_multi_mod),
    TEST_FUNCTION(fmpq_mat_sub),
    TEST_FUNCTION(fmpq_mat_trace),
    TEST_FUNCTION(fmpq_mat_transpose),
    TEST_FUNCTION(fmpq_mat_window_init_clear)
};

/* main function *************************************************************/

TEST_MAIN(tests)
