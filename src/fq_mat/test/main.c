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

#include "t-add_sub.c"
#include "t-can_solve.c"
#include "t-charpoly.c"
#include "t-concat_horizontal.c"
#include "t-concat_vertical.c"
#include "t-equal.c"
#include "t-inv.c"
#include "t-invert_rows_cols.c"
#include "t-is_zero.c"
#include "t-lu_classical.c"
#include "t-lu_recursive.c"
#include "t-minpoly.c"
#include "t-mul.c"
#include "t-mul_KS.c"
#include "t-mul_vec.c"
#include "t-nullspace.c"
#include "t-one.c"
#include "t-rank.c"
#include "t-rref.c"
#include "t-set_fmpz_mod_mat.c"
#include "t-set_nmod_mat.c"
#include "t-solve.c"
#include "t-solve_tril.c"
#include "t-solve_tril_classical.c"
#include "t-solve_tril_recursive.c"
#include "t-solve_triu.c"
#include "t-solve_triu_classical.c"
#include "t-solve_triu_recursive.c"
#include "t-submul.c"
#include "t-vec_mul.c"
#include "t-window_init_clear.c"
#include "t-zero.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(fq_mat_add_sub),
    TEST_FUNCTION(fq_mat_can_solve),
    TEST_FUNCTION(fq_mat_charpoly),
    TEST_FUNCTION(fq_mat_concat_horizontal),
    TEST_FUNCTION(fq_mat_concat_vertical),
    TEST_FUNCTION(fq_mat_equal),
    TEST_FUNCTION(fq_mat_inv),
    TEST_FUNCTION(fq_mat_invert_rows_cols),
    TEST_FUNCTION(fq_mat_is_zero),
    TEST_FUNCTION(fq_mat_lu_classical),
    TEST_FUNCTION(fq_mat_lu_recursive),
    TEST_FUNCTION(fq_mat_minpoly),
    TEST_FUNCTION(fq_mat_mul),
    TEST_FUNCTION(fq_mat_mul_KS),
    TEST_FUNCTION(fq_mat_mul_vec),
    TEST_FUNCTION(fq_mat_nullspace),
    TEST_FUNCTION(fq_mat_one),
    TEST_FUNCTION(fq_mat_rank),
    TEST_FUNCTION(fq_mat_rref),
    TEST_FUNCTION(fq_mat_set_fmpz_mod_mat),
    TEST_FUNCTION(fq_mat_set_nmod_mat),
    TEST_FUNCTION(fq_mat_solve),
    TEST_FUNCTION(fq_mat_solve_tril),
    TEST_FUNCTION(fq_mat_solve_tril_classical),
    TEST_FUNCTION(fq_mat_solve_tril_recursive),
    TEST_FUNCTION(fq_mat_solve_triu),
    TEST_FUNCTION(fq_mat_solve_triu_classical),
    TEST_FUNCTION(fq_mat_solve_triu_recursive),
    TEST_FUNCTION(fq_mat_submul),
    TEST_FUNCTION(fq_mat_vec_mul),
    TEST_FUNCTION(fq_mat_window_init_clear),
    TEST_FUNCTION(fq_mat_zero)
};

/* main function *************************************************************/

TEST_MAIN(tests)
