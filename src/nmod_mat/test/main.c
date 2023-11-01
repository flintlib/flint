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
#include "t-can_solve.c"
#include "t-charpoly_berkowitz.c"
#include "t-charpoly.c"
#include "t-charpoly_danilevsky.c"
#include "t-concat_horizontal.c"
#include "t-concat_vertical.c"
#include "t-det.c"
#include "t-det_howell.c"
#include "t-howell_form.c"
#include "t-init_clear.c"
#include "t-inv.c"
#include "t-invert_rows_cols.c"
#include "t-lu_classical.c"
#include "t-lu_classical_delayed.c"
#include "t-lu_recursive.c"
#include "t-minpoly.c"
#include "t-mul_blas.c"
#include "t-mul.c"
#include "t-mul_classical_threaded.c"
#include "t-mul_nmod_vec.c"
#include "t-mul_strassen.c"
#include "t-neg.c"
#include "t-nmod_vec_mul.c"
#include "t-nullspace.c"
#include "t-permute_rows.c"
#include "t-pow.c"
#include "t-rank.c"
#include "t-rref.c"
#include "t-scalar_addmul_ui.c"
#include "t-scalar_mul.c"
#include "t-solve.c"
#include "t-solve_tril.c"
#include "t-solve_tril_classical.c"
#include "t-solve_tril_recursive.c"
#include "t-solve_triu.c"
#include "t-solve_triu_classical.c"
#include "t-solve_triu_recursive.c"
#include "t-solve_vec.c"
#include "t-submul.c"
#include "t-trace.c"
#include "t-transpose.c"
#include "t-window_init_clear.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(nmod_mat_add),
    TEST_FUNCTION(nmod_mat_addmul),
    TEST_FUNCTION(nmod_mat_can_solve),
    TEST_FUNCTION(nmod_mat_charpoly_berkowitz),
    TEST_FUNCTION(nmod_mat_charpoly),
    TEST_FUNCTION(nmod_mat_charpoly_danilevsky),
    TEST_FUNCTION(nmod_mat_concat_horizontal),
    TEST_FUNCTION(nmod_mat_concat_vertical),
    TEST_FUNCTION(nmod_mat_det),
    TEST_FUNCTION(nmod_mat_det_howell),
    TEST_FUNCTION(nmod_mat_howell_form),
    TEST_FUNCTION(nmod_mat_init_clear),
    TEST_FUNCTION(nmod_mat_inv),
    TEST_FUNCTION(nmod_mat_invert_rows_cols),
    TEST_FUNCTION(nmod_mat_lu_classical),
    TEST_FUNCTION(nmod_mat_lu_classical_delayed),
    TEST_FUNCTION(nmod_mat_lu_recursive),
    TEST_FUNCTION(nmod_mat_minpoly),
    TEST_FUNCTION(nmod_mat_mul_blas),
    TEST_FUNCTION(nmod_mat_mul),
    TEST_FUNCTION(nmod_mat_mul_classical_threaded),
    TEST_FUNCTION(nmod_mat_mul_nmod_vec),
    TEST_FUNCTION(nmod_mat_mul_strassen),
    TEST_FUNCTION(nmod_mat_neg),
    TEST_FUNCTION(nmod_mat_nmod_vec_mul),
    TEST_FUNCTION(nmod_mat_nullspace),
    TEST_FUNCTION(nmod_mat_permute_rows),
    TEST_FUNCTION(nmod_mat_pow),
    TEST_FUNCTION(nmod_mat_rank),
    TEST_FUNCTION(nmod_mat_rref),
    TEST_FUNCTION(nmod_mat_scalar_addmul_ui),
    TEST_FUNCTION(nmod_mat_scalar_mul),
    TEST_FUNCTION(nmod_mat_solve),
    TEST_FUNCTION(nmod_mat_solve_tril),
    TEST_FUNCTION(nmod_mat_solve_tril_classical),
    TEST_FUNCTION(nmod_mat_solve_tril_recursive),
    TEST_FUNCTION(nmod_mat_solve_triu),
    TEST_FUNCTION(nmod_mat_solve_triu_classical),
    TEST_FUNCTION(nmod_mat_solve_triu_recursive),
    TEST_FUNCTION(nmod_mat_solve_vec),
    TEST_FUNCTION(nmod_mat_submul),
    TEST_FUNCTION(nmod_mat_trace),
    TEST_FUNCTION(nmod_mat_transpose),
    TEST_FUNCTION(nmod_mat_window_init_clear)
};

/* main function *************************************************************/

TEST_MAIN(tests)
