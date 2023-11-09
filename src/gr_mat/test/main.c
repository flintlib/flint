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

#include "t-adjugate.c"
#include "t-charpoly_danilevsky.c"
#include "t-charpoly_faddeev_bsgs.c"
#include "t-charpoly_faddeev.c"
#include "t-charpoly_gauss.c"
#include "t-charpoly_householder.c"
#include "t-concat_horizontal.c"
#include "t-concat_vertical.c"
#include "t-det_berkowitz.c"
#include "t-det_cofactor.c"
#include "t-det_fflu.c"
#include "t-det_lu.c"
#include "t-diagonalization.c"
#include "t-hadamard.c"
#include "t-hessenberg.c"
#include "t-hessenberg_gauss.c"
#include "t-hessenberg_householder.c"
#include "t-inv.c"
#include "t-invert_rows_cols.c"
#include "t-lu.c"
#include "t-lu_classical.c"
#include "t-lu_recursive.c"
#include "t-minpoly_field.c"
#include "t-mul_strassen.c"
#include "t-nullspace.c"
#include "t-properties.c"
#include "t-randrank.c"
#include "t-rank.c"
#include "t-rank_fflu.c"
#include "t-rank_lu.c"
#include "t-rref_den_fflu.c"
#include "t-rref_fflu.c"
#include "t-rref_lu.c"
#include "t-solve.c"
#include "t-solve_den.c"
#include "t-solve_den_fflu.c"
#include "t-solve_fflu.c"
#include "t-solve_field.c"
#include "t-solve_lu.c"
#include "t-solve_tril.c"
#include "t-solve_triu.c"
#include "t-window_init_clear.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(gr_mat_adjugate),
    TEST_FUNCTION(gr_mat_charpoly_danilevsky),
    TEST_FUNCTION(gr_mat_charpoly_faddeev_bsgs),
    TEST_FUNCTION(gr_mat_charpoly_faddeev),
    TEST_FUNCTION(gr_mat_charpoly_gauss),
    TEST_FUNCTION(gr_mat_charpoly_householder),
    TEST_FUNCTION(gr_mat_concat_horizontal),
    TEST_FUNCTION(gr_mat_concat_vertical),
    TEST_FUNCTION(gr_mat_det_berkowitz),
    TEST_FUNCTION(gr_mat_det_cofactor),
    TEST_FUNCTION(gr_mat_det_fflu),
    TEST_FUNCTION(gr_mat_det_lu),
    TEST_FUNCTION(gr_mat_diagonalization),
    TEST_FUNCTION(gr_mat_hadamard),
    TEST_FUNCTION(gr_mat_hessenberg),
    TEST_FUNCTION(gr_mat_hessenberg_gauss),
    TEST_FUNCTION(gr_mat_hessenberg_householder),
    TEST_FUNCTION(gr_mat_inv),
    TEST_FUNCTION(gr_mat_invert_rows_cols),
    TEST_FUNCTION(gr_mat_lu),
    TEST_FUNCTION(gr_mat_lu_classical),
    TEST_FUNCTION(gr_mat_lu_recursive),
    TEST_FUNCTION(gr_mat_minpoly_field),
    TEST_FUNCTION(gr_mat_mul_strassen),
    TEST_FUNCTION(gr_mat_nullspace),
    TEST_FUNCTION(gr_mat_properties),
    TEST_FUNCTION(gr_mat_randrank),
    TEST_FUNCTION(gr_mat_rank),
    TEST_FUNCTION(gr_mat_rank_fflu),
    TEST_FUNCTION(gr_mat_rank_lu),
    TEST_FUNCTION(gr_mat_rref_den_fflu),
    TEST_FUNCTION(gr_mat_rref_fflu),
    TEST_FUNCTION(gr_mat_rref_lu),
    TEST_FUNCTION(gr_mat_solve),
    TEST_FUNCTION(gr_mat_solve_den),
    TEST_FUNCTION(gr_mat_solve_den_fflu),
    TEST_FUNCTION(gr_mat_solve_fflu),
    TEST_FUNCTION(gr_mat_solve_field),
    TEST_FUNCTION(gr_mat_solve_lu),
    TEST_FUNCTION(gr_mat_solve_tril),
    TEST_FUNCTION(gr_mat_solve_triu),
    TEST_FUNCTION(gr_mat_window_init_clear)
};

/* main function *************************************************************/

TEST_MAIN(tests)
