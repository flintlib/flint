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

#include "t-approx_eig_qr.c"
#include "t-charpoly.c"
#include "t-companion.c"
#include "t-det.c"
#include "t-det_precond.c"
#include "t-dft.c"
#include "t-eig_enclosure_rump.c"
#include "t-eig_global_enclosure.c"
#include "t-eig_multiple.c"
#include "t-eig_simple.c"
#include "t-exp.c"
#include "t-exp_taylor_sum.c"
#include "t-frobenius_norm.c"
#include "t-inv.c"
#include "t-lu.c"
#include "t-lu_recursive.c"
#include "t-mul.c"
#include "t-mul_entrywise.c"
#include "t-mul_reorder.c"
#include "t-mul_threaded.c"
#include "t-set_real_imag.c"
#include "t-solve.c"
#include "t-solve_lu.c"
#include "t-solve_precond.c"
#include "t-solve_tril.c"
#include "t-solve_triu.c"
#include "t-sqr.c"
#include "t-trace.c"
#include "t-transpose.c"
#include "t-vector_mul.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(acb_mat_approx_eig_qr),
    TEST_FUNCTION(acb_mat_charpoly),
    TEST_FUNCTION(acb_mat_companion),
    TEST_FUNCTION(acb_mat_det),
    TEST_FUNCTION(acb_mat_det_precond),
    TEST_FUNCTION(acb_mat_dft),
    TEST_FUNCTION(acb_mat_eig_enclosure_rump),
    TEST_FUNCTION(acb_mat_eig_global_enclosure),
    TEST_FUNCTION(acb_mat_eig_multiple),
    TEST_FUNCTION(acb_mat_eig_simple),
    TEST_FUNCTION(acb_mat_exp),
    TEST_FUNCTION(acb_mat_exp_taylor_sum),
    TEST_FUNCTION(acb_mat_frobenius_norm),
    TEST_FUNCTION(acb_mat_inv),
    TEST_FUNCTION(acb_mat_lu),
    TEST_FUNCTION(acb_mat_lu_recursive),
    TEST_FUNCTION(acb_mat_mul),
    TEST_FUNCTION(acb_mat_mul_entrywise),
    TEST_FUNCTION(acb_mat_mul_reorder),
    TEST_FUNCTION(acb_mat_mul_threaded),
    TEST_FUNCTION(acb_mat_set_real_imag),
    TEST_FUNCTION(acb_mat_solve),
    TEST_FUNCTION(acb_mat_solve_lu),
    TEST_FUNCTION(acb_mat_solve_precond),
    TEST_FUNCTION(acb_mat_solve_tril),
    TEST_FUNCTION(acb_mat_solve_triu),
    TEST_FUNCTION(acb_mat_sqr),
    TEST_FUNCTION(acb_mat_trace),
    TEST_FUNCTION(acb_mat_transpose),
    TEST_FUNCTION(acb_mat_vector_mul)
};

/* main function *************************************************************/

TEST_MAIN(tests)
