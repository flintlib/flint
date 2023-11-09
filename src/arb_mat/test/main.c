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

#include "t-addmul_rad_mag_fast.c"
#include "t-charpoly.c"
#include "t-cho.c"
#include "t-companion.c"
#include "t-dct.c"
#include "t-det.c"
#include "t-det_precond.c"
#include "t-exp.c"
#include "t-exp_taylor_sum.c"
#include "t-frobenius_norm.c"
#include "t-inv.c"
#include "t-inv_cho_precomp.c"
#include "t-inv_ldl_precomp.c"
#include "t-ldl.c"
#include "t-lu.c"
#include "t-lu_recursive.c"
#include "t-mul_block.c"
#include "t-mul.c"
#include "t-mul_entrywise.c"
#include "t-mul_threaded.c"
#include "t-pascal.c"
#include "t-solve.c"
#include "t-solve_cho_precomp.c"
#include "t-solve_ldl_precomp.c"
#include "t-solve_lu.c"
#include "t-solve_preapprox.c"
#include "t-solve_precond.c"
#include "t-solve_tril.c"
#include "t-solve_triu.c"
#include "t-spd_get_fmpz_mat.c"
#include "t-spd_inv.c"
#include "t-spd_lll_reduce.c"
#include "t-spd_solve.c"
#include "t-sqr.c"
#include "t-stirling.c"
#include "t-trace.c"
#include "t-transpose.c"
#include "t-vector_mul.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(arb_mat_addmul_rad_mag_fast),
    TEST_FUNCTION(arb_mat_charpoly),
    TEST_FUNCTION(arb_mat_cho),
    TEST_FUNCTION(arb_mat_companion),
    TEST_FUNCTION(arb_mat_dct),
    TEST_FUNCTION(arb_mat_det),
    TEST_FUNCTION(arb_mat_det_precond),
    TEST_FUNCTION(arb_mat_exp),
    TEST_FUNCTION(arb_mat_exp_taylor_sum),
    TEST_FUNCTION(arb_mat_frobenius_norm),
    TEST_FUNCTION(arb_mat_inv),
    TEST_FUNCTION(arb_mat_inv_cho_precomp),
    TEST_FUNCTION(arb_mat_inv_ldl_precomp),
    TEST_FUNCTION(arb_mat_ldl),
    TEST_FUNCTION(arb_mat_lu),
    TEST_FUNCTION(arb_mat_lu_recursive),
    TEST_FUNCTION(arb_mat_mul_block),
    TEST_FUNCTION(arb_mat_mul),
    TEST_FUNCTION(arb_mat_mul_entrywise),
    TEST_FUNCTION(arb_mat_mul_threaded),
    TEST_FUNCTION(arb_mat_pascal),
    TEST_FUNCTION(arb_mat_solve),
    TEST_FUNCTION(arb_mat_solve_cho_precomp),
    TEST_FUNCTION(arb_mat_solve_ldl_precomp),
    TEST_FUNCTION(arb_mat_solve_lu),
    TEST_FUNCTION(arb_mat_solve_preapprox),
    TEST_FUNCTION(arb_mat_solve_precond),
    TEST_FUNCTION(arb_mat_solve_tril),
    TEST_FUNCTION(arb_mat_solve_triu),
    TEST_FUNCTION(arb_mat_spd_get_fmpz_mat),
    TEST_FUNCTION(arb_mat_spd_inv),
    TEST_FUNCTION(arb_mat_spd_lll_reduce),
    TEST_FUNCTION(arb_mat_spd_solve),
    TEST_FUNCTION(arb_mat_sqr),
    TEST_FUNCTION(arb_mat_stirling),
    TEST_FUNCTION(arb_mat_trace),
    TEST_FUNCTION(arb_mat_transpose),
    TEST_FUNCTION(arb_mat_vector_mul)
};

/* main function *************************************************************/

TEST_MAIN(tests)
