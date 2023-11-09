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
#include "t-ca_poly_evaluate.c"
#include "t-charpoly.c"
#include "t-charpoly_danilevsky.c"
#include "t-companion.c"
#include "t-det.c"
#include "t-dft.c"
#include "t-diagonalization.c"
#include "t-exp.c"
#include "t-inv.c"
#include "t-jordan_blocks.c"
#include "t-jordan_form.c"
#include "t-lu.c"
#include "t-lu_classical.c"
#include "t-lu_recursive.c"
#include "t-mul.c"
#include "t-mul_same_nf.c"
#include "t-nonsingular_solve_adjugate.c"
#include "t-nonsingular_solve.c"
#include "t-nonsingular_solve_fflu.c"
#include "t-nonsingular_solve_lu.c"
#include "t-rank.c"
#include "t-right_kernel.c"
#include "t-rref.c"
#include "t-rref_fflu.c"
#include "t-rref_lu.c"
#include "t-solve_tril.c"
#include "t-solve_triu.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(ca_mat_adjugate),
    TEST_FUNCTION(ca_mat_ca_poly_evaluate),
    TEST_FUNCTION(ca_mat_charpoly),
    TEST_FUNCTION(ca_mat_charpoly_danilevsky),
    TEST_FUNCTION(ca_mat_companion),
    TEST_FUNCTION(ca_mat_det),
    TEST_FUNCTION(ca_mat_dft),
    TEST_FUNCTION(ca_mat_diagonalization),
    TEST_FUNCTION(ca_mat_exp),
    TEST_FUNCTION(ca_mat_inv),
    TEST_FUNCTION(ca_mat_jordan_blocks),
    TEST_FUNCTION(ca_mat_jordan_form),
    TEST_FUNCTION(ca_mat_lu),
    TEST_FUNCTION(ca_mat_lu_classical),
    TEST_FUNCTION(ca_mat_lu_recursive),
    TEST_FUNCTION(ca_mat_mul),
    TEST_FUNCTION(ca_mat_mul_same_nf),
    TEST_FUNCTION(ca_mat_nonsingular_solve_adjugate),
    TEST_FUNCTION(ca_mat_nonsingular_solve),
    TEST_FUNCTION(ca_mat_nonsingular_solve_fflu),
    TEST_FUNCTION(ca_mat_nonsingular_solve_lu),
    TEST_FUNCTION(ca_mat_rank),
    TEST_FUNCTION(ca_mat_right_kernel),
    TEST_FUNCTION(ca_mat_rref),
    TEST_FUNCTION(ca_mat_rref_fflu),
    TEST_FUNCTION(ca_mat_rref_lu),
    TEST_FUNCTION(ca_mat_solve_tril),
    TEST_FUNCTION(ca_mat_solve_triu)
};

/* main function *************************************************************/

TEST_MAIN(tests)
