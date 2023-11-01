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

#include "t-add_sub_neg.c"
#include "t-can_solve.c"
#include "t-charpoly.c"
#include "t-fmpz_vec_mul.c"
#include "t-get_set_fmpz_mat.c"
#include "t-howell_form.c"
#include "t-init_clear.c"
#include "t-inv.c"
#include "t-lu.c"
#include "t-minpoly.c"
#include "t-mul.c"
#include "t-mul_classical_threaded.c"
#include "t-mul_fmpz_vec.c"
#include "t-nullspace.c"
#include "t-rank.c"
#include "t-rref.c"
#include "t-scalar_mul_fmpz.c"
#include "t-scalar_mul_si.c"
#include "t-scalar_mul_ui.c"
#include "t-solve.c"
#include "t-solve_tril.c"
#include "t-solve_triu.c"
#include "t-sqr.c"
#include "t-trace.c"
#include "t-window_init_clear.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(fmpz_mod_mat_add_sub_neg),
    TEST_FUNCTION(fmpz_mod_mat_can_solve),
    TEST_FUNCTION(fmpz_mod_mat_charpoly),
    TEST_FUNCTION(fmpz_mod_mat_fmpz_vec_mul),
    TEST_FUNCTION(fmpz_mod_mat_get_set_fmpz_mat),
    TEST_FUNCTION(fmpz_mod_mat_howell_form),
    TEST_FUNCTION(fmpz_mod_mat_init_clear),
    TEST_FUNCTION(fmpz_mod_mat_inv),
    TEST_FUNCTION(fmpz_mod_mat_lu),
    TEST_FUNCTION(fmpz_mod_mat_minpoly),
    TEST_FUNCTION(fmpz_mod_mat_mul),
    TEST_FUNCTION(fmpz_mod_mat_mul_classical_threaded),
    TEST_FUNCTION(fmpz_mod_mat_mul_fmpz_vec),
    TEST_FUNCTION(fmpz_mod_mat_nullspace),
    TEST_FUNCTION(fmpz_mod_mat_rank),
    TEST_FUNCTION(fmpz_mod_mat_rref),
    TEST_FUNCTION(fmpz_mod_mat_scalar_mul_fmpz),
    TEST_FUNCTION(fmpz_mod_mat_scalar_mul_si),
    TEST_FUNCTION(fmpz_mod_mat_scalar_mul_ui),
    TEST_FUNCTION(fmpz_mod_mat_solve),
    TEST_FUNCTION(fmpz_mod_mat_solve_tril),
    TEST_FUNCTION(fmpz_mod_mat_solve_triu),
    TEST_FUNCTION(fmpz_mod_mat_sqr),
    TEST_FUNCTION(fmpz_mod_mat_trace),
    TEST_FUNCTION(fmpz_mod_mat_window_init_clear)
};

/* main function *************************************************************/

TEST_MAIN(tests)
