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
#include "t-div.c"
#include "t-equal_fmpz_fmpq.c"
#include "t-get_fmpz_mod_poly.c"
#include "t-get_nmod_poly.c"
#include "t-get_set_den.c"
#include "t-get_set_fmpq_poly.c"
#include "t-get_set_fmpz_mat_row.c"
#include "t-init_clear.c"
#include "t-inv.c"
#include "t-is_rational_integer.c"
#include "t-mod_fmpz.c"
#include "t-mul.c"
#include "t-mul_div_fmpq.c"
#include "t-mul_gen.c"
#include "t-norm.c"
#include "t-norm_div.c"
#include "t-pow.c"
#include "t-rep_mat.c"
#include "t-rep_mat_fmpz_mat_den.c"
#include "t-set_coeff_num_fmpz.c"
#include "t-set_equal.c"
#include "t-set_equal_si_ui.c"
#include "t-trace.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(nf_elem_add_sub),
    TEST_FUNCTION(nf_elem_div),
    TEST_FUNCTION(nf_elem_equal_fmpz_fmpq),
    TEST_FUNCTION(nf_elem_get_fmpz_mod_poly),
    TEST_FUNCTION(nf_elem_get_nmod_poly),
    TEST_FUNCTION(nf_elem_get_set_den),
    TEST_FUNCTION(nf_elem_get_set_fmpq_poly),
    TEST_FUNCTION(nf_elem_get_set_fmpz_mat_row),
    TEST_FUNCTION(nf_elem_init_clear),
    TEST_FUNCTION(nf_elem_inv),
    TEST_FUNCTION(nf_elem_is_rational_integer),
    TEST_FUNCTION(nf_elem_mod_fmpz),
    TEST_FUNCTION(nf_elem_mul),
    TEST_FUNCTION(nf_elem_mul_div_fmpq),
    TEST_FUNCTION(nf_elem_mul_gen),
    TEST_FUNCTION(nf_elem_norm),
    TEST_FUNCTION(nf_elem_norm_div),
    TEST_FUNCTION(nf_elem_pow),
    TEST_FUNCTION(nf_elem_rep_mat),
    TEST_FUNCTION(nf_elem_rep_mat_fmpz_mat_den),
    TEST_FUNCTION(nf_elem_set_coeff_num_fmpz),
    TEST_FUNCTION(nf_elem_set_equal),
    TEST_FUNCTION(nf_elem_set_equal_si_ui),
    TEST_FUNCTION(nf_elem_trace)
};

/* main function *************************************************************/

TEST_MAIN(tests)
