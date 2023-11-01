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
#include "t-content.c"
#include "t-dot.c"
#include "t-get_d_vec_2exp.c"
#include "t-get_set_fft.c"
#include "t-get_set_nmod_vec.c"
#include "t-height.c"
#include "t-height_index.c"
#include "t-init_clear.c"
#include "t-is_zero.c"
#include "t-lcm.c"
#include "t-max_bits.c"
#include "t-max_limbs.c"
#include "t-neg.c"
#include "t-prod.c"
#include "t-scalar_abs.c"
#include "t-scalar_addmul_fmpz.c"
#include "t-scalar_addmul_si_2exp.c"
#include "t-scalar_addmul_si.c"
#include "t-scalar_addmul_ui.c"
#include "t-scalar_divexact_fmpz.c"
#include "t-scalar_divexact_si.c"
#include "t-scalar_divexact_ui.c"
#include "t-scalar_fdiv_q_fmpz.c"
#include "t-scalar_mod_fmpz.c"
#include "t-scalar_mul_2exp.c"
#include "t-scalar_mul_fmpz.c"
#include "t-scalar_mul_si.c"
#include "t-scalar_mul_ui.c"
#include "t-scalar_smod_fmpz.c"
#include "t-scalar_submul_fmpz.c"
#include "t-scalar_submul_si_2exp.c"
#include "t-scalar_submul_si.c"
#include "t-set_equal.c"
#include "t-sub.c"
#include "t-sum.c"
#include "t-sum_max_bits.c"
#include "t-swap.c"
#include "t-zero.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(fmpz_vec_add),
    TEST_FUNCTION(fmpz_vec_content),
    TEST_FUNCTION(fmpz_vec_dot),
    TEST_FUNCTION(fmpz_vec_get_d_vec_2exp),
    TEST_FUNCTION(fmpz_vec_get_set_fft),
    TEST_FUNCTION(fmpz_vec_get_set_nmod_vec),
    TEST_FUNCTION(fmpz_vec_height),
    TEST_FUNCTION(fmpz_vec_height_index),
    TEST_FUNCTION(fmpz_vec_init_clear),
    TEST_FUNCTION(fmpz_vec_is_zero),
    TEST_FUNCTION(fmpz_vec_lcm),
    TEST_FUNCTION(fmpz_vec_max_bits),
    TEST_FUNCTION(fmpz_vec_max_limbs),
    TEST_FUNCTION(fmpz_vec_neg),
    TEST_FUNCTION(fmpz_vec_prod),
    TEST_FUNCTION(fmpz_vec_scalar_abs),
    TEST_FUNCTION(fmpz_vec_scalar_addmul_fmpz),
    TEST_FUNCTION(fmpz_vec_scalar_addmul_si_2exp),
    TEST_FUNCTION(fmpz_vec_scalar_addmul_si),
    TEST_FUNCTION(fmpz_vec_scalar_addmul_ui),
    TEST_FUNCTION(fmpz_vec_scalar_divexact_fmpz),
    TEST_FUNCTION(fmpz_vec_scalar_divexact_si),
    TEST_FUNCTION(fmpz_vec_scalar_divexact_ui),
    TEST_FUNCTION(fmpz_vec_scalar_fdiv_q_fmpz),
    TEST_FUNCTION(fmpz_vec_scalar_mod_fmpz),
    TEST_FUNCTION(fmpz_vec_scalar_mul_2exp),
    TEST_FUNCTION(fmpz_vec_scalar_mul_fmpz),
    TEST_FUNCTION(fmpz_vec_scalar_mul_si),
    TEST_FUNCTION(fmpz_vec_scalar_mul_ui),
    TEST_FUNCTION(fmpz_vec_scalar_smod_fmpz),
    TEST_FUNCTION(fmpz_vec_scalar_submul_fmpz),
    TEST_FUNCTION(fmpz_vec_scalar_submul_si_2exp),
    TEST_FUNCTION(fmpz_vec_scalar_submul_si),
    TEST_FUNCTION(fmpz_vec_set_equal),
    TEST_FUNCTION(fmpz_vec_sub),
    TEST_FUNCTION(fmpz_vec_sum),
    TEST_FUNCTION(fmpz_vec_sum_max_bits),
    TEST_FUNCTION(fmpz_vec_swap),
    TEST_FUNCTION(fmpz_vec_zero)
};

/* main function *************************************************************/

TEST_MAIN(tests)
