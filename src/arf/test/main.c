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
#include <mpfr.h>

/* Include functions *********************************************************/

#include "t-abs_bound_le_2exp_fmpz.c"
#include "t-abs_bound_lt_2exp_fmpz.c"
#include "t-abs_bound_lt_2exp_si.c"
#include "t-add.c"
#include "t-add_fmpz_2exp.c"
#include "t-add_fmpz.c"
#include "t-addmul.c"
#include "t-addmul_fmpz.c"
#include "t-addmul_si.c"
#include "t-addmul_ui.c"
#include "t-add_si.c"
#include "t-add_ui.c"
#include "t-approx_dot.c"
#include "t-ceil.c"
#include "t-cmp_2exp_si.c"
#include "t-cmpabs_2exp_si.c"
#include "t-cmpabs.c"
#include "t-cmp.c"
#include "t-complex_mul.c"
#include "t-complex_sqr.c"
#include "t-div.c"
#include "t-dump_file.c"
#include "t-dump_str.c"
#include "t-floor.c"
#include "t-fma.c"
#include "t-frexp.c"
#include "t-get_d.c"
#include "t-get_fmpz.c"
#include "t-get_mpfr.c"
#include "t-get_str.c"
#include "t-is_int_2exp_si.c"
#include "t-mul.c"
#include "t-mul_fmpz.c"
#include "t-mul_si.c"
#include "t-mul_ui.c"
#include "t-mul_via_mpfr.c"
#include "t-neg_round.c"
#include "t-root.c"
#include "t-rsqrt.c"
#include "t-set_d.c"
#include "t-set_fmpq.c"
#include "t-set_fmpz_2exp.c"
#include "t-set_round.c"
#include "t-set_round_fmpz.c"
#include "t-set_round_mpz.c"
#include "t-set_round_ui.c"
#include "t-set_round_uiui.c"
#include "t-sgn.c"
#include "t-sosq.c"
#include "t-sqrt.c"
#include "t-sub.c"
#include "t-sub_fmpz.c"
#include "t-submul.c"
#include "t-submul_fmpz.c"
#include "t-submul_si.c"
#include "t-submul_ui.c"
#include "t-sub_si.c"
#include "t-sub_ui.c"
#include "t-sum.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(arf_abs_bound_le_2exp_fmpz),
    TEST_FUNCTION(arf_abs_bound_lt_2exp_fmpz),
    TEST_FUNCTION(arf_abs_bound_lt_2exp_si),
    TEST_FUNCTION(arf_add),
    TEST_FUNCTION(arf_add_fmpz_2exp),
    TEST_FUNCTION(arf_add_fmpz),
    TEST_FUNCTION(arf_addmul),
    TEST_FUNCTION(arf_addmul_fmpz),
    TEST_FUNCTION(arf_addmul_si),
    TEST_FUNCTION(arf_addmul_ui),
    TEST_FUNCTION(arf_add_si),
    TEST_FUNCTION(arf_add_ui),
    TEST_FUNCTION(arf_approx_dot),
    TEST_FUNCTION(arf_ceil),
    TEST_FUNCTION(arf_cmp_2exp_si),
    TEST_FUNCTION(arf_cmpabs_2exp_si),
    TEST_FUNCTION(arf_cmpabs),
    TEST_FUNCTION(arf_cmp),
    TEST_FUNCTION(arf_complex_mul),
    TEST_FUNCTION(arf_complex_sqr),
    TEST_FUNCTION(arf_div),
    TEST_FUNCTION(arf_dump_file),
    TEST_FUNCTION(arf_dump_str),
    TEST_FUNCTION(arf_floor),
    TEST_FUNCTION(arf_fma),
    TEST_FUNCTION(arf_frexp),
    TEST_FUNCTION(arf_get_d),
    TEST_FUNCTION(arf_get_fmpz),
    TEST_FUNCTION(arf_get_mpfr),
    TEST_FUNCTION(arf_get_str),
    TEST_FUNCTION(arf_is_int_2exp_si),
    TEST_FUNCTION(arf_mul),
    TEST_FUNCTION(arf_mul_fmpz),
    TEST_FUNCTION(arf_mul_si),
    TEST_FUNCTION(arf_mul_ui),
    TEST_FUNCTION(arf_mul_via_mpfr),
    TEST_FUNCTION(arf_neg_round),
    TEST_FUNCTION(arf_root),
    TEST_FUNCTION(arf_rsqrt),
    TEST_FUNCTION(arf_set_d),
    TEST_FUNCTION(arf_set_fmpq),
    TEST_FUNCTION(arf_set_fmpz_2exp),
    TEST_FUNCTION(arf_set_round),
    TEST_FUNCTION(arf_set_round_fmpz),
    TEST_FUNCTION(arf_set_round_mpz),
    TEST_FUNCTION(arf_set_round_ui),
    TEST_FUNCTION(arf_set_round_uiui),
    TEST_FUNCTION(arf_sgn),
    TEST_FUNCTION(arf_sosq),
    TEST_FUNCTION(arf_sqrt),
    TEST_FUNCTION(arf_sub),
    TEST_FUNCTION(arf_sub_fmpz),
    TEST_FUNCTION(arf_submul),
    TEST_FUNCTION(arf_submul_fmpz),
    TEST_FUNCTION(arf_submul_si),
    TEST_FUNCTION(arf_submul_ui),
    TEST_FUNCTION(arf_sub_si),
    TEST_FUNCTION(arf_sub_ui),
    TEST_FUNCTION(arf_sum)
};

/* main function *************************************************************/

TEST_MAIN(tests)
