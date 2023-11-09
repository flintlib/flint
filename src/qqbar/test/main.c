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
#include "fexpr.h"

/* Include functions *********************************************************/

#include "t-abs2.c"
#include "t-abs.c"
#include "t-acos_pi.c"
#include "t-acot_pi.c"
#include "t-acsc_pi.c"
#include "t-add.c"
#include "t-asec_pi.c"
#include "t-asin_pi.c"
#include "t-atan_pi.c"
#include "t-ceil.c"
#include "t-cmpabs.c"
#include "t-cmpabs_im.c"
#include "t-cmpabs_re.c"
#include "t-cmp_im.c"
#include "t-cmp_re.c"
#include "t-conjugates.c"
#include "t-cos_pi.c"
#include "t-cot_pi.c"
#include "t-csc_pi.c"
#include "t-csgn.c"
#include "t-div.c"
#include "t-equal_fmpq_poly_val.c"
#include "t-evaluate_fmpq_poly.c"
#include "t-evaluate_fmpz_mpoly.c"
#include "t-exp_pi_i.c"
#include "t-express_in_field.c"
#include "t-floor.c"
#include "t-fmpz_poly_composed_op.c"
#include "t-get_acb.c"
#include "t-get_fexpr.c"
#include "t-get_fexpr_formula.c"
#include "t-get_quadratic.c"
#include "t-guess.c"
#include "t-inv.c"
#include "t-log_pi_i.c"
#include "t-mul_2exp_si.c"
#include "t-mul.c"
#include "t-pow.c"
#include "t-pow_fmpq.c"
#include "t-pow_fmpz.c"
#include "t-pow_si.c"
#include "t-pow_ui.c"
#include "t-randtest.c"
#include "t-re_im.c"
#include "t-root_of_unity.c"
#include "t-roots_fmpz_poly.c"
#include "t-root_ui.c"
#include "t-sec_pi.c"
#include "t-set_d.c"
#include "t-set_re_im_d.c"
#include "t-sgn.c"
#include "t-sgn_re.c"
#include "t-sin_pi.c"
#include "t-sub.c"
#include "t-tan_pi.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(qqbar_abs2),
    TEST_FUNCTION(qqbar_abs),
    TEST_FUNCTION(qqbar_acos_pi),
    TEST_FUNCTION(qqbar_acot_pi),
    TEST_FUNCTION(qqbar_acsc_pi),
    TEST_FUNCTION(qqbar_add),
    TEST_FUNCTION(qqbar_asec_pi),
    TEST_FUNCTION(qqbar_asin_pi),
    TEST_FUNCTION(qqbar_atan_pi),
    TEST_FUNCTION(qqbar_ceil),
    TEST_FUNCTION(qqbar_cmpabs),
    TEST_FUNCTION(qqbar_cmpabs_im),
    TEST_FUNCTION(qqbar_cmpabs_re),
    TEST_FUNCTION(qqbar_cmp_im),
    TEST_FUNCTION(qqbar_cmp_re),
    TEST_FUNCTION(qqbar_conjugates),
    TEST_FUNCTION(qqbar_cos_pi),
    TEST_FUNCTION(qqbar_cot_pi),
    TEST_FUNCTION(qqbar_csc_pi),
    TEST_FUNCTION(qqbar_csgn),
    TEST_FUNCTION(qqbar_div),
    TEST_FUNCTION(qqbar_equal_fmpq_poly_val),
    TEST_FUNCTION(qqbar_evaluate_fmpq_poly),
    TEST_FUNCTION(qqbar_evaluate_fmpz_mpoly),
    TEST_FUNCTION(qqbar_exp_pi_i),
    TEST_FUNCTION(qqbar_express_in_field),
    TEST_FUNCTION(qqbar_floor),
    TEST_FUNCTION(qqbar_fmpz_poly_composed_op),
    TEST_FUNCTION(qqbar_get_acb),
    TEST_FUNCTION(qqbar_get_fexpr),
    TEST_FUNCTION(qqbar_get_fexpr_formula),
    TEST_FUNCTION(qqbar_get_quadratic),
    TEST_FUNCTION(qqbar_guess),
    TEST_FUNCTION(qqbar_inv),
    TEST_FUNCTION(qqbar_log_pi_i),
    TEST_FUNCTION(qqbar_mul_2exp_si),
    TEST_FUNCTION(qqbar_mul),
    TEST_FUNCTION(qqbar_pow),
    TEST_FUNCTION(qqbar_pow_fmpq),
    TEST_FUNCTION(qqbar_pow_fmpz),
    TEST_FUNCTION(qqbar_pow_si),
    TEST_FUNCTION(qqbar_pow_ui),
    TEST_FUNCTION(qqbar_randtest),
    TEST_FUNCTION(qqbar_re_im),
    TEST_FUNCTION(qqbar_root_of_unity),
    TEST_FUNCTION(qqbar_roots_fmpz_poly),
    TEST_FUNCTION(qqbar_root_ui),
    TEST_FUNCTION(qqbar_sec_pi),
    TEST_FUNCTION(qqbar_set_d),
    TEST_FUNCTION(qqbar_set_re_im_d),
    TEST_FUNCTION(qqbar_sgn),
    TEST_FUNCTION(qqbar_sgn_re),
    TEST_FUNCTION(qqbar_sin_pi),
    TEST_FUNCTION(qqbar_sub),
    TEST_FUNCTION(qqbar_tan_pi)
};

/* main function *************************************************************/

TEST_MAIN(tests)
