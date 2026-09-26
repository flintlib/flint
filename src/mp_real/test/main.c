/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Include functions *********************************************************/

#include "t-atan_bitwise_rs.c"
#include "t-atan_frac_bsplit.c"
#include "t-const.c"
#include "t-misc.c"
#include "t-prec.c"
#include "t-bitwise_rs_stress.c"
#include "t-const_euler.c"
#include "t-const_log2.c"
#include "t-const_misc.c"
#include "t-const_pi4.c"
#include "t-exp_bitwise_rs.c"
#include "t-exp_rs.c"
#include "t-log1p_bitwise_rs.c"
#include "t-machin_caches.c"
#include "t-sin_cos_bitwise_rs.c"
#include "t-sin_cos_bits.c"
#include "t-exp_log_atan_bits.c"
#include "t-rel_tab.c"
#include "t-sin_cos_diophantine.c"
#include "t-sin_cos_reduced.c"
#include "t-sin_cos_sum_bs.c"
#include "t-div_newton.c"
#include "t-sqrt_newton.c"
#include "t-tab_bsplit.c"
#include "t-arith.c"
#include "t-api.c"
#include "t-submul_bounded.c"
#include "t-mul_complex.c"
#include "t-agm.c"
#include "t-const_series.c"
#include "t-hypgeom_series.c"
#include "t-newton.c"
#include "t-exp_diophantine.c"
#include "t-exp_notab.c"
#include "t-exp_reduced.c"
#include "t-exp_sum_bs.c"
#include "t-tab_floor.c"
#include "t-tab_static.c"
#include "t-tan_bitwise_rs.c"
#include "t-series_rs.c"
#include "t-series_tapered.c"
#include "t-trig_rs.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(mp_real_atan_bitwise_rs),
    TEST_FUNCTION(mp_real_atan_frac_bsplit),
    TEST_FUNCTION(mp_real_const),
    TEST_FUNCTION(mp_real_misc),
    TEST_FUNCTION(mp_real_prec),
    TEST_FUNCTION(mp_real_bitwise_rs_stress),
    TEST_FUNCTION(mp_real_const_euler),
    TEST_FUNCTION(mp_real_const_log2),
    TEST_FUNCTION(mp_real_const_misc),
    TEST_FUNCTION(mp_real_const_pi4),
    TEST_FUNCTION(mp_real_exp_bitwise_rs),
    TEST_FUNCTION(mp_real_exp_rs),
    TEST_FUNCTION(mp_real_log1p_bitwise_rs),
    TEST_FUNCTION(mp_real_machin_caches),
    TEST_FUNCTION(mp_real_sin_cos_bitwise_rs),
    TEST_FUNCTION(mp_real_sin_cos_bits),
    TEST_FUNCTION(mp_real_exp_log_atan_bits),
    TEST_FUNCTION(mp_real_rel_tab),
    TEST_FUNCTION(mp_real_sin_cos_diophantine),
    TEST_FUNCTION(mp_real_sin_cos_reduced),
    TEST_FUNCTION(mp_real_sin_cos_sum_bs),
    TEST_FUNCTION(mp_real_div_newton),
    TEST_FUNCTION(mp_real_sqrt_newton),
    TEST_FUNCTION(mp_real_tab_bsplit),
    TEST_FUNCTION(mp_real_arith),
    TEST_FUNCTION(mp_real_api),
    TEST_FUNCTION(mp_real_submul_bounded),
    TEST_FUNCTION(mp_real_mul_complex),
    TEST_FUNCTION(mp_real_agm),
    TEST_FUNCTION(mp_real_const_series),
    TEST_FUNCTION(mp_real_hypgeom_series),
    TEST_FUNCTION(mp_real_newton),
    TEST_FUNCTION(mp_real_exp_diophantine),
    TEST_FUNCTION(mp_real_exp_notab),
    TEST_FUNCTION(mp_real_exp_reduced),
    TEST_FUNCTION(mp_real_exp_sum_bs),
    TEST_FUNCTION(mp_real_tab_floor),
    TEST_FUNCTION(mp_real_tab_static),
    TEST_FUNCTION(mp_real_tan_bitwise_rs),
    TEST_FUNCTION(mp_real_series_rs),
    TEST_FUNCTION(mp_real_series_tapered),
    TEST_FUNCTION(mp_real_trig_rs)
};

/* main function *************************************************************/

TEST_MAIN(tests)
