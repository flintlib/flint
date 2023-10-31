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

#include "t-add_2exp_fmpz.c"
#include "t-add.c"
#include "t-addmul.c"
#include "t-atan.c"
#include "t-binpow_uiui.c"
#include "t-bin_uiui.c"
#include "t-cmp_2exp_si.c"
#include "t-cmp.c"
#include "t-cosh.c"
#include "t-div.c"
#include "t-div_lower.c"
#include "t-d_log_lower_bound.c"
#include "t-d_log_upper_bound.c"
#include "t-dump_file.c"
#include "t-dump_str.c"
#include "t-exp.c"
#include "t-expinv.c"
#include "t-expm1.c"
#include "t-exp_tail.c"
#include "t-fac_ui.c"
#include "t-fast_add_2exp_si.c"
#include "t-fast_addmul.c"
#include "t-fast_mul_2exp_si.c"
#include "t-fast_mul.c"
#include "t-geom_series.c"
#include "t-get_d.c"
#include "t-hurwitz_zeta_uiui.c"
#include "t-log1p.c"
#include "t-log.c"
#include "t-mul_2exp_fmpz.c"
#include "t-mul_2exp_si.c"
#include "t-mul.c"
#include "t-mul_lower.c"
#include "t-neg_log.c"
#include "t-polylog_tail.c"
#include "t-pow_fmpz.c"
#include "t-pow_ui.c"
#include "t-rfac_ui.c"
#include "t-root.c"
#include "t-rsqrt.c"
#include "t-rsqrt_lower.c"
#include "t-set_d_2exp_fmpz.c"
#include "t-set_d.c"
#include "t-set_ui.c"
#include "t-set_ui_lower.c"
#include "t-sinh.c"
#include "t-sqrt.c"
#include "t-sqrt_lower.c"
#include "t-sub.c"
#include "t-sub_lower.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(mag_add_2exp_fmpz),
    TEST_FUNCTION(mag_add),
    TEST_FUNCTION(mag_addmul),
    TEST_FUNCTION(mag_atan),
    TEST_FUNCTION(mag_binpow_uiui),
    TEST_FUNCTION(mag_bin_uiui),
    TEST_FUNCTION(mag_cmp_2exp_si),
    TEST_FUNCTION(mag_cmp),
    TEST_FUNCTION(mag_cosh),
    TEST_FUNCTION(mag_div),
    TEST_FUNCTION(mag_div_lower),
    TEST_FUNCTION(mag_d_log_lower_bound),
    TEST_FUNCTION(mag_d_log_upper_bound),
    TEST_FUNCTION(mag_dump_file),
    TEST_FUNCTION(mag_dump_str),
    TEST_FUNCTION(mag_exp),
    TEST_FUNCTION(mag_expinv),
    TEST_FUNCTION(mag_expm1),
    TEST_FUNCTION(mag_exp_tail),
    TEST_FUNCTION(mag_fac_ui),
    TEST_FUNCTION(mag_fast_add_2exp_si),
    TEST_FUNCTION(mag_fast_addmul),
    TEST_FUNCTION(mag_fast_mul_2exp_si),
    TEST_FUNCTION(mag_fast_mul),
    TEST_FUNCTION(mag_geom_series),
    TEST_FUNCTION(mag_get_d),
    TEST_FUNCTION(mag_hurwitz_zeta_uiui),
    TEST_FUNCTION(mag_log1p),
    TEST_FUNCTION(mag_log),
    TEST_FUNCTION(mag_mul_2exp_fmpz),
    TEST_FUNCTION(mag_mul_2exp_si),
    TEST_FUNCTION(mag_mul),
    TEST_FUNCTION(mag_mul_lower),
    TEST_FUNCTION(mag_neg_log),
    TEST_FUNCTION(mag_polylog_tail),
    TEST_FUNCTION(mag_pow_fmpz),
    TEST_FUNCTION(mag_pow_ui),
    TEST_FUNCTION(mag_rfac_ui),
    TEST_FUNCTION(mag_root),
    TEST_FUNCTION(mag_rsqrt),
    TEST_FUNCTION(mag_rsqrt_lower),
    TEST_FUNCTION(mag_set_d_2exp_fmpz),
    TEST_FUNCTION(mag_set_d),
    TEST_FUNCTION(mag_set_ui),
    TEST_FUNCTION(mag_set_ui_lower),
    TEST_FUNCTION(mag_sinh),
    TEST_FUNCTION(mag_sqrt),
    TEST_FUNCTION(mag_sqrt_lower),
    TEST_FUNCTION(mag_sub),
    TEST_FUNCTION(mag_sub_lower)
};

/* main function *************************************************************/

TEST_MAIN(tests)
