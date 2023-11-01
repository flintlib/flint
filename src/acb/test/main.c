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

#include "t-acos.c"
#include "t-acosh.c"
#include "t-agm1.c"
#include "t-agm.c"
#include "t-approx_dot.c"
#include "t-asin.c"
#include "t-asinh.c"
#include "t-atan.c"
#include "t-atanh.c"
#include "t-barnes_g.c"
#include "t-bernoulli_poly_ui.c"
#include "t-chebyshev_t_ui.c"
#include "t-chebyshev_u_ui.c"
#include "t-cos_pi.c"
#include "t-cot.c"
#include "t-coth.c"
#include "t-cot_pi.c"
#include "t-csc.c"
#include "t-csch.c"
#include "t-csc_pi.c"
#include "t-csgn.c"
#include "t-digamma.c"
#include "t-div.c"
#include "t-dot.c"
#include "t-dot_fmpz.c"
#include "t-dot_si.c"
#include "t-dot_siui.c"
#include "t-dot_ui.c"
#include "t-dot_uiui.c"
#include "t-exp.c"
#include "t-exp_invexp.c"
#include "t-expm1.c"
#include "t-exp_pi_i.c"
#include "t-gamma.c"
#include "t-get_abs_lbound_arf.c"
#include "t-get_abs_ubound_arf.c"
#include "t-get_mag.c"
#include "t-get_mag_lower.c"
#include "t-inv.c"
#include "t-lambertw.c"
#include "t-lgamma.c"
#include "t-log1p.c"
#include "t-log.c"
#include "t-log_sin_pi.c"
#include "t-mul.c"
#include "t-mul_naive.c"
#include "t-polygamma.c"
#include "t-pow.c"
#include "t-pow_fmpz.c"
#include "t-quadratic_roots_fmpz.c"
#include "t-rel_accuracy_bits.c"
#include "t-rgamma.c"
#include "t-rising2_ui.c"
#include "t-rising_ui.c"
#include "t-rising_ui_get_mag.c"
#include "t-root_ui.c"
#include "t-rsqrt.c"
#include "t-sec.c"
#include "t-sech.c"
#include "t-sgn.c"
#include "t-sinc.c"
#include "t-sin_cos.c"
#include "t-sinc_pi.c"
#include "t-sinh_cosh.c"
#include "t-sin_pi.c"
#include "t-sqrt.c"
#include "t-sqrts.c"
#include "t-tan.c"
#include "t-tanh.c"
#include "t-tan_pi.c"
#include "t-urandom.c"
#include "t-vec_set_real_imag.c"
#include "t-vec_unit_roots.c"
#include "t-zeta.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(acb_acos),
    TEST_FUNCTION(acb_acosh),
    TEST_FUNCTION(acb_agm1),
    TEST_FUNCTION(acb_agm),
    TEST_FUNCTION(acb_approx_dot),
    TEST_FUNCTION(acb_asin),
    TEST_FUNCTION(acb_asinh),
    TEST_FUNCTION(acb_atan),
    TEST_FUNCTION(acb_atanh),
    TEST_FUNCTION(acb_barnes_g),
    TEST_FUNCTION(acb_bernoulli_poly_ui),
    TEST_FUNCTION(acb_chebyshev_t_ui),
    TEST_FUNCTION(acb_chebyshev_u_ui),
    TEST_FUNCTION(acb_cos_pi),
    TEST_FUNCTION(acb_cot),
    TEST_FUNCTION(acb_coth),
    TEST_FUNCTION(acb_cot_pi),
    TEST_FUNCTION(acb_csc),
    TEST_FUNCTION(acb_csch),
    TEST_FUNCTION(acb_csc_pi),
    TEST_FUNCTION(acb_csgn),
    TEST_FUNCTION(acb_digamma),
    TEST_FUNCTION(acb_div),
    TEST_FUNCTION(acb_dot),
    TEST_FUNCTION(acb_dot_fmpz),
    TEST_FUNCTION(acb_dot_si),
    TEST_FUNCTION(acb_dot_siui),
    TEST_FUNCTION(acb_dot_ui),
    TEST_FUNCTION(acb_dot_uiui),
    TEST_FUNCTION(acb_exp),
    TEST_FUNCTION(acb_exp_invexp),
    TEST_FUNCTION(acb_expm1),
    TEST_FUNCTION(acb_exp_pi_i),
    TEST_FUNCTION(acb_gamma),
    TEST_FUNCTION(acb_get_abs_lbound_arf),
    TEST_FUNCTION(acb_get_abs_ubound_arf),
    TEST_FUNCTION(acb_get_mag),
    TEST_FUNCTION(acb_get_mag_lower),
    TEST_FUNCTION(acb_inv),
    TEST_FUNCTION(acb_lambertw),
    TEST_FUNCTION(acb_lgamma),
    TEST_FUNCTION(acb_log1p),
    TEST_FUNCTION(acb_log),
    TEST_FUNCTION(acb_log_sin_pi),
    TEST_FUNCTION(acb_mul),
    TEST_FUNCTION(acb_mul_naive),
    TEST_FUNCTION(acb_polygamma),
    TEST_FUNCTION(acb_pow),
    TEST_FUNCTION(acb_pow_fmpz),
    TEST_FUNCTION(acb_quadratic_roots_fmpz),
    TEST_FUNCTION(acb_rel_accuracy_bits),
    TEST_FUNCTION(acb_rgamma),
    TEST_FUNCTION(acb_rising2_ui),
    TEST_FUNCTION(acb_rising_ui),
    TEST_FUNCTION(acb_rising_ui_get_mag),
    TEST_FUNCTION(acb_root_ui),
    TEST_FUNCTION(acb_rsqrt),
    TEST_FUNCTION(acb_sec),
    TEST_FUNCTION(acb_sech),
    TEST_FUNCTION(acb_sgn),
    TEST_FUNCTION(acb_sinc),
    TEST_FUNCTION(acb_sin_cos),
    TEST_FUNCTION(acb_sinc_pi),
    TEST_FUNCTION(acb_sinh_cosh),
    TEST_FUNCTION(acb_sin_pi),
    TEST_FUNCTION(acb_sqrt),
    TEST_FUNCTION(acb_sqrts),
    TEST_FUNCTION(acb_tan),
    TEST_FUNCTION(acb_tanh),
    TEST_FUNCTION(acb_tan_pi),
    TEST_FUNCTION(acb_urandom),
    TEST_FUNCTION(acb_vec_set_real_imag),
    TEST_FUNCTION(acb_vec_unit_roots),
    TEST_FUNCTION(acb_zeta)
};

/* main function *************************************************************/

TEST_MAIN(tests)
