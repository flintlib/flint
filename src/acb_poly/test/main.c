/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <stdlib.h>

/* Include functions *********************************************************/

#include "t-add.c"
#include "t-add_series.c"
#include "t-add_si.c"
#include "t-agm1_series.c"
#include "t-atan_series.c"
#include "t-binomial_transform_basecase.c"
#include "t-binomial_transform.c"
#include "t-binomial_transform_convolution.c"
#include "t-borel_transform.c"
#include "t-compose.c"
#include "t-compose_series.c"
#include "t-cos_pi_series.c"
#include "t-cot_pi_series.c"
#include "t-digamma_series.c"
#include "t-divrem.c"
#include "t-div_series.c"
#include "t-elliptic_k_series.c"
#include "t-elliptic_p_series.c"
#include "t-evaluate2.c"
#include "t-evaluate2_horner.c"
#include "t-evaluate2_rectangular.c"
#include "t-evaluate.c"
#include "t-evaluate_horner.c"
#include "t-evaluate_rectangular.c"
#include "t-evaluate_vec_fast.c"
#include "t-evaluate_vec_iter.c"
#include "t-exp_pi_i_series.c"
#include "t-exp_series_basecase.c"
#include "t-exp_series.c"
#include "t-find_roots.c"
#include "t-gamma_series.c"
#include "t-get_coeff_ptr.c"
#include "t-get_set_coeff_acb.c"
#include "t-get_unique_fmpz_poly.c"
#include "t-graeffe_transform.c"
#include "t-interpolate_barycentric.c"
#include "t-interpolate_fast.c"
#include "t-interpolate_newton.c"
#include "t-inv_series.c"
#include "t-lambertw_series.c"
#include "t-lgamma_series.c"
#include "t-log1p_series.c"
#include "t-log_series.c"
#include "t-mullow.c"
#include "t-mullow_classical.c"
#include "t-mullow_transpose.c"
#include "t-mullow_transpose_gauss.c"
#include "t-polylog_cpx.c"
#include "t-polylog_series.c"
#include "t-pow_acb_series.c"
#include "t-pow_series.c"
#include "t-powsum_one_series_sieved.c"
#include "t-powsum_series_naive_threaded.c"
#include "t-pow_ui.c"
#include "t-pow_ui_trunc_binexp.c"
#include "t-product_roots.c"
#include "t-revert_series.c"
#include "t-rgamma_series.c"
#include "t-rising_ui_series.c"
#include "t-root_bound_fujiwara.c"
#include "t-rsqrt_series.c"
#include "t-set_trunc_round.c"
#include "t-shift_left_right.c"
#include "t-sin_cos_pi_series.c"
#include "t-sin_cos_series.c"
#include "t-sinc_series.c"
#include "t-sinh_cosh_series.c"
#include "t-sin_pi_series.c"
#include "t-sin_series_cos_series.c"
#include "t-sqrt_series.c"
#include "t-sub.c"
#include "t-sub_series.c"
#include "t-tan_series.c"
#include "t-taylor_shift.c"
#include "t-zeta_cpx_series.c"
#include "t-zeta_em_tail_bsplit.c"
#include "t-zeta_series.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(acb_poly_add),
    TEST_FUNCTION(acb_poly_add_series),
    TEST_FUNCTION(acb_poly_add_si),
    TEST_FUNCTION(acb_poly_agm1_series),
    TEST_FUNCTION(acb_poly_atan_series),
    TEST_FUNCTION(acb_poly_binomial_transform_basecase),
    TEST_FUNCTION(acb_poly_binomial_transform),
    TEST_FUNCTION(acb_poly_binomial_transform_convolution),
    TEST_FUNCTION(acb_poly_borel_transform),
    TEST_FUNCTION(acb_poly_compose),
    TEST_FUNCTION(acb_poly_compose_series),
    TEST_FUNCTION(acb_poly_cos_pi_series),
    TEST_FUNCTION(acb_poly_cot_pi_series),
    TEST_FUNCTION(acb_poly_digamma_series),
    TEST_FUNCTION(acb_poly_divrem),
    TEST_FUNCTION(acb_poly_div_series),
    TEST_FUNCTION(acb_poly_elliptic_k_series),
    TEST_FUNCTION(acb_poly_elliptic_p_series),
    TEST_FUNCTION(acb_poly_evaluate2),
    TEST_FUNCTION(acb_poly_evaluate2_horner),
    TEST_FUNCTION(acb_poly_evaluate2_rectangular),
    TEST_FUNCTION(acb_poly_evaluate),
    TEST_FUNCTION(acb_poly_evaluate_horner),
    TEST_FUNCTION(acb_poly_evaluate_rectangular),
    TEST_FUNCTION(acb_poly_evaluate_vec_fast),
    TEST_FUNCTION(acb_poly_evaluate_vec_iter),
    TEST_FUNCTION(acb_poly_exp_pi_i_series),
    TEST_FUNCTION(acb_poly_exp_series_basecase),
    TEST_FUNCTION(acb_poly_exp_series),
    TEST_FUNCTION(acb_poly_find_roots),
    TEST_FUNCTION(acb_poly_gamma_series),
    TEST_FUNCTION(acb_poly_get_coeff_ptr),
    TEST_FUNCTION(acb_poly_get_set_coeff_acb),
    TEST_FUNCTION(acb_poly_get_unique_fmpz_poly),
    TEST_FUNCTION(acb_poly_graeffe_transform),
    TEST_FUNCTION(acb_poly_interpolate_barycentric),
    TEST_FUNCTION(acb_poly_interpolate_fast),
    TEST_FUNCTION(acb_poly_interpolate_newton),
    TEST_FUNCTION(acb_poly_inv_series),
    TEST_FUNCTION(acb_poly_lambertw_series),
    TEST_FUNCTION(acb_poly_lgamma_series),
    TEST_FUNCTION(acb_poly_log1p_series),
    TEST_FUNCTION(acb_poly_log_series),
    TEST_FUNCTION(acb_poly_mullow),
    TEST_FUNCTION(acb_poly_mullow_classical),
    TEST_FUNCTION(acb_poly_mullow_transpose),
    TEST_FUNCTION(acb_poly_mullow_transpose_gauss),
    TEST_FUNCTION(acb_poly_polylog_cpx),
    TEST_FUNCTION(acb_poly_polylog_series),
    TEST_FUNCTION(acb_poly_pow_acb_series),
    TEST_FUNCTION(acb_poly_pow_series),
    TEST_FUNCTION(acb_poly_powsum_one_series_sieved),
    TEST_FUNCTION(acb_poly_powsum_series_naive_threaded),
    TEST_FUNCTION(acb_poly_pow_ui),
    TEST_FUNCTION(acb_poly_pow_ui_trunc_binexp),
    TEST_FUNCTION(acb_poly_product_roots),
    TEST_FUNCTION(acb_poly_revert_series),
    TEST_FUNCTION(acb_poly_rgamma_series),
    TEST_FUNCTION(acb_poly_rising_ui_series),
    TEST_FUNCTION(acb_poly_root_bound_fujiwara),
    TEST_FUNCTION(acb_poly_rsqrt_series),
    TEST_FUNCTION(acb_poly_set_trunc_round),
    TEST_FUNCTION(acb_poly_shift_left_right),
    TEST_FUNCTION(acb_poly_sin_cos_pi_series),
    TEST_FUNCTION(acb_poly_sin_cos_series),
    TEST_FUNCTION(acb_poly_sinc_series),
    TEST_FUNCTION(acb_poly_sinh_cosh_series),
    TEST_FUNCTION(acb_poly_sin_pi_series),
    TEST_FUNCTION(acb_poly_sin_series_cos_series),
    TEST_FUNCTION(acb_poly_sqrt_series),
    TEST_FUNCTION(acb_poly_sub),
    TEST_FUNCTION(acb_poly_sub_series),
    TEST_FUNCTION(acb_poly_tan_series),
    TEST_FUNCTION(acb_poly_taylor_shift),
    TEST_FUNCTION(acb_poly_zeta_cpx_series),
    TEST_FUNCTION(acb_poly_zeta_em_tail_bsplit),
    TEST_FUNCTION(acb_poly_zeta_series)
};

/* main function *************************************************************/

TEST_MAIN(tests)
