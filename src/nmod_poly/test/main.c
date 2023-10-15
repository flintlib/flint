/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>

/* Include functions *********************************************************/

#include "t-add.c"
#include "t-add_series.c"
#include "t-asinh_series.c"
#include "t-asin_series.c"
#include "t-atanh_series.c"
#include "t-atan_series.c"
#include "t-berlekamp_massey.c"
#include "t-bit_pack.c"
#include "t-compose.c"
#include "t-compose_horner.c"
#include "t-compose_mod_brent_kung.c"
#include "t-compose_mod_brent_kung_precomp_preinv.c"
#include "t-compose_mod_brent_kung_precomp_preinv_threaded.c"
#include "t-compose_mod_brent_kung_preinv.c"
#include "t-compose_mod_brent_kung_vec_preinv.c"
#include "t-compose_mod_brent_kung_vec_preinv_threaded.c"
#include "t-compose_mod.c"
#include "t-compose_mod_horner.c"
#include "t-compose_series.c"
#include "t-cosh_series.c"
#include "t-cos_series.c"
#include "t-deflate.c"
#include "t-derivative.c"
#include "t-discriminant.c"
#include "t-div.c"
#include "t-divides.c"
#include "t-divides_classical.c"
#include "t-div_newton_n_preinv.c"
#include "t-divrem_basecase.c"
#include "t-divrem.c"
#include "t-divrem_newton_n_preinv.c"
#include "t-div_root.c"
#include "t-div_series_basecase.c"
#include "t-div_series.c"
#include "t-equal_trunc.c"
#include "t-evaluate_mat_horner.c"
#include "t-evaluate_mat_paterson_stockmeyer.c"
#include "t-evaluate_nmod.c"
#include "t-evaluate_nmod_vec_fast.c"
#include "t-exp_series.c"
#include "t-find_distinct_nonzero_roots.c"
#include "t-fread_print.c"
#include "t-gcd.c"
#include "t-gcd_euclidean.c"
#include "t-gcd_hgcd.c"
#include "t-gcdinv.c"
#include "t-get_set_coeff_ui.c"
#include "t-get_set_str.c"
#include "t-hgcd.c"
#include "t-inflate.c"
#include "t-init_realloc_clear.c"
#include "t-integral.c"
#include "t-interpolate_nmod_vec_barycentric.c"
#include "t-interpolate_nmod_vec.c"
#include "t-interpolate_nmod_vec_fast.c"
#include "t-interpolate_nmod_vec_newton.c"
#include "t-invmod.c"
#include "t-inv_series_basecase.c"
#include "t-inv_series_newton.c"
#include "t-invsqrt_series.c"
#include "t-log_series.c"
#include "t-make_monic.c"
#include "t-mul.c"
#include "t-mul_classical.c"
#include "t-mulhigh.c"
#include "t-mulhigh_classical.c"
#include "t-mul_KS2.c"
#include "t-mul_KS4.c"
#include "t-mul_KS.c"
#include "t-mullow.c"
#include "t-mullow_classical.c"
#include "t-mullow_KS.c"
#include "t-mulmod.c"
#include "t-mulmod_preinv.c"
#include "t-multi_crt.c"
#include "t-neg.c"
#include "t-pow_binexp.c"
#include "t-pow.c"
#include "t-powers_mod_bsgs.c"
#include "t-powers_mod_naive.c"
#include "t-power_sums.c"
#include "t-power_sums_naive.c"
#include "t-power_sums_schoenhage.c"
#include "t-powmod_fmpz_binexp.c"
#include "t-powmod_fmpz_binexp_preinv.c"
#include "t-powmod_ui_binexp.c"
#include "t-powmod_ui_binexp_preinv.c"
#include "t-powmod_x_fmpz_preinv.c"
#include "t-powmod_x_ui_preinv.c"
#include "t-pow_trunc_binexp.c"
#include "t-pow_trunc.c"
#include "t-product_roots_nmod_vec.c"
#include "t-rem.c"
#include "t-resultant.c"
#include "t-resultant_euclidean.c"
#include "t-resultant_hgcd.c"
#include "t-reverse.c"
#include "t-revert_series.c"
#include "t-revert_series_lagrange.c"
#include "t-revert_series_lagrange_fast.c"
#include "t-revert_series_newton.c"
#include "t-scalar_addmul_nmod.c"
#include "t-scalar_mul_nmod.c"
#include "t-set_trunc.c"
#include "t-shift_left_right.c"
#include "t-sinh_series.c"
#include "t-sin_series.c"
#include "t-sqrt.c"
#include "t-sqrt_series.c"
#include "t-sub.c"
#include "t-sub_series.c"
#include "t-tanh_series.c"
#include "t-tan_series.c"
#include "t-taylor_shift.c"
#include "t-taylor_shift_convolution.c"
#include "t-taylor_shift_horner.c"
#include "t-xgcd.c"
#include "t-xgcd_euclidean.c"
#include "t-xgcd_hgcd.c"

/* Array of test functions ***************************************************/

int (*test_functions[])(void) =
{
    TEST_FUNCTION(nmod_poly_add),
    TEST_FUNCTION(nmod_poly_add_series),
    TEST_FUNCTION(nmod_poly_asinh_series),
    TEST_FUNCTION(nmod_poly_asin_series),
    TEST_FUNCTION(nmod_poly_atanh_series),
    TEST_FUNCTION(nmod_poly_atan_series),
    TEST_FUNCTION(nmod_poly_berlekamp_massey),
    TEST_FUNCTION(nmod_poly_bit_pack),
    TEST_FUNCTION(nmod_poly_compose),
    TEST_FUNCTION(nmod_poly_compose_horner),
    TEST_FUNCTION(nmod_poly_compose_mod_brent_kung),
    TEST_FUNCTION(nmod_poly_compose_mod_brent_kung_precomp_preinv),
    TEST_FUNCTION(nmod_poly_compose_mod_brent_kung_precomp_preinv_threaded),
    TEST_FUNCTION(nmod_poly_compose_mod_brent_kung_preinv),
    TEST_FUNCTION(nmod_poly_compose_mod_brent_kung_vec_preinv),
    TEST_FUNCTION(nmod_poly_compose_mod_brent_kung_vec_preinv_threaded),
    TEST_FUNCTION(nmod_poly_compose_mod),
    TEST_FUNCTION(nmod_poly_compose_mod_horner),
    TEST_FUNCTION(nmod_poly_compose_series),
    TEST_FUNCTION(nmod_poly_cosh_series),
    TEST_FUNCTION(nmod_poly_cos_series),
    TEST_FUNCTION(nmod_poly_deflate),
    TEST_FUNCTION(nmod_poly_derivative),
    TEST_FUNCTION(nmod_poly_discriminant),
    TEST_FUNCTION(nmod_poly_div),
    TEST_FUNCTION(nmod_poly_divides),
    TEST_FUNCTION(nmod_poly_divides_classical),
    TEST_FUNCTION(nmod_poly_div_newton_n_preinv),
    TEST_FUNCTION(nmod_poly_divrem_basecase),
    TEST_FUNCTION(nmod_poly_divrem),
    TEST_FUNCTION(nmod_poly_divrem_newton_n_preinv),
    TEST_FUNCTION(nmod_poly_div_root),
    TEST_FUNCTION(nmod_poly_div_series_basecase),
    TEST_FUNCTION(nmod_poly_div_series),
    TEST_FUNCTION(nmod_poly_equal_trunc),
    TEST_FUNCTION(nmod_poly_evaluate_mat_horner),
    TEST_FUNCTION(nmod_poly_evaluate_mat_paterson_stockmeyer),
    TEST_FUNCTION(nmod_poly_evaluate_nmod),
    TEST_FUNCTION(nmod_poly_evaluate_nmod_vec_fast),
    TEST_FUNCTION(nmod_poly_exp_series),
    TEST_FUNCTION(nmod_poly_find_distinct_nonzero_roots),
    TEST_FUNCTION(nmod_poly_fread_print),
    TEST_FUNCTION(nmod_poly_gcd),
    TEST_FUNCTION(nmod_poly_gcd_euclidean),
    TEST_FUNCTION(nmod_poly_gcd_hgcd),
    TEST_FUNCTION(nmod_poly_gcdinv),
    TEST_FUNCTION(nmod_poly_get_set_coeff_ui),
    TEST_FUNCTION(nmod_poly_get_set_str),
    TEST_FUNCTION(nmod_poly_hgcd),
    TEST_FUNCTION(nmod_poly_inflate),
    TEST_FUNCTION(nmod_poly_init_realloc_clear),
    TEST_FUNCTION(nmod_poly_integral),
    TEST_FUNCTION(nmod_poly_interpolate_nmod_vec_barycentric),
    TEST_FUNCTION(nmod_poly_interpolate_nmod_vec),
    TEST_FUNCTION(nmod_poly_interpolate_nmod_vec_fast),
    TEST_FUNCTION(nmod_poly_interpolate_nmod_vec_newton),
    TEST_FUNCTION(nmod_poly_invmod),
    TEST_FUNCTION(nmod_poly_inv_series_basecase),
    TEST_FUNCTION(nmod_poly_inv_series_newton),
    TEST_FUNCTION(nmod_poly_invsqrt_series),
    TEST_FUNCTION(nmod_poly_log_series),
    TEST_FUNCTION(nmod_poly_make_monic),
    TEST_FUNCTION(nmod_poly_mul),
    TEST_FUNCTION(nmod_poly_mul_classical),
    TEST_FUNCTION(nmod_poly_mulhigh),
    TEST_FUNCTION(nmod_poly_mulhigh_classical),
    TEST_FUNCTION(nmod_poly_mul_KS2),
    TEST_FUNCTION(nmod_poly_mul_KS4),
    TEST_FUNCTION(nmod_poly_mul_KS),
    TEST_FUNCTION(nmod_poly_mullow),
    TEST_FUNCTION(nmod_poly_mullow_classical),
    TEST_FUNCTION(nmod_poly_mullow_KS),
    TEST_FUNCTION(nmod_poly_mulmod),
    TEST_FUNCTION(nmod_poly_mulmod_preinv),
    TEST_FUNCTION(nmod_poly_multi_crt),
    TEST_FUNCTION(nmod_poly_neg),
    TEST_FUNCTION(nmod_poly_pow_binexp),
    TEST_FUNCTION(nmod_poly_pow),
    TEST_FUNCTION(nmod_poly_powers_mod_bsgs),
    TEST_FUNCTION(nmod_poly_powers_mod_naive),
    TEST_FUNCTION(nmod_poly_power_sums),
    TEST_FUNCTION(nmod_poly_power_sums_naive),
    TEST_FUNCTION(nmod_poly_power_sums_schoenhage),
    TEST_FUNCTION(nmod_poly_powmod_fmpz_binexp),
    TEST_FUNCTION(nmod_poly_powmod_fmpz_binexp_preinv),
    TEST_FUNCTION(nmod_poly_powmod_ui_binexp),
    TEST_FUNCTION(nmod_poly_powmod_ui_binexp_preinv),
    TEST_FUNCTION(nmod_poly_powmod_x_fmpz_preinv),
    TEST_FUNCTION(nmod_poly_powmod_x_ui_preinv),
    TEST_FUNCTION(nmod_poly_pow_trunc_binexp),
    TEST_FUNCTION(nmod_poly_pow_trunc),
    TEST_FUNCTION(nmod_poly_product_roots_nmod_vec),
    TEST_FUNCTION(nmod_poly_rem),
    TEST_FUNCTION(nmod_poly_resultant),
    TEST_FUNCTION(nmod_poly_resultant_euclidean),
    TEST_FUNCTION(nmod_poly_resultant_hgcd),
    TEST_FUNCTION(nmod_poly_reverse),
    TEST_FUNCTION(nmod_poly_revert_series),
    TEST_FUNCTION(nmod_poly_revert_series_lagrange),
    TEST_FUNCTION(nmod_poly_revert_series_lagrange_fast),
    TEST_FUNCTION(nmod_poly_revert_series_newton),
    TEST_FUNCTION(nmod_poly_scalar_addmul_nmod),
    TEST_FUNCTION(nmod_poly_scalar_mul_nmod),
    TEST_FUNCTION(nmod_poly_set_trunc),
    TEST_FUNCTION(nmod_poly_shift_left_right),
    TEST_FUNCTION(nmod_poly_sinh_series),
    TEST_FUNCTION(nmod_poly_sin_series),
    TEST_FUNCTION(nmod_poly_sqrt),
    TEST_FUNCTION(nmod_poly_sqrt_series),
    TEST_FUNCTION(nmod_poly_sub),
    TEST_FUNCTION(nmod_poly_sub_series),
    TEST_FUNCTION(nmod_poly_tanh_series),
    TEST_FUNCTION(nmod_poly_tan_series),
    TEST_FUNCTION(nmod_poly_taylor_shift),
    TEST_FUNCTION(nmod_poly_taylor_shift_convolution),
    TEST_FUNCTION(nmod_poly_taylor_shift_horner),
    TEST_FUNCTION(nmod_poly_xgcd),
    TEST_FUNCTION(nmod_poly_xgcd_euclidean),
    TEST_FUNCTION(nmod_poly_xgcd_hgcd)
};

char nmod_poly_add_name[] = "nmod_poly_add";
char nmod_poly_add_series_name[] = "nmod_poly_add_series";
char nmod_poly_asinh_series_name[] = "nmod_poly_asinh_series";
char nmod_poly_asin_series_name[] = "nmod_poly_asin_series";
char nmod_poly_atanh_series_name[] = "nmod_poly_atanh_series";
char nmod_poly_atan_series_name[] = "nmod_poly_atan_series";
char nmod_poly_berlekamp_massey_name[] = "nmod_poly_berlekamp_massey";
char nmod_poly_bit_pack_name[] = "nmod_poly_bit_pack";
char nmod_poly_compose_name[] = "nmod_poly_compose";
char nmod_poly_compose_horner_name[] = "nmod_poly_compose_horner";
char nmod_poly_compose_mod_brent_kung_name[] = "nmod_poly_compose_mod_brent_kung";
char nmod_poly_compose_mod_brent_kung_precomp_preinv_name[] = "nmod_poly_compose_mod_brent_kung_precomp_preinv";
char nmod_poly_compose_mod_brent_kung_precomp_preinv_threaded_name[] = "nmod_poly_compose_mod_brent_kung_precomp_preinv_threaded";
char nmod_poly_compose_mod_brent_kung_preinv_name[] = "nmod_poly_compose_mod_brent_kung_preinv";
char nmod_poly_compose_mod_brent_kung_vec_preinv_name[] = "nmod_poly_compose_mod_brent_kung_vec_preinv";
char nmod_poly_compose_mod_brent_kung_vec_preinv_threaded_name[] = "nmod_poly_compose_mod_brent_kung_vec_preinv_threaded";
char nmod_poly_compose_mod_name[] = "nmod_poly_compose_mod";
char nmod_poly_compose_mod_horner_name[] = "nmod_poly_compose_mod_horner";
char nmod_poly_compose_series_name[] = "nmod_poly_compose_series";
char nmod_poly_cosh_series_name[] = "nmod_poly_cosh_series";
char nmod_poly_cos_series_name[] = "nmod_poly_cos_series";
char nmod_poly_deflate_name[] = "nmod_poly_deflate";
char nmod_poly_derivative_name[] = "nmod_poly_derivative";
char nmod_poly_discriminant_name[] = "nmod_poly_discriminant";
char nmod_poly_div_name[] = "nmod_poly_div";
char nmod_poly_divides_name[] = "nmod_poly_divides";
char nmod_poly_divides_classical_name[] = "nmod_poly_divides_classical";
char nmod_poly_div_newton_n_preinv_name[] = "nmod_poly_div_newton_n_preinv";
char nmod_poly_divrem_basecase_name[] = "nmod_poly_divrem_basecase";
char nmod_poly_divrem_name[] = "nmod_poly_divrem";
char nmod_poly_divrem_newton_n_preinv_name[] = "nmod_poly_divrem_newton_n_preinv";
char nmod_poly_div_root_name[] = "nmod_poly_div_root";
char nmod_poly_div_series_basecase_name[] = "nmod_poly_div_series_basecase";
char nmod_poly_div_series_name[] = "nmod_poly_div_series";
char nmod_poly_equal_trunc_name[] = "nmod_poly_equal_trunc";
char nmod_poly_evaluate_mat_horner_name[] = "nmod_poly_evaluate_mat_horner";
char nmod_poly_evaluate_mat_paterson_stockmeyer_name[] = "nmod_poly_evaluate_mat_paterson_stockmeyer";
char nmod_poly_evaluate_nmod_name[] = "nmod_poly_evaluate_nmod";
char nmod_poly_evaluate_nmod_vec_fast_name[] = "nmod_poly_evaluate_nmod_vec_fast";
char nmod_poly_exp_series_name[] = "nmod_poly_exp_series";
char nmod_poly_find_distinct_nonzero_roots_name[] = "nmod_poly_find_distinct_nonzero_roots";
char nmod_poly_fread_print_name[] = "nmod_poly_fread_print";
char nmod_poly_gcd_name[] = "nmod_poly_gcd";
char nmod_poly_gcd_euclidean_name[] = "nmod_poly_gcd_euclidean";
char nmod_poly_gcd_hgcd_name[] = "nmod_poly_gcd_hgcd";
char nmod_poly_gcdinv_name[] = "nmod_poly_gcdinv";
char nmod_poly_get_set_coeff_ui_name[] = "nmod_poly_get_set_coeff_ui";
char nmod_poly_get_set_str_name[] = "nmod_poly_get_set_str";
char nmod_poly_hgcd_name[] = "nmod_poly_hgcd";
char nmod_poly_inflate_name[] = "nmod_poly_inflate";
char nmod_poly_init_realloc_clear_name[] = "nmod_poly_init_realloc_clear";
char nmod_poly_integral_name[] = "nmod_poly_integral";
char nmod_poly_interpolate_nmod_vec_barycentric_name[] = "nmod_poly_interpolate_nmod_vec_barycentric";
char nmod_poly_interpolate_nmod_vec_name[] = "nmod_poly_interpolate_nmod_vec";
char nmod_poly_interpolate_nmod_vec_fast_name[] = "nmod_poly_interpolate_nmod_vec_fast";
char nmod_poly_interpolate_nmod_vec_newton_name[] = "nmod_poly_interpolate_nmod_vec_newton";
char nmod_poly_invmod_name[] = "nmod_poly_invmod";
char nmod_poly_inv_series_basecase_name[] = "nmod_poly_inv_series_basecase";
char nmod_poly_inv_series_newton_name[] = "nmod_poly_inv_series_newton";
char nmod_poly_invsqrt_series_name[] = "nmod_poly_invsqrt_series";
char nmod_poly_log_series_name[] = "nmod_poly_log_series";
char nmod_poly_make_monic_name[] = "nmod_poly_make_monic";
char nmod_poly_mul_name[] = "nmod_poly_mul";
char nmod_poly_mul_classical_name[] = "nmod_poly_mul_classical";
char nmod_poly_mulhigh_name[] = "nmod_poly_mulhigh";
char nmod_poly_mulhigh_classical_name[] = "nmod_poly_mulhigh_classical";
char nmod_poly_mul_KS2_name[] = "nmod_poly_mul_KS2";
char nmod_poly_mul_KS4_name[] = "nmod_poly_mul_KS4";
char nmod_poly_mul_KS_name[] = "nmod_poly_mul_KS";
char nmod_poly_mullow_name[] = "nmod_poly_mullow";
char nmod_poly_mullow_classical_name[] = "nmod_poly_mullow_classical";
char nmod_poly_mullow_KS_name[] = "nmod_poly_mullow_KS";
char nmod_poly_mulmod_name[] = "nmod_poly_mulmod";
char nmod_poly_mulmod_preinv_name[] = "nmod_poly_mulmod_preinv";
char nmod_poly_multi_crt_name[] = "nmod_poly_multi_crt";
char nmod_poly_neg_name[] = "nmod_poly_neg";
char nmod_poly_pow_binexp_name[] = "nmod_poly_pow_binexp";
char nmod_poly_pow_name[] = "nmod_poly_pow";
char nmod_poly_powers_mod_bsgs_name[] = "nmod_poly_powers_mod_bsgs";
char nmod_poly_powers_mod_naive_name[] = "nmod_poly_powers_mod_naive";
char nmod_poly_power_sums_name[] = "nmod_poly_power_sums";
char nmod_poly_power_sums_naive_name[] = "nmod_poly_power_sums_naive";
char nmod_poly_power_sums_schoenhage_name[] = "nmod_poly_power_sums_schoenhage";
char nmod_poly_powmod_fmpz_binexp_name[] = "nmod_poly_powmod_fmpz_binexp";
char nmod_poly_powmod_fmpz_binexp_preinv_name[] = "nmod_poly_powmod_fmpz_binexp_preinv";
char nmod_poly_powmod_ui_binexp_name[] = "nmod_poly_powmod_ui_binexp";
char nmod_poly_powmod_ui_binexp_preinv_name[] = "nmod_poly_powmod_ui_binexp_preinv";
char nmod_poly_powmod_x_fmpz_preinv_name[] = "nmod_poly_powmod_x_fmpz_preinv";
char nmod_poly_powmod_x_ui_preinv_name[] = "nmod_poly_powmod_x_ui_preinv";
char nmod_poly_pow_trunc_binexp_name[] = "nmod_poly_pow_trunc_binexp";
char nmod_poly_pow_trunc_name[] = "nmod_poly_pow_trunc";
char nmod_poly_product_roots_nmod_vec_name[] = "nmod_poly_product_roots_nmod_vec";
char nmod_poly_rem_name[] = "nmod_poly_rem";
char nmod_poly_resultant_name[] = "nmod_poly_resultant";
char nmod_poly_resultant_euclidean_name[] = "nmod_poly_resultant_euclidean";
char nmod_poly_resultant_hgcd_name[] = "nmod_poly_resultant_hgcd";
char nmod_poly_reverse_name[] = "nmod_poly_reverse";
char nmod_poly_revert_series_name[] = "nmod_poly_revert_series";
char nmod_poly_revert_series_lagrange_name[] = "nmod_poly_revert_series_lagrange";
char nmod_poly_revert_series_lagrange_fast_name[] = "nmod_poly_revert_series_lagrange_fast";
char nmod_poly_revert_series_newton_name[] = "nmod_poly_revert_series_newton";
char nmod_poly_scalar_addmul_nmod_name[] = "nmod_poly_scalar_addmul_nmod";
char nmod_poly_scalar_mul_nmod_name[] = "nmod_poly_scalar_mul_nmod";
char nmod_poly_set_trunc_name[] = "nmod_poly_set_trunc";
char nmod_poly_shift_left_right_name[] = "nmod_poly_shift_left_right";
char nmod_poly_sinh_series_name[] = "nmod_poly_sinh_series";
char nmod_poly_sin_series_name[] = "nmod_poly_sin_series";
char nmod_poly_sqrt_name[] = "nmod_poly_sqrt";
char nmod_poly_sqrt_series_name[] = "nmod_poly_sqrt_series";
char nmod_poly_sub_name[] = "nmod_poly_sub";
char nmod_poly_sub_series_name[] = "nmod_poly_sub_series";
char nmod_poly_tanh_series_name[] = "nmod_poly_tanh_series";
char nmod_poly_tan_series_name[] = "nmod_poly_tan_series";
char nmod_poly_taylor_shift_name[] = "nmod_poly_taylor_shift";
char nmod_poly_taylor_shift_convolution_name[] = "nmod_poly_taylor_shift_convolution";
char nmod_poly_taylor_shift_horner_name[] = "nmod_poly_taylor_shift_horner";
char nmod_poly_xgcd_name[] = "nmod_poly_xgcd";
char nmod_poly_xgcd_euclidean_name[] = "nmod_poly_xgcd_euclidean";
char nmod_poly_xgcd_hgcd_name[] = "nmod_poly_xgcd_hgcd";

char * test_names[] =
{
    nmod_poly_add_name,
    nmod_poly_add_series_name,
    nmod_poly_asinh_series_name,
    nmod_poly_asin_series_name,
    nmod_poly_atanh_series_name,
    nmod_poly_atan_series_name,
    nmod_poly_berlekamp_massey_name,
    nmod_poly_bit_pack_name,
    nmod_poly_compose_name,
    nmod_poly_compose_horner_name,
    nmod_poly_compose_mod_brent_kung_name,
    nmod_poly_compose_mod_brent_kung_precomp_preinv_name,
    nmod_poly_compose_mod_brent_kung_precomp_preinv_threaded_name,
    nmod_poly_compose_mod_brent_kung_preinv_name,
    nmod_poly_compose_mod_brent_kung_vec_preinv_name,
    nmod_poly_compose_mod_brent_kung_vec_preinv_threaded_name,
    nmod_poly_compose_mod_name,
    nmod_poly_compose_mod_horner_name,
    nmod_poly_compose_series_name,
    nmod_poly_cosh_series_name,
    nmod_poly_cos_series_name,
    nmod_poly_deflate_name,
    nmod_poly_derivative_name,
    nmod_poly_discriminant_name,
    nmod_poly_div_name,
    nmod_poly_divides_name,
    nmod_poly_divides_classical_name,
    nmod_poly_div_newton_n_preinv_name,
    nmod_poly_divrem_basecase_name,
    nmod_poly_divrem_name,
    nmod_poly_divrem_newton_n_preinv_name,
    nmod_poly_div_root_name,
    nmod_poly_div_series_basecase_name,
    nmod_poly_div_series_name,
    nmod_poly_equal_trunc_name,
    nmod_poly_evaluate_mat_horner_name,
    nmod_poly_evaluate_mat_paterson_stockmeyer_name,
    nmod_poly_evaluate_nmod_name,
    nmod_poly_evaluate_nmod_vec_fast_name,
    nmod_poly_exp_series_name,
    nmod_poly_find_distinct_nonzero_roots_name,
    nmod_poly_fread_print_name,
    nmod_poly_gcd_name,
    nmod_poly_gcd_euclidean_name,
    nmod_poly_gcd_hgcd_name,
    nmod_poly_gcdinv_name,
    nmod_poly_get_set_coeff_ui_name,
    nmod_poly_get_set_str_name,
    nmod_poly_hgcd_name,
    nmod_poly_inflate_name,
    nmod_poly_init_realloc_clear_name,
    nmod_poly_integral_name,
    nmod_poly_interpolate_nmod_vec_barycentric_name,
    nmod_poly_interpolate_nmod_vec_name,
    nmod_poly_interpolate_nmod_vec_fast_name,
    nmod_poly_interpolate_nmod_vec_newton_name,
    nmod_poly_invmod_name,
    nmod_poly_inv_series_basecase_name,
    nmod_poly_inv_series_newton_name,
    nmod_poly_invsqrt_series_name,
    nmod_poly_log_series_name,
    nmod_poly_make_monic_name,
    nmod_poly_mul_name,
    nmod_poly_mul_classical_name,
    nmod_poly_mulhigh_name,
    nmod_poly_mulhigh_classical_name,
    nmod_poly_mul_KS2_name,
    nmod_poly_mul_KS4_name,
    nmod_poly_mul_KS_name,
    nmod_poly_mullow_name,
    nmod_poly_mullow_classical_name,
    nmod_poly_mullow_KS_name,
    nmod_poly_mulmod_name,
    nmod_poly_mulmod_preinv_name,
    nmod_poly_multi_crt_name,
    nmod_poly_neg_name,
    nmod_poly_pow_binexp_name,
    nmod_poly_pow_name,
    nmod_poly_powers_mod_bsgs_name,
    nmod_poly_powers_mod_naive_name,
    nmod_poly_power_sums_name,
    nmod_poly_power_sums_naive_name,
    nmod_poly_power_sums_schoenhage_name,
    nmod_poly_powmod_fmpz_binexp_name,
    nmod_poly_powmod_fmpz_binexp_preinv_name,
    nmod_poly_powmod_ui_binexp_name,
    nmod_poly_powmod_ui_binexp_preinv_name,
    nmod_poly_powmod_x_fmpz_preinv_name,
    nmod_poly_powmod_x_ui_preinv_name,
    nmod_poly_pow_trunc_binexp_name,
    nmod_poly_pow_trunc_name,
    nmod_poly_product_roots_nmod_vec_name,
    nmod_poly_rem_name,
    nmod_poly_resultant_name,
    nmod_poly_resultant_euclidean_name,
    nmod_poly_resultant_hgcd_name,
    nmod_poly_reverse_name,
    nmod_poly_revert_series_name,
    nmod_poly_revert_series_lagrange_name,
    nmod_poly_revert_series_lagrange_fast_name,
    nmod_poly_revert_series_newton_name,
    nmod_poly_scalar_addmul_nmod_name,
    nmod_poly_scalar_mul_nmod_name,
    nmod_poly_set_trunc_name,
    nmod_poly_shift_left_right_name,
    nmod_poly_sinh_series_name,
    nmod_poly_sin_series_name,
    nmod_poly_sqrt_name,
    nmod_poly_sqrt_series_name,
    nmod_poly_sub_name,
    nmod_poly_sub_series_name,
    nmod_poly_tanh_series_name,
    nmod_poly_tan_series_name,
    nmod_poly_taylor_shift_name,
    nmod_poly_taylor_shift_convolution_name,
    nmod_poly_taylor_shift_horner_name,
    nmod_poly_xgcd_name,
    nmod_poly_xgcd_euclidean_name,
    nmod_poly_xgcd_hgcd_name
};

/* main function *************************************************************/

#define NUMBER_OF_TESTS (sizeof(test_functions) / sizeof(int (*)(void)))

int
main(int argc, char ** argv)
{
    int ix, jx;

    if (argc < 2)
    {
        for (ix = 0; ix < NUMBER_OF_TESTS; ix++)
            if (test_functions[ix]())
                flint_abort();
    }
    else
    {
        for (ix = 1; ix < argc; ix++)
        {
            for (jx = 0; jx < NUMBER_OF_TESTS; jx++)
            {
                /* If argument equals to test name, run it */
                if (strcmp(argv[ix], test_names[jx]) == 0)
                {
                    if (test_functions[jx]())
                        flint_abort();
                    break;
                }
            }

            if (jx == NUMBER_OF_TESTS)
            {
                fprintf(stderr, "Error: Could not find test function for %s\n", argv[ix]);
                flint_abort();
            }
        }
    }

    return 0;
}
