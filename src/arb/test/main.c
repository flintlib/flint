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

#include "t-acos.c"
#include "t-acosh.c"
#include "t-add_arf.c"
#include "t-add.c"
#include "t-add_error.c"
#include "t-add_fmpz_2exp.c"
#include "t-add_fmpz.c"
#include "t-addmul_arf.c"
#include "t-addmul.c"
#include "t-addmul_fmpz.c"
#include "t-addmul_si.c"
#include "t-addmul_ui.c"
#include "t-add_si.c"
#include "t-add_ui.c"
#include "t-agm.c"
#include "t-approx_dot.c"
#include "t-asin.c"
#include "t-asinh.c"
#include "t-atan2.c"
#include "t-atan_arf_bb.c"
#include "t-atan_arf.c"
#include "t-atan_arf_newton.c"
#include "t-atan.c"
#include "t-atan_frac_bsplit.c"
#include "t-atan_gauss_primes_vec_bsplit.c"
#include "t-atanh.c"
#include "t-atan_newton.c"
#include "t-atan_sum_bs_powtab.c"
#include "t-atan_tab.c"
#include "t-atan_taylor_rs.c"
#include "t-bell_fmpz.c"
#include "t-bell_sum_taylor.c"
#include "t-bernoulli_poly_ui.c"
#include "t-bernoulli_ui.c"
#include "t-can_round_mpfr.c"
#include "t-ceil.c"
#include "t-chebyshev_t_ui.c"
#include "t-chebyshev_u_ui.c"
#include "t-const_apery.c"
#include "t-const_catalan.c"
#include "t-const_e.c"
#include "t-const_euler.c"
#include "t-const_glaisher.c"
#include "t-const_khinchin.c"
#include "t-const_log10.c"
#include "t-const_log2.c"
#include "t-const_pi.c"
#include "t-contains_arf.c"
#include "t-contains.c"
#include "t-contains_fmpq.c"
#include "t-contains_int.c"
#include "t-contains_interior.c"
#include "t-cos.c"
#include "t-cosh.c"
#include "t-cos_pi.c"
#include "t-cos_pi_fmpq_algebraic.c"
#include "t-cos_pi_fmpq.c"
#include "t-coth.c"
#include "t-cot_pi.c"
#include "t-csc.c"
#include "t-csch.c"
#include "t-csc_pi.c"
#include "t-digamma.c"
#include "t-digits_round_inplace.c"
#include "t-div_2expm1_ui.c"
#include "t-div_arf.c"
#include "t-div.c"
#include "t-div_fmpz.c"
#include "t-div_newton.c"
#include "t-div_si.c"
#include "t-div_ui.c"
#include "t-dot.c"
#include "t-dot_fmpz.c"
#include "t-dot_si.c"
#include "t-dot_siui.c"
#include "t-dot_ui.c"
#include "t-dot_uiui.c"
#include "t-doublefac_ui.c"
#include "t-dump_file.c"
#include "t-dump_str.c"
#include "t-euler_number_fmpz.c"
#include "t-euler_number_ui.c"
#include "t-exp_arf_bb.c"
#include "t-exp_arf_rs_generic.c"
#include "t-exp.c"
#include "t-exp_invexp.c"
#include "t-expm1.c"
#include "t-exp_sum_bs_powtab.c"
#include "t-exp_tab.c"
#include "t-exp_taylor_rs.c"
#include "t-fac_ui.c"
#include "t-fib.c"
#include "t-floor.c"
#include "t-fma.c"
#include "t-gamma.c"
#include "t-gamma_fmpq.c"
#include "t-get_abs_lbound_arf.c"
#include "t-get_fmpz_mid_rad_10exp.c"
#include "t-get_interval_arf.c"
#include "t-get_interval_fmpz_2exp.c"
#include "t-get_interval_mpfr.c"
#include "t-get_lbound_arf.c"
#include "t-get_mag.c"
#include "t-get_mag_lower.c"
#include "t-get_mag_lower_nonnegative.c"
#include "t-get_mpn_fixed_mod_log2.c"
#include "t-get_mpn_fixed_mod_pi4.c"
#include "t-get_rand_fmpq.c"
#include "t-get_str.c"
#include "t-get_unique_fmpz.c"
#include "t-hurwitz_zeta.c"
#include "t-intersection.c"
#include "t-lambertw.c"
#include "t-lgamma.c"
#include "t-log1p.c"
#include "t-log_arf.c"
#include "t-log_base_ui.c"
#include "t-log.c"
#include "t-log_hypot.c"
#include "t-log_newton.c"
#include "t-log_primes_vec_bsplit.c"
#include "t-log_tab.c"
#include "t-log_ui_from_prev.c"
#include "t-max.c"
#include "t-min.c"
#include "t-minmax.c"
#include "t-mul_arf.c"
#include "t-mul.c"
#include "t-mul_fmpz.c"
#include "t-mul_more.c"
#include "t-mul_si.c"
#include "t-mul_ui.c"
#include "t-nonnegative_abs.c"
#include "t-overlaps.c"
#include "t-partitions_fmpz.c"
#include "t-pos_times_posinf.c"
#include "t-pow.c"
#include "t-power_sum_vec.c"
#include "t-pow_fmpq.c"
#include "t-pow_fmpz.c"
#include "t-primorial.c"
#include "t-rgamma.c"
#include "t-richcmp.c"
#include "t-rising2_ui.c"
#include "t-rising_ui.c"
#include "t-root_ui.c"
#include "t-rsqrt.c"
#include "t-sec.c"
#include "t-sech.c"
#include "t-set_interval_arf.c"
#include "t-set_interval_mag.c"
#include "t-set_interval_mpfr.c"
#include "t-set_interval_neg_pos_mag.c"
#include "t-set_str.c"
#include "t-sgn.c"
#include "t-sin.c"
#include "t-sinc.c"
#include "t-sin_cos_arf_bb.c"
#include "t-sin_cos_arf_generic.c"
#include "t-sin_cos.c"
#include "t-sin_cos_generic.c"
#include "t-sin_cos_pi.c"
#include "t-sin_cos_pi_fmpq_algebraic.c"
#include "t-sin_cos_pi_fmpq.c"
#include "t-sin_cos_tab.c"
#include "t-sin_cos_taylor_rs.c"
#include "t-sinc_pi.c"
#include "t-sinh.c"
#include "t-sinh_cosh.c"
#include "t-sin_pi.c"
#include "t-sin_pi_fmpq_algebraic.c"
#include "t-sin_pi_fmpq.c"
#include "t-special.c"
#include "t-sqrt1pm1.c"
#include "t-sqrt.c"
#include "t-sqrt_newton.c"
#include "t-sqrtpos.c"
#include "t-sub_arf.c"
#include "t-sub.c"
#include "t-sub_fmpz.c"
#include "t-submul_arf.c"
#include "t-submul.c"
#include "t-submul_fmpz.c"
#include "t-submul_si.c"
#include "t-submul_ui.c"
#include "t-sub_si.c"
#include "t-sub_ui.c"
#include "t-tanh.c"
#include "t-tan_pi.c"
#include "t-trim.c"
#include "t-ui_pow_ui.c"
#include "t-union.c"
#include "t-urandom.c"
#include "t-zeta.c"
#include "t-zeta_ui_asymp.c"
#include "t-zeta_ui_bernoulli.c"
#include "t-zeta_ui_borwein_bsplit.c"
#include "t-zeta_ui.c"
#include "t-zeta_ui_euler_product.c"
#include "t-zeta_ui_vec_borwein.c"
#include "t-zeta_ui_vec.c"

/* Array of test functions ***************************************************/

int (*test_functions[])(void) =
{
    TEST_FUNCTION(arb_acos),
    TEST_FUNCTION(arb_acosh),
    TEST_FUNCTION(arb_add_arf),
    TEST_FUNCTION(arb_add),
    TEST_FUNCTION(arb_add_error),
    TEST_FUNCTION(arb_add_fmpz_2exp),
    TEST_FUNCTION(arb_add_fmpz),
    TEST_FUNCTION(arb_addmul_arf),
    TEST_FUNCTION(arb_addmul),
    TEST_FUNCTION(arb_addmul_fmpz),
    TEST_FUNCTION(arb_addmul_si),
    TEST_FUNCTION(arb_addmul_ui),
    TEST_FUNCTION(arb_add_si),
    TEST_FUNCTION(arb_add_ui),
    TEST_FUNCTION(arb_agm),
    TEST_FUNCTION(arb_approx_dot),
    TEST_FUNCTION(arb_asin),
    TEST_FUNCTION(arb_asinh),
    TEST_FUNCTION(arb_atan2),
    TEST_FUNCTION(arb_atan_arf_bb),
    TEST_FUNCTION(arb_atan_arf),
    TEST_FUNCTION(arb_atan_arf_newton),
    TEST_FUNCTION(arb_atan),
    TEST_FUNCTION(arb_atan_frac_bsplit),
    TEST_FUNCTION(arb_atan_gauss_primes_vec_bsplit),
    TEST_FUNCTION(arb_atanh),
    TEST_FUNCTION(arb_atan_newton),
    TEST_FUNCTION(arb_atan_sum_bs_powtab),
    TEST_FUNCTION(arb_atan_tab),
    TEST_FUNCTION(arb_atan_taylor_rs),
    TEST_FUNCTION(arb_bell_fmpz),
    TEST_FUNCTION(arb_bell_sum_taylor),
    TEST_FUNCTION(arb_bernoulli_poly_ui),
    TEST_FUNCTION(arb_bernoulli_ui),
    TEST_FUNCTION(arb_can_round_mpfr),
    TEST_FUNCTION(arb_ceil),
    TEST_FUNCTION(arb_chebyshev_t_ui),
    TEST_FUNCTION(arb_chebyshev_u_ui),
    TEST_FUNCTION(arb_const_apery),
    TEST_FUNCTION(arb_const_catalan),
    TEST_FUNCTION(arb_const_e),
    TEST_FUNCTION(arb_const_euler),
    TEST_FUNCTION(arb_const_glaisher),
    TEST_FUNCTION(arb_const_khinchin),
    TEST_FUNCTION(arb_const_log10),
    TEST_FUNCTION(arb_const_log2),
    TEST_FUNCTION(arb_const_pi),
    TEST_FUNCTION(arb_contains_arf),
    TEST_FUNCTION(arb_contains),
    TEST_FUNCTION(arb_contains_fmpq),
    TEST_FUNCTION(arb_contains_int),
    TEST_FUNCTION(arb_contains_interior),
    TEST_FUNCTION(arb_cos),
    TEST_FUNCTION(arb_cosh),
    TEST_FUNCTION(arb_cos_pi),
    TEST_FUNCTION(arb_cos_pi_fmpq_algebraic),
    TEST_FUNCTION(arb_cos_pi_fmpq),
    TEST_FUNCTION(arb_coth),
    TEST_FUNCTION(arb_cot_pi),
    TEST_FUNCTION(arb_csc),
    TEST_FUNCTION(arb_csch),
    TEST_FUNCTION(arb_csc_pi),
    TEST_FUNCTION(arb_digamma),
    TEST_FUNCTION(arb_digits_round_inplace),
    TEST_FUNCTION(arb_div_2expm1_ui),
    TEST_FUNCTION(arb_div_arf),
    TEST_FUNCTION(arb_div),
    TEST_FUNCTION(arb_div_fmpz),
    TEST_FUNCTION(arb_div_newton),
    TEST_FUNCTION(arb_div_si),
    TEST_FUNCTION(arb_div_ui),
    TEST_FUNCTION(arb_dot),
    TEST_FUNCTION(arb_dot_fmpz),
    TEST_FUNCTION(arb_dot_si),
    TEST_FUNCTION(arb_dot_siui),
    TEST_FUNCTION(arb_dot_ui),
    TEST_FUNCTION(arb_dot_uiui),
    TEST_FUNCTION(arb_doublefac_ui),
    TEST_FUNCTION(arb_dump_file),
    TEST_FUNCTION(arb_dump_str),
    TEST_FUNCTION(arb_euler_number_fmpz),
    TEST_FUNCTION(arb_euler_number_ui),
    TEST_FUNCTION(arb_exp_arf_bb),
    TEST_FUNCTION(arb_exp_arf_rs_generic),
    TEST_FUNCTION(arb_exp),
    TEST_FUNCTION(arb_exp_invexp),
    TEST_FUNCTION(arb_expm1),
    TEST_FUNCTION(arb_exp_sum_bs_powtab),
    TEST_FUNCTION(arb_exp_tab),
    TEST_FUNCTION(arb_exp_taylor_rs),
    TEST_FUNCTION(arb_fac_ui),
    TEST_FUNCTION(arb_fib),
    TEST_FUNCTION(arb_floor),
    TEST_FUNCTION(arb_fma),
    TEST_FUNCTION(arb_gamma),
    TEST_FUNCTION(arb_gamma_fmpq),
    TEST_FUNCTION(arb_get_abs_lbound_arf),
    TEST_FUNCTION(arb_get_fmpz_mid_rad_10exp),
    TEST_FUNCTION(arb_get_interval_arf),
    TEST_FUNCTION(arb_get_interval_fmpz_2exp),
    TEST_FUNCTION(arb_get_interval_mpfr),
    TEST_FUNCTION(arb_get_lbound_arf),
    TEST_FUNCTION(arb_get_mag),
    TEST_FUNCTION(arb_get_mag_lower),
    TEST_FUNCTION(arb_get_mag_lower_nonnegative),
    TEST_FUNCTION(arb_get_mpn_fixed_mod_log2),
    TEST_FUNCTION(arb_get_mpn_fixed_mod_pi4),
    TEST_FUNCTION(arb_get_rand_fmpq),
    TEST_FUNCTION(arb_get_str),
    TEST_FUNCTION(arb_get_unique_fmpz),
    TEST_FUNCTION(arb_hurwitz_zeta),
    TEST_FUNCTION(arb_intersection),
    TEST_FUNCTION(arb_lambertw),
    TEST_FUNCTION(arb_lgamma),
    TEST_FUNCTION(arb_log1p),
    TEST_FUNCTION(arb_log_arf),
    TEST_FUNCTION(arb_log_base_ui),
    TEST_FUNCTION(arb_log),
    TEST_FUNCTION(arb_log_hypot),
    TEST_FUNCTION(arb_log_newton),
    TEST_FUNCTION(arb_log_primes_vec_bsplit),
    TEST_FUNCTION(arb_log_tab),
    TEST_FUNCTION(arb_log_ui_from_prev),
    TEST_FUNCTION(arb_max),
    TEST_FUNCTION(arb_min),
    TEST_FUNCTION(arb_minmax),
    TEST_FUNCTION(arb_mul_arf),
    TEST_FUNCTION(arb_mul),
    TEST_FUNCTION(arb_mul_fmpz),
    TEST_FUNCTION(arb_mul_more),
    TEST_FUNCTION(arb_mul_si),
    TEST_FUNCTION(arb_mul_ui),
    TEST_FUNCTION(arb_nonnegative_abs),
    TEST_FUNCTION(arb_overlaps),
    TEST_FUNCTION(arb_partitions_fmpz),
    TEST_FUNCTION(arb_pos_times_posinf),
    TEST_FUNCTION(arb_pow),
    TEST_FUNCTION(arb_power_sum_vec),
    TEST_FUNCTION(arb_pow_fmpq),
    TEST_FUNCTION(arb_pow_fmpz),
    TEST_FUNCTION(arb_primorial),
    TEST_FUNCTION(arb_rgamma),
    TEST_FUNCTION(arb_richcmp),
    TEST_FUNCTION(arb_rising2_ui),
    TEST_FUNCTION(arb_rising_ui),
    TEST_FUNCTION(arb_root_ui),
    TEST_FUNCTION(arb_rsqrt),
    TEST_FUNCTION(arb_sec),
    TEST_FUNCTION(arb_sech),
    TEST_FUNCTION(arb_set_interval_arf),
    TEST_FUNCTION(arb_set_interval_mag),
    TEST_FUNCTION(arb_set_interval_mpfr),
    TEST_FUNCTION(arb_set_interval_neg_pos_mag),
    TEST_FUNCTION(arb_set_str),
    TEST_FUNCTION(arb_sgn),
    TEST_FUNCTION(arb_sin),
    TEST_FUNCTION(arb_sinc),
    TEST_FUNCTION(arb_sin_cos_arf_bb),
    TEST_FUNCTION(arb_sin_cos_arf_generic),
    TEST_FUNCTION(arb_sin_cos),
    TEST_FUNCTION(arb_sin_cos_generic),
    TEST_FUNCTION(arb_sin_cos_pi),
    TEST_FUNCTION(arb_sin_cos_pi_fmpq_algebraic),
    TEST_FUNCTION(arb_sin_cos_pi_fmpq),
    TEST_FUNCTION(arb_sin_cos_tab),
    TEST_FUNCTION(arb_sin_cos_taylor_rs),
    TEST_FUNCTION(arb_sinc_pi),
    TEST_FUNCTION(arb_sinh),
    TEST_FUNCTION(arb_sinh_cosh),
    TEST_FUNCTION(arb_sin_pi),
    TEST_FUNCTION(arb_sin_pi_fmpq_algebraic),
    TEST_FUNCTION(arb_sin_pi_fmpq),
    TEST_FUNCTION(arb_special),
    TEST_FUNCTION(arb_sqrt1pm1),
    TEST_FUNCTION(arb_sqrt),
    TEST_FUNCTION(arb_sqrt_newton),
    TEST_FUNCTION(arb_sqrtpos),
    TEST_FUNCTION(arb_sub_arf),
    TEST_FUNCTION(arb_sub),
    TEST_FUNCTION(arb_sub_fmpz),
    TEST_FUNCTION(arb_submul_arf),
    TEST_FUNCTION(arb_submul),
    TEST_FUNCTION(arb_submul_fmpz),
    TEST_FUNCTION(arb_submul_si),
    TEST_FUNCTION(arb_submul_ui),
    TEST_FUNCTION(arb_sub_si),
    TEST_FUNCTION(arb_sub_ui),
    TEST_FUNCTION(arb_tanh),
    TEST_FUNCTION(arb_tan_pi),
    TEST_FUNCTION(arb_trim),
    TEST_FUNCTION(arb_ui_pow_ui),
    TEST_FUNCTION(arb_union),
    TEST_FUNCTION(arb_urandom),
    TEST_FUNCTION(arb_zeta),
    TEST_FUNCTION(arb_zeta_ui_asymp),
    TEST_FUNCTION(arb_zeta_ui_bernoulli),
    TEST_FUNCTION(arb_zeta_ui_borwein_bsplit),
    TEST_FUNCTION(arb_zeta_ui),
    TEST_FUNCTION(arb_zeta_ui_euler_product),
    TEST_FUNCTION(arb_zeta_ui_vec_borwein),
    TEST_FUNCTION(arb_zeta_ui_vec)
};

char arb_acos_name[] = "arb_acos";
char arb_acosh_name[] = "arb_acosh";
char arb_add_arf_name[] = "arb_add_arf";
char arb_add_name[] = "arb_add";
char arb_add_error_name[] = "arb_add_error";
char arb_add_fmpz_2exp_name[] = "arb_add_fmpz_2exp";
char arb_add_fmpz_name[] = "arb_add_fmpz";
char arb_addmul_arf_name[] = "arb_addmul_arf";
char arb_addmul_name[] = "arb_addmul";
char arb_addmul_fmpz_name[] = "arb_addmul_fmpz";
char arb_addmul_si_name[] = "arb_addmul_si";
char arb_addmul_ui_name[] = "arb_addmul_ui";
char arb_add_si_name[] = "arb_add_si";
char arb_add_ui_name[] = "arb_add_ui";
char arb_agm_name[] = "arb_agm";
char arb_approx_dot_name[] = "arb_approx_dot";
char arb_asin_name[] = "arb_asin";
char arb_asinh_name[] = "arb_asinh";
char arb_atan2_name[] = "arb_atan2";
char arb_atan_arf_bb_name[] = "arb_atan_arf_bb";
char arb_atan_arf_name[] = "arb_atan_arf";
char arb_atan_arf_newton_name[] = "arb_atan_arf_newton";
char arb_atan_name[] = "arb_atan";
char arb_atan_frac_bsplit_name[] = "arb_atan_frac_bsplit";
char arb_atan_gauss_primes_vec_bsplit_name[] = "arb_atan_gauss_primes_vec_bsplit";
char arb_atanh_name[] = "arb_atanh";
char arb_atan_newton_name[] = "arb_atan_newton";
char arb_atan_sum_bs_powtab_name[] = "arb_atan_sum_bs_powtab";
char arb_atan_tab_name[] = "arb_atan_tab";
char arb_atan_taylor_rs_name[] = "arb_atan_taylor_rs";
char arb_bell_fmpz_name[] = "arb_bell_fmpz";
char arb_bell_sum_taylor_name[] = "arb_bell_sum_taylor";
char arb_bernoulli_poly_ui_name[] = "arb_bernoulli_poly_ui";
char arb_bernoulli_ui_name[] = "arb_bernoulli_ui";
char arb_can_round_mpfr_name[] = "arb_can_round_mpfr";
char arb_ceil_name[] = "arb_ceil";
char arb_chebyshev_t_ui_name[] = "arb_chebyshev_t_ui";
char arb_chebyshev_u_ui_name[] = "arb_chebyshev_u_ui";
char arb_const_apery_name[] = "arb_const_apery";
char arb_const_catalan_name[] = "arb_const_catalan";
char arb_const_e_name[] = "arb_const_e";
char arb_const_euler_name[] = "arb_const_euler";
char arb_const_glaisher_name[] = "arb_const_glaisher";
char arb_const_khinchin_name[] = "arb_const_khinchin";
char arb_const_log10_name[] = "arb_const_log10";
char arb_const_log2_name[] = "arb_const_log2";
char arb_const_pi_name[] = "arb_const_pi";
char arb_contains_arf_name[] = "arb_contains_arf";
char arb_contains_name[] = "arb_contains";
char arb_contains_fmpq_name[] = "arb_contains_fmpq";
char arb_contains_int_name[] = "arb_contains_int";
char arb_contains_interior_name[] = "arb_contains_interior";
char arb_cos_name[] = "arb_cos";
char arb_cosh_name[] = "arb_cosh";
char arb_cos_pi_name[] = "arb_cos_pi";
char arb_cos_pi_fmpq_algebraic_name[] = "arb_cos_pi_fmpq_algebraic";
char arb_cos_pi_fmpq_name[] = "arb_cos_pi_fmpq";
char arb_coth_name[] = "arb_coth";
char arb_cot_pi_name[] = "arb_cot_pi";
char arb_csc_name[] = "arb_csc";
char arb_csch_name[] = "arb_csch";
char arb_csc_pi_name[] = "arb_csc_pi";
char arb_digamma_name[] = "arb_digamma";
char arb_digits_round_inplace_name[] = "arb_digits_round_inplace";
char arb_div_2expm1_ui_name[] = "arb_div_2expm1_ui";
char arb_div_arf_name[] = "arb_div_arf";
char arb_div_name[] = "arb_div";
char arb_div_fmpz_name[] = "arb_div_fmpz";
char arb_div_newton_name[] = "arb_div_newton";
char arb_div_si_name[] = "arb_div_si";
char arb_div_ui_name[] = "arb_div_ui";
char arb_dot_name[] = "arb_dot";
char arb_dot_fmpz_name[] = "arb_dot_fmpz";
char arb_dot_si_name[] = "arb_dot_si";
char arb_dot_siui_name[] = "arb_dot_siui";
char arb_dot_ui_name[] = "arb_dot_ui";
char arb_dot_uiui_name[] = "arb_dot_uiui";
char arb_doublefac_ui_name[] = "arb_doublefac_ui";
char arb_dump_file_name[] = "arb_dump_file";
char arb_dump_str_name[] = "arb_dump_str";
char arb_euler_number_fmpz_name[] = "arb_euler_number_fmpz";
char arb_euler_number_ui_name[] = "arb_euler_number_ui";
char arb_exp_arf_bb_name[] = "arb_exp_arf_bb";
char arb_exp_arf_rs_generic_name[] = "arb_exp_arf_rs_generic";
char arb_exp_name[] = "arb_exp";
char arb_exp_invexp_name[] = "arb_exp_invexp";
char arb_expm1_name[] = "arb_expm1";
char arb_exp_sum_bs_powtab_name[] = "arb_exp_sum_bs_powtab";
char arb_exp_tab_name[] = "arb_exp_tab";
char arb_exp_taylor_rs_name[] = "arb_exp_taylor_rs";
char arb_fac_ui_name[] = "arb_fac_ui";
char arb_fib_name[] = "arb_fib";
char arb_floor_name[] = "arb_floor";
char arb_fma_name[] = "arb_fma";
char arb_gamma_name[] = "arb_gamma";
char arb_gamma_fmpq_name[] = "arb_gamma_fmpq";
char arb_get_abs_lbound_arf_name[] = "arb_get_abs_lbound_arf";
char arb_get_fmpz_mid_rad_10exp_name[] = "arb_get_fmpz_mid_rad_10exp";
char arb_get_interval_arf_name[] = "arb_get_interval_arf";
char arb_get_interval_fmpz_2exp_name[] = "arb_get_interval_fmpz_2exp";
char arb_get_interval_mpfr_name[] = "arb_get_interval_mpfr";
char arb_get_lbound_arf_name[] = "arb_get_lbound_arf";
char arb_get_mag_name[] = "arb_get_mag";
char arb_get_mag_lower_name[] = "arb_get_mag_lower";
char arb_get_mag_lower_nonnegative_name[] = "arb_get_mag_lower_nonnegative";
char arb_get_mpn_fixed_mod_log2_name[] = "arb_get_mpn_fixed_mod_log2";
char arb_get_mpn_fixed_mod_pi4_name[] = "arb_get_mpn_fixed_mod_pi4";
char arb_get_rand_fmpq_name[] = "arb_get_rand_fmpq";
char arb_get_str_name[] = "arb_get_str";
char arb_get_unique_fmpz_name[] = "arb_get_unique_fmpz";
char arb_hurwitz_zeta_name[] = "arb_hurwitz_zeta";
char arb_intersection_name[] = "arb_intersection";
char arb_lambertw_name[] = "arb_lambertw";
char arb_lgamma_name[] = "arb_lgamma";
char arb_log1p_name[] = "arb_log1p";
char arb_log_arf_name[] = "arb_log_arf";
char arb_log_base_ui_name[] = "arb_log_base_ui";
char arb_log_name[] = "arb_log";
char arb_log_hypot_name[] = "arb_log_hypot";
char arb_log_newton_name[] = "arb_log_newton";
char arb_log_primes_vec_bsplit_name[] = "arb_log_primes_vec_bsplit";
char arb_log_tab_name[] = "arb_log_tab";
char arb_log_ui_from_prev_name[] = "arb_log_ui_from_prev";
char arb_max_name[] = "arb_max";
char arb_min_name[] = "arb_min";
char arb_minmax_name[] = "arb_minmax";
char arb_mul_arf_name[] = "arb_mul_arf";
char arb_mul_name[] = "arb_mul";
char arb_mul_fmpz_name[] = "arb_mul_fmpz";
char arb_mul_more_name[] = "arb_mul_more";
char arb_mul_si_name[] = "arb_mul_si";
char arb_mul_ui_name[] = "arb_mul_ui";
char arb_nonnegative_abs_name[] = "arb_nonnegative_abs";
char arb_overlaps_name[] = "arb_overlaps";
char arb_partitions_fmpz_name[] = "arb_partitions_fmpz";
char arb_pos_times_posinf_name[] = "arb_pos_times_posinf";
char arb_pow_name[] = "arb_pow";
char arb_power_sum_vec_name[] = "arb_power_sum_vec";
char arb_pow_fmpq_name[] = "arb_pow_fmpq";
char arb_pow_fmpz_name[] = "arb_pow_fmpz";
char arb_primorial_name[] = "arb_primorial";
char arb_rgamma_name[] = "arb_rgamma";
char arb_richcmp_name[] = "arb_richcmp";
char arb_rising2_ui_name[] = "arb_rising2_ui";
char arb_rising_ui_name[] = "arb_rising_ui";
char arb_root_ui_name[] = "arb_root_ui";
char arb_rsqrt_name[] = "arb_rsqrt";
char arb_sec_name[] = "arb_sec";
char arb_sech_name[] = "arb_sech";
char arb_set_interval_arf_name[] = "arb_set_interval_arf";
char arb_set_interval_mag_name[] = "arb_set_interval_mag";
char arb_set_interval_mpfr_name[] = "arb_set_interval_mpfr";
char arb_set_interval_neg_pos_mag_name[] = "arb_set_interval_neg_pos_mag";
char arb_set_str_name[] = "arb_set_str";
char arb_sgn_name[] = "arb_sgn";
char arb_sin_name[] = "arb_sin";
char arb_sinc_name[] = "arb_sinc";
char arb_sin_cos_arf_bb_name[] = "arb_sin_cos_arf_bb";
char arb_sin_cos_arf_generic_name[] = "arb_sin_cos_arf_generic";
char arb_sin_cos_name[] = "arb_sin_cos";
char arb_sin_cos_generic_name[] = "arb_sin_cos_generic";
char arb_sin_cos_pi_name[] = "arb_sin_cos_pi";
char arb_sin_cos_pi_fmpq_algebraic_name[] = "arb_sin_cos_pi_fmpq_algebraic";
char arb_sin_cos_pi_fmpq_name[] = "arb_sin_cos_pi_fmpq";
char arb_sin_cos_tab_name[] = "arb_sin_cos_tab";
char arb_sin_cos_taylor_rs_name[] = "arb_sin_cos_taylor_rs";
char arb_sinc_pi_name[] = "arb_sinc_pi";
char arb_sinh_name[] = "arb_sinh";
char arb_sinh_cosh_name[] = "arb_sinh_cosh";
char arb_sin_pi_name[] = "arb_sin_pi";
char arb_sin_pi_fmpq_algebraic_name[] = "arb_sin_pi_fmpq_algebraic";
char arb_sin_pi_fmpq_name[] = "arb_sin_pi_fmpq";
char arb_special_name[] = "arb_special";
char arb_sqrt1pm1_name[] = "arb_sqrt1pm1";
char arb_sqrt_name[] = "arb_sqrt";
char arb_sqrt_newton_name[] = "arb_sqrt_newton";
char arb_sqrtpos_name[] = "arb_sqrtpos";
char arb_sub_arf_name[] = "arb_sub_arf";
char arb_sub_name[] = "arb_sub";
char arb_sub_fmpz_name[] = "arb_sub_fmpz";
char arb_submul_arf_name[] = "arb_submul_arf";
char arb_submul_name[] = "arb_submul";
char arb_submul_fmpz_name[] = "arb_submul_fmpz";
char arb_submul_si_name[] = "arb_submul_si";
char arb_submul_ui_name[] = "arb_submul_ui";
char arb_sub_si_name[] = "arb_sub_si";
char arb_sub_ui_name[] = "arb_sub_ui";
char arb_tanh_name[] = "arb_tanh";
char arb_tan_pi_name[] = "arb_tan_pi";
char arb_trim_name[] = "arb_trim";
char arb_ui_pow_ui_name[] = "arb_ui_pow_ui";
char arb_union_name[] = "arb_union";
char arb_urandom_name[] = "arb_urandom";
char arb_zeta_name[] = "arb_zeta";
char arb_zeta_ui_asymp_name[] = "arb_zeta_ui_asymp";
char arb_zeta_ui_bernoulli_name[] = "arb_zeta_ui_bernoulli";
char arb_zeta_ui_borwein_bsplit_name[] = "arb_zeta_ui_borwein_bsplit";
char arb_zeta_ui_name[] = "arb_zeta_ui";
char arb_zeta_ui_euler_product_name[] = "arb_zeta_ui_euler_product";
char arb_zeta_ui_vec_borwein_name[] = "arb_zeta_ui_vec_borwein";
char arb_zeta_ui_vec_name[] = "arb_zeta_ui_vec";

char * test_names[] =
{
    arb_acos_name,
    arb_acosh_name,
    arb_add_arf_name,
    arb_add_name,
    arb_add_error_name,
    arb_add_fmpz_2exp_name,
    arb_add_fmpz_name,
    arb_addmul_arf_name,
    arb_addmul_name,
    arb_addmul_fmpz_name,
    arb_addmul_si_name,
    arb_addmul_ui_name,
    arb_add_si_name,
    arb_add_ui_name,
    arb_agm_name,
    arb_approx_dot_name,
    arb_asin_name,
    arb_asinh_name,
    arb_atan2_name,
    arb_atan_arf_bb_name,
    arb_atan_arf_name,
    arb_atan_arf_newton_name,
    arb_atan_name,
    arb_atan_frac_bsplit_name,
    arb_atan_gauss_primes_vec_bsplit_name,
    arb_atanh_name,
    arb_atan_newton_name,
    arb_atan_sum_bs_powtab_name,
    arb_atan_tab_name,
    arb_atan_taylor_rs_name,
    arb_bell_fmpz_name,
    arb_bell_sum_taylor_name,
    arb_bernoulli_poly_ui_name,
    arb_bernoulli_ui_name,
    arb_can_round_mpfr_name,
    arb_ceil_name,
    arb_chebyshev_t_ui_name,
    arb_chebyshev_u_ui_name,
    arb_const_apery_name,
    arb_const_catalan_name,
    arb_const_e_name,
    arb_const_euler_name,
    arb_const_glaisher_name,
    arb_const_khinchin_name,
    arb_const_log10_name,
    arb_const_log2_name,
    arb_const_pi_name,
    arb_contains_arf_name,
    arb_contains_name,
    arb_contains_fmpq_name,
    arb_contains_int_name,
    arb_contains_interior_name,
    arb_cos_name,
    arb_cosh_name,
    arb_cos_pi_name,
    arb_cos_pi_fmpq_algebraic_name,
    arb_cos_pi_fmpq_name,
    arb_coth_name,
    arb_cot_pi_name,
    arb_csc_name,
    arb_csch_name,
    arb_csc_pi_name,
    arb_digamma_name,
    arb_digits_round_inplace_name,
    arb_div_2expm1_ui_name,
    arb_div_arf_name,
    arb_div_name,
    arb_div_fmpz_name,
    arb_div_newton_name,
    arb_div_si_name,
    arb_div_ui_name,
    arb_dot_name,
    arb_dot_fmpz_name,
    arb_dot_si_name,
    arb_dot_siui_name,
    arb_dot_ui_name,
    arb_dot_uiui_name,
    arb_doublefac_ui_name,
    arb_dump_file_name,
    arb_dump_str_name,
    arb_euler_number_fmpz_name,
    arb_euler_number_ui_name,
    arb_exp_arf_bb_name,
    arb_exp_arf_rs_generic_name,
    arb_exp_name,
    arb_exp_invexp_name,
    arb_expm1_name,
    arb_exp_sum_bs_powtab_name,
    arb_exp_tab_name,
    arb_exp_taylor_rs_name,
    arb_fac_ui_name,
    arb_fib_name,
    arb_floor_name,
    arb_fma_name,
    arb_gamma_name,
    arb_gamma_fmpq_name,
    arb_get_abs_lbound_arf_name,
    arb_get_fmpz_mid_rad_10exp_name,
    arb_get_interval_arf_name,
    arb_get_interval_fmpz_2exp_name,
    arb_get_interval_mpfr_name,
    arb_get_lbound_arf_name,
    arb_get_mag_name,
    arb_get_mag_lower_name,
    arb_get_mag_lower_nonnegative_name,
    arb_get_mpn_fixed_mod_log2_name,
    arb_get_mpn_fixed_mod_pi4_name,
    arb_get_rand_fmpq_name,
    arb_get_str_name,
    arb_get_unique_fmpz_name,
    arb_hurwitz_zeta_name,
    arb_intersection_name,
    arb_lambertw_name,
    arb_lgamma_name,
    arb_log1p_name,
    arb_log_arf_name,
    arb_log_base_ui_name,
    arb_log_name,
    arb_log_hypot_name,
    arb_log_newton_name,
    arb_log_primes_vec_bsplit_name,
    arb_log_tab_name,
    arb_log_ui_from_prev_name,
    arb_max_name,
    arb_min_name,
    arb_minmax_name,
    arb_mul_arf_name,
    arb_mul_name,
    arb_mul_fmpz_name,
    arb_mul_more_name,
    arb_mul_si_name,
    arb_mul_ui_name,
    arb_nonnegative_abs_name,
    arb_overlaps_name,
    arb_partitions_fmpz_name,
    arb_pos_times_posinf_name,
    arb_pow_name,
    arb_power_sum_vec_name,
    arb_pow_fmpq_name,
    arb_pow_fmpz_name,
    arb_primorial_name,
    arb_rgamma_name,
    arb_richcmp_name,
    arb_rising2_ui_name,
    arb_rising_ui_name,
    arb_root_ui_name,
    arb_rsqrt_name,
    arb_sec_name,
    arb_sech_name,
    arb_set_interval_arf_name,
    arb_set_interval_mag_name,
    arb_set_interval_mpfr_name,
    arb_set_interval_neg_pos_mag_name,
    arb_set_str_name,
    arb_sgn_name,
    arb_sin_name,
    arb_sinc_name,
    arb_sin_cos_arf_bb_name,
    arb_sin_cos_arf_generic_name,
    arb_sin_cos_name,
    arb_sin_cos_generic_name,
    arb_sin_cos_pi_name,
    arb_sin_cos_pi_fmpq_algebraic_name,
    arb_sin_cos_pi_fmpq_name,
    arb_sin_cos_tab_name,
    arb_sin_cos_taylor_rs_name,
    arb_sinc_pi_name,
    arb_sinh_name,
    arb_sinh_cosh_name,
    arb_sin_pi_name,
    arb_sin_pi_fmpq_algebraic_name,
    arb_sin_pi_fmpq_name,
    arb_special_name,
    arb_sqrt1pm1_name,
    arb_sqrt_name,
    arb_sqrt_newton_name,
    arb_sqrtpos_name,
    arb_sub_arf_name,
    arb_sub_name,
    arb_sub_fmpz_name,
    arb_submul_arf_name,
    arb_submul_name,
    arb_submul_fmpz_name,
    arb_submul_si_name,
    arb_submul_ui_name,
    arb_sub_si_name,
    arb_sub_ui_name,
    arb_tanh_name,
    arb_tan_pi_name,
    arb_trim_name,
    arb_ui_pow_ui_name,
    arb_union_name,
    arb_urandom_name,
    arb_zeta_name,
    arb_zeta_ui_asymp_name,
    arb_zeta_ui_bernoulli_name,
    arb_zeta_ui_borwein_bsplit_name,
    arb_zeta_ui_name,
    arb_zeta_ui_euler_product_name,
    arb_zeta_ui_vec_borwein_name,
    arb_zeta_ui_vec_name
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
