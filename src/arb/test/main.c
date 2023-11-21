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
#include "t-const_reciprocal_fibonacci.c"
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

test_struct tests[] =
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
    TEST_FUNCTION(arb_const_reciprocal_fibonacci),
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

/* main function *************************************************************/

TEST_MAIN(tests)
