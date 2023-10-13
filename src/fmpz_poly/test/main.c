/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Try to get fdopen declared for fmpz_poly_[print/read] */
#if defined __STRICT_ANSI__
# undef __STRICT_ANSI__
#endif

#include <string.h>
#include <stdio.h>
#include <mpfr.h>

/* Include functions *********************************************************/

#include "t-2norm_normalised_bits.c"
#include "t-add.c"
#include "t-add_series.c"
#include "t-add_sub_fmpz.c"
#include "t-add_sub_si.c"
#include "t-bit_pack.c"
#include "t-bound_roots.c"
#include "t-chebyshev_t.c"
#include "t-chebyshev_u.c"
#include "t-CLD_bound.c"
#include "t-compose.c"
#include "t-compose_divconquer.c"
#include "t-compose_horner.c"
#include "t-compose_series_brent_kung.c"
#include "t-compose_series.c"
#include "t-compose_series_horner.c"
#include "t-content.c"
#include "t-cos_minpoly.c"
#include "t-CRT_ui.c"
#include "t-CRT_ui_unsigned.c"
#include "t-cyclotomic.c"
#include "t-deflate.c"
#include "t-derivative.c"
#include "t-discriminant.c"
#include "t-div_basecase.c"
#include "t-div_divconquer.c"
#include "t-divhigh_smodp.c"
#include "t-divides.c"
#include "t-divlow_smodp.c"
#include "t-div_preinv.c"
#include "t-divrem_basecase.c"
#include "t-divrem.c"
#include "t-divrem_divconquer.c"
#include "t-divrem_preinv.c"
#include "t-div_root.c"
#include "t-div_series_basecase.c"
#include "t-div_series.c"
#include "t-div_series_divconquer.c"
#include "t-equal_fmpz.c"
#include "t-equal_trunc.c"
#include "t-eta_qexp.c"
#include "t-eulerian_polynomial.c"
#include "t-evaluate_divconquer_fmpq.c"
#include "t-evaluate_divconquer_fmpz.c"
#include "t-evaluate_fmpq.c"
#include "t-evaluate_fmpz.c"
#include "t-evaluate_horner_d_2exp.c"
#include "t-evaluate_horner_fmpq.c"
#include "t-evaluate_horner_fmpz.c"
#include "t-evaluate_mod.c"
#include "t-fibonacci.c"
#include "t-gcd.c"
#include "t-gcd_heuristic.c"
#include "t-gcd_modular.c"
#include "t-gcd_subresultant.c"
#include "t-get_coeff_ptr.c"
#include "t-get_nmod_poly.c"
#include "t-get_set_coeff_fmpz.c"
#include "t-get_set_coeff_si.c"
#include "t-get_set_coeff_ui.c"
#include "t-get_set_str.c"
#include "t-get_str.c"
#include "t-get_str_pretty.c"
#include "t-hensel_lift.c"
#include "t-hensel_lift_once.c"
#include "t-hensel_lift_without_only_inverse.c"
#include "t-hensel_start_continue_lift.c"
#include "t-hermite_h.c"
#include "t-hermite_he.c"
#include "t-inflate.c"
#include "t-init_realloc_clear.c"
#include "t-interpolate_fmpz_vec.c"
#include "t-inv_series_basecase.c"
#include "t-inv_series.c"
#include "t-inv_series_newton.c"
#include "t-is_cyclotomic.c"
#include "t-is_squarefree.c"
#include "t-lcm.c"
#include "t-legendre_pt.c"
#include "t-mul.c"
#include "t-mul_classical.c"
#include "t-mulhigh_classical.c"
#include "t-mulhigh_karatsuba_n.c"
#include "t-mulhigh_n.c"
#include "t-mul_karatsuba.c"
#include "t-mul_KS.c"
#include "t-mullow.c"
#include "t-mullow_classical.c"
#include "t-mullow_karatsuba_n.c"
#include "t-mullow_KS.c"
#include "t-mullow_SS.c"
#include "t-mullow_SS_precache.c"
#include "t-mulmid_classical.c"
#include "t-mul_SS.c"
#include "t-mul_SS_precache.c"
#include "t-neg.c"
#include "t-newton_to_monomial.c"
#include "t-nth_derivative.c"
#include "t-num_real_roots.c"
#include "t-num_real_roots_sturm.c"
#include "t-pow_addchains.c"
#include "t-pow_binexp.c"
#include "t-pow_binomial.c"
#include "t-pow.c"
#include "t-power_sums.c"
#include "t-pow_multinomial.c"
#include "t-pow_trunc.c"
#include "t-primitive_part.c"
#include "t-print_read.c"
#include "t-print_read_pretty.c"
#include "t-product_roots_fmpq_vec.c"
#include "t-product_roots_fmpz_vec.c"
#include "t-pseudo_div.c"
#include "t-pseudo_divrem_basecase.c"
#include "t-pseudo_divrem_cohen.c"
#include "t-pseudo_divrem_divconquer.c"
#include "t-pseudo_rem.c"
#include "t-pseudo_rem_cohen.c"
#include "t-randtest_no_real_root.c"
#include "t-rem_basecase.c"
#include "t-remove.c"
#include "t-remove_content_2exp.c"
#include "t-rem_powers_precomp.c"
#include "t-resultant.c"
#include "t-resultant_euclidean.c"
#include "t-resultant_modular.c"
#include "t-resultant_modular_div.c"
#include "t-reverse.c"
#include "t-revert_series.c"
#include "t-revert_series_lagrange.c"
#include "t-revert_series_lagrange_fast.c"
#include "t-revert_series_newton.c"
#include "t-scalar_abs.c"
#include "t-scalar_addmul_fmpz.c"
#include "t-scalar_addmul_si.c"
#include "t-scalar_addmul_ui.c"
#include "t-scalar_mul_fmpz.c"
#include "t-scalar_mul_si.c"
#include "t-scalar_mul_ui.c"
#include "t-scalar_submul_fmpz.c"
#include "t-scale_2exp.c"
#include "t-set_equal.c"
#include "t-set_fmpz_equal.c"
#include "t-set_si_equal.c"
#include "t-set_trunc.c"
#include "t-set_ui_equal.c"
#include "t-shift_left_right.c"
#include "t-signature.c"
#include "t-sqr.c"
#include "t-sqr_classical.c"
#include "t-sqr_karatsuba.c"
#include "t-sqr_KS.c"
#include "t-sqrlow.c"
#include "t-sqrlow_classical.c"
#include "t-sqrlow_karatsuba_n.c"
#include "t-sqrlow_KS.c"
#include "t-sqrt.c"
#include "t-sqrt_classical.c"
#include "t-sqrt_divconquer.c"
#include "t-sqrt_KS.c"
#include "t-sqrtrem_classical.c"
#include "t-sqrtrem_divconquer.c"
#include "t-sqrt_series.c"
#include "t-sub.c"
#include "t-sub_series.c"
#include "t-swap.c"
#include "t-swinnerton_dyer.c"
#include "t-taylor_shift.c"
#include "t-taylor_shift_divconquer.c"
#include "t-taylor_shift_horner.c"
#include "t-taylor_shift_multi_mod_threaded.c"
#include "t-theta_qexp.c"
#include "t-xgcd_modular.c"
#include "t-zero.c"
#include "t-zero_coeffs.c"

/* Array of test functions ***************************************************/

int (*test_functions[])(void) =
{
    TEST_FUNCTION(fmpz_poly_2norm_normalised_bits),
    TEST_FUNCTION(fmpz_poly_add),
    TEST_FUNCTION(fmpz_poly_add_series),
    TEST_FUNCTION(fmpz_poly_add_sub_fmpz),
    TEST_FUNCTION(fmpz_poly_add_sub_si),
    TEST_FUNCTION(fmpz_poly_bit_pack),
    TEST_FUNCTION(fmpz_poly_bound_roots),
    TEST_FUNCTION(fmpz_poly_chebyshev_t),
    TEST_FUNCTION(fmpz_poly_chebyshev_u),
    TEST_FUNCTION(fmpz_poly_CLD_bound),
    TEST_FUNCTION(fmpz_poly_compose),
    TEST_FUNCTION(fmpz_poly_compose_divconquer),
    TEST_FUNCTION(fmpz_poly_compose_horner),
    TEST_FUNCTION(fmpz_poly_compose_series_brent_kung),
    TEST_FUNCTION(fmpz_poly_compose_series),
    TEST_FUNCTION(fmpz_poly_compose_series_horner),
    TEST_FUNCTION(fmpz_poly_content),
    TEST_FUNCTION(fmpz_poly_cos_minpoly),
    TEST_FUNCTION(fmpz_poly_CRT_ui),
    TEST_FUNCTION(fmpz_poly_CRT_ui_unsigned),
    TEST_FUNCTION(fmpz_poly_cyclotomic),
    TEST_FUNCTION(fmpz_poly_deflate),
    TEST_FUNCTION(fmpz_poly_derivative),
    TEST_FUNCTION(fmpz_poly_discriminant),
    TEST_FUNCTION(fmpz_poly_div_basecase),
    TEST_FUNCTION(fmpz_poly_div_divconquer),
    TEST_FUNCTION(fmpz_poly_divhigh_smodp),
    TEST_FUNCTION(fmpz_poly_divides),
    TEST_FUNCTION(fmpz_poly_divlow_smodp),
    TEST_FUNCTION(fmpz_poly_div_preinv),
    TEST_FUNCTION(fmpz_poly_divrem_basecase),
    TEST_FUNCTION(fmpz_poly_divrem),
    TEST_FUNCTION(fmpz_poly_divrem_divconquer),
    TEST_FUNCTION(fmpz_poly_divrem_preinv),
    TEST_FUNCTION(fmpz_poly_div_root),
    TEST_FUNCTION(fmpz_poly_div_series_basecase),
    TEST_FUNCTION(fmpz_poly_div_series),
    TEST_FUNCTION(fmpz_poly_div_series_divconquer),
    TEST_FUNCTION(fmpz_poly_equal_fmpz),
    TEST_FUNCTION(fmpz_poly_equal_trunc),
    TEST_FUNCTION(fmpz_poly_eta_qexp),
    TEST_FUNCTION(fmpz_poly_eulerian_polynomial),
    TEST_FUNCTION(fmpz_poly_evaluate_divconquer_fmpq),
    TEST_FUNCTION(fmpz_poly_evaluate_divconquer_fmpz),
    TEST_FUNCTION(fmpz_poly_evaluate_fmpq),
    TEST_FUNCTION(fmpz_poly_evaluate_fmpz),
    TEST_FUNCTION(fmpz_poly_evaluate_horner_d_2exp),
    TEST_FUNCTION(fmpz_poly_evaluate_horner_fmpq),
    TEST_FUNCTION(fmpz_poly_evaluate_horner_fmpz),
    TEST_FUNCTION(fmpz_poly_evaluate_mod),
    TEST_FUNCTION(fmpz_poly_fibonacci),
    TEST_FUNCTION(fmpz_poly_gcd),
    TEST_FUNCTION(fmpz_poly_gcd_heuristic),
    TEST_FUNCTION(fmpz_poly_gcd_modular),
    TEST_FUNCTION(fmpz_poly_gcd_subresultant),
    TEST_FUNCTION(fmpz_poly_get_coeff_ptr),
    TEST_FUNCTION(fmpz_poly_get_nmod_poly),
    TEST_FUNCTION(fmpz_poly_get_set_coeff_fmpz),
    TEST_FUNCTION(fmpz_poly_get_set_coeff_si),
    TEST_FUNCTION(fmpz_poly_get_set_coeff_ui),
    TEST_FUNCTION(fmpz_poly_get_set_str),
    TEST_FUNCTION(fmpz_poly_get_str),
    TEST_FUNCTION(fmpz_poly_get_str_pretty),
    TEST_FUNCTION(fmpz_poly_hensel_lift),
    TEST_FUNCTION(fmpz_poly_hensel_lift_once),
    TEST_FUNCTION(fmpz_poly_hensel_lift_without_only_inverse),
    TEST_FUNCTION(fmpz_poly_hensel_start_continue_lift),
    TEST_FUNCTION(fmpz_poly_hermite_h),
    TEST_FUNCTION(fmpz_poly_hermite_he),
    TEST_FUNCTION(fmpz_poly_inflate),
    TEST_FUNCTION(fmpz_poly_init_realloc_clear),
    TEST_FUNCTION(fmpz_poly_interpolate_fmpz_vec),
    TEST_FUNCTION(fmpz_poly_inv_series_basecase),
    TEST_FUNCTION(fmpz_poly_inv_series),
    TEST_FUNCTION(fmpz_poly_inv_series_newton),
    TEST_FUNCTION(fmpz_poly_is_cyclotomic),
    TEST_FUNCTION(fmpz_poly_is_squarefree),
    TEST_FUNCTION(fmpz_poly_lcm),
    TEST_FUNCTION(fmpz_poly_legendre_pt),
    TEST_FUNCTION(fmpz_poly_mul),
    TEST_FUNCTION(fmpz_poly_mul_classical),
    TEST_FUNCTION(fmpz_poly_mulhigh_classical),
    TEST_FUNCTION(fmpz_poly_mulhigh_karatsuba_n),
    TEST_FUNCTION(fmpz_poly_mulhigh_n),
    TEST_FUNCTION(fmpz_poly_mul_karatsuba),
    TEST_FUNCTION(fmpz_poly_mul_KS),
    TEST_FUNCTION(fmpz_poly_mullow),
    TEST_FUNCTION(fmpz_poly_mullow_classical),
    TEST_FUNCTION(fmpz_poly_mullow_karatsuba_n),
    TEST_FUNCTION(fmpz_poly_mullow_KS),
    TEST_FUNCTION(fmpz_poly_mullow_SS),
    TEST_FUNCTION(fmpz_poly_mullow_SS_precache),
    TEST_FUNCTION(fmpz_poly_mulmid_classical),
    TEST_FUNCTION(fmpz_poly_mul_SS),
    TEST_FUNCTION(fmpz_poly_mul_SS_precache),
    TEST_FUNCTION(fmpz_poly_neg),
    TEST_FUNCTION(fmpz_poly_newton_to_monomial),
    TEST_FUNCTION(fmpz_poly_nth_derivative),
    TEST_FUNCTION(fmpz_poly_num_real_roots),
    TEST_FUNCTION(fmpz_poly_num_real_roots_sturm),
    TEST_FUNCTION(fmpz_poly_pow_addchains),
    TEST_FUNCTION(fmpz_poly_pow_binexp),
    TEST_FUNCTION(fmpz_poly_pow_binomial),
    TEST_FUNCTION(fmpz_poly_pow),
    TEST_FUNCTION(fmpz_poly_power_sums),
    TEST_FUNCTION(fmpz_poly_pow_multinomial),
    TEST_FUNCTION(fmpz_poly_pow_trunc),
    TEST_FUNCTION(fmpz_poly_primitive_part),
    TEST_FUNCTION(fmpz_poly_print_read),
    TEST_FUNCTION(fmpz_poly_print_read_pretty),
    TEST_FUNCTION(fmpz_poly_product_roots_fmpq_vec),
    TEST_FUNCTION(fmpz_poly_product_roots_fmpz_vec),
    TEST_FUNCTION(fmpz_poly_pseudo_div),
    TEST_FUNCTION(fmpz_poly_pseudo_divrem_basecase),
    TEST_FUNCTION(fmpz_poly_pseudo_divrem_cohen),
    TEST_FUNCTION(fmpz_poly_pseudo_divrem_divconquer),
    TEST_FUNCTION(fmpz_poly_pseudo_rem),
    TEST_FUNCTION(fmpz_poly_pseudo_rem_cohen),
    TEST_FUNCTION(fmpz_poly_randtest_no_real_root),
    TEST_FUNCTION(fmpz_poly_rem_basecase),
    TEST_FUNCTION(fmpz_poly_remove),
    TEST_FUNCTION(fmpz_poly_remove_content_2exp),
    TEST_FUNCTION(fmpz_poly_rem_powers_precomp),
    TEST_FUNCTION(fmpz_poly_resultant),
    TEST_FUNCTION(fmpz_poly_resultant_euclidean),
    TEST_FUNCTION(fmpz_poly_resultant_modular),
    TEST_FUNCTION(fmpz_poly_resultant_modular_div),
    TEST_FUNCTION(fmpz_poly_reverse),
    TEST_FUNCTION(fmpz_poly_revert_series),
    TEST_FUNCTION(fmpz_poly_revert_series_lagrange),
    TEST_FUNCTION(fmpz_poly_revert_series_lagrange_fast),
    TEST_FUNCTION(fmpz_poly_revert_series_newton),
    TEST_FUNCTION(fmpz_poly_scalar_abs),
    TEST_FUNCTION(fmpz_poly_scalar_addmul_fmpz),
    TEST_FUNCTION(fmpz_poly_scalar_addmul_si),
    TEST_FUNCTION(fmpz_poly_scalar_addmul_ui),
    TEST_FUNCTION(fmpz_poly_scalar_mul_fmpz),
    TEST_FUNCTION(fmpz_poly_scalar_mul_si),
    TEST_FUNCTION(fmpz_poly_scalar_mul_ui),
    TEST_FUNCTION(fmpz_poly_scalar_submul_fmpz),
    TEST_FUNCTION(fmpz_poly_scale_2exp),
    TEST_FUNCTION(fmpz_poly_set_equal),
    TEST_FUNCTION(fmpz_poly_set_fmpz_equal),
    TEST_FUNCTION(fmpz_poly_set_si_equal),
    TEST_FUNCTION(fmpz_poly_set_trunc),
    TEST_FUNCTION(fmpz_poly_set_ui_equal),
    TEST_FUNCTION(fmpz_poly_shift_left_right),
    TEST_FUNCTION(fmpz_poly_signature),
    TEST_FUNCTION(fmpz_poly_sqr),
    TEST_FUNCTION(fmpz_poly_sqr_classical),
    TEST_FUNCTION(fmpz_poly_sqr_karatsuba),
    TEST_FUNCTION(fmpz_poly_sqr_KS),
    TEST_FUNCTION(fmpz_poly_sqrlow),
    TEST_FUNCTION(fmpz_poly_sqrlow_classical),
    TEST_FUNCTION(fmpz_poly_sqrlow_karatsuba_n),
    TEST_FUNCTION(fmpz_poly_sqrlow_KS),
    TEST_FUNCTION(fmpz_poly_sqrt),
    TEST_FUNCTION(fmpz_poly_sqrt_classical),
    TEST_FUNCTION(fmpz_poly_sqrt_divconquer),
    TEST_FUNCTION(fmpz_poly_sqrt_KS),
    TEST_FUNCTION(fmpz_poly_sqrtrem_classical),
    TEST_FUNCTION(fmpz_poly_sqrtrem_divconquer),
    TEST_FUNCTION(fmpz_poly_sqrt_series),
    TEST_FUNCTION(fmpz_poly_sub),
    TEST_FUNCTION(fmpz_poly_sub_series),
    TEST_FUNCTION(fmpz_poly_swap),
    TEST_FUNCTION(fmpz_poly_swinnerton_dyer),
    TEST_FUNCTION(fmpz_poly_taylor_shift),
    TEST_FUNCTION(fmpz_poly_taylor_shift_divconquer),
    TEST_FUNCTION(fmpz_poly_taylor_shift_horner),
    TEST_FUNCTION(fmpz_poly_taylor_shift_multi_mod_threaded),
    TEST_FUNCTION(fmpz_poly_theta_qexp),
    TEST_FUNCTION(fmpz_poly_xgcd_modular),
    TEST_FUNCTION(fmpz_poly_zero),
    TEST_FUNCTION(fmpz_poly_zero_coeffs)
};

char fmpz_poly_2norm_normalised_bits_name[] = "fmpz_poly_2norm_normalised_bits";
char fmpz_poly_add_name[] = "fmpz_poly_add";
char fmpz_poly_add_series_name[] = "fmpz_poly_add_series";
char fmpz_poly_add_sub_fmpz_name[] = "fmpz_poly_add_sub_fmpz";
char fmpz_poly_add_sub_si_name[] = "fmpz_poly_add_sub_si";
char fmpz_poly_bit_pack_name[] = "fmpz_poly_bit_pack";
char fmpz_poly_bound_roots_name[] = "fmpz_poly_bound_roots";
char fmpz_poly_chebyshev_t_name[] = "fmpz_poly_chebyshev_t";
char fmpz_poly_chebyshev_u_name[] = "fmpz_poly_chebyshev_u";
char fmpz_poly_CLD_bound_name[] = "fmpz_poly_CLD_bound";
char fmpz_poly_compose_name[] = "fmpz_poly_compose";
char fmpz_poly_compose_divconquer_name[] = "fmpz_poly_compose_divconquer";
char fmpz_poly_compose_horner_name[] = "fmpz_poly_compose_horner";
char fmpz_poly_compose_series_brent_kung_name[] = "fmpz_poly_compose_series_brent_kung";
char fmpz_poly_compose_series_name[] = "fmpz_poly_compose_series";
char fmpz_poly_compose_series_horner_name[] = "fmpz_poly_compose_series_horner";
char fmpz_poly_content_name[] = "fmpz_poly_content";
char fmpz_poly_cos_minpoly_name[] = "fmpz_poly_cos_minpoly";
char fmpz_poly_CRT_ui_name[] = "fmpz_poly_CRT_ui";
char fmpz_poly_CRT_ui_unsigned_name[] = "fmpz_poly_CRT_ui_unsigned";
char fmpz_poly_cyclotomic_name[] = "fmpz_poly_cyclotomic";
char fmpz_poly_deflate_name[] = "fmpz_poly_deflate";
char fmpz_poly_derivative_name[] = "fmpz_poly_derivative";
char fmpz_poly_discriminant_name[] = "fmpz_poly_discriminant";
char fmpz_poly_div_basecase_name[] = "fmpz_poly_div_basecase";
char fmpz_poly_div_divconquer_name[] = "fmpz_poly_div_divconquer";
char fmpz_poly_divhigh_smodp_name[] = "fmpz_poly_divhigh_smodp";
char fmpz_poly_divides_name[] = "fmpz_poly_divides";
char fmpz_poly_divlow_smodp_name[] = "fmpz_poly_divlow_smodp";
char fmpz_poly_div_preinv_name[] = "fmpz_poly_div_preinv";
char fmpz_poly_divrem_basecase_name[] = "fmpz_poly_divrem_basecase";
char fmpz_poly_divrem_name[] = "fmpz_poly_divrem";
char fmpz_poly_divrem_divconquer_name[] = "fmpz_poly_divrem_divconquer";
char fmpz_poly_divrem_preinv_name[] = "fmpz_poly_divrem_preinv";
char fmpz_poly_div_root_name[] = "fmpz_poly_div_root";
char fmpz_poly_div_series_basecase_name[] = "fmpz_poly_div_series_basecase";
char fmpz_poly_div_series_name[] = "fmpz_poly_div_series";
char fmpz_poly_div_series_divconquer_name[] = "fmpz_poly_div_series_divconquer";
char fmpz_poly_equal_fmpz_name[] = "fmpz_poly_equal_fmpz";
char fmpz_poly_equal_trunc_name[] = "fmpz_poly_equal_trunc";
char fmpz_poly_eta_qexp_name[] = "fmpz_poly_eta_qexp";
char fmpz_poly_eulerian_polynomial_name[] = "fmpz_poly_eulerian_polynomial";
char fmpz_poly_evaluate_divconquer_fmpq_name[] = "fmpz_poly_evaluate_divconquer_fmpq";
char fmpz_poly_evaluate_divconquer_fmpz_name[] = "fmpz_poly_evaluate_divconquer_fmpz";
char fmpz_poly_evaluate_fmpq_name[] = "fmpz_poly_evaluate_fmpq";
char fmpz_poly_evaluate_fmpz_name[] = "fmpz_poly_evaluate_fmpz";
char fmpz_poly_evaluate_horner_d_2exp_name[] = "fmpz_poly_evaluate_horner_d_2exp";
char fmpz_poly_evaluate_horner_fmpq_name[] = "fmpz_poly_evaluate_horner_fmpq";
char fmpz_poly_evaluate_horner_fmpz_name[] = "fmpz_poly_evaluate_horner_fmpz";
char fmpz_poly_evaluate_mod_name[] = "fmpz_poly_evaluate_mod";
char fmpz_poly_fibonacci_name[] = "fmpz_poly_fibonacci";
char fmpz_poly_gcd_name[] = "fmpz_poly_gcd";
char fmpz_poly_gcd_heuristic_name[] = "fmpz_poly_gcd_heuristic";
char fmpz_poly_gcd_modular_name[] = "fmpz_poly_gcd_modular";
char fmpz_poly_gcd_subresultant_name[] = "fmpz_poly_gcd_subresultant";
char fmpz_poly_get_coeff_ptr_name[] = "fmpz_poly_get_coeff_ptr";
char fmpz_poly_get_nmod_poly_name[] = "fmpz_poly_get_nmod_poly";
char fmpz_poly_get_set_coeff_fmpz_name[] = "fmpz_poly_get_set_coeff_fmpz";
char fmpz_poly_get_set_coeff_si_name[] = "fmpz_poly_get_set_coeff_si";
char fmpz_poly_get_set_coeff_ui_name[] = "fmpz_poly_get_set_coeff_ui";
char fmpz_poly_get_set_str_name[] = "fmpz_poly_get_set_str";
char fmpz_poly_get_str_name[] = "fmpz_poly_get_str";
char fmpz_poly_get_str_pretty_name[] = "fmpz_poly_get_str_pretty";
char fmpz_poly_hensel_lift_name[] = "fmpz_poly_hensel_lift";
char fmpz_poly_hensel_lift_once_name[] = "fmpz_poly_hensel_lift_once";
char fmpz_poly_hensel_lift_without_only_inverse_name[] = "fmpz_poly_hensel_lift_without_only_inverse";
char fmpz_poly_hensel_start_continue_lift_name[] = "fmpz_poly_hensel_start_continue_lift";
char fmpz_poly_hermite_h_name[] = "fmpz_poly_hermite_h";
char fmpz_poly_hermite_he_name[] = "fmpz_poly_hermite_he";
char fmpz_poly_inflate_name[] = "fmpz_poly_inflate";
char fmpz_poly_init_realloc_clear_name[] = "fmpz_poly_init_realloc_clear";
char fmpz_poly_interpolate_fmpz_vec_name[] = "fmpz_poly_interpolate_fmpz_vec";
char fmpz_poly_inv_series_basecase_name[] = "fmpz_poly_inv_series_basecase";
char fmpz_poly_inv_series_name[] = "fmpz_poly_inv_series";
char fmpz_poly_inv_series_newton_name[] = "fmpz_poly_inv_series_newton";
char fmpz_poly_is_cyclotomic_name[] = "fmpz_poly_is_cyclotomic";
char fmpz_poly_is_squarefree_name[] = "fmpz_poly_is_squarefree";
char fmpz_poly_lcm_name[] = "fmpz_poly_lcm";
char fmpz_poly_legendre_pt_name[] = "fmpz_poly_legendre_pt";
char fmpz_poly_mul_name[] = "fmpz_poly_mul";
char fmpz_poly_mul_classical_name[] = "fmpz_poly_mul_classical";
char fmpz_poly_mulhigh_classical_name[] = "fmpz_poly_mulhigh_classical";
char fmpz_poly_mulhigh_karatsuba_n_name[] = "fmpz_poly_mulhigh_karatsuba_n";
char fmpz_poly_mulhigh_n_name[] = "fmpz_poly_mulhigh_n";
char fmpz_poly_mul_karatsuba_name[] = "fmpz_poly_mul_karatsuba";
char fmpz_poly_mul_KS_name[] = "fmpz_poly_mul_KS";
char fmpz_poly_mullow_name[] = "fmpz_poly_mullow";
char fmpz_poly_mullow_classical_name[] = "fmpz_poly_mullow_classical";
char fmpz_poly_mullow_karatsuba_n_name[] = "fmpz_poly_mullow_karatsuba_n";
char fmpz_poly_mullow_KS_name[] = "fmpz_poly_mullow_KS";
char fmpz_poly_mullow_SS_name[] = "fmpz_poly_mullow_SS";
char fmpz_poly_mullow_SS_precache_name[] = "fmpz_poly_mullow_SS_precache";
char fmpz_poly_mulmid_classical_name[] = "fmpz_poly_mulmid_classical";
char fmpz_poly_mul_SS_name[] = "fmpz_poly_mul_SS";
char fmpz_poly_mul_SS_precache_name[] = "fmpz_poly_mul_SS_precache";
char fmpz_poly_neg_name[] = "fmpz_poly_neg";
char fmpz_poly_newton_to_monomial_name[] = "fmpz_poly_newton_to_monomial";
char fmpz_poly_nth_derivative_name[] = "fmpz_poly_nth_derivative";
char fmpz_poly_num_real_roots_name[] = "fmpz_poly_num_real_roots";
char fmpz_poly_num_real_roots_sturm_name[] = "fmpz_poly_num_real_roots_sturm";
char fmpz_poly_pow_addchains_name[] = "fmpz_poly_pow_addchains";
char fmpz_poly_pow_binexp_name[] = "fmpz_poly_pow_binexp";
char fmpz_poly_pow_binomial_name[] = "fmpz_poly_pow_binomial";
char fmpz_poly_pow_name[] = "fmpz_poly_pow";
char fmpz_poly_power_sums_name[] = "fmpz_poly_power_sums";
char fmpz_poly_pow_multinomial_name[] = "fmpz_poly_pow_multinomial";
char fmpz_poly_pow_trunc_name[] = "fmpz_poly_pow_trunc";
char fmpz_poly_primitive_part_name[] = "fmpz_poly_primitive_part";
char fmpz_poly_print_read_name[] = "fmpz_poly_print_read";
char fmpz_poly_print_read_pretty_name[] = "fmpz_poly_print_read_pretty";
char fmpz_poly_product_roots_fmpq_vec_name[] = "fmpz_poly_product_roots_fmpq_vec";
char fmpz_poly_product_roots_fmpz_vec_name[] = "fmpz_poly_product_roots_fmpz_vec";
char fmpz_poly_pseudo_div_name[] = "fmpz_poly_pseudo_div";
char fmpz_poly_pseudo_divrem_basecase_name[] = "fmpz_poly_pseudo_divrem_basecase";
char fmpz_poly_pseudo_divrem_cohen_name[] = "fmpz_poly_pseudo_divrem_cohen";
char fmpz_poly_pseudo_divrem_divconquer_name[] = "fmpz_poly_pseudo_divrem_divconquer";
char fmpz_poly_pseudo_rem_name[] = "fmpz_poly_pseudo_rem";
char fmpz_poly_pseudo_rem_cohen_name[] = "fmpz_poly_pseudo_rem_cohen";
char fmpz_poly_randtest_no_real_root_name[] = "fmpz_poly_randtest_no_real_root";
char fmpz_poly_rem_basecase_name[] = "fmpz_poly_rem_basecase";
char fmpz_poly_remove_name[] = "fmpz_poly_remove";
char fmpz_poly_remove_content_2exp_name[] = "fmpz_poly_remove_content_2exp";
char fmpz_poly_rem_powers_precomp_name[] = "fmpz_poly_rem_powers_precomp";
char fmpz_poly_resultant_name[] = "fmpz_poly_resultant";
char fmpz_poly_resultant_euclidean_name[] = "fmpz_poly_resultant_euclidean";
char fmpz_poly_resultant_modular_name[] = "fmpz_poly_resultant_modular";
char fmpz_poly_resultant_modular_div_name[] = "fmpz_poly_resultant_modular_div";
char fmpz_poly_reverse_name[] = "fmpz_poly_reverse";
char fmpz_poly_revert_series_name[] = "fmpz_poly_revert_series";
char fmpz_poly_revert_series_lagrange_name[] = "fmpz_poly_revert_series_lagrange";
char fmpz_poly_revert_series_lagrange_fast_name[] = "fmpz_poly_revert_series_lagrange_fast";
char fmpz_poly_revert_series_newton_name[] = "fmpz_poly_revert_series_newton";
char fmpz_poly_scalar_abs_name[] = "fmpz_poly_scalar_abs";
char fmpz_poly_scalar_addmul_fmpz_name[] = "fmpz_poly_scalar_addmul_fmpz";
char fmpz_poly_scalar_addmul_si_name[] = "fmpz_poly_scalar_addmul_si";
char fmpz_poly_scalar_addmul_ui_name[] = "fmpz_poly_scalar_addmul_ui";
char fmpz_poly_scalar_mul_fmpz_name[] = "fmpz_poly_scalar_mul_fmpz";
char fmpz_poly_scalar_mul_si_name[] = "fmpz_poly_scalar_mul_si";
char fmpz_poly_scalar_mul_ui_name[] = "fmpz_poly_scalar_mul_ui";
char fmpz_poly_scalar_submul_fmpz_name[] = "fmpz_poly_scalar_submul_fmpz";
char fmpz_poly_scale_2exp_name[] = "fmpz_poly_scale_2exp";
char fmpz_poly_set_equal_name[] = "fmpz_poly_set_equal";
char fmpz_poly_set_fmpz_equal_name[] = "fmpz_poly_set_fmpz_equal";
char fmpz_poly_set_si_equal_name[] = "fmpz_poly_set_si_equal";
char fmpz_poly_set_trunc_name[] = "fmpz_poly_set_trunc";
char fmpz_poly_set_ui_equal_name[] = "fmpz_poly_set_ui_equal";
char fmpz_poly_shift_left_right_name[] = "fmpz_poly_shift_left_right";
char fmpz_poly_signature_name[] = "fmpz_poly_signature";
char fmpz_poly_sqr_name[] = "fmpz_poly_sqr";
char fmpz_poly_sqr_classical_name[] = "fmpz_poly_sqr_classical";
char fmpz_poly_sqr_karatsuba_name[] = "fmpz_poly_sqr_karatsuba";
char fmpz_poly_sqr_KS_name[] = "fmpz_poly_sqr_KS";
char fmpz_poly_sqrlow_name[] = "fmpz_poly_sqrlow";
char fmpz_poly_sqrlow_classical_name[] = "fmpz_poly_sqrlow_classical";
char fmpz_poly_sqrlow_karatsuba_n_name[] = "fmpz_poly_sqrlow_karatsuba_n";
char fmpz_poly_sqrlow_KS_name[] = "fmpz_poly_sqrlow_KS";
char fmpz_poly_sqrt_name[] = "fmpz_poly_sqrt";
char fmpz_poly_sqrt_classical_name[] = "fmpz_poly_sqrt_classical";
char fmpz_poly_sqrt_divconquer_name[] = "fmpz_poly_sqrt_divconquer";
char fmpz_poly_sqrt_KS_name[] = "fmpz_poly_sqrt_KS";
char fmpz_poly_sqrtrem_classical_name[] = "fmpz_poly_sqrtrem_classical";
char fmpz_poly_sqrtrem_divconquer_name[] = "fmpz_poly_sqrtrem_divconquer";
char fmpz_poly_sqrt_series_name[] = "fmpz_poly_sqrt_series";
char fmpz_poly_sub_name[] = "fmpz_poly_sub";
char fmpz_poly_sub_series_name[] = "fmpz_poly_sub_series";
char fmpz_poly_swap_name[] = "fmpz_poly_swap";
char fmpz_poly_swinnerton_dyer_name[] = "fmpz_poly_swinnerton_dyer";
char fmpz_poly_taylor_shift_name[] = "fmpz_poly_taylor_shift";
char fmpz_poly_taylor_shift_divconquer_name[] = "fmpz_poly_taylor_shift_divconquer";
char fmpz_poly_taylor_shift_horner_name[] = "fmpz_poly_taylor_shift_horner";
char fmpz_poly_taylor_shift_multi_mod_threaded_name[] = "fmpz_poly_taylor_shift_multi_mod_threaded";
char fmpz_poly_theta_qexp_name[] = "fmpz_poly_theta_qexp";
char fmpz_poly_xgcd_modular_name[] = "fmpz_poly_xgcd_modular";
char fmpz_poly_zero_name[] = "fmpz_poly_zero";
char fmpz_poly_zero_coeffs_name[] = "fmpz_poly_zero_coeffs";

char * test_names[] =
{
    fmpz_poly_2norm_normalised_bits_name,
    fmpz_poly_add_name,
    fmpz_poly_add_series_name,
    fmpz_poly_add_sub_fmpz_name,
    fmpz_poly_add_sub_si_name,
    fmpz_poly_bit_pack_name,
    fmpz_poly_bound_roots_name,
    fmpz_poly_chebyshev_t_name,
    fmpz_poly_chebyshev_u_name,
    fmpz_poly_CLD_bound_name,
    fmpz_poly_compose_name,
    fmpz_poly_compose_divconquer_name,
    fmpz_poly_compose_horner_name,
    fmpz_poly_compose_series_brent_kung_name,
    fmpz_poly_compose_series_name,
    fmpz_poly_compose_series_horner_name,
    fmpz_poly_content_name,
    fmpz_poly_cos_minpoly_name,
    fmpz_poly_CRT_ui_name,
    fmpz_poly_CRT_ui_unsigned_name,
    fmpz_poly_cyclotomic_name,
    fmpz_poly_deflate_name,
    fmpz_poly_derivative_name,
    fmpz_poly_discriminant_name,
    fmpz_poly_div_basecase_name,
    fmpz_poly_div_divconquer_name,
    fmpz_poly_divhigh_smodp_name,
    fmpz_poly_divides_name,
    fmpz_poly_divlow_smodp_name,
    fmpz_poly_div_preinv_name,
    fmpz_poly_divrem_basecase_name,
    fmpz_poly_divrem_name,
    fmpz_poly_divrem_divconquer_name,
    fmpz_poly_divrem_preinv_name,
    fmpz_poly_div_root_name,
    fmpz_poly_div_series_basecase_name,
    fmpz_poly_div_series_name,
    fmpz_poly_div_series_divconquer_name,
    fmpz_poly_equal_fmpz_name,
    fmpz_poly_equal_trunc_name,
    fmpz_poly_eta_qexp_name,
    fmpz_poly_eulerian_polynomial_name,
    fmpz_poly_evaluate_divconquer_fmpq_name,
    fmpz_poly_evaluate_divconquer_fmpz_name,
    fmpz_poly_evaluate_fmpq_name,
    fmpz_poly_evaluate_fmpz_name,
    fmpz_poly_evaluate_horner_d_2exp_name,
    fmpz_poly_evaluate_horner_fmpq_name,
    fmpz_poly_evaluate_horner_fmpz_name,
    fmpz_poly_evaluate_mod_name,
    fmpz_poly_fibonacci_name,
    fmpz_poly_gcd_name,
    fmpz_poly_gcd_heuristic_name,
    fmpz_poly_gcd_modular_name,
    fmpz_poly_gcd_subresultant_name,
    fmpz_poly_get_coeff_ptr_name,
    fmpz_poly_get_nmod_poly_name,
    fmpz_poly_get_set_coeff_fmpz_name,
    fmpz_poly_get_set_coeff_si_name,
    fmpz_poly_get_set_coeff_ui_name,
    fmpz_poly_get_set_str_name,
    fmpz_poly_get_str_name,
    fmpz_poly_get_str_pretty_name,
    fmpz_poly_hensel_lift_name,
    fmpz_poly_hensel_lift_once_name,
    fmpz_poly_hensel_lift_without_only_inverse_name,
    fmpz_poly_hensel_start_continue_lift_name,
    fmpz_poly_hermite_h_name,
    fmpz_poly_hermite_he_name,
    fmpz_poly_inflate_name,
    fmpz_poly_init_realloc_clear_name,
    fmpz_poly_interpolate_fmpz_vec_name,
    fmpz_poly_inv_series_basecase_name,
    fmpz_poly_inv_series_name,
    fmpz_poly_inv_series_newton_name,
    fmpz_poly_is_cyclotomic_name,
    fmpz_poly_is_squarefree_name,
    fmpz_poly_lcm_name,
    fmpz_poly_legendre_pt_name,
    fmpz_poly_mul_name,
    fmpz_poly_mul_classical_name,
    fmpz_poly_mulhigh_classical_name,
    fmpz_poly_mulhigh_karatsuba_n_name,
    fmpz_poly_mulhigh_n_name,
    fmpz_poly_mul_karatsuba_name,
    fmpz_poly_mul_KS_name,
    fmpz_poly_mullow_name,
    fmpz_poly_mullow_classical_name,
    fmpz_poly_mullow_karatsuba_n_name,
    fmpz_poly_mullow_KS_name,
    fmpz_poly_mullow_SS_name,
    fmpz_poly_mullow_SS_precache_name,
    fmpz_poly_mulmid_classical_name,
    fmpz_poly_mul_SS_name,
    fmpz_poly_mul_SS_precache_name,
    fmpz_poly_neg_name,
    fmpz_poly_newton_to_monomial_name,
    fmpz_poly_nth_derivative_name,
    fmpz_poly_num_real_roots_name,
    fmpz_poly_num_real_roots_sturm_name,
    fmpz_poly_pow_addchains_name,
    fmpz_poly_pow_binexp_name,
    fmpz_poly_pow_binomial_name,
    fmpz_poly_pow_name,
    fmpz_poly_power_sums_name,
    fmpz_poly_pow_multinomial_name,
    fmpz_poly_pow_trunc_name,
    fmpz_poly_primitive_part_name,
    fmpz_poly_print_read_name,
    fmpz_poly_print_read_pretty_name,
    fmpz_poly_product_roots_fmpq_vec_name,
    fmpz_poly_product_roots_fmpz_vec_name,
    fmpz_poly_pseudo_div_name,
    fmpz_poly_pseudo_divrem_basecase_name,
    fmpz_poly_pseudo_divrem_cohen_name,
    fmpz_poly_pseudo_divrem_divconquer_name,
    fmpz_poly_pseudo_rem_name,
    fmpz_poly_pseudo_rem_cohen_name,
    fmpz_poly_randtest_no_real_root_name,
    fmpz_poly_rem_basecase_name,
    fmpz_poly_remove_name,
    fmpz_poly_remove_content_2exp_name,
    fmpz_poly_rem_powers_precomp_name,
    fmpz_poly_resultant_name,
    fmpz_poly_resultant_euclidean_name,
    fmpz_poly_resultant_modular_name,
    fmpz_poly_resultant_modular_div_name,
    fmpz_poly_reverse_name,
    fmpz_poly_revert_series_name,
    fmpz_poly_revert_series_lagrange_name,
    fmpz_poly_revert_series_lagrange_fast_name,
    fmpz_poly_revert_series_newton_name,
    fmpz_poly_scalar_abs_name,
    fmpz_poly_scalar_addmul_fmpz_name,
    fmpz_poly_scalar_addmul_si_name,
    fmpz_poly_scalar_addmul_ui_name,
    fmpz_poly_scalar_mul_fmpz_name,
    fmpz_poly_scalar_mul_si_name,
    fmpz_poly_scalar_mul_ui_name,
    fmpz_poly_scalar_submul_fmpz_name,
    fmpz_poly_scale_2exp_name,
    fmpz_poly_set_equal_name,
    fmpz_poly_set_fmpz_equal_name,
    fmpz_poly_set_si_equal_name,
    fmpz_poly_set_trunc_name,
    fmpz_poly_set_ui_equal_name,
    fmpz_poly_shift_left_right_name,
    fmpz_poly_signature_name,
    fmpz_poly_sqr_name,
    fmpz_poly_sqr_classical_name,
    fmpz_poly_sqr_karatsuba_name,
    fmpz_poly_sqr_KS_name,
    fmpz_poly_sqrlow_name,
    fmpz_poly_sqrlow_classical_name,
    fmpz_poly_sqrlow_karatsuba_n_name,
    fmpz_poly_sqrlow_KS_name,
    fmpz_poly_sqrt_name,
    fmpz_poly_sqrt_classical_name,
    fmpz_poly_sqrt_divconquer_name,
    fmpz_poly_sqrt_KS_name,
    fmpz_poly_sqrtrem_classical_name,
    fmpz_poly_sqrtrem_divconquer_name,
    fmpz_poly_sqrt_series_name,
    fmpz_poly_sub_name,
    fmpz_poly_sub_series_name,
    fmpz_poly_swap_name,
    fmpz_poly_swinnerton_dyer_name,
    fmpz_poly_taylor_shift_name,
    fmpz_poly_taylor_shift_divconquer_name,
    fmpz_poly_taylor_shift_horner_name,
    fmpz_poly_taylor_shift_multi_mod_threaded_name,
    fmpz_poly_theta_qexp_name,
    fmpz_poly_xgcd_modular_name,
    fmpz_poly_zero_name,
    fmpz_poly_zero_coeffs_name
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
