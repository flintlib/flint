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
#include <stdlib.h>
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

test_struct tests[] =
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

/* main function *************************************************************/

TEST_MAIN(tests)
