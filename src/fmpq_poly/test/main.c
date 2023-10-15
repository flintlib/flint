/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Try to get fdopen declared for fmpq_poly_[print/read] */
#if defined __STRICT_ANSI__
# undef __STRICT_ANSI__
#endif

#include <string.h>

/* Include functions *********************************************************/

#include "t-add.c"
#include "t-add_series.c"
#include "t-add_sub_fmpq.c"
#include "t-add_sub_fmpz.c"
#include "t-add_sub_si.c"
#include "t-asinh_series.c"
#include "t-asin_series.c"
#include "t-atanh_series.c"
#include "t-atan_series.c"
#include "t-cmp.c"
#include "t-compose.c"
#include "t-compose_series_brent_kung.c"
#include "t-compose_series.c"
#include "t-compose_series_horner.c"
#include "t-content.c"
#include "t-cosh_series.c"
#include "t-cos_series.c"
#include "t-derivative.c"
#include "t-div.c"
#include "t-divides.c"
#include "t-divrem.c"
#include "t-div_series.c"
#include "t-equal_trunc.c"
#include "t-evaluate_fmpq.c"
#include "t-evaluate_fmpz.c"
#include "t-exp_expinv_series.c"
#include "t-exp_series.c"
#include "t-gcd.c"
#include "t-gegenbauer_c.c"
#include "t-get_nmod_poly.c"
#include "t-get_set_coeff_fmpq.c"
#include "t-get_set_coeff_fmpz.c"
#include "t-get_set_coeff_si.c"
#include "t-get_set_coeff_ui.c"
#include "t-get_set_str.c"
#include "t-get_slice.c"
#include "t-init_realloc_clear.c"
#include "t-integral.c"
#include "t-interpolate_fmpz_vec.c"
#include "t-inv.c"
#include "t-inv_series_newton.c"
#include "t-invsqrt_series.c"
#include "t-is_squarefree.c"
#include "t-laguerre_l.c"
#include "t-lcm.c"
#include "t-legendre_p.c"
#include "t-log_series.c"
#include "t-make_monic.c"
#include "t-mul.c"
#include "t-mullow.c"
#include "t-neg.c"
#include "t-nth_derivative.c"
#include "t-pow.c"
#include "t-power_sums.c"
#include "t-pow_trunc.c"
#include "t-primitive_part.c"
#include "t-print_read.c"
#include "t-rem.c"
#include "t-remove.c"
#include "t-rem_powers_precomp.c"
#include "t-rescale.c"
#include "t-resultant.c"
#include "t-resultant_div.c"
#include "t-reverse.c"
#include "t-revert_series.c"
#include "t-revert_series_lagrange.c"
#include "t-revert_series_lagrange_fast.c"
#include "t-revert_series_newton.c"
#include "t-scalar_div_fmpq.c"
#include "t-scalar_div_fmpz.c"
#include "t-scalar_div_si.c"
#include "t-scalar_div_ui.c"
#include "t-scalar_mul_fmpq.c"
#include "t-scalar_mul_fmpz.c"
#include "t-scalar_mul_si.c"
#include "t-scalar_mul_ui.c"
#include "t-set_equal.c"
#include "t-set_trunc.c"
#include "t-shift_left_right.c"
#include "t-sin_cos_series.c"
#include "t-sinh_cosh_series.c"
#include "t-sinh_series.c"
#include "t-sin_series.c"
#include "t-sqrt_series.c"
#include "t-sub.c"
#include "t-sub_series.c"
#include "t-swap.c"
#include "t-tanh_series.c"
#include "t-tan_series.c"
#include "t-xgcd.c"
#include "t-zero.c"

/* Array of test functions ***************************************************/

int (*test_functions[])(void) =
{
    TEST_FUNCTION(fmpq_poly_add),
    TEST_FUNCTION(fmpq_poly_add_series),
    TEST_FUNCTION(fmpq_poly_add_sub_fmpq),
    TEST_FUNCTION(fmpq_poly_add_sub_fmpz),
    TEST_FUNCTION(fmpq_poly_add_sub_si),
    TEST_FUNCTION(fmpq_poly_asinh_series),
    TEST_FUNCTION(fmpq_poly_asin_series),
    TEST_FUNCTION(fmpq_poly_atanh_series),
    TEST_FUNCTION(fmpq_poly_atan_series),
    TEST_FUNCTION(fmpq_poly_cmp),
    TEST_FUNCTION(fmpq_poly_compose),
    TEST_FUNCTION(fmpq_poly_compose_series_brent_kung),
    TEST_FUNCTION(fmpq_poly_compose_series),
    TEST_FUNCTION(fmpq_poly_compose_series_horner),
    TEST_FUNCTION(fmpq_poly_content),
    TEST_FUNCTION(fmpq_poly_cosh_series),
    TEST_FUNCTION(fmpq_poly_cos_series),
    TEST_FUNCTION(fmpq_poly_derivative),
    TEST_FUNCTION(fmpq_poly_div),
    TEST_FUNCTION(fmpq_poly_divides),
    TEST_FUNCTION(fmpq_poly_divrem),
    TEST_FUNCTION(fmpq_poly_div_series),
    TEST_FUNCTION(fmpq_poly_equal_trunc),
    TEST_FUNCTION(fmpq_poly_evaluate_fmpq),
    TEST_FUNCTION(fmpq_poly_evaluate_fmpz),
    TEST_FUNCTION(fmpq_poly_exp_expinv_series),
    TEST_FUNCTION(fmpq_poly_exp_series),
    TEST_FUNCTION(fmpq_poly_gcd),
    TEST_FUNCTION(fmpq_poly_gegenbauer_c),
    TEST_FUNCTION(fmpq_poly_get_nmod_poly),
    TEST_FUNCTION(fmpq_poly_get_set_coeff_fmpq),
    TEST_FUNCTION(fmpq_poly_get_set_coeff_fmpz),
    TEST_FUNCTION(fmpq_poly_get_set_coeff_si),
    TEST_FUNCTION(fmpq_poly_get_set_coeff_ui),
    TEST_FUNCTION(fmpq_poly_get_set_str),
    TEST_FUNCTION(fmpq_poly_get_slice),
    TEST_FUNCTION(fmpq_poly_init_realloc_clear),
    TEST_FUNCTION(fmpq_poly_integral),
    TEST_FUNCTION(fmpq_poly_interpolate_fmpz_vec),
    TEST_FUNCTION(fmpq_poly_inv),
    TEST_FUNCTION(fmpq_poly_inv_series_newton),
    TEST_FUNCTION(fmpq_poly_invsqrt_series),
    TEST_FUNCTION(fmpq_poly_is_squarefree),
    TEST_FUNCTION(fmpq_poly_laguerre_l),
    TEST_FUNCTION(fmpq_poly_lcm),
    TEST_FUNCTION(fmpq_poly_legendre_p),
    TEST_FUNCTION(fmpq_poly_log_series),
    TEST_FUNCTION(fmpq_poly_make_monic),
    TEST_FUNCTION(fmpq_poly_mul),
    TEST_FUNCTION(fmpq_poly_mullow),
    TEST_FUNCTION(fmpq_poly_neg),
    TEST_FUNCTION(fmpq_poly_nth_derivative),
    TEST_FUNCTION(fmpq_poly_pow),
    TEST_FUNCTION(fmpq_poly_power_sums),
    TEST_FUNCTION(fmpq_poly_pow_trunc),
    TEST_FUNCTION(fmpq_poly_primitive_part),
    TEST_FUNCTION(fmpq_poly_print_read),
    TEST_FUNCTION(fmpq_poly_rem),
    TEST_FUNCTION(fmpq_poly_remove),
    TEST_FUNCTION(fmpq_poly_rem_powers_precomp),
    TEST_FUNCTION(fmpq_poly_rescale),
    TEST_FUNCTION(fmpq_poly_resultant),
    TEST_FUNCTION(fmpq_poly_resultant_div),
    TEST_FUNCTION(fmpq_poly_reverse),
    TEST_FUNCTION(fmpq_poly_revert_series),
    TEST_FUNCTION(fmpq_poly_revert_series_lagrange),
    TEST_FUNCTION(fmpq_poly_revert_series_lagrange_fast),
    TEST_FUNCTION(fmpq_poly_revert_series_newton),
    TEST_FUNCTION(fmpq_poly_scalar_div_fmpq),
    TEST_FUNCTION(fmpq_poly_scalar_div_fmpz),
    TEST_FUNCTION(fmpq_poly_scalar_div_si),
    TEST_FUNCTION(fmpq_poly_scalar_div_ui),
    TEST_FUNCTION(fmpq_poly_scalar_mul_fmpq),
    TEST_FUNCTION(fmpq_poly_scalar_mul_fmpz),
    TEST_FUNCTION(fmpq_poly_scalar_mul_si),
    TEST_FUNCTION(fmpq_poly_scalar_mul_ui),
    TEST_FUNCTION(fmpq_poly_set_equal),
    TEST_FUNCTION(fmpq_poly_set_trunc),
    TEST_FUNCTION(fmpq_poly_shift_left_right),
    TEST_FUNCTION(fmpq_poly_sin_cos_series),
    TEST_FUNCTION(fmpq_poly_sinh_cosh_series),
    TEST_FUNCTION(fmpq_poly_sinh_series),
    TEST_FUNCTION(fmpq_poly_sin_series),
    TEST_FUNCTION(fmpq_poly_sqrt_series),
    TEST_FUNCTION(fmpq_poly_sub),
    TEST_FUNCTION(fmpq_poly_sub_series),
    TEST_FUNCTION(fmpq_poly_swap),
    TEST_FUNCTION(fmpq_poly_tanh_series),
    TEST_FUNCTION(fmpq_poly_tan_series),
    TEST_FUNCTION(fmpq_poly_xgcd),
    TEST_FUNCTION(fmpq_poly_zero)
};

char fmpq_poly_add_name[] = "fmpq_poly_add";
char fmpq_poly_add_series_name[] = "fmpq_poly_add_series";
char fmpq_poly_add_sub_fmpq_name[] = "fmpq_poly_add_sub_fmpq";
char fmpq_poly_add_sub_fmpz_name[] = "fmpq_poly_add_sub_fmpz";
char fmpq_poly_add_sub_si_name[] = "fmpq_poly_add_sub_si";
char fmpq_poly_asinh_series_name[] = "fmpq_poly_asinh_series";
char fmpq_poly_asin_series_name[] = "fmpq_poly_asin_series";
char fmpq_poly_atanh_series_name[] = "fmpq_poly_atanh_series";
char fmpq_poly_atan_series_name[] = "fmpq_poly_atan_series";
char fmpq_poly_cmp_name[] = "fmpq_poly_cmp";
char fmpq_poly_compose_name[] = "fmpq_poly_compose";
char fmpq_poly_compose_series_brent_kung_name[] = "fmpq_poly_compose_series_brent_kung";
char fmpq_poly_compose_series_name[] = "fmpq_poly_compose_series";
char fmpq_poly_compose_series_horner_name[] = "fmpq_poly_compose_series_horner";
char fmpq_poly_content_name[] = "fmpq_poly_content";
char fmpq_poly_cosh_series_name[] = "fmpq_poly_cosh_series";
char fmpq_poly_cos_series_name[] = "fmpq_poly_cos_series";
char fmpq_poly_derivative_name[] = "fmpq_poly_derivative";
char fmpq_poly_div_name[] = "fmpq_poly_div";
char fmpq_poly_divides_name[] = "fmpq_poly_divides";
char fmpq_poly_divrem_name[] = "fmpq_poly_divrem";
char fmpq_poly_div_series_name[] = "fmpq_poly_div_series";
char fmpq_poly_equal_trunc_name[] = "fmpq_poly_equal_trunc";
char fmpq_poly_evaluate_fmpq_name[] = "fmpq_poly_evaluate_fmpq";
char fmpq_poly_evaluate_fmpz_name[] = "fmpq_poly_evaluate_fmpz";
char fmpq_poly_exp_expinv_series_name[] = "fmpq_poly_exp_expinv_series";
char fmpq_poly_exp_series_name[] = "fmpq_poly_exp_series";
char fmpq_poly_gcd_name[] = "fmpq_poly_gcd";
char fmpq_poly_gegenbauer_c_name[] = "fmpq_poly_gegenbauer_c";
char fmpq_poly_get_nmod_poly_name[] = "fmpq_poly_get_nmod_poly";
char fmpq_poly_get_set_coeff_fmpq_name[] = "fmpq_poly_get_set_coeff_fmpq";
char fmpq_poly_get_set_coeff_fmpz_name[] = "fmpq_poly_get_set_coeff_fmpz";
char fmpq_poly_get_set_coeff_si_name[] = "fmpq_poly_get_set_coeff_si";
char fmpq_poly_get_set_coeff_ui_name[] = "fmpq_poly_get_set_coeff_ui";
char fmpq_poly_get_set_str_name[] = "fmpq_poly_get_set_str";
char fmpq_poly_get_slice_name[] = "fmpq_poly_get_slice";
char fmpq_poly_init_realloc_clear_name[] = "fmpq_poly_init_realloc_clear";
char fmpq_poly_integral_name[] = "fmpq_poly_integral";
char fmpq_poly_interpolate_fmpz_vec_name[] = "fmpq_poly_interpolate_fmpz_vec";
char fmpq_poly_inv_name[] = "fmpq_poly_inv";
char fmpq_poly_inv_series_newton_name[] = "fmpq_poly_inv_series_newton";
char fmpq_poly_invsqrt_series_name[] = "fmpq_poly_invsqrt_series";
char fmpq_poly_is_squarefree_name[] = "fmpq_poly_is_squarefree";
char fmpq_poly_laguerre_l_name[] = "fmpq_poly_laguerre_l";
char fmpq_poly_lcm_name[] = "fmpq_poly_lcm";
char fmpq_poly_legendre_p_name[] = "fmpq_poly_legendre_p";
char fmpq_poly_log_series_name[] = "fmpq_poly_log_series";
char fmpq_poly_make_monic_name[] = "fmpq_poly_make_monic";
char fmpq_poly_mul_name[] = "fmpq_poly_mul";
char fmpq_poly_mullow_name[] = "fmpq_poly_mullow";
char fmpq_poly_neg_name[] = "fmpq_poly_neg";
char fmpq_poly_nth_derivative_name[] = "fmpq_poly_nth_derivative";
char fmpq_poly_pow_name[] = "fmpq_poly_pow";
char fmpq_poly_power_sums_name[] = "fmpq_poly_power_sums";
char fmpq_poly_pow_trunc_name[] = "fmpq_poly_pow_trunc";
char fmpq_poly_primitive_part_name[] = "fmpq_poly_primitive_part";
char fmpq_poly_print_read_name[] = "fmpq_poly_print_read";
char fmpq_poly_rem_name[] = "fmpq_poly_rem";
char fmpq_poly_remove_name[] = "fmpq_poly_remove";
char fmpq_poly_rem_powers_precomp_name[] = "fmpq_poly_rem_powers_precomp";
char fmpq_poly_rescale_name[] = "fmpq_poly_rescale";
char fmpq_poly_resultant_name[] = "fmpq_poly_resultant";
char fmpq_poly_resultant_div_name[] = "fmpq_poly_resultant_div";
char fmpq_poly_reverse_name[] = "fmpq_poly_reverse";
char fmpq_poly_revert_series_name[] = "fmpq_poly_revert_series";
char fmpq_poly_revert_series_lagrange_name[] = "fmpq_poly_revert_series_lagrange";
char fmpq_poly_revert_series_lagrange_fast_name[] = "fmpq_poly_revert_series_lagrange_fast";
char fmpq_poly_revert_series_newton_name[] = "fmpq_poly_revert_series_newton";
char fmpq_poly_scalar_div_fmpq_name[] = "fmpq_poly_scalar_div_fmpq";
char fmpq_poly_scalar_div_fmpz_name[] = "fmpq_poly_scalar_div_fmpz";
char fmpq_poly_scalar_div_si_name[] = "fmpq_poly_scalar_div_si";
char fmpq_poly_scalar_div_ui_name[] = "fmpq_poly_scalar_div_ui";
char fmpq_poly_scalar_mul_fmpq_name[] = "fmpq_poly_scalar_mul_fmpq";
char fmpq_poly_scalar_mul_fmpz_name[] = "fmpq_poly_scalar_mul_fmpz";
char fmpq_poly_scalar_mul_si_name[] = "fmpq_poly_scalar_mul_si";
char fmpq_poly_scalar_mul_ui_name[] = "fmpq_poly_scalar_mul_ui";
char fmpq_poly_set_equal_name[] = "fmpq_poly_set_equal";
char fmpq_poly_set_trunc_name[] = "fmpq_poly_set_trunc";
char fmpq_poly_shift_left_right_name[] = "fmpq_poly_shift_left_right";
char fmpq_poly_sin_cos_series_name[] = "fmpq_poly_sin_cos_series";
char fmpq_poly_sinh_cosh_series_name[] = "fmpq_poly_sinh_cosh_series";
char fmpq_poly_sinh_series_name[] = "fmpq_poly_sinh_series";
char fmpq_poly_sin_series_name[] = "fmpq_poly_sin_series";
char fmpq_poly_sqrt_series_name[] = "fmpq_poly_sqrt_series";
char fmpq_poly_sub_name[] = "fmpq_poly_sub";
char fmpq_poly_sub_series_name[] = "fmpq_poly_sub_series";
char fmpq_poly_swap_name[] = "fmpq_poly_swap";
char fmpq_poly_tanh_series_name[] = "fmpq_poly_tanh_series";
char fmpq_poly_tan_series_name[] = "fmpq_poly_tan_series";
char fmpq_poly_xgcd_name[] = "fmpq_poly_xgcd";
char fmpq_poly_zero_name[] = "fmpq_poly_zero";

char * test_names[] =
{
    fmpq_poly_add_name,
    fmpq_poly_add_series_name,
    fmpq_poly_add_sub_fmpq_name,
    fmpq_poly_add_sub_fmpz_name,
    fmpq_poly_add_sub_si_name,
    fmpq_poly_asinh_series_name,
    fmpq_poly_asin_series_name,
    fmpq_poly_atanh_series_name,
    fmpq_poly_atan_series_name,
    fmpq_poly_cmp_name,
    fmpq_poly_compose_name,
    fmpq_poly_compose_series_brent_kung_name,
    fmpq_poly_compose_series_name,
    fmpq_poly_compose_series_horner_name,
    fmpq_poly_content_name,
    fmpq_poly_cosh_series_name,
    fmpq_poly_cos_series_name,
    fmpq_poly_derivative_name,
    fmpq_poly_div_name,
    fmpq_poly_divides_name,
    fmpq_poly_divrem_name,
    fmpq_poly_div_series_name,
    fmpq_poly_equal_trunc_name,
    fmpq_poly_evaluate_fmpq_name,
    fmpq_poly_evaluate_fmpz_name,
    fmpq_poly_exp_expinv_series_name,
    fmpq_poly_exp_series_name,
    fmpq_poly_gcd_name,
    fmpq_poly_gegenbauer_c_name,
    fmpq_poly_get_nmod_poly_name,
    fmpq_poly_get_set_coeff_fmpq_name,
    fmpq_poly_get_set_coeff_fmpz_name,
    fmpq_poly_get_set_coeff_si_name,
    fmpq_poly_get_set_coeff_ui_name,
    fmpq_poly_get_set_str_name,
    fmpq_poly_get_slice_name,
    fmpq_poly_init_realloc_clear_name,
    fmpq_poly_integral_name,
    fmpq_poly_interpolate_fmpz_vec_name,
    fmpq_poly_inv_name,
    fmpq_poly_inv_series_newton_name,
    fmpq_poly_invsqrt_series_name,
    fmpq_poly_is_squarefree_name,
    fmpq_poly_laguerre_l_name,
    fmpq_poly_lcm_name,
    fmpq_poly_legendre_p_name,
    fmpq_poly_log_series_name,
    fmpq_poly_make_monic_name,
    fmpq_poly_mul_name,
    fmpq_poly_mullow_name,
    fmpq_poly_neg_name,
    fmpq_poly_nth_derivative_name,
    fmpq_poly_pow_name,
    fmpq_poly_power_sums_name,
    fmpq_poly_pow_trunc_name,
    fmpq_poly_primitive_part_name,
    fmpq_poly_print_read_name,
    fmpq_poly_rem_name,
    fmpq_poly_remove_name,
    fmpq_poly_rem_powers_precomp_name,
    fmpq_poly_rescale_name,
    fmpq_poly_resultant_name,
    fmpq_poly_resultant_div_name,
    fmpq_poly_reverse_name,
    fmpq_poly_revert_series_name,
    fmpq_poly_revert_series_lagrange_name,
    fmpq_poly_revert_series_lagrange_fast_name,
    fmpq_poly_revert_series_newton_name,
    fmpq_poly_scalar_div_fmpq_name,
    fmpq_poly_scalar_div_fmpz_name,
    fmpq_poly_scalar_div_si_name,
    fmpq_poly_scalar_div_ui_name,
    fmpq_poly_scalar_mul_fmpq_name,
    fmpq_poly_scalar_mul_fmpz_name,
    fmpq_poly_scalar_mul_si_name,
    fmpq_poly_scalar_mul_ui_name,
    fmpq_poly_set_equal_name,
    fmpq_poly_set_trunc_name,
    fmpq_poly_shift_left_right_name,
    fmpq_poly_sin_cos_series_name,
    fmpq_poly_sinh_cosh_series_name,
    fmpq_poly_sinh_series_name,
    fmpq_poly_sin_series_name,
    fmpq_poly_sqrt_series_name,
    fmpq_poly_sub_name,
    fmpq_poly_sub_series_name,
    fmpq_poly_swap_name,
    fmpq_poly_tanh_series_name,
    fmpq_poly_tan_series_name,
    fmpq_poly_xgcd_name,
    fmpq_poly_zero_name
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
