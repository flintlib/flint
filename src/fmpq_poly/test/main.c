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
#include <stdlib.h>

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

test_struct tests[] =
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

/* main function *************************************************************/

TEST_MAIN(tests)
