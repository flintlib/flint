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

#include "t-atan_series.c"
#include "t-compose.c"
#include "t-compose_divconquer.c"
#include "t-compose_horner.c"
#include "t-compose_series.c"
#include "t-div_basecase.c"
#include "t-div.c"
#include "t-div_divconquer.c"
#include "t-divexact.c"
#include "t-div_newton.c"
#include "t-divrem_basecase.c"
#include "t-divrem.c"
#include "t-divrem_divconquer.c"
#include "t-divrem_newton.c"
#include "t-div_series.c"
#include "t-evaluate.c"
#include "t-evaluate_horner.c"
#include "t-evaluate_modular.c"
#include "t-evaluate_other.c"
#include "t-evaluate_other_rectangular.c"
#include "t-evaluate_rectangular.c"
#include "t-evaluate_vec_fast.c"
#include "t-exp_series.c"
#include "t-factor_squarefree.c"
#include "t-gcd.c"
#include "t-gcd_euclidean.c"
#include "t-gcd_hgcd.c"
#include "t-hgcd.c"
#include "t-integral.c"
#include "t-inv_series.c"
#include "t-log_series.c"
#include "t-make_monic.c"
#include "t-nth_derivative.c"
#include "t-pow_series_fmpq.c"
#include "t-pow_series_ui.c"
#include "t-pow_ui.c"
#include "t-rem.c"
#include "t-resultant.c"
#include "t-resultant_euclidean.c"
#include "t-resultant_hgcd.c"
#include "t-resultant_sylvester.c"
#include "t-revert_series.c"
#include "t-roots.c"
#include "t-roots_other.c"
#include "t-rsqrt_series.c"
#include "t-shift_left_right.c"
#include "t-sqrt_series.c"
#include "t-squarefree_part.c"
#include "t-taylor_shift.c"
#include "t-taylor_shift_convolution.c"
#include "t-taylor_shift_divconquer.c"
#include "t-taylor_shift_horner.c"
#include "t-xgcd_euclidean.c"
#include "t-xgcd_hgcd.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(gr_poly_atan_series),
    TEST_FUNCTION(gr_poly_compose),
    TEST_FUNCTION(gr_poly_compose_divconquer),
    TEST_FUNCTION(gr_poly_compose_horner),
    TEST_FUNCTION(gr_poly_compose_series),
    TEST_FUNCTION(gr_poly_div_basecase),
    TEST_FUNCTION(gr_poly_div),
    TEST_FUNCTION(gr_poly_div_divconquer),
    TEST_FUNCTION(gr_poly_divexact),
    TEST_FUNCTION(gr_poly_div_newton),
    TEST_FUNCTION(gr_poly_divrem_basecase),
    TEST_FUNCTION(gr_poly_divrem),
    TEST_FUNCTION(gr_poly_divrem_divconquer),
    TEST_FUNCTION(gr_poly_divrem_newton),
    TEST_FUNCTION(gr_poly_div_series),
    TEST_FUNCTION(gr_poly_evaluate),
    TEST_FUNCTION(gr_poly_evaluate_horner),
    TEST_FUNCTION(gr_poly_evaluate_modular),
    TEST_FUNCTION(gr_poly_evaluate_other),
    TEST_FUNCTION(gr_poly_evaluate_other_rectangular),
    TEST_FUNCTION(gr_poly_evaluate_rectangular),
    TEST_FUNCTION(gr_poly_evaluate_vec_fast),
    TEST_FUNCTION(gr_poly_exp_series),
    TEST_FUNCTION(gr_poly_factor_squarefree),
    TEST_FUNCTION(gr_poly_gcd),
    TEST_FUNCTION(gr_poly_gcd_euclidean),
    TEST_FUNCTION(gr_poly_gcd_hgcd),
    TEST_FUNCTION(gr_poly_hgcd),
    TEST_FUNCTION(gr_poly_integral),
    TEST_FUNCTION(gr_poly_inv_series),
    TEST_FUNCTION(gr_poly_log_series),
    TEST_FUNCTION(gr_poly_make_monic),
    TEST_FUNCTION(gr_poly_nth_derivative),
    TEST_FUNCTION(gr_poly_pow_series_fmpq),
    TEST_FUNCTION(gr_poly_pow_series_ui),
    TEST_FUNCTION(gr_poly_pow_ui),
    TEST_FUNCTION(gr_poly_rem),
    TEST_FUNCTION(gr_poly_resultant),
    TEST_FUNCTION(gr_poly_resultant_euclidean),
    TEST_FUNCTION(gr_poly_resultant_hgcd),
    TEST_FUNCTION(gr_poly_resultant_sylvester),
    TEST_FUNCTION(gr_poly_revert_series),
    TEST_FUNCTION(gr_poly_roots),
    TEST_FUNCTION(gr_poly_roots_other),
    TEST_FUNCTION(gr_poly_rsqrt_series),
    TEST_FUNCTION(gr_poly_shift_left_right),
    TEST_FUNCTION(gr_poly_sqrt_series),
    TEST_FUNCTION(gr_poly_squarefree_part),
    TEST_FUNCTION(gr_poly_taylor_shift),
    TEST_FUNCTION(gr_poly_taylor_shift_convolution),
    TEST_FUNCTION(gr_poly_taylor_shift_divconquer),
    TEST_FUNCTION(gr_poly_taylor_shift_horner),
    TEST_FUNCTION(gr_poly_xgcd_euclidean),
    TEST_FUNCTION(gr_poly_xgcd_hgcd)
};

/* main function *************************************************************/

TEST_MAIN(tests)
