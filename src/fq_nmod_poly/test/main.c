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

#include "t-add.c"
#include "t-add_series.c"
#include "t-compose.c"
#include "t-compose_mod_brent_kung.c"
#include "t-compose_mod_brent_kung_preinv.c"
#include "t-compose_mod.c"
#include "t-compose_mod_horner.c"
#include "t-compose_mod_horner_preinv.c"
#include "t-compose_mod_preinv.c"
#include "t-deflate.c"
#include "t-derivative.c"
#include "t-div.c"
#include "t-divides.c"
#include "t-div_newton_n_preinv.c"
#include "t-divrem.c"
#include "t-divrem_newton_n_preinv.c"
#include "t-div_series.c"
#include "t-equal_trunc.c"
#include "t-evaluate_fq.c"
#include "t-evaluate_fq_vec_fast.c"
#include "t-gcd.c"
#include "t-gcd_euclidean_f.c"
#include "t-get_str.c"
#include "t-get_str_pretty.c"
#include "t-hamming_weight.c"
#include "t-inflate.c"
#include "t-inv_series_newton.c"
#include "t-invsqrt_series.c"
#include "t-make_monic.c"
#include "t-mul.c"
#include "t-mul_classical.c"
#include "t-mulhigh.c"
#include "t-mulhigh_classical.c"
#include "t-mul_KS.c"
#include "t-mullow.c"
#include "t-mullow_classical.c"
#include "t-mullow_KS.c"
#include "t-mullow_univariate.c"
#include "t-mulmod.c"
#include "t-mulmod_preinv.c"
#include "t-mul_univariate.c"
#include "t-neg.c"
#include "t-pow.c"
#include "t-powmod_fmpz_binexp.c"
#include "t-powmod_fmpz_binexp_preinv.c"
#include "t-powmod_fmpz_sliding_preinv.c"
#include "t-powmod_ui_binexp.c"
#include "t-powmod_ui_binexp_preinv.c"
#include "t-powmod_x_fmpz_preinv.c"
#include "t-pow_trunc_binexp.c"
#include "t-pow_trunc.c"
#include "t-randtest_irreducible.c"
#include "t-scalar_addmul_fq.c"
#include "t-scalar_div_fq.c"
#include "t-scalar_mul_fq.c"
#include "t-scalar_submul_fq.c"
#include "t-set_fmpz_mod_poly.c"
#include "t-set_nmod_poly.c"
#include "t-set_trunc.c"
#include "t-shift_left_right.c"
#include "t-sqr.c"
#include "t-sqr_classical.c"
#include "t-sqr_KS.c"
#include "t-sqrt.c"
#include "t-sqrt_series.c"
#include "t-sub.c"
#include "t-sub_series.c"
#include "t-xgcd.c"
#include "t-xgcd_euclidean_f.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(fq_nmod_poly_add),
    TEST_FUNCTION(fq_nmod_poly_add_series),
    TEST_FUNCTION(fq_nmod_poly_compose),
    TEST_FUNCTION(fq_nmod_poly_compose_mod_brent_kung),
    TEST_FUNCTION(fq_nmod_poly_compose_mod_brent_kung_preinv),
    TEST_FUNCTION(fq_nmod_poly_compose_mod),
    TEST_FUNCTION(fq_nmod_poly_compose_mod_horner),
    TEST_FUNCTION(fq_nmod_poly_compose_mod_horner_preinv),
    TEST_FUNCTION(fq_nmod_poly_compose_mod_preinv),
    TEST_FUNCTION(fq_nmod_poly_deflate),
    TEST_FUNCTION(fq_nmod_poly_derivative),
    TEST_FUNCTION(fq_nmod_poly_div),
    TEST_FUNCTION(fq_nmod_poly_divides),
    TEST_FUNCTION(fq_nmod_poly_div_newton_n_preinv),
    TEST_FUNCTION(fq_nmod_poly_divrem),
    TEST_FUNCTION(fq_nmod_poly_divrem_newton_n_preinv),
    TEST_FUNCTION(fq_nmod_poly_div_series),
    TEST_FUNCTION(fq_nmod_poly_equal_trunc),
    TEST_FUNCTION(fq_nmod_poly_evaluate_fq),
    TEST_FUNCTION(fq_nmod_poly_evaluate_fq_nmod_vec_fast),
    TEST_FUNCTION(fq_nmod_poly_gcd),
    TEST_FUNCTION(fq_nmod_poly_gcd_euclidean_f),
    TEST_FUNCTION(fq_nmod_poly_get_str),
    TEST_FUNCTION(fq_nmod_poly_get_str_pretty),
    TEST_FUNCTION(fq_nmod_poly_hamming_weight),
    TEST_FUNCTION(fq_nmod_poly_inflate),
    TEST_FUNCTION(fq_nmod_poly_inv_series_newton),
    TEST_FUNCTION(fq_nmod_poly_invsqrt_series),
    TEST_FUNCTION(fq_nmod_poly_make_monic),
    TEST_FUNCTION(fq_nmod_poly_mul),
    TEST_FUNCTION(fq_nmod_poly_mul_classical),
    TEST_FUNCTION(fq_nmod_poly_mulhigh),
    TEST_FUNCTION(fq_nmod_poly_mulhigh_classical),
    TEST_FUNCTION(fq_nmod_poly_mul_KS),
    TEST_FUNCTION(fq_nmod_poly_mullow),
    TEST_FUNCTION(fq_nmod_poly_mullow_classical),
    TEST_FUNCTION(fq_nmod_poly_mullow_KS),
    TEST_FUNCTION(fq_nmod_poly_mullow_univariate),
    TEST_FUNCTION(fq_nmod_poly_mulmod),
    TEST_FUNCTION(fq_nmod_poly_mulmod_preinv),
    TEST_FUNCTION(fq_nmod_poly_mul_univariate),
    TEST_FUNCTION(fq_nmod_poly_neg),
    TEST_FUNCTION(fq_nmod_poly_pow),
    TEST_FUNCTION(fq_nmod_poly_powmod_fmpz_binexp),
    TEST_FUNCTION(fq_nmod_poly_powmod_fmpz_binexp_preinv),
    TEST_FUNCTION(fq_nmod_poly_powmod_fmpz_sliding_preinv),
    TEST_FUNCTION(fq_nmod_poly_powmod_ui_binexp),
    TEST_FUNCTION(fq_nmod_poly_powmod_ui_binexp_preinv),
    TEST_FUNCTION(fq_nmod_poly_powmod_x_fmpz_preinv),
    TEST_FUNCTION(fq_nmod_poly_pow_trunc_binexp),
    TEST_FUNCTION(fq_nmod_poly_pow_trunc),
    TEST_FUNCTION(fq_nmod_poly_randtest_irreducible),
    TEST_FUNCTION(fq_nmod_poly_scalar_addmul_fq),
    TEST_FUNCTION(fq_nmod_poly_scalar_div_fq),
    TEST_FUNCTION(fq_nmod_poly_scalar_mul_fq),
    TEST_FUNCTION(fq_nmod_poly_scalar_submul_fq),
    TEST_FUNCTION(fq_nmod_poly_set_fmpz_mod_poly),
    TEST_FUNCTION(fq_nmod_poly_set_nmod_poly),
    TEST_FUNCTION(fq_nmod_poly_set_trunc),
    TEST_FUNCTION(fq_nmod_poly_shift_left_right),
    TEST_FUNCTION(fq_nmod_poly_sqr),
    TEST_FUNCTION(fq_nmod_poly_sqr_classical),
    TEST_FUNCTION(fq_nmod_poly_sqr_KS),
    TEST_FUNCTION(fq_nmod_poly_sqrt),
    TEST_FUNCTION(fq_nmod_poly_sqrt_series),
    TEST_FUNCTION(fq_nmod_poly_sub),
    TEST_FUNCTION(fq_nmod_poly_sub_series),
    TEST_FUNCTION(fq_nmod_poly_xgcd),
    TEST_FUNCTION(fq_nmod_poly_xgcd_euclidean_f)
};

/* main function *************************************************************/

TEST_MAIN(tests)
