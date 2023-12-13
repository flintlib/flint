/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Try to get fdopen declared for t-print_read.c */
#if defined __STRICT_ANSI__
# undef __STRICT_ANSI__
#endif

#include <string.h>
#include <stdlib.h>

/* Include functions *********************************************************/

#include "t-add.c"
#include "t-add_series.c"
#include "t-berlekamp_massey.c"
#include "t-compose.c"
#include "t-compose_mod_brent_kung.c"
#include "t-compose_mod_brent_kung_precomp_preinv.c"
#include "t-compose_mod_brent_kung_precomp_preinv_threaded.c"
#include "t-compose_mod_brent_kung_preinv.c"
#include "t-compose_mod_brent_kung_vec_preinv.c"
#include "t-compose_mod_brent_kung_vec_preinv_threaded.c"
#include "t-compose_mod.c"
#include "t-compose_mod_horner.c"
#include "t-deflate_deflation_inflate.c"
#include "t-derivative.c"
#include "t-discriminant.c"
#include "t-div.c"
#include "t-divides.c"
#include "t-divides_classical.c"
#include "t-div_newton_n_preinv.c"
#include "t-divrem_basecase.c"
#include "t-divrem_f.c"
#include "t-divrem_newton_n_preinv.c"
#include "t-div_series.c"
#include "t-evaluate_fmpz.c"
#include "t-evaluate_fmpz_vec.c"
#include "t-find_distinct_nonzero_roots.c"
#include "t-frobenius_powers_precomp.c"
#include "t-gcd.c"
#include "t-gcd_euclidean_f.c"
#include "t-gcdinv.c"
#include "t-gcdinv_euclidean.c"
#include "t-get_set_fmpz_poly.c"
#include "t-get_set_nmod_poly.c"
#include "t-init_realloc_clear.c"
#include "t-invmod.c"
#include "t-inv_series.c"
#include "t-invsqrt_series.c"
#include "t-minpoly.c"
#include "t-mul.c"
#include "t-mulhigh.c"
#include "t-mullow.c"
#include "t-mulmod.c"
#include "t-mulmod_preinv.c"
#include "t-neg.c"
#include "t-powers_mod_bsgs.c"
#include "t-powers_mod_naive.c"
#include "t-powmod_fmpz_binexp.c"
#include "t-powmod_fmpz_binexp_preinv.c"
#include "t-powmod_ui_binexp.c"
#include "t-powmod_ui_binexp_preinv.c"
#include "t-powmod_x_fmpz_preinv.c"
#include "t-pow_trunc_binexp.c"
#include "t-pow_trunc.c"
#include "t-print_read.c"
#include "t-product_roots_fmpz_vec.c"
#include "t-radix.c"
#include "t-randtest_monic_primitive.c"
#include "t-rem_basecase.c"
#include "t-resultant.c"
#include "t-resultant_euclidean.c"
#include "t-resultant_hgcd.c"
#include "t-scalar_div_fmpz.c"
#include "t-scalar_mul_fmpz.c"
#include "t-set_equal.c"
#include "t-set_trunc.c"
#include "t-shift_left_right.c"
#include "t-sqrt.c"
#include "t-sqrt_series.c"
#include "t-sub.c"
#include "t-sub_series.c"
#include "t-swap.c"
#include "t-xgcd.c"
#include "t-zero.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(fmpz_mod_poly_add),
    TEST_FUNCTION(fmpz_mod_poly_add_series),
    TEST_FUNCTION(fmpz_mod_poly_berlekamp_massey),
    TEST_FUNCTION(fmpz_mod_poly_compose),
    TEST_FUNCTION(fmpz_mod_poly_compose_mod_brent_kung),
    TEST_FUNCTION(fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv),
    TEST_FUNCTION(fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv_threaded),
    TEST_FUNCTION(fmpz_mod_poly_compose_mod_brent_kung_preinv),
    TEST_FUNCTION(fmpz_mod_poly_compose_mod_brent_kung_vec_preinv),
    TEST_FUNCTION(fmpz_mod_poly_compose_mod_brent_kung_vec_preinv_threaded),
    TEST_FUNCTION(fmpz_mod_poly_compose_mod),
    TEST_FUNCTION(fmpz_mod_poly_compose_mod_horner),
    TEST_FUNCTION(fmpz_mod_poly_deflate_deflation_inflate),
    TEST_FUNCTION(fmpz_mod_poly_derivative),
    TEST_FUNCTION(fmpz_mod_poly_discriminant),
    TEST_FUNCTION(fmpz_mod_poly_div),
    TEST_FUNCTION(fmpz_mod_poly_divides),
    TEST_FUNCTION(fmpz_mod_poly_divides_classical),
    TEST_FUNCTION(fmpz_mod_poly_div_newton_n_preinv),
    TEST_FUNCTION(fmpz_mod_poly_divrem_basecase),
    TEST_FUNCTION(fmpz_mod_poly_divrem_f),
    TEST_FUNCTION(fmpz_mod_poly_divrem_newton_n_preinv),
    TEST_FUNCTION(fmpz_mod_poly_div_series),
    TEST_FUNCTION(fmpz_mod_poly_evaluate_fmpz),
    TEST_FUNCTION(fmpz_mod_poly_evaluate_fmpz_vec),
    TEST_FUNCTION(fmpz_mod_poly_find_distinct_nonzero_roots),
    TEST_FUNCTION(fmpz_mod_poly_frobenius_powers_precomp),
    TEST_FUNCTION(fmpz_mod_poly_gcd),
    TEST_FUNCTION(fmpz_mod_poly_gcd_euclidean_f),
    TEST_FUNCTION(fmpz_mod_poly_gcdinv),
    TEST_FUNCTION(fmpz_mod_poly_gcdinv_euclidean),
    TEST_FUNCTION(fmpz_mod_poly_get_set_fmpz_poly),
    TEST_FUNCTION(fmpz_mod_poly_get_set_nmod_poly),
    TEST_FUNCTION(fmpz_mod_poly_init_realloc_clear),
    TEST_FUNCTION(fmpz_mod_poly_invmod),
    TEST_FUNCTION(fmpz_mod_poly_inv_series),
    TEST_FUNCTION(fmpz_mod_poly_invsqrt_series),
    TEST_FUNCTION(fmpz_mod_poly_minpoly),
    TEST_FUNCTION(fmpz_mod_poly_mul),
    TEST_FUNCTION(fmpz_mod_poly_mulhigh),
    TEST_FUNCTION(fmpz_mod_poly_mullow),
    TEST_FUNCTION(fmpz_mod_poly_mulmod),
    TEST_FUNCTION(fmpz_mod_poly_mulmod_preinv),
    TEST_FUNCTION(fmpz_mod_poly_neg),
    TEST_FUNCTION(fmpz_mod_poly_powers_mod_bsgs),
    TEST_FUNCTION(fmpz_mod_poly_powers_mod_naive),
    TEST_FUNCTION(fmpz_mod_poly_powmod_fmpz_binexp),
    TEST_FUNCTION(fmpz_mod_poly_powmod_fmpz_binexp_preinv),
    TEST_FUNCTION(fmpz_mod_poly_powmod_ui_binexp),
    TEST_FUNCTION(fmpz_mod_poly_powmod_ui_binexp_preinv),
    TEST_FUNCTION(fmpz_mod_poly_powmod_x_fmpz_preinv),
    TEST_FUNCTION(fmpz_mod_poly_pow_trunc_binexp),
    TEST_FUNCTION(fmpz_mod_poly_pow_trunc),
    TEST_FUNCTION(fmpz_mod_poly_print_read),
    TEST_FUNCTION(fmpz_mod_poly_product_roots_fmpz_vec),
    TEST_FUNCTION(fmpz_mod_poly_radix),
    TEST_FUNCTION(fmpz_mod_poly_randtest_monic_primitive),
    TEST_FUNCTION(fmpz_mod_poly_rem_basecase),
    TEST_FUNCTION(fmpz_mod_poly_resultant),
    TEST_FUNCTION(fmpz_mod_poly_resultant_euclidean),
    TEST_FUNCTION(fmpz_mod_poly_resultant_hgcd),
    TEST_FUNCTION(fmpz_mod_poly_scalar_div_fmpz),
    TEST_FUNCTION(fmpz_mod_poly_scalar_mul_fmpz),
    TEST_FUNCTION(fmpz_mod_poly_set_equal),
    TEST_FUNCTION(fmpz_mod_poly_set_trunc),
    TEST_FUNCTION(fmpz_mod_poly_shift_left_right),
    TEST_FUNCTION(fmpz_mod_poly_sqrt),
    TEST_FUNCTION(fmpz_mod_poly_sqrt_series),
    TEST_FUNCTION(fmpz_mod_poly_sub),
    TEST_FUNCTION(fmpz_mod_poly_sub_series),
    TEST_FUNCTION(fmpz_mod_poly_swap),
    TEST_FUNCTION(fmpz_mod_poly_xgcd),
    TEST_FUNCTION(fmpz_mod_poly_zero)
};

/* main function *************************************************************/

TEST_MAIN(tests)
