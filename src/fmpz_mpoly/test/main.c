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

#include "t-add_sub.c"
#include "t-add_sub_fmpz.c"
#include "t-add_sub_si.c"
#include "t-add_sub_ui.c"
#include "t-buchberger_naive.c"
#include "t-cmp.c"
#include "t-compose_fmpz_mpoly.c"
#include "t-compose_fmpz_poly.c"
#include "t-content_vars.c"
#include "t-degree.c"
#include "t-degrees_term_exp_fits_ui_si.c"
#include "t-derivative_integral.c"
#include "t-divides_array.c"
#include "t-divides.c"
#include "t-divides_heap_threaded.c"
#include "t-divides_monagan_pearce.c"
#include "t-div_monagan_pearce.c"
#include "t-divrem_array.c"
#include "t-divrem_ideal_monagan_pearce.c"
#include "t-divrem_monagan_pearce.c"
#include "t-evaluate.c"
#include "t-gcd_brown.c"
#include "t-gcd.c"
#include "t-gcd_cofactors.c"
#include "t-gcd_hensel.c"
#include "t-gcd_subresultant.c"
#include "t-gcd_zippel2.c"
#include "t-gcd_zippel.c"
#include "t-gen.c"
#include "t-get_coeff_vars_ui.c"
#include "t-get_set_coeff.c"
#include "t-get_set_coeff_fmpz_monomial.c"
#include "t-get_set_is_fmpz.c"
#include "t-get_set_str_pretty.c"
#include "t-get_set_term_coeff_fmpz.c"
#include "t-get_set_term_exp_fmpz.c"
#include "t-get_set_term_exp_si.c"
#include "t-get_set_term_exp_ui.c"
#include "t-get_term.c"
#include "t-get_term_monomial.c"
#include "t-inflate_deflate.c"
#include "t-init.c"
#include "t-mul_array.c"
#include "t-mul_array_threaded.c"
#include "t-mul.c"
#include "t-mul_dense.c"
#include "t-mul_heap_threaded.c"
#include "t-mul_johnson.c"
#include "t-mul_monomial.c"
#include "t-neg.c"
#include "t-pow_fps.c"
#include "t-pow_ui.c"
#include "t-push_term_fmpz_fmpz.c"
#include "t-push_term_fmpz_ui.c"
#include "t-quasidiv_heap.c"
#include "t-quasidivrem_heap.c"
#include "t-quasidivrem_ideal_heap.c"
#include "t-repack_bits.c"
#include "t-resize.c"
#include "t-resultant_discriminant.c"
#include "t-reverse.c"
#include "t-scalar_divexact_fmpz.c"
#include "t-scalar_divexact_si.c"
#include "t-scalar_divexact_ui.c"
#include "t-scalar_divides_fmpz.c"
#include "t-scalar_fmma.c"
#include "t-scalar_mul_fmpz.c"
#include "t-scalar_mul_si.c"
#include "t-scalar_mul_ui.c"
#include "t-set_equal.c"
#include "t-set_fmpz.c"
#include "t-set_gen_fmpz_poly.c"
#include "t-set_si.c"
#include "t-set_ui.c"
#include "t-sort_terms.c"
#include "t-sqrt_heap.c"
#include "t-symmetric.c"
#include "t-term_content.c"
#include "t-total_degree.c"
#include "t-univar.c"
#include "t-univar_resultant.c"
#include "t-used_vars.c"
#include "t-vec_autoreduction.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(fmpz_mpoly_add_sub),
    TEST_FUNCTION(fmpz_mpoly_add_sub_fmpz),
    TEST_FUNCTION(fmpz_mpoly_add_sub_si),
    TEST_FUNCTION(fmpz_mpoly_add_sub_ui),
    TEST_FUNCTION(fmpz_mpoly_buchberger_naive),
    TEST_FUNCTION(fmpz_mpoly_cmp),
    TEST_FUNCTION(fmpz_mpoly_compose_fmpz_mpoly),
    TEST_FUNCTION(fmpz_mpoly_compose_fmpz_poly),
    TEST_FUNCTION(fmpz_mpoly_content_vars),
    TEST_FUNCTION(fmpz_mpoly_degree),
    TEST_FUNCTION(fmpz_mpoly_degrees_term_exp_fits_ui_si),
    TEST_FUNCTION(fmpz_mpoly_derivative_integral),
    TEST_FUNCTION(fmpz_mpoly_divides_array),
    TEST_FUNCTION(fmpz_mpoly_divides),
    TEST_FUNCTION(fmpz_mpoly_divides_heap_threaded),
    TEST_FUNCTION(fmpz_mpoly_divides_monagan_pearce),
    TEST_FUNCTION(fmpz_mpoly_div_monagan_pearce),
    TEST_FUNCTION(fmpz_mpoly_divrem_array),
    TEST_FUNCTION(fmpz_mpoly_divrem_ideal_monagan_pearce),
    TEST_FUNCTION(fmpz_mpoly_divrem_monagan_pearce),
    TEST_FUNCTION(fmpz_mpoly_evaluate),
    TEST_FUNCTION(fmpz_mpoly_gcd_brown),
    TEST_FUNCTION(fmpz_mpoly_gcd),
    TEST_FUNCTION(fmpz_mpoly_gcd_cofactors),
    TEST_FUNCTION(fmpz_mpoly_gcd_hensel),
    TEST_FUNCTION(fmpz_mpoly_gcd_subresultant),
    TEST_FUNCTION(fmpz_mpoly_gcd_zippel2),
    TEST_FUNCTION(fmpz_mpoly_gcd_zippel),
    TEST_FUNCTION(fmpz_mpoly_gen),
    TEST_FUNCTION(fmpz_mpoly_get_coeff_vars_ui),
    TEST_FUNCTION(fmpz_mpoly_get_set_coeff),
    TEST_FUNCTION(fmpz_mpoly_get_set_coeff_fmpz_monomial),
    TEST_FUNCTION(fmpz_mpoly_get_set_is_fmpz),
    TEST_FUNCTION(fmpz_mpoly_get_set_str_pretty),
    TEST_FUNCTION(fmpz_mpoly_get_set_term_coeff_fmpz),
    TEST_FUNCTION(fmpz_mpoly_get_set_term_exp_fmpz),
    TEST_FUNCTION(fmpz_mpoly_get_set_term_exp_si),
    TEST_FUNCTION(fmpz_mpoly_get_set_term_exp_ui),
    TEST_FUNCTION(fmpz_mpoly_get_term),
    TEST_FUNCTION(fmpz_mpoly_get_term_monomial),
    TEST_FUNCTION(fmpz_mpoly_inflate_deflate),
    TEST_FUNCTION(fmpz_mpoly_init),
    TEST_FUNCTION(fmpz_mpoly_mul_array),
    TEST_FUNCTION(fmpz_mpoly_mul_array_threaded),
    TEST_FUNCTION(fmpz_mpoly_mul),
    TEST_FUNCTION(fmpz_mpoly_mul_dense),
    TEST_FUNCTION(fmpz_mpoly_mul_heap_threaded),
    TEST_FUNCTION(fmpz_mpoly_mul_johnson),
    TEST_FUNCTION(fmpz_mpoly_mul_monomial),
    TEST_FUNCTION(fmpz_mpoly_neg),
    TEST_FUNCTION(fmpz_mpoly_pow_fps),
    TEST_FUNCTION(fmpz_mpoly_pow_ui),
    TEST_FUNCTION(fmpz_mpoly_push_term_fmpz_fmpz),
    TEST_FUNCTION(fmpz_mpoly_push_term_fmpz_ui),
    TEST_FUNCTION(fmpz_mpoly_quasidiv_heap),
    TEST_FUNCTION(fmpz_mpoly_quasidivrem_heap),
    TEST_FUNCTION(fmpz_mpoly_quasidivrem_ideal_heap),
    TEST_FUNCTION(fmpz_mpoly_repack_bits),
    TEST_FUNCTION(fmpz_mpoly_resize),
    TEST_FUNCTION(fmpz_mpoly_resultant_discriminant),
    TEST_FUNCTION(fmpz_mpoly_reverse),
    TEST_FUNCTION(fmpz_mpoly_scalar_divexact_fmpz),
    TEST_FUNCTION(fmpz_mpoly_scalar_divexact_si),
    TEST_FUNCTION(fmpz_mpoly_scalar_divexact_ui),
    TEST_FUNCTION(fmpz_mpoly_scalar_divides_fmpz),
    TEST_FUNCTION(fmpz_mpoly_scalar_fmma),
    TEST_FUNCTION(fmpz_mpoly_scalar_mul_fmpz),
    TEST_FUNCTION(fmpz_mpoly_scalar_mul_si),
    TEST_FUNCTION(fmpz_mpoly_scalar_mul_ui),
    TEST_FUNCTION(fmpz_mpoly_set_equal),
    TEST_FUNCTION(fmpz_mpoly_set_fmpz),
    TEST_FUNCTION(fmpz_mpoly_set_gen_fmpz_poly),
    TEST_FUNCTION(fmpz_mpoly_set_si),
    TEST_FUNCTION(fmpz_mpoly_set_ui),
    TEST_FUNCTION(fmpz_mpoly_sort_terms),
    TEST_FUNCTION(fmpz_mpoly_sqrt_heap),
    TEST_FUNCTION(fmpz_mpoly_symmetric),
    TEST_FUNCTION(fmpz_mpoly_term_content),
    TEST_FUNCTION(fmpz_mpoly_total_degree),
    TEST_FUNCTION(fmpz_mpoly_univar),
    TEST_FUNCTION(fmpz_mpoly_univar_resultant),
    TEST_FUNCTION(fmpz_mpoly_used_vars),
    TEST_FUNCTION(fmpz_mpoly_vec_autoreduction)
};

/* main function *************************************************************/

TEST_MAIN(tests)
