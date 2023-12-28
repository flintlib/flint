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
#include "t-add_sub_ui.c"
#include "t-cmp.c"
#include "t-compose_nmod_mpoly.c"
#include "t-compose_nmod_poly.c"
#include "t-content_vars.c"
#include "t-degree.c"
#include "t-derivative.c"
#include "t-div.c"
#include "t-divides.c"
#include "t-divides_dense.c"
#include "t-divides_heap_threaded.c"
#include "t-divides_monagan_pearce.c"
#include "t-div_monagan_pearce.c"
#include "t-divrem.c"
#include "t-divrem_ideal.c"
#include "t-divrem_ideal_monagan_pearce.c"
#include "t-divrem_monagan_pearce.c"
#include "t-evaluate.c"
#include "t-gcd_brown.c"
#include "t-gcd.c"
#include "t-gcd_cofactors.c"
#include "t-gcd_hensel.c"
#include "t-gcd_zippel2.c"
#include "t-gcd_zippel.c"
#include "t-gen.c"
#include "t-get_coeff_vars_ui.c"
#include "t-get_set_string_pretty.c"
#include "t-get_set_term_coeff_ui.c"
#include "t-get_set_term_exp_si.c"
#include "t-get_set_term_exp_ui.c"
#include "t-get_term.c"
#include "t-get_term_monomial.c"
#include "t-inflate_deflate.c"
#include "t-mpolyn_divides_threaded_pool.c"
#include "t-mpolyuu_divides.c"
#include "t-mul_array.c"
#include "t-mul_array_threaded.c"
#include "t-mul.c"
#include "t-mul_dense.c"
#include "t-mul_heap_threaded.c"
#include "t-mul_johnson.c"
#include "t-pow_rmul.c"
#include "t-pow_ui.c"
#include "t-push_term_ui_fmpz.c"
#include "t-push_term_ui_ui.c"
#include "t-quadratic_root.c"
#include "t-repack_bits.c"
#include "t-resize.c"
#include "t-resultant_discriminant.c"
#include "t-scalar_addmul_ui.c"
#include "t-scalar_mul_ui.c"
#include "t-sqrt.c"
#include "t-term_content.c"
#include "t-total_degree.c"
#include "t-univar.c"
#include "t-univar_resultant.c"
#include "t-used_vars.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(nmod_mpoly_add_sub),
    TEST_FUNCTION(nmod_mpoly_add_sub_ui),
    TEST_FUNCTION(nmod_mpoly_cmp),
    TEST_FUNCTION(nmod_mpoly_compose_nmod_mpoly),
    TEST_FUNCTION(nmod_mpoly_compose_nmod_poly),
    TEST_FUNCTION(nmod_mpoly_content_vars),
    TEST_FUNCTION(nmod_mpoly_degree),
    TEST_FUNCTION(nmod_mpoly_derivative),
    TEST_FUNCTION(nmod_mpoly_div),
    TEST_FUNCTION(nmod_mpoly_divides),
    TEST_FUNCTION(nmod_mpoly_divides_dense),
    TEST_FUNCTION(nmod_mpoly_divides_heap_threaded),
    TEST_FUNCTION(nmod_mpoly_divides_monagan_pearce),
    TEST_FUNCTION(nmod_mpoly_div_monagan_pearce),
    TEST_FUNCTION(nmod_mpoly_divrem),
    TEST_FUNCTION(nmod_mpoly_divrem_ideal),
    TEST_FUNCTION(nmod_mpoly_divrem_ideal_monagan_pearce),
    TEST_FUNCTION(nmod_mpoly_divrem_monagan_pearce),
    TEST_FUNCTION(nmod_mpoly_evaluate),
    TEST_FUNCTION(nmod_mpoly_gcd_brown),
    TEST_FUNCTION(nmod_mpoly_gcd),
    TEST_FUNCTION(nmod_mpoly_gcd_cofactors),
    TEST_FUNCTION(nmod_mpoly_gcd_hensel),
    TEST_FUNCTION(nmod_mpoly_gcd_zippel2),
    TEST_FUNCTION(nmod_mpoly_gcd_zippel),
    TEST_FUNCTION(nmod_mpoly_gen),
    TEST_FUNCTION(nmod_mpoly_get_coeff_vars_ui),
    TEST_FUNCTION(nmod_mpoly_get_set_string_pretty),
    TEST_FUNCTION(nmod_mpoly_get_set_term_coeff_ui),
    TEST_FUNCTION(nmod_mpoly_get_set_term_exp_si),
    TEST_FUNCTION(nmod_mpoly_get_set_term_exp_ui),
    TEST_FUNCTION(nmod_mpoly_get_term),
    TEST_FUNCTION(nmod_mpoly_get_term_monomial),
    TEST_FUNCTION(nmod_mpoly_inflate_deflate),
    TEST_FUNCTION(nmod_mpolyn_divides_threaded_pool),
    TEST_FUNCTION(nmod_mpoly_mpolyuu_divides),
    TEST_FUNCTION(nmod_mpoly_mul_array),
    TEST_FUNCTION(nmod_mpoly_mul_array_threaded),
    TEST_FUNCTION(nmod_mpoly_mul),
    TEST_FUNCTION(nmod_mpoly_mul_dense),
    TEST_FUNCTION(nmod_mpoly_mul_heap_threaded),
    TEST_FUNCTION(nmod_mpoly_mul_johnson),
    TEST_FUNCTION(nmod_mpoly_pow_rmul),
    TEST_FUNCTION(nmod_mpoly_pow_ui),
    TEST_FUNCTION(nmod_mpoly_push_term_ui_fmpz),
    TEST_FUNCTION(nmod_mpoly_push_term_ui_ui),
    TEST_FUNCTION(nmod_mpoly_quadratic_root),
    TEST_FUNCTION(nmod_mpoly_repack_bits),
    TEST_FUNCTION(nmod_mpoly_resize),
    TEST_FUNCTION(nmod_mpoly_resultant_discriminant),
    TEST_FUNCTION(nmod_mpoly_scalar_addmul_ui),
    TEST_FUNCTION(nmod_mpoly_scalar_mul_ui),
    TEST_FUNCTION(nmod_mpoly_sqrt),
    TEST_FUNCTION(nmod_mpoly_term_content),
    TEST_FUNCTION(nmod_mpoly_total_degree),
    TEST_FUNCTION(nmod_mpoly_univar),
    TEST_FUNCTION(nmod_mpoly_univar_resultant),
    TEST_FUNCTION(nmod_mpoly_used_vars)
};

/* main function *************************************************************/

TEST_MAIN(tests)
