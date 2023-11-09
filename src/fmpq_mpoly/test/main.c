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
#include "t-add_sub_fmpq.c"
#include "t-add_sub_fmpz.c"
#include "t-cmp.c"
#include "t-compose_fmpq_mpoly.c"
#include "t-compose_fmpq_poly.c"
#include "t-content.c"
#include "t-content_vars.c"
#include "t-degree.c"
#include "t-degrees_term_exp_fits_ui_si.c"
#include "t-derivative_integral.c"
#include "t-div.c"
#include "t-divides.c"
#include "t-divrem.c"
#include "t-divrem_ideal.c"
#include "t-equal_is_fmpq.c"
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
#include "t-get_set_coeff_fmpq_fmpz.c"
#include "t-get_set_coeff_fmpq_monomial.c"
#include "t-get_set_coeff_fmpq_ui.c"
#include "t-get_set_is_fmpq.c"
#include "t-get_set_str_pretty.c"
#include "t-get_set_term_coeff_fmpq.c"
#include "t-get_set_term_exp_fmpz.c"
#include "t-get_set_term_exp_si.c"
#include "t-get_set_term_exp_ui.c"
#include "t-get_term.c"
#include "t-get_term_monomial.c"
#include "t-mul.c"
#include "t-pow_fmpz.c"
#include "t-pow_ui.c"
#include "t-push_term_fmpq_fmpz.c"
#include "t-push_term_fmpq_ui.c"
#include "t-resultant_discriminant.c"
#include "t-scalar_mul_div_fmpq.c"
#include "t-scalar_mul_div_fmpz.c"
#include "t-scalar_mul_fmpq.c"
#include "t-sqrt.c"
#include "t-term_content.c"
#include "t-total_degree.c"
#include "t-univar.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(fmpq_mpoly_add_sub),
    TEST_FUNCTION(fmpq_mpoly_add_sub_fmpq),
    TEST_FUNCTION(fmpq_mpoly_add_sub_fmpz),
    TEST_FUNCTION(fmpq_mpoly_cmp),
    TEST_FUNCTION(fmpq_mpoly_compose_fmpq_mpoly),
    TEST_FUNCTION(fmpq_mpoly_compose_fmpq_poly),
    TEST_FUNCTION(fmpq_mpoly_content),
    TEST_FUNCTION(fmpq_mpoly_content_vars),
    TEST_FUNCTION(fmpq_mpoly_degree),
    TEST_FUNCTION(fmpq_mpoly_degrees_term_exp_fits_ui_si),
    TEST_FUNCTION(fmpq_mpoly_derivative_integral),
    TEST_FUNCTION(fmpq_mpoly_div),
    TEST_FUNCTION(fmpq_mpoly_divides),
    TEST_FUNCTION(fmpq_mpoly_divrem),
    TEST_FUNCTION(fmpq_mpoly_divrem_ideal),
    TEST_FUNCTION(fmpq_mpoly_equal_is_fmpq),
    TEST_FUNCTION(fmpq_mpoly_evaluate),
    TEST_FUNCTION(fmpq_mpoly_gcd_brown),
    TEST_FUNCTION(fmpq_mpoly_gcd),
    TEST_FUNCTION(fmpq_mpoly_gcd_cofactors),
    TEST_FUNCTION(fmpq_mpoly_gcd_hensel),
    TEST_FUNCTION(fmpq_mpoly_gcd_subresultant),
    TEST_FUNCTION(fmpq_mpoly_gcd_zippel2),
    TEST_FUNCTION(fmpq_mpoly_gcd_zippel),
    TEST_FUNCTION(fmpq_mpoly_gen),
    TEST_FUNCTION(fmpq_mpoly_get_coeff_vars_ui),
    TEST_FUNCTION(fmpq_mpoly_get_set_coeff_fmpq_fmpz),
    TEST_FUNCTION(fmpq_mpoly_get_set_coeff_fmpq_monomial),
    TEST_FUNCTION(fmpq_mpoly_get_set_coeff_fmpq_ui),
    TEST_FUNCTION(fmpq_mpoly_get_set_is_fmpq),
    TEST_FUNCTION(fmpq_mpoly_get_set_str_pretty),
    TEST_FUNCTION(fmpq_mpoly_get_set_term_coeff_fmpq),
    TEST_FUNCTION(fmpq_mpoly_get_set_term_exp_fmpz),
    TEST_FUNCTION(fmpq_mpoly_get_set_term_exp_si),
    TEST_FUNCTION(fmpq_mpoly_get_set_term_exp_ui),
    TEST_FUNCTION(fmpq_mpoly_get_term),
    TEST_FUNCTION(fmpq_mpoly_get_term_monomial),
    TEST_FUNCTION(fmpq_mpoly_mul),
    TEST_FUNCTION(fmpq_mpoly_pow_fmpz),
    TEST_FUNCTION(fmpq_mpoly_pow_ui),
    TEST_FUNCTION(fmpq_mpoly_push_term_fmpq_fmpz),
    TEST_FUNCTION(fmpq_mpoly_push_term_fmpq_ui),
    TEST_FUNCTION(fmpq_mpoly_resultant_discriminant),
    TEST_FUNCTION(fmpq_mpoly_scalar_mul_div_fmpq),
    TEST_FUNCTION(fmpq_mpoly_scalar_mul_div_fmpz),
    TEST_FUNCTION(fmpq_mpoly_scalar_mul_fmpq),
    TEST_FUNCTION(fmpq_mpoly_sqrt),
    TEST_FUNCTION(fmpq_mpoly_term_content),
    TEST_FUNCTION(fmpq_mpoly_total_degree),
    TEST_FUNCTION(fmpq_mpoly_univar)
};

/* main function *************************************************************/

TEST_MAIN(tests)
