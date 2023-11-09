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
#include "t-cmp.c"
#include "t-degree.c"
#include "t-degrees_term_exp_fits_ui_si.c"
#include "t-derivative.c"
#include "t-divides.c"
#include "t-divides_dense.c"
#include "t-divides_monagan_pearce.c"
#include "t-div_monagan_pearce.c"
#include "t-divrem.c"
#include "t-divrem_ideal_monagan_pearce.c"
#include "t-evaluate.c"
#include "t-gcd_brown.c"
#include "t-gcd_cofactors.c"
#include "t-gcd_hensel.c"
#include "t-gcd_subresultant.c"
#include "t-gcd_zippel2.c"
#include "t-gcd_zippel.c"
#include "t-gen.c"
#include "t-get_coeff_vars_ui.c"
#include "t-get_set_coeff_fmpz_fmpz.c"
#include "t-get_set_coeff_fmpz_monomial.c"
#include "t-get_set_coeff_fmpz_ui.c"
#include "t-get_set_is_fmpz.c"
#include "t-get_set_str_pretty.c"
#include "t-get_set_term_coeff_fmpz.c"
#include "t-get_set_term_exp_fmpz.c"
#include "t-get_set_term_exp_si.c"
#include "t-get_set_term_exp_ui.c"
#include "t-get_term.c"
#include "t-get_term_monomial.c"
#include "t-mul.c"
#include "t-mul_dense.c"
#include "t-mul_johnson.c"
#include "t-push_term_fmpz_fmpz.c"
#include "t-push_term_fmpz_ui.c"
#include "t-quadratic_root.c"
#include "t-resultant_discriminant.c"
#include "t-scalar_addmul_fmpz.c"
#include "t-scalar_mul_fmpz.c"
#include "t-sqrt.c"
#include "t-total_degree.c"
#include "t-univar_resultant.c"
#include "t-used_vars.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(fmpz_mod_mpoly_add_sub),
    TEST_FUNCTION(fmpz_mod_mpoly_add_sub_fmpz),
    TEST_FUNCTION(fmpz_mod_mpoly_add_sub_si),
    TEST_FUNCTION(fmpz_mod_mpoly_cmp),
    TEST_FUNCTION(fmpz_mod_mpoly_degree),
    TEST_FUNCTION(fmpz_mod_mpoly_degrees_term_exp_fits_ui_si),
    TEST_FUNCTION(fmpz_mod_mpoly_derivative),
    TEST_FUNCTION(fmpz_mod_mpoly_divides),
    TEST_FUNCTION(fmpz_mod_mpoly_divides_dense),
    TEST_FUNCTION(fmpz_mod_mpoly_divides_monagan_pearce),
    TEST_FUNCTION(fmpz_mod_mpoly_div_monagan_pearce),
    TEST_FUNCTION(fmpz_mod_mpoly_divrem),
    TEST_FUNCTION(fmpz_mod_mpoly_divrem_ideal_monagan_pearce),
    TEST_FUNCTION(fmpz_mod_mpoly_evaluate),
    TEST_FUNCTION(fmpz_mod_mpoly_gcd_brown),
    TEST_FUNCTION(fmpz_mod_mpoly_gcd_cofactors),
    TEST_FUNCTION(fmpz_mod_mpoly_gcd_hensel),
    TEST_FUNCTION(fmpz_mod_mpoly_gcd_subresultant),
    TEST_FUNCTION(fmpz_mod_mpoly_gcd_zippel2),
    TEST_FUNCTION(fmpz_mod_mpoly_gcd_zippel),
    TEST_FUNCTION(fmpz_mod_mpoly_gen),
    TEST_FUNCTION(fmpz_mod_mpoly_get_coeff_vars_ui),
    TEST_FUNCTION(fmpz_mod_mpoly_get_set_coeff_fmpz_fmpz),
    TEST_FUNCTION(fmpz_mod_mpoly_get_set_coeff_fmpz_monomial),
    TEST_FUNCTION(fmpz_mod_mpoly_get_set_coeff_fmpz_ui),
    TEST_FUNCTION(fmpz_mod_mpoly_get_set_is_fmpz),
    TEST_FUNCTION(fmpz_mod_mpoly_get_set_str_pretty),
    TEST_FUNCTION(fmpz_mod_mpoly_get_set_term_coeff_fmpz),
    TEST_FUNCTION(fmpz_mod_mpoly_get_set_term_exp_fmpz),
    TEST_FUNCTION(fmpz_mod_mpoly_get_set_term_exp_si),
    TEST_FUNCTION(fmpz_mod_mpoly_get_set_term_exp_ui),
    TEST_FUNCTION(fmpz_mod_mpoly_get_term),
    TEST_FUNCTION(fmpz_mod_mpoly_get_term_monomial),
    TEST_FUNCTION(fmpz_mod_mpoly_mul),
    TEST_FUNCTION(fmpz_mod_mpoly_mul_dense),
    TEST_FUNCTION(fmpz_mod_mpoly_mul_johnson),
    TEST_FUNCTION(fmpz_mod_mpoly_push_term_fmpz_fmpz),
    TEST_FUNCTION(fmpz_mod_mpoly_push_term_fmpz_ui),
    TEST_FUNCTION(fmpz_mod_mpoly_quadratic_root),
    TEST_FUNCTION(fmpz_mod_mpoly_resultant_discriminant),
    TEST_FUNCTION(fmpz_mod_mpoly_scalar_addmul_fmpz),
    TEST_FUNCTION(fmpz_mod_mpoly_scalar_mul_fmpz),
    TEST_FUNCTION(fmpz_mod_mpoly_sqrt),
    TEST_FUNCTION(fmpz_mod_mpoly_total_degree),
    TEST_FUNCTION(fmpz_mod_mpoly_univar_resultant),
    TEST_FUNCTION(fmpz_mod_mpoly_used_vars)
};

/* main function *************************************************************/

TEST_MAIN(tests)
